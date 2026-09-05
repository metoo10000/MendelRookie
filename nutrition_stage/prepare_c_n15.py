from __future__ import annotations

import io
import json
from pathlib import Path

import numpy as np
import pandas as pd
import requests

from prepare_c_n13 import CYCLES, merge_cycle, weighted_mean

OUT = Path("nutrition_stage/output_c_n15")
OUT.mkdir(parents=True, exist_ok=True)


def load_extra(cycle: str, stem: str, required: bool = True) -> tuple[pd.DataFrame, dict[str, object]]:
    year = CYCLES[cycle]["year"]
    url = f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{year}/DataFiles/{stem}_{cycle}.xpt"
    try:
        response = requests.get(url, timeout=180)
        if response.status_code != 200:
            raise RuntimeError(f"HTTP {response.status_code}")
        frame = pd.read_sas(io.BytesIO(response.content), format="xport", encoding="latin1")
        frame.columns = [str(c).upper() for c in frame.columns]
        return frame, {
            "cycle": cycle, "component": stem, "file": f"{stem}_{cycle}.xpt", "url": url,
            "required": required, "available": True, "rows": int(len(frame)),
            "columns": int(frame.shape[1]), "column_names": "|".join(frame.columns), "error": None,
        }
    except Exception as exc:  # noqa: BLE001
        if required:
            raise RuntimeError(f"Unable to load {stem}_{cycle}: {exc}") from exc
        return pd.DataFrame(columns=["SEQN"]), {
            "cycle": cycle, "component": stem, "file": f"{stem}_{cycle}.xpt", "url": url,
            "required": required, "available": False, "rows": 0, "columns": 0,
            "column_names": "SEQN", "error": f"{type(exc).__name__}: {exc}",
        }


def keep(frame: pd.DataFrame, names: list[str]) -> pd.DataFrame:
    cols = [name for name in names if name in frame.columns]
    if "SEQN" in frame.columns and "SEQN" not in cols:
        cols.insert(0, "SEQN")
    return frame.loc[:, cols].copy()


def first_numeric(frame: pd.DataFrame, names: list[str]) -> pd.Series:
    out = pd.Series(np.nan, index=frame.index, dtype="float64")
    for name in names:
        if name in frame.columns:
            out = out.fillna(pd.to_numeric(frame[name], errors="coerce"))
    return out


def smoking_status(frame: pd.DataFrame) -> pd.Series:
    ever = first_numeric(frame, ["SMQ020"])
    now = first_numeric(frame, ["SMQ040"])
    result = pd.Series(pd.NA, index=frame.index, dtype="string")
    result.loc[ever.eq(2)] = "never"
    result.loc[ever.eq(1) & now.eq(3)] = "former"
    result.loc[ever.eq(1) & now.isin([1, 2])] = "current"
    return result


def mean_sbp(frame: pd.DataFrame) -> pd.Series:
    cols = [c for c in ["BPXSY1", "BPXSY2", "BPXSY3", "BPXSY4"] if c in frame.columns]
    if not cols:
        return pd.Series(np.nan, index=frame.index, dtype="float64")
    values = frame[cols].apply(pd.to_numeric, errors="coerce")
    values = values.where(values.between(60, 260))
    return values.mean(axis=1, skipna=True)


def main() -> None:
    frames: list[pd.DataFrame] = []
    audits: list[dict[str, object]] = []

    for cycle in CYCLES:
        base, base_audit = merge_cycle(cycle)
        audits.extend(base_audit)
        urine, a1 = load_extra(cycle, "ALB_CR")
        ghb, a2 = load_extra(cycle, "GHB")
        bpx, a3 = load_extra(cycle, "BPX")
        smq, a4 = load_extra(cycle, "SMQ")
        audits.extend([a1, a2, a3, a4])

        extra = keep(urine, ["SEQN", "URDACT", "URXUMA", "URXUCR"])
        extra = extra.merge(keep(ghb, ["SEQN", "LBXGH"]), on="SEQN", how="outer")
        bpx_small = keep(bpx, ["SEQN", "BPXSY1", "BPXSY2", "BPXSY3", "BPXSY4"])
        bpx_small["sbp_mean"] = mean_sbp(bpx_small)
        extra = extra.merge(bpx_small[["SEQN", "sbp_mean"]], on="SEQN", how="outer")
        smq_small = keep(smq, ["SEQN", "SMQ020", "SMQ040"])
        smq_small["smoking"] = smoking_status(smq_small)
        extra = extra.merge(smq_small[["SEQN", "smoking"]], on="SEQN", how="outer")
        base = base.merge(extra.drop_duplicates("SEQN"), on="SEQN", how="left")
        base["uacr"] = first_numeric(base, ["URDACT"])
        base["hba1c"] = first_numeric(base, ["LBXGH"])
        base["sbp"] = first_numeric(base, ["sbp_mean"])
        frames.append(base)

    raw = pd.concat(frames, ignore_index=True, sort=False)
    flow: list[dict[str, object]] = []

    def record(step: str, frame: pd.DataFrame) -> None:
        flow.append({"step": step, "n": int(len(frame))})

    record("all_selected_cycle_records", raw)
    cohort = raw.loc[raw["age"].between(40, 79, inclusive="both")].copy()
    record("age_40_79", cohort)
    cohort = cohort.loc[cohort["pregnancy"].ne(1) | cohort["pregnancy"].isna()].copy()
    record("exclude_known_pregnancy", cohort)
    cohort = cohort.loc[cohort["raas"]].copy()
    record("current_acei_or_arb_users", cohort)
    cohort = cohort.loc[cohort["diet1_valid"] & cohort["diet2_valid"]].copy()
    record("two_valid_dietary_recalls", cohort)
    cohort = cohort.loc[
        cohort["mean_kcal"].between(600, 6000, inclusive="both")
        & cohort["mean_sodium"].between(300, 10000, inclusive="both")
        & cohort["mean_potassium"].between(300, 7000, inclusive="both")
        & cohort["mean_protein"].between(10, 300, inclusive="both")
        & cohort["uacr"].between(0.1, 5000, inclusive="both")
        & cohort["egfr"].ge(30)
        & cohort["bmi"].between(15, 70, inclusive="both")
        & cohort["weight"].gt(0)
        & cohort["psu_u"].notna()
        & cohort["strata_u"].notna()
    ].copy()
    record("valid_nutrition_kidney_and_survey_data", cohort)

    cohort["exposure"] = np.select(
        [cohort["mean_sodium"].lt(2300), cohort["mean_sodium"].ge(3000)],
        ["lower", "high"], default="intermediate"
    )
    cohort = cohort.loc[cohort["exposure"].isin(["lower", "high"])].copy()
    record("frozen_lower_vs_high_sodium_comparison", cohort)

    cohort["albuminuria"] = cohort["uacr"].ge(30).astype(int)
    cohort["severe_albuminuria"] = cohort["uacr"].ge(300).astype(int)
    cohort["log_uacr"] = np.log(cohort["uacr"])
    cohort["log_sodium"] = np.log(cohort["mean_sodium"])
    cohort["sodium_density"] = cohort["mean_sodium"] / cohort["mean_kcal"] * 1000.0
    cohort["any_diuretic"] = cohort["thiazide"] | cohort["loop"]

    covars = [
        "age", "sex", "race", "bmi", "egfr", "diabetes", "hypertension", "cvd", "sbp", "hba1c",
        "mean_kcal", "mean_protein", "mean_potassium", "any_diuretic", "smoking", "pir",
    ]
    cohort["complete_case"] = cohort[covars].notna().all(axis=1)
    record("complete_case_full_adjustment", cohort.loc[cohort["complete_case"]])

    group_rows: list[dict[str, object]] = []
    for scope, frame in [
        ("overall", cohort),
        ("discovery", cohort.loc[cohort["period"].eq("discovery")]),
        ("validation", cohort.loc[cohort["period"].eq("validation")]),
    ]:
        for exposure, group in frame.groupby("exposure", observed=True):
            group_rows.append({
                "scope": scope, "exposure": exposure, "n": int(len(group)),
                "albuminuria_events": int(group["albuminuria"].sum()),
                "severe_albuminuria_events": int(group["severe_albuminuria"].sum()),
                "weighted_albuminuria_pct": 100 * weighted_mean(group["albuminuria"], group["weight"]),
                "weighted_geometric_mean_uacr": float(np.exp(weighted_mean(group["log_uacr"], group["weight"]))),
                "weighted_mean_sodium": weighted_mean(group["mean_sodium"], group["weight"]),
                "weighted_mean_sbp": weighted_mean(group["sbp"], group["weight"]),
            })

    missing_rows: list[dict[str, object]] = []
    for scope, frame in [("overall", cohort)] + [(f"cycle_{c}", cohort.loc[cohort["cycle"].eq(c)]) for c in CYCLES]:
        for variable in ["mean_sodium", "mean_kcal", "mean_protein", "mean_potassium", "uacr", "egfr", "bmi", "diabetes", "hypertension", "cvd", "sbp", "hba1c", "smoking", "pir"]:
            missing_rows.append({
                "scope": scope, "variable": variable, "n": int(len(frame)),
                "missing_n": int(frame[variable].isna().sum()),
                "missing_pct": float(frame[variable].isna().mean() * 100) if len(frame) else np.nan,
            })

    semantic = pd.DataFrame([
        {"concept": "dietary_sodium", "fields": "DR1TSODI and DR2TSODI", "meaning": "mean sodium from two complete 24-hour recalls", "unit": "mg/day", "rule": "lower <2300; high >=3000; middle excluded"},
        {"concept": "raas_therapy", "fields": "RXDDRUG generic names", "meaning": "current ACE inhibitor or ARB prescription use in past 30 days", "unit": "person-level class", "rule": "generic name ending in -pril or -sartan plus sacubitril/aliskiren mapping"},
        {"concept": "kidney_outcome", "fields": "URDACT", "meaning": "spot urine albumin-creatinine ratio", "unit": "mg/g", "rule": "primary >=30; severe >=300; log-UACR continuous"},
        {"concept": "blood_pressure", "fields": "BPXSY1-BPXSY4", "meaning": "mean of valid measured systolic blood pressure readings", "unit": "mmHg", "rule": "valid 60-260"},
        {"concept": "survey_design", "fields": "WTDR2D, SDMVPSU, SDMVSTRA", "meaning": "day-2 dietary weight and NHANES complex design", "unit": "survey", "rule": "pooled weight divided by six; cycle-unique PSU/strata"},
    ])

    columns = [
        "SEQN", "cycle", "period", "age", "sex", "race", "education", "pir", "weight", "psu_u", "strata_u",
        "mean_kcal", "mean_protein", "mean_potassium", "mean_sodium", "mean_magnesium", "day1_sodium", "day2_sodium",
        "sodium_density", "uacr", "albuminuria", "severe_albuminuria", "log_uacr", "log_sodium", "creatinine", "egfr",
        "albumin", "bmi", "body_weight", "height", "diabetes", "hypertension", "cvd", "heart_failure", "sbp", "hba1c",
        "thiazide", "loop", "any_diuretic", "raas", "ppi", "smoking", "exposure", "complete_case",
    ]
    cohort.loc[:, columns].to_csv(OUT / "c_n15_analysis_core.csv", index=False)
    pd.DataFrame(audits).to_csv(OUT / "c_n15_source_audit.csv", index=False)
    pd.DataFrame(flow).to_csv(OUT / "c_n15_flow.csv", index=False)
    pd.DataFrame(group_rows).to_csv(OUT / "c_n15_group_counts.csv", index=False)
    pd.DataFrame(missing_rows).to_csv(OUT / "c_n15_missingness.csv", index=False)
    semantic.to_csv(OUT / "c_n15_semantic_audit.csv", index=False)

    key = cohort.loc[cohort["exposure"].eq("high")]
    val_key = key.loc[key["period"].eq("validation")]
    prep = {
        "candidate_code": "C-N15", "actual_n": int(len(cohort)),
        "complete_case_n": int(cohort["complete_case"].sum()),
        "complete_case_retention": float(cohort["complete_case"].mean()),
        "primary_events": int(cohort["albuminuria"].sum()),
        "severe_events": int(cohort["severe_albuminuria"].sum()),
        "key_group_n": int(len(key)), "key_group_events": int(key["albuminuria"].sum()),
        "validation_key_group_n": int(len(val_key)), "validation_key_group_events": int(val_key["albuminuria"].sum()),
        "cycles": sorted(cohort["cycle"].unique().tolist()),
    }
    (OUT / "c_n15_prep_status.json").write_text(json.dumps(prep, indent=2), encoding="utf-8")
    print(json.dumps(prep, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
