from __future__ import annotations

import io
import json
from pathlib import Path

import numpy as np
import pandas as pd
import requests

OUT = Path("nutrition_stage/output_c_n11")
OUT.mkdir(parents=True, exist_ok=True)

CYCLES = {
    "C": {"year": "2003", "period": "discovery", "index": 1},
    "D": {"year": "2005", "period": "validation", "index": 2},
}

FILES = {
    "C": {
        "demo": "DEMO_C", "vitc": "L06VIT_C", "iron": "L06TFR_C", "cbc": "L25_C",
        "crp": "L11_C", "biochem": "L40_C", "body": "BMX_C", "smoking": "SMQ_C",
        "diabetes": "DIQ_C", "medical": "MCQ_C",
    },
    "D": {
        "demo": "DEMO_D", "vitc": "VIC_D", "iron": "FERTIN_D", "cbc": "CBC_D",
        "crp": "CRP_D", "biochem": "BIOPRO_D", "body": "BMX_D", "smoking": "SMQ_D",
        "diabetes": "DIQ_D", "medical": "MCQ_D",
    },
}


def load_xpt(cycle: str, component: str, required: bool = True) -> tuple[pd.DataFrame, dict[str, object]]:
    stem = FILES[cycle][component]
    url = f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{CYCLES[cycle]['year']}/DataFiles/{stem}.xpt"
    try:
        response = requests.get(url, timeout=180)
        if response.status_code != 200:
            raise RuntimeError(f"HTTP {response.status_code}")
        frame = pd.read_sas(io.BytesIO(response.content), format="xport", encoding="latin1")
        frame.columns = [str(c).upper() for c in frame.columns]
        return frame, {
            "cycle": cycle, "component": component, "file": f"{stem}.xpt", "url": url,
            "required": required, "available": True, "rows": int(len(frame)),
            "columns": int(frame.shape[1]), "column_names": "|".join(frame.columns), "error": None,
        }
    except Exception as exc:  # noqa: BLE001
        if required:
            raise RuntimeError(f"Unable to load {component} in cycle {cycle}: {exc}") from exc
        return pd.DataFrame(columns=["SEQN"]), {
            "cycle": cycle, "component": component, "file": f"{stem}.xpt", "url": url,
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


def egfr_2021(creatinine: pd.Series, age: pd.Series, sex: pd.Series) -> pd.Series:
    female = sex.eq(2)
    kappa = np.where(female, 0.7, 0.9)
    alpha = np.where(female, -0.241, -0.302)
    ratio = creatinine / kappa
    return 142.0 * np.minimum(ratio, 1.0) ** alpha * np.maximum(ratio, 1.0) ** -1.2 * 0.9938 ** age * np.where(female, 1.012, 1.0)


def weighted_mean(x: pd.Series, w: pd.Series) -> float:
    x = pd.to_numeric(x, errors="coerce")
    ok = x.notna() & w.notna() & w.gt(0)
    return float(np.average(x.loc[ok], weights=w.loc[ok])) if ok.any() else float("nan")


def merge_cycle(cycle: str) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    frames: dict[str, pd.DataFrame] = {}
    audits: list[dict[str, object]] = []
    for component in FILES[cycle]:
        frame, audit = load_xpt(cycle, component, required=component != "medical")
        frames[component] = frame
        audits.append(audit)

    pieces = [
        keep(frames["demo"], ["SEQN", "RIAGENDR", "RIDAGEYR", "RIDRETH1", "RIDRETH3", "RIDEXPRG", "WTMEC2YR", "SDMVPSU", "SDMVSTRA", "INDFMPIR", "DMDEDUC2"]),
        keep(frames["vitc"], ["SEQN", "LBXVIC", "LBDVICSI"]),
        keep(frames["iron"], ["SEQN", "LBXFER", "LBDFER", "LBDFERSI", "LBXTFR"]),
        keep(frames["cbc"], ["SEQN", "LBXHGB", "LBXMCVSI", "LBXMCV", "LBXRDW", "LBXMCHSI"]),
        keep(frames["crp"], ["SEQN", "LBXCRP"]),
        keep(frames["biochem"], ["SEQN", "LBXSCR", "LBDSCR", "LBXSAL", "LBXSTP", "LBXSATSI", "LBXAST"]),
        keep(frames["body"], ["SEQN", "BMXBMI", "BMXWT", "BMXHT"]),
        keep(frames["smoking"], ["SEQN", "SMQ020", "SMQ040"]),
        keep(frames["diabetes"], ["SEQN", "DIQ010"]),
        keep(frames["medical"], ["SEQN", "MCQ220", "MCQ160L"]),
    ]
    merged = pieces[0]
    for piece in pieces[1:]:
        if not piece.empty and "SEQN" in piece.columns:
            merged = merged.merge(piece.drop_duplicates("SEQN"), on="SEQN", how="left")

    merged["cycle"] = cycle
    merged["period"] = CYCLES[cycle]["period"]
    merged["cycle_index"] = CYCLES[cycle]["index"]
    merged["age"] = first_numeric(merged, ["RIDAGEYR"])
    merged["sex"] = first_numeric(merged, ["RIAGENDR"])
    merged["race"] = first_numeric(merged, ["RIDRETH3", "RIDRETH1"])
    merged["pregnancy"] = first_numeric(merged, ["RIDEXPRG"])
    merged["education"] = first_numeric(merged, ["DMDEDUC2"])
    merged["pir"] = first_numeric(merged, ["INDFMPIR"])
    merged["weight"] = first_numeric(merged, ["WTMEC2YR"]) / 2.0
    merged["psu_u"] = merged["cycle_index"] * 100 + first_numeric(merged, ["SDMVPSU"])
    merged["strata_u"] = merged["cycle_index"] * 1000 + first_numeric(merged, ["SDMVSTRA"])
    merged["vitamin_c_umol_l"] = first_numeric(merged, ["LBDVICSI"])
    merged.loc[merged["vitamin_c_umol_l"].isna(), "vitamin_c_umol_l"] = first_numeric(merged, ["LBXVIC"]) * 56.78
    merged["ferritin"] = first_numeric(merged, ["LBXFER", "LBDFER", "LBDFERSI"])
    merged["stfr"] = first_numeric(merged, ["LBXTFR"])
    merged["hgb"] = first_numeric(merged, ["LBXHGB"])
    merged["mcv"] = first_numeric(merged, ["LBXMCVSI", "LBXMCV"])
    merged["rdw"] = first_numeric(merged, ["LBXRDW"])
    merged["mch"] = first_numeric(merged, ["LBXMCHSI"])
    merged["crp_mg_l"] = first_numeric(merged, ["LBXCRP"]) * 10.0
    merged["creatinine"] = first_numeric(merged, ["LBXSCR", "LBDSCR"])
    merged["albumin"] = first_numeric(merged, ["LBXSAL"])
    merged["total_protein"] = first_numeric(merged, ["LBXSTP"])
    merged["alt"] = first_numeric(merged, ["LBXSATSI"])
    merged["ast"] = first_numeric(merged, ["LBXAST"])
    merged["bmi"] = first_numeric(merged, ["BMXBMI"])
    merged["body_weight"] = first_numeric(merged, ["BMXWT"])
    merged["height"] = first_numeric(merged, ["BMXHT"])
    merged["smoked100"] = first_numeric(merged, ["SMQ020"])
    merged["current_smoke"] = first_numeric(merged, ["SMQ040"])
    merged["diabetes_q"] = first_numeric(merged, ["DIQ010"])
    merged["cancer_q"] = first_numeric(merged, ["MCQ220"])
    merged["liver_q"] = first_numeric(merged, ["MCQ160L"])
    merged["egfr"] = egfr_2021(merged["creatinine"], merged["age"], merged["sex"])
    return merged, audits


def main() -> None:
    frames, audits = [], []
    for cycle in CYCLES:
        frame, audit = merge_cycle(cycle)
        frames.append(frame)
        audits.extend(audit)
    raw = pd.concat(frames, ignore_index=True, sort=False)

    raw["smoking"] = pd.Series(pd.NA, index=raw.index, dtype="string")
    raw.loc[raw["smoked100"].eq(2), "smoking"] = "never"
    raw.loc[raw["smoked100"].eq(1) & raw["current_smoke"].eq(3), "smoking"] = "former"
    raw.loc[raw["smoked100"].eq(1) & raw["current_smoke"].isin([1, 2]), "smoking"] = "current"
    raw["diabetes"] = np.where(raw["diabetes_q"].eq(1), 1.0, np.where(raw["diabetes_q"].isin([2, 3]), 0.0, np.nan))
    raw["known_cancer"] = np.where(raw["cancer_q"].eq(1), 1.0, np.where(raw["cancer_q"].eq(2), 0.0, np.nan))
    raw["known_liver"] = np.where(raw["liver_q"].eq(1), 1.0, np.where(raw["liver_q"].eq(2), 0.0, np.nan))

    flow = []
    def record(step: str, frame: pd.DataFrame) -> None:
        flow.append({"step": step, "n": int(len(frame))})

    record("all_selected_cycle_records", raw)
    cohort = raw.loc[raw["age"].between(20, 79, inclusive="both")].copy()
    record("age_20_79", cohort)
    cohort = cohort.loc[~cohort["pregnancy"].eq(1)].copy()
    record("exclude_known_pregnancy", cohort)
    cohort = cohort.loc[
        cohort["vitamin_c_umol_l"].between(2, 200, inclusive="both")
        & cohort["ferritin"].between(30, 500, inclusive="both")
        & cohort["hgb"].between(5, 20, inclusive="both")
        & cohort["mcv"].between(55, 125, inclusive="both")
        & cohort["rdw"].between(8, 30, inclusive="both")
    ].copy()
    record("valid_nutrition_iron_and_hematology", cohort)
    cohort = cohort.loc[
        cohort["crp_mg_l"].between(0, 10, inclusive="both")
        & cohort["egfr"].ge(60)
        & cohort["albumin"].ge(3.5)
        & cohort["bmi"].between(15, 70, inclusive="both")
    ].copy()
    record("exclude_overt_inflammation_ckd_hypoalbuminemia", cohort)
    cohort = cohort.loc[~cohort["known_cancer"].eq(1) & ~cohort["known_liver"].eq(1)].copy()
    record("exclude_known_cancer_or_liver_disease", cohort)
    cohort = cohort.loc[cohort["weight"].gt(0) & cohort["psu_u"].notna() & cohort["strata_u"].notna()].copy()
    record("valid_survey_design", cohort)

    cohort["exposure"] = np.select(
        [cohort["vitamin_c_umol_l"].lt(23), cohort["vitamin_c_umol_l"].ge(50)],
        ["depleted", "adequate"], default="intermediate"
    )
    cohort = cohort.loc[cohort["exposure"].isin(["depleted", "adequate"])].copy()
    record("frozen_exposure_contrast", cohort)

    cohort["anemia"] = np.where(cohort["sex"].eq(1), cohort["hgb"].lt(13), cohort["hgb"].lt(12)).astype(int)
    cohort["high_rdw"] = cohort["rdw"].gt(14.5).astype(int)
    cohort["microcytosis"] = cohort["mcv"].lt(80).astype(int)
    cohort["log_vitamin_c"] = np.log(cohort["vitamin_c_umol_l"])
    cohort["log_crp"] = np.log(cohort["crp_mg_l"].clip(lower=0.05))
    cohort["log_ferritin"] = np.log(cohort["ferritin"])
    full_covars = ["age", "sex", "race", "bmi", "log_crp", "log_ferritin", "albumin", "smoking", "diabetes", "pir"]
    cohort["complete_case"] = cohort[full_covars].notna().all(axis=1)
    record("complete_case_full_adjustment", cohort.loc[cohort["complete_case"]])

    group_rows = []
    for scope, frame in [("overall", cohort), ("discovery", cohort.loc[cohort["period"].eq("discovery")]), ("validation", cohort.loc[cohort["period"].eq("validation")])]:
        for exposure, group in frame.groupby("exposure", observed=True):
            group_rows.append({
                "scope": scope, "exposure": exposure, "n": int(len(group)),
                "anemia_events": int(group["anemia"].sum()), "high_rdw_events": int(group["high_rdw"].sum()),
                "microcytosis_events": int(group["microcytosis"].sum()),
                "weighted_anemia_pct": 100 * weighted_mean(group["anemia"], group["weight"]),
                "weighted_mean_hgb": weighted_mean(group["hgb"], group["weight"]),
                "weighted_mean_rdw": weighted_mean(group["rdw"], group["weight"]),
                "weighted_mean_vitamin_c": weighted_mean(group["vitamin_c_umol_l"], group["weight"]),
            })

    missing_rows = []
    for scope, frame in [("overall", cohort)] + [(f"cycle_{c}", cohort.loc[cohort["cycle"].eq(c)]) for c in CYCLES]:
        for variable in ["vitamin_c_umol_l", "ferritin", "hgb", "mcv", "rdw", "crp_mg_l", "egfr", "albumin", "bmi", "smoking", "diabetes", "pir"]:
            missing_rows.append({"scope": scope, "variable": variable, "n": int(len(frame)), "missing_n": int(frame[variable].isna().sum()), "missing_pct": float(frame[variable].isna().mean() * 100) if len(frame) else np.nan})

    semantic = pd.DataFrame([
        {"concept": "vitamin_c", "fields": "LBDVICSI or LBXVIC converted", "unit": "umol/L", "meaning": "measured serum vitamin C", "rule": "valid 2-200; depleted <23; adequate >=50; intermediate excluded"},
        {"concept": "iron_replete_restriction", "fields": "LBXFER/LBDFER/LBDFERSI", "unit": "ng/mL", "meaning": "measured serum ferritin", "rule": "30-500"},
        {"concept": "hematologic_outcomes", "fields": "LBXHGB/LBXMCVSI/LBXRDW", "unit": "g/dL, fL, percent", "meaning": "measured hemoglobin, MCV and RDW", "rule": "sex-specific anemia; RDW >14.5 supportive"},
        {"concept": "survey_design", "fields": "WTMEC2YR/SDMVPSU/SDMVSTRA", "unit": "survey", "meaning": "MEC weight, PSU and strata", "rule": "pooled weight divided by two; cycle-unique IDs"},
    ])

    columns = ["SEQN", "cycle", "period", "age", "sex", "race", "education", "pir", "pregnancy", "weight", "psu_u", "strata_u", "vitamin_c_umol_l", "ferritin", "stfr", "hgb", "mcv", "rdw", "mch", "crp_mg_l", "creatinine", "egfr", "albumin", "total_protein", "alt", "ast", "bmi", "body_weight", "height", "smoking", "diabetes", "known_cancer", "known_liver", "exposure", "anemia", "high_rdw", "microcytosis", "log_vitamin_c", "log_crp", "log_ferritin", "complete_case"]
    cohort.loc[:, columns].to_csv(OUT / "c_n11_analysis_core.csv", index=False)
    pd.DataFrame(audits).to_csv(OUT / "c_n11_source_audit.csv", index=False)
    pd.DataFrame(flow).to_csv(OUT / "c_n11_flow.csv", index=False)
    pd.DataFrame(group_rows).to_csv(OUT / "c_n11_group_counts.csv", index=False)
    pd.DataFrame(missing_rows).to_csv(OUT / "c_n11_missingness.csv", index=False)
    semantic.to_csv(OUT / "c_n11_semantic_audit.csv", index=False)

    key = cohort.loc[cohort["exposure"].eq("depleted")]
    val_key = key.loc[key["period"].eq("validation")]
    status = {
        "candidate_code": "C-N11", "actual_n": int(len(cohort)),
        "complete_case_n": int(cohort["complete_case"].sum()), "complete_case_retention": float(cohort["complete_case"].mean()),
        "anemia_events": int(cohort["anemia"].sum()), "key_group_n": int(len(key)),
        "key_group_anemia_events": int(key["anemia"].sum()), "validation_key_group_n": int(len(val_key)),
        "validation_key_group_anemia_events": int(val_key["anemia"].sum()), "cycles": sorted(cohort["cycle"].unique().tolist()),
    }
    (OUT / "c_n11_prep_status.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
    print(json.dumps(status, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
