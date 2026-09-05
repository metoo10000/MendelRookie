from __future__ import annotations

import io
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import requests

OUT = Path("nutrition_stage/output_c_n13")
OUT.mkdir(parents=True, exist_ok=True)

CYCLES = {
    "E": {"year": "2007", "period": "discovery", "index": 1},
    "F": {"year": "2009", "period": "discovery", "index": 2},
    "G": {"year": "2011", "period": "discovery", "index": 3},
    "H": {"year": "2013", "period": "validation", "index": 4},
    "I": {"year": "2015", "period": "validation", "index": 5},
    "J": {"year": "2017", "period": "validation", "index": 6},
}

COMPONENTS = {
    "demo": "DEMO",
    "diet1": "DR1TOT",
    "diet2": "DR2TOT",
    "biochem": "BIOPRO",
    "rx": "RXQ_RX",
    "body": "BMX",
    "diabetes": "DIQ",
    "medical": "MCQ",
    "bpq": "BPQ",
}

THIAZIDE = re.compile(r"hydrochlorothiazide|chlorthalidone|chlorothiazide|indapamide|metolazone|bendroflumethiazide|methyclothiazide|polythiazide|trichlormethiazide", re.I)
LOOP = re.compile(r"furosemide|bumetanide|torsemide|ethacrynic", re.I)
K_SPARING = re.compile(r"spironolactone|eplerenone|amiloride|triamterene", re.I)
K_SUPPLEMENT = re.compile(r"potassium chloride|potassium citrate|potassium gluconate|potassium bicarbonate|potassium acetate", re.I)
RAAS = re.compile(r"\b[a-z]+pril\b|\b[a-z]+sartan\b|sacubitril|aliskiren", re.I)
PPI = re.compile(r"omeprazole|esomeprazole|lansoprazole|dexlansoprazole|pantoprazole|rabeprazole", re.I)


def xpt_url(cycle: str, stem: str) -> str:
    return f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{CYCLES[cycle]['year']}/DataFiles/{stem}_{cycle}.xpt"


def load_xpt(cycle: str, component: str, required: bool = True) -> tuple[pd.DataFrame, dict[str, object]]:
    stem = COMPONENTS[component]
    url = xpt_url(cycle, stem)
    try:
        response = requests.get(url, timeout=180)
        if response.status_code != 200:
            raise RuntimeError(f"HTTP {response.status_code}")
        frame = pd.read_sas(io.BytesIO(response.content), format="xport", encoding="latin1")
        frame.columns = [str(c).upper() for c in frame.columns]
        return frame, {
            "cycle": cycle, "component": component, "file": f"{stem}_{cycle}.xpt", "url": url,
            "required": required, "available": True, "rows": int(len(frame)), "columns": int(frame.shape[1]),
            "column_names": "|".join(frame.columns), "error": None,
        }
    except Exception as exc:  # noqa: BLE001
        if required:
            raise RuntimeError(f"Unable to load {component} in {cycle}: {exc}") from exc
        return pd.DataFrame(columns=["SEQN"]), {
            "cycle": cycle, "component": component, "file": f"{stem}_{cycle}.xpt", "url": url,
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


def rx_person_level(rx: pd.DataFrame) -> pd.DataFrame:
    if rx.empty or "SEQN" not in rx.columns:
        return pd.DataFrame(columns=["SEQN"])
    name_col = next((c for c in ["RXDDRUG", "RXD240B"] if c in rx.columns), None)
    if name_col is None:
        return pd.DataFrame(columns=["SEQN"])
    tmp = rx[["SEQN", name_col] + (["RXDDAYS"] if "RXDDAYS" in rx.columns else [])].copy()
    tmp["drug"] = tmp[name_col].fillna("").astype(str).str.lower()
    tmp["is_thiazide"] = tmp["drug"].str.contains(THIAZIDE, na=False)
    tmp["is_loop"] = tmp["drug"].str.contains(LOOP, na=False)
    tmp["is_k_sparing"] = tmp["drug"].str.contains(K_SPARING, na=False)
    tmp["is_k_supplement"] = tmp["drug"].str.contains(K_SUPPLEMENT, na=False)
    tmp["is_raas"] = tmp["drug"].str.contains(RAAS, na=False)
    tmp["is_ppi"] = tmp["drug"].str.contains(PPI, na=False)
    if "RXDDAYS" in tmp.columns:
        tmp["days"] = pd.to_numeric(tmp["RXDDAYS"], errors="coerce")
    else:
        tmp["days"] = np.nan
    tmp["wasting_days"] = np.where(tmp["is_thiazide"] | tmp["is_loop"], tmp["days"], np.nan)
    result = tmp.groupby("SEQN", as_index=False).agg(
        thiazide=("is_thiazide", "max"),
        loop=("is_loop", "max"),
        k_sparing=("is_k_sparing", "max"),
        k_supplement=("is_k_supplement", "max"),
        raas=("is_raas", "max"),
        ppi=("is_ppi", "max"),
        wasting_days=("wasting_days", "max"),
        drug_count=("drug", lambda x: int((x != "").sum())),
    )
    result["potassium_wasting_diuretic"] = result["thiazide"] | result["loop"]
    result["diuretic_class"] = np.select(
        [result["thiazide"] & ~result["loop"], result["loop"] & ~result["thiazide"], result["thiazide"] & result["loop"]],
        ["thiazide", "loop", "both"], default="none"
    )
    return result


def merge_cycle(cycle: str) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    frames: dict[str, pd.DataFrame] = {}
    audits: list[dict[str, object]] = []
    for component in COMPONENTS:
        frame, audit = load_xpt(cycle, component, required=component not in {"medical", "bpq"})
        frames[component] = frame
        audits.append(audit)

    rxp = rx_person_level(frames["rx"])
    pieces = [
        keep(frames["demo"], ["SEQN", "RIAGENDR", "RIDAGEYR", "RIDRETH1", "RIDRETH3", "RIDEXPRG", "SDMVPSU", "SDMVSTRA", "INDFMPIR", "DMDEDUC2"]),
        keep(frames["diet1"], ["SEQN", "WTDRD1", "DR1DRSTZ", "DR1TKCAL", "DR1TPROT", "DR1TPOTA", "DR1TSODI", "DR1TMAGN"]),
        keep(frames["diet2"], ["SEQN", "WTDR2D", "DR2DRSTZ", "DR2TKCAL", "DR2TPROT", "DR2TPOTA", "DR2TSODI", "DR2TMAGN"]),
        keep(frames["biochem"], ["SEQN", "LBXSKSI", "LBXSNASI", "LBXSCLSI", "LBXSCR", "LBDSCR", "LBXSAL", "LBXSC3SI"]),
        keep(frames["body"], ["SEQN", "BMXBMI", "BMXWT", "BMXHT"]),
        keep(frames["diabetes"], ["SEQN", "DIQ010"]),
        keep(frames["medical"], ["SEQN", "MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E", "MCQ160F"]),
        keep(frames["bpq"], ["SEQN", "BPQ020"]),
        rxp,
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
    merged["weight"] = first_numeric(merged, ["WTDR2D"]) / len(CYCLES)
    merged["psu_u"] = merged["cycle_index"] * 100 + first_numeric(merged, ["SDMVPSU"])
    merged["strata_u"] = merged["cycle_index"] * 1000 + first_numeric(merged, ["SDMVSTRA"])
    merged["diet1_valid"] = first_numeric(merged, ["DR1DRSTZ"]).eq(1)
    merged["diet2_valid"] = first_numeric(merged, ["DR2DRSTZ"]).eq(1)
    for nutrient in ["kcal", "protein", "potassium", "sodium", "magnesium"]:
        day1_field = {"kcal": "DR1TKCAL", "protein": "DR1TPROT", "potassium": "DR1TPOTA", "sodium": "DR1TSODI", "magnesium": "DR1TMAGN"}[nutrient]
        day2_field = day1_field.replace("DR1", "DR2")
        d1 = first_numeric(merged, [day1_field])
        d2 = first_numeric(merged, [day2_field])
        merged[f"day1_{nutrient}"] = d1
        merged[f"day2_{nutrient}"] = d2
        merged[f"mean_{nutrient}"] = (d1 + d2) / 2.0
    merged["serum_potassium"] = first_numeric(merged, ["LBXSKSI"])
    merged["serum_sodium"] = first_numeric(merged, ["LBXSNASI"])
    merged["serum_chloride"] = first_numeric(merged, ["LBXSCLSI"])
    merged["serum_bicarbonate"] = first_numeric(merged, ["LBXSC3SI"])
    merged["creatinine"] = first_numeric(merged, ["LBXSCR", "LBDSCR"])
    merged["albumin"] = first_numeric(merged, ["LBXSAL"])
    merged["bmi"] = first_numeric(merged, ["BMXBMI"])
    merged["body_weight"] = first_numeric(merged, ["BMXWT"])
    merged["height"] = first_numeric(merged, ["BMXHT"])
    merged["egfr"] = egfr_2021(merged["creatinine"], merged["age"], merged["sex"])
    merged["diabetes"] = np.where(first_numeric(merged, ["DIQ010"]).eq(1), 1.0, np.where(first_numeric(merged, ["DIQ010"]).isin([2, 3]), 0.0, np.nan))
    merged["hypertension"] = np.where(first_numeric(merged, ["BPQ020"]).eq(1), 1.0, np.where(first_numeric(merged, ["BPQ020"]).eq(2), 0.0, np.nan))
    cvd_fields = [c for c in ["MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E", "MCQ160F"] if c in merged.columns]
    if cvd_fields:
        cvd_matrix = merged[cvd_fields].apply(pd.to_numeric, errors="coerce")
        merged["cvd"] = np.where(cvd_matrix.eq(1).any(axis=1), 1.0, np.where(cvd_matrix.eq(2).all(axis=1), 0.0, np.nan))
        merged["heart_failure"] = np.where(first_numeric(merged, ["MCQ160B"]).eq(1), 1.0, np.where(first_numeric(merged, ["MCQ160B"]).eq(2), 0.0, np.nan))
    else:
        merged["cvd"] = np.nan
        merged["heart_failure"] = np.nan
    for col in ["thiazide", "loop", "k_sparing", "k_supplement", "raas", "ppi", "potassium_wasting_diuretic"]:
        if col not in merged.columns:
            merged[col] = False
        merged[col] = merged[col].fillna(False).astype(bool)
    if "diuretic_class" not in merged.columns:
        merged["diuretic_class"] = "none"
    if "wasting_days" not in merged.columns:
        merged["wasting_days"] = np.nan
    return merged, audits


def weighted_mean(x: pd.Series, w: pd.Series) -> float:
    x = pd.to_numeric(x, errors="coerce")
    ok = x.notna() & w.notna() & w.gt(0)
    return float(np.average(x.loc[ok], weights=w.loc[ok])) if ok.any() else float("nan")


def main() -> None:
    frames, audits = [], []
    for cycle in CYCLES:
        frame, audit = merge_cycle(cycle)
        frames.append(frame)
        audits.extend(audit)
    raw = pd.concat(frames, ignore_index=True, sort=False)

    flow: list[dict[str, object]] = []
    def record(step: str, frame: pd.DataFrame) -> None:
        flow.append({"step": step, "n": int(len(frame))})

    record("all_selected_cycle_records", raw)
    cohort = raw.loc[raw["age"].between(40, 79, inclusive="both")].copy()
    record("age_40_79", cohort)
    cohort = cohort.loc[cohort["pregnancy"].ne(1) | cohort["pregnancy"].isna()].copy()
    record("exclude_known_pregnancy", cohort)
    cohort = cohort.loc[cohort["potassium_wasting_diuretic"]].copy()
    record("current_potassium_wasting_diuretic_users", cohort)
    cohort = cohort.loc[~cohort["k_sparing"] & ~cohort["k_supplement"]].copy()
    record("exclude_potassium_sparing_or_prescription_potassium", cohort)
    cohort = cohort.loc[cohort["diet1_valid"] & cohort["diet2_valid"]].copy()
    record("two_valid_dietary_recalls", cohort)
    cohort = cohort.loc[
        cohort["mean_kcal"].between(600, 6000, inclusive="both")
        & cohort["mean_potassium"].between(300, 7000, inclusive="both")
        & cohort["mean_sodium"].between(300, 10000, inclusive="both")
        & cohort["mean_magnesium"].between(30, 1000, inclusive="both")
        & cohort["serum_potassium"].between(2.0, 6.5, inclusive="both")
        & cohort["egfr"].ge(30)
        & cohort["bmi"].between(15, 70, inclusive="both")
        & cohort["weight"].gt(0)
        & cohort["psu_u"].notna()
        & cohort["strata_u"].notna()
    ].copy()
    record("valid_nutrition_electrolyte_and_survey_data", cohort)

    cohort["exposure"] = np.select(
        [cohort["mean_potassium"].lt(2000), cohort["mean_potassium"].ge(2500)],
        ["low", "higher"], default="intermediate"
    )
    cohort = cohort.loc[cohort["exposure"].isin(["low", "higher"])].copy()
    record("frozen_low_vs_higher_intake_comparison", cohort)

    cohort["low_k_3_6"] = cohort["serum_potassium"].lt(3.6).astype(int)
    cohort["hypokalemia_3_5"] = cohort["serum_potassium"].lt(3.5).astype(int)
    cohort["low_normal_k_4_0"] = cohort["serum_potassium"].lt(4.0).astype(int)
    cohort["log_potassium_intake"] = np.log(cohort["mean_potassium"])
    cohort["potassium_density"] = cohort["mean_potassium"] / cohort["mean_kcal"] * 1000.0
    cohort["sodium_potassium_ratio"] = cohort["mean_sodium"] / cohort["mean_potassium"]
    cohort["chronic_diuretic"] = cohort["wasting_days"].ge(30)

    covars = ["age", "sex", "race", "bmi", "egfr", "diabetes", "hypertension", "cvd", "mean_kcal", "mean_sodium", "mean_magnesium", "diuretic_class", "raas", "pir"]
    cohort["complete_case"] = cohort[covars].notna().all(axis=1)
    record("complete_case_full_adjustment", cohort.loc[cohort["complete_case"]])

    group_rows = []
    for scope, frame in [("overall", cohort), ("discovery", cohort.loc[cohort["period"].eq("discovery")]), ("validation", cohort.loc[cohort["period"].eq("validation")])]:
        for exposure, group in frame.groupby("exposure", observed=True):
            group_rows.append({
                "scope": scope, "exposure": exposure, "n": int(len(group)),
                "low_k_3_6_events": int(group["low_k_3_6"].sum()),
                "hypokalemia_3_5_events": int(group["hypokalemia_3_5"].sum()),
                "low_normal_k_4_0_events": int(group["low_normal_k_4_0"].sum()),
                "weighted_low_k_3_6_pct": 100 * weighted_mean(group["low_k_3_6"], group["weight"]),
                "weighted_mean_serum_potassium": weighted_mean(group["serum_potassium"], group["weight"]),
                "weighted_mean_dietary_potassium": weighted_mean(group["mean_potassium"], group["weight"]),
                "thiazide_n": int(group["thiazide"].sum()), "loop_n": int(group["loop"].sum()),
            })

    missing_rows = []
    for scope, frame in [("overall", cohort)] + [(f"cycle_{c}", cohort.loc[cohort["cycle"].eq(c)]) for c in CYCLES]:
        for variable in ["mean_potassium", "mean_sodium", "mean_magnesium", "mean_kcal", "serum_potassium", "egfr", "bmi", "diabetes", "hypertension", "cvd", "pir", "wasting_days"]:
            missing_rows.append({"scope": scope, "variable": variable, "n": int(len(frame)), "missing_n": int(frame[variable].isna().sum()), "missing_pct": float(frame[variable].isna().mean() * 100) if len(frame) else np.nan})

    semantic = pd.DataFrame([
        {"concept": "dietary_potassium", "fields": "DR1TPOTA and DR2TPOTA", "meaning": "mean potassium from two complete 24-hour recalls", "unit": "mg/day", "rule": "low <2000; higher >=2500; 2000-2499 excluded"},
        {"concept": "drug_exposure", "fields": "RXDDRUG in RXQ_RX", "meaning": "generic names of prescription medicines used in past 30 days", "unit": "person-level class", "rule": "thiazide or loop; exclude potassium-sparing and prescription potassium"},
        {"concept": "outcome", "fields": "LBXSKSI", "meaning": "measured serum potassium", "unit": "mmol/L", "rule": "primary <3.6; strict <3.5; continuous secondary"},
        {"concept": "survey_design", "fields": "WTDR2D, SDMVPSU, SDMVSTRA", "meaning": "day-2 dietary weight and NHANES design", "unit": "survey", "rule": "pooled weight divided by six; cycle-unique PSU/strata"},
    ])

    columns = ["SEQN", "cycle", "period", "age", "sex", "race", "education", "pir", "weight", "psu_u", "strata_u", "mean_kcal", "mean_protein", "mean_potassium", "mean_sodium", "mean_magnesium", "day1_potassium", "day2_potassium", "potassium_density", "sodium_potassium_ratio", "serum_potassium", "serum_sodium", "serum_chloride", "serum_bicarbonate", "creatinine", "egfr", "albumin", "bmi", "body_weight", "height", "diabetes", "hypertension", "cvd", "heart_failure", "thiazide", "loop", "diuretic_class", "wasting_days", "chronic_diuretic", "raas", "ppi", "exposure", "low_k_3_6", "hypokalemia_3_5", "low_normal_k_4_0", "log_potassium_intake", "complete_case"]
    cohort.loc[:, columns].to_csv(OUT / "c_n13_analysis_core.csv", index=False)
    pd.DataFrame(audits).to_csv(OUT / "c_n13_source_audit.csv", index=False)
    pd.DataFrame(flow).to_csv(OUT / "c_n13_flow.csv", index=False)
    pd.DataFrame(group_rows).to_csv(OUT / "c_n13_group_counts.csv", index=False)
    pd.DataFrame(missing_rows).to_csv(OUT / "c_n13_missingness.csv", index=False)
    semantic.to_csv(OUT / "c_n13_semantic_audit.csv", index=False)

    key = cohort.loc[cohort["exposure"].eq("low")]
    val_key = key.loc[key["period"].eq("validation")]
    prep = {
        "candidate_code": "C-N13", "actual_n": int(len(cohort)),
        "complete_case_n": int(cohort["complete_case"].sum()), "complete_case_retention": float(cohort["complete_case"].mean()),
        "primary_events": int(cohort["low_k_3_6"].sum()), "strict_events": int(cohort["hypokalemia_3_5"].sum()),
        "key_group_n": int(len(key)), "key_group_events": int(key["low_k_3_6"].sum()),
        "validation_key_group_n": int(len(val_key)), "validation_key_group_events": int(val_key["low_k_3_6"].sum()),
        "cycles": sorted(cohort["cycle"].unique().tolist()),
    }
    (OUT / "c_n13_prep_status.json").write_text(json.dumps(prep, indent=2), encoding="utf-8")
    print(json.dumps(prep, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
