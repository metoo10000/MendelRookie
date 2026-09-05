from __future__ import annotations

import io
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import requests

OUT = Path("nutrition_stage/output_c_n09")
OUT.mkdir(parents=True, exist_ok=True)

CYCLES = {
    "D": {"year": "2005", "years": "2005-2006", "index": 1, "period": "discovery"},
    "E": {"year": "2007", "years": "2007-2008", "index": 2, "period": "discovery"},
    "F": {"year": "2009", "years": "2009-2010", "index": 3, "period": "validation"},
}

COMPONENTS = {
    "demo": "DEMO",
    "b6": "VIT_B6",
    "ferritin": "FERTIN",
    "cbc": "CBC",
    "crp": "CRP",
    "biochem": "BIOPRO",
    "body": "BMX",
    "smoking": "SMQ",
    "diabetes": "DIQ",
    "medical": "MCQ",
}


def xpt_url(cycle: str, stem: str) -> str:
    year = CYCLES[cycle]["year"]
    return f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{year}/DataFiles/{stem}_{cycle}.xpt"


def load_xpt(cycle: str, component: str, required: bool = True) -> tuple[pd.DataFrame, dict[str, object]]:
    stem = COMPONENTS[component]
    url = xpt_url(cycle, stem)
    try:
        response = requests.get(url, timeout=180)
        if response.status_code != 200:
            raise RuntimeError(f"HTTP {response.status_code}")
        frame = pd.read_sas(io.BytesIO(response.content), format="xport", encoding="latin1")
        frame.columns = [str(c).upper() for c in frame.columns]
        audit = {
            "cycle": cycle,
            "component": component,
            "file": f"{stem}_{cycle}.xpt",
            "url": url,
            "required": required,
            "available": True,
            "rows": int(frame.shape[0]),
            "columns": int(frame.shape[1]),
            "column_names": "|".join(frame.columns),
            "error": None,
        }
        return frame, audit
    except Exception as exc:  # noqa: BLE001
        if required:
            raise RuntimeError(f"Unable to load {component} for cycle {cycle}: {exc}") from exc
        return pd.DataFrame(columns=["SEQN"]), {
            "cycle": cycle,
            "component": component,
            "file": f"{stem}_{cycle}.xpt",
            "url": url,
            "required": required,
            "available": False,
            "rows": 0,
            "columns": 0,
            "column_names": "SEQN",
            "error": f"{type(exc).__name__}: {exc}",
        }


def keep(frame: pd.DataFrame, columns: Iterable[str]) -> pd.DataFrame:
    selected = [c for c in columns if c in frame.columns]
    if "SEQN" in frame.columns and "SEQN" not in selected:
        selected.insert(0, "SEQN")
    return frame.loc[:, selected].copy()


def numeric(frame: pd.DataFrame, names: Iterable[str]) -> pd.Series:
    result = pd.Series(np.nan, index=frame.index, dtype="float64")
    for name in names:
        if name in frame.columns:
            result = result.fillna(pd.to_numeric(frame[name], errors="coerce"))
    return result


def egfr_2021(creatinine: pd.Series, age: pd.Series, sex: pd.Series) -> pd.Series:
    female = sex.eq(2)
    kappa = np.where(female, 0.7, 0.9)
    alpha = np.where(female, -0.241, -0.302)
    ratio = creatinine / kappa
    return (
        142.0
        * np.minimum(ratio, 1.0) ** alpha
        * np.maximum(ratio, 1.0) ** -1.200
        * 0.9938 ** age
        * np.where(female, 1.012, 1.0)
    )


def weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    values = pd.to_numeric(values, errors="coerce")
    valid = values.notna() & weights.notna() & weights.gt(0)
    if not valid.any():
        return float("nan")
    return float(np.average(values.loc[valid], weights=weights.loc[valid]))


def merge_cycle(cycle: str) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    frames: dict[str, pd.DataFrame] = {}
    audits: list[dict[str, object]] = []
    for component in COMPONENTS:
        frame, audit = load_xpt(cycle, component, required=component not in {"medical"})
        frames[component] = frame
        audits.append(audit)

    pieces = [
        keep(frames["demo"], [
            "SEQN", "RIAGENDR", "RIDAGEYR", "RIDRETH1", "RIDRETH3", "RIDEXPRG",
            "WTMEC2YR", "SDMVPSU", "SDMVSTRA", "INDFMPIR", "DMDEDUC2",
        ]),
        keep(frames["b6"], ["SEQN", "LBXPLP", "LBX4PA"]),
        keep(frames["ferritin"], ["SEQN", "LBXFER", "LBDFERSI"]),
        keep(frames["cbc"], ["SEQN", "LBXHGB", "LBXMCVSI", "LBXMCV", "LBXMCHSI", "LBXRDW"]),
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
    merged["cycle_years"] = CYCLES[cycle]["years"]
    merged["period"] = CYCLES[cycle]["period"]
    merged["cycle_index"] = CYCLES[cycle]["index"]
    merged["age"] = numeric(merged, ["RIDAGEYR"])
    merged["sex"] = numeric(merged, ["RIAGENDR"])
    merged["race"] = numeric(merged, ["RIDRETH3", "RIDRETH1"])
    merged["pregnancy"] = numeric(merged, ["RIDEXPRG"])
    merged["education"] = numeric(merged, ["DMDEDUC2"])
    merged["pir"] = numeric(merged, ["INDFMPIR"])
    merged["weight"] = numeric(merged, ["WTMEC2YR"]) / len(CYCLES)
    merged["psu"] = numeric(merged, ["SDMVPSU"])
    merged["strata"] = numeric(merged, ["SDMVSTRA"])
    merged["psu_u"] = merged["cycle_index"] * 100 + merged["psu"]
    merged["strata_u"] = merged["cycle_index"] * 1000 + merged["strata"]
    merged["plp"] = numeric(merged, ["LBXPLP"])
    merged["four_pa"] = numeric(merged, ["LBX4PA"])
    merged["ferritin"] = numeric(merged, ["LBXFER", "LBDFERSI"])
    merged["hgb"] = numeric(merged, ["LBXHGB"])
    merged["mcv"] = numeric(merged, ["LBXMCVSI", "LBXMCV"])
    merged["mch"] = numeric(merged, ["LBXMCHSI"])
    merged["rdw"] = numeric(merged, ["LBXRDW"])
    merged["crp_mg_l"] = numeric(merged, ["LBXCRP"]) * 10.0
    merged["creatinine_mg_dl"] = numeric(merged, ["LBXSCR", "LBDSCR"])
    merged["albumin_g_dl"] = numeric(merged, ["LBXSAL"])
    merged["total_protein_g_dl"] = numeric(merged, ["LBXSTP"])
    merged["alt_u_l"] = numeric(merged, ["LBXSATSI"])
    merged["ast_u_l"] = numeric(merged, ["LBXAST"])
    merged["bmi"] = numeric(merged, ["BMXBMI"])
    merged["body_weight_kg"] = numeric(merged, ["BMXWT"])
    merged["height_cm"] = numeric(merged, ["BMXHT"])
    merged["smoked100"] = numeric(merged, ["SMQ020"])
    merged["current_smoke"] = numeric(merged, ["SMQ040"])
    merged["diabetes_q"] = numeric(merged, ["DIQ010"])
    merged["cancer_q"] = numeric(merged, ["MCQ220"])
    merged["liver_q"] = numeric(merged, ["MCQ160L"])
    merged["egfr"] = egfr_2021(merged["creatinine_mg_dl"], merged["age"], merged["sex"])
    return merged, audits


def main() -> None:
    frames: list[pd.DataFrame] = []
    audits: list[dict[str, object]] = []
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
    raw["known_liver_disease"] = np.where(raw["liver_q"].eq(1), 1.0, np.where(raw["liver_q"].eq(2), 0.0, np.nan))

    flow: list[dict[str, object]] = []
    def record(step: str, frame: pd.DataFrame) -> None:
        flow.append({"step": step, "n": int(len(frame))})

    record("all_selected_cycle_records", raw)
    cohort = raw.loc[raw["age"].between(20, 79, inclusive="both")].copy()
    record("age_20_79", cohort)
    cohort = cohort.loc[~cohort["pregnancy"].eq(1)].copy()
    record("exclude_known_pregnancy", cohort)
    cohort = cohort.loc[
        cohort["plp"].between(1, 1000, inclusive="both")
        & cohort["ferritin"].between(30, 500, inclusive="both")
        & cohort["hgb"].between(5, 20, inclusive="both")
        & cohort["mcv"].between(55, 125, inclusive="both")
    ].copy()
    record("valid_core_biomarkers_and_iron_replete", cohort)
    cohort = cohort.loc[
        cohort["crp_mg_l"].between(0, 10, inclusive="both")
        & cohort["egfr"].ge(60)
        & cohort["albumin_g_dl"].ge(3.5)
        & cohort["bmi"].between(15, 70, inclusive="both")
    ].copy()
    record("exclude_overt_inflammation_ckd_hypoalbuminemia", cohort)
    cohort = cohort.loc[~cohort["known_cancer"].eq(1) & ~cohort["known_liver_disease"].eq(1)].copy()
    record("exclude_known_cancer_or_liver_disease", cohort)
    cohort = cohort.loc[cohort["weight"].gt(0) & cohort["psu_u"].notna() & cohort["strata_u"].notna()].copy()
    record("valid_survey_design", cohort)

    cohort["exposure"] = np.select(
        [cohort["plp"].lt(20), cohort["plp"].ge(30)],
        ["deficient", "adequate"],
        default="marginal",
    )
    cohort = cohort.loc[cohort["exposure"].isin(["deficient", "adequate"])].copy()
    record("clean_frozen_exposure_comparison", cohort)

    cohort["anemia"] = np.where(cohort["sex"].eq(1), cohort["hgb"].lt(13.0), cohort["hgb"].lt(12.0)).astype(int)
    cohort["microcytosis"] = cohort["mcv"].lt(80.0).astype(int)
    cohort["nonmacrocytic_anemia"] = (cohort["anemia"].eq(1) & cohort["mcv"].lt(100.0)).astype(int)
    cohort["log_plp"] = np.log(cohort["plp"])
    cohort["log_crp"] = np.log(cohort["crp_mg_l"].clip(lower=0.05))
    cohort["log_ferritin"] = np.log(cohort["ferritin"])

    full_covariates = [
        "age", "sex", "race", "bmi", "log_crp", "log_ferritin", "albumin_g_dl",
        "smoking", "diabetes", "pir",
    ]
    cohort["complete_case"] = cohort[full_covariates].notna().all(axis=1)
    record("complete_case_full_adjustment", cohort.loc[cohort["complete_case"]])

    group_rows: list[dict[str, object]] = []
    for scope, frame in [
        ("overall", cohort),
        ("discovery", cohort.loc[cohort["period"].eq("discovery")]),
        ("validation", cohort.loc[cohort["period"].eq("validation")]),
    ]:
        for exposure, group in frame.groupby("exposure", observed=True):
            group_rows.append({
                "scope": scope,
                "exposure": exposure,
                "n": int(len(group)),
                "anemia_events": int(group["anemia"].sum()),
                "nonmacrocytic_anemia_events": int(group["nonmacrocytic_anemia"].sum()),
                "microcytosis_events": int(group["microcytosis"].sum()),
                "weighted_anemia_pct": 100.0 * weighted_mean(group["anemia"], group["weight"]),
                "weighted_mean_hgb": weighted_mean(group["hgb"], group["weight"]),
                "weighted_mean_mcv": weighted_mean(group["mcv"], group["weight"]),
                "weighted_mean_plp": weighted_mean(group["plp"], group["weight"]),
            })

    missing_rows: list[dict[str, object]] = []
    missing_vars = [
        "plp", "ferritin", "hgb", "mcv", "crp_mg_l", "creatinine_mg_dl", "egfr",
        "albumin_g_dl", "bmi", "smoking", "diabetes", "pir", "known_cancer", "known_liver_disease",
    ]
    for scope, frame in [("overall", cohort)] + [(f"cycle_{c}", cohort.loc[cohort["cycle"].eq(c)]) for c in CYCLES]:
        for variable in missing_vars:
            missing_rows.append({
                "scope": scope,
                "variable": variable,
                "n": int(len(frame)),
                "missing_n": int(frame[variable].isna().sum()),
                "missing_pct": float(frame[variable].isna().mean() * 100) if len(frame) else np.nan,
            })

    semantic = pd.DataFrame([
        {"concept": "biochemical_exposure", "fields": "LBXPLP", "meaning": "serum pyridoxal 5-prime-phosphate", "unit": "nmol/L", "quality_rule": "1-1000; deficient <20; adequate >=30; marginal excluded"},
        {"concept": "iron_replete_restriction", "fields": "LBXFER/LBDFERSI", "meaning": "serum ferritin", "unit": "ng/mL", "quality_rule": "30-500"},
        {"concept": "hematologic_outcomes", "fields": "LBXHGB/LBXMCVSI/LBXMCV", "meaning": "measured hemoglobin and mean corpuscular volume", "unit": "g/dL and fL", "quality_rule": "sex-specific anemia; measured continuous outcomes"},
        {"concept": "inflammation", "fields": "LBXCRP", "meaning": "C-reactive protein", "unit": "released mg/dL converted x10 to mg/L", "quality_rule": "0-10 mg/L; log transformed"},
        {"concept": "kidney_function", "fields": "LBXSCR/LBDSCR", "meaning": "serum creatinine and 2021 CKD-EPI eGFR", "unit": "mg/dL and mL/min/1.73m2", "quality_rule": "eGFR >=60"},
        {"concept": "survey_design", "fields": "WTMEC2YR/SDMVPSU/SDMVSTRA", "meaning": "MEC examination weight, PSU and stratum", "unit": "survey", "quality_rule": "weight divided by three pooled cycles; cycle-unique PSU/strata"},
    ])

    columns = [
        "SEQN", "cycle", "cycle_years", "period", "age", "sex", "race", "education", "pir",
        "weight", "psu_u", "strata_u", "plp", "four_pa", "ferritin", "hgb", "mcv", "mch", "rdw",
        "crp_mg_l", "creatinine_mg_dl", "egfr", "albumin_g_dl", "total_protein_g_dl", "alt_u_l", "ast_u_l",
        "bmi", "body_weight_kg", "height_cm", "smoking", "diabetes", "known_cancer", "known_liver_disease",
        "exposure", "anemia", "nonmacrocytic_anemia", "microcytosis", "log_plp", "log_crp", "log_ferritin",
        "complete_case",
    ]
    cohort.loc[:, columns].to_csv(OUT / "c_n09_analysis_core.csv", index=False)
    pd.DataFrame(audits).to_csv(OUT / "c_n09_source_audit.csv", index=False)
    pd.DataFrame(flow).to_csv(OUT / "c_n09_flow.csv", index=False)
    pd.DataFrame(group_rows).to_csv(OUT / "c_n09_group_counts.csv", index=False)
    pd.DataFrame(missing_rows).to_csv(OUT / "c_n09_missingness.csv", index=False)
    semantic.to_csv(OUT / "c_n09_semantic_audit.csv", index=False)

    key = cohort.loc[cohort["exposure"].eq("deficient")]
    val_key = key.loc[key["period"].eq("validation")]
    status = {
        "candidate_code": "C-N09",
        "actual_n": int(len(cohort)),
        "complete_case_n": int(cohort["complete_case"].sum()),
        "complete_case_retention": float(cohort["complete_case"].mean()),
        "anemia_events": int(cohort["anemia"].sum()),
        "key_group_n": int(len(key)),
        "key_group_anemia_events": int(key["anemia"].sum()),
        "validation_key_group_n": int(len(val_key)),
        "validation_key_group_anemia_events": int(val_key["anemia"].sum()),
        "cycles": sorted(cohort["cycle"].unique().tolist()),
    }
    (OUT / "c_n09_prep_status.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
    print(json.dumps(status, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
