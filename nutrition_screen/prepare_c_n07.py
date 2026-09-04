from __future__ import annotations

import io
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import requests

OUT = Path("nutrition_screen/output_c_n07")
OUT.mkdir(parents=True, exist_ok=True)

CYCLES = {
    "G": {"year": "2011", "years": "2011-2012", "index": 1, "period": "discovery"},
    "H": {"year": "2013", "years": "2013-2014", "index": 2, "period": "validation"},
}

COMPONENTS = {
    "demo": "DEMO",
    "day1": "DR1TOT",
    "day2": "DR2TOT",
    "grip": "MGX",
    "body": "BMX",
    "activity": "PAQ",
    "smoking": "SMQ",
    "diabetes": "DIQ",
    "hba1c": "GHB",
    "biochem": "BIOPRO",
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
        frame.columns = [str(column).upper() for column in frame.columns]
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
        audit = {
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
        return pd.DataFrame(columns=["SEQN"]), audit


def keep(frame: pd.DataFrame, columns: Iterable[str]) -> pd.DataFrame:
    selected = [column for column in columns if column in frame.columns]
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


def merge_cycle(cycle: str) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    frames: dict[str, pd.DataFrame] = {}
    audits: list[dict[str, object]] = []
    for component in COMPONENTS:
        frame, audit = load_xpt(cycle, component, required=True)
        frames[component] = frame
        audits.append(audit)

    pieces = [
        keep(
            frames["demo"],
            [
                "SEQN",
                "RIAGENDR",
                "RIDAGEYR",
                "RIDRETH1",
                "RIDRETH3",
                "WTMEC2YR",
                "SDMVPSU",
                "SDMVSTRA",
                "INDFMPIR",
                "DMDEDUC2",
            ],
        ),
        keep(
            frames["day1"],
            [
                "SEQN",
                "WTDR2D",
                "DR1DRSTZ",
                "DR1TPROT",
                "DR1TKCAL",
                "DR1_300",
                "DRQSDIET",
            ],
        ),
        keep(
            frames["day2"],
            ["SEQN", "DR2DRSTZ", "DR2TPROT", "DR2TKCAL", "DR2_300"],
        ),
        keep(
            frames["grip"],
            [
                "SEQN",
                "MGDEXSTS",
                "MGXH1T1",
                "MGXH1T2",
                "MGXH1T3",
                "MGXH2T1",
                "MGXH2T2",
                "MGXH2T3",
                "MGDCGSZ",
            ],
        ),
        keep(frames["body"], ["SEQN", "BMXWT", "BMXHT", "BMXBMI"]),
        keep(frames["activity"], ["SEQN", "PAQ650", "PAQ665", "PAD680"]),
        keep(frames["smoking"], ["SEQN", "SMQ020", "SMQ040"]),
        keep(frames["diabetes"], ["SEQN", "DIQ010"]),
        keep(frames["hba1c"], ["SEQN", "LBXGH"]),
        keep(frames["biochem"], ["SEQN", "LBXSCR", "LBXSAL", "LBXSTP"]),
    ]

    merged = pieces[0]
    for piece in pieces[1:]:
        merged = merged.merge(piece.drop_duplicates("SEQN"), on="SEQN", how="left")

    merged["cycle"] = cycle
    merged["cycle_years"] = CYCLES[cycle]["years"]
    merged["period"] = CYCLES[cycle]["period"]
    merged["cycle_index"] = CYCLES[cycle]["index"]
    merged["age"] = numeric(merged, ["RIDAGEYR"])
    merged["sex"] = numeric(merged, ["RIAGENDR"])
    merged["race"] = numeric(merged, ["RIDRETH3", "RIDRETH1"])
    merged["education"] = numeric(merged, ["DMDEDUC2"])
    merged["pir"] = numeric(merged, ["INDFMPIR"])
    merged["weight_2day"] = numeric(merged, ["WTDR2D"]) / len(CYCLES)
    merged["psu"] = numeric(merged, ["SDMVPSU"])
    merged["strata"] = numeric(merged, ["SDMVSTRA"])
    merged["psu_u"] = merged["cycle_index"] * 100 + merged["psu"]
    merged["strata_u"] = merged["cycle_index"] * 1000 + merged["strata"]
    merged["recall1_status"] = numeric(merged, ["DR1DRSTZ"])
    merged["recall2_status"] = numeric(merged, ["DR2DRSTZ"])
    merged["protein1_g"] = numeric(merged, ["DR1TPROT"])
    merged["protein2_g"] = numeric(merged, ["DR2TPROT"])
    merged["energy1_kcal"] = numeric(merged, ["DR1TKCAL"])
    merged["energy2_kcal"] = numeric(merged, ["DR2TKCAL"])
    merged["special_diet"] = numeric(merged, ["DRQSDIET"])
    merged["body_weight_kg"] = numeric(merged, ["BMXWT"])
    merged["height_cm"] = numeric(merged, ["BMXHT"])
    merged["bmi"] = numeric(merged, ["BMXBMI"])
    merged["grip_status"] = numeric(merged, ["MGDEXSTS"])
    grip_columns = [
        column
        for column in ["MGXH1T1", "MGXH1T2", "MGXH1T3", "MGXH2T1", "MGXH2T2", "MGXH2T3"]
        if column in merged.columns
    ]
    merged["grip_max_kg"] = merged[grip_columns].apply(pd.to_numeric, errors="coerce").max(axis=1)
    merged["grip_combined_kg"] = numeric(merged, ["MGDCGSZ"])
    merged["vigorous_recreation"] = numeric(merged, ["PAQ650"])
    merged["moderate_recreation"] = numeric(merged, ["PAQ665"])
    merged["sedentary_minutes"] = numeric(merged, ["PAD680"])
    merged["smoked100"] = numeric(merged, ["SMQ020"])
    merged["current_smoke"] = numeric(merged, ["SMQ040"])
    merged["diabetes_q"] = numeric(merged, ["DIQ010"])
    merged["hba1c"] = numeric(merged, ["LBXGH"])
    merged["creatinine_mg_dl"] = numeric(merged, ["LBXSCR"])
    merged["albumin_g_dl"] = numeric(merged, ["LBXSAL"])
    merged["serum_total_protein_g_dl"] = numeric(merged, ["LBXSTP"])
    merged["egfr"] = egfr_2021(merged["creatinine_mg_dl"], merged["age"], merged["sex"])
    return merged, audits


def weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    values = pd.to_numeric(values, errors="coerce")
    ok = values.notna() & weights.notna() & weights.gt(0)
    if not ok.any():
        return float("nan")
    return float(np.average(values.loc[ok], weights=weights.loc[ok]))


def main() -> None:
    all_frames: list[pd.DataFrame] = []
    all_audits: list[dict[str, object]] = []
    for cycle in CYCLES:
        frame, audits = merge_cycle(cycle)
        all_frames.append(frame)
        all_audits.extend(audits)
    raw = pd.concat(all_frames, ignore_index=True, sort=False)

    raw["active_recreation"] = np.where(
        raw["vigorous_recreation"].eq(1) | raw["moderate_recreation"].eq(1),
        1.0,
        np.where(
            raw["vigorous_recreation"].eq(2) & raw["moderate_recreation"].eq(2),
            0.0,
            np.nan,
        ),
    )
    raw["smoking"] = pd.Series(pd.NA, index=raw.index, dtype="string")
    raw.loc[raw["smoked100"].eq(2), "smoking"] = "never"
    raw.loc[raw["smoked100"].eq(1) & raw["current_smoke"].eq(3), "smoking"] = "former"
    raw.loc[raw["smoked100"].eq(1) & raw["current_smoke"].isin([1, 2]), "smoking"] = "current"
    raw["diabetes"] = np.where(
        raw["diabetes_q"].eq(1) | raw["hba1c"].ge(6.5),
        1.0,
        np.where(raw["diabetes_q"].isin([2, 3]) & raw["hba1c"].lt(6.5), 0.0, np.nan),
    )

    flow: list[dict[str, object]] = []

    def record(step: str, frame: pd.DataFrame) -> None:
        flow.append({"step": step, "n": int(len(frame))})

    record("all_selected_cycle_records", raw)
    cohort = raw.loc[raw["age"].ge(60)].copy()
    record("age_60_plus", cohort)
    cohort = cohort.loc[cohort["recall1_status"].eq(1) & cohort["recall2_status"].eq(1)].copy()
    record("two_reliable_dietary_recalls", cohort)
    cohort = cohort.loc[
        cohort["weight_2day"].gt(0)
        & cohort["psu_u"].notna()
        & cohort["strata_u"].notna()
    ].copy()
    record("valid_two_day_survey_design", cohort)
    cohort = cohort.loc[
        cohort["body_weight_kg"].between(30, 250, inclusive="both")
        & cohort["height_cm"].between(120, 210, inclusive="both")
        & cohort["bmi"].between(12, 70, inclusive="both")
    ].copy()
    record("valid_anthropometrics", cohort)
    cohort = cohort.loc[
        cohort["grip_status"].eq(1)
        & cohort["grip_max_kg"].between(2, 100, inclusive="both")
    ].copy()
    record("valid_bilateral_objective_grip", cohort)
    cohort = cohort.loc[
        cohort["protein1_g"].between(1, 400, inclusive="both")
        & cohort["protein2_g"].between(1, 400, inclusive="both")
        & cohort["energy1_kcal"].between(100, 10000, inclusive="both")
        & cohort["energy2_kcal"].between(100, 10000, inclusive="both")
    ].copy()
    record("valid_two_day_nutrient_intakes", cohort)

    cohort["protein1_gkg"] = cohort["protein1_g"] / cohort["body_weight_kg"]
    cohort["protein2_gkg"] = cohort["protein2_g"] / cohort["body_weight_kg"]
    cohort["mean_protein_gkg"] = (cohort["protein1_gkg"] + cohort["protein2_gkg"]) / 2.0
    cohort["absolute_protein_difference_gkg"] = (cohort["protein1_gkg"] - cohort["protein2_gkg"]).abs()
    cohort["protein_cv"] = cohort["absolute_protein_difference_gkg"] / cohort["mean_protein_gkg"].clip(lower=0.01)
    cohort["mean_energy_kcal"] = (cohort["energy1_kcal"] + cohort["energy2_kcal"]) / 2.0
    cohort["absolute_energy_difference_kcal"] = (cohort["energy1_kcal"] - cohort["energy2_kcal"]).abs()

    height_m = cohort["height_cm"] / 100.0
    ideal_weight_kg = 25.0 * height_m.pow(2)
    cohort["adjusted_weight_kg"] = np.where(
        cohort["bmi"].ge(30),
        ideal_weight_kg + 0.4 * (cohort["body_weight_kg"] - ideal_weight_kg),
        cohort["body_weight_kg"],
    )
    cohort["protein1_gkg_adjw"] = cohort["protein1_g"] / cohort["adjusted_weight_kg"]
    cohort["protein2_gkg_adjw"] = cohort["protein2_g"] / cohort["adjusted_weight_kg"]
    cohort["mean_protein_gkg_adjw"] = (cohort["protein1_gkg_adjw"] + cohort["protein2_gkg_adjw"]) / 2.0

    cohort = cohort.loc[cohort["mean_protein_gkg"].ge(0.8) & cohort["mean_protein_gkg"].le(3.0)].copy()
    record("mean_two_day_protein_at_least_rda", cohort)
    cohort["pattern"] = np.where(
        cohort["protein1_gkg"].ge(0.8) & cohort["protein2_gkg"].ge(0.8),
        "consistent_adequate",
        "unstable_adequacy",
    )
    cohort["pattern_adjw"] = np.where(
        cohort["mean_protein_gkg_adjw"].ge(0.8)
        & cohort["protein1_gkg_adjw"].ge(0.8)
        & cohort["protein2_gkg_adjw"].ge(0.8),
        "consistent_adequate",
        np.where(
            cohort["mean_protein_gkg_adjw"].ge(0.8),
            "unstable_adequacy",
            "below_mean_requirement",
        ),
    )
    record("classified_primary_exposure", cohort)

    cohort["weak_fnih"] = np.where(
        cohort["sex"].eq(1),
        (cohort["grip_max_kg"] < 26).astype(int),
        (cohort["grip_max_kg"] < 16).astype(int),
    )
    cohort["weak_ewgsop2"] = np.where(
        cohort["sex"].eq(1),
        (cohort["grip_max_kg"] < 27).astype(int),
        (cohort["grip_max_kg"] < 16).astype(int),
    )
    cohort["grip_bmi_ratio"] = cohort["grip_max_kg"] / cohort["bmi"]
    cohort["weak_fnih_bmi"] = np.where(
        cohort["sex"].eq(1),
        (cohort["grip_bmi_ratio"] < 1.0).astype(int),
        (cohort["grip_bmi_ratio"] < 0.56).astype(int),
    )

    full_covariates = [
        "age",
        "sex",
        "race",
        "bmi",
        "height_cm",
        "mean_protein_gkg",
        "mean_energy_kcal",
        "active_recreation",
        "sedentary_minutes",
        "smoking",
        "diabetes",
        "egfr",
        "pir",
    ]
    cohort["complete_case"] = cohort[full_covariates].notna().all(axis=1)
    record("complete_case_full_adjustment", cohort.loc[cohort["complete_case"]])

    group_rows: list[dict[str, object]] = []
    scopes = [
        ("overall", cohort),
        ("discovery", cohort.loc[cohort["period"].eq("discovery")]),
        ("validation", cohort.loc[cohort["period"].eq("validation")]),
    ]
    for scope, frame in scopes:
        for pattern, group in frame.groupby("pattern", observed=True):
            group_rows.append(
                {
                    "scope": scope,
                    "pattern": pattern,
                    "n": int(len(group)),
                    "weak_fnih_events": int(group["weak_fnih"].sum()),
                    "weak_ewgsop2_events": int(group["weak_ewgsop2"].sum()),
                    "weighted_weak_fnih_pct": 100.0 * weighted_mean(group["weak_fnih"], group["weight_2day"]),
                    "mean_protein_gkg": float(group["mean_protein_gkg"].mean()),
                    "mean_absolute_difference_gkg": float(group["absolute_protein_difference_gkg"].mean()),
                    "mean_grip_kg": float(group["grip_max_kg"].mean()),
                }
            )

    missing_rows: list[dict[str, object]] = []
    variables = [
        "protein1_g",
        "protein2_g",
        "grip_max_kg",
        "body_weight_kg",
        "height_cm",
        "bmi",
        "mean_energy_kcal",
        "active_recreation",
        "sedentary_minutes",
        "smoking",
        "diabetes",
        "egfr",
        "pir",
        "albumin_g_dl",
    ]
    missing_scopes = [("overall", cohort)] + [
        (f"cycle_{cycle}", cohort.loc[cohort["cycle"].eq(cycle)]) for cycle in CYCLES
    ]
    for scope, frame in missing_scopes:
        for variable in variables:
            missing_rows.append(
                {
                    "scope": scope,
                    "variable": variable,
                    "n": int(len(frame)),
                    "missing_n": int(frame[variable].isna().sum()),
                    "missing_pct": float(frame[variable].isna().mean() * 100) if len(frame) else np.nan,
                }
            )

    semantic = pd.DataFrame(
        [
            {
                "concept": "two_day_protein_intake",
                "fields": "DR1TPROT, DR2TPROT",
                "semantics": "grams of protein from foods and beverages on each 24-hour recall day",
                "unit": "g/day and g/kg/day using measured body weight",
                "quality_rule": "both recalls reliable; 1-400 g/day",
            },
            {
                "concept": "two_day_diet_weight",
                "fields": "WTDR2D",
                "semantics": "NHANES two-day dietary sample weight",
                "unit": "survey weight",
                "quality_rule": "positive weight divided by two pooled cycles",
            },
            {
                "concept": "objective_grip_strength",
                "fields": "MGXH1T1-3 and MGXH2T1-3; MGDEXSTS",
                "semantics": "maximum measured grip across six trials among participants completing both hands",
                "unit": "kg",
                "quality_rule": "status=1; 2-100 kg",
            },
            {
                "concept": "primary_weakness",
                "fields": "derived from grip_max_kg",
                "semantics": "FNIH clinically relevant weakness",
                "unit": "binary",
                "quality_rule": "men <26 kg; women <16 kg",
            },
            {
                "concept": "renal_function",
                "fields": "LBXSCR, age, sex",
                "semantics": "2021 CKD-EPI race-free eGFR",
                "unit": "mL/min/1.73m2",
                "quality_rule": "adjusted in main model; eGFR >=60 sensitivity",
            },
            {
                "concept": "exposure_pattern",
                "fields": "two daily protein g/kg values",
                "semantics": "among mean intake >=0.8 g/kg/day, both days adequate versus one day below 0.8",
                "unit": "binary",
                "quality_rule": "thresholds frozen before outcome modeling",
            },
        ]
    )

    output_columns = [
        "SEQN",
        "cycle",
        "cycle_years",
        "period",
        "age",
        "sex",
        "race",
        "education",
        "pir",
        "weight_2day",
        "psu_u",
        "strata_u",
        "protein1_g",
        "protein2_g",
        "protein1_gkg",
        "protein2_gkg",
        "mean_protein_gkg",
        "absolute_protein_difference_gkg",
        "protein_cv",
        "energy1_kcal",
        "energy2_kcal",
        "mean_energy_kcal",
        "absolute_energy_difference_kcal",
        "body_weight_kg",
        "adjusted_weight_kg",
        "height_cm",
        "bmi",
        "protein1_gkg_adjw",
        "protein2_gkg_adjw",
        "mean_protein_gkg_adjw",
        "pattern",
        "pattern_adjw",
        "grip_max_kg",
        "grip_combined_kg",
        "weak_fnih",
        "weak_ewgsop2",
        "grip_bmi_ratio",
        "weak_fnih_bmi",
        "active_recreation",
        "sedentary_minutes",
        "smoking",
        "diabetes",
        "hba1c",
        "creatinine_mg_dl",
        "egfr",
        "albumin_g_dl",
        "serum_total_protein_g_dl",
        "special_diet",
        "complete_case",
    ]
    cohort.loc[:, output_columns].to_csv(OUT / "c_n07_analysis_core.csv", index=False)
    pd.DataFrame(all_audits).to_csv(OUT / "c_n07_source_audit.csv", index=False)
    pd.DataFrame(flow).to_csv(OUT / "c_n07_flow.csv", index=False)
    pd.DataFrame(group_rows).to_csv(OUT / "c_n07_group_counts.csv", index=False)
    pd.DataFrame(missing_rows).to_csv(OUT / "c_n07_missingness.csv", index=False)
    semantic.to_csv(OUT / "c_n07_semantic_audit.csv", index=False)

    key = cohort.loc[cohort["pattern"].eq("unstable_adequacy")]
    validation_key = key.loc[key["period"].eq("validation")]
    status = {
        "candidate_code": "C-N07",
        "actual_n": int(len(cohort)),
        "complete_case_n": int(cohort["complete_case"].sum()),
        "complete_case_retention": float(cohort["complete_case"].mean()),
        "primary_events": int(cohort["weak_fnih"].sum()),
        "key_group_n": int(len(key)),
        "key_group_events": int(key["weak_fnih"].sum()),
        "validation_key_group_n": int(len(validation_key)),
        "validation_key_group_events": int(validation_key["weak_fnih"].sum()),
        "cycles": sorted(cohort["cycle"].unique().tolist()),
    }
    (OUT / "c_n07_prep_status.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
    print(json.dumps(status, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
