from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from prepare_c_n13 import CYCLES, merge_cycle, weighted_mean

OUT = Path("nutrition_stage/output_c_n14")
OUT.mkdir(parents=True, exist_ok=True)


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
    cohort = raw.loc[raw["age"].between(50, 79, inclusive="both")].copy()
    record("age_50_79", cohort)
    cohort = cohort.loc[cohort["pregnancy"].ne(1) | cohort["pregnancy"].isna()].copy()
    record("exclude_known_pregnancy", cohort)
    cohort = cohort.loc[cohort["thiazide"] & ~cohort["loop"] & ~cohort["k_sparing"]].copy()
    record("current_thiazide_only_without_potassium_sparing", cohort)
    cohort = cohort.loc[cohort["heart_failure"].ne(1) | cohort["heart_failure"].isna()].copy()
    record("exclude_known_heart_failure", cohort)
    cohort = cohort.loc[cohort["diet1_valid"] & cohort["diet2_valid"]].copy()
    record("two_valid_dietary_recalls", cohort)
    cohort = cohort.loc[
        cohort["mean_kcal"].between(600, 6000, inclusive="both")
        & cohort["mean_sodium"].between(300, 10000, inclusive="both")
        & cohort["mean_potassium"].between(300, 7000, inclusive="both")
        & cohort["mean_protein"].between(10, 300, inclusive="both")
        & cohort["serum_sodium"].between(120, 155, inclusive="both")
        & cohort["egfr"].ge(30)
        & cohort["bmi"].between(15, 70, inclusive="both")
        & cohort["weight"].gt(0)
        & cohort["psu_u"].notna()
        & cohort["strata_u"].notna()
    ].copy()
    record("valid_nutrition_electrolyte_and_survey_data", cohort)

    cohort["exposure"] = np.select(
        [cohort["mean_sodium"].lt(2000), cohort["mean_sodium"].ge(2500)],
        ["low", "higher"], default="intermediate"
    )
    cohort = cohort.loc[cohort["exposure"].isin(["low", "higher"])].copy()
    record("frozen_low_vs_higher_sodium_comparison", cohort)

    cohort["hyponatremia_135"] = cohort["serum_sodium"].lt(136).astype(int)
    cohort["strict_hyponatremia_134"] = cohort["serum_sodium"].lt(135).astype(int)
    cohort["low_normal_na_138"] = cohort["serum_sodium"].lt(138).astype(int)
    cohort["log_sodium_intake"] = np.log(cohort["mean_sodium"])
    cohort["sodium_density"] = cohort["mean_sodium"] / cohort["mean_kcal"] * 1000.0
    cohort["chronic_thiazide"] = cohort["wasting_days"].ge(30)

    covars = ["age", "sex", "race", "bmi", "egfr", "diabetes", "hypertension", "cvd", "mean_kcal", "mean_protein", "mean_potassium", "mean_magnesium", "raas", "ppi", "pir"]
    cohort["complete_case"] = cohort[covars].notna().all(axis=1)
    record("complete_case_full_adjustment", cohort.loc[cohort["complete_case"]])

    group_rows = []
    for scope, frame in [("overall", cohort), ("discovery", cohort.loc[cohort["period"].eq("discovery")]), ("validation", cohort.loc[cohort["period"].eq("validation")])]:
        for exposure, group in frame.groupby("exposure", observed=True):
            group_rows.append({
                "scope": scope, "exposure": exposure, "n": int(len(group)),
                "hyponatremia_135_events": int(group["hyponatremia_135"].sum()),
                "strict_hyponatremia_134_events": int(group["strict_hyponatremia_134"].sum()),
                "low_normal_na_138_events": int(group["low_normal_na_138"].sum()),
                "weighted_hyponatremia_135_pct": 100 * weighted_mean(group["hyponatremia_135"], group["weight"]),
                "weighted_mean_serum_sodium": weighted_mean(group["serum_sodium"], group["weight"]),
                "weighted_mean_dietary_sodium": weighted_mean(group["mean_sodium"], group["weight"]),
            })

    missing_rows = []
    for scope, frame in [("overall", cohort)] + [(f"cycle_{c}", cohort.loc[cohort["cycle"].eq(c)]) for c in CYCLES]:
        for variable in ["mean_sodium", "mean_potassium", "mean_protein", "mean_magnesium", "mean_kcal", "serum_sodium", "egfr", "bmi", "diabetes", "hypertension", "cvd", "pir", "wasting_days"]:
            missing_rows.append({"scope": scope, "variable": variable, "n": int(len(frame)), "missing_n": int(frame[variable].isna().sum()), "missing_pct": float(frame[variable].isna().mean() * 100) if len(frame) else np.nan})

    semantic = pd.DataFrame([
        {"concept": "dietary_sodium", "fields": "DR1TSODI and DR2TSODI", "meaning": "mean sodium from two complete 24-hour recalls", "unit": "mg/day", "rule": "low <2000; higher >=2500; intermediate excluded"},
        {"concept": "drug_exposure", "fields": "RXDDRUG in RXQ_RX", "meaning": "generic-name thiazide use in past 30 days", "unit": "person-level class", "rule": "thiazide only; loop and potassium-sparing agents excluded"},
        {"concept": "outcome", "fields": "LBXSNASI", "meaning": "measured serum sodium", "unit": "mmol/L", "rule": "primary <=135; strict <135; continuous secondary"},
        {"concept": "survey_design", "fields": "WTDR2D, SDMVPSU, SDMVSTRA", "meaning": "day-2 dietary weight and NHANES design", "unit": "survey", "rule": "pooled weight divided by six; cycle-unique PSU/strata"},
    ])

    columns = ["SEQN", "cycle", "period", "age", "sex", "race", "education", "pir", "weight", "psu_u", "strata_u", "mean_kcal", "mean_protein", "mean_potassium", "mean_sodium", "mean_magnesium", "day1_sodium", "day2_sodium", "sodium_density", "serum_sodium", "serum_potassium", "serum_chloride", "serum_bicarbonate", "creatinine", "egfr", "albumin", "bmi", "body_weight", "height", "diabetes", "hypertension", "cvd", "heart_failure", "wasting_days", "chronic_thiazide", "raas", "ppi", "exposure", "hyponatremia_135", "strict_hyponatremia_134", "low_normal_na_138", "log_sodium_intake", "complete_case"]
    cohort.loc[:, columns].to_csv(OUT / "c_n14_analysis_core.csv", index=False)
    pd.DataFrame(audits).to_csv(OUT / "c_n14_source_audit.csv", index=False)
    pd.DataFrame(flow).to_csv(OUT / "c_n14_flow.csv", index=False)
    pd.DataFrame(group_rows).to_csv(OUT / "c_n14_group_counts.csv", index=False)
    pd.DataFrame(missing_rows).to_csv(OUT / "c_n14_missingness.csv", index=False)
    semantic.to_csv(OUT / "c_n14_semantic_audit.csv", index=False)

    key = cohort.loc[cohort["exposure"].eq("low")]
    val_key = key.loc[key["period"].eq("validation")]
    prep = {
        "candidate_code": "C-N14", "actual_n": int(len(cohort)),
        "complete_case_n": int(cohort["complete_case"].sum()), "complete_case_retention": float(cohort["complete_case"].mean()),
        "primary_events": int(cohort["hyponatremia_135"].sum()), "strict_events": int(cohort["strict_hyponatremia_134"].sum()),
        "key_group_n": int(len(key)), "key_group_events": int(key["hyponatremia_135"].sum()),
        "validation_key_group_n": int(len(val_key)), "validation_key_group_events": int(val_key["hyponatremia_135"].sum()),
        "cycles": sorted(cohort["cycle"].unique().tolist()),
    }
    (OUT / "c_n14_prep_status.json").write_text(json.dumps(prep, indent=2), encoding="utf-8")
    print(json.dumps(prep, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
