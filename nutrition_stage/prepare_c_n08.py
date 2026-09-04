from __future__ import annotations

import io
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import requests

OUT = Path("nutrition_stage/output_c_n08")
OUT.mkdir(parents=True, exist_ok=True)

CYCLES = {
    "E": {"year": "2007", "years": "2007-2008", "mort": "2007_2008", "index": 1, "period": "discovery"},
    "F": {"year": "2009", "years": "2009-2010", "mort": "2009_2010", "index": 2, "period": "discovery"},
    "G": {"year": "2011", "years": "2011-2012", "mort": "2011_2012", "index": 3, "period": "discovery"},
    "H": {"year": "2013", "years": "2013-2014", "mort": "2013_2014", "index": 4, "period": "validation"},
    "I": {"year": "2015", "years": "2015-2016", "mort": "2015_2016", "index": 5, "period": "validation"},
    "J": {"year": "2017", "years": "2017-2018", "mort": "2017_2018", "index": 6, "period": "validation"},
}

COMPONENTS = {
    "demo": "DEMO",
    "day1": "DR1TOT",
    "day2": "DR2TOT",
    "body": "BMX",
    "activity": "PAQ",
    "smoking": "SMQ",
    "diabetes": "DIQ",
    "hba1c": "GHB",
    "biochem": "BIOPRO",
    "medical": "MCQ",
}


def xpt_url(cycle: str, stem: str) -> str:
    year = CYCLES[cycle]["year"]
    return f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{year}/DataFiles/{stem}_{cycle}.xpt"


def mortality_url(cycle: str) -> str:
    mort = CYCLES[cycle]["mort"]
    return (
        "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/"
        f"NHANES_{mort}_MORT_2019_PUBLIC.dat"
    )


def load_xpt(cycle: str, component: str) -> tuple[pd.DataFrame, dict[str, object]]:
    stem = COMPONENTS[component]
    url = xpt_url(cycle, stem)
    response = requests.get(url, timeout=180)
    if response.status_code != 200:
        raise RuntimeError(f"Unable to load {component} {cycle}: HTTP {response.status_code}")
    frame = pd.read_sas(io.BytesIO(response.content), format="xport", encoding="latin1")
    frame.columns = [str(column).upper() for column in frame.columns]
    audit = {
        "cycle": cycle,
        "component": component,
        "file": f"{stem}_{cycle}.xpt",
        "url": url,
        "rows": int(frame.shape[0]),
        "columns": int(frame.shape[1]),
        "column_names": "|".join(frame.columns),
    }
    return frame, audit


def load_mortality(cycle: str) -> tuple[pd.DataFrame, dict[str, object]]:
    url = mortality_url(cycle)
    response = requests.get(url, timeout=180)
    if response.status_code != 200:
        raise RuntimeError(f"Unable to load mortality {cycle}: HTTP {response.status_code}")
    lines = response.text.splitlines()
    records: list[dict[str, object]] = []
    for line in lines:
        if len(line) < 48:
            continue
        try:
            records.append(
                {
                    "SEQN": int(line[0:14].strip()),
                    "ELIGSTAT": int(line[14:15].strip()) if line[14:15].strip() else np.nan,
                    "MORTSTAT": int(line[15:16].strip()) if line[15:16].strip() else np.nan,
                    "UCOD_LEADING": line[16:19].strip(),
                    "DIABETES_MORT": int(line[19:20].strip()) if line[19:20].strip() else np.nan,
                    "HYPERTEN_MORT": int(line[20:21].strip()) if line[20:21].strip() else np.nan,
                    "PERMTH_INT": int(line[42:45].strip()) if line[42:45].strip() else np.nan,
                    "PERMTH_EXM": int(line[45:48].strip()) if line[45:48].strip() else np.nan,
                }
            )
        except ValueError:
            continue
    frame = pd.DataFrame(records)
    audit = {
        "cycle": cycle,
        "component": "mortality",
        "file": url.rsplit("/", 1)[-1],
        "url": url,
        "rows": int(frame.shape[0]),
        "columns": int(frame.shape[1]),
        "column_names": "|".join(frame.columns),
    }
    return frame, audit


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
        frame, audit = load_xpt(cycle, component)
        frames[component] = frame
        audits.append(audit)
    mortality, mortality_audit = load_mortality(cycle)
    audits.append(mortality_audit)

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
        keep(frames["day1"], ["SEQN", "WTDR2D", "DR1DRSTZ", "DR1TPROT", "DR1TKCAL", "DRQSDIET"]),
        keep(frames["day2"], ["SEQN", "DR2DRSTZ", "DR2TPROT", "DR2TKCAL"]),
        keep(frames["body"], ["SEQN", "BMXWT", "BMXHT", "BMXBMI"]),
        keep(frames["activity"], ["SEQN", "PAQ650", "PAQ665", "PAD680"]),
        keep(frames["smoking"], ["SEQN", "SMQ020", "SMQ040"]),
        keep(frames["diabetes"], ["SEQN", "DIQ010"]),
        keep(frames["hba1c"], ["SEQN", "LBXGH"]),
        keep(frames["biochem"], ["SEQN", "LBXSCR", "LBXSAL"]),
        keep(
            frames["medical"],
            ["SEQN", "MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E", "MCQ160F", "MCQ220"],
        ),
        mortality,
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
    merged["vigorous_recreation"] = numeric(merged, ["PAQ650"])
    merged["moderate_recreation"] = numeric(merged, ["PAQ665"])
    merged["sedentary_minutes"] = numeric(merged, ["PAD680"])
    merged["smoked100"] = numeric(merged, ["SMQ020"])
    merged["current_smoke"] = numeric(merged, ["SMQ040"])
    merged["diabetes_q"] = numeric(merged, ["DIQ010"])
    merged["hba1c"] = numeric(merged, ["LBXGH"])
    merged["creatinine_mg_dl"] = numeric(merged, ["LBXSCR"])
    merged["albumin_g_dl"] = numeric(merged, ["LBXSAL"])
    merged["egfr"] = egfr_2021(merged["creatinine_mg_dl"], merged["age"], merged["sex"])
    cvd_columns = [column for column in ["MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E"] if column in merged.columns]
    merged["cvd"] = np.where(merged[cvd_columns].eq(1).any(axis=1), 1.0, np.where(merged[cvd_columns].isin([2]).all(axis=1), 0.0, np.nan))
    merged["cancer"] = np.where(numeric(merged, ["MCQ220"]).eq(1), 1.0, np.where(numeric(merged, ["MCQ220"]).eq(2), 0.0, np.nan))
    merged["mortality_eligible"] = numeric(merged, ["ELIGSTAT"])
    merged["death"] = numeric(merged, ["MORTSTAT"])
    merged["follow_up_months"] = numeric(merged, ["PERMTH_EXM"])
    return merged, audits


def weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    values = pd.to_numeric(values, errors="coerce")
    valid = values.notna() & weights.notna() & weights.gt(0)
    if not valid.any():
        return float("nan")
    return float(np.average(values.loc[valid], weights=weights.loc[valid]))


def main() -> None:
    frames: list[pd.DataFrame] = []
    audits: list[dict[str, object]] = []
    for cycle in CYCLES:
        frame, audit = merge_cycle(cycle)
        frames.append(frame)
        audits.extend(audit)
    raw = pd.concat(frames, ignore_index=True, sort=False)

    raw["active_recreation"] = np.where(
        raw["vigorous_recreation"].eq(1) | raw["moderate_recreation"].eq(1),
        1.0,
        np.where(raw["vigorous_recreation"].eq(2) & raw["moderate_recreation"].eq(2), 0.0, np.nan),
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
    cohort = raw.loc[raw["age"].ge(50)].copy()
    record("age_50_plus", cohort)
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
        cohort["protein1_g"].between(1, 400, inclusive="both")
        & cohort["protein2_g"].between(1, 400, inclusive="both")
        & cohort["energy1_kcal"].between(100, 10000, inclusive="both")
        & cohort["energy2_kcal"].between(100, 10000, inclusive="both")
    ].copy()
    record("valid_two_day_nutrient_intakes", cohort)
    cohort = cohort.loc[
        cohort["mortality_eligible"].eq(1)
        & cohort["death"].isin([0, 1])
        & cohort["follow_up_months"].gt(0)
    ].copy()
    record("valid_linked_mortality_follow_up", cohort)

    cohort["protein1_gkg"] = cohort["protein1_g"] / cohort["body_weight_kg"]
    cohort["protein2_gkg"] = cohort["protein2_g"] / cohort["body_weight_kg"]
    cohort["mean_protein_gkg"] = (cohort["protein1_gkg"] + cohort["protein2_gkg"]) / 2.0
    cohort["minimum_protein_gkg"] = cohort[["protein1_gkg", "protein2_gkg"]].min(axis=1)
    cohort["absolute_protein_difference_gkg"] = (cohort["protein1_gkg"] - cohort["protein2_gkg"]).abs()
    cohort["mean_energy_kcal"] = (cohort["energy1_kcal"] + cohort["energy2_kcal"]) / 2.0
    cohort["absolute_energy_difference_kcal"] = (cohort["energy1_kcal"] - cohort["energy2_kcal"]).abs()
    cohort["follow_up_years"] = cohort["follow_up_months"] / 12.0

    height_m = cohort["height_cm"] / 100.0
    ideal_weight = 25.0 * height_m.pow(2)
    cohort["adjusted_weight_kg"] = np.where(
        cohort["bmi"].ge(30),
        ideal_weight + 0.4 * (cohort["body_weight_kg"] - ideal_weight),
        cohort["body_weight_kg"],
    )
    cohort["protein1_gkg_adjw"] = cohort["protein1_g"] / cohort["adjusted_weight_kg"]
    cohort["protein2_gkg_adjw"] = cohort["protein2_g"] / cohort["adjusted_weight_kg"]
    cohort["mean_protein_gkg_adjw"] = (cohort["protein1_gkg_adjw"] + cohort["protein2_gkg_adjw"]) / 2.0

    cohort = cohort.loc[cohort["mean_protein_gkg"].between(0.8, 3.0, inclusive="both")].copy()
    record("mean_two_day_protein_at_least_rda", cohort)
    cohort["pattern"] = np.select(
        [
            cohort["protein1_gkg"].ge(0.8) & cohort["protein2_gkg"].ge(0.8),
            cohort["minimum_protein_gkg"].lt(0.66),
        ],
        ["consistent_rda", "episodic_sub_ear"],
        default="intermediate",
    )
    cohort = cohort.loc[cohort["pattern"].isin(["consistent_rda", "episodic_sub_ear"])].copy()
    record("clean_frozen_exposure_comparison", cohort)

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
        "cvd",
        "cancer",
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
                    "deaths": int(group["death"].sum()),
                    "weighted_death_pct": 100.0 * weighted_mean(group["death"], group["weight_2day"]),
                    "mean_follow_up_years": float(group["follow_up_years"].mean()),
                    "mean_protein_gkg": float(group["mean_protein_gkg"].mean()),
                    "mean_minimum_day_gkg": float(group["minimum_protein_gkg"].mean()),
                }
            )

    missing_rows: list[dict[str, object]] = []
    variables = [
        "protein1_g",
        "protein2_g",
        "body_weight_kg",
        "height_cm",
        "bmi",
        "mean_energy_kcal",
        "active_recreation",
        "sedentary_minutes",
        "smoking",
        "diabetes",
        "egfr",
        "cvd",
        "cancer",
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
                "meaning": "two reliable 24-hour food-and-beverage protein intakes",
                "unit": "g/kg measured body weight/day",
                "quality_rule": "both recalls reliable; 1-400 g/day",
            },
            {
                "concept": "episodic_sub_ear",
                "fields": "minimum of two protein g/kg days",
                "meaning": "mean intake >=0.8 g/kg/day but one day <0.66 g/kg/day",
                "unit": "binary",
                "quality_rule": "clean comparison excludes intermediate 0.66-0.79 day",
            },
            {
                "concept": "mortality",
                "fields": "MORTSTAT and PERMTH_EXM",
                "meaning": "public-use NDI-linked all-cause mortality through 2019",
                "unit": "event and months from MEC examination",
                "quality_rule": "ELIGSTAT=1; positive follow-up",
            },
            {
                "concept": "survey_design",
                "fields": "WTDR2D, SDMVPSU, SDMVSTRA",
                "meaning": "two-day dietary weights and design variables",
                "unit": "survey",
                "quality_rule": "weights divided by six pooled cycles; cycle-unique PSU/strata",
            },
        ]
    )

    columns = [
        "SEQN", "cycle", "cycle_years", "period", "age", "sex", "race", "education", "pir",
        "weight_2day", "psu_u", "strata_u", "protein1_g", "protein2_g", "protein1_gkg",
        "protein2_gkg", "mean_protein_gkg", "minimum_protein_gkg", "absolute_protein_difference_gkg",
        "energy1_kcal", "energy2_kcal", "mean_energy_kcal", "absolute_energy_difference_kcal",
        "body_weight_kg", "adjusted_weight_kg", "height_cm", "bmi", "protein1_gkg_adjw",
        "protein2_gkg_adjw", "mean_protein_gkg_adjw", "pattern", "active_recreation",
        "sedentary_minutes", "smoking", "diabetes", "hba1c", "creatinine_mg_dl", "egfr",
        "cvd", "cancer", "albumin_g_dl", "special_diet", "death", "follow_up_months",
        "follow_up_years", "complete_case",
    ]
    cohort.loc[:, columns].to_csv(OUT / "c_n08_analysis_core.csv", index=False)
    pd.DataFrame(audits).to_csv(OUT / "c_n08_source_audit.csv", index=False)
    pd.DataFrame(flow).to_csv(OUT / "c_n08_flow.csv", index=False)
    pd.DataFrame(group_rows).to_csv(OUT / "c_n08_group_counts.csv", index=False)
    pd.DataFrame(missing_rows).to_csv(OUT / "c_n08_missingness.csv", index=False)
    semantic.to_csv(OUT / "c_n08_semantic_audit.csv", index=False)

    key = cohort.loc[cohort["pattern"].eq("episodic_sub_ear")]
    validation_key = key.loc[key["period"].eq("validation")]
    status = {
        "candidate_code": "C-N08",
        "actual_n": int(len(cohort)),
        "complete_case_n": int(cohort["complete_case"].sum()),
        "complete_case_retention": float(cohort["complete_case"].mean()),
        "deaths": int(cohort["death"].sum()),
        "key_group_n": int(len(key)),
        "key_group_deaths": int(key["death"].sum()),
        "validation_key_group_n": int(len(validation_key)),
        "validation_key_group_deaths": int(validation_key["death"].sum()),
        "cycles": sorted(cohort["cycle"].unique().tolist()),
    }
    (OUT / "c_n08_prep_status.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
    print(json.dumps(status, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
