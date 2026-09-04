from __future__ import annotations

import io
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import requests

OUT = Path("nutrition_screen/output_c_n06")
OUT.mkdir(parents=True, exist_ok=True)

CYCLES = {
    "A": {"year": "1999", "years": "1999-2000", "index": 1},
    "B": {"year": "2001", "years": "2001-2002", "index": 2},
    "C": {"year": "2003", "years": "2003-2004", "index": 3},
}

FILES = {
    "A": {
        "demo": ["DEMO"],
        "pn": ["LEXPN"],
        "nutrition": ["LAB06"],
        "hba1c": ["LAB10"],
        "biochem": ["LAB18"],
        "bmx": ["BMX"],
        "diq": ["DIQ"],
        "smq": ["SMQ"],
    },
    "B": {
        "demo": ["DEMO_B"],
        "pn": ["LEXPN_B"],
        "nutrition": ["L06_B"],
        "hba1c": ["L10_B"],
        "biochem": ["L40_B"],
        "bmx": ["BMX_B"],
        "diq": ["DIQ_B"],
        "smq": ["SMQ_B"],
    },
    "C": {
        "demo": ["DEMO_C"],
        "pn": ["LEXPN_C"],
        "b12": ["L06NB_C"],
        "mma": ["L06MH_C"],
        "hba1c": ["L10_C"],
        "biochem": ["L40_C"],
        "bmx": ["BMX_C"],
        "diq": ["DIQ_C"],
        "smq": ["SMQ_C"],
    },
}


def xpt_url(cycle: str, stem: str) -> str:
    year = CYCLES[cycle]["year"]
    return f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{year}/DataFiles/{stem}.xpt"


def read_xpt(cycle: str, stems: list[str], required: bool = True) -> tuple[pd.DataFrame, str, str | None]:
    errors: list[str] = []
    for stem in stems:
        url = xpt_url(cycle, stem)
        try:
            response = requests.get(url, timeout=120)
            if response.status_code != 200:
                errors.append(f"{stem}:HTTP{response.status_code}")
                continue
            frame = pd.read_sas(io.BytesIO(response.content), format="xport", encoding="latin1")
            frame.columns = [str(c).upper() for c in frame.columns]
            return frame, stem, None
        except Exception as exc:  # noqa: BLE001
            errors.append(f"{stem}:{type(exc).__name__}:{exc}")
    if required:
        raise RuntimeError(f"Unable to load required component for cycle {cycle}: {' | '.join(errors)}")
    return pd.DataFrame(columns=["SEQN"]), "", " | ".join(errors)


def keep(frame: pd.DataFrame, columns: Iterable[str]) -> pd.DataFrame:
    wanted = [c for c in columns if c in frame.columns]
    if "SEQN" in frame.columns and "SEQN" not in wanted:
        wanted.insert(0, "SEQN")
    return frame.loc[:, wanted].copy()


def first_numeric(frame: pd.DataFrame, names: Iterable[str]) -> pd.Series:
    result = pd.Series(np.nan, index=frame.index, dtype="float64")
    for name in names:
        if name in frame.columns:
            result = result.fillna(pd.to_numeric(frame[name], errors="coerce"))
    return result


def merge_cycle(cycle: str) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    audit: list[dict[str, object]] = []

    def load(key: str, required: bool = True) -> pd.DataFrame:
        frame, selected, error = read_xpt(cycle, FILES[cycle][key], required=required)
        audit.append(
            {
                "cycle": cycle,
                "component": key,
                "selected_file": selected,
                "required": required,
                "rows": int(frame.shape[0]),
                "columns": int(frame.shape[1]),
                "error": error,
                "column_names": "|".join(frame.columns),
            }
        )
        return frame

    demo = load("demo")
    pn = load("pn")
    hba1c = load("hba1c")
    biochem = load("biochem")
    bmx = load("bmx")
    diq = load("diq", required=False)
    smq = load("smq", required=False)

    if cycle in ("A", "B"):
        nutrition = load("nutrition")
        nutrient_piece = keep(
            nutrition,
            [
                "SEQN",
                "LBXB12",
                "LBDB12SI",
                "LBXMMA",
                "LBDFOLSI",
                "LBXFOL",
                "LBXHCY",
                "LBDHCY",
            ],
        )
    else:
        b12 = load("b12")
        mma = load("mma")
        nutrient_piece = keep(
            b12,
            ["SEQN", "LBXB12", "LBDB12SI", "LBDFOLSI", "LBXFOL"],
        ).merge(
            keep(mma, ["SEQN", "LBXMMA", "LBXHCY", "LBDHCY"]),
            on="SEQN",
            how="outer",
        )

    pieces = [
        keep(
            demo,
            [
                "SEQN",
                "RIAGENDR",
                "RIDAGEYR",
                "RIDRETH1",
                "RIDRETH2",
                "RIDRETH3",
                "RIDEXPRG",
                "WTMEC2YR",
                "WTMEC4YR",
                "SDMVPSU",
                "SDMVSTRA",
                "INDFMPIR",
            ],
        ),
        keep(pn, ["SEQN", "LEALPN", "LEARPN"]),
        nutrient_piece,
        keep(hba1c, ["SEQN", "LBXGH"]),
        keep(
            biochem,
            [
                "SEQN",
                "LBXSCR",
                "LBDSCR",
                "LBDSCRSI",
                "LBXSCRSI",
                "LBXSC3SI",
            ],
        ),
        keep(bmx, ["SEQN", "BMXBMI"]),
        keep(diq, ["SEQN", "DIQ010"]),
        keep(smq, ["SEQN", "SMQ020", "SMQ040"]),
    ]

    merged = pieces[0]
    for piece in pieces[1:]:
        if not piece.empty and "SEQN" in piece.columns:
            merged = merged.merge(piece.drop_duplicates("SEQN"), on="SEQN", how="left")

    merged["cycle"] = cycle
    merged["cycle_years"] = CYCLES[cycle]["years"]
    merged["cycle_index"] = CYCLES[cycle]["index"]
    merged["age"] = first_numeric(merged, ["RIDAGEYR"])
    merged["sex"] = first_numeric(merged, ["RIAGENDR"])
    merged["race"] = first_numeric(merged, ["RIDRETH1", "RIDRETH2", "RIDRETH3"])
    merged["pregnancy"] = first_numeric(merged, ["RIDEXPRG"])
    merged["weight_raw"] = first_numeric(merged, ["WTMEC4YR", "WTMEC2YR"])
    merged["weight"] = merged["weight_raw"] / 3.0
    merged["psu"] = first_numeric(merged, ["SDMVPSU"])
    merged["strata"] = first_numeric(merged, ["SDMVSTRA"])
    merged["psu_u"] = merged["cycle_index"] * 100 + merged["psu"]
    merged["strata_u"] = merged["cycle_index"] * 1000 + merged["strata"]
    merged["pir"] = first_numeric(merged, ["INDFMPIR"])
    merged["bmi"] = first_numeric(merged, ["BMXBMI"])
    merged["hba1c"] = first_numeric(merged, ["LBXGH"])
    merged["diabetes_q"] = first_numeric(merged, ["DIQ010"])
    merged["smoked100"] = first_numeric(merged, ["SMQ020"])
    merged["current_smoke"] = first_numeric(merged, ["SMQ040"])
    merged["b12_pg_ml"] = first_numeric(merged, ["LBXB12"])
    merged["b12_pmol_l"] = first_numeric(merged, ["LBDB12SI"])
    merged["mma_umol_l"] = first_numeric(merged, ["LBXMMA"])
    merged["folate_nmol_l"] = first_numeric(merged, ["LBDFOLSI"])
    merged["folate_ng_ml"] = first_numeric(merged, ["LBXFOL"])
    merged["homocysteine"] = first_numeric(merged, ["LBXHCY", "LBDHCY"])
    merged["left_insensate"] = first_numeric(merged, ["LEALPN"])
    merged["right_insensate"] = first_numeric(merged, ["LEARPN"])
    merged["creatinine_mg_dl"] = first_numeric(merged, ["LBXSCR", "LBDSCR"])
    if cycle == "A":
        merged["creatinine_mg_dl"] = 1.013 * merged["creatinine_mg_dl"] + 0.147
    return merged, audit


def egfr_2021(creatinine: pd.Series, age: pd.Series, sex: pd.Series) -> pd.Series:
    female = sex.eq(2)
    k = np.where(female, 0.7, 0.9)
    alpha = np.where(female, -0.241, -0.302)
    ratio = creatinine / k
    return 142.0 * np.minimum(ratio, 1.0) ** alpha * np.maximum(ratio, 1.0) ** -1.200 * 0.9938 ** age * np.where(female, 1.012, 1.0)


def weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    ok = values.notna() & weights.notna() & weights.gt(0)
    if not ok.any():
        return float("nan")
    return float(np.average(values.loc[ok].astype(float), weights=weights.loc[ok]))


def main() -> None:
    frames: list[pd.DataFrame] = []
    audits: list[dict[str, object]] = []
    for cycle in CYCLES:
        frame, audit = merge_cycle(cycle)
        frames.append(frame)
        audits.extend(audit)

    raw = pd.concat(frames, ignore_index=True, sort=False)
    raw["egfr"] = egfr_2021(raw["creatinine_mg_dl"], raw["age"], raw["sex"])
    raw["diabetes"] = np.where(
        raw["diabetes_q"].eq(1) | raw["hba1c"].ge(6.5),
        1,
        np.where(raw["diabetes_q"].isin([2, 3]) & raw["hba1c"].lt(6.5), 0, np.nan),
    )
    raw["smoking"] = pd.Series(pd.NA, index=raw.index, dtype="string")
    raw.loc[raw["smoked100"].eq(2), "smoking"] = "never"
    raw.loc[raw["smoked100"].eq(1) & raw["current_smoke"].isin([1, 2]), "smoking"] = "current"
    raw.loc[raw["smoked100"].eq(1) & raw["current_smoke"].eq(3), "smoking"] = "former"
    raw["pn_sites"] = raw["left_insensate"] + raw["right_insensate"]
    raw["pn_any"] = (raw["pn_sites"] >= 1).astype(float)
    raw["pn_two_plus"] = (raw["pn_sites"] >= 2).astype(float)
    raw["log_mma"] = np.log(raw["mma_umol_l"])
    raw["period"] = np.where(raw["cycle"].isin(["A", "B"]), "discovery", "validation")

    flow: list[dict[str, object]] = []

    def record(step: str, frame: pd.DataFrame) -> None:
        flow.append({"step": step, "n": int(len(frame))})

    record("all_selected_cycle_records", raw)
    cohort = raw.loc[raw["age"].ge(40)].copy()
    record("age_40_plus", cohort)
    cohort = cohort.loc[cohort["pregnancy"].ne(1) | cohort["pregnancy"].isna()].copy()
    record("exclude_known_pregnancy", cohort)
    cohort = cohort.loc[
        cohort["left_insensate"].between(0, 3, inclusive="both")
        & cohort["right_insensate"].between(0, 3, inclusive="both")
    ].copy()
    record("valid_bilateral_monofilament_exam", cohort)
    cohort = cohort.loc[
        cohort["b12_pg_ml"].between(50, 3000, inclusive="both")
        & cohort["mma_umol_l"].between(0.02, 20, inclusive="both")
        & cohort["creatinine_mg_dl"].between(0.2, 15, inclusive="both")
        & cohort["egfr"].between(5, 200, inclusive="both")
    ].copy()
    record("valid_nutrition_and_kidney_biomarkers", cohort)
    cohort = cohort.loc[cohort["egfr"].ge(60)].copy()
    record("primary_renal_function_restriction", cohort)
    cohort = cohort.loc[
        cohort["weight"].gt(0)
        & cohort["psu_u"].notna()
        & cohort["strata_u"].notna()
    ].copy()
    record("valid_survey_design", cohort)

    cohort["group"] = np.select(
        [
            cohort["b12_pg_ml"].lt(300),
            cohort["b12_pg_ml"].ge(300) & cohort["mma_umol_l"].gt(0.40),
            cohort["b12_pg_ml"].ge(300) & cohort["mma_umol_l"].le(0.40),
        ],
        ["low_b12", "discordant", "adequate"],
        default="unclassified",
    )
    cohort = cohort.loc[cohort["group"].ne("unclassified")].copy()
    record("three_group_classification", cohort)

    full_covariates = ["age", "sex", "race", "bmi", "pir", "smoking", "diabetes", "hba1c", "egfr"]
    cohort["complete_case"] = cohort[full_covariates].notna().all(axis=1)
    record("complete_case_full_model", cohort.loc[cohort["complete_case"]])

    group_rows: list[dict[str, object]] = []
    scopes = [
        ("overall", cohort),
        ("discovery", cohort.loc[cohort["period"].eq("discovery")]),
        ("validation", cohort.loc[cohort["period"].eq("validation")]),
    ]
    for scope, frame in scopes:
        for group, g in frame.groupby("group", observed=True):
            group_rows.append(
                {
                    "scope": scope,
                    "group": group,
                    "n": int(len(g)),
                    "pn_any_events": int(g["pn_any"].sum()),
                    "pn_two_plus_events": int(g["pn_two_plus"].sum()),
                    "pn_any_pct_unweighted": float(g["pn_any"].mean() * 100),
                    "pn_any_pct_weighted": weighted_mean(g["pn_any"], g["weight"]) * 100,
                    "mean_mma": float(g["mma_umol_l"].mean()),
                    "mean_b12": float(g["b12_pg_ml"].mean()),
                }
            )

    missing_rows: list[dict[str, object]] = []
    variables = [
        "b12_pg_ml", "mma_umol_l", "left_insensate", "right_insensate", "creatinine_mg_dl",
        "age", "sex", "race", "bmi", "pir", "smoking", "diabetes", "hba1c", "egfr", "folate_nmol_l",
    ]
    cycle_scopes = [("overall", cohort)] + [(f"cycle_{c}", cohort.loc[cohort["cycle"].eq(c)]) for c in CYCLES]
    for scope, frame in cycle_scopes:
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
            {"concept": "serum_b12", "field": "LBXB12", "unit": "pg/mL", "meaning": "measured serum vitamin B12", "rule": "valid 50-3000"},
            {"concept": "methylmalonic_acid", "field": "LBXMMA", "unit": "umol/L", "meaning": "measured serum methylmalonic acid", "rule": "valid 0.02-20"},
            {"concept": "peripheral_neuropathy", "field": "LEALPN+LEARPN", "unit": "0-6 insensate plantar sites", "meaning": "objective 10-g monofilament summary", "rule": "both feet valid; primary >=1 site"},
            {"concept": "creatinine", "field": "LBXSCR/LBDSCR", "unit": "mg/dL", "meaning": "serum creatinine", "rule": "1999-2000 corrected as 1.013*x+0.147"},
            {"concept": "egfr", "field": "derived", "unit": "mL/min/1.73m2", "meaning": "2021 CKD-EPI race-free eGFR", "rule": "primary >=60"},
            {"concept": "survey_design", "field": "WTMEC2YR/WTMEC4YR, SDMVPSU, SDMVSTRA", "unit": "survey", "meaning": "MEC weights and design", "rule": "pooled weights divided by 3; cycle-unique PSU/strata"},
        ]
    )

    columns = [
        "SEQN", "cycle", "cycle_years", "period", "age", "sex", "race", "pregnancy",
        "weight", "psu_u", "strata_u", "pir", "bmi", "smoking", "diabetes", "hba1c",
        "b12_pg_ml", "b12_pmol_l", "mma_umol_l", "log_mma", "folate_nmol_l", "folate_ng_ml",
        "homocysteine", "creatinine_mg_dl", "egfr", "left_insensate", "right_insensate",
        "pn_sites", "pn_any", "pn_two_plus", "group", "complete_case",
    ]
    cohort.loc[:, columns].to_csv(OUT / "c_n06_analysis_core.csv", index=False)
    pd.DataFrame(audits).to_csv(OUT / "c_n06_source_audit.csv", index=False)
    pd.DataFrame(flow).to_csv(OUT / "c_n06_flow.csv", index=False)
    pd.DataFrame(group_rows).to_csv(OUT / "c_n06_group_counts.csv", index=False)
    pd.DataFrame(missing_rows).to_csv(OUT / "c_n06_missingness.csv", index=False)
    semantic.to_csv(OUT / "c_n06_semantic_audit.csv", index=False)

    key = cohort.loc[cohort["group"].eq("discordant")]
    validation_key = key.loc[key["period"].eq("validation")]
    status = {
        "candidate_code": "C-N06",
        "actual_n": int(len(cohort)),
        "complete_case_n": int(cohort["complete_case"].sum()),
        "complete_case_retention": float(cohort["complete_case"].mean()),
        "primary_events": int(cohort["pn_any"].sum()),
        "key_group_n": int(len(key)),
        "key_group_events": int(key["pn_any"].sum()),
        "validation_key_group_n": int(len(validation_key)),
        "validation_key_group_events": int(validation_key["pn_any"].sum()),
        "cycles": sorted(cohort["cycle"].unique().tolist()),
    }
    (OUT / "c_n06_prep_status.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
    print(json.dumps(status, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
