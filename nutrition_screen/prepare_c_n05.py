from __future__ import annotations

import io
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import requests

BASE = "https://raw.githubusercontent.com/protobi/nhanes-continuous/main"
OUT = Path("nutrition_screen/output")
OUT.mkdir(parents=True, exist_ok=True)

CYCLES = {
    "C": {"years": "2003-2004", "index": 1},
    "D": {"years": "2005-2006", "index": 2},
    "E": {"years": "2007-2008", "index": 3},
    "F": {"years": "2009-2010", "index": 4},
}


def fetch_csv(path: str, *, required: bool = True) -> pd.DataFrame:
    url = f"{BASE}/{path}"
    response = requests.get(url, timeout=120)
    if response.status_code != 200:
        if required:
            raise RuntimeError(f"Required source unavailable: {path}; HTTP {response.status_code}")
        return pd.DataFrame(columns=["SEQN"])
    try:
        return pd.read_csv(io.BytesIO(response.content), low_memory=False)
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(f"Cannot parse {path}: {type(exc).__name__}: {exc}") from exc


def keep(df: pd.DataFrame, columns: Iterable[str]) -> pd.DataFrame:
    cols = [c for c in columns if c in df.columns]
    if "SEQN" not in cols and "SEQN" in df.columns:
        cols.insert(0, "SEQN")
    return df.loc[:, cols].copy()


def first_existing(df: pd.DataFrame, names: Iterable[str]) -> pd.Series:
    out = pd.Series(np.nan, index=df.index, dtype="float64")
    for name in names:
        if name in df.columns:
            out = out.fillna(pd.to_numeric(df[name], errors="coerce"))
    return out


def merge_one(cycle: str) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    sources: list[dict[str, object]] = []

    def load(path: str, required: bool = True) -> pd.DataFrame:
        df = fetch_csv(path, required=required)
        sources.append(
            {
                "cycle": cycle,
                "path": path,
                "required": required,
                "available": not df.empty,
                "rows": int(len(df)),
                "columns": int(df.shape[1]),
            }
        )
        return df

    demo = load(f"Demographics/DEMO_{cycle}.csv")
    bmx = load(f"Examination/BMX_{cycle}.csv")
    smq = load(f"Questionnaire/SMQ_{cycle}.csv", required=False)
    diq = load(f"Questionnaire/DIQ_{cycle}.csv", required=False)

    if cycle == "C":
        iron = load("Laboratory/L06TFR_C.csv")
        cbc = load("Laboratory/L25_C.csv")
        crp = load("Laboratory/L11_C.csv")
        lab = keep(iron, ["SEQN", "LBXTFR", "LBDFER", "LBDFERSI"])
    else:
        tfr = load(f"Laboratory/TFR_{cycle}.csv")
        fer = load(f"Laboratory/FERTIN_{cycle}.csv")
        cbc = load(f"Laboratory/CBC_{cycle}.csv")
        crp = load(f"Laboratory/CRP_{cycle}.csv")
        lab = keep(tfr, ["SEQN", "LBXTFR", "LBDTFRSI"]).merge(
            keep(fer, ["SEQN", "LBXFER", "LBDFERSI"]),
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
                "RIDRETH3",
                "RIDEXPRG",
                "WTMEC2YR",
                "SDMVPSU",
                "SDMVSTRA",
                "INDFMPIR",
            ],
        ),
        lab,
        keep(cbc, ["SEQN", "LBXHGB", "LBXMCVSI", "LBXMCV", "LBXMCHSI", "LBXRDW"]),
        keep(crp, ["SEQN", "LBXCRP", "LBDHRPLC"]),
        keep(bmx, ["SEQN", "BMXBMI"]),
        keep(smq, ["SEQN", "SMQ020", "SMQ040"]),
        keep(diq, ["SEQN", "DIQ010"]),
    ]

    merged = pieces[0]
    for piece in pieces[1:]:
        if not piece.empty and "SEQN" in piece.columns:
            merged = merged.merge(piece.drop_duplicates("SEQN"), on="SEQN", how="left")

    merged["cycle"] = cycle
    merged["cycle_years"] = CYCLES[cycle]["years"]
    merged["cycle_index"] = CYCLES[cycle]["index"]
    merged["stfr"] = first_existing(merged, ["LBXTFR"])
    merged["ferritin"] = first_existing(merged, ["LBXFER", "LBDFER", "LBDFERSI"])
    merged["hgb"] = first_existing(merged, ["LBXHGB"])
    merged["mcv"] = first_existing(merged, ["LBXMCVSI", "LBXMCV"])
    merged["crp_mg_dl"] = first_existing(merged, ["LBXCRP"])
    merged["crp_mg_l"] = merged["crp_mg_dl"] * 10.0
    merged["bmi"] = first_existing(merged, ["BMXBMI"])
    merged["age"] = first_existing(merged, ["RIDAGEYR"])
    merged["sex"] = first_existing(merged, ["RIAGENDR"])
    merged["race"] = first_existing(merged, ["RIDRETH1", "RIDRETH3"])
    merged["pregnancy"] = first_existing(merged, ["RIDEXPRG"])
    merged["weight"] = first_existing(merged, ["WTMEC2YR"]) / len(CYCLES)
    merged["psu"] = first_existing(merged, ["SDMVPSU"])
    merged["strata"] = first_existing(merged, ["SDMVSTRA"])
    merged["pir"] = first_existing(merged, ["INDFMPIR"])
    merged["smoked100"] = first_existing(merged, ["SMQ020"])
    merged["diabetes_q"] = first_existing(merged, ["DIQ010"])
    merged["psu_u"] = merged["cycle_index"] * 100 + merged["psu"]
    merged["strata_u"] = merged["cycle_index"] * 1000 + merged["strata"]
    return merged, sources


def weighted_percent(values: pd.Series, weight: pd.Series) -> float:
    values_numeric = pd.to_numeric(values, errors="coerce")
    ok = values_numeric.notna() & weight.notna() & (weight > 0)
    if not ok.any():
        return float("nan")
    return float(np.average(values_numeric.loc[ok], weights=weight.loc[ok]) * 100.0)


def main() -> None:
    all_frames: list[pd.DataFrame] = []
    all_sources: list[dict[str, object]] = []
    for cycle in CYCLES:
        frame, sources = merge_one(cycle)
        all_frames.append(frame)
        all_sources.extend(sources)

    raw = pd.concat(all_frames, ignore_index=True, sort=False)
    flow: list[dict[str, object]] = []

    def add_flow(step: str, frame: pd.DataFrame) -> None:
        flow.append({"step": step, "n": int(len(frame))})

    add_flow("all_examined_records_in_selected_cycles", raw)
    cohort = raw.loc[(raw["sex"] == 2) & raw["age"].between(20, 49, inclusive="both")].copy()
    add_flow("women_age_20_49", cohort)
    cohort = cohort.loc[cohort["pregnancy"].ne(1) | cohort["pregnancy"].isna()].copy()
    add_flow("exclude_known_pregnancy", cohort)
    cohort = cohort.loc[
        cohort["stfr"].between(0.2, 30, inclusive="both")
        & cohort["ferritin"].between(1, 2000, inclusive="both")
        & cohort["hgb"].between(5, 20, inclusive="both")
    ].copy()
    add_flow("valid_core_iron_and_hemoglobin", cohort)
    cohort = cohort.loc[
        cohort["weight"].gt(0)
        & cohort["psu_u"].notna()
        & cohort["strata_u"].notna()
    ].copy()
    add_flow("valid_survey_design", cohort)

    cohort["group"] = np.select(
        [
            cohort["ferritin"] < 15,
            (cohort["ferritin"] >= 15) & (cohort["stfr"] > 5.33),
            (cohort["ferritin"] >= 15) & (cohort["stfr"] <= 5.33),
        ],
        ["ferritin_low", "discordant", "adequate"],
        default="unclassified",
    )
    cohort = cohort.loc[cohort["group"].ne("unclassified")].copy()
    add_flow("classified_three_group_cohort", cohort)

    cohort["anemia"] = (cohort["hgb"] < 12.0).astype(int)
    cohort["low_mcv"] = np.where(
        cohort["mcv"].notna(),
        (cohort["mcv"] < 80.0).astype(int),
        np.nan,
    )
    cohort["log_crp"] = np.log(cohort["crp_mg_l"].clip(lower=0.01))
    cohort["log_stfr"] = np.log(cohort["stfr"])
    cohort["smoking"] = pd.Series(pd.NA, index=cohort.index, dtype="string")
    cohort.loc[cohort["smoked100"].eq(1), "smoking"] = "ever"
    cohort.loc[cohort["smoked100"].eq(2), "smoking"] = "never"
    cohort["diabetes"] = np.where(
        cohort["diabetes_q"].eq(1),
        1,
        np.where(cohort["diabetes_q"].isin([2, 3]), 0, np.nan),
    )
    cohort["period"] = np.where(
        cohort["cycle"].isin(["C", "D"]),
        "discovery",
        "validation",
    )

    covariates = ["bmi", "crp_mg_l", "pir", "smoking", "race"]
    cohort["complete_case"] = cohort[covariates].notna().all(axis=1)
    add_flow("complete_case_full_adjustment", cohort.loc[cohort["complete_case"]])

    missing_rows: list[dict[str, object]] = []
    scopes = [("overall", cohort)] + [
        (f"cycle_{cycle}", cohort.loc[cohort["cycle"].eq(cycle)])
        for cycle in CYCLES
    ]
    for variable in [
        "stfr",
        "ferritin",
        "hgb",
        "mcv",
        "bmi",
        "crp_mg_l",
        "pir",
        "smoking",
        "diabetes",
    ]:
        for scope_name, frame in scopes:
            missing_rows.append(
                {
                    "scope": scope_name,
                    "variable": variable,
                    "n": int(len(frame)),
                    "missing_n": int(frame[variable].isna().sum()),
                    "missing_pct": float(frame[variable].isna().mean() * 100)
                    if len(frame)
                    else np.nan,
                }
            )

    group_rows: list[dict[str, object]] = []
    group_scopes = [
        ("overall", cohort),
        ("discovery", cohort.loc[cohort["period"].eq("discovery")]),
        ("validation", cohort.loc[cohort["period"].eq("validation")]),
    ]
    for scope_name, frame in group_scopes:
        for group_name, group_frame in frame.groupby("group", observed=True):
            group_rows.append(
                {
                    "scope": scope_name,
                    "group": group_name,
                    "n": int(len(group_frame)),
                    "anemia_events": int(group_frame["anemia"].sum()),
                    "anemia_pct_unweighted": float(group_frame["anemia"].mean() * 100),
                    "anemia_pct_weighted": weighted_percent(
                        group_frame["anemia"],
                        group_frame["weight"],
                    ),
                    "mean_hgb_unweighted": float(group_frame["hgb"].mean()),
                    "mean_mcv_unweighted": float(group_frame["mcv"].mean()),
                    "weight_sum": float(group_frame["weight"].sum()),
                }
            )

    source_df = pd.DataFrame(all_sources)
    flow_df = pd.DataFrame(flow)
    missing_df = pd.DataFrame(missing_rows)
    group_df = pd.DataFrame(group_rows)

    semantic = pd.DataFrame(
        [
            {
                "variable": "stfr",
                "source": "L06TFR_C or TFR_D/E/F",
                "released_field": "LBXTFR",
                "semantics": "serum soluble transferrin receptor",
                "unit": "mg/L",
                "range_rule": "0.2-30",
            },
            {
                "variable": "ferritin",
                "source": "L06TFR_C or FERTIN_D/E/F",
                "released_field": "LBDFER/LBXFER",
                "semantics": "serum ferritin",
                "unit": "ng/mL",
                "range_rule": "1-2000",
            },
            {
                "variable": "hgb",
                "source": "L25_C or CBC_D/E/F",
                "released_field": "LBXHGB",
                "semantics": "measured hemoglobin",
                "unit": "g/dL",
                "range_rule": "5-20",
            },
            {
                "variable": "mcv",
                "source": "L25_C or CBC_D/E/F",
                "released_field": "LBXMCVSI/LBXMCV",
                "semantics": "mean corpuscular volume",
                "unit": "fL",
                "range_rule": "no post-hoc truncation",
            },
            {
                "variable": "crp_mg_l",
                "source": "L11_C or CRP_D/E/F",
                "released_field": "LBXCRP",
                "semantics": "C-reactive protein",
                "unit": "released mg/dL converted x10 to mg/L",
                "range_rule": "log transform",
            },
            {
                "variable": "survey_design",
                "source": "DEMO_C/D/E/F",
                "released_field": "WTMEC2YR/SDMVPSU/SDMVSTRA",
                "semantics": "MEC weight, PSU and strata",
                "unit": "pooled weight divided by four",
                "range_rule": "cycle-unique PSU and stratum IDs",
            },
        ]
    )

    out_cols = [
        "SEQN",
        "cycle",
        "cycle_years",
        "period",
        "age",
        "race",
        "pregnancy",
        "weight",
        "psu_u",
        "strata_u",
        "stfr",
        "ferritin",
        "hgb",
        "mcv",
        "crp_mg_l",
        "log_crp",
        "bmi",
        "pir",
        "smoking",
        "diabetes",
        "group",
        "anemia",
        "low_mcv",
        "log_stfr",
        "complete_case",
    ]
    cohort.loc[:, out_cols].to_csv(OUT / "c_n05_analysis_core.csv", index=False)
    source_df.to_csv(OUT / "c_n05_source_audit.csv", index=False)
    flow_df.to_csv(OUT / "c_n05_flow.csv", index=False)
    missing_df.to_csv(OUT / "c_n05_missingness.csv", index=False)
    group_df.to_csv(OUT / "c_n05_group_counts.csv", index=False)
    semantic.to_csv(OUT / "c_n05_semantic_audit.csv", index=False)

    required_failures = source_df.loc[
        source_df["required"].astype(bool) & ~source_df["available"].astype(bool),
        "path",
    ].tolist()
    prep_status = {
        "candidate_code": "C-N05",
        "core_n": int(len(cohort)),
        "complete_case_n": int(cohort["complete_case"].sum()),
        "complete_case_retention": float(cohort["complete_case"].mean()),
        "key_group_n": int(cohort["group"].eq("discordant").sum()),
        "key_group_anemia_events": int(
            cohort.loc[cohort["group"].eq("discordant"), "anemia"].sum()
        ),
        "cycles": sorted(cohort["cycle"].unique().tolist()),
        "source_failures": required_failures,
    }
    (OUT / "c_n05_prep_status.json").write_text(
        json.dumps(prep_status, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(prep_status, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
