from __future__ import annotations

import io
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import requests

OUT = Path("nutrition_stage/output_c_n12")
OUT.mkdir(parents=True, exist_ok=True)

CYCLES = {
    "A": {"year": "1999", "years": "1999-2000", "index": 1, "period": "discovery"},
    "B": {"year": "2001", "years": "2001-2002", "index": 2, "period": "discovery"},
    "C": {"year": "2003", "years": "2003-2004", "index": 3, "period": "validation"},
}

FILES = {
    "A": {
        "demo": ["DEMO"],
        "pn": ["LEXPN"],
        "b12": ["LAB06"],
        "hba1c": ["LAB10"],
        "biochem": ["LAB18"],
        "body": ["BMX"],
        "diabetes": ["DIQ"],
        "smoking": ["SMQ"],
    },
    "B": {
        "demo": ["DEMO_B"],
        "pn": ["LEXPN_B"],
        "b12": ["L06_B"],
        "hba1c": ["L10_B"],
        "biochem": ["L40_B"],
        "body": ["BMX_B"],
        "diabetes": ["DIQ_B"],
        "smoking": ["SMQ_B"],
    },
    "C": {
        "demo": ["DEMO_C"],
        "pn": ["LEXPN_C"],
        "b12": ["L06NB_C"],
        "mma": ["L06MH_C"],
        "hba1c": ["L10_C"],
        "biochem": ["L40_C"],
        "body": ["BMX_C"],
        "diabetes": ["DIQ_C"],
        "smoking": ["SMQ_C"],
    },
}


def xpt_url(cycle: str, stem: str) -> str:
    return f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{CYCLES[cycle]['year']}/DataFiles/{stem}.xpt"


def read_xpt(cycle: str, stems: list[str], required: bool = True) -> tuple[pd.DataFrame, str, str | None]:
    errors: list[str] = []
    for stem in stems:
        url = xpt_url(cycle, stem)
        try:
            response = requests.get(url, timeout=180)
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


def pooled_weight(cycle: str, frame: pd.DataFrame) -> pd.Series:
    """Correct six-year combination: 1999-2002 four-year weights contribute 4/6; 2003-04 contributes 2/6."""
    if cycle in {"A", "B"} and "WTMEC4YR" in frame.columns:
        return pd.to_numeric(frame["WTMEC4YR"], errors="coerce") * (4.0 / 6.0)
    return numeric(frame, ["WTMEC2YR"]) * (2.0 / 6.0)


def weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    values = pd.to_numeric(values, errors="coerce")
    valid = values.notna() & weights.notna() & weights.gt(0)
    if not valid.any():
        return float("nan")
    return float(np.average(values.loc[valid], weights=weights.loc[valid]))


def merge_cycle(cycle: str) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    audits: list[dict[str, object]] = []
    frames: dict[str, pd.DataFrame] = {}

    for component, stems in FILES[cycle].items():
        frame, selected, error = read_xpt(cycle, stems, required=True)
        frames[component] = frame
        audits.append({
            "cycle": cycle,
            "component": component,
            "selected_file": selected,
            "required": True,
            "rows": int(frame.shape[0]),
            "columns": int(frame.shape[1]),
            "column_names": "|".join(frame.columns),
            "error": error,
        })

    if cycle == "C":
        b12_piece = keep(frames["b12"], ["SEQN", "LBXB12", "LBDB12SI", "LBDFOLSI", "LBXFOL"])
        b12_piece = b12_piece.merge(keep(frames["mma"], ["SEQN", "LBXMMA"]), on="SEQN", how="outer")
    else:
        b12_piece = keep(frames["b12"], ["SEQN", "LBXB12", "LBDB12SI", "LBXMMA", "LBDFOLSI", "LBXFOL"])

    pieces = [
        keep(frames["demo"], [
            "SEQN", "RIAGENDR", "RIDAGEYR", "RIDRETH1", "RIDRETH2", "RIDRETH3", "RIDEXPRG",
            "WTMEC2YR", "WTMEC4YR", "SDMVPSU", "SDMVSTRA", "INDFMPIR", "DMDEDUC2",
        ]),
        keep(frames["pn"], ["SEQN", "LEALPN", "LEARPN"]),
        b12_piece,
        keep(frames["hba1c"], ["SEQN", "LBXGH"]),
        keep(frames["biochem"], ["SEQN", "LBXSCR", "LBDSCR", "LBDSCRSI", "LBXSAL", "LBXSTP"]),
        keep(frames["body"], ["SEQN", "BMXBMI", "BMXWT", "BMXHT"]),
        keep(frames["diabetes"], ["SEQN", "DIQ010"]),
        keep(frames["smoking"], ["SEQN", "SMQ020", "SMQ040"]),
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
    merged["race"] = numeric(merged, ["RIDRETH1", "RIDRETH2", "RIDRETH3"])
    merged["pregnancy"] = numeric(merged, ["RIDEXPRG"])
    merged["education"] = numeric(merged, ["DMDEDUC2"])
    merged["pir"] = numeric(merged, ["INDFMPIR"])
    merged["weight"] = pooled_weight(cycle, merged)
    merged["psu"] = numeric(merged, ["SDMVPSU"])
    merged["strata"] = numeric(merged, ["SDMVSTRA"])
    merged["psu_u"] = merged["cycle_index"] * 100 + merged["psu"]
    merged["strata_u"] = merged["cycle_index"] * 1000 + merged["strata"]
    merged["b12_pg_ml"] = numeric(merged, ["LBXB12"])
    merged["mma_umol_l"] = numeric(merged, ["LBXMMA"])
    merged["folate_nmol_l"] = numeric(merged, ["LBDFOLSI"])
    merged["folate_ng_ml"] = numeric(merged, ["LBXFOL"])
    merged["hba1c"] = numeric(merged, ["LBXGH"])
    merged["creatinine_mg_dl"] = numeric(merged, ["LBXSCR", "LBDSCR"])
    if cycle == "A":
        merged["creatinine_mg_dl"] = 1.013 * merged["creatinine_mg_dl"] + 0.147
    merged["egfr"] = egfr_2021(merged["creatinine_mg_dl"], merged["age"], merged["sex"])
    merged["albumin_g_dl"] = numeric(merged, ["LBXSAL"])
    merged["total_protein_g_dl"] = numeric(merged, ["LBXSTP"])
    merged["bmi"] = numeric(merged, ["BMXBMI"])
    merged["body_weight_kg"] = numeric(merged, ["BMXWT"])
    merged["height_cm"] = numeric(merged, ["BMXHT"])
    merged["diabetes_q"] = numeric(merged, ["DIQ010"])
    merged["smoked100"] = numeric(merged, ["SMQ020"])
    merged["current_smoke"] = numeric(merged, ["SMQ040"])
    merged["left_insensate"] = numeric(merged, ["LEALPN"])
    merged["right_insensate"] = numeric(merged, ["LEARPN"])
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
    raw["diabetes"] = np.where(
        raw["diabetes_q"].eq(1) | raw["hba1c"].ge(6.5),
        1.0,
        np.where(raw["diabetes_q"].isin([2, 3]) & raw["hba1c"].lt(6.5), 0.0, np.nan),
    )
    raw["pn_sites"] = raw["left_insensate"] + raw["right_insensate"]
    raw["pn_any"] = (raw["pn_sites"] >= 1).astype(float)
    raw["pn_two_plus"] = (raw["pn_sites"] >= 2).astype(float)
    raw["log_mma"] = np.log(raw["mma_umol_l"])
    raw["log_b12"] = np.log(raw["b12_pg_ml"])

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
        cohort["b12_pg_ml"].between(100, 3000, inclusive="both")
        & cohort["mma_umol_l"].between(0.02, 2.0, inclusive="both")
        & cohort["folate_nmol_l"].ge(7.0)
    ].copy()
    record("valid_b12_mma_and_nondeficient_folate", cohort)
    cohort = cohort.loc[
        cohort["egfr"].ge(60)
        & cohort["albumin_g_dl"].ge(3.5)
        & cohort["bmi"].between(15, 70, inclusive="both")
    ].copy()
    record("exclude_ckd_hypoalbuminemia_and_extreme_bmi", cohort)
    cohort = cohort.loc[cohort["diabetes"].eq(0)].copy()
    record("exclude_diabetes", cohort)
    cohort = cohort.loc[
        cohort["weight"].gt(0)
        & cohort["psu_u"].notna()
        & cohort["strata_u"].notna()
    ].copy()
    record("valid_survey_design", cohort)

    cohort["group"] = np.select(
        [
            cohort["b12_pg_ml"].ge(200) & cohort["mma_umol_l"].gt(0.26),
            cohort["b12_pg_ml"].ge(300) & cohort["mma_umol_l"].le(0.26),
            cohort["b12_pg_ml"].lt(200),
        ],
        ["metabolic_discordance", "concordant_adequate", "serum_low"],
        default="intermediate",
    )
    cohort = cohort.loc[cohort["group"].isin(["metabolic_discordance", "concordant_adequate"])].copy()
    record("frozen_discordant_vs_concordant_comparison", cohort)

    full_covariates = [
        "age", "sex", "race", "bmi", "smoking", "hba1c", "egfr",
        "folate_nmol_l", "pir", "albumin_g_dl",
    ]
    cohort["complete_case"] = cohort[full_covariates].notna().all(axis=1)
    record("complete_case_full_adjustment", cohort.loc[cohort["complete_case"]])

    group_rows: list[dict[str, object]] = []
    for scope, frame in [
        ("overall", cohort),
        ("discovery", cohort.loc[cohort["period"].eq("discovery")]),
        ("validation", cohort.loc[cohort["period"].eq("validation")]),
    ]:
        for group_name, group in frame.groupby("group", observed=True):
            group_rows.append({
                "scope": scope,
                "group": group_name,
                "n": int(len(group)),
                "pn_any_events": int(group["pn_any"].sum()),
                "pn_two_plus_events": int(group["pn_two_plus"].sum()),
                "weighted_pn_any_pct": 100.0 * weighted_mean(group["pn_any"], group["weight"]),
                "weighted_pn_two_plus_pct": 100.0 * weighted_mean(group["pn_two_plus"], group["weight"]),
                "mean_b12_pg_ml": float(group["b12_pg_ml"].mean()),
                "mean_mma_umol_l": float(group["mma_umol_l"].mean()),
                "mean_egfr": float(group["egfr"].mean()),
            })

    missing_rows: list[dict[str, object]] = []
    missing_vars = [
        "b12_pg_ml", "mma_umol_l", "folate_nmol_l", "left_insensate", "right_insensate",
        "hba1c", "egfr", "albumin_g_dl", "bmi", "smoking", "pir",
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
        {"concept": "serum_b12", "fields": "LBXB12", "meaning": "measured serum vitamin B12", "unit": "pg/mL", "rule": "valid 100-3000; non-low >=200"},
        {"concept": "functional_b12_marker", "fields": "LBXMMA", "meaning": "measured serum methylmalonic acid", "unit": "umol/L", "rule": "elevated >0.26; valid 0.02-2.0"},
        {"concept": "peripheral_neuropathy", "fields": "LEALPN+LEARPN", "meaning": "objective 10-g monofilament insensate plantar sites", "unit": "0-6", "rule": "both feet valid; primary >=1; strict >=2"},
        {"concept": "major_alternative_causes", "fields": "LBXSCR/LBDSCR, LBXGH, DIQ010, folate", "meaning": "kidney function, diabetes and folate status", "unit": "mixed", "rule": "eGFR>=60; diabetes excluded; folate>=7 nmol/L"},
        {"concept": "survey_design", "fields": "WTMEC4YR/WTMEC2YR, SDMVPSU, SDMVSTRA", "meaning": "MEC weights and design", "unit": "survey", "rule": "A/B weight x4/6; C weight x2/6; cycle-unique PSU/strata"},
    ])

    columns = [
        "SEQN", "cycle", "cycle_years", "period", "age", "sex", "race", "education", "pir",
        "weight", "psu_u", "strata_u", "b12_pg_ml", "mma_umol_l", "log_b12", "log_mma",
        "folate_nmol_l", "folate_ng_ml", "hba1c", "creatinine_mg_dl", "egfr", "albumin_g_dl",
        "total_protein_g_dl", "bmi", "body_weight_kg", "height_cm", "smoking", "diabetes",
        "left_insensate", "right_insensate", "pn_sites", "pn_any", "pn_two_plus", "group", "complete_case",
    ]
    cohort.loc[:, columns].to_csv(OUT / "c_n12_analysis_core.csv", index=False)
    pd.DataFrame(audits).to_csv(OUT / "c_n12_source_audit.csv", index=False)
    pd.DataFrame(flow).to_csv(OUT / "c_n12_flow.csv", index=False)
    pd.DataFrame(group_rows).to_csv(OUT / "c_n12_group_counts.csv", index=False)
    pd.DataFrame(missing_rows).to_csv(OUT / "c_n12_missingness.csv", index=False)
    semantic.to_csv(OUT / "c_n12_semantic_audit.csv", index=False)

    key = cohort.loc[cohort["group"].eq("metabolic_discordance")]
    val_key = key.loc[key["period"].eq("validation")]
    prep = {
        "candidate_code": "C-N12",
        "actual_n": int(len(cohort)),
        "complete_case_n": int(cohort["complete_case"].sum()),
        "complete_case_retention": float(cohort["complete_case"].mean()),
        "primary_events": int(cohort["pn_any"].sum()),
        "key_group_n": int(len(key)),
        "key_group_events": int(key["pn_any"].sum()),
        "validation_key_group_n": int(len(val_key)),
        "validation_key_group_events": int(val_key["pn_any"].sum()),
        "cycles": sorted(cohort["cycle"].unique().tolist()),
    }
    (OUT / "c_n12_prep_status.json").write_text(json.dumps(prep, indent=2), encoding="utf-8")
    print(json.dumps(prep, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
