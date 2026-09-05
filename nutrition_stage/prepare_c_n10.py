from __future__ import annotations

import io
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import requests

OUT = Path("nutrition_stage/output_c_n10")
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
        "vit_e": ["LAB06"],
        "lipids": ["LAB13"],
        "hba1c": ["LAB10"],
        "biochem": ["LAB18"],
        "body": ["BMX"],
        "diabetes": ["DIQ"],
        "smoking": ["SMQ"],
    },
    "B": {
        "demo": ["DEMO_B"],
        "pn": ["LEXPN_B"],
        "vit_e": ["L06VIT_B"],
        "b12": ["L06_B"],
        "lipids": ["L13_B"],
        "hba1c": ["L10_B"],
        "biochem": ["L40_B"],
        "body": ["BMX_B"],
        "diabetes": ["DIQ_B"],
        "smoking": ["SMQ_B"],
    },
    "C": {
        "demo": ["DEMO_C"],
        "pn": ["LEXPN_C"],
        "vit_e": ["L45VIT_C"],
        "b12": ["L06NB_C"],
        "mma": ["L06MH_C"],
        "lipids": ["L13_C"],
        "hba1c": ["L10_C"],
        "biochem": ["L40_C"],
        "body": ["BMX_C"],
        "diabetes": ["DIQ_C"],
        "smoking": ["SMQ_C"],
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
    return 142.0 * np.minimum(ratio, 1.0) ** alpha * np.maximum(ratio, 1.0) ** -1.200 * 0.9938 ** age * np.where(female, 1.012, 1.0)


def weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    values = pd.to_numeric(values, errors="coerce")
    valid = values.notna() & weights.notna() & weights.gt(0)
    if not valid.any():
        return float("nan")
    return float(np.average(values.loc[valid], weights=weights.loc[valid]))


def merge_cycle(cycle: str) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    audit: list[dict[str, object]] = []

    def load(key: str, required: bool = True) -> pd.DataFrame:
        frame, selected, error = read_xpt(cycle, FILES[cycle][key], required=required)
        audit.append({
            "cycle": cycle,
            "component": key,
            "selected_file": selected,
            "required": required,
            "rows": int(frame.shape[0]),
            "columns": int(frame.shape[1]),
            "error": error,
            "column_names": "|".join(frame.columns),
        })
        return frame

    demo = load("demo")
    pn = load("pn")
    vit_e = load("vit_e")
    lipids = load("lipids")
    hba1c = load("hba1c")
    biochem = load("biochem")
    body = load("body")
    diabetes = load("diabetes", required=False)
    smoking = load("smoking", required=False)

    if cycle == "A":
        b12_piece = keep(vit_e, ["SEQN", "LBXB12", "LBDB12SI", "LBXMMA", "LBDFOLSI", "LBXFOL"])
    elif cycle == "B":
        b12 = load("b12")
        b12_piece = keep(b12, ["SEQN", "LBXB12", "LBDB12SI", "LBXMMA", "LBDFOLSI", "LBXFOL"])
    else:
        b12 = load("b12")
        mma = load("mma")
        b12_piece = keep(b12, ["SEQN", "LBXB12", "LBDB12SI", "LBDFOLSI", "LBXFOL"]).merge(
            keep(mma, ["SEQN", "LBXMMA"]), on="SEQN", how="outer"
        )

    pieces = [
        keep(demo, [
            "SEQN", "RIAGENDR", "RIDAGEYR", "RIDRETH1", "RIDRETH2", "RIDRETH3", "RIDEXPRG",
            "WTMEC2YR", "WTMEC4YR", "SDMVPSU", "SDMVSTRA", "INDFMPIR", "DMDEDUC2",
        ]),
        keep(pn, ["SEQN", "LEALPN", "LEARPN"]),
        keep(vit_e, ["SEQN", "LBXVIE", "LBDVIESI", "LBXATC", "LBDATCSI"]),
        b12_piece,
        keep(lipids, ["SEQN", "LBXTC", "LBDTCSI", "LBXHDD", "LBDHDL"]),
        keep(hba1c, ["SEQN", "LBXGH"]),
        keep(biochem, ["SEQN", "LBXSCR", "LBDSCR", "LBDSCRSI", "LBXSAL", "LBXSTP"]),
        keep(body, ["SEQN", "BMXBMI", "BMXWT", "BMXHT"]),
        keep(diabetes, ["SEQN", "DIQ010"]),
        keep(smoking, ["SEQN", "SMQ020", "SMQ040"]),
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
    merged["weight_raw"] = numeric(merged, ["WTMEC4YR", "WTMEC2YR"])
    merged["weight"] = merged["weight_raw"] / len(CYCLES)
    merged["psu"] = numeric(merged, ["SDMVPSU"])
    merged["strata"] = numeric(merged, ["SDMVSTRA"])
    merged["psu_u"] = merged["cycle_index"] * 100 + merged["psu"]
    merged["strata_u"] = merged["cycle_index"] * 1000 + merged["strata"]
    merged["alpha_tocopherol_umol_l"] = numeric(merged, ["LBDVIESI", "LBDATCSI"])
    merged["total_cholesterol_mg_dl"] = numeric(merged, ["LBXTC"])
    merged["total_cholesterol_mmol_l"] = numeric(merged, ["LBDTCSI"])
    merged.loc[merged["total_cholesterol_mmol_l"].isna(), "total_cholesterol_mmol_l"] = (
        merged.loc[merged["total_cholesterol_mmol_l"].isna(), "total_cholesterol_mg_dl"] / 38.67
    )
    merged["hdl_mg_dl"] = numeric(merged, ["LBXHDD", "LBDHDL"])
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
    return merged, audit


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
    raw["vit_e_chol_ratio"] = raw["alpha_tocopherol_umol_l"] / raw["total_cholesterol_mmol_l"]
    raw["log_ratio"] = np.log(raw["vit_e_chol_ratio"])
    raw["log_alpha_tocopherol"] = np.log(raw["alpha_tocopherol_umol_l"])
    raw["log_mma"] = np.log(raw["mma_umol_l"])

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
        cohort["alpha_tocopherol_umol_l"].between(2, 100, inclusive="both")
        & cohort["total_cholesterol_mmol_l"].between(1.5, 15, inclusive="both")
        & cohort["vit_e_chol_ratio"].between(0.2, 15, inclusive="both")
    ].copy()
    record("valid_vitamin_e_and_cholesterol", cohort)
    cohort = cohort.loc[
        cohort["b12_pg_ml"].between(200, 3000, inclusive="both")
        & cohort["mma_umol_l"].between(0.02, 0.40, inclusive="both")
        & cohort["egfr"].ge(60)
    ].copy()
    record("exclude_b12_functional_deficiency_and_ckd", cohort)
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
            cohort["alpha_tocopherol_umol_l"].lt(12),
            cohort["alpha_tocopherol_umol_l"].ge(12) & cohort["vit_e_chol_ratio"].lt(2.2),
            cohort["alpha_tocopherol_umol_l"].ge(12) & cohort["vit_e_chol_ratio"].ge(2.2),
        ],
        ["absolute_low", "ratio_low_absolute_normal", "concordant_adequate"],
        default="unclassified",
    )
    cohort = cohort.loc[cohort["group"].ne("unclassified")].copy()
    record("three_group_classification", cohort)

    full_covariates = [
        "age", "sex", "race", "bmi", "smoking", "hba1c", "egfr", "b12_pg_ml",
        "log_mma", "folate_nmol_l", "pir", "albumin_g_dl",
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
                "mean_alpha_tocopherol_umol_l": float(group["alpha_tocopherol_umol_l"].mean()),
                "mean_vit_e_chol_ratio": float(group["vit_e_chol_ratio"].mean()),
                "mean_total_cholesterol_mmol_l": float(group["total_cholesterol_mmol_l"].mean()),
            })

    missing_rows: list[dict[str, object]] = []
    missing_vars = [
        "alpha_tocopherol_umol_l", "total_cholesterol_mmol_l", "vit_e_chol_ratio", "left_insensate",
        "right_insensate", "b12_pg_ml", "mma_umol_l", "folate_nmol_l", "egfr", "bmi", "smoking",
        "hba1c", "pir", "albumin_g_dl",
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
        {"concept": "absolute_vitamin_e", "fields": "LBDVIESI/LBDATCSI", "meaning": "serum alpha-tocopherol", "unit": "umol/L", "quality_rule": "2-100; conventional low <12"},
        {"concept": "lipid_adjusted_vitamin_e", "fields": "alpha-tocopherol divided by total cholesterol", "meaning": "alpha-tocopherol-to-total-cholesterol ratio", "unit": "umol/mmol", "quality_rule": "low <2.2"},
        {"concept": "peripheral_neuropathy", "fields": "LEALPN+LEARPN", "meaning": "objective 10-g monofilament insensate plantar sites", "unit": "0-6", "quality_rule": "both feet valid; primary >=1; sensitivity >=2"},
        {"concept": "confounder_exclusion", "fields": "LBXB12/LBXMMA/creatinine/DIQ010/LBXGH", "meaning": "B12 status, kidney function and diabetes", "unit": "mixed", "quality_rule": "B12 >=200; MMA <=0.40; eGFR >=60; diabetes excluded"},
        {"concept": "survey_design", "fields": "WTMEC2YR/WTMEC4YR, SDMVPSU, SDMVSTRA", "meaning": "MEC weights and design", "unit": "survey", "quality_rule": "pooled weights divided by 3; cycle-unique PSU/strata"},
    ])

    columns = [
        "SEQN", "cycle", "cycle_years", "period", "age", "sex", "race", "education", "pir", "pregnancy",
        "weight", "psu_u", "strata_u", "alpha_tocopherol_umol_l", "total_cholesterol_mg_dl",
        "total_cholesterol_mmol_l", "hdl_mg_dl", "vit_e_chol_ratio", "log_ratio", "log_alpha_tocopherol",
        "b12_pg_ml", "mma_umol_l", "log_mma", "folate_nmol_l", "folate_ng_ml", "hba1c",
        "creatinine_mg_dl", "egfr", "albumin_g_dl", "total_protein_g_dl", "bmi", "body_weight_kg",
        "height_cm", "smoking", "diabetes", "left_insensate", "right_insensate", "pn_sites", "pn_any",
        "pn_two_plus", "group", "complete_case",
    ]
    cohort.loc[:, columns].to_csv(OUT / "c_n10_analysis_core.csv", index=False)
    pd.DataFrame(audits).to_csv(OUT / "c_n10_source_audit.csv", index=False)
    pd.DataFrame(flow).to_csv(OUT / "c_n10_flow.csv", index=False)
    pd.DataFrame(group_rows).to_csv(OUT / "c_n10_group_counts.csv", index=False)
    pd.DataFrame(missing_rows).to_csv(OUT / "c_n10_missingness.csv", index=False)
    semantic.to_csv(OUT / "c_n10_semantic_audit.csv", index=False)

    key = cohort.loc[cohort["group"].eq("ratio_low_absolute_normal")]
    val_key = key.loc[key["period"].eq("validation")]
    status = {
        "candidate_code": "C-N10",
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
    (OUT / "c_n10_prep_status.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
    print(json.dumps(status, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
