from __future__ import annotations

import io
import json
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import requests

BASE = "https://raw.githubusercontent.com/protobi/nhanes-continuous/main"
OUT = Path("nutrition_stage/output_batch_n11")
OUT.mkdir(parents=True, exist_ok=True)

CYCLES = {
    "D": ("2005-2006", 1),
    "E": ("2007-2008", 2),
    "F": ("2009-2010", 3),
    "G": ("2011-2012", 4),
    "H": ("2013-2014", 5),
    "I": ("2015-2016", 6),
    "J": ("2017-2018", 7),
}

THIAZIDE_LOOP = re.compile(
    r"hydrochlorothiazide|chlorthalidone|indapamide|metolazone|furosemide|bumetanide|torsemide|ethacrynic|chlorothiazide|bendroflumethiazide|trichlormethiazide|methyclothiazide|polythiazide",
    re.I,
)
THIAZIDE = re.compile(
    r"hydrochlorothiazide|chlorthalidone|indapamide|metolazone|chlorothiazide|bendroflumethiazide|trichlormethiazide|methyclothiazide|polythiazide",
    re.I,
)
LOOP = re.compile(r"furosemide|bumetanide|torsemide|ethacrynic", re.I)
K_SPARING = re.compile(r"spironolactone|eplerenone|amiloride|triamterene", re.I)
POTASSIUM_RX = re.compile(r"potassium chloride|potassium citrate|potassium bicarbonate|potassium gluconate", re.I)
ACE_ARB = re.compile(
    r"lisinopril|enalapril|ramipril|benazepril|captopril|quinapril|fosinopril|moexipril|perindopril|trandolapril|losartan|valsartan|irbesartan|candesartan|telmisartan|olmesartan|eprosartan|azilsartan",
    re.I,
)


def fetch_csv(paths: list[str], required: bool = True) -> tuple[pd.DataFrame, str]:
    errors: list[str] = []
    for path in paths:
        url = f"{BASE}/{path}"
        r = requests.get(url, timeout=120)
        if r.status_code != 200:
            errors.append(f"{path}:HTTP{r.status_code}")
            continue
        try:
            df = pd.read_csv(io.BytesIO(r.content), low_memory=False)
            return df, path
        except Exception as exc:  # noqa: BLE001
            errors.append(f"{path}:{type(exc).__name__}")
    if required:
        raise RuntimeError("No usable source: " + " | ".join(errors))
    return pd.DataFrame(columns=["SEQN"]), ""


def keep(df: pd.DataFrame, cols: Iterable[str]) -> pd.DataFrame:
    use = [c for c in cols if c in df.columns]
    if "SEQN" in df.columns and "SEQN" not in use:
        use.insert(0, "SEQN")
    return df.loc[:, use].copy()


def first_num(df: pd.DataFrame, names: Iterable[str]) -> pd.Series:
    out = pd.Series(np.nan, index=df.index, dtype=float)
    for name in names:
        if name in df.columns:
            out = out.fillna(pd.to_numeric(df[name], errors="coerce"))
    return out


def egfr_2021(scr: pd.Series, age: pd.Series, sex: pd.Series) -> pd.Series:
    female = sex.eq(2)
    k = np.where(female, 0.7, 0.9)
    alpha = np.where(female, -0.241, -0.302)
    ratio = scr / k
    return 142 * np.minimum(ratio, 1) ** alpha * np.maximum(ratio, 1) ** -1.2 * 0.9938 ** age * np.where(female, 1.012, 1)


def med_flags(rx: pd.DataFrame) -> pd.DataFrame:
    if rx.empty or "SEQN" not in rx.columns:
        return pd.DataFrame(columns=["SEQN", "diuretic", "thiazide", "loop", "k_sparing", "potassium_rx", "ace_arb"])
    text_cols = [c for c in ["RXDDRUG", "RXDDRGID", "RXDINGFL", "RXDRSC1"] if c in rx.columns]
    if not text_cols:
        text = pd.Series("", index=rx.index)
    else:
        text = rx[text_cols].fillna("").astype(str).agg(" ".join, axis=1)
    tmp = pd.DataFrame({"SEQN": rx["SEQN"]})
    tmp["diuretic"] = text.str.contains(THIAZIDE_LOOP)
    tmp["thiazide"] = text.str.contains(THIAZIDE)
    tmp["loop"] = text.str.contains(LOOP)
    tmp["k_sparing"] = text.str.contains(K_SPARING)
    tmp["potassium_rx"] = text.str.contains(POTASSIUM_RX)
    tmp["ace_arb"] = text.str.contains(ACE_ARB)
    return tmp.groupby("SEQN", as_index=False).max()


def merge_cycle(cycle: str) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    years, idx = CYCLES[cycle]
    audit: list[dict[str, object]] = []

    def load(paths: list[str], required: bool = True) -> pd.DataFrame:
        df, used = fetch_csv(paths, required=required)
        audit.append({"cycle": cycle, "years": years, "requested": "|".join(paths), "used": used, "required": required, "rows": len(df), "columns": df.shape[1]})
        return df

    demo = load([f"Demographics/DEMO_{cycle}.csv"])
    bmx = load([f"Examination/BMX_{cycle}.csv"])
    d1 = load([f"Dietary/DR1TOT_{cycle}.csv"])
    d2 = load([f"Dietary/DR2TOT_{cycle}.csv"], required=False)
    bio = load([f"Laboratory/BIOPRO_{cycle}.csv", f"Laboratory/L40_{cycle}.csv", f"Laboratory/LAB18_{cycle}.csv"])
    rx = load([f"Questionnaire/RXQ_RX_{cycle}.csv"])
    smq = load([f"Questionnaire/SMQ_{cycle}.csv"], required=False)
    diq = load([f"Questionnaire/DIQ_{cycle}.csv"], required=False)
    bpq = load([f"Questionnaire/BPQ_{cycle}.csv"], required=False)
    mcq = load([f"Questionnaire/MCQ_{cycle}.csv"], required=False)

    m = keep(demo, ["SEQN", "RIAGENDR", "RIDAGEYR", "RIDRETH1", "RIDRETH3", "RIDEXPRG", "WTMEC2YR", "SDMVPSU", "SDMVSTRA", "INDFMPIR"])
    pieces = [
        keep(bmx, ["SEQN", "BMXBMI"]),
        keep(d1, ["SEQN", "WTDRD1", "DR1TKCAL", "DR1TPOTA", "DR1TSODI", "DR1TMAGN", "DR1TPROT"]),
        keep(d2, ["SEQN", "WTDR2D", "DR2TKCAL", "DR2TPOTA", "DR2TSODI", "DR2TMAGN", "DR2TPROT"]),
        keep(bio, ["SEQN", "LBXSKSI", "LBXSNA", "LBXSCR", "LBXSAL", "LBXSMGSI", "LBXSMSI"]),
        med_flags(rx),
        keep(smq, ["SEQN", "SMQ020", "SMQ040"]),
        keep(diq, ["SEQN", "DIQ010"]),
        keep(bpq, ["SEQN", "BPQ020", "BPQ040A"]),
        keep(mcq, ["SEQN", "MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E", "MCQ160F"]),
    ]
    for piece in pieces:
        if not piece.empty:
            m = m.merge(piece.drop_duplicates("SEQN"), on="SEQN", how="left")

    m["cycle"] = cycle
    m["cycle_index"] = idx
    m["period"] = "discovery" if cycle in ["D", "E", "F", "G"] else "validation"
    m["age"] = first_num(m, ["RIDAGEYR"])
    m["sex"] = first_num(m, ["RIAGENDR"])
    m["race"] = first_num(m, ["RIDRETH1", "RIDRETH3"])
    m["pregnancy"] = first_num(m, ["RIDEXPRG"])
    m["bmi"] = first_num(m, ["BMXBMI"])
    m["pir"] = first_num(m, ["INDFMPIR"])
    m["weight"] = first_num(m, ["WTDRD1"]) / len(CYCLES)
    m["psu"] = first_num(m, ["SDMVPSU"])
    m["strata"] = first_num(m, ["SDMVSTRA"])
    m["psu_u"] = idx * 100 + m["psu"]
    m["strata_u"] = idx * 1000 + m["strata"]
    m["k_serum"] = first_num(m, ["LBXSKSI"])
    m["na_serum"] = first_num(m, ["LBXSNA"])
    m["scr"] = first_num(m, ["LBXSCR"])
    m["albumin"] = first_num(m, ["LBXSAL"])
    m["energy1"] = first_num(m, ["DR1TKCAL"])
    m["potassium1"] = first_num(m, ["DR1TPOTA"])
    m["sodium1"] = first_num(m, ["DR1TSODI"])
    m["magnesium1"] = first_num(m, ["DR1TMAGN"])
    m["protein1"] = first_num(m, ["DR1TPROT"])
    m["energy2"] = first_num(m, ["DR2TKCAL"])
    m["potassium2"] = first_num(m, ["DR2TPOTA"])
    m["sodium2"] = first_num(m, ["DR2TSODI"])
    m["potassium_mean2"] = m[["potassium1", "potassium2"]].mean(axis=1, skipna=False)
    m["sodium_mean2"] = m[["sodium1", "sodium2"]].mean(axis=1, skipna=False)
    m["egfr"] = egfr_2021(m["scr"], m["age"], m["sex"])
    m["smoking"] = np.select(
        [first_num(m, ["SMQ020"]).eq(2), first_num(m, ["SMQ020"]).eq(1) & first_num(m, ["SMQ040"]).isin([1, 2]), first_num(m, ["SMQ020"]).eq(1) & first_num(m, ["SMQ040"]).eq(3)],
        ["never", "current", "former"], default=np.nan,
    )
    m["diabetes"] = np.where(first_num(m, ["DIQ010"]).eq(1), 1, np.where(first_num(m, ["DIQ010"]).isin([2, 3]), 0, np.nan))
    m["hypertension"] = np.where(first_num(m, ["BPQ020"]).eq(1), 1, np.where(first_num(m, ["BPQ020"]).eq(2), 0, np.nan))
    cv = pd.concat([first_num(m, [x]).eq(1) for x in ["MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E", "MCQ160F"]], axis=1)
    m["cvd"] = cv.any(axis=1).astype(int)
    for flag in ["diuretic", "thiazide", "loop", "k_sparing", "potassium_rx", "ace_arb"]:
        if flag not in m:
            m[flag] = False
        m[flag] = m[flag].fillna(False).astype(bool)
    return m, audit


def main() -> None:
    frames, audits = [], []
    for cycle in CYCLES:
        frame, audit = merge_cycle(cycle)
        frames.append(frame)
        audits.extend(audit)
    raw = pd.concat(frames, ignore_index=True, sort=False)
    flow = []
    def record(step: str, df: pd.DataFrame) -> None:
        flow.append({"step": step, "n": int(len(df))})
    record("all_selected_cycles", raw)
    c = raw.loc[raw["age"].ge(40)].copy(); record("age_40_plus", c)
    c = c.loc[c["pregnancy"].ne(1) | c["pregnancy"].isna()].copy(); record("exclude_known_pregnancy", c)
    c = c.loc[c["energy1"].between(500, 5000) & c["potassium1"].between(100, 10000) & c["sodium1"].between(100, 15000)].copy(); record("valid_day1_diet", c)
    c = c.loc[c["k_serum"].between(2.0, 7.0) & c["na_serum"].between(115, 160) & c["scr"].between(0.2, 15)].copy(); record("valid_serum_electrolytes_and_creatinine", c)
    c = c.loc[c["weight"].gt(0) & c["psu_u"].notna() & c["strata_u"].notna()].copy(); record("valid_survey_design", c)
    c["hypokalemia"] = (c["k_serum"] < 3.5).astype(int)
    c["hyponatremia"] = (c["na_serum"] < 135).astype(int)
    c["hyperkalemia"] = (c["k_serum"] > 5.0).astype(int)
    c["low_k_intake"] = (c["potassium1"] < 2000).astype(int)
    c["low_na_intake"] = (c["sodium1"] < 2000).astype(int)
    c["high_k_intake"] = (c["potassium1"] >= 3500).astype(int)
    c["k_per_1000kcal"] = c["potassium1"] / c["energy1"] * 1000
    c["na_per_1000kcal"] = c["sodium1"] / c["energy1"] * 1000
    full_cov = ["age", "sex", "race", "bmi", "egfr", "diabetes", "hypertension", "cvd", "smoking", "energy1", "sodium1"]
    c["complete_case_k"] = c[full_cov].notna().all(axis=1)
    full_cov_na = ["age", "sex", "race", "bmi", "egfr", "diabetes", "hypertension", "cvd", "smoking", "energy1", "potassium1"]
    c["complete_case_na"] = c[full_cov_na].notna().all(axis=1)
    record("complete_case_k_model", c.loc[c["complete_case_k"]])
    record("complete_case_na_model", c.loc[c["complete_case_na"]])

    cols = ["SEQN", "cycle", "period", "age", "sex", "race", "pregnancy", "weight", "psu_u", "strata_u", "pir", "bmi", "smoking", "diabetes", "hypertension", "cvd", "egfr", "albumin", "energy1", "potassium1", "sodium1", "magnesium1", "protein1", "potassium_mean2", "sodium_mean2", "k_serum", "na_serum", "diuretic", "thiazide", "loop", "k_sparing", "potassium_rx", "ace_arb", "hypokalemia", "hyponatremia", "hyperkalemia", "low_k_intake", "low_na_intake", "high_k_intake", "k_per_1000kcal", "na_per_1000kcal", "complete_case_k", "complete_case_na"]
    c[cols].to_csv(OUT / "batch_n11_core.csv", index=False)
    pd.DataFrame(audits).to_csv(OUT / "source_audit.csv", index=False)
    pd.DataFrame(flow).to_csv(OUT / "flow.csv", index=False)

    miss = []
    for scope, df in [("overall", c), ("discovery", c[c.period == "discovery"]), ("validation", c[c.period == "validation"])]:
        for var in full_cov + ["potassium1", "sodium1", "k_serum", "na_serum"]:
            miss.append({"scope": scope, "variable": var, "n": len(df), "missing_n": int(df[var].isna().sum()), "missing_pct": float(df[var].isna().mean()*100)})
    pd.DataFrame(miss).to_csv(OUT / "missingness.csv", index=False)

    groups = []
    candidate_specs = [
        ("C-N11", c[c.diuretic & ~c.k_sparing], "low_k_intake", "hypokalemia"),
        ("C-N12", c[c.thiazide & ~c.k_sparing], "low_na_intake", "hyponatremia"),
        ("C-N13", c[c.ace_arb & ~c.diuretic], "high_k_intake", "hyperkalemia"),
    ]
    for code, df, exposure, outcome in candidate_specs:
        for scope, z in [("overall", df), ("discovery", df[df.period == "discovery"]), ("validation", df[df.period == "validation"])]:
            for x, g in z.groupby(exposure):
                groups.append({"candidate_code": code, "scope": scope, "exposed": int(x), "n": len(g), "events": int(g[outcome].sum()), "event_pct": float(g[outcome].mean()*100) if len(g) else np.nan, "mean_k": float(g.k_serum.mean()), "mean_na": float(g.na_serum.mean())})
    pd.DataFrame(groups).to_csv(OUT / "group_counts.csv", index=False)

    semantic = pd.DataFrame([
        {"concept": "dietary_potassium", "field": "DR1TPOTA", "unit": "mg/day", "semantics": "day-1 24-hour dietary recall total", "main_rule": "low <2000 mg/day"},
        {"concept": "dietary_sodium", "field": "DR1TSODI", "unit": "mg/day", "semantics": "day-1 24-hour dietary recall total", "main_rule": "low <2000 mg/day"},
        {"concept": "serum_potassium", "field": "LBXSKSI", "unit": "mmol/L", "semantics": "measured serum potassium", "main_rule": "hypokalemia <3.5; hyperkalemia >5.0"},
        {"concept": "serum_sodium", "field": "LBXSNA", "unit": "mmol/L", "semantics": "measured serum sodium", "main_rule": "hyponatremia <135"},
        {"concept": "medication_exposure", "field": "RXDDRUG plus ingredient strings", "unit": "current prescription drug", "semantics": "participant-level any current named ingredient", "main_rule": "frozen regex classes; potassium-sparing co-use excluded"},
        {"concept": "survey_design", "field": "WTDRD1/SDMVPSU/SDMVSTRA", "unit": "survey", "semantics": "day-1 dietary weight and NHANES design", "main_rule": "weight divided by seven cycles; cycle-unique PSU/strata"},
    ])
    semantic.to_csv(OUT / "semantic_audit.csv", index=False)
    print(json.dumps({"rows": len(c), "cycles": sorted(c.cycle.unique()), "candidate_codes": [x[0] for x in candidate_specs]}, indent=2))

if __name__ == "__main__":
    main()
