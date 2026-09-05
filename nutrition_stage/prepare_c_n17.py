from __future__ import annotations

import io
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import requests

BASE = "https://raw.githubusercontent.com/protobi/nhanes-continuous/main"
TREE_URL = "https://api.github.com/repos/protobi/nhanes-continuous/git/trees/main?recursive=1"
OUT = Path("nutrition_stage/output_c_n17")
OUT.mkdir(parents=True, exist_ok=True)
CYCLES = {"C": ("2003-2004", 1, "discovery"), "D": ("2005-2006", 2, "validation")}


def get_tree() -> list[str]:
    r = requests.get(TREE_URL, timeout=120)
    r.raise_for_status()
    return [x["path"] for x in r.json()["tree"] if x.get("type") == "blob"]


def get_csv(path: str) -> pd.DataFrame:
    r = requests.get(f"{BASE}/{path}", timeout=120)
    r.raise_for_status()
    return pd.read_csv(io.BytesIO(r.content), low_memory=False)


def load_exact(tree: list[str], paths: list[str], required: bool = True) -> tuple[pd.DataFrame, str]:
    for path in paths:
        if path not in tree:
            continue
        try:
            return get_csv(path), path
        except Exception:
            pass
    if required:
        raise RuntimeError(f"No usable exact source: {paths}")
    return pd.DataFrame(columns=["SEQN"]), ""


def find_lab(tree: list[str], cycle: str, target: str) -> tuple[pd.DataFrame, str, str]:
    if target == "pth":
        wanted = {"LBXPTH", "LBDPTHSI", "LBXPTHSI"}
        name_re = re.compile(r"PTH|PARATH|L06", re.I)
    else:
        wanted = {"LBXVIDMS", "LBDVIDMS", "LBXVID", "LBDVID", "LBXVD2MS", "LBXVD3MS"}
        name_re = re.compile(r"VID|VITD|DXX|25OH", re.I)
    exact = [
        f"Laboratory/PTH_{cycle}.csv",
        f"Laboratory/VID_{cycle}.csv",
        f"Laboratory/VIDMS_{cycle}.csv",
        f"Laboratory/L06_{cycle}.csv",
        f"Laboratory/L06_2_{cycle}.csv",
        f"Laboratory/L06PTH_{cycle}.csv",
        f"Laboratory/L06VID_{cycle}.csv",
    ]
    candidates = exact + [p for p in tree if p.startswith("Laboratory/") and p.endswith(f"_{cycle}.csv") and name_re.search(p)]
    seen: set[str] = set()
    for path in candidates:
        if path in seen or path not in tree:
            continue
        seen.add(path)
        try:
            df = get_csv(path)
        except Exception:
            continue
        matches = [c for c in df.columns if c.upper() in wanted]
        if not matches:
            if target == "pth":
                matches = [c for c in df.columns if "PTH" in c.upper() and c.upper().startswith(("LBX", "LBD"))]
            else:
                matches = [c for c in df.columns if re.search(r"VID|VD[23]|25OH", c.upper()) and c.upper().startswith(("LBX", "LBD"))]
        for c in matches:
            x = pd.to_numeric(df[c], errors="coerce")
            med = x.median()
            if target == "pth" and pd.notna(med) and 5 < med < 300:
                return df, path, c
            if target == "vitd" and pd.notna(med) and 5 < med < 250:
                return df, path, c
    raise RuntimeError(f"No {target} field found for cycle {cycle}")


def keep(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    use = [c for c in cols if c in df.columns]
    if "SEQN" in df.columns and "SEQN" not in use:
        use.insert(0, "SEQN")
    return df[use].copy()


def num(df: pd.DataFrame, names: list[str]) -> pd.Series:
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


def merge_cycle(tree: list[str], cycle: str) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    years, idx, period = CYCLES[cycle]
    audit: list[dict[str, object]] = []

    def load(paths: list[str], required: bool = True) -> pd.DataFrame:
        df, used = load_exact(tree, paths, required=required)
        audit.append({"cycle": cycle, "years": years, "used": used, "required": required, "rows": len(df), "columns": df.shape[1]})
        return df

    demo = load([f"Demographics/DEMO_{cycle}.csv"])
    bmx = load([f"Examination/BMX_{cycle}.csv"])
    d1 = load([f"Dietary/DR1TOT_{cycle}.csv"])
    d2 = load([f"Dietary/DR2TOT_{cycle}.csv"], required=False)
    bio = load([f"Laboratory/BIOPRO_{cycle}.csv", f"Laboratory/L40_{cycle}.csv"])
    smq = load([f"Questionnaire/SMQ_{cycle}.csv"], required=False)
    diq = load([f"Questionnaire/DIQ_{cycle}.csv"], required=False)
    pth, pth_path, pth_col = find_lab(tree, cycle, "pth")
    vitd, vitd_path, vitd_col = find_lab(tree, cycle, "vitd")
    audit += [
        {"cycle": cycle, "years": years, "used": pth_path, "rows": len(pth), "columns": pth.shape[1], "selected_field": pth_col, "concept": "PTH"},
        {"cycle": cycle, "years": years, "used": vitd_path, "rows": len(vitd), "columns": vitd.shape[1], "selected_field": vitd_col, "concept": "25OHD"},
    ]

    m = keep(demo, ["SEQN", "RIAGENDR", "RIDAGEYR", "RIDRETH1", "RIDRETH3", "WTMEC2YR", "SDMVPSU", "SDMVSTRA", "INDFMPIR"])
    pieces = [
        keep(bmx, ["SEQN", "BMXBMI"]),
        keep(d1, ["SEQN", "WTDRD1", "DR1TKCAL", "DR1TCALC", "DR1TPHOS", "DR1TPROT"]),
        keep(d2, ["SEQN", "DR2TKCAL", "DR2TCALC", "DR2TPHOS"]),
        keep(bio, ["SEQN", "LBXSCR", "LBXSCA", "LBXSAL", "LBXSPH"]),
        keep(smq, ["SEQN", "SMQ020", "SMQ040"]),
        keep(diq, ["SEQN", "DIQ010"]),
        keep(pth, ["SEQN", pth_col]),
        keep(vitd, ["SEQN", vitd_col, "LBXVD2MS", "LBXVD3MS"]),
    ]
    for piece in pieces:
        m = m.merge(piece.drop_duplicates("SEQN"), on="SEQN", how="left")

    m["cycle"] = cycle
    m["period"] = period
    m["age"] = num(m, ["RIDAGEYR"])
    m["sex"] = num(m, ["RIAGENDR"])
    m["race"] = num(m, ["RIDRETH1", "RIDRETH3"])
    m["bmi"] = num(m, ["BMXBMI"])
    m["pir"] = num(m, ["INDFMPIR"])
    m["weight"] = num(m, ["WTDRD1"]) / 2
    m["psu_u"] = idx * 100 + num(m, ["SDMVPSU"])
    m["strata_u"] = idx * 1000 + num(m, ["SDMVSTRA"])
    m["energy1"] = num(m, ["DR1TKCAL"])
    m["calcium1"] = num(m, ["DR1TCALC"])
    m["phosphorus1"] = num(m, ["DR1TPHOS"])
    m["calcium2"] = num(m, ["DR2TCALC"])
    m["calcium_mean2"] = pd.concat([m["calcium1"], m["calcium2"]], axis=1).mean(axis=1, skipna=False)
    m["scr"] = num(m, ["LBXSCR"])
    m["serum_calcium"] = num(m, ["LBXSCA"])
    m["albumin"] = num(m, ["LBXSAL"])
    m["serum_phosphorus"] = num(m, ["LBXSPH"])
    m["egfr"] = egfr_2021(m["scr"], m["age"], m["sex"])
    m["pth_pg_ml"] = num(m, [pth_col])
    vit = num(m, [vitd_col])
    med = vit.median()
    m["vitd_ng_ml"] = vit / 2.496 if pd.notna(med) and med > 80 else vit
    m["smoking"] = np.select(
        [num(m, ["SMQ020"]).eq(2), num(m, ["SMQ020"]).eq(1) & num(m, ["SMQ040"]).eq(3), num(m, ["SMQ020"]).eq(1) & num(m, ["SMQ040"]).isin([1, 2])],
        ["never", "former", "current"], default=np.nan,
    )
    m["diabetes"] = np.where(num(m, ["DIQ010"]).eq(1), 1, np.where(num(m, ["DIQ010"]).isin([2, 3]), 0, np.nan))
    return m, audit


def main() -> None:
    tree = get_tree()
    frames, audits = [], []
    for cycle in CYCLES:
        x, a = merge_cycle(tree, cycle)
        frames.append(x)
        audits.extend(a)
    raw = pd.concat(frames, ignore_index=True, sort=False)
    flow: list[dict[str, object]] = []
    def rec(step: str, frame: pd.DataFrame) -> None:
        flow.append({"step": step, "n": len(frame)})
    rec("all_linked", raw)
    c = raw.loc[raw["age"].between(20, 79)].copy(); rec("age_20_79", c)
    c = c.loc[c["vitd_ng_ml"].between(3, 100) & c["pth_pg_ml"].between(5, 500)].copy(); rec("valid_25ohd_and_pth", c)
    c = c.loc[c["energy1"].between(500, 5000) & c["calcium1"].between(50, 4000)].copy(); rec("valid_day1_diet", c)
    c = c.loc[c["egfr"].ge(60) & c["serum_calcium"].between(8.4, 10.2)].copy(); rec("exclude_ckd_and_abnormal_serum_calcium", c)
    c = c.loc[c["vitd_ng_ml"].lt(20)].copy(); rec("vitamin_d_insufficient_primary_population", c)
    c = c.loc[c["weight"].gt(0) & c["psu_u"].notna() & c["strata_u"].notna()].copy(); rec("valid_survey_design", c)
    c["low_calcium"] = (c["calcium1"] < 800).astype(int)
    c["elevated_pth"] = (c["pth_pg_ml"] > 65).astype(int)
    cov = ["age", "sex", "race", "bmi", "smoking", "diabetes", "energy1", "phosphorus1", "vitd_ng_ml", "egfr", "serum_calcium", "albumin"]
    c["complete_case"] = c[cov].notna().all(axis=1); rec("complete_case_main", c.loc[c["complete_case"]])
    outcols = ["SEQN", "cycle", "period", "age", "sex", "race", "weight", "psu_u", "strata_u", "pir", "bmi", "smoking", "diabetes", "energy1", "calcium1", "calcium_mean2", "phosphorus1", "vitd_ng_ml", "pth_pg_ml", "egfr", "serum_calcium", "albumin", "serum_phosphorus", "low_calcium", "elevated_pth", "complete_case"]
    c[outcols].to_csv(OUT / "core.csv", index=False)
    pd.DataFrame(flow).to_csv(OUT / "flow.csv", index=False)
    pd.DataFrame(audits).to_csv(OUT / "source_audit.csv", index=False)
    miss = []
    for scope, x in [("overall", c), ("discovery", c[c.period == "discovery"]), ("validation", c[c.period == "validation"])]:
        for v in cov + ["pir", "calcium_mean2", "serum_phosphorus"]:
            miss.append({"scope": scope, "variable": v, "n": len(x), "missing_n": int(x[v].isna().sum()), "missing_pct": float(x[v].isna().mean() * 100)})
    pd.DataFrame(miss).to_csv(OUT / "missingness.csv", index=False)
    groups = []
    for scope, x in [("overall", c), ("discovery", c[c.period == "discovery"]), ("validation", c[c.period == "validation"])]:
        for z, g in x.groupby("low_calcium"):
            groups.append({"scope": scope, "low_calcium": int(z), "n": len(g), "events": int(g.elevated_pth.sum()), "event_pct": float(g.elevated_pth.mean() * 100), "mean_pth": float(g.pth_pg_ml.mean())})
    pd.DataFrame(groups).to_csv(OUT / "group_counts.csv", index=False)
    pd.DataFrame([
        {"concept": "vitamin_D_insufficiency", "field": "cycle-specific 25OHD field", "unit": "ng/mL", "semantics": "measured total 25-hydroxyvitamin D", "rule": "primary population <20; plausible 3-100; nmol/L converted if required"},
        {"concept": "dietary_calcium", "field": "DR1TCALC", "unit": "mg/day", "semantics": "day-1 24-hour dietary recall", "rule": "low <800 mg/day"},
        {"concept": "PTH", "field": "cycle-specific LBXPTH-equivalent", "unit": "pg/mL", "semantics": "measured intact parathyroid hormone", "rule": "elevated >65; plausible 5-500"},
        {"concept": "survey_design", "field": "WTDRD1/SDMVPSU/SDMVSTRA", "unit": "survey", "semantics": "day-1 dietary design", "rule": "pooled weight divided by two; cycle-unique PSU/strata"},
    ]).to_csv(OUT / "semantic_audit.csv", index=False)
    prep = {"candidate_code": "C-N17", "actual_n": len(c), "complete_case_n": int(c.complete_case.sum()), "retention": float(c.complete_case.mean()), "events": int(c.elevated_pth.sum()), "key_group_n": int(c.low_calcium.sum()), "key_group_events": int(c.loc[c.low_calcium.eq(1), "elevated_pth"].sum()), "validation_key_group_n": int(c.loc[(c.period == "validation") & c.low_calcium.eq(1)].shape[0]), "validation_key_group_events": int(c.loc[(c.period == "validation") & c.low_calcium.eq(1), "elevated_pth"].sum())}
    (OUT / "prep_status.json").write_text(json.dumps(prep, indent=2))
    print(json.dumps(prep, indent=2))

if __name__ == "__main__":
    main()
