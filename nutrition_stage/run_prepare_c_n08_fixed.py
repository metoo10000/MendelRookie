from __future__ import annotations

import re
import traceback
from pathlib import Path

source_path = Path(__file__).with_name("prepare_c_n08.py")
source = source_path.read_text(encoding="utf-8")

mortality_function = r'''
def load_mortality(cycle: str) -> tuple[pd.DataFrame, dict[str, object]]:
    url = mortality_url(cycle)
    response = requests.get(url, timeout=180)
    if response.status_code != 200:
        raise RuntimeError(f"Unable to load mortality {cycle}: HTTP {response.status_code}")
    names = ["SEQN", "ELIGSTAT", "MORTSTAT", "UCOD_LEADING", "DIABETES_MORT", "HYPERTEN_MORT", "PERMTH_INT", "PERMTH_EXM"]
    colspecs = [(0, 6), (14, 15), (15, 16), (16, 19), (19, 20), (20, 21), (42, 45), (45, 48)]
    frame = pd.read_fwf(io.StringIO(response.text), colspecs=colspecs, names=names, dtype=str)
    for column in ["SEQN", "ELIGSTAT", "MORTSTAT", "DIABETES_MORT", "HYPERTEN_MORT", "PERMTH_INT", "PERMTH_EXM"]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame = frame.loc[frame["SEQN"].notna()].copy()
    frame["SEQN"] = frame["SEQN"].astype(int)
    audit = {
        "cycle": cycle,
        "component": "mortality",
        "file": url.rsplit("/", 1)[-1],
        "url": url,
        "rows": int(frame.shape[0]),
        "columns": int(frame.shape[1]),
        "column_names": "|".join(frame.columns),
    }
    if frame.empty:
        raise RuntimeError(f"Mortality parser produced zero rows for cycle {cycle}")
    return frame, audit
'''

patched, n = re.subn(
    r"def load_mortality\(cycle: str\).*?(?=\ndef keep\()",
    mortality_function + "\n",
    source,
    count=1,
    flags=re.S,
)
if n != 1:
    raise RuntimeError("Could not replace the frozen mortality parser")

old_merge = '''    merged = pieces[0]
    for piece in pieces[1:]:
        merged = merged.merge(piece.drop_duplicates("SEQN"), on="SEQN", how="left")
'''
new_merge = '''    merged = pieces[0]
    piece_names = ["demo", "day1", "day2", "body", "activity", "smoking", "diabetes", "hba1c", "biochem", "medical", "mortality"]
    for piece_name, piece in zip(piece_names[1:], pieces[1:]):
        if "SEQN" not in piece.columns:
            raise RuntimeError(f"Merge component {piece_name} lacks SEQN; columns={list(piece.columns)}; shape={piece.shape}")
        merged = merged.merge(piece.drop_duplicates("SEQN"), on="SEQN", how="left")
'''
if old_merge not in patched:
    raise RuntimeError("Expected frozen merge block was not found")
patched = patched.replace(old_merge, new_merge, 1)

namespace = {"__name__": "__main__", "__file__": str(source_path)}
try:
    exec(compile(patched, str(source_path), "exec"), namespace, namespace)
except Exception:
    out = Path("nutrition_stage/output_c_n08")
    out.mkdir(parents=True, exist_ok=True)
    text = traceback.format_exc()
    (out / "prepare_error.txt").write_text(text, encoding="utf-8")
    print("BEGIN_C_N08_PREPARE_ERROR")
    print(text)
    print("END_C_N08_PREPARE_ERROR")
    raise
