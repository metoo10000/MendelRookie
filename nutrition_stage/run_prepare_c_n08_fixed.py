from __future__ import annotations

import traceback
from pathlib import Path

source_path = Path(__file__).with_name("prepare_c_n08.py")
source = source_path.read_text(encoding="utf-8")
old = '"SEQN": int(line[0:14].strip()),'
new = '"SEQN": int(line[0:6].strip()),'
if old not in source:
    raise RuntimeError("Expected frozen mortality parser line was not found")
patched = source.replace(old, new, 1)
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
