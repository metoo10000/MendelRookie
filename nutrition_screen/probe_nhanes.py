from __future__ import annotations

import io
import json
import sys
from pathlib import Path

import pandas as pd
import requests

CYCLES = {
    "2005-2006": "D",
    "2007-2008": "E",
    "2009-2010": "F",
    "2011-2012": "G",
    "2013-2014": "H",
    "2015-2016": "I",
    "2017-2018": "J",
}

MODULES = [
    "DEMO",
    "DR1TOT",
    "DR2TOT",
    "BIOPRO",
    "RXQ_RX",
    "CBC",
    "MCQ",
    "BPQ",
    "DIQ",
    "KIQ_U",
    "ALQ",
    "SMQ",
]

OUT = Path("nutrition_screen/output")
OUT.mkdir(parents=True, exist_ok=True)


def read_xpt(url: str) -> pd.DataFrame:
    response = requests.get(url, timeout=90)
    response.raise_for_status()
    return pd.read_sas(io.BytesIO(response.content), format="xport", encoding="latin1")


def main() -> int:
    results: list[dict[str, object]] = []
    for cycle, suffix in CYCLES.items():
        for module in MODULES:
            filename = f"{module}_{suffix}.XPT"
            url = f"https://wwwn.cdc.gov/Nchs/Nhanes/{cycle}/{filename}"
            row: dict[str, object] = {
                "cycle": cycle,
                "suffix": suffix,
                "module": module,
                "filename": filename,
                "url": url,
            }
            try:
                df = read_xpt(url)
                row.update(
                    ok=True,
                    rows=int(df.shape[0]),
                    columns=int(df.shape[1]),
                    column_names=list(map(str, df.columns)),
                )
            except Exception as exc:  # noqa: BLE001
                row.update(ok=False, error=f"{type(exc).__name__}: {exc}")
            results.append(row)
            print(json.dumps(row, ensure_ascii=False), flush=True)

    output_path = OUT / "nhanes_module_probe.json"
    output_path.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")

    ok_rows = [r for r in results if r.get("ok")]
    print(f"PROBE_COMPLETE modules_ok={len(ok_rows)} modules_total={len(results)} output={output_path}")
    return 0 if ok_rows else 1


if __name__ == "__main__":
    sys.exit(main())
