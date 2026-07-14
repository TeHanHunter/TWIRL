#!/usr/bin/env python3
"""Verify the blinded Franklin current-A2v1 re-review handoff."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queue", type=Path, required=True)
    parser.add_argument("--private-manifest", type=Path, required=True)
    parser.add_argument("--render-summary", type=Path, required=True)
    parser.add_argument("--sheet-dir", type=Path, required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    args = parser.parse_args()

    queue = pd.read_csv(args.queue, low_memory=False)
    private = pd.read_parquet(args.private_manifest)
    render = json.loads(args.render_summary.read_text())
    failures: list[str] = []

    def require(condition: bool, message: str) -> None:
        if not condition:
            failures.append(message)

    require(len(queue) == 134, f"public queue has {len(queue)} rows, expected 134")
    require(len(private) == 134, f"private manifest has {len(private)} rows, expected 134")
    require(queue["review_id"].is_unique, "public review IDs are duplicated")
    require(private["review_id"].is_unique, "private review IDs are duplicated")
    require(
        set(queue["review_id"].astype(str)) == set(private["review_id"].astype(str)),
        "public/private review identities differ",
    )
    hidden = ("label", "franklin", "legacy", "prior", "truth", "recovery")
    exposed = [
        column for column in queue.columns if any(token in column.lower() for token in hidden)
    ]
    require(not exposed, f"public queue exposes hidden columns: {exposed}")
    signals = private["label"].isin(
        {
            "planet_like",
            "eclipsing_binary_or_pceb",
            "stellar_variability",
            "wide_transit_like",
        }
    )
    require(int(signals.sum()) == 34, "private manifest does not contain 34 signal rows")
    require(int((~signals).sum()) == 100, "private manifest does not contain 100 controls")
    pngs = list(args.sheet_dir.glob("*.png"))
    pdfs = list(args.sheet_dir.glob("*.pdf"))
    require(len(pngs) == 134, f"sheet directory has {len(pngs)} PNGs, expected 134")
    require(not pdfs, f"sheet directory unexpectedly has {len(pdfs)} PDFs")
    status_counts = render.get("status_counts", {})
    n_success = int(status_counts.get("ok", 0)) + int(
        status_counts.get("reused", 0)
    )
    require(n_success == 134, "render summary does not report 134 successful sheets")
    result = {
        "passed": not failures,
        "n_rows": int(len(queue)),
        "n_signal_rows": int(signals.sum()),
        "n_control_rows": int((~signals).sum()),
        "n_png": int(len(pngs)),
        "n_pdf": int(len(pdfs)),
        "failures": failures,
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
