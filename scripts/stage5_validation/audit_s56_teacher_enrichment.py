#!/usr/bin/env python3
"""Measure batch-one Planet enrichment and gate teacher inference on S57+."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.vetting.teacher_active_learning import sector_rollout_readiness  # noqa: E402


def read(path: Path) -> pd.DataFrame:
    return pd.read_parquet(path) if path.suffix.lower() == ".parquet" else pd.read_csv(path, low_memory=False)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hidden-provenance", type=Path, required=True)
    parser.add_argument("--labels", type=Path, required=True)
    parser.add_argument("--transfer-summary", type=Path, required=True)
    parser.add_argument("--product-qa-summary", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    hidden = read(args.hidden_provenance)
    labels = read(args.labels)
    labels["row_id"] = labels["row_id"].astype(str)
    hidden["row_id"] = hidden["row_id"].astype(str)
    latest = labels.sort_values("updated_utc", kind="stable").drop_duplicates("row_id", keep="last")
    joined = hidden.merge(
        latest.loc[:, ["row_id", "candidate_key", "label", "labeler", "updated_utc"]].rename(
            columns={"candidate_key": "human_candidate_key", "label": "human_label"}
        ),
        on="row_id",
        how="left",
        validate="one_to_one",
    )
    mismatch = joined["human_candidate_key"].fillna("").astype(str).ne(
        joined["candidate_key"].fillna("").astype(str)
    ) & joined["human_label"].fillna("").astype(str).ne("")
    if mismatch.any():
        raise ValueError(f"{int(mismatch.sum())} batch labels have candidate-key mismatches")
    transfer = json.loads(args.transfer_summary.read_text())
    product_qa = json.loads(args.product_qa_summary.read_text())
    transfer_passed = bool(
        transfer.get("scored_transfer", {}).get("transfer_gate", {}).get("passed", False)
    )
    summary = sector_rollout_readiness(
        joined,
        transfer_gate_passed=transfer_passed,
        product_qa_passed=bool(product_qa.get("passed", False)),
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=True))


if __name__ == "__main__":
    main()
