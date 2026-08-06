#!/usr/bin/env python3
"""Build the prospective S63 rank-one ADP-small candidate table."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.vetting.teacher_v3_prospective import (  # noqa: E402
    build_s63_rank_one_candidates,
    file_sha256,
    load_prospective_plan,
    read_tic_inventory,
    tic_inventory_sha256,
    validate_git_sha,
    validate_s63_bls_summary,
    validate_sha256,
    write_s63_rank_one_candidates,
)


def _read_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".parquet":
        return pd.read_parquet(path)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path, low_memory=False)
    raise ValueError(f"unsupported BLS table format: {path}")


def _bound_json(path: Path, expected_sha256: str) -> dict[str, object]:
    expected = validate_sha256(expected_sha256, context=f"expected hash for {path}")
    before = file_sha256(path)
    if before != expected:
        raise ValueError(f"SHA-256 mismatch for {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if file_sha256(path) != before:
        raise RuntimeError(f"artifact changed while it was read: {path}")
    if not isinstance(payload, dict):
        raise ValueError(f"JSON artifact must be an object: {path}")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--expected-plan-sha256", required=True)
    parser.add_argument("--model-ready-allowlist", type=Path, required=True)
    parser.add_argument("--expected-model-ready-allowlist-sha256", required=True)
    parser.add_argument("--bls-peaks", type=Path, required=True)
    parser.add_argument("--expected-bls-peaks-sha256", required=True)
    parser.add_argument("--bls-summary", type=Path, required=True)
    parser.add_argument("--expected-bls-summary-sha256", required=True)
    parser.add_argument("--producer-git-sha", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    producer_git_sha = validate_git_sha(
        args.producer_git_sha, context="candidate producer Git SHA"
    )

    _, plan_sha256 = load_prospective_plan(
        args.plan,
        expected_sha256=args.expected_plan_sha256,
    )
    model_tics = read_tic_inventory(
        args.model_ready_allowlist,
        expected_sha256=args.expected_model_ready_allowlist_sha256,
    )
    expected_peaks_sha = validate_sha256(
        args.expected_bls_peaks_sha256,
        context="expected S63 BLS peak table",
    )
    before = file_sha256(args.bls_peaks)
    if before != expected_peaks_sha:
        raise ValueError("S63 BLS peak-table SHA-256 mismatch")
    peaks = _read_table(args.bls_peaks)
    if file_sha256(args.bls_peaks) != before:
        raise RuntimeError("S63 BLS peak table changed while it was read")
    bls_summary = _bound_json(args.bls_summary, args.expected_bls_summary_sha256)
    validate_s63_bls_summary(
        bls_summary,
        peak_table_sha256=expected_peaks_sha,
        model_ready_allowlist_sha256=args.expected_model_ready_allowlist_sha256,
        model_ready_tics_sha256=tic_inventory_sha256(set(model_tics)),
        n_model_ready_tics=len(model_tics),
    )
    candidates, summary = build_s63_rank_one_candidates(
        peaks,
        model_ready_tics=model_tics,
        artifact_hashes={
            "prospective_plan_sha256": plan_sha256,
            "model_ready_allowlist_sha256": args.expected_model_ready_allowlist_sha256,
            "bls_merged_table_sha256": expected_peaks_sha,
        },
    )
    output = write_s63_rank_one_candidates(
        candidates,
        summary,
        out_path=args.out,
        producer_git_sha=producer_git_sha,
    )
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
