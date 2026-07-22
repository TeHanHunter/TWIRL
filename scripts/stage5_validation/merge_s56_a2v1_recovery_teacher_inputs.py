#!/usr/bin/env python3
"""Merge all fresh-recovery candidate and native-input shards."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import pandas as pd

from twirl.vetting.harmonic_export import merge_raw_pair_shards
from twirl.vetting.harmonic_inputs import native_group_path, verify_raw_pair_contract


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-shards", type=Path, nargs="+", required=True)
    parser.add_argument("--native-shards", type=Path, nargs="+", required=True)
    parser.add_argument("--out-candidates", type=Path, required=True)
    parser.add_argument("--out-native-h5", type=Path, required=True)
    parser.add_argument("--maximum-injections", type=int, default=20_000)
    args = parser.parse_args()

    candidate_shards = tuple(sorted(path.resolve() for path in args.candidate_shards))
    native_shards = tuple(sorted(path.resolve() for path in args.native_shards))
    if len(candidate_shards) != len(native_shards):
        raise ValueError("candidate/native shard counts differ")
    candidate_shard_sha256 = {
        str(path): _sha256(path) for path in candidate_shards
    }
    native_shard_sha256 = {str(path): _sha256(path) for path in native_shards}
    native_shard_training: list[dict[str, str]] = []
    for candidate_path, native_path in zip(
        candidate_shards, native_shards, strict=True
    ):
        candidate_tag = candidate_path.stem.removeprefix("candidates_")
        native_tag = native_path.stem.removeprefix("native_")
        if candidate_tag != native_tag:
            raise ValueError(
                f"candidate/native shard tags differ: {candidate_path}, {native_path}"
            )
        with h5py.File(native_path, "r") as h5:
            declared_table = Path(str(h5.attrs.get("training_table", ""))).resolve()
            declared_sha256 = str(h5.attrs.get("training_table_sha256", ""))
        if declared_table != candidate_path:
            raise ValueError(
                f"native shard {native_path} is not bound to {candidate_path}"
            )
        if declared_sha256 != candidate_shard_sha256[str(candidate_path)]:
            raise ValueError(
                f"native shard {native_path} candidate-table SHA256 mismatch"
            )
        native_shard_training.append(
            {
                "candidate_shard": str(candidate_path),
                "candidate_shard_sha256": declared_sha256,
                "native_shard": str(native_path),
            }
        )
    if native_shard_sha256 != {
        str(path): _sha256(path) for path in native_shards
    }:
        raise RuntimeError("native shards changed while their table binding was read")

    candidates = pd.concat(
        [pd.read_csv(path, low_memory=False) for path in candidate_shards],
        ignore_index=True,
    )
    if candidates["review_id"].duplicated().any():
        raise ValueError("candidate shards contain duplicate review_id values")
    n_candidate_injections = int(candidates["injection_id"].nunique())
    if not 0 < n_candidate_injections <= args.maximum_injections:
        raise ValueError(
            f"candidate shards cover {n_candidate_injections:,} injections; "
            f"expected between 1 and {args.maximum_injections:,}"
        )
    args.out_candidates.parent.mkdir(parents=True, exist_ok=True)
    candidate_temporary = args.out_candidates.with_suffix(
        args.out_candidates.suffix + ".tmp"
    )
    candidate_temporary.unlink(missing_ok=True)
    candidates.to_csv(candidate_temporary, index=False)
    candidate_temporary.replace(args.out_candidates)
    merge = merge_raw_pair_shards(
        shard_paths=native_shards,
        out_h5=args.out_native_h5,
        merged_training_table=args.out_candidates,
    )
    verification = verify_raw_pair_contract(
        args.out_native_h5,
        require_errors=True,
        require_periodograms=True,
    )
    expected_counts = {"targets": 0, "injections": n_candidate_injections}
    expected_groups = {
        native_group_path(row) for row in candidates.to_dict("records")
    }
    with h5py.File(args.out_native_h5, "r") as h5:
        observed_groups = {
            f"{root}/{key}"
            for root in ("targets", "injections")
            if root in h5
            for key in h5[root]
        }
        merged_training_sha256 = str(h5.attrs.get("training_table_sha256", ""))
    exact_group_match = observed_groups == expected_groups
    merged_candidate_sha256 = _sha256(args.out_candidates)
    source_candidates_unchanged = candidate_shard_sha256 == {
        str(path): _sha256(path) for path in candidate_shards
    }
    source_native_unchanged = native_shard_sha256 == {
        str(path): _sha256(path) for path in native_shards
    }
    passed = (
        verification["passed"]
        and merge["counts"] == expected_counts
        and exact_group_match
        and source_candidates_unchanged
        and source_native_unchanged
        and merge["shard_sha256"] == native_shard_sha256
        and merged_training_sha256 == merged_candidate_sha256
    )
    summary = {
        "n_candidate_rows": int(len(candidates)),
        "n_injections": n_candidate_injections,
        "merge": merge,
        "verification": verification,
        "expected_counts": expected_counts,
        "exact_group_match": exact_group_match,
        "candidate_shard_sha256": candidate_shard_sha256,
        "native_shard_sha256": native_shard_sha256,
        "native_shard_training": native_shard_training,
        "source_candidates_unchanged": source_candidates_unchanged,
        "source_native_unchanged": source_native_unchanged,
        "merged_candidate_sha256": merged_candidate_sha256,
        "merged_native_training_table_sha256": merged_training_sha256,
        "passed": bool(passed),
    }
    args.out_native_h5.with_suffix(".merge_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if passed else 3


if __name__ == "__main__":
    raise SystemExit(main())
