#!/usr/bin/env python3
"""Score fresh A2v1 injection candidates with the frozen Teacher-v1 ensemble."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from twirl.vetting.harmonic_inference import score_harmonic_teacher_ensemble


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--native-h5", type=Path, required=True)
    parser.add_argument("--checkpoints", type=Path, nargs=5, required=True)
    parser.add_argument("--out-scores", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--filter-manifest", type=Path)
    parser.add_argument("--filter-key", default="injection_id")
    parser.add_argument("--allow-cpu", action="store_true")
    args = parser.parse_args()

    candidates = (
        pd.read_parquet(args.candidates)
        if args.candidates.suffix.lower() == ".parquet"
        else pd.read_csv(args.candidates, low_memory=False)
    )
    if args.filter_manifest:
        manifest = (
            pd.read_parquet(args.filter_manifest)
            if args.filter_manifest.suffix.lower() == ".parquet"
            else pd.read_csv(args.filter_manifest, low_memory=False)
        )
        if args.filter_key not in candidates or args.filter_key not in manifest:
            raise KeyError(
                f"filter key {args.filter_key!r} is absent from candidates or manifest"
            )
        allowed = frozenset(manifest[args.filter_key].fillna("").astype(str))
        candidates = candidates.loc[
            candidates[args.filter_key].fillna("").astype(str).isin(allowed)
        ].copy()
    if args.limit is not None:
        candidates = candidates.head(int(args.limit)).copy()
    scores, summary = score_harmonic_teacher_ensemble(
        candidates=candidates,
        native_h5=args.native_h5,
        checkpoint_paths=args.checkpoints,
        batch_size=args.batch_size,
        workers=args.workers,
        require_cuda=not args.allow_cpu,
        allow_injections=True,
    )
    args.out_scores.parent.mkdir(parents=True, exist_ok=True)
    scores.to_parquet(args.out_scores, compression="zstd", index=False)
    summary["out_scores"] = str(args.out_scores)
    args.out_scores.with_suffix(".summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
