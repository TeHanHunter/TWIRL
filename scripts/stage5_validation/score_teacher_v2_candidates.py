#!/usr/bin/env python3
"""Score real or injected candidates with a frozen Teacher-v2 ensemble."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from twirl.vetting.teacher_v2_inference import score_teacher_v2_ensemble


def _read(path: Path) -> pd.DataFrame:
    return pd.read_parquet(path) if path.suffix.lower() == ".parquet" else pd.read_csv(path, low_memory=False)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--native-h5", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, action="append", required=True)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--out-path", type=Path, required=True)
    parser.add_argument("--filter-manifest", type=Path)
    parser.add_argument("--filter-key", default="injection_id")
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--allow-cpu", action="store_true")
    args = parser.parse_args()

    candidates = _read(args.candidates)
    if args.filter_manifest:
        manifest = _read(args.filter_manifest)
        if args.filter_key not in candidates or args.filter_key not in manifest:
            raise KeyError(f"filter key {args.filter_key!r} is absent from candidates or manifest")
        allowed = frozenset(manifest[args.filter_key].fillna("").astype(str))
        candidates = candidates.loc[
            candidates[args.filter_key].fillna("").astype(str).isin(allowed)
        ].copy()
    if args.limit is not None:
        candidates = candidates.head(max(0, int(args.limit))).copy()
    scored, summary = score_teacher_v2_ensemble(
        candidates=candidates,
        native_h5=args.native_h5,
        checkpoint_paths=args.checkpoint,
        profile=args.profile,
        batch_size=args.batch_size,
        workers=args.workers,
        require_cuda=not args.allow_cpu,
    )
    args.out_path.parent.mkdir(parents=True, exist_ok=True)
    if args.out_path.suffix.lower() == ".parquet":
        scored.to_parquet(args.out_path, compression="zstd", index=False)
    else:
        scored.to_csv(args.out_path, index=False)
    summary["out_path"] = str(args.out_path)
    args.out_path.with_suffix(".summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
