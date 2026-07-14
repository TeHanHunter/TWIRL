#!/usr/bin/env python3
"""Merge and verify resumable ADP-only real-BLS shards."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import pandas as pd

from twirl.vetting.adp_only import ADP_ONLY_APERTURES


def merge_shards(
    *,
    shard_dir: Path,
    out_path: Path,
    n_shards: int,
) -> dict[str, object]:
    paths = [shard_dir / f"real_adp_bls_peaks_{index:03d}.parquet" for index in range(n_shards)]
    summaries = [shard_dir / f"summary_{index:03d}.json" for index in range(n_shards)]
    missing = [str(path) for path in [*paths, *summaries] if not path.exists()]
    if missing:
        raise FileNotFoundError(f"missing {len(missing)} real-BLS shard products; first={missing[:5]}")
    summary_payloads = [json.loads(path.read_text()) for path in summaries]
    for index, payload in enumerate(summary_payloads):
        if int(payload.get("shard_index", -1)) != index:
            raise ValueError(f"real-BLS shard {index} has the wrong shard_index")
        if int(payload.get("n_shards", -1)) != int(n_shards):
            raise ValueError(f"real-BLS shard {index} has the wrong n_shards")
    frame = pd.concat([pd.read_parquet(path) for path in paths], ignore_index=True)
    apertures = set(frame["aperture"].dropna().astype(str))
    if not apertures.issubset(set(ADP_ONLY_APERTURES)):
        raise ValueError(f"non-ADP apertures found in real-BLS shards: {sorted(apertures)}")
    key = ["tic", "aperture", "peak_rank"]
    valid_key = frame[key].notna().all(axis=1)
    if frame.loc[valid_key, key].duplicated().any():
        raise ValueError("real-BLS shards contain duplicate TIC/aperture/peak-rank rows")
    expected_total = {int(payload["n_targets_total"]) for payload in summary_payloads}
    if len(expected_total) != 1:
        raise ValueError("real-BLS shards disagree on the total target count")
    expected_targets = expected_total.pop()
    observed_targets = int(pd.to_numeric(frame["tic"], errors="coerce").nunique())
    if observed_targets != expected_targets:
        raise ValueError(
            f"real-BLS target coverage is incomplete: {observed_targets} != {expected_targets}"
        )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    frame = frame.sort_values(key, kind="stable").reset_index(drop=True)
    frame.to_parquet(out_path, compression="zstd", index=False)
    status = frame["status"].fillna("").astype(str)
    summary: dict[str, object] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "n_shards": int(n_shards),
        "n_rows": int(len(frame)),
        "n_targets": observed_targets,
        "apertures": sorted(apertures),
        "status_counts": {str(key): int(value) for key, value in status.value_counts().items()},
        "out_path": str(out_path),
        "passed": True,
    }
    out_path.with_suffix(".summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shard-dir", type=Path, required=True)
    parser.add_argument("--out-path", type=Path, required=True)
    parser.add_argument("--n-shards", type=int, required=True)
    args = parser.parse_args()
    summary = merge_shards(
        shard_dir=args.shard_dir,
        out_path=args.out_path,
        n_shards=args.n_shards,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
