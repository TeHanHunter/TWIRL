#!/usr/bin/env python3
"""Audit the exact 34-observation detector-consistency BLS v4 canary."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from twirl.vetting.ssl_full_pool_native import shard_for_tic


EXPECTED_COUNTS = {56: 7, 58: 3, 61: 16, 62: 8}
EXPECTED_APERTURES = {"DET_FLUX_ADP_SML", "DET_FLUX_ADP"}
INPUT_MODE = "immutable_raw_v1_detector_consistent"
INPUT_CONTRACT = "twirl_teacher_ssl_fullpool_raw_v1_detector_consistent_bls_v4"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _expected(config_root: Path) -> dict[tuple[int, int], list[int]]:
    grouped: dict[tuple[int, int], list[int]] = {}
    observed_counts: dict[int, int] = {}
    for path in sorted(config_root.glob("s*_affected.csv")):
        sector = int(path.stem[1:3])
        tics = pd.read_csv(path)["tic"].astype(np.int64).tolist()
        if len(tics) != len(set(tics)) or any(tic <= 0 for tic in tics):
            raise ValueError(f"invalid canary TIC inventory in {path}")
        observed_counts[sector] = len(tics)
        for tic in tics:
            shard = shard_for_tic(sector=sector, tic=int(tic), n_shards=16)
            grouped.setdefault((sector, shard), []).append(int(tic))
    if observed_counts != EXPECTED_COUNTS:
        raise ValueError(
            f"canary sector counts differ: {observed_counts} != {EXPECTED_COUNTS}"
        )
    return {key: sorted(values) for key, values in sorted(grouped.items())}


def audit(*, canary_root: Path, config_root: Path) -> dict[str, object]:
    expected = _expected(config_root)
    records: list[dict[str, object]] = []
    observed: set[tuple[int, int]] = set()
    max_time_delta = 0.0
    for (sector, shard), expected_tics in expected.items():
        shard_root = canary_root / f"s{sector}" / "shards"
        table_path = shard_root / f"real_adp_bls_peaks_{shard:03d}.parquet"
        summary_path = shard_root / f"summary_{shard:03d}.json"
        if not table_path.is_file() or not summary_path.is_file():
            raise FileNotFoundError(
                f"missing S{sector} shard {shard} canary products"
            )
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        expected_summary = {
            "sector": sector,
            "shard_index": shard,
            "n_shards": 16,
            "n_targets": len(expected_tics),
            "n_targets_total": EXPECTED_COUNTS[sector],
            "n_periods": 50_000,
            "n_peaks": 10,
            "input_mode": INPUT_MODE,
            "raw_v4_input_contract_version": INPUT_CONTRACT,
            "raw_compact_cadence_inventory_passed": True,
            "n_raw_compact_cadence_inventories_verified": len(expected_tics),
            "source_product_tag": "A2v1-fullpool-v1",
        }
        mismatches = {
            name: {"expected": value, "observed": summary.get(name)}
            for name, value in expected_summary.items()
            if summary.get(name) != value
        }
        if mismatches:
            raise ValueError(
                f"S{sector} shard {shard} summary mismatch: "
                + json.dumps(mismatches, sort_keys=True)
            )
        if summary.get("peak_table_sha256") != _sha256(table_path):
            raise ValueError(f"S{sector} shard {shard} peak hash mismatch")
        delta = float(summary["raw_compact_time_delta_max_s"])
        if not np.isfinite(delta) or delta > 2.0:
            raise ValueError(f"S{sector} shard {shard} time audit failed")
        max_time_delta = max(max_time_delta, delta)

        frame = pd.read_parquet(table_path)
        frame_tics = set(frame["tic"].astype(np.int64).tolist())
        if frame_tics != set(expected_tics):
            raise ValueError(f"S{sector} shard {shard} target coverage differs")
        if set(frame["aperture"].astype(str)) != EXPECTED_APERTURES:
            raise ValueError(f"S{sector} shard {shard} aperture coverage differs")
        if not frame["status"].astype(str).eq("ok").all():
            raise ValueError(f"S{sector} shard {shard} has non-ok BLS status")
        for (tic, aperture), group in frame.groupby(
            ["tic", "aperture"], sort=False
        ):
            ranks = sorted(group["peak_rank"].astype(int).tolist())
            if ranks != list(range(1, 11)):
                raise ValueError(
                    f"S{sector} TIC {tic} {aperture} peak ranks differ"
                )
            observed.add((sector, int(tic)))
        records.append(
            {
                "sector": sector,
                "shard_index": shard,
                "tics": expected_tics,
                "peak_table_sha256": _sha256(table_path),
                "summary_sha256": _sha256(summary_path),
            }
        )
    expected_keys = {
        (sector, tic)
        for (sector, _shard), tics in expected.items()
        for tic in tics
    }
    if observed != expected_keys or len(observed) != 34:
        raise ValueError("canary union does not equal the exact affected set")
    return {
        "passed": True,
        "schema_version": "twirl_teacher_ssl_fullpool_v4_bls_canary_audit_v1",
        "n_observations": len(observed),
        "n_source_shards": len(expected),
        "by_sector": EXPECTED_COUNTS,
        "raw_compact_time_delta_max_s": max_time_delta,
        "records": records,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--canary-root", type=Path, required=True)
    parser.add_argument("--config-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = audit(
        canary_root=args.canary_root.resolve(),
        config_root=args.config_root.resolve(),
    )
    payload = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(payload, encoding="utf-8")
    print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
