#!/usr/bin/env python3
"""Merge and provenance-verify resumable ADP-only real-BLS shards."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from twirl.vetting.adp_only import ADP_ONLY_APERTURES


_INVARIANT_SUMMARY_FIELDS = (
    "sector",
    "contract_version",
    "bls_search_contract_version",
    "bls_config_sha256",
    "external_quality_policy_contract",
    "compact_lc",
    "compact_lc_sha256",
    "cadence_reference",
    "cadence_reference_sha256",
    "cadence_reference_manifest",
    "cadence_reference_manifest_sha256",
    "cadence_reference_contract_version",
    "cadence_reference_cadence_authority",
    "cadence_reference_quality_authority",
    "cadence_reference_source_hashes_sha256",
    "apertures",
    "n_targets_total",
    "n_periods",
    "n_peaks",
    "source_product_tag",
    "config",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def merge_shards(
    *,
    shard_dir: Path,
    out_path: Path,
    n_shards: int,
) -> dict[str, object]:
    if int(n_shards) < 1:
        raise ValueError("n_shards must be positive")
    shard_dir = Path(shard_dir)
    out_path = Path(out_path)
    paths = [
        shard_dir / f"real_adp_bls_peaks_{index:03d}.parquet"
        for index in range(n_shards)
    ]
    summaries = [
        shard_dir / f"summary_{index:03d}.json" for index in range(n_shards)
    ]
    missing = [str(path) for path in [*paths, *summaries] if not path.exists()]
    if missing:
        raise FileNotFoundError(
            f"missing {len(missing)} real-BLS shard products; first={missing[:5]}"
        )
    initial_hashes = {
        str(path.resolve()): _sha256(path) for path in [*paths, *summaries]
    }
    summary_payloads = [json.loads(path.read_text()) for path in summaries]
    baseline = summary_payloads[0]
    for field in _INVARIANT_SUMMARY_FIELDS:
        if field not in baseline:
            raise ValueError(f"real-BLS shard summary is missing {field}")
    frames: list[pd.DataFrame] = []
    target_sets: list[set[int]] = []
    expected_columns: tuple[str, ...] | None = None
    for index, (path, summary_path, payload) in enumerate(
        zip(paths, summaries, summary_payloads, strict=True)
    ):
        if int(payload.get("shard_index", -1)) != index:
            raise ValueError(f"real-BLS shard {index} has the wrong shard_index")
        if int(payload.get("n_shards", -1)) != int(n_shards):
            raise ValueError(f"real-BLS shard {index} has the wrong n_shards")
        for field in _INVARIANT_SUMMARY_FIELDS:
            if payload.get(field) != baseline[field]:
                raise ValueError(
                    f"real-BLS shard {index} disagrees on provenance field {field}"
                )
        if payload.get("peak_table_sha256") != initial_hashes[str(path.resolve())]:
            raise ValueError(f"real-BLS shard {index} peak-table SHA256 mismatch")
        outputs = payload.get("outputs", {})
        if not isinstance(outputs, dict) or Path(
            str(outputs.get("peak_table", ""))
        ).resolve() != path.resolve():
            raise ValueError(f"real-BLS shard {index} output path mismatch")
        frame = pd.read_parquet(path)
        columns = tuple(str(value) for value in frame.columns)
        if expected_columns is None:
            expected_columns = columns
        elif columns != expected_columns:
            raise ValueError("real-BLS shards have different table schemas")
        required_columns = {"tic", "aperture", "peak_rank", "status"}
        missing_columns = sorted(required_columns - set(frame.columns))
        if missing_columns:
            raise ValueError(
                f"real-BLS shard {index} lacks required columns: {missing_columns}"
            )
        if int(payload.get("n_rows", -1)) != len(frame):
            raise ValueError(f"real-BLS shard {index} row-count mismatch")
        tic = pd.to_numeric(frame.get("tic"), errors="coerce")
        if tic.isna().any() or (tic <= 0).any() or (tic != np.floor(tic)).any():
            raise ValueError(f"real-BLS shard {index} has invalid TIC values")
        targets = set(tic.astype(np.int64).tolist())
        if len(targets) != int(payload.get("n_targets", -1)):
            raise ValueError(f"real-BLS shard {index} target-count mismatch")
        if any(targets & previous for previous in target_sets):
            raise ValueError("real-BLS shard target sets overlap")
        target_sets.append(targets)
        frames.append(frame)

    frame = pd.concat(frames, ignore_index=True)
    apertures = set(frame["aperture"].dropna().astype(str))
    if not apertures.issubset(set(ADP_ONLY_APERTURES)):
        raise ValueError(
            f"non-ADP apertures found in real-BLS shards: {sorted(apertures)}"
        )
    if sorted(apertures) != sorted(str(value) for value in baseline["apertures"]):
        raise ValueError("merged real-BLS apertures disagree with shard summaries")
    key = ["tic", "aperture", "peak_rank"]
    valid_key = frame[key].notna().all(axis=1)
    if frame.loc[valid_key, key].duplicated().any():
        raise ValueError("real-BLS shards contain duplicate TIC/aperture/peak-rank rows")
    expected_targets = int(baseline["n_targets_total"])
    observed_targets = len(set().union(*target_sets))
    if observed_targets != expected_targets:
        raise ValueError(
            f"real-BLS target coverage is incomplete: {observed_targets} != {expected_targets}"
        )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    frame = frame.sort_values(key, kind="stable").reset_index(drop=True)
    temporary = out_path.with_suffix(out_path.suffix + ".tmp.parquet")
    temporary.unlink(missing_ok=True)
    frame.to_parquet(temporary, compression="zstd", index=False)
    final_hashes = {
        str(path.resolve()): _sha256(path) for path in [*paths, *summaries]
    }
    if final_hashes != initial_hashes:
        temporary.unlink(missing_ok=True)
        raise RuntimeError("real-BLS shard inputs changed during merge")
    temporary.replace(out_path)
    status = frame["status"].fillna("").astype(str)
    valid = status.eq("ok") & pd.to_numeric(
        frame.get("peak_rank"), errors="coerce"
    ).gt(0)
    summary: dict[str, Any] = {
        field: baseline[field] for field in _INVARIANT_SUMMARY_FIELDS
    }
    summary.update(
        {
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "out_dir": str(out_path.parent),
            "n_shards": 1,
            "shard_index": 0,
            "n_source_shards": int(n_shards),
            "n_targets": observed_targets,
            "n_targets_total": expected_targets,
            "n_rows": int(len(frame)),
            "n_unique_tics": observed_targets,
            "n_valid_peak_rows": int(valid.sum()),
            "workers": 0,
            "peak_table_sha256": _sha256(out_path),
            "status_counts": {
                str(key): int(value)
                for key, value in status.value_counts().sort_index().items()
            },
            "aperture_counts": {
                str(key): int(value)
                for key, value in frame["aperture"]
                .fillna("")
                .astype(str)
                .value_counts()
                .sort_index()
                .items()
            },
            "source_shards": [
                {
                    "shard_index": index,
                    "peak_table": str(path),
                    "peak_table_sha256": initial_hashes[str(path.resolve())],
                    "summary": str(summary_path),
                    "summary_sha256": initial_hashes[
                        str(summary_path.resolve())
                    ],
                    "n_targets": int(payload["n_targets"]),
                    "n_rows": int(payload["n_rows"]),
                }
                for index, (path, summary_path, payload) in enumerate(
                    zip(paths, summaries, summary_payloads, strict=True)
                )
            ],
            "outputs": {
                "peak_table": str(out_path),
                "summary": str(out_path.with_suffix(".summary.json")),
            },
            "passed": True,
        }
    )
    summary_path = out_path.with_suffix(".summary.json")
    summary_tmp = summary_path.with_suffix(summary_path.suffix + ".tmp")
    summary_tmp.write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    summary_tmp.replace(summary_path)
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
