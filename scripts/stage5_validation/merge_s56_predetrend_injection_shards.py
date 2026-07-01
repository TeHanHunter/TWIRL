#!/usr/bin/env python3
"""Merge metadata tables from sharded S56 pre-detrend injection outputs."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any


DEFAULT_CHUNK_ROOT = (
    Path("data_local/stage3_injections/s56_twirlfs_v2_injection_training/")
    / "pdo_allhost_predetrend_batman_periodradius_grid_sharded/chunks"
)
DEFAULT_OUT_DIR = DEFAULT_CHUNK_ROOT.parent


def _json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                fieldnames.append(key)
                seen.add(key)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def merge_shard_metadata(*, chunk_root: Path, out_dir: Path) -> dict[str, Any]:
    chunk_dirs = sorted(path for path in chunk_root.glob("chunk_*") if path.is_dir())
    if not chunk_dirs:
        raise FileNotFoundError(f"no chunk_* directories under {chunk_root}")
    manifest_rows: list[dict[str, Any]] = []
    label_rows: list[dict[str, Any]] = []
    chunk_summaries: list[dict[str, Any]] = []
    skipped_totals: dict[str, int] = {}
    source_targets_total = 0

    for chunk_dir in chunk_dirs:
        manifest_csv = chunk_dir / "injection_manifest.csv"
        labels_csv = chunk_dir / "injection_labels.csv"
        summary_json = chunk_dir / "summary.json"
        if not manifest_csv.exists() or not labels_csv.exists() or not summary_json.exists():
            continue
        summary = _load_json(summary_json)
        chunk_name = chunk_dir.name
        for row in _read_csv(manifest_csv):
            row["source_chunk"] = chunk_name
            row["source_h5"] = str(chunk_dir / "injected_lightcurves.h5")
            manifest_rows.append(row)
        for row in _read_csv(labels_csv):
            row["source_chunk"] = chunk_name
            row["source_h5"] = str(chunk_dir / "injected_lightcurves.h5")
            label_rows.append(row)
        chunk_summaries.append(
            {
                "chunk": chunk_name,
                "n_source_targets": int(summary.get("n_source_targets", 0)),
                "n_injections": int(summary.get("n_injections", 0)),
                "n_unique_injected_tics": int(summary.get("n_unique_injected_tics", 0)),
                "all_source_targets_injected": bool(summary.get("all_source_targets_injected", False)),
                "summary_json": str(summary_json),
                "h5": str(chunk_dir / "injected_lightcurves.h5"),
            }
        )
        source_targets_total += int(summary.get("n_source_targets", 0))
        skipped = summary.get("skipped", {})
        if isinstance(skipped, dict):
            for key, value in skipped.items():
                skipped_totals[str(key)] = skipped_totals.get(str(key), 0) + int(value)

    if not manifest_rows:
        raise FileNotFoundError(f"no completed shard metadata under {chunk_root}")
    unique_tics = {int(row["tic"]) for row in manifest_rows if str(row.get("tic", "")).strip()}
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_out = out_dir / "injection_manifest.csv"
    labels_out = out_dir / "injection_labels.csv"
    _write_csv(manifest_out, manifest_rows)
    _write_csv(labels_out, label_rows)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "chunk_root": str(chunk_root),
        "manifest_csv": str(manifest_out),
        "labels_csv": str(labels_out),
        "n_completed_shards": int(len(chunk_summaries)),
        "n_source_targets": int(source_targets_total),
        "n_injections": int(len(manifest_rows)),
        "n_unique_injected_tics": int(len(unique_tics)),
        "unique_injected_fraction_of_source_targets": (
            float(len(unique_tics) / source_targets_total) if source_targets_total else None
        ),
        "all_source_targets_injected": bool(source_targets_total > 0 and len(unique_tics) >= source_targets_total),
        "skipped": skipped_totals,
        "chunk_summaries": chunk_summaries,
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n")
    return summary


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chunk-root", type=Path, default=DEFAULT_CHUNK_ROOT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    summary = merge_shard_metadata(chunk_root=args.chunk_root, out_dir=args.out_dir)
    print("[merge-injection-shards] complete")
    print(f"  completed shards: {summary['n_completed_shards']:,}")
    print(f"  injections: {summary['n_injections']:,}")
    print(f"  unique TICs: {summary['n_unique_injected_tics']:,}")
    print(f"  manifest: {summary['manifest_csv']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
