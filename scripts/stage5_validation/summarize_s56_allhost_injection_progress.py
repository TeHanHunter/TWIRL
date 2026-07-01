#!/usr/bin/env python3
"""Summarize the restartable S56 all-host injection build progress."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any


DEFAULT_OUT_ROOT = (
    Path("data_local/stage3_injections/s56_twirlfs_v2_injection_training")
    / "pdo_allhost_predetrend_batman_periodradius_grid_sharded"
)


def _load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text())


def _chunk_name(path: Path) -> str:
    return path.name


def summarize_progress(out_root: Path) -> dict[str, Any]:
    shard_manifest_path = out_root / "target_shards" / "manifest.json"
    shard_manifest = _load_json(shard_manifest_path)
    shards = shard_manifest.get("shards", [])
    shard_names = [str(row.get("shard", "")) for row in shards if row.get("shard")]
    shard_name_set = set(shard_names)

    chunk_root = out_root / "chunks"
    chunk_dirs = sorted(path for path in chunk_root.glob("chunk_*") if path.is_dir())
    completed: list[str] = []
    completed_with_h5: list[str] = []
    started: list[str] = []
    skipped_totals: dict[str, int] = {}
    n_injections = 0
    n_unique_tics = 0

    for chunk_dir in chunk_dirs:
        started.append(_chunk_name(chunk_dir))
        summary_path = chunk_dir / "summary.json"
        if not summary_path.exists():
            continue
        completed.append(_chunk_name(chunk_dir))
        if (chunk_dir / "injected_lightcurves.h5").exists():
            completed_with_h5.append(_chunk_name(chunk_dir))
        summary = _load_json(summary_path)
        n_injections += int(summary.get("n_injections", 0) or 0)
        n_unique_tics += int(summary.get("n_unique_injected_tics", 0) or 0)
        for key, value in (summary.get("skipped") or {}).items():
            skipped_totals[key] = skipped_totals.get(key, 0) + int(value or 0)

    completed_set = set(completed)
    started_set = set(started)
    pending = [name for name in shard_names if name not in started_set]
    active_or_incomplete = [name for name in started if name not in completed_set]
    unknown_started = [name for name in started if shard_name_set and name not in shard_name_set]
    n_shards = int(shard_manifest.get("n_shards", len(shard_names)) or 0)
    n_completed = len(completed)

    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "out_root": str(out_root),
        "shard_manifest": str(shard_manifest_path),
        "n_targets": int(shard_manifest.get("n_targets", 0) or 0),
        "n_paths": int(shard_manifest.get("n_paths", 0) or 0),
        "n_shards": n_shards,
        "n_started_shards": len(started),
        "n_completed_shards": n_completed,
        "n_completed_h5_shards": len(completed_with_h5),
        "n_active_or_incomplete_shards": len(active_or_incomplete),
        "n_pending_shards": len(pending),
        "completed_fraction": (float(n_completed) / float(n_shards)) if n_shards else 0.0,
        "n_injections_completed": n_injections,
        "n_unique_injected_tics_completed": n_unique_tics,
        "skipped_totals_completed": skipped_totals,
        "active_or_incomplete_shards": active_or_incomplete[:20],
        "pending_shards_head": pending[:20],
        "unknown_started_shards": unknown_started[:20],
        "combined_summary_exists": (out_root / "summary.json").exists(),
        "combined_manifest_exists": (out_root / "injection_manifest.csv").exists(),
    }
    return summary


def write_outputs(summary: dict[str, Any], out_dir: Path | None) -> None:
    if out_dir is None:
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "progress_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    lines = [
        "# S56 All-Host Injection Progress",
        "",
        f"- Out root: `{summary['out_root']}`",
        f"- Shards: `{summary['n_completed_shards']} / {summary['n_shards']}` complete",
        f"- Completed fraction: `{summary['completed_fraction']:.1%}`",
        f"- Started shards: `{summary['n_started_shards']}`",
        f"- Active/incomplete shards: `{summary['n_active_or_incomplete_shards']}`",
        f"- Pending shards: `{summary['n_pending_shards']}`",
        f"- Completed injections: `{summary['n_injections_completed']}`",
        f"- Completed unique TICs: `{summary['n_unique_injected_tics_completed']}`",
        f"- Combined summary exists: `{summary['combined_summary_exists']}`",
        "",
    ]
    if summary["active_or_incomplete_shards"]:
        lines.append("Active/incomplete shard head:")
        lines.extend(f"- `{name}`" for name in summary["active_or_incomplete_shards"])
        lines.append("")
    if summary["skipped_totals_completed"]:
        lines.append("Completed-shard skip totals:")
        for key, value in sorted(summary["skipped_totals_completed"].items()):
            lines.append(f"- `{key}`: `{value}`")
        lines.append("")
    (out_dir / "progress_summary.md").write_text("\n".join(lines) + "\n")


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-root", type=Path, default=DEFAULT_OUT_ROOT)
    parser.add_argument("--out-dir", type=Path, default=None)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    summary = summarize_progress(args.out_root)
    print(json.dumps(summary, indent=2, sort_keys=True))
    write_outputs(summary, args.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
