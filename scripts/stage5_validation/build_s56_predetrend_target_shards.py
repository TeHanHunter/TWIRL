#!/usr/bin/env python3
"""Build TIC-grouped raw-HDF5 target shards for S56 pre-detrend injections."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any


DEFAULT_ORBIT_ROOTS = (
    Path("/pdo/users/tehan/tglc-gpu-production/orbit-119/ffi"),
    Path("/pdo/users/tehan/tglc-gpu-production/orbit-120/ffi"),
)
DEFAULT_OUT_DIR = (
    Path("data_local/stage3_injections/s56_twirlfs_v2_injection_training/")
    / "pdo_allhost_predetrend_batman_periodradius_grid_sharded/target_shards"
)


def _json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _tic_from_path(path: Path) -> int | None:
    try:
        return int(Path(path).stem)
    except ValueError:
        return None


def _read_path_list(path: Path | None) -> list[Path]:
    if path is None:
        return []
    rows: list[Path] = []
    for raw in path.read_text().splitlines():
        text = raw.strip()
        if not text or text.startswith("#"):
            continue
        rows.append(Path(text))
    return rows


def discover_raw_h5_targets(
    *,
    orbit_roots: tuple[Path, ...],
    target_h5_paths: tuple[Path, ...] = (),
    limit_targets: int | None = None,
) -> dict[int, list[Path]]:
    """Return raw TGLC HDF5 paths grouped by TIC."""
    grouped: dict[int, list[Path]] = {}
    if target_h5_paths:
        candidates = list(target_h5_paths)
    else:
        candidates = []
        for root in orbit_roots:
            candidates.extend(sorted(Path(root).glob("cam*/ccd*/LC/*.h5")))

    for path in candidates:
        tic = _tic_from_path(path)
        if tic is None:
            continue
        grouped.setdefault(int(tic), []).append(Path(path))

    ordered = {tic: sorted(paths) for tic, paths in sorted(grouped.items())}
    if limit_targets is not None:
        keep = set(list(ordered)[: int(limit_targets)])
        ordered = {tic: paths for tic, paths in ordered.items() if tic in keep}
    return ordered


def write_target_shards(
    *,
    targets: dict[int, list[Path]],
    out_dir: Path,
    shard_size: int,
    overwrite: bool = False,
) -> dict[str, Any]:
    if shard_size <= 0:
        raise ValueError("shard_size must be positive")
    if not targets:
        raise ValueError("no targets to shard")
    out_dir.mkdir(parents=True, exist_ok=True)
    existing = sorted(out_dir.glob("chunk_*.txt")) + sorted(out_dir.glob("chunk_*.meta.json"))
    if existing and not overwrite:
        raise FileExistsError(f"{out_dir} already has shard files; pass --overwrite")
    if overwrite:
        for path in existing:
            path.unlink()

    tic_rows = [(tic, targets[tic]) for tic in sorted(targets)]
    shard_rows: list[dict[str, Any]] = []
    target_index_path = out_dir / "target_index.csv"
    with target_index_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=("tic", "n_paths", "paths"))
        writer.writeheader()
        for tic, paths in tic_rows:
            writer.writerow({"tic": tic, "n_paths": len(paths), "paths": "|".join(str(path) for path in paths)})

    for shard_index, start in enumerate(range(0, len(tic_rows), shard_size)):
        chunk = tic_rows[start : start + shard_size]
        shard_name = f"chunk_{shard_index:03d}"
        list_path = out_dir / f"{shard_name}.txt"
        meta_path = out_dir / f"{shard_name}.meta.json"
        lines = [str(path) for _, paths in chunk for path in paths]
        list_path.write_text("\n".join(lines) + "\n")
        tics = [tic for tic, _ in chunk]
        meta = {
            "shard": shard_name,
            "shard_index": int(shard_index),
            "start_index": int(start),
            "n_tics": int(len(tics)),
            "n_paths": int(len(lines)),
            "tic_min": int(min(tics)),
            "tic_max": int(max(tics)),
            "list_path": str(list_path),
        }
        meta_path.write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
        shard_rows.append(meta)

    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "out_dir": str(out_dir),
        "target_index_csv": str(target_index_path),
        "shard_size": int(shard_size),
        "n_targets": int(len(tic_rows)),
        "n_paths": int(sum(len(paths) for _, paths in tic_rows)),
        "n_shards": int(len(shard_rows)),
        "shards": shard_rows,
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True, default=_json_default) + "\n")
    return manifest


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orbit-roots", type=Path, nargs="+", default=list(DEFAULT_ORBIT_ROOTS))
    parser.add_argument("--target-h5-paths", type=Path, nargs="+", default=None)
    parser.add_argument("--target-h5-list", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--shard-size", type=int, default=250)
    parser.add_argument("--limit-targets", type=int, default=None)
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    target_h5_paths = tuple(args.target_h5_paths or ()) + tuple(_read_path_list(args.target_h5_list))
    targets = discover_raw_h5_targets(
        orbit_roots=tuple(args.orbit_roots),
        target_h5_paths=target_h5_paths,
        limit_targets=args.limit_targets,
    )
    manifest = write_target_shards(
        targets=targets,
        out_dir=args.out_dir,
        shard_size=args.shard_size,
        overwrite=args.overwrite,
    )
    print("[target-shards] complete")
    print(f"  targets: {manifest['n_targets']:,}")
    print(f"  paths: {manifest['n_paths']:,}")
    print(f"  shards: {manifest['n_shards']:,}")
    print(f"  manifest: {Path(manifest['out_dir']) / 'manifest.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
