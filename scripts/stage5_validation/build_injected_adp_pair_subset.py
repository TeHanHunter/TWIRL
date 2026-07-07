#!/usr/bin/env python3
"""Build a queue-specific injected ADP 1x1 + 3x3 HDF5 subset."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))


def _target_key(tic: int) -> str:
    return f"{int(tic):016d}"


def _copy_attrs(src: Any, dst: Any) -> None:
    for key, value in src.attrs.items():
        dst.attrs[key] = value


def _dataset(group: Any, name: str, data: np.ndarray, compression: str | None) -> None:
    if name in group:
        del group[name]
    group.create_dataset(name, data=data, compression=compression)


def build_subset(
    *,
    queue_csv: Path,
    old_injection_h5: Path,
    compact_export_h5: Path,
    out_h5: Path,
    apertures: tuple[str, str],
    compression: str | None,
    overwrite: bool,
) -> dict[str, Any]:
    import h5py

    if out_h5.exists() and not overwrite:
        raise FileExistsError(f"output exists; pass --overwrite: {out_h5}")
    out_h5.parent.mkdir(parents=True, exist_ok=True)
    tmp_h5 = out_h5.with_suffix(out_h5.suffix + ".tmp")
    if tmp_h5.exists():
        tmp_h5.unlink()

    queue = pd.read_csv(queue_csv)
    if "h5_group" not in queue.columns:
        raise ValueError(f"queue lacks h5_group column: {queue_csv}")
    queue = queue[queue["h5_group"].notna()].copy()
    created_utc = datetime.now(timezone.utc).isoformat()
    rows: list[dict[str, Any]] = []
    missing_target = 0
    missing_aperture = 0
    time_mismatch = 0

    with h5py.File(old_injection_h5, "r") as old, h5py.File(compact_export_h5, "r") as export, h5py.File(tmp_h5, "w") as out:
        out.attrs["created_utc"] = created_utc
        out.attrs["source_injection_h5"] = str(old_injection_h5)
        out.attrs["source_compact_export_h5"] = str(compact_export_h5)
        out.attrs["queue_csv"] = str(queue_csv)
        out.attrs["aperture"] = apertures[0]
        out.attrs["apertures"] = json.dumps(list(apertures))
        out.attrs["product_note"] = "Queue subset rebuilt to carry injected ADP 1x1 and ADP 3x3 flux columns."
        out.require_group("injections")
        targets = export.get("targets")
        if targets is None:
            raise ValueError(f"compact export lacks /targets: {compact_export_h5}")

        for _, row in queue.iterrows():
            group_path = str(row["h5_group"])
            if group_path not in old:
                rows.append({"h5_group": group_path, "status": "missing_old_group"})
                continue
            old_group = old[group_path]
            tic = int(float(row.get("tic", old_group.attrs.get("tic", -1))))
            target_key = _target_key(tic)
            if target_key not in targets:
                missing_target += 1
                rows.append({"h5_group": group_path, "tic": tic, "status": "missing_target"})
                continue
            target = targets[target_key]
            missing = [ap for ap in apertures if ap not in target]
            if missing:
                missing_aperture += 1
                rows.append({"h5_group": group_path, "tic": tic, "status": "missing_aperture", "missing": "|".join(missing)})
                continue
            if "transit_model" not in old_group:
                rows.append({"h5_group": group_path, "tic": tic, "status": "missing_transit_model"})
                continue

            old_time = np.asarray(old_group["time"], dtype=np.float64)
            target_time = np.asarray(target["time"], dtype=np.float64)
            if len(old_time) != len(target_time) or np.nanmax(np.abs(old_time - target_time)) > 1.0e-8:
                time_mismatch += 1
                rows.append({"h5_group": group_path, "tic": tic, "status": "time_mismatch"})
                continue

            quality = np.asarray(target["quality"], dtype=np.int32)
            model = np.asarray(old_group["transit_model"], dtype=np.float64)
            out_group = out.create_group(group_path)
            _copy_attrs(old_group, out_group)
            out_group.attrs["aperture"] = apertures[0]
            out_group.attrs["apertures"] = json.dumps(list(apertures))
            out_group.attrs["source_injection_h5"] = str(old_injection_h5)
            out_group.attrs["source_compact_target"] = f"/targets/{target_key}"
            out_group.attrs["rebuilt_adp_pair"] = True

            _dataset(out_group, "time", target_time.astype(np.float64), compression)
            _dataset(out_group, "quality", quality.astype(np.int32), compression)
            orbitid = np.asarray(target["orbitid"] if "orbitid" in target else old_group["orbitid"], dtype=np.int16)
            _dataset(out_group, "orbitid", orbitid.astype(np.int16), compression)
            if "cadenceno" in target:
                _dataset(out_group, "cadenceno", np.asarray(target["cadenceno"], dtype=np.int32), compression)
            elif "cadenceno" in old_group:
                _dataset(out_group, "cadenceno", np.asarray(old_group["cadenceno"], dtype=np.int32), compression)
            if "in_transit" in old_group:
                _dataset(out_group, "in_transit", np.asarray(old_group["in_transit"], dtype=np.bool_), compression)
            _dataset(out_group, "transit_model", model.astype(np.float32), compression)

            injected_by_aperture: dict[str, np.ndarray] = {}
            baseline_by_aperture: dict[str, float] = {}
            for ap in apertures:
                original = np.asarray(target[ap], dtype=np.float64)
                good = (quality == 0) & np.isfinite(original)
                baseline = float(np.nanmedian(original[good])) if np.any(good) else float(np.nanmedian(original[np.isfinite(original)]))
                if not np.isfinite(baseline):
                    baseline = 1.0
                injected = original + baseline * (model - 1.0)
                injected_by_aperture[ap] = injected
                baseline_by_aperture[ap] = baseline
                out_group.attrs[f"baseline_{ap}"] = float(baseline)
                _dataset(out_group, f"{ap}_original", original.astype(np.float32), compression)
                _dataset(out_group, f"{ap}_injected", injected.astype(np.float32), compression)

            primary = apertures[0]
            _dataset(out_group, "flux_original", np.asarray(target[primary], dtype=np.float32), compression)
            _dataset(out_group, "flux_injected", injected_by_aperture[primary].astype(np.float32), compression)
            rows.append({"h5_group": group_path, "tic": tic, "status": "ok"})

    tmp_h5.replace(out_h5)
    status_counts = pd.Series([r["status"] for r in rows], dtype="object").value_counts().sort_index().to_dict()
    summary = {
        "created_utc": created_utc,
        "queue_csv": str(queue_csv),
        "old_injection_h5": str(old_injection_h5),
        "compact_export_h5": str(compact_export_h5),
        "out_h5": str(out_h5),
        "apertures": list(apertures),
        "n_rows": int(len(rows)),
        "status_counts": {str(k): int(v) for k, v in status_counts.items()},
        "missing_target": int(missing_target),
        "missing_aperture": int(missing_aperture),
        "time_mismatch": int(time_mismatch),
    }
    summary_path = out_h5.with_suffix(".summary.json")
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    rows_path = out_h5.with_suffix(".rows.csv")
    pd.DataFrame(rows).to_csv(rows_path, index=False)
    return summary


def _build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--queue-csv", type=Path, required=True)
    ap.add_argument("--old-injection-h5", type=Path, required=True)
    ap.add_argument("--compact-export-h5", type=Path, required=True)
    ap.add_argument("--out-h5", type=Path, required=True)
    ap.add_argument("--apertures", nargs=2, default=["DET_FLUX_ADP_SML", "DET_FLUX_ADP"])
    ap.add_argument("--compression", default="gzip")
    ap.add_argument("--overwrite", action="store_true")
    return ap


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    summary = build_subset(
        queue_csv=args.queue_csv,
        old_injection_h5=args.old_injection_h5,
        compact_export_h5=args.compact_export_h5,
        out_h5=args.out_h5,
        apertures=tuple(args.apertures),
        compression=args.compression or None,
        overwrite=bool(args.overwrite),
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
