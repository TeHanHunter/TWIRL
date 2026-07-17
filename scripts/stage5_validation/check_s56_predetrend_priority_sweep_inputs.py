#!/usr/bin/env python3
"""Preflight-check inputs for the S56 pre-detrend priority recovery sweep."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import pandas as pd


REQUIRED_PRIORITY_COLUMNS = {"injection_id", "priority_reason", "recommended_aperture", "failure_class"}
REQUIRED_SURVIVAL_COLUMNS = {"injection_id", "aperture", "delta_signal_multi_snr_mad"}
REQUIRED_GROUP_DATASETS = ("time", "quality", "in_transit")


def _json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    try:
        import numpy as np

        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, np.ndarray):
            return value.tolist()
    except Exception:
        pass
    return str(value)


def load_injection_ids(path: Path) -> tuple[str, ...]:
    ids: list[str] = []
    seen: set[str] = set()
    for line in path.read_text().splitlines():
        item = line.strip()
        if not item or item.startswith("#"):
            continue
        item = item.split(",", 1)[0].strip()
        if item and item not in seen:
            ids.append(item)
            seen.add(item)
    return tuple(ids)


def _parse_aperture_attr(value: Any) -> tuple[str, ...]:
    if value is None:
        return ()
    raw = str(value)
    if not raw:
        return ()
    try:
        decoded = json.loads(raw)
    except json.JSONDecodeError:
        return tuple(part for part in raw.split("|") if part)
    if isinstance(decoded, list):
        return tuple(str(part) for part in decoded if str(part))
    return (str(decoded),)


def _available_apertures(group: Any) -> set[str]:
    out = set(_parse_aperture_attr(group.attrs.get("apertures")))
    single = str(group.attrs.get("aperture", ""))
    if single:
        out.add(single)
    for name in group.keys():
        if name.endswith("_injected"):
            out.add(name[: -len("_injected")])
    return out


def _has_injected_flux(group: Any, aperture: str) -> bool:
    if f"{aperture}_injected" in group:
        return True
    return "flux_injected" in group and str(group.attrs.get("aperture", "")) == aperture


def _check_h5(
    injection_h5: Path,
    injection_ids: tuple[str, ...],
    required_apertures: tuple[str, ...],
) -> dict[str, Any]:
    import h5py

    summary: dict[str, Any] = {
        "path": str(injection_h5),
        "exists": injection_h5.exists(),
        "n_groups": 0,
        "n_requested_ids": len(injection_ids),
        "missing_ids": [],
        "missing_required_datasets": {},
        "missing_required_apertures": {},
        "length_mismatch_ids": [],
        "available_aperture_counts": {},
    }
    if not injection_h5.exists():
        return summary

    with h5py.File(injection_h5, "r") as h5:
        if "injections" not in h5:
            summary["missing_injections_group"] = True
            return summary
        root = h5["injections"]
        all_ids = set(root.keys())
        summary["n_groups"] = len(all_ids)
        missing_ids = [inj_id for inj_id in injection_ids if inj_id not in all_ids]
        summary["missing_ids"] = missing_ids
        for inj_id in injection_ids:
            if inj_id not in all_ids:
                continue
            group = root[inj_id]
            missing_datasets = [name for name in REQUIRED_GROUP_DATASETS if name not in group]
            if missing_datasets:
                summary["missing_required_datasets"][inj_id] = missing_datasets

            available = _available_apertures(group)
            for aperture in available:
                summary["available_aperture_counts"][aperture] = (
                    int(summary["available_aperture_counts"].get(aperture, 0)) + 1
                )
            missing_apertures = [ap for ap in required_apertures if not _has_injected_flux(group, ap)]
            if missing_apertures:
                summary["missing_required_apertures"][inj_id] = missing_apertures

            lengths = []
            for name in ("time", "quality", "in_transit"):
                if name in group:
                    lengths.append(tuple(group[name].shape))
            for aperture in required_apertures:
                dataset_name = f"{aperture}_injected"
                if dataset_name in group:
                    lengths.append(tuple(group[dataset_name].shape))
                elif "flux_injected" in group and str(group.attrs.get("aperture", "")) == aperture:
                    lengths.append(tuple(group["flux_injected"].shape))
            if lengths and len(set(lengths)) > 1:
                summary["length_mismatch_ids"].append(inj_id)
    return summary


def _check_priority_csv(priority_csv: Path, injection_ids: tuple[str, ...]) -> dict[str, Any]:
    summary: dict[str, Any] = {"path": str(priority_csv), "exists": priority_csv.exists()}
    if not priority_csv.exists():
        return summary
    df = pd.read_csv(priority_csv)
    missing_cols = sorted(REQUIRED_PRIORITY_COLUMNS - set(df.columns))
    ids_in_csv = set(df["injection_id"].dropna().astype(str)) if "injection_id" in df else set()
    summary.update(
        {
            "n_rows": int(len(df)),
            "missing_columns": missing_cols,
            "n_unique_ids": int(len(ids_in_csv)),
            "ids_missing_from_csv": [inj_id for inj_id in injection_ids if inj_id not in ids_in_csv],
            "extra_ids_in_csv": sorted(ids_in_csv - set(injection_ids)),
            "priority_reason_counts": df["priority_reason"].value_counts().to_dict()
            if "priority_reason" in df
            else {},
            "recommended_aperture_counts": df["recommended_aperture"].value_counts().to_dict()
            if "recommended_aperture" in df
            else {},
        }
    )
    return summary


def _check_survival_csv(
    survival_csv: Path,
    injection_ids: tuple[str, ...],
    required_apertures: tuple[str, ...],
) -> dict[str, Any]:
    summary: dict[str, Any] = {"path": str(survival_csv), "exists": survival_csv.exists()}
    if not survival_csv.exists():
        return summary
    df = pd.read_csv(survival_csv)
    missing_cols = sorted(REQUIRED_SURVIVAL_COLUMNS - set(df.columns))
    summary["n_rows"] = int(len(df))
    summary["missing_columns"] = missing_cols
    if missing_cols:
        return summary
    df["injection_id"] = df["injection_id"].astype(str)
    df["aperture"] = df["aperture"].astype(str)
    id_set = set(injection_ids)
    sub = df[df["injection_id"].isin(id_set)]
    summary.update(
        {
            "n_priority_rows": int(len(sub)),
            "n_priority_ids": int(sub["injection_id"].nunique()),
            "ids_missing_from_survival": [inj_id for inj_id in injection_ids if inj_id not in set(sub["injection_id"])],
            "aperture_counts_for_priority_ids": sub["aperture"].value_counts().to_dict(),
            "missing_survival_by_aperture": {},
        }
    )
    for aperture in required_apertures:
        ids_for_ap = set(sub.loc[sub["aperture"].eq(aperture), "injection_id"])
        missing = [inj_id for inj_id in injection_ids if inj_id not in ids_for_ap]
        if missing:
            summary["missing_survival_by_aperture"][aperture] = missing
    return summary


def preflight(
    *,
    injection_h5: Path,
    injection_id_file: Path,
    priority_queue_csv: Path,
    survival_csv: Path,
    required_apertures: tuple[str, ...],
    expected_id_count: int = 0,
) -> dict[str, Any]:
    injection_ids = load_injection_ids(injection_id_file)
    summary: dict[str, Any] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "injection_id_file": str(injection_id_file),
        "n_injection_ids": len(injection_ids),
        "expected_id_count": int(expected_id_count),
        "required_apertures": list(required_apertures),
        "id_count_matches_expected": bool(expected_id_count <= 0 or len(injection_ids) == expected_id_count),
        "h5": _check_h5(injection_h5, injection_ids, required_apertures),
        "priority_csv": _check_priority_csv(priority_queue_csv, injection_ids),
        "survival_csv": _check_survival_csv(survival_csv, injection_ids, required_apertures),
    }

    errors: list[str] = []
    if not injection_id_file.exists():
        errors.append(f"missing injection ID file: {injection_id_file}")
    if expected_id_count > 0 and len(injection_ids) != expected_id_count:
        errors.append(f"expected {expected_id_count} injection IDs, found {len(injection_ids)}")

    h5 = summary["h5"]
    if not h5.get("exists"):
        errors.append(f"missing injection HDF5: {injection_h5}")
    if h5.get("missing_injections_group"):
        errors.append("injection HDF5 is missing /injections group")
    if h5.get("missing_ids"):
        errors.append(f"HDF5 missing {len(h5['missing_ids'])} requested IDs")
    if h5.get("missing_required_datasets"):
        errors.append(f"HDF5 missing required datasets for {len(h5['missing_required_datasets'])} IDs")
    if h5.get("missing_required_apertures"):
        errors.append(f"HDF5 missing required apertures for {len(h5['missing_required_apertures'])} IDs")
    if h5.get("length_mismatch_ids"):
        errors.append(f"HDF5 has length mismatches for {len(h5['length_mismatch_ids'])} IDs")

    priority = summary["priority_csv"]
    if not priority.get("exists"):
        errors.append(f"missing priority queue CSV: {priority_queue_csv}")
    if priority.get("missing_columns"):
        errors.append(f"priority queue missing columns: {priority['missing_columns']}")
    if priority.get("ids_missing_from_csv"):
        errors.append(f"priority queue missing {len(priority['ids_missing_from_csv'])} requested IDs")

    survival = summary["survival_csv"]
    if not survival.get("exists"):
        errors.append(f"missing survival CSV: {survival_csv}")
    if survival.get("missing_columns"):
        errors.append(f"survival CSV missing columns: {survival['missing_columns']}")
    if survival.get("ids_missing_from_survival"):
        errors.append(f"survival CSV missing {len(survival['ids_missing_from_survival'])} requested IDs")
    if survival.get("missing_survival_by_aperture"):
        counts = {ap: len(ids) for ap, ids in survival["missing_survival_by_aperture"].items()}
        errors.append(f"survival CSV missing required aperture rows: {counts}")

    summary["ok"] = not errors
    summary["errors"] = errors
    return summary


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--injection-h5", type=Path, required=True)
    parser.add_argument("--injection-id-file", type=Path, required=True)
    parser.add_argument("--priority-queue-csv", type=Path, required=True)
    parser.add_argument("--survival-csv", type=Path, required=True)
    parser.add_argument("--required-aperture", action="append", default=[])
    parser.add_argument("--expected-id-count", type=int, default=0)
    parser.add_argument("--out-json", type=Path, default=None)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    summary = preflight(
        injection_h5=args.injection_h5,
        injection_id_file=args.injection_id_file,
        priority_queue_csv=args.priority_queue_csv,
        survival_csv=args.survival_csv,
        required_apertures=tuple(args.required_aperture),
        expected_id_count=args.expected_id_count,
    )
    text = json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n"
    if args.out_json:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(text)
    print(text)
    return 0 if summary["ok"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
