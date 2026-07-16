#!/usr/bin/env python3
"""Validate an A2v1 Stage 1 product tree.

Checks two product contracts:

1. TGLC HDF5 coverage under ``orbit-*/ffi/cam*/ccd*/LC/*.h5`` against the
   TWIRL observation table for the requested sector/orbits, with an optional
   file-openability audit.
2. A2v1 HLSP FITS schema: ADP and ADP015 branch columns are present, while the
   canonical/default ``DET_FLUX*`` and ``SYS_RM_FLUX`` columns are absent.

The script is intentionally read-only and writes a compact JSON report.
"""
from __future__ import annotations

import argparse
import json
import re
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
import h5py
from astropy.io import fits
from astropy.table import Table


DEFAULT_A2V1_ROOT = Path("/pdo/users/tehan/tglc-gpu-production-A2v1")
DEFAULT_OBSERVATIONS = Path(
    "data_local/catalogs/twirl_master_catalog/twirl_wd_tess_observations_v0.fits"
)

REQUIRED_BASE_COLUMNS = (
    "TIME",
    "CADENCENO",
    "SAP_FLUX",
    "QUALITY",
    "ORBITID",
    "SAP_X",
    "SAP_Y",
    "SAP_BKG",
    "SAP_BKG_ERR",
)
REQUIRED_A2V1_COLUMNS = (
    "DET_FLUX_ADP_SML",
    "DET_FLUX_ADP",
    "DET_FLUX_ADP_LAG",
    "DET_FLUX_ADP_ERR",
    "DET_FLUX_ADP015_SML",
    "DET_FLUX_ADP015",
    "DET_FLUX_ADP015_LAG",
    "DET_FLUX_ADP015_ERR",
)
FORBIDDEN_COLUMNS = (
    "DET_FLUX",
    "DET_FLUX_ERR",
    "DET_FLUX_SML",
    "DET_FLUX_LAG",
    "SYS_RM_FLUX",
)
FITS_RE = re.compile(r"hlsp_twirlfs_tess_ffi_s(?P<sector>\d{4})-(?P<tic>\d{16})_")


@dataclass(frozen=True)
class ExpectedH5:
    orbit: int
    camera: int
    ccd: int
    tic: int
    edge_warn: bool = False

    @property
    def rel_path(self) -> Path:
        return (
            Path(f"orbit-{self.orbit}")
            / "ffi"
            / f"cam{self.camera}"
            / f"ccd{self.ccd}"
            / "LC"
            / f"{self.tic}.h5"
        )

    @property
    def ccd_key(self) -> str:
        return f"o{self.orbit}_cam{self.camera}_ccd{self.ccd}"


def _json_default(value: object) -> object:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def _read_expected_h5(
    observations_path: Path,
    *,
    sector: int,
    orbits: set[int],
) -> list[ExpectedH5]:
    table = Table.read(observations_path, memmap=True)
    mask = (table["sector"] == sector) & np.isin(table["orbit"], sorted(orbits))
    mask &= table["tic_id"] > 0
    has_edge_warn = "edge_warn" in table.colnames

    expected: list[ExpectedH5] = []
    seen: set[tuple[int, int, int, int]] = set()
    for row in table[mask]:
        item = ExpectedH5(
            orbit=int(row["orbit"]),
            camera=int(row["camera"]),
            ccd=int(row["ccd"]),
            tic=int(row["tic_id"]),
            edge_warn=bool(row["edge_warn"]) if has_edge_warn else False,
        )
        key = (item.orbit, item.camera, item.ccd, item.tic)
        if key in seen:
            continue
        seen.add(key)
        expected.append(item)
    return expected


def _summarize_expected(expected: Iterable[ExpectedH5]) -> dict[str, object]:
    expected = list(expected)
    by_orbit = Counter(str(item.orbit) for item in expected)
    by_ccd = Counter(item.ccd_key for item in expected)
    unique_tics = {item.tic for item in expected}
    edge_warn_tics = {item.tic for item in expected if item.edge_warn}
    return {
        "n_expected_h5": len(expected),
        "n_expected_unique_tics": len(unique_tics),
        "n_expected_h5_edge_warn": sum(item.edge_warn for item in expected),
        "n_expected_unique_tics_with_any_edge_warn": len(edge_warn_tics),
        "expected_h5_by_orbit": dict(sorted(by_orbit.items())),
        "expected_h5_by_ccd": dict(sorted(by_ccd.items())),
    }


def _audit_h5(
    root: Path,
    expected: list[ExpectedH5],
    *,
    max_examples: int,
    check_open: bool,
) -> dict[str, object]:
    missing: list[ExpectedH5] = []
    present_by_orbit: Counter[str] = Counter()
    present_by_ccd: Counter[str] = Counter()
    zero_byte: list[str] = []
    unreadable: list[dict[str, str]] = []

    for item in expected:
        path = root / item.rel_path
        if not path.exists():
            missing.append(item)
            continue
        present_by_orbit[str(item.orbit)] += 1
        present_by_ccd[item.ccd_key] += 1
        try:
            is_zero_byte = path.stat().st_size == 0
        except OSError:
            is_zero_byte = True
        if is_zero_byte:
            zero_byte.append(str(path))
            continue
        if check_open:
            try:
                with h5py.File(path, "r") as handle:
                    # Accessing the top-level keys forces HDF5 to read the file
                    # signature and root object rather than merely opening a path.
                    tuple(handle.keys())
            except OSError as exc:
                unreadable.append(
                    {"path": str(path), "error": f"{type(exc).__name__}: {exc}"}
                )

    missing_by_orbit = Counter(str(item.orbit) for item in missing)
    missing_by_ccd = Counter(item.ccd_key for item in missing)
    missing_edge_warn = [item for item in missing if item.edge_warn]
    missing_non_edge = [item for item in missing if not item.edge_warn]
    return {
        "n_present_h5": len(expected) - len(missing),
        "n_missing_h5": len(missing),
        "n_zero_byte_h5": len(zero_byte),
        "n_unreadable_h5": len(unreadable),
        "present_h5_by_orbit": dict(sorted(present_by_orbit.items())),
        "present_h5_by_ccd": dict(sorted(present_by_ccd.items())),
        "missing_h5_by_orbit": dict(sorted(missing_by_orbit.items())),
        "missing_h5_by_ccd": dict(sorted(missing_by_ccd.items())),
        "n_missing_h5_edge_warn": len(missing_edge_warn),
        "n_missing_h5_non_edge": len(missing_non_edge),
        "missing_h5_examples": [asdict(item) | {"path": str(item.rel_path)} for item in missing[:max_examples]],
        "missing_h5_edge_warn_examples": [
            asdict(item) | {"path": str(item.rel_path)} for item in missing_edge_warn[:max_examples]
        ],
        "missing_h5_non_edge_examples": [
            asdict(item) | {"path": str(item.rel_path)} for item in missing_non_edge[:max_examples]
        ],
        "zero_byte_h5_examples": zero_byte[:max_examples],
        "unreadable_h5_examples": unreadable[:max_examples],
    }


def _discover_fits(hlsp_root: Path, sector: int) -> list[Path]:
    pattern = f"hlsp_twirlfs_tess_ffi_s{sector:04d}-*.fits"
    return sorted(hlsp_root.rglob(pattern))


def _tic_from_fits_name(path: Path, sector: int) -> int | None:
    match = FITS_RE.search(path.name)
    if match is None:
        return None
    if int(match.group("sector")) != sector:
        return None
    return int(match.group("tic"))


def _audit_one_fits(path: Path, *, check_finite_values: bool) -> dict[str, object]:
    row: dict[str, object] = {"path": str(path), "ok": False}
    try:
        with fits.open(path, memmap=False) as hdul:
            header = hdul[0].header
            lightcurve_hdu = hdul["LIGHTCURVE"]
            columns = set(lightcurve_hdu.columns.names)
            missing_required = [
                col for col in REQUIRED_BASE_COLUMNS + REQUIRED_A2V1_COLUMNS if col not in columns
            ]
            forbidden_present = [col for col in FORBIDDEN_COLUMNS if col in columns]
            n_cadences = int(lightcurve_hdu.header.get("NAXIS2", 0))
            n_quality0: int | None = None
            finite_counts = {}
            if check_finite_values:
                data = lightcurve_hdu.data
                q0 = (
                    np.asarray(data["QUALITY"] == 0)
                    if "QUALITY" in columns
                    else np.zeros(n_cadences, dtype=bool)
                )
                n_quality0 = int(q0.sum())
                for col in REQUIRED_A2V1_COLUMNS:
                    if col in columns:
                        values = np.asarray(data[col], dtype=float)
                        finite_counts[col] = int(np.isfinite(values[q0]).sum())
            row.update(
                ok=not missing_required
                and not forbidden_present
                and str(header.get("METHOD", "")) == "A2v1"
                and str(header.get("PRODTAG", "")) == "A2v1"
                and bool(header.get("A2V1", False)),
                tic=int(header.get("TICID", -1)),
                sector=int(header.get("SECTOR", -1)),
                method=str(header.get("METHOD", "")),
                prodtag=str(header.get("PRODTAG", "")),
                a2v1=bool(header.get("A2V1", False)),
                branches=str(header.get("BRANCHES", "")),
                n_cadences=int(n_cadences),
                n_quality0=n_quality0,
                missing_required_columns=missing_required,
                forbidden_columns_present=forbidden_present,
                finite_quality0_counts=finite_counts,
            )
    except Exception as exc:  # noqa: BLE001 - validator should report all bad files.
        row.update(error=str(exc))
    return row


def _audit_fits(
    hlsp_root: Path,
    *,
    sector: int,
    expected_tic_edge_warn: dict[int, bool],
    max_open: int,
    max_examples: int,
    check_finite_values: bool,
    fits_workers: int,
    progress_every: int,
) -> dict[str, object]:
    fits_paths = _discover_fits(hlsp_root, sector)
    found_tics = {
        tic for path in fits_paths if (tic := _tic_from_fits_name(path, sector)) is not None
    }
    expected_unique_tics = set(expected_tic_edge_warn)
    missing_tics = sorted(expected_unique_tics - found_tics)
    missing_edge_warn_tics = [tic for tic in missing_tics if expected_tic_edge_warn[tic]]
    missing_non_edge_tics = [tic for tic in missing_tics if not expected_tic_edge_warn[tic]]
    extra_tics = sorted(found_tics - expected_unique_tics)

    if max_open == 0:
        paths_to_open = fits_paths
    else:
        paths_to_open = fits_paths[:max_open]

    checked: list[dict[str, object]] = []
    total_to_check = len(paths_to_open)
    if fits_workers == 1:
        for index, path in enumerate(paths_to_open, start=1):
            checked.append(_audit_one_fits(path, check_finite_values=check_finite_values))
            if progress_every > 0 and (index % progress_every == 0 or index == total_to_check):
                print(f"[validate-a2v1] FITS schema {index:,}/{total_to_check:,}", flush=True)
    else:
        with ThreadPoolExecutor(max_workers=fits_workers) as executor:
            futures = [
                executor.submit(_audit_one_fits, path, check_finite_values=check_finite_values)
                for path in paths_to_open
            ]
            for index, future in enumerate(as_completed(futures), start=1):
                checked.append(future.result())
                if progress_every > 0 and (index % progress_every == 0 or index == total_to_check):
                    print(f"[validate-a2v1] FITS schema {index:,}/{total_to_check:,}", flush=True)
    checked.sort(key=lambda row: str(row.get("path", "")))
    bad = [row for row in checked if not row.get("ok", False)]
    methods = sorted({row.get("method", "") for row in checked if row.get("method")})
    prodtags = sorted({row.get("prodtag", "") for row in checked if row.get("prodtag")})
    branches = sorted({row.get("branches", "") for row in checked if row.get("branches")})

    return {
        "hlsp_root": str(hlsp_root),
        "n_fits": len(fits_paths),
        "n_found_unique_tics": len(found_tics),
        "n_expected_unique_tics": len(expected_unique_tics),
        "n_missing_fits_tics": len(missing_tics),
        "n_missing_fits_edge_warn_tics": len(missing_edge_warn_tics),
        "n_missing_fits_non_edge_tics": len(missing_non_edge_tics),
        "n_extra_fits_tics": len(extra_tics),
        "n_checked_fits": len(checked),
        "n_bad_checked_fits": len(bad),
        "checked_methods": methods,
        "checked_prodtags": prodtags,
        "checked_branches": branches,
        "missing_fits_tic_examples": missing_tics[:max_examples],
        "missing_fits_edge_warn_tic_examples": missing_edge_warn_tics[:max_examples],
        "missing_fits_non_edge_tic_examples": missing_non_edge_tics[:max_examples],
        "extra_fits_tic_examples": extra_tics[:max_examples],
        "bad_fits_examples": bad[:max_examples],
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--a2v1-root", type=Path, default=DEFAULT_A2V1_ROOT)
    parser.add_argument("--hlsp-root", type=Path)
    parser.add_argument("--observations", type=Path, default=DEFAULT_OBSERVATIONS)
    parser.add_argument("--sector", type=int, default=56)
    parser.add_argument("--orbits", type=int, nargs="+", default=[119, 120])
    parser.add_argument(
        "--max-open-fits",
        type=int,
        default=0,
        help="Maximum FITS files to open for schema checks; 0 means all.",
    )
    parser.add_argument("--max-examples", type=int, default=20)
    parser.add_argument("--summary-json", type=Path)
    parser.add_argument(
        "--skip-fits",
        action="store_true",
        help="Audit only the HDF5 checkpoint and do not discover or open FITS files.",
    )
    parser.add_argument(
        "--check-h5-open",
        action="store_true",
        help="Open every present, nonzero HDF5 file and fail unreadable files.",
    )
    parser.add_argument(
        "--allow-missing-h5",
        action="store_true",
        help="Return success even when requested HDF5 files are missing.",
    )
    parser.add_argument(
        "--allow-missing-fits",
        action="store_true",
        help="Return success even when expected target FITS files are missing.",
    )
    parser.add_argument(
        "--allow-edge-warn-missing",
        action="store_true",
        help=(
            "Treat missing targets flagged edge_warn in the observation table as expected "
            "TGLC aperture-edge exclusions. Non-edge gaps still fail validation."
        ),
    )
    parser.add_argument(
        "--schema-only",
        action="store_true",
        help="Validate FITS headers and column schema without reading cadence arrays.",
    )
    parser.add_argument(
        "--fits-workers",
        type=int,
        default=1,
        help="Concurrent FITS readers for schema checks (default: 1).",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=500,
        help="Print a FITS schema checkpoint every N files; 0 disables checkpoints.",
    )
    args = parser.parse_args(argv)
    if args.fits_workers < 1:
        parser.error("--fits-workers must be >= 1")
    if args.progress_every < 0:
        parser.error("--progress-every must be >= 0")

    a2v1_root = args.a2v1_root.expanduser()
    hlsp_root = args.hlsp_root or (a2v1_root / f"hlsp_s{args.sector:04d}_A2v1")
    expected = _read_expected_h5(
        args.observations.expanduser(),
        sector=args.sector,
        orbits=set(args.orbits),
    )
    expected_tic_edge_warn: dict[int, bool] = {}
    for item in expected:
        expected_tic_edge_warn[item.tic] = expected_tic_edge_warn.get(item.tic, True) and item.edge_warn

    requested_orbits = {int(orbit) for orbit in args.orbits}
    observed_orbits = {item.orbit for item in expected}
    missing_requested_orbits = sorted(requested_orbits - observed_orbits)
    expected_contract_ok = bool(expected_tic_edge_warn) and not missing_requested_orbits

    summary: dict[str, object] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "sector": int(args.sector),
        "orbits": [int(orbit) for orbit in args.orbits],
        "a2v1_root": str(a2v1_root),
        "observations": str(args.observations),
        "schema": {
            "required_base_columns": REQUIRED_BASE_COLUMNS,
            "required_a2v1_columns": REQUIRED_A2V1_COLUMNS,
            "forbidden_columns": FORBIDDEN_COLUMNS,
            "expected_method": "A2v1",
            "expected_prodtag": "A2v1",
            "schema_only": bool(args.schema_only),
            "check_h5_open": bool(args.check_h5_open),
            "allow_edge_warn_missing": bool(args.allow_edge_warn_missing),
            "fits_workers": int(args.fits_workers),
            "progress_every": int(args.progress_every),
        },
    }
    summary.update(_summarize_expected(expected))
    summary["expected_contract"] = {
        "ok": bool(expected_contract_ok),
        "requested_orbits": sorted(requested_orbits),
        "observed_orbits": sorted(observed_orbits),
        "missing_requested_orbits": missing_requested_orbits,
        "has_expected_rows": bool(expected),
        "has_expected_unique_tics": bool(expected_tic_edge_warn),
    }
    summary["h5"] = _audit_h5(
        a2v1_root,
        expected,
        max_examples=args.max_examples,
        check_open=args.check_h5_open,
    )
    if args.skip_fits:
        summary["fits"] = {"skipped": True}
    else:
        summary["fits"] = _audit_fits(
            hlsp_root,
            sector=args.sector,
            expected_tic_edge_warn=expected_tic_edge_warn,
            max_open=args.max_open_fits,
            max_examples=args.max_examples,
            check_finite_values=not args.schema_only,
            fits_workers=args.fits_workers,
            progress_every=args.progress_every,
        )

    h5_missing = (
        summary["h5"]["n_missing_h5_non_edge"]
        if args.allow_edge_warn_missing
        else summary["h5"]["n_missing_h5"]
    )
    h5_ok = bool(
        expected_contract_ok
        and (args.allow_missing_h5 or h5_missing == 0)
        and summary["h5"]["n_zero_byte_h5"] == 0
        and summary["h5"]["n_unreadable_h5"] == 0
    )
    if args.skip_fits:
        fits_ok = True
    else:
        fits_missing = (
            summary["fits"]["n_missing_fits_non_edge_tics"]
            if args.allow_edge_warn_missing
            else summary["fits"]["n_missing_fits_tics"]
        )
        fits_ok = bool(
            (args.allow_missing_fits or fits_missing == 0)
            and summary["fits"]["n_bad_checked_fits"] == 0
        )
    summary["ok"] = bool(h5_ok and fits_ok)
    summary["ok_h5"] = h5_ok
    summary["ok_fits"] = fits_ok

    text = json.dumps(summary, indent=2, sort_keys=True, default=_json_default)
    if args.summary_json:
        args.summary_json.parent.mkdir(parents=True, exist_ok=True)
        args.summary_json.write_text(text + "\n")
    print(text)
    return 0 if summary["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
