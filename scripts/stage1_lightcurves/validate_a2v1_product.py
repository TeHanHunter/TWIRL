#!/usr/bin/env python3
"""Validate an A2v1 Stage 1 product tree.

Checks two product contracts:

1. TGLC HDF5 coverage under ``orbit-*/ffi/cam*/ccd*/LC/*.h5`` against the
   TWIRL observation table for the requested sector/orbits.
2. A2v1 HLSP FITS schema: ADP and ADP015 branch columns are present, while the
   canonical/default ``DET_FLUX*`` and ``SYS_RM_FLUX`` columns are absent.

The script is intentionally read-only and writes a compact JSON report.
"""
from __future__ import annotations

import argparse
import json
import re
from collections import Counter
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
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

    expected: list[ExpectedH5] = []
    seen: set[tuple[int, int, int, int]] = set()
    for row in table[mask]:
        item = ExpectedH5(
            orbit=int(row["orbit"]),
            camera=int(row["camera"]),
            ccd=int(row["ccd"]),
            tic=int(row["tic_id"]),
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
    return {
        "n_expected_h5": len(expected),
        "n_expected_unique_tics": len(unique_tics),
        "expected_h5_by_orbit": dict(sorted(by_orbit.items())),
        "expected_h5_by_ccd": dict(sorted(by_ccd.items())),
    }


def _audit_h5(
    root: Path,
    expected: list[ExpectedH5],
    *,
    max_examples: int,
) -> dict[str, object]:
    missing: list[ExpectedH5] = []
    present_by_orbit: Counter[str] = Counter()
    present_by_ccd: Counter[str] = Counter()
    zero_byte: list[str] = []

    for item in expected:
        path = root / item.rel_path
        if not path.exists():
            missing.append(item)
            continue
        present_by_orbit[str(item.orbit)] += 1
        present_by_ccd[item.ccd_key] += 1
        try:
            if path.stat().st_size == 0:
                zero_byte.append(str(path))
        except OSError:
            zero_byte.append(str(path))

    missing_by_orbit = Counter(str(item.orbit) for item in missing)
    missing_by_ccd = Counter(item.ccd_key for item in missing)
    return {
        "n_present_h5": len(expected) - len(missing),
        "n_missing_h5": len(missing),
        "n_zero_byte_h5": len(zero_byte),
        "present_h5_by_orbit": dict(sorted(present_by_orbit.items())),
        "present_h5_by_ccd": dict(sorted(present_by_ccd.items())),
        "missing_h5_by_orbit": dict(sorted(missing_by_orbit.items())),
        "missing_h5_by_ccd": dict(sorted(missing_by_ccd.items())),
        "missing_h5_examples": [asdict(item) | {"path": str(item.rel_path)} for item in missing[:max_examples]],
        "zero_byte_h5_examples": zero_byte[:max_examples],
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


def _audit_one_fits(path: Path) -> dict[str, object]:
    row: dict[str, object] = {"path": str(path), "ok": False}
    try:
        with fits.open(path, memmap=False) as hdul:
            header = hdul[0].header
            data = hdul["LIGHTCURVE"].data
            columns = set(hdul["LIGHTCURVE"].columns.names)
            missing_required = [
                col for col in REQUIRED_BASE_COLUMNS + REQUIRED_A2V1_COLUMNS if col not in columns
            ]
            forbidden_present = [col for col in FORBIDDEN_COLUMNS if col in columns]
            n_cadences = len(data)
            q0 = np.asarray(data["QUALITY"] == 0) if "QUALITY" in columns else np.zeros(n_cadences, dtype=bool)
            finite_counts = {}
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
                n_quality0=int(q0.sum()),
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
    expected_unique_tics: set[int],
    max_open: int,
    max_examples: int,
) -> dict[str, object]:
    fits_paths = _discover_fits(hlsp_root, sector)
    found_tics = {
        tic for path in fits_paths if (tic := _tic_from_fits_name(path, sector)) is not None
    }
    missing_tics = sorted(expected_unique_tics - found_tics)
    extra_tics = sorted(found_tics - expected_unique_tics)

    if max_open == 0:
        paths_to_open = fits_paths
    else:
        paths_to_open = fits_paths[:max_open]

    checked = [_audit_one_fits(path) for path in paths_to_open]
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
        "n_extra_fits_tics": len(extra_tics),
        "n_checked_fits": len(checked),
        "n_bad_checked_fits": len(bad),
        "checked_methods": methods,
        "checked_prodtags": prodtags,
        "checked_branches": branches,
        "missing_fits_tic_examples": missing_tics[:max_examples],
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
        "--allow-missing-h5",
        action="store_true",
        help="Return success even when requested HDF5 files are missing.",
    )
    parser.add_argument(
        "--allow-missing-fits",
        action="store_true",
        help="Return success even when expected target FITS files are missing.",
    )
    args = parser.parse_args(argv)

    a2v1_root = args.a2v1_root.expanduser()
    hlsp_root = args.hlsp_root or (a2v1_root / f"hlsp_s{args.sector:04d}_A2v1")
    expected = _read_expected_h5(
        args.observations.expanduser(),
        sector=args.sector,
        orbits=set(args.orbits),
    )
    expected_unique_tics = {item.tic for item in expected}

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
        },
    }
    summary.update(_summarize_expected(expected))
    summary["h5"] = _audit_h5(a2v1_root, expected, max_examples=args.max_examples)
    summary["fits"] = _audit_fits(
        hlsp_root,
        sector=args.sector,
        expected_unique_tics=expected_unique_tics,
        max_open=args.max_open_fits,
        max_examples=args.max_examples,
    )

    h5_ok = bool(args.allow_missing_h5 or summary["h5"]["n_missing_h5"] == 0)
    fits_ok = bool(
        args.allow_missing_fits
        or (
            summary["fits"]["n_missing_fits_tics"] == 0
            and summary["fits"]["n_bad_checked_fits"] == 0
        )
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
