#!/usr/bin/env python3
"""Freeze the host-only S63 reservation used by full-pool SSL.

This reads target identities from the authoritative TWIRL observation table.
It does not discover, open, or otherwise inspect any S63 light curve.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
from astropy.table import Table


SCHEMA_VERSION = "twirl_teacher_ssl_reserved_sector_tics_v1"
DEFAULT_SECTOR = 63
DEFAULT_ORBITS = (133, 134)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def build_reserved_tics(
    *,
    observations: Path,
    out_tics: Path,
    sector: int = DEFAULT_SECTOR,
    orbits: tuple[int, ...] = DEFAULT_ORBITS,
) -> dict[str, Any]:
    """Write a deterministic one-TIC-per-line whole-host reservation."""

    observations = Path(observations).resolve()
    out_tics = Path(out_tics).resolve()
    if int(sector) != DEFAULT_SECTOR:
        raise ValueError("Teacher SSL prospective reservation is bounded to S63")
    normalized_orbits = tuple(sorted({int(value) for value in orbits}))
    if normalized_orbits != DEFAULT_ORBITS:
        raise ValueError(
            f"S63 reservation requires exact orbits {DEFAULT_ORBITS}"
        )
    if not observations.is_file():
        raise FileNotFoundError(observations)

    input_sha256 = _sha256(observations)
    table = Table.read(observations, memmap=True)
    required = {"sector", "orbit", "tic_id"}
    missing = sorted(required - set(table.colnames))
    if missing:
        raise KeyError(f"observation table lacks columns: {missing}")
    mask = (
        np.asarray(table["sector"], dtype=np.int64) == int(sector)
    ) & np.isin(
        np.asarray(table["orbit"], dtype=np.int64),
        np.asarray(normalized_orbits, dtype=np.int64),
    )
    tic_values = np.asarray(table["tic_id"], dtype=np.int64)
    mask &= tic_values > 0
    selected = tic_values[mask]
    tics = np.unique(selected)
    if tics.size == 0:
        raise ValueError("S63 reservation contains no positive TIC IDs")

    payload = "".join(f"{int(tic)}\n" for tic in tics)
    out_tics.parent.mkdir(parents=True, exist_ok=True)
    temporary = out_tics.with_suffix(out_tics.suffix + ".tmp")
    temporary.write_text(payload)
    temporary.replace(out_tics)
    output_sha256 = _sha256(out_tics)
    summary = {
        "schema_version": SCHEMA_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "identity_only": True,
        "light_curves_opened": False,
        "sector": int(sector),
        "orbits": list(normalized_orbits),
        "selection": (
            "sector == 63; orbit in {133,134}; tic_id > 0; unique TIC"
        ),
        "observations": str(observations),
        "observations_sha256": input_sha256,
        "n_selected_observation_rows": int(selected.size),
        "n_reserved_tics": int(tics.size),
        "reserved_tics": str(out_tics),
        "reserved_tics_sha256": output_sha256,
    }
    summary_path = out_tics.with_suffix(".summary.json")
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--observations", type=Path, required=True)
    parser.add_argument("--out-tics", type=Path, required=True)
    parser.add_argument("--sector", type=int, default=DEFAULT_SECTOR)
    parser.add_argument(
        "--orbits",
        type=int,
        nargs="+",
        default=list(DEFAULT_ORBITS),
    )
    args = parser.parse_args()
    summary = build_reserved_tics(
        observations=args.observations,
        out_tics=args.out_tics,
        sector=args.sector,
        orbits=tuple(args.orbits),
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
