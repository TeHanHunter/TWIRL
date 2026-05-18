#!/usr/bin/env python3
"""Scan /pdo/qlp-data/tica-delivery/s{NNNN}/ for TICA FFI filenames and emit a
JSON file in the schema that detrend_wrapper.py and hlsp_wrapper.py expect via
the TWIRL_SECTOR_ORBIT_JSON env var.

TICA FFI filename pattern (from tglc.utils.manifest):
    hlsp_tica_tess_ffi_s{sector:04d}-o{orbit_within_sector}-{cadence:08d}-cam{C}-ccd{D}_tess_v01_img.fits

We extract sector / orbit-within-sector / cadence from each filename and group
into per-(sector,orbit_within_sector) cadence ranges. The absolute orbit_number
in the OrbitalInfo schema is the orbit-within-sector tag mapped to the survey's
running orbit count; for sectors >= 56 the convention TWIRL has been using is
contiguous orbit numbers per sector (e.g. S56 -> 119, 120; S57 -> 121, 122).

The mapping orbit-within-sector -> absolute orbit_number is provided either
via --base-orbit (one orbit per (sector, ow) pair, ow=1 maps to base, ow=2 to
base+1, ...) or via --orbit-map (explicit "sector:ow:absolute" entries; useful
for sectors with non-contiguous orbits).

Output JSON schema (keyed by sector as string for dict round-trip in JSON):
    {
        "57": [
            {"orbit_number": 121, "cadence_range": [702400, 708100], "mjd_range": [0.0, 0.0]},
            {"orbit_number": 122, "cadence_range": [708200, 714000], "mjd_range": [0.0, 0.0]}
        ],
        ...
    }

Usage:
    python build_sector_orbit_json.py --sector 57 58 59 \\
        --base-orbit 57:121 \\
        --output /pdo/users/tehan/tglc-gpu-production/sector_orbit.json

The base-orbit list lets you anchor each sector's orbit numbering. For
sectors where you've already produced data, you can also pass
--inspect-only (just report what was found, don't write).
"""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

DEFAULT_TICA_ROOT = Path("/pdo/qlp-data/tica-delivery")
FNAME_RE = re.compile(
    r"hlsp_tica_tess_ffi_s(?P<sector>\d{4})-o(?P<ow>\d+)-(?P<cad>\d{8})-cam(?P<cam>\d)-ccd(?P<ccd>\d)_tess_v\d+_img\.fits"
)


def _parse_base_orbit(values: list[str]) -> dict[int, int]:
    out: dict[int, int] = {}
    for v in values or []:
        if ":" not in v:
            raise SystemExit(f"--base-orbit entry must look like S:N, got {v!r}")
        s, n = v.split(":", 1)
        out[int(s)] = int(n)
    return out


def _parse_orbit_map(values: list[str]) -> dict[tuple[int, int], int]:
    out: dict[tuple[int, int], int] = {}
    for v in values or []:
        parts = v.split(":")
        if len(parts) != 3:
            raise SystemExit(f"--orbit-map entry must look like S:OW:ABS, got {v!r}")
        out[(int(parts[0]), int(parts[1]))] = int(parts[2])
    return out


def _scan_sector(sector_dir: Path) -> dict[int, tuple[int, int]]:
    """Return {orbit_within_sector: (cad_min, cad_max)} from one sector dir."""
    cads_by_ow: dict[int, list[int]] = defaultdict(list)
    for ccd_dir in sector_dir.iterdir():
        if not ccd_dir.is_dir():
            continue
        for fits in ccd_dir.iterdir():
            m = FNAME_RE.match(fits.name)
            if not m:
                continue
            cads_by_ow[int(m["ow"])].append(int(m["cad"]))
    return {ow: (min(c), max(c)) for ow, c in cads_by_ow.items()}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--sector", type=int, nargs="+", required=True,
        help="Sectors to scan (e.g. 57 58 59)."
    )
    ap.add_argument(
        "--tica-root", type=Path, default=DEFAULT_TICA_ROOT,
        help="Root of TICA delivery tree."
    )
    ap.add_argument(
        "--base-orbit", action="append", default=[],
        help="Map sector to first absolute orbit (e.g. 57:121). Orbit ow=1 "
             "becomes base, ow=2 becomes base+1, etc."
    )
    ap.add_argument(
        "--orbit-map", action="append", default=[],
        help='Explicit "sector:orbit_within_sector:absolute_orbit" entries. '
             "Use for sectors with non-contiguous orbits. Overrides --base-orbit."
    )
    ap.add_argument("--output", type=Path, help="Write JSON here. If omitted, prints to stdout.")
    ap.add_argument("--inspect-only", action="store_true", help="Only report findings, do not write.")
    args = ap.parse_args()

    base = _parse_base_orbit(args.base_orbit)
    explicit = _parse_orbit_map(args.orbit_map)

    out: dict[str, list[dict[str, object]]] = {}
    for sector in args.sector:
        sector_dir = args.tica_root / f"s{sector:04d}"
        if not sector_dir.is_dir():
            print(f"[build_sector_orbit_json] sector {sector}: NO DIR {sector_dir}, skipping")
            continue
        ranges = _scan_sector(sector_dir)
        if not ranges:
            print(f"[build_sector_orbit_json] sector {sector}: no FFIs matched")
            continue
        entries = []
        for ow in sorted(ranges):
            cad_min, cad_max = ranges[ow]
            if (sector, ow) in explicit:
                abs_orbit = explicit[(sector, ow)]
            elif sector in base:
                abs_orbit = base[sector] + (ow - 1)
            else:
                raise SystemExit(
                    f"sector {sector}: missing --base-orbit or --orbit-map for ow={ow}; "
                    f"saw cadence range [{cad_min}, {cad_max}]"
                )
            entries.append({
                "orbit_number": abs_orbit,
                "cadence_range": [cad_min, cad_max],
                "mjd_range": [0.0, 0.0],
            })
            print(
                f"[build_sector_orbit_json] sector {sector} ow{ow} -> orbit {abs_orbit} "
                f"cadence_range=[{cad_min}, {cad_max}]"
            )
        out[str(sector)] = entries

    if args.inspect_only or args.output is None:
        print(json.dumps(out, indent=2))
        return 0

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"[build_sector_orbit_json] wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
