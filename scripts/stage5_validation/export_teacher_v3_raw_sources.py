#!/usr/bin/env python3
"""Export one sector's compact raw/error host subset for Teacher v3."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import h5py


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.vetting.harmonic_export import export_tglc_raw_sources  # noqa: E402
from twirl.vetting.s63_preprocessing import (  # noqa: E402
    file_sha256,
    validate_producer_git_sha,
)

S63_PROSPECTIVE_CONTRACT = "s63_teacher_v3_prospective_v1"


def _authorize_sector(*, sector: int, orbits: list[int], contract: str | None) -> str:
    if 56 <= int(sector) <= 62:
        if contract is not None:
            raise ValueError("legacy Teacher-v3 raw export must not set a prospective contract")
        return "teacher_v3_s56_s62_legacy"
    if int(sector) == 63:
        if contract != S63_PROSPECTIVE_CONTRACT:
            raise ValueError(
                "S63 raw export requires --prospective-contract "
                f"{S63_PROSPECTIVE_CONTRACT}"
            )
        if sorted(int(value) for value in orbits) != [133, 134]:
            raise ValueError("prospective S63 raw export is locked to orbits 133/134")
        return S63_PROSPECTIVE_CONTRACT
    raise ValueError(
        "Teacher-v3 raw export is bounded to legacy S56-S62 plus explicitly "
        "authorized prospective S63"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sector", type=int, required=True)
    parser.add_argument("--training-table", type=Path, required=True)
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--compact-adp-h5", type=Path, required=True)
    parser.add_argument("--out-h5", type=Path, required=True)
    parser.add_argument("--orbits", type=int, nargs="+", required=True)
    parser.add_argument(
        "--prospective-contract",
        choices=(S63_PROSPECTIVE_CONTRACT,),
        help="Required only for the sealed prospective S63 inference lane.",
    )
    parser.add_argument("--producer-git-sha")
    args = parser.parse_args()
    if len(set(args.orbits)) != len(args.orbits) or any(
        orbit <= 0 for orbit in args.orbits
    ):
        raise SystemExit("--orbits must contain distinct positive orbit IDs")
    try:
        authorization = _authorize_sector(
            sector=args.sector,
            orbits=args.orbits,
            contract=args.prospective_contract,
        )
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    producer_git_sha: str | None = None
    if args.sector == 63:
        if args.producer_git_sha is None:
            raise SystemExit("prospective S63 raw export requires --producer-git-sha")
        producer_git_sha = validate_producer_git_sha(args.producer_git_sha)
    elif args.producer_git_sha is not None:
        raise SystemExit("legacy Teacher-v3 raw export must not set --producer-git-sha")
    training_table_sha256 = file_sha256(args.training_table)
    compact_h5_sha256 = file_sha256(args.compact_adp_h5)
    summary = export_tglc_raw_sources(
        training_table=args.training_table,
        raw_root=args.raw_root,
        out_h5=args.out_h5,
        orbits=tuple(args.orbits),
        compact_adp_h5=args.compact_adp_h5,
    )
    summary["sector"] = int(args.sector)
    summary["orbits"] = [int(value) for value in args.orbits]
    summary["authorization_contract"] = authorization
    if producer_git_sha is not None:
        if file_sha256(args.training_table) != training_table_sha256 or file_sha256(
            args.compact_adp_h5
        ) != compact_h5_sha256:
            raise RuntimeError("S63 raw-export inputs changed during production")
        with h5py.File(args.out_h5, "r+") as h5:
            h5.attrs["producer_git_sha"] = producer_git_sha
            h5.attrs["training_table_sha256"] = training_table_sha256
            h5.attrs["compact_adp_h5_sha256"] = compact_h5_sha256
            h5.flush()
        summary["producer_git_sha"] = producer_git_sha
        summary["training_table_sha256"] = training_table_sha256
        summary["compact_adp_h5_sha256"] = compact_h5_sha256
        summary["out_h5_sha256"] = file_sha256(args.out_h5)
    summary_path = args.out_h5.with_suffix(".summary.json")
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
