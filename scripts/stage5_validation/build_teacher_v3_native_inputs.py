#!/usr/bin/env python3
"""Build one quality-aware Teacher-v3 native-input HDF5 shard."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import h5py


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.vetting.harmonic_export import build_raw_pair_export  # noqa: E402
from twirl.lightcurves.external_quality import ORBITID_POLICIES  # noqa: E402
from twirl.vetting.harmonic_inputs import verify_raw_pair_contract  # noqa: E402
from twirl.vetting.s63_preprocessing import (  # noqa: E402
    file_sha256,
    validate_producer_git_sha,
)

S63_PROSPECTIVE_CONTRACT = "s63_teacher_v3_prospective_v1"


def _authorize_sector(*, sector: int, prospective_contract: str | None) -> str:
    """Keep the frozen S56--S62 lane unchanged and authorize only sealed S63."""

    if 56 <= int(sector) <= 62:
        if prospective_contract is not None:
            raise ValueError(
                "legacy Teacher-v3 sectors must not set --prospective-contract"
            )
        return "teacher_v3_s56_s62_legacy"
    if int(sector) == 63:
        if prospective_contract != S63_PROSPECTIVE_CONTRACT:
            raise ValueError(
                "S63 native input requires --prospective-contract "
                f"{S63_PROSPECTIVE_CONTRACT}"
            )
        return S63_PROSPECTIVE_CONTRACT
    raise ValueError(
        "Teacher-v3 native input is bounded to legacy S56-S62 plus "
        "explicitly authorized prospective S63"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sector", type=int, required=True)
    parser.add_argument("--training-table", type=Path, required=True)
    parser.add_argument("--raw-source-h5", type=Path, required=True)
    parser.add_argument("--compact-adp-h5", type=Path, required=True)
    parser.add_argument("--injection-pair-h5", type=Path)
    parser.add_argument("--cadence-reference-table", type=Path, required=True)
    parser.add_argument("--cadence-reference-manifest", type=Path, required=True)
    parser.add_argument("--out-h5", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    parser.add_argument("--n-periods", type=int, default=4096)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--n-shards", type=int, default=1)
    parser.add_argument(
        "--orbitid-policy",
        choices=ORBITID_POLICIES,
        default="strict",
    )
    parser.add_argument("--producer-git-sha")
    parser.add_argument(
        "--prospective-contract",
        choices=(S63_PROSPECTIVE_CONTRACT,),
        help="Required only for the sealed prospective S63 inference lane.",
    )
    args = parser.parse_args()
    try:
        authorization = _authorize_sector(
            sector=args.sector,
            prospective_contract=args.prospective_contract,
        )
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    if args.sector == 63 and args.orbitid_policy != "strict":
        raise SystemExit("prospective S63 native input requires strict orbit IDs")
    if args.sector == 63 and args.injection_pair_h5 is not None:
        raise SystemExit("prospective S63 native input is real-only")
    producer_git_sha: str | None = None
    if args.sector == 63:
        if args.producer_git_sha is None:
            raise SystemExit("prospective S63 native input requires --producer-git-sha")
        producer_git_sha = validate_producer_git_sha(args.producer_git_sha)
    elif args.producer_git_sha is not None:
        raise SystemExit("legacy Teacher-v3 native input must not set --producer-git-sha")
    build = build_raw_pair_export(
        training_table=args.training_table,
        raw_source_h5=args.raw_source_h5,
        compact_adp_h5=args.compact_adp_h5,
        injection_pair_h5=args.injection_pair_h5,
        cadence_reference_table=args.cadence_reference_table,
        cadence_reference_manifest=args.cadence_reference_manifest,
        out_h5=args.out_h5,
        repo_root=args.repo_root,
        sector=args.sector,
        n_periods=args.n_periods,
        shard_index=args.shard_index,
        n_shards=args.n_shards,
        orbitid_policy=args.orbitid_policy,
    )
    if producer_git_sha is not None:
        with h5py.File(args.out_h5, "r+") as h5:
            h5.attrs["producer_git_sha"] = producer_git_sha
            h5.flush()
    verification = verify_raw_pair_contract(
        args.out_h5,
        require_errors=True,
        require_periodograms=True,
    )
    summary = {
        "authorization_contract": authorization,
        "build": build,
        "verification": verification,
    }
    if producer_git_sha is not None:
        summary.update(
            {
                "sector": int(args.sector),
                "producer_git_sha": producer_git_sha,
                "out_h5_sha256": file_sha256(args.out_h5),
            }
        )
    args.out_h5.with_suffix(".summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))
    if not verification["passed"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
