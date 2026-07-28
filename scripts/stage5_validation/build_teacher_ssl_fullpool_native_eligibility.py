#!/usr/bin/env python3
"""Freeze the BLS-derived full-pool native/model eligibility authority."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from twirl.vetting.ssl_full_pool_eligibility import (  # noqa: E402
    NATIVE_MODEL_ELIGIBILITY_ANCHOR_APERTURE,
    write_native_model_eligibility,
)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--frozen-pool",
        type=Path,
        required=True,
        help="Frozen leakage-safe S56--S62 observation table.",
    )
    parser.add_argument(
        "--frozen-pool-summary",
        type=Path,
        required=True,
        help="Checksum-bound summary for --frozen-pool.",
    )
    parser.add_argument(
        "--bls-summary",
        type=Path,
        required=True,
        help="Passed checksum-bound global full-pool BLS summary.",
    )
    parser.add_argument(
        "--bls-peaks",
        type=Path,
        help=(
            "Optional byte-identical staged copy of the summary-authorized "
            "global BLS Parquet."
        ),
    )
    parser.add_argument(
        "--exclusions-out",
        type=Path,
        required=True,
        help="Immutable sorted CSV containing only native/model exclusions.",
    )
    parser.add_argument(
        "--summary-out",
        type=Path,
        required=True,
        help="Immutable checksum-bound eligibility summary JSON.",
    )
    parser.add_argument(
        "--anchor-aperture",
        default=NATIVE_MODEL_ELIGIBILITY_ANCHOR_APERTURE,
    )
    parser.add_argument(
        "--production-lock",
        action="store_true",
        help=(
            "Require the preregistered 175366/175347/19 production counts "
            "and exact full/eligible/excluded identity hashes."
        ),
    )
    return parser


def _authorized_bls_path(
    summary_path: Path,
    override: Path | None,
) -> Path:
    if override is not None:
        return override.expanduser().resolve()
    try:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        value = summary["output"]["path"]
    except (
        OSError,
        UnicodeDecodeError,
        json.JSONDecodeError,
        KeyError,
        TypeError,
    ) as exc:
        raise ValueError(
            "unable to resolve global BLS artifact from --bls-summary"
        ) from exc
    if not isinstance(value, str) or not value.strip():
        raise ValueError("global BLS summary output path is blank")
    return Path(value).expanduser().resolve()


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    pool_path = args.frozen_pool.expanduser().resolve()
    pool_summary_path = args.frozen_pool_summary.expanduser().resolve()
    bls_summary_path = args.bls_summary.expanduser().resolve()
    bls_path = _authorized_bls_path(bls_summary_path, args.bls_peaks)
    for name, path in (
        ("frozen pool", pool_path),
        ("frozen pool summary", pool_summary_path),
        ("global BLS", bls_path),
        ("global BLS summary", bls_summary_path),
    ):
        if not path.is_file() or path.stat().st_size <= 0:
            raise FileNotFoundError(f"{name} is missing or empty: {path}")
    authority = write_native_model_eligibility(
        pool_path=pool_path,
        pool_summary_path=pool_summary_path,
        bls_path=bls_path,
        bls_summary_path=bls_summary_path,
        exclusions_path=args.exclusions_out,
        summary_path=args.summary_out,
        anchor_aperture=str(args.anchor_aperture),
        production_lock=bool(args.production_lock),
    )
    print(
        json.dumps(
            {
                "passed": True,
                "contract_version": authority.contract_version,
                "release_binding": authority.release_binding,
                "n_full": authority.n_full,
                "n_eligible": authority.n_eligible,
                "n_excluded": authority.n_excluded,
                "full_observation_identity_sha256": (
                    authority.full_observation_identity_sha256
                ),
                "eligible_observation_identity_sha256": (
                    authority.eligible_observation_identity_sha256
                ),
                "excluded_observation_identity_sha256": (
                    authority.excluded_observation_identity_sha256
                ),
                "exclusions_path": str(authority.bindings["exclusions"].path),
                "summary_path": str(authority.bindings["eligibility_summary"].path),
            },
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
