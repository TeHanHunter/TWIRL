#!/usr/bin/env python3
"""Freeze the label-free S56--S62 full-pool SSL observation registry."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from twirl.vetting.teacher_native_registry import (  # noqa: E402
    file_sha256,
    read_table,
)
from twirl.vetting.ssl_full_pool_eligibility import (  # noqa: E402
    load_native_model_eligibility,
)
from twirl.vetting.ssl_full_pool_native import (  # noqa: E402
    load_full_pool_native_registry_release,
)
from twirl.vetting.teacher_ssl_fullpool import (  # noqa: E402
    FULLPOOL_SSL_ANCHOR_APERTURE,
    FULLPOOL_SSL_SECTORS,
    build_fullpool_ssl_registry,
    load_frozen_ssl_full_pool,
    load_global_full_pool_bls,
    read_tic_inventory,
    write_fullpool_ssl_registry,
)


def _sectors(value: str) -> tuple[int, ...]:
    try:
        sectors = tuple(sorted({int(item) for item in value.split(",")}))
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "sectors must be comma-separated integers"
        ) from exc
    if not sectors or any(sector <= 0 for sector in sectors):
        raise argparse.ArgumentTypeError("sectors must be positive")
    return sectors


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--frozen-pool",
        type=Path,
        required=True,
        help="Preregistered retained S56--S62 observation table.",
    )
    parser.add_argument(
        "--frozen-pool-summary",
        type=Path,
        required=True,
        help="Checksum-bound summary for --frozen-pool and its allowlists.",
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
            "Optional explicitly staged copy of the summary-authorized BLS "
            "Parquet; byte count and SHA-256 must match exactly."
        ),
    )
    parser.add_argument(
        "--native-registry",
        type=Path,
        required=True,
        help="Separate unlabeled observation-keyed native-input registry.",
    )
    parser.add_argument(
        "--native-registry-summary",
        type=Path,
        required=True,
        help="Generic checksum-bound native-input registry summary.",
    )
    parser.add_argument(
        "--native-release-summary",
        type=Path,
        required=True,
        help="Full-pool native-v2 eligibility/coverage release summary.",
    )
    parser.add_argument(
        "--eligibility-exclusions",
        type=Path,
        required=True,
        help="Frozen native/model exclusions CSV.",
    )
    parser.add_argument(
        "--eligibility-summary",
        type=Path,
        required=True,
        help="Checksum-bound native/model eligibility summary.",
    )
    parser.add_argument(
        "--frozen-split-registry",
        type=Path,
        required=True,
        help="Frozen Teacher-v3 TIC split registry used only for host exclusion/folds.",
    )
    parser.add_argument(
        "--reserved-hosts",
        type=Path,
        required=True,
        help=(
            "Accepted prospective S63 host inventory: a strict sorted "
            "one-TIC-per-line .txt file, or a table with a tic column."
        ),
    )
    parser.add_argument("--registry-out", type=Path, required=True)
    parser.add_argument("--summary-out", type=Path, required=True)
    parser.add_argument(
        "--sectors",
        type=_sectors,
        default=FULLPOOL_SSL_SECTORS,
    )
    parser.add_argument(
        "--anchor-aperture",
        default=FULLPOOL_SSL_ANCHOR_APERTURE,
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    inputs = {
        "frozen_pool": args.frozen_pool.expanduser().resolve(),
        "frozen_pool_summary": (
            args.frozen_pool_summary.expanduser().resolve()
        ),
        "bls_summary": args.bls_summary.expanduser().resolve(),
        "native_registry": args.native_registry.expanduser().resolve(),
        "native_registry_summary": (
            args.native_registry_summary.expanduser().resolve()
        ),
        "native_release_summary": (
            args.native_release_summary.expanduser().resolve()
        ),
        "eligibility_exclusions": (
            args.eligibility_exclusions.expanduser().resolve()
        ),
        "eligibility_summary": (
            args.eligibility_summary.expanduser().resolve()
        ),
        "frozen_split_registry": (
            args.frozen_split_registry.expanduser().resolve()
        ),
        "reserved_hosts": args.reserved_hosts.expanduser().resolve(),
    }
    bls_override = (
        args.bls_peaks.expanduser().resolve()
        if args.bls_peaks is not None
        else None
    )
    if bls_override is not None:
        inputs["bls_peaks_override"] = bls_override
    for name, path in inputs.items():
        if not path.is_file() or path.stat().st_size <= 0:
            raise FileNotFoundError(f"{name} is missing or empty: {path}")

    frozen_pool, frozen_pool_summary = load_frozen_ssl_full_pool(
        pool_path=inputs["frozen_pool"],
        summary_path=inputs["frozen_pool_summary"],
        validate_allowlists=True,
    )
    split_sha256 = file_sha256(inputs["frozen_split_registry"])
    reserved_sha256 = file_sha256(inputs["reserved_hosts"])
    declared_inputs = frozen_pool_summary.get("inputs", {})
    declared_split_sha256 = declared_inputs.get(
        "tic_split_registry", {}
    ).get("sha256")
    declared_reserved_sha256 = declared_inputs.get(
        "s63_reserved_tics", {}
    ).get("sha256")
    if declared_split_sha256 != split_sha256:
        raise ValueError(
            "--frozen-split-registry does not match the frozen pool summary"
        )
    if declared_reserved_sha256 != reserved_sha256:
        raise ValueError(
            "--reserved-hosts does not match the frozen pool summary"
        )
    bls_rows, bls_summary, bls_path = load_global_full_pool_bls(
        summary_path=inputs["bls_summary"],
        frozen_pool=frozen_pool,
        frozen_pool_summary_path=inputs["frozen_pool_summary"],
        output_path_override=bls_override,
    )
    eligibility = load_native_model_eligibility(
        inputs["eligibility_exclusions"],
        inputs["eligibility_summary"],
        pool_path=inputs["frozen_pool"],
        pool_summary_path=inputs["frozen_pool_summary"],
        bls_path=bls_path,
        bls_summary_path=inputs["bls_summary"],
        production_lock=True,
    )
    native_registry, native_release = (
        load_full_pool_native_registry_release(
            registry_path=inputs["native_registry"],
            registry_summary_path=inputs["native_registry_summary"],
            release_summary_path=inputs["native_release_summary"],
            eligibility=eligibility,
        )
    )
    registry, audit = build_fullpool_ssl_registry(
        bls_rows,
        native_registry,
        read_table(inputs["frozen_split_registry"]),
        read_tic_inventory(inputs["reserved_hosts"]),
        frozen_pool=frozen_pool,
        eligibility=eligibility,
        sectors=args.sectors,
        anchor_aperture=str(args.anchor_aperture),
    )
    source_provenance = {
        name: {
            "path": str(path),
            "sha256": file_sha256(path),
            "size_bytes": int(path.stat().st_size),
        }
        for name, path in inputs.items()
    }
    source_provenance["frozen_pool_authority_bindings"] = {
        "split_registry_sha256_equal": True,
        "reserved_hosts_sha256_equal": True,
        "split_registry_sha256": split_sha256,
        "reserved_hosts_sha256": reserved_sha256,
    }
    source_provenance["global_bls_authority_binding"] = {
        "summary_path": str(inputs["bls_summary"]),
        "summary_sha256": file_sha256(inputs["bls_summary"]),
        "artifact_path": str(bls_path),
        "artifact_sha256": file_sha256(bls_path),
        "artifact_matches_summary": True,
        "summary_schema_version": bls_summary["schema_version"],
        "summary_contract_version": bls_summary["contract_version"],
    }
    source_provenance["native_model_eligibility_binding"] = {
        "contract_version": eligibility.contract_version,
        "release_binding": eligibility.release_binding,
        "full_observations": eligibility.n_full,
        "eligible_observations": eligibility.n_eligible,
        "excluded_observations": eligibility.n_excluded,
        "full_observation_identity_sha256": (
            eligibility.full_observation_identity_sha256
        ),
        "eligible_observation_identity_sha256": (
            eligibility.eligible_observation_identity_sha256
        ),
        "excluded_observation_identity_sha256": (
            eligibility.excluded_observation_identity_sha256
        ),
    }
    source_provenance["native_release_binding"] = {
        "schema_version": native_release["schema_version"],
        "release_binding": native_release["release_binding"],
        "native_contract_version": native_release[
            "native_contract_version"
        ],
        "release_summary_sha256": file_sha256(
            inputs["native_release_summary"]
        ),
    }
    result = write_fullpool_ssl_registry(
        registry,
        audit,
        registry_path=args.registry_out,
        summary_path=args.summary_out,
        source_provenance=source_provenance,
    )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
