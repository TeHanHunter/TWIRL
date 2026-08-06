#!/usr/bin/env python3
"""Bind the final S63 compact ADP pair to its Stage-1 receipt and producer."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

import h5py


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.vetting.s63_preprocessing import (  # noqa: E402
    S63_COMPACT_ATTESTATION_CONTRACT,
    file_sha256,
    require_producer_git_sha,
    validate_producer_git_sha,
    validate_s63_stage1_compact,
    validate_sha256,
)


def attest_compact(
    *,
    accepted_validation: Path,
    expected_accepted_validation_sha256: str,
    compact_h5: Path,
    expected_compact_h5_sha256: str,
    compact_manifest: Path,
    expected_compact_manifest_sha256: str,
    producer_git_sha: str,
    published_h5_path: Path | None = None,
) -> dict[str, object]:
    accepted_validation = Path(accepted_validation).resolve()
    compact_h5 = Path(compact_h5).resolve()
    compact_manifest = Path(compact_manifest).resolve()
    expected_inputs = {
        accepted_validation: validate_sha256(
            expected_accepted_validation_sha256,
            label="expected accepted validation SHA-256",
        ),
        compact_h5: validate_sha256(
            expected_compact_h5_sha256, label="expected compact HDF5 SHA-256"
        ),
        compact_manifest: validate_sha256(
            expected_compact_manifest_sha256,
            label="expected compact manifest SHA-256",
        ),
    }
    for path, expected in expected_inputs.items():
        if file_sha256(path) != expected:
            raise ValueError(f"pre-attestation SHA-256 mismatch: {path}")
    producer = validate_producer_git_sha(producer_git_sha)
    accepted = json.loads(accepted_validation.read_text(encoding="utf-8"))
    require_producer_git_sha(accepted, producer, label="Stage-1 receipt")
    validate_s63_stage1_compact(
        accepted_validation_path=accepted_validation,
        compact_h5_path=compact_h5,
        compact_manifest_path=compact_manifest,
        expected_producer_git_sha=None,
        require_attestation=False,
    )
    manifest = json.loads(compact_manifest.read_text(encoding="utf-8"))
    if not isinstance(manifest, dict):
        raise ValueError("compact manifest must contain one object")
    source_host_out_h5 = str(manifest["out_h5"])
    if published_h5_path is not None:
        published_h5 = Path(published_h5_path).expanduser().resolve()
        if not published_h5.is_absolute() or published_h5.exists():
            raise ValueError("published compact HDF5 path must be fresh and absolute")
        source_host_out_h5 = str(published_h5)
        manifest["out_h5"] = source_host_out_h5
    accepted_sha256 = file_sha256(accepted_validation)
    with h5py.File(compact_h5, "r+") as h5:
        existing = h5.attrs.get("producer_git_sha")
        if existing is not None and validate_producer_git_sha(existing) != producer:
            raise ValueError("compact HDF5 already declares another producer Git SHA")
        h5.attrs["producer_git_sha"] = producer
        h5.attrs["attestation_contract_version"] = (
            S63_COMPACT_ATTESTATION_CONTRACT
        )
        h5.attrs["accepted_stage1_validation_sha256"] = accepted_sha256
        h5.attrs["source_host_out_h5"] = source_host_out_h5
        h5.flush()
    manifest.update(
        {
            "producer_git_sha": producer,
            "attestation_contract_version": S63_COMPACT_ATTESTATION_CONTRACT,
            "accepted_stage1_validation_sha256": accepted_sha256,
            "out_h5_sha256": file_sha256(compact_h5),
        }
    )
    pending = compact_manifest.with_suffix(
        compact_manifest.suffix + ".attest_pending"
    )
    if pending.exists():
        raise FileExistsError(f"stale compact attestation pending file: {pending}")
    pending.write_text(
        json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    os.replace(pending, compact_manifest)
    audit = validate_s63_stage1_compact(
        accepted_validation_path=accepted_validation,
        compact_h5_path=compact_h5,
        compact_manifest_path=compact_manifest,
        expected_producer_git_sha=producer,
    )
    return {key: value for key, value in audit.items() if key != "target_tics"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accepted-validation", type=Path, required=True)
    parser.add_argument("--expected-accepted-validation-sha256", required=True)
    parser.add_argument("--compact-h5", type=Path, required=True)
    parser.add_argument("--expected-compact-h5-sha256", required=True)
    parser.add_argument("--compact-manifest", type=Path, required=True)
    parser.add_argument("--expected-compact-manifest-sha256", required=True)
    parser.add_argument("--producer-git-sha", required=True)
    parser.add_argument("--published-h5-path", type=Path)
    args = parser.parse_args()
    result = attest_compact(
        accepted_validation=args.accepted_validation,
        expected_accepted_validation_sha256=(
            args.expected_accepted_validation_sha256
        ),
        compact_h5=args.compact_h5,
        expected_compact_h5_sha256=args.expected_compact_h5_sha256,
        compact_manifest=args.compact_manifest,
        expected_compact_manifest_sha256=args.expected_compact_manifest_sha256,
        producer_git_sha=args.producer_git_sha,
        published_h5_path=args.published_h5_path,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
