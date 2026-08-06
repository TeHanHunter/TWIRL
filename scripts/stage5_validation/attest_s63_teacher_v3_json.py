#!/usr/bin/env python3
"""Atomically add an exact producer-Git attestation to S63 JSON outputs."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.vetting.s63_preprocessing import (  # noqa: E402
    S63_STAGE1_RECEIPT_ATTESTATION_CONTRACT,
    file_sha256,
    validate_producer_git_sha,
    validate_sha256,
)


def attest_json(
    path: Path,
    *,
    expected_sha256: str,
    producer_git_sha: str,
    source_json: Path | None = None,
    expected_source_sha256: str | None = None,
    bls_published_peak_table: Path | None = None,
    bls_published_summary: Path | None = None,
) -> dict[str, object]:
    path = Path(path).resolve()
    expected = validate_sha256(expected_sha256, label="expected JSON SHA-256")
    producer = validate_producer_git_sha(producer_git_sha)
    if file_sha256(path) != expected:
        raise ValueError(f"pre-attestation JSON SHA-256 mismatch: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("attested JSON must contain one object")
    existing = payload.get("producer_git_sha")
    if existing is not None and validate_producer_git_sha(existing) != producer:
        raise ValueError("JSON already declares a different producer_git_sha")
    if (source_json is None) != (expected_source_sha256 is None):
        raise ValueError(
            "source_json and expected_source_sha256 must be supplied together"
        )
    if source_json is not None:
        source = Path(source_json).expanduser().resolve(strict=True)
        source_sha256 = validate_sha256(
            expected_source_sha256, label="expected source JSON SHA-256"
        )
        if source == path:
            raise ValueError("source JSON and attested copy must be distinct files")
        if file_sha256(source) != source_sha256:
            raise ValueError("source JSON SHA-256 mismatch")
        if expected != source_sha256:
            raise ValueError("pre-attestation copy is not byte-identical to source")
        payload.update(
            {
                "attestation_contract_version": (
                    S63_STAGE1_RECEIPT_ATTESTATION_CONTRACT
                ),
                "source_validation_path": str(source),
                "source_validation_sha256": source_sha256,
                "pre_attestation_sha256": expected,
            }
        )
    if (bls_published_peak_table is None) != (bls_published_summary is None):
        raise ValueError("both published BLS paths must be supplied together")
    if bls_published_peak_table is not None:
        peak_path = Path(bls_published_peak_table).expanduser().resolve()
        summary_path = Path(bls_published_summary).expanduser().resolve()
        if not peak_path.is_absolute() or not summary_path.is_absolute():
            raise ValueError("published BLS paths must be absolute")
        outputs = payload.get("outputs")
        if not isinstance(outputs, dict) or "peak_table" not in outputs:
            raise ValueError("BLS summary lacks outputs.peak_table")
        outputs["peak_table"] = str(peak_path)
        outputs["summary"] = str(summary_path)
        payload["out_dir"] = str(peak_path.parent)
    if file_sha256(path) != expected:
        raise RuntimeError("JSON changed while its attestation was prepared")
    payload["producer_git_sha"] = producer
    temporary = path.with_suffix(path.suffix + ".attest_pending")
    if temporary.exists():
        raise FileExistsError(f"stale attestation pending file exists: {temporary}")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)
    return {
        "path": str(path),
        "producer_git_sha": producer,
        "sha256": file_sha256(path),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--expected-sha256", required=True)
    parser.add_argument("--producer-git-sha", required=True)
    parser.add_argument("--source-json", type=Path)
    parser.add_argument("--expected-source-sha256")
    parser.add_argument("--bls-published-peak-table", type=Path)
    parser.add_argument("--bls-published-summary", type=Path)
    args = parser.parse_args()
    result = attest_json(
        args.json,
        expected_sha256=args.expected_sha256,
        producer_git_sha=args.producer_git_sha,
        source_json=args.source_json,
        expected_source_sha256=args.expected_source_sha256,
        bls_published_peak_table=args.bls_published_peak_table,
        bls_published_summary=args.bls_published_summary,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
