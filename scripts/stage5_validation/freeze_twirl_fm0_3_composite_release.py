#!/usr/bin/env python3
"""Freeze the identity-only FM0.3 S56--S64 + S66--S77 release."""

from __future__ import annotations

import argparse
import json
from collections.abc import Mapping
from pathlib import Path

from twirl.models.fm0.composite_release import (
    LATER_SECTORS,
    SourceManifestBinding,
    freeze_composite_release,
)
from twirl.models.fm0.registry import FM0ContractError, sha256_file

AUTHORITY_SCHEMA_VERSION = "twirl_fm0_3_composite_freeze_authority_v1"


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--authority", required=True, type=Path)
    parser.add_argument("--authority-sha256", required=True)
    parser.add_argument("--producer-git-sha", required=True)
    parser.add_argument("--output", required=True, type=Path)
    return parser


def _object(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise FM0ContractError(f"{label} must be an object")
    return value


def main() -> None:
    args = _parser().parse_args()
    authority = args.authority.expanduser().resolve(strict=True)
    if sha256_file(authority) != args.authority_sha256:
        raise FM0ContractError("composite authority hash drifted")
    payload = _object(json.loads(authority.read_text(encoding="utf-8")), "authority")
    if (
        set(payload)
        != {
            "schema_version",
            "legacy",
            "later",
            "admission_receipt",
        }
        or payload["schema_version"] != AUTHORITY_SCHEMA_VERSION
    ):
        raise FM0ContractError("composite authority schema/fields drifted")

    legacy = _object(payload["legacy"], "legacy authority")
    if set(legacy) != {
        "release_root",
        "manifest_sha256",
        "receipt_path",
        "receipt_sha256",
    }:
        raise FM0ContractError("legacy authority fields drifted")
    legacy_binding = SourceManifestBinding(
        source_kind="legacy",
        sector_min=56,
        sector_max=64,
        release_root=Path(str(legacy["release_root"])),
        manifest_sha256=str(legacy["manifest_sha256"]),
        receipt_path=Path(str(legacy["receipt_path"])),
        receipt_sha256=str(legacy["receipt_sha256"]),
    )

    later_raw = payload["later"]
    if not isinstance(later_raw, list) or len(later_raw) != len(LATER_SECTORS):
        raise FM0ContractError("later authority must contain chronological S66--S77")
    later: list[SourceManifestBinding] = []
    for sector, raw_value in zip(LATER_SECTORS, later_raw, strict=True):
        raw = _object(raw_value, f"S{sector} authority")
        if (
            set(raw)
            != {
                "sector",
                "release_root",
                "manifest_sha256",
                "receipt_sha256",
            }
            or raw["sector"] != sector
        ):
            raise FM0ContractError(f"S{sector} authority fields/order drifted")
        root = Path(str(raw["release_root"]))
        later.append(
            SourceManifestBinding(
                source_kind="later",
                sector_min=sector,
                sector_max=sector,
                release_root=root,
                manifest_sha256=str(raw["manifest_sha256"]),
                receipt_path=root / "receipt.json",
                receipt_sha256=str(raw["receipt_sha256"]),
            )
        )

    admission = _object(payload["admission_receipt"], "admission authority")
    if set(admission) != {"path", "sha256"}:
        raise FM0ContractError("admission authority fields drifted")
    result = freeze_composite_release(
        args.output,
        legacy_binding=legacy_binding,
        later_bindings=later,
        admission_receipt_path=Path(str(admission["path"])),
        admission_receipt_sha256=str(admission["sha256"]),
        producer_git_sha=args.producer_git_sha,
    )
    print(
        json.dumps(
            {
                "root": str(result.root),
                "receipt_sha256": result.receipt_sha256,
                "source_bindings_sha256": result.source_bindings_sha256,
                "role_index_sha256": result.role_index_sha256,
                "selection": result.receipt["selection"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
