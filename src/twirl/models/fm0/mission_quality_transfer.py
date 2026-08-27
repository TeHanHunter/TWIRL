"""Package receipt-bound FM mission-quality authorities for ORCD transfer."""

from __future__ import annotations

import json
import os
import re
import shutil
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .mission_quality_admission import (
    MISSION_QUALITY_READY_STATE,
    MISSION_QUALITY_RECEIPT_SCHEMA_VERSION,
)
from .registry import FM0ContractError, sha256_file

MISSION_QUALITY_TRANSFER_SCHEMA_VERSION = "twirl_fm0_mission_quality_transfer_v1"
_SHA40 = re.compile(r"^[0-9a-f]{40}$")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")


def _load_receipt(path: Path) -> Mapping[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise FM0ContractError(f"invalid mission-quality receipt: {path}") from exc
    if not isinstance(payload, Mapping):
        raise FM0ContractError(f"mission-quality receipt is not an object: {path}")
    if (
        payload.get("schema_version") != MISSION_QUALITY_RECEIPT_SCHEMA_VERSION
        or payload.get("quality_state") != MISSION_QUALITY_READY_STATE
        or payload.get("passed") is not True
        or payload.get("panel_admission_authorized") is not False
    ):
        raise FM0ContractError(f"mission-quality receipt is not transferable: {path}")
    return payload


def _binding_destination(*, sector: int, binding: Mapping[str, Any]) -> Path:
    role = str(binding.get("role", ""))
    source_name = Path(str(binding.get("path", ""))).name
    if role == "qlp_quaternion":
        orbit = int(binding.get("orbit", -1))
        return Path(f"sources/s{sector:04d}/qlp/orbit-{orbit}/ffi/run") / source_name
    if role == "qlp_detector_qflag":
        orbit = int(binding.get("orbit", -1))
        return Path(f"sources/s{sector:04d}/qlp/orbit-{orbit}/ffi/run") / source_name
    if role in {
        "spoc_mission_quality",
        "tica_mission_quality",
        "tica_materialization_summary",
        "tica_materialization_ready",
    }:
        return Path(f"sources/s{sector:04d}/mission") / source_name
    raise FM0ContractError(f"unsupported mission-quality source role: {role!r}")


def package_mission_quality_transfer(
    *,
    receipt_paths: Sequence[str | Path],
    output_dir: str | Path,
    producer_git_sha: str,
) -> dict[str, Any]:
    """Copy verified compact authorities into one immutable transfer release."""

    if _SHA40.fullmatch(producer_git_sha) is None:
        raise FM0ContractError("producer_git_sha must be a full lowercase Git SHA")
    receipts = tuple(Path(path).expanduser().resolve() for path in receipt_paths)
    if not receipts or len(receipts) != len(set(receipts)):
        raise FM0ContractError("receipt_paths must be nonempty and unique")
    final = Path(output_dir).expanduser().resolve()
    partial = final.with_name(final.name + ".partial")
    if final.exists() or partial.exists():
        raise FM0ContractError(f"refusing to overwrite transfer release: {final}")
    partial.mkdir(parents=True)

    try:
        receipt_rows: list[dict[str, Any]] = []
        source_rows: list[dict[str, Any]] = []
        sectors: list[int] = []
        destinations: set[Path] = set()
        for receipt_path in receipts:
            receipt = _load_receipt(receipt_path)
            sector = int(receipt.get("sector", -1))
            if sector < 56 or sector in sectors:
                raise FM0ContractError("transfer receipts have invalid/duplicate sectors")
            sectors.append(sector)
            receipt_hash = sha256_file(receipt_path)
            staged_receipt = Path(f"receipts/s{sector:04d}.mission_quality.json")
            target = partial / staged_receipt
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(receipt_path, target)
            if sha256_file(target) != receipt_hash:
                raise FM0ContractError("staged mission-quality receipt hash mismatch")
            receipt_rows.append(
                {
                    "sector": sector,
                    "original_path": str(receipt_path),
                    "staged_path": str(staged_receipt),
                    "sha256": receipt_hash,
                    "bytes": target.stat().st_size,
                }
            )

            bindings = receipt.get("source_bindings")
            if not isinstance(bindings, list) or not bindings:
                raise FM0ContractError(
                    f"S{sector} mission-quality receipt lacks source bindings"
                )
            for binding in bindings:
                if not isinstance(binding, Mapping):
                    raise FM0ContractError(
                        f"S{sector} mission-quality source binding is invalid"
                    )
                source = Path(str(binding.get("path", ""))).expanduser().resolve()
                declared_hash = str(binding.get("sha256", ""))
                if (
                    not source.is_file()
                    or _SHA256.fullmatch(declared_hash) is None
                    or sha256_file(source) != declared_hash
                ):
                    raise FM0ContractError(
                        f"S{sector} source binding is missing or hash-mismatched: {source}"
                    )
                staged = _binding_destination(sector=sector, binding=binding)
                if staged in destinations:
                    raise FM0ContractError(f"duplicate staged source path: {staged}")
                destinations.add(staged)
                target = partial / staged
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source, target)
                if sha256_file(target) != declared_hash:
                    raise FM0ContractError(f"staged source hash mismatch: {staged}")
                source_rows.append(
                    {
                        "sector": sector,
                        "role": str(binding["role"]),
                        "original_path": str(source),
                        "staged_path": str(staged),
                        "sha256": declared_hash,
                        "bytes": target.stat().st_size,
                    }
                )

        if sectors != sorted(sectors):
            raise FM0ContractError("mission-quality receipts must be sector-sorted")
        manifest = {
            "schema_version": MISSION_QUALITY_TRANSFER_SCHEMA_VERSION,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "producer_git_sha": producer_git_sha,
            "sectors": sectors,
            "n_sectors": len(sectors),
            "n_receipts": len(receipt_rows),
            "n_source_files": len(source_rows),
            "n_bytes": sum(row["bytes"] for row in receipt_rows + source_rows),
            "receipts": receipt_rows,
            "sources": source_rows,
            "claim_limit": (
                "compact mission-quality transfer only; HDF5 joins, six-view "
                "shards, panel admission, and A2v1 acceptance remain unverified"
            ),
        }
        (partial / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        (partial / "READY").write_text(producer_git_sha + "\n", encoding="utf-8")
        os.replace(partial, final)
        return manifest
    except Exception:
        shutil.rmtree(partial, ignore_errors=True)
        raise


__all__ = [
    "MISSION_QUALITY_TRANSFER_SCHEMA_VERSION",
    "package_mission_quality_transfer",
]
