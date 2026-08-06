#!/usr/bin/env python3
"""Derive the authoritative S63 model-ready TICs from compact ``/targets``."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import sys
from typing import Any

import h5py


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.vetting.adp_only import ADP_ONLY_APERTURES  # noqa: E402
from twirl.vetting.teacher_v3_prospective import (  # noqa: E402
    file_sha256,
    load_prospective_plan,
    read_tic_inventory,
    tic_inventory_sha256,
    validate_sha256,
)
from twirl.vetting.s63_preprocessing import (  # noqa: E402
    validate_producer_git_sha,
    validate_s63_stage1_compact,
)


CONTRACT_VERSION = "twirl_teacher_v3_s63_model_ready_allowlist_v1"


def _bound_json(path: Path, expected_sha256: str, *, label: str) -> dict[str, Any]:
    expected = validate_sha256(expected_sha256, context=f"expected {label}")
    before = file_sha256(path)
    if before != expected:
        raise ValueError(f"{label} SHA-256 mismatch")
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if file_sha256(path) != before:
        raise RuntimeError(f"{label} changed while it was read")
    if not isinstance(payload, dict):
        raise ValueError(f"{label} must be a JSON object")
    return payload


def _validate_accepted_stage1(payload: dict[str, Any]) -> None:
    if payload.get("sector") != 63 or sorted(payload.get("orbits", [])) != [133, 134]:
        raise ValueError("accepted Stage-1 validation is not S63 orbits 133/134")
    if any(payload.get(field) is not True for field in ("ok", "ok_h5", "ok_fits")):
        raise ValueError("accepted Stage-1 validation did not pass")
    h5 = payload.get("h5", {})
    fits = payload.get("fits", {})
    if any(
        h5.get(field) != 0
        for field in ("n_missing_h5_non_edge", "n_zero_byte_h5", "n_unreadable_h5")
    ):
        raise ValueError("accepted Stage-1 HDF5 validation is blocking")
    if any(
        fits.get(field) != 0
        for field in ("n_missing_fits_non_edge_tics", "n_bad_checked_fits")
    ):
        raise ValueError("accepted Stage-1 FITS validation is blocking")


def derive_model_ready_tics(compact_h5: Path, *, expected_sha256: str) -> tuple[int, ...]:
    """Read only compact group identities/attrs; no S63 light-curve tensors."""

    expected = validate_sha256(expected_sha256, context="expected compact HDF5")
    before = file_sha256(compact_h5)
    if before != expected:
        raise ValueError("compact HDF5 SHA-256 mismatch")
    with h5py.File(compact_h5, "r") as h5:
        if int(h5.attrs.get("sector", -1)) != 63:
            raise ValueError("compact HDF5 sector is not 63")
        try:
            columns = json.loads(str(h5.attrs["flux_columns"]))
        except (KeyError, TypeError, json.JSONDecodeError) as exc:
            raise ValueError("compact HDF5 lacks valid flux_columns") from exc
        if columns != list(ADP_ONLY_APERTURES):
            raise ValueError("compact HDF5 is not exactly the locked ADP pair")
        if "targets" not in h5:
            raise ValueError("compact HDF5 lacks /targets")
        try:
            tics = tuple(sorted(int(key) for key in h5["targets"].keys()))
        except ValueError as exc:
            raise ValueError("compact HDF5 contains a non-TIC target group") from exc
        if not tics or any(tic <= 0 for tic in tics) or len(tics) != len(set(tics)):
            raise ValueError("compact HDF5 target identities are invalid")
        if int(h5.attrs.get("n_targets", -1)) != len(tics):
            raise ValueError("compact HDF5 n_targets disagrees with /targets")
        required_datasets = {
            "time",
            "cadenceno",
            "quality",
            "orbitid",
            *ADP_ONLY_APERTURES,
        }
        for key, group in h5["targets"].items():
            missing = sorted(required_datasets - set(group.keys()))
            if missing:
                raise ValueError(f"compact target {key} lacks datasets: {missing}")
            shapes = {name: group[name].shape for name in required_datasets}
            if any(len(shape) != 1 for shape in shapes.values()):
                raise ValueError(f"compact target {key} has non-vector datasets")
            lengths = {shape[0] for shape in shapes.values()}
            if len(lengths) != 1 or next(iter(lengths), 0) <= 0:
                raise ValueError(
                    f"compact target {key} has empty/misaligned dataset lengths"
                )
            tic = int(key)
            if int(group.attrs.get("tic", -1)) != tic or int(
                group.attrs.get("sector", -1)
            ) != 63:
                raise ValueError(f"compact target {key} has invalid TIC/sector attrs")
            for field in ("camera", "ccd"):
                if int(group.attrs.get(field, -1)) not in range(1, 5):
                    raise ValueError(f"compact target {key} has invalid {field}")
    if file_sha256(compact_h5) != before:
        raise RuntimeError("compact HDF5 changed while target identities were read")
    return tics


def _publish_pair(
    *,
    tics: tuple[int, ...],
    summary: dict[str, Any],
    output: Path,
    summary_path: Path,
) -> None:
    output = Path(output).resolve()
    summary_path = Path(summary_path).resolve()
    if output.exists() or summary_path.exists():
        raise FileExistsError("refusing to overwrite model-ready outputs")
    output.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    claim = summary_path.with_suffix(summary_path.suffix + ".claim")
    output_tmp = output.with_suffix(output.suffix + ".tmp")
    summary_tmp = summary_path.with_suffix(summary_path.suffix + ".tmp")
    descriptor: int | None = None
    published_output = False
    try:
        descriptor = os.open(claim, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        os.write(descriptor, f"pid={os.getpid()}\n".encode("ascii"))
        os.close(descriptor)
        descriptor = None
        output_tmp.write_text("".join(f"{tic}\n" for tic in tics), encoding="ascii")
        summary["outputs"] = {
            "model_ready_allowlist": {
                "path": str(output),
                "sha256": file_sha256(output_tmp),
                "size_bytes": int(output_tmp.stat().st_size),
            },
            "summary": {"path": str(summary_path)},
        }
        summary_tmp.write_text(
            json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        os.replace(output_tmp, output)
        published_output = True
        os.replace(summary_tmp, summary_path)
    except Exception:
        if published_output:
            output.unlink(missing_ok=True)
        raise
    finally:
        if descriptor is not None:
            os.close(descriptor)
        output_tmp.unlink(missing_ok=True)
        summary_tmp.unlink(missing_ok=True)
        claim.unlink(missing_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--expected-plan-sha256", required=True)
    parser.add_argument("--accepted-validation", type=Path, required=True)
    parser.add_argument("--expected-accepted-validation-sha256", required=True)
    parser.add_argument("--compact-h5", type=Path, required=True)
    parser.add_argument("--expected-compact-h5-sha256", required=True)
    parser.add_argument("--compact-manifest", type=Path, required=True)
    parser.add_argument("--expected-compact-manifest-sha256", required=True)
    parser.add_argument("--reserved-tics", type=Path, required=True)
    parser.add_argument("--expected-reserved-tics-sha256", required=True)
    parser.add_argument("--producer-git-sha", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()
    producer_git_sha = validate_producer_git_sha(args.producer_git_sha)

    plan, plan_sha = load_prospective_plan(
        args.plan,
        expected_sha256=args.expected_plan_sha256,
    )
    validation = _bound_json(
        args.accepted_validation,
        args.expected_accepted_validation_sha256,
        label="accepted Stage-1 validation",
    )
    _validate_accepted_stage1(validation)
    compact_manifest = _bound_json(
        args.compact_manifest,
        args.expected_compact_manifest_sha256,
        label="compact manifest",
    )
    if compact_manifest.get("sector") != 63 or compact_manifest.get(
        "requested_columns"
    ) != list(ADP_ONLY_APERTURES):
        raise ValueError("compact manifest is not the S63 locked ADP pair")
    compact_receipt = validate_s63_stage1_compact(
        accepted_validation_path=args.accepted_validation,
        compact_h5_path=args.compact_h5,
        compact_manifest_path=args.compact_manifest,
        expected_producer_git_sha=producer_git_sha,
    )
    reserved = read_tic_inventory(
        args.reserved_tics,
        expected_sha256=args.expected_reserved_tics_sha256,
    )
    reservation = plan["s63_identity_reservation"]
    if args.expected_reserved_tics_sha256 != reservation["reserved_tics_sha256"]:
        raise ValueError("reserved-TIC hash differs from prospective plan")
    if len(reserved) != int(reservation["n_requested_tics"]):
        raise ValueError("reserved-TIC count differs from prospective plan")
    model_ready = derive_model_ready_tics(
        args.compact_h5,
        expected_sha256=args.expected_compact_h5_sha256,
    )
    if list(model_ready) != compact_receipt["target_tics"]:
        raise RuntimeError("compact receipt target identities changed during derivation")
    if int(compact_manifest.get("n_exported_targets", -1)) != len(model_ready):
        raise ValueError("compact manifest target count disagrees with HDF5")
    missing_from_seal = sorted(set(model_ready) - set(reserved))
    if missing_from_seal:
        raise ValueError(
            "compact /targets contains TICs outside the sealed S63 reservation: "
            f"{missing_from_seal[:10]}"
        )
    summary = {
        "schema_version": CONTRACT_VERSION,
        "producer_git_sha": producer_git_sha,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "sector": 63,
        "orbits": [133, 134],
        "derivation": "sorted unique integer keys under validated compact /targets",
        "light_curve_tensors_read": False,
        "source_hashes": {
            "prospective_plan_sha256": plan_sha,
            "accepted_stage1_validation_sha256": args.expected_accepted_validation_sha256,
            "compact_h5_sha256": args.expected_compact_h5_sha256,
            "compact_manifest_sha256": args.expected_compact_manifest_sha256,
            "reserved_tics_sha256": args.expected_reserved_tics_sha256,
        },
        "counts": {
            "n_reserved_tics": len(reserved),
            "n_model_ready_tics": len(model_ready),
            "n_reserved_not_model_ready": len(set(reserved) - set(model_ready)),
        },
        "model_ready_tics_sha256": tic_inventory_sha256(model_ready),
        "partition_checks": {
            "model_ready_subset_of_reservation": True,
            "compact_target_count_exact": True,
        },
        "passed": True,
    }
    _publish_pair(
        tics=model_ready,
        summary=summary,
        output=args.out,
        summary_path=args.summary,
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
