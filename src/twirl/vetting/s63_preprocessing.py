"""Fail-closed provenance checks for prospective S63 preprocessing products."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
from typing import Any, Mapping

import h5py

from twirl.vetting.adp_only import ADP_ONLY_APERTURES


S63 = 63
S63_ORBITS = (133, 134)
S63_STAGE1_RECEIPT_ATTESTATION_CONTRACT = (
    "twirl_teacher_v3_s63_stage1_receipt_attestation_v1"
)
S63_COMPACT_ATTESTATION_CONTRACT = "twirl_teacher_v3_s63_compact_attestation_v1"
_GIT_SHA = re.compile(r"[0-9a-f]{40}")
_SHA256 = re.compile(r"[0-9a-f]{64}")
_COMPACT_SKIP_KEYS = {
    "read_failed",
    "tic_filter",
    "duplicate_tic",
    "no_flux_columns",
}


def file_sha256(path: Path, *, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(chunk_size):
            digest.update(block)
    return digest.hexdigest()


def validate_producer_git_sha(value: Any, *, label: str = "producer Git SHA") -> str:
    normalized = str(value).strip().lower()
    if _GIT_SHA.fullmatch(normalized) is None:
        raise ValueError(f"{label} must be a full lowercase 40-character Git SHA")
    return normalized


def validate_sha256(value: Any, *, label: str) -> str:
    normalized = str(value).strip().lower()
    if _SHA256.fullmatch(normalized) is None:
        raise ValueError(f"{label} must be a resolved SHA-256 digest")
    return normalized


def require_producer_git_sha(
    payload: Mapping[str, Any],
    expected_git_sha: str,
    *,
    label: str,
) -> str:
    expected = validate_producer_git_sha(expected_git_sha)
    observed = validate_producer_git_sha(
        payload.get("producer_git_sha"), label=f"{label} producer_git_sha"
    )
    if observed != expected:
        raise ValueError(
            f"{label} producer_git_sha={observed} differs from launch Git SHA {expected}"
        )
    return observed


def _load_json(path: Path, *, label: str) -> dict[str, Any]:
    path = Path(path)
    if not path.is_file() or path.stat().st_size <= 0:
        raise ValueError(f"{label} must be a nonempty regular file: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{label} must contain one JSON object")
    return payload


def _absolute_normalized(value: Any, *, label: str) -> Path:
    text = str(value).strip()
    path = Path(text).expanduser()
    if not text or not path.is_absolute():
        raise ValueError(f"{label} must be an absolute path")
    return Path(os.path.normpath(str(path)))


def _under_root(path: Path, root: Path) -> bool:
    return path != root and root in path.parents


def _int(value: Any, *, label: str, minimum: int | None = None) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{label} must be an integer")
    try:
        integer = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be an integer") from exc
    if isinstance(value, float) and not value.is_integer():
        raise ValueError(f"{label} must be an integer")
    if isinstance(value, str) and str(integer) != value.strip():
        raise ValueError(f"{label} must be an integer")
    if minimum is not None and integer < minimum:
        raise ValueError(f"{label} must be >= {minimum}")
    return integer


def validate_s63_stage1_receipt_attestation(
    path: Path,
    *,
    expected_producer_git_sha: str,
) -> dict[str, Any]:
    """Validate the immutable-source/copy attestation for the S63 receipt.

    The accepted A2v1 validation in the production log tree is immutable.  A
    prospective run consumes an exact copy carrying this attestation instead.
    Cross-host consumers may not have the source path mounted; when it is
    mounted, its current bytes must still match the recorded original digest.
    """

    attested_path = Path(path).expanduser().resolve(strict=True)
    payload = _load_json(attested_path, label="attested S63 Stage-1 receipt")
    require_producer_git_sha(
        payload, expected_producer_git_sha, label="attested Stage-1 receipt"
    )
    if payload.get("attestation_contract_version") != (
        S63_STAGE1_RECEIPT_ATTESTATION_CONTRACT
    ):
        raise ValueError("Stage-1 receipt attestation contract is missing")
    source_path = _absolute_normalized(
        payload.get("source_validation_path"),
        label="Stage-1 receipt source_validation_path",
    )
    source_sha256 = validate_sha256(
        payload.get("source_validation_sha256"),
        label="Stage-1 receipt source_validation_sha256",
    )
    pre_attestation_sha256 = validate_sha256(
        payload.get("pre_attestation_sha256"),
        label="Stage-1 receipt pre_attestation_sha256",
    )
    if source_sha256 != pre_attestation_sha256:
        raise ValueError("attested Stage-1 copy was not byte-identical to its source")
    if source_path.resolve(strict=False) == attested_path:
        raise ValueError("Stage-1 receipt source and attested copy must be distinct")
    if source_path.exists():
        if not source_path.is_file() or file_sha256(source_path) != source_sha256:
            raise ValueError("original Stage-1 validation no longer matches its receipt")
        source_verified_here = True
    else:
        source_verified_here = False
    return {
        "source_validation_path": str(source_path),
        "source_validation_sha256": source_sha256,
        "source_verified_on_this_host": source_verified_here,
        "attested_validation_sha256": file_sha256(attested_path),
        "producer_git_sha": validate_producer_git_sha(expected_producer_git_sha),
        "passed": True,
    }


def validate_s63_stage1_compact(
    *,
    accepted_validation_path: Path,
    compact_h5_path: Path,
    compact_manifest_path: Path,
    expected_producer_git_sha: str | None,
    require_attestation: bool = True,
) -> dict[str, Any]:
    """Prove that a compact ADP pair is the exact accepted S63 FITS receipt.

    Source-host paths remain authoritative strings.  A hash-identical compact
    HDF5 may be relocated to another host only when the manifest's original
    ``out_h5`` path is absent there; if that path exists, it must resolve to the
    supplied file.
    """

    accepted_validation_path = Path(accepted_validation_path).resolve()
    compact_h5_path = Path(compact_h5_path).resolve()
    compact_manifest_path = Path(compact_manifest_path).resolve()
    accepted = _load_json(
        accepted_validation_path, label="accepted S63 Stage-1 validation"
    )
    manifest = _load_json(compact_manifest_path, label="S63 compact manifest")
    if accepted.get("sector") != S63 or sorted(accepted.get("orbits", [])) != list(
        S63_ORBITS
    ):
        raise ValueError("accepted Stage-1 receipt is not S63 orbits 133/134")
    if any(accepted.get(field) is not True for field in ("ok", "ok_h5", "ok_fits")):
        raise ValueError("accepted Stage-1 receipt did not pass")
    fits = accepted.get("fits")
    if not isinstance(fits, Mapping) or fits.get("skipped") is True:
        raise ValueError("accepted Stage-1 receipt lacks a FITS audit")
    for field in (
        "n_bad_checked_fits",
        "n_extra_fits_tics",
        "n_missing_fits_non_edge_tics",
    ):
        if _int(fits.get(field, -1), label=f"accepted fits.{field}") != 0:
            raise ValueError(f"accepted Stage-1 fits.{field} must be zero")
    receipt_root = _absolute_normalized(
        fits.get("hlsp_root"), label="accepted Stage-1 fits.hlsp_root"
    )
    receipt_n_fits = _int(
        fits.get("n_fits"), label="accepted Stage-1 fits.n_fits", minimum=1
    )
    receipt_n_tics = _int(
        fits.get("n_found_unique_tics"),
        label="accepted Stage-1 fits.n_found_unique_tics",
        minimum=1,
    )

    if manifest.get("sector") != S63:
        raise ValueError("compact manifest is not S63")
    if manifest.get("requested_columns") != list(ADP_ONLY_APERTURES):
        raise ValueError("compact manifest is not exactly the locked ADP pair")
    manifest_root = _absolute_normalized(
        manifest.get("hlsp_root"), label="compact manifest hlsp_root"
    )
    if manifest_root != receipt_root:
        raise ValueError(
            "compact manifest hlsp_root differs from accepted Stage-1 fits.hlsp_root"
        )
    discovered = _int(
        manifest.get("n_discovered_files"),
        label="compact manifest n_discovered_files",
        minimum=1,
    )
    exported = _int(
        manifest.get("n_exported_targets"),
        label="compact manifest n_exported_targets",
        minimum=1,
    )
    if discovered != receipt_n_fits or exported != receipt_n_tics:
        raise ValueError(
            "compact discovered/exported counts differ from the accepted FITS receipt"
        )
    skipped = manifest.get("skipped")
    if not isinstance(skipped, Mapping) or not _COMPACT_SKIP_KEYS.issubset(skipped):
        raise ValueError("compact manifest lacks the complete skip-count mapping")
    invalid_skip = {
        str(key): value
        for key, value in skipped.items()
        if _int(value, label=f"compact skipped.{key}", minimum=0) != 0
    }
    if invalid_skip:
        raise ValueError(f"accepted S63 compact skip counts must all be zero: {invalid_skip}")
    if discovered != exported:
        raise ValueError("zero-skip compact receipt must discover and export equal counts")

    records = manifest.get("records")
    if not isinstance(records, list) or len(records) != exported:
        raise ValueError("compact manifest records do not match its exported count")
    record_by_tic: dict[int, Mapping[str, Any]] = {}
    for index, record in enumerate(records):
        if not isinstance(record, Mapping):
            raise ValueError(f"compact manifest record {index} is not an object")
        tic = _int(record.get("tic"), label=f"compact record {index} TIC", minimum=1)
        if tic in record_by_tic:
            raise ValueError(f"compact manifest contains duplicate TIC {tic}")
        if _int(record.get("sector"), label=f"compact record {tic} sector") != S63:
            raise ValueError(f"compact manifest record {tic} is not S63")
        for field in ("camera", "ccd"):
            detector = _int(
                record.get(field), label=f"compact record {tic} {field}"
            )
            if detector not in range(1, 5):
                raise ValueError(f"compact manifest record {tic} has invalid {field}")
        if record.get("flux_columns") != list(ADP_ONLY_APERTURES):
            raise ValueError(f"compact manifest record {tic} lacks the exact ADP pair")
        source_fits = _absolute_normalized(
            record.get("source_fits"), label=f"compact record {tic} source_fits"
        )
        if not _under_root(source_fits, receipt_root):
            raise ValueError(
                f"compact manifest record {tic} source_fits is outside the accepted root"
            )
        record_by_tic[tic] = record

    if not compact_h5_path.is_file() or compact_h5_path.stat().st_size <= 0:
        raise ValueError("compact HDF5 must be a nonempty regular file")
    compact_sha256 = file_sha256(compact_h5_path)
    declared_out = _absolute_normalized(
        manifest.get("out_h5"), label="compact manifest out_h5"
    )
    canonical_match = declared_out.resolve(strict=False) == compact_h5_path
    if declared_out.exists() and not canonical_match:
        raise ValueError(
            "compact manifest source-host out_h5 exists but is not the supplied product"
        )
    relocated_by_hash = not canonical_match
    accepted_sha256 = file_sha256(accepted_validation_path)

    with h5py.File(compact_h5_path, "r") as h5:
        if int(h5.attrs.get("sector", -1)) != S63:
            raise ValueError("compact HDF5 is not S63")
        try:
            columns = json.loads(str(h5.attrs["flux_columns"]))
        except (KeyError, TypeError, json.JSONDecodeError) as exc:
            raise ValueError("compact HDF5 lacks valid flux_columns") from exc
        if columns != list(ADP_ONLY_APERTURES):
            raise ValueError("compact HDF5 is not exactly the locked ADP pair")
        h5_root = _absolute_normalized(
            h5.attrs.get("hlsp_root"), label="compact HDF5 hlsp_root"
        )
        if h5_root != receipt_root:
            raise ValueError("compact HDF5 hlsp_root differs from the accepted root")
        if "targets" not in h5:
            raise ValueError("compact HDF5 lacks /targets")
        target_keys = list(h5["targets"].keys())
        try:
            target_tics = {int(key) for key in target_keys}
        except ValueError as exc:
            raise ValueError("compact HDF5 has a non-TIC target key") from exc
        if len(target_tics) != len(target_keys) or target_tics != set(record_by_tic):
            raise ValueError("compact manifest records do not exactly equal /targets")
        if int(h5.attrs.get("n_targets", -1)) != exported:
            raise ValueError("compact HDF5 n_targets differs from the FITS receipt")
        required_datasets = {
            "time",
            "cadenceno",
            "quality",
            "orbitid",
            *ADP_ONLY_APERTURES,
        }
        for key, group in h5["targets"].items():
            tic = int(key)
            record = record_by_tic[tic]
            missing = sorted(required_datasets - set(group.keys()))
            if missing:
                raise ValueError(f"compact target {tic} lacks datasets: {missing}")
            lengths = {
                int(group[name].shape[0])
                for name in required_datasets
                if len(group[name].shape) == 1
            }
            if len(lengths) != 1 or next(iter(lengths), 0) <= 0 or any(
                len(group[name].shape) != 1 for name in required_datasets
            ):
                raise ValueError(f"compact target {tic} has misaligned datasets")
            n_cadences = next(iter(lengths))
            expected_identity = {
                "tic": tic,
                "sector": S63,
                "camera": _int(record["camera"], label=f"record {tic} camera"),
                "ccd": _int(record["ccd"], label=f"record {tic} ccd"),
            }
            for field, expected in expected_identity.items():
                if _int(group.attrs.get(field, -1), label=f"target {tic} {field}") != expected:
                    raise ValueError(f"compact target {tic} {field} differs from manifest")
            if _int(record.get("n_cadences"), label=f"record {tic} n_cadences") != n_cadences:
                raise ValueError(f"compact target {tic} cadence count differs from manifest")
            record_source = _absolute_normalized(
                record.get("source_fits"), label=f"record {tic} source_fits"
            )
            group_source = _absolute_normalized(
                group.attrs.get("source_fits"), label=f"target {tic} source_fits"
            )
            if group_source != record_source or not _under_root(group_source, receipt_root):
                raise ValueError(f"compact target {tic} source_fits provenance mismatch")

        if require_attestation:
            expected_git_sha = validate_producer_git_sha(expected_producer_git_sha)
            validate_s63_stage1_receipt_attestation(
                accepted_validation_path,
                expected_producer_git_sha=expected_git_sha,
            )
            require_producer_git_sha(manifest, expected_git_sha, label="compact manifest")
            h5_git_sha = validate_producer_git_sha(
                h5.attrs.get("producer_git_sha"),
                label="compact HDF5 producer_git_sha",
            )
            if h5_git_sha != expected_git_sha:
                raise ValueError("compact HDF5 producer Git SHA differs from launch")
            if str(manifest.get("attestation_contract_version", "")) != (
                S63_COMPACT_ATTESTATION_CONTRACT
            ) or str(h5.attrs.get("attestation_contract_version", "")) != (
                S63_COMPACT_ATTESTATION_CONTRACT
            ):
                raise ValueError("compact attestation contract is missing or incompatible")
            declared_h5_hash = validate_sha256(
                manifest.get("out_h5_sha256"), label="compact manifest out_h5_sha256"
            )
            if declared_h5_hash != compact_sha256:
                raise ValueError("compact manifest does not hash-bind the supplied HDF5")
            for label, observed in (
                (
                    "compact manifest accepted receipt",
                    manifest.get("accepted_stage1_validation_sha256"),
                ),
                (
                    "compact HDF5 accepted receipt",
                    h5.attrs.get("accepted_stage1_validation_sha256"),
                ),
            ):
                if validate_sha256(observed, label=label) != accepted_sha256:
                    raise ValueError(f"{label} hash binding mismatch")
            h5_source_out = _absolute_normalized(
                h5.attrs.get("source_host_out_h5"),
                label="compact HDF5 source_host_out_h5",
            )
            if h5_source_out != declared_out:
                raise ValueError("compact HDF5/manifest source-host out_h5 disagree")

    return {
        "accepted_stage1_validation_sha256": accepted_sha256,
        "hlsp_root": str(receipt_root),
        "n_discovered_files": discovered,
        "n_exported_targets": exported,
        "n_manifest_records": len(record_by_tic),
        "compact_h5_sha256": compact_sha256,
        "source_host_out_h5": str(declared_out),
        "canonical_source_path_match": canonical_match,
        "relocated_by_hash": relocated_by_hash,
        "target_tics": sorted(record_by_tic),
        "passed": True,
    }


__all__ = [
    "S63_COMPACT_ATTESTATION_CONTRACT",
    "S63_STAGE1_RECEIPT_ATTESTATION_CONTRACT",
    "file_sha256",
    "require_producer_git_sha",
    "validate_producer_git_sha",
    "validate_s63_stage1_compact",
    "validate_s63_stage1_receipt_attestation",
    "validate_sha256",
]
