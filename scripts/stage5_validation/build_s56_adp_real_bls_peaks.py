#!/usr/bin/env python3
"""Build a full-S56 real-candidate BLS peak table from the ADP pair only."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Mapping

import h5py
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.io.compact_export import read_compact_lc_export  # noqa: E402
from twirl.io.hlsp import BJDREFI, HLSPLightCurve  # noqa: E402
from twirl.lightcurves.a2v1_cadence_reference import (  # noqa: E402
    AUTHORITY_EXCLUSION_POLICY,
    CADENCE_REFERENCE_BUILDER_VERSION,
)
from twirl.lightcurves.external_quality import (  # noqa: E402
    AUTHORITY_EXCLUSION_EXTERNAL_BIT,
    AUTHORITY_EXCLUSION_POLICY_CONTRACT,
    authority_exclusions_sha256,
)
from twirl.search.a2v1_bls_contract import (  # noqa: E402
    A2V1_TEACHER_BLS_SEARCH_CONTRACT,
    approved_a2v1_teacher_bls_config,
    bls_config_sha256,
)
from twirl.search.bls import BLSConfig, run_bls_on_lc  # noqa: E402
from twirl.vetting.adp_only import (  # noqa: E402
    ADP_ONLY_APERTURES,
    ADP_ONLY_CONTRACT_VERSION,
    validate_adp_only_apertures,
)
from twirl.vetting.recovery50_teacher import json_default, write_table  # noqa: E402
from twirl.vetting.two_aperture import (  # noqa: E402
    measure_two_aperture_candidate_metadata,
)
from twirl.vetting.ssl_full_pool_native import (  # noqa: E402
    FULL_POOL_NATIVE_DETREND_CONFIG_SHA256,
    FULL_POOL_NATIVE_DETREND_CONTRACT_VERSION,
    FULL_POOL_NATIVE_DETREND_TIME_CONTRACT_VERSION,
    _effective_quality_adp03q,
    load_production_raw_source_release,
    load_sector_pool_authority,
    select_sector_shard,
)


DEFAULT_COMPACT_LC = (
    REPO_ROOT / "data_local/stage3_injections/s56_twirlfs_v2_lc_export/"
    "s56_twirlfs_v2_adp_lc_export_pdo.h5"
)
DEFAULT_OUT_DIR = REPO_ROOT / "reports/stage5_validation/s56_adp_real_bls_peaks"

EXTERNAL_QUALITY_POLICY_CONTRACT = (
    "a2v1_bls_internal_or_authoritative_external_quality_v2"
)
CADENCE_REFERENCE_CONTRACT_TEMPLATE = "s{sector}_a2v1_cadence_reference_v1"
CADENCE_REFERENCE_COLUMNS = (
    "sector",
    "orbitid",
    "camera",
    "ccd",
    "cadenceno",
    "spoc_quality",
    "qlp_quality",
    "external_quality",
)
EXPECTED_CADENCE_AUTHORITY = "qlp_cam_quat"
EXPECTED_QUALITY_AUTHORITY = "spoc_and_qlp_quality_flags"
EXPECTED_QUALITY_COMPOSITION = {
    "external_quality": "spoc_quality | (qlp_quality << 30)",
    "qlp_quality_raw_values": [0, 1],
    "qlp_quality_external_bit": 30,
}
ORBITID_POLICIES = ("strict", "reference_by_cadence")
ORBITID_RECONCILIATION_CONTRACT_VERSION = "a2v1_compact_orbitid_reconciliation_v1"
TARGET_SELECTION_CONTRACT_VERSION = "a2v1_bls_target_allowlist_v1"
RAW_V4_INPUT_CONTRACT_VERSION = (
    "twirl_teacher_ssl_fullpool_raw_v1_detector_consistent_bls_v4"
)
S63_TARGET_METADATA_CONTRACT_VERSION = (
    "twirl_teacher_v3_s63_rank1_candidate_metadata_v1"
)
S63_TARGET_METADATA_COLUMNS = (
    "adp_sml_own_even_depth",
    "adp_sml_own_odd_depth",
    "adp_sml_own_even_odd_depth_delta",
    "adp_sml_own_even_odd_sigma_delta",
    "adp_sml_trend_ptp",
    "adp_own_even_depth",
    "adp_own_odd_depth",
    "adp_own_even_odd_depth_delta",
    "adp_own_even_odd_sigma_delta",
    "adp_trend_ptp",
)

# Populated once per worker by ``_initialize_external_quality_worker``.  The
# authoritative cadence map is deliberately not repeated in every task
# payload: a full S56 run has tens of thousands of targets but only sixteen
# detector maps.
_EXTERNAL_QUALITY_BY_DETECTOR: (
    dict[tuple[int, int, int], tuple[np.ndarray, np.ndarray, np.ndarray]] | None
) = None
_EXTERNAL_QUALITY_AUTHORITY_EXCLUSIONS: (
    dict[tuple[int, int, int], np.ndarray] | None
) = None
_EXTERNAL_QUALITY_PROVENANCE: dict[str, Any] | None = None
_ORBITID_POLICY = "strict"


def _sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def _read_table(path: Path) -> pd.DataFrame:
    suffix = Path(path).suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix == ".parquet":
        return pd.read_parquet(path)
    raise ValueError("cadence reference must end in .csv or .parquet")


def _integer_column(frame: pd.DataFrame, column: str) -> pd.Series:
    values = pd.to_numeric(frame[column], errors="raise")
    array = values.to_numpy(dtype=np.float64)
    if not np.all(np.isfinite(array)) or not np.all(array == np.floor(array)):
        raise ValueError(f"cadence-reference {column} must contain finite integers")
    return values.astype(np.int64)


def _canonical_json_sha256(payload: Mapping[str, Any]) -> str:
    text = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _validate_orbitid_policy(value: str) -> str:
    policy = str(value)
    if policy not in ORBITID_POLICIES:
        raise ValueError(
            f"orbitid_policy must be one of {ORBITID_POLICIES}; got {policy!r}"
        )
    return policy


def _validate_expected_sha256(
    value: str | None,
    *,
    context: str,
) -> str | None:
    if value is None:
        return None
    normalized = str(value).strip().lower()
    if len(normalized) != 64 or any(
        character not in "0123456789abcdef" for character in normalized
    ):
        raise ValueError(f"{context} must be a lowercase SHA-256 digest")
    return normalized


def _tic_inventory_sha256(tics: list[int] | tuple[int, ...]) -> str:
    """Hash a canonical, sorted, unique positive-TIC inventory."""

    normalized = sorted(int(value) for value in tics)
    if any(value <= 0 for value in normalized) or len(normalized) != len(
        set(normalized)
    ):
        raise ValueError("TIC inventory must contain unique positive integers")
    payload = "".join(f"{value}\n" for value in normalized)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def load_target_allowlist(path: Path) -> tuple[list[int], dict[str, Any]]:
    """Load a checksum-bound CSV or one-TIC-per-line text allowlist."""

    path = Path(path)
    digest = _sha256(path)
    suffix = path.suffix.lower()
    if suffix == ".csv":
        frame = pd.read_csv(path, dtype=str, keep_default_na=False)
        if "tic" not in frame.columns:
            raise ValueError("target allowlist CSV must contain a 'tic' column")
        raw_values = frame["tic"].astype(str).str.strip().tolist()
    elif suffix in {".txt", ".list"}:
        raw_values = []
        for raw_line in path.read_text(encoding="utf-8").splitlines():
            value = raw_line.strip()
            if not value or value.startswith("#"):
                continue
            if not raw_values and value.lower() == "tic":
                continue
            raw_values.append(value)
    else:
        raise ValueError("target allowlist must end in .csv, .txt, or .list")
    if not raw_values:
        raise ValueError("target allowlist is empty")
    invalid = [
        value
        for value in raw_values
        if not value or not value.isascii() or not value.isdigit()
    ]
    if invalid:
        raise ValueError(
            "target allowlist TIC values must be positive base-10 integers; "
            f"examples={invalid[:10]}"
        )
    tics = [int(value) for value in raw_values]
    if any(value <= 0 or value > np.iinfo(np.int64).max for value in tics):
        raise ValueError("target allowlist TIC values must be positive int64 integers")
    duplicate = pd.Series(tics, dtype=np.int64).duplicated(keep=False)
    if duplicate.any():
        examples = sorted(
            set(pd.Series(tics, dtype=np.int64).loc[duplicate].astype(int))
        )[:10]
        raise ValueError(f"target allowlist contains duplicate TICs: {examples}")
    tics = sorted(tics)
    return tics, {
        "target_selection_contract_version": TARGET_SELECTION_CONTRACT_VERSION,
        "target_allowlist": str(path),
        "target_allowlist_sha256": digest,
        "target_allowlist_count": len(tics),
        "target_allowlist_tics_sha256": _tic_inventory_sha256(tics),
    }


def _manifest_integer(value: Any, *, context: str) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{context} must be an integer")
    numeric = int(value)
    if numeric < 0 or numeric > np.iinfo(np.int64).max:
        raise ValueError(f"{context} must be a nonnegative int64 integer")
    return numeric


def _validate_authority_exclusions(
    *,
    manifest: Mapping[str, Any],
    frame: pd.DataFrame,
    sector: int,
    detector_names: list[str],
) -> tuple[dict[tuple[int, int, int], np.ndarray], str, int]:
    """Validate the exact cadence exclusions copied from SPOC provenance."""

    declaration = manifest.get("authority_exclusions")
    if not isinstance(declaration, Mapping):
        raise ValueError("cadence-reference authority_exclusions must be an object")
    required_fields = {
        "contract_version",
        "policy",
        "external_bit",
        "n_rows",
        "by_detector",
    }
    if set(declaration) != required_fields:
        raise ValueError("cadence-reference authority_exclusions has the wrong fields")
    if str(declaration["contract_version"]) != AUTHORITY_EXCLUSION_POLICY_CONTRACT:
        raise ValueError("cadence-reference authority-exclusion contract mismatch")
    if str(declaration["policy"]) != AUTHORITY_EXCLUSION_POLICY:
        raise ValueError("cadence-reference authority-exclusion policy mismatch")
    external_bit = _manifest_integer(
        declaration["external_bit"],
        context="cadence-reference authority-exclusion external_bit",
    )
    if external_bit != AUTHORITY_EXCLUSION_EXTERNAL_BIT:
        raise ValueError("cadence-reference authority-exclusion bit mismatch")

    declared_sha256 = str(manifest.get("authority_exclusions_sha256", "")).lower()
    observed_sha256 = authority_exclusions_sha256(declaration)
    if declared_sha256 != observed_sha256:
        raise ValueError("cadence-reference authority-exclusions hash mismatch")

    declared_total = _manifest_integer(
        declaration["n_rows"],
        context="cadence-reference authority-exclusion n_rows",
    )
    aggregate_total = _manifest_integer(
        manifest.get("n_spoc_rows_excluded_by_quat"),
        context="cadence-reference n_spoc_rows_excluded_by_quat",
    )
    if aggregate_total != declared_total:
        raise ValueError(
            "cadence-reference authority-exclusion count disagrees with "
            "n_spoc_rows_excluded_by_quat"
        )
    by_detector = declaration["by_detector"]
    if not isinstance(by_detector, Mapping):
        raise ValueError(
            "cadence-reference authority-exclusion by_detector must be an object"
        )
    if set(str(value) for value in by_detector) != set(detector_names):
        raise ValueError(
            "cadence-reference authority-exclusion detector inventory mismatch"
        )

    output: dict[tuple[int, int, int], np.ndarray] = {}
    observed_total = 0
    for camera, ccd in sorted(
        {
            (int(camera), int(ccd))
            for camera, ccd in zip(frame["camera"], frame["ccd"], strict=True)
        }
    ):
        name = f"cam{camera}_ccd{ccd}"
        detector = by_detector[name]
        if not isinstance(detector, Mapping):
            raise ValueError(
                f"cadence-reference authority exclusion {name} must be an object"
            )
        if set(detector) != {"n_rows", "rows"}:
            raise ValueError(
                f"cadence-reference authority exclusion {name} has the wrong fields"
            )
        rows = detector["rows"]
        if not isinstance(rows, list):
            raise ValueError(
                f"cadence-reference authority exclusion {name} rows must be a list"
            )
        declared_detector_rows = _manifest_integer(
            detector["n_rows"],
            context=f"cadence-reference authority exclusion {name} n_rows",
        )
        if declared_detector_rows != len(rows):
            raise ValueError(
                f"cadence-reference authority exclusion {name} row-count mismatch"
            )
        cadences: list[int] = []
        for index, row in enumerate(rows):
            if not isinstance(row, Mapping):
                raise ValueError(
                    f"cadence-reference authority exclusion {name} row {index} "
                    "must be an object"
                )
            if set(row) != {"cadenceno", "spoc_quality"}:
                raise ValueError(
                    f"cadence-reference authority exclusion {name} row {index} "
                    "has the wrong fields"
                )
            cadences.append(
                _manifest_integer(
                    row["cadenceno"],
                    context=(
                        f"cadence-reference authority exclusion {name} "
                        f"row {index} cadenceno"
                    ),
                )
            )
            _manifest_integer(
                row["spoc_quality"],
                context=(
                    f"cadence-reference authority exclusion {name} "
                    f"row {index} spoc_quality"
                ),
            )
        if len(cadences) != len(set(cadences)):
            raise ValueError(
                f"cadence-reference authority exclusion {name} has duplicate cadences"
            )
        declared_cadences = np.asarray(sorted(cadences), dtype=np.int64)
        retained_cadences = frame.loc[
            (frame["camera"] == camera) & (frame["ccd"] == ccd),
            "cadenceno",
        ].to_numpy(dtype=np.int64)
        overlap = np.intersect1d(
            declared_cadences,
            retained_cadences,
            assume_unique=True,
        )
        if len(overlap):
            raise ValueError(
                f"cadence-reference authority exclusion {name} is still present "
                f"in the table: {overlap[:20].astype(int).tolist()}"
            )
        output[(int(sector), camera, ccd)] = declared_cadences
        observed_total += len(declared_cadences)
    if observed_total != declared_total:
        raise ValueError(
            "cadence-reference authority-exclusion total disagrees with detector rows"
        )
    return output, observed_sha256, declared_total


def _validate_source_hash_declarations(manifest: Mapping[str, Any]) -> str:
    """Validate and fingerprint the authority hashes declared by the manifest.

    The manifest itself and the cadence table are rehashed before and after a
    BLS run.  This canonical fingerprint additionally binds the complete map
    of source-file hashes that the cadence-reference builder verified when it
    published the immutable evidence pair.
    """

    declared = manifest.get("source_file_sha256")
    if not isinstance(declared, Mapping) or not declared:
        raise ValueError("cadence-reference source_file_sha256 is missing or empty")
    normalized: dict[str, str] = {}
    for raw_path, raw_digest in declared.items():
        path = str(raw_path).strip()
        digest = str(raw_digest).lower()
        if not path:
            raise ValueError("cadence-reference source hash has an empty path")
        if len(digest) != 64 or any(char not in "0123456789abcdef" for char in digest):
            raise ValueError(f"cadence-reference source hash is invalid for {path!r}")
        if path in normalized:
            raise ValueError(f"duplicate cadence-reference source path: {path}")
        normalized[path] = digest

    sources = manifest.get("sources")
    if not isinstance(sources, list) or not sources:
        raise ValueError("cadence-reference sources inventory is missing or empty")
    observed_paths: set[str] = set()
    for index, source in enumerate(sources):
        if not isinstance(source, Mapping):
            raise ValueError(f"cadence-reference source {index} is not an object")
        path = str(source.get("path", "")).strip()
        digest = str(source.get("sha256", "")).lower()
        if not path or path in observed_paths:
            raise ValueError(
                "cadence-reference sources contain an empty/duplicate path"
            )
        observed_paths.add(path)
        if normalized.get(path) != digest:
            raise ValueError(
                "cadence-reference sources inventory disagrees with "
                f"source_file_sha256 for {path!r}"
            )
    if observed_paths != set(normalized):
        raise ValueError(
            "cadence-reference sources inventory does not exactly cover "
            "source_file_sha256"
        )
    return _canonical_json_sha256(normalized)


def load_external_quality_reference(
    *,
    table_path: Path,
    manifest_path: Path,
    sector: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Load and validate one authoritative detector-cadence quality overlay."""

    table_path = Path(table_path)
    manifest_path = Path(manifest_path)
    table_sha256 = _sha256(table_path)
    manifest_sha256 = _sha256(manifest_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(manifest, dict):
        raise ValueError("cadence-reference manifest must contain a JSON object")

    required_manifest = {
        "contract_version",
        "builder_version",
        "sector",
        "cadence_authority",
        "quality_authority",
        "quality_composition",
        "table_sha256",
        "n_rows",
        "detectors",
        "orbits",
        "source_file_sha256",
        "sources",
        "n_nonzero_spoc_quality",
        "n_nonzero_qlp_quality",
        "n_nonzero_external_quality",
        "n_spoc_rows_excluded_by_quat",
        "authority_exclusions",
        "authority_exclusions_sha256",
    }
    missing_manifest = sorted(required_manifest - set(manifest))
    if missing_manifest:
        raise ValueError(
            f"cadence-reference manifest is missing fields: {missing_manifest}"
        )
    expected_contract = CADENCE_REFERENCE_CONTRACT_TEMPLATE.format(sector=int(sector))
    if str(manifest["contract_version"]) != expected_contract:
        raise ValueError("cadence-reference manifest contract mismatch")
    if str(manifest["builder_version"]) != CADENCE_REFERENCE_BUILDER_VERSION:
        raise ValueError("cadence-reference manifest builder version mismatch")
    if int(manifest["sector"]) != int(sector):
        raise ValueError("cadence-reference manifest sector mismatch")
    if str(manifest["cadence_authority"]) != EXPECTED_CADENCE_AUTHORITY:
        raise ValueError("cadence-reference cadence authority mismatch")
    if str(manifest["quality_authority"]) != EXPECTED_QUALITY_AUTHORITY:
        raise ValueError("cadence-reference quality authority mismatch")
    if manifest["quality_composition"] != EXPECTED_QUALITY_COMPOSITION:
        raise ValueError("cadence-reference external-quality composition mismatch")
    if str(manifest["table_sha256"]) != table_sha256:
        raise ValueError("cadence-reference table hash mismatch")

    frame = _read_table(table_path)
    if tuple(frame.columns) != CADENCE_REFERENCE_COLUMNS:
        raise ValueError(
            "cadence-reference columns must be exactly "
            f"{CADENCE_REFERENCE_COLUMNS}; observed {tuple(frame.columns)}"
        )
    if len(frame) != int(manifest["n_rows"]):
        raise ValueError("cadence-reference row count disagrees with manifest")
    if frame.empty:
        raise ValueError("cadence-reference table is empty")
    for column in CADENCE_REFERENCE_COLUMNS:
        frame[column] = _integer_column(frame, column)
    if set(frame["sector"].astype(int)) != {int(sector)}:
        raise ValueError("cadence-reference table sector mismatch")
    if not frame["camera"].between(1, 4).all() or not frame["ccd"].between(1, 4).all():
        raise ValueError("cadence-reference camera/ccd values must be in 1..4")
    if (frame["orbitid"] <= 0).any() or (frame["cadenceno"] < 0).any():
        raise ValueError("cadence-reference orbit/cadence values are invalid")
    if (frame[["spoc_quality", "qlp_quality", "external_quality"]] < 0).any().any():
        raise ValueError("cadence-reference quality values must be nonnegative")
    if not set(frame["qlp_quality"].astype(int)).issubset({0, 1}):
        raise ValueError("cadence-reference qlp_quality must contain only 0/1")
    expected_external = np.bitwise_or(
        frame["spoc_quality"].to_numpy(dtype=np.int64),
        np.left_shift(frame["qlp_quality"].to_numpy(dtype=np.int64), 30),
    )
    if not np.array_equal(
        expected_external, frame["external_quality"].to_numpy(dtype=np.int64)
    ):
        raise ValueError("cadence-reference external_quality composition is invalid")

    key = ["sector", "camera", "ccd", "cadenceno"]
    duplicate = frame.duplicated(key, keep=False)
    if duplicate.any():
        examples = frame.loc[duplicate, [*key, "orbitid"]].head(10).to_dict("records")
        raise ValueError(
            f"cadence-reference has duplicate detector-cadence mappings: {examples}"
        )
    observed_detectors = sorted(
        {
            f"cam{int(camera)}_ccd{int(ccd)}"
            for camera, ccd in zip(frame["camera"], frame["ccd"], strict=True)
        }
    )
    if sorted(str(value) for value in manifest["detectors"]) != observed_detectors:
        raise ValueError("cadence-reference detector inventory mismatch")
    observed_orbits = sorted(set(frame["orbitid"].astype(int)))
    if sorted(int(value) for value in manifest["orbits"]) != observed_orbits:
        raise ValueError("cadence-reference orbit inventory mismatch")
    count_fields = {
        "n_nonzero_spoc_quality": "spoc_quality",
        "n_nonzero_qlp_quality": "qlp_quality",
        "n_nonzero_external_quality": "external_quality",
    }
    for manifest_field, column in count_fields.items():
        observed = int(np.count_nonzero(frame[column].to_numpy(dtype=np.int64)))
        if int(manifest[manifest_field]) != observed:
            raise ValueError(
                f"cadence-reference {manifest_field} disagrees with the table"
            )

    authority_exclusions, authority_exclusions_sha256, n_authority_exclusions = (
        _validate_authority_exclusions(
            manifest=manifest,
            frame=frame,
            sector=int(sector),
            detector_names=observed_detectors,
        )
    )
    source_hashes_sha256 = _validate_source_hash_declarations(manifest)
    frame = frame.sort_values(
        ["sector", "camera", "ccd", "cadenceno"], kind="stable"
    ).reset_index(drop=True)
    return frame, {
        "contract_version": expected_contract,
        "cadence_authority": EXPECTED_CADENCE_AUTHORITY,
        "quality_authority": EXPECTED_QUALITY_AUTHORITY,
        "table_sha256": table_sha256,
        "manifest_sha256": manifest_sha256,
        "source_file_sha256_declaration_sha256": source_hashes_sha256,
        "authority_exclusion_policy_contract": (AUTHORITY_EXCLUSION_POLICY_CONTRACT),
        "authority_exclusion_external_bit": AUTHORITY_EXCLUSION_EXTERNAL_BIT,
        "authority_exclusions_sha256": authority_exclusions_sha256,
        "n_authority_exclusions": n_authority_exclusions,
        "_authority_exclusions_by_detector": authority_exclusions,
    }


def _reference_worker_payload(
    frame: pd.DataFrame,
) -> dict[tuple[int, int, int], tuple[np.ndarray, np.ndarray, np.ndarray]]:
    payload: dict[tuple[int, int, int], tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for (sector, camera, ccd), detector in frame.groupby(
        ["sector", "camera", "ccd"], sort=True
    ):
        ordered = detector.sort_values("cadenceno", kind="stable")
        payload[(int(sector), int(camera), int(ccd))] = (
            ordered["cadenceno"].to_numpy(dtype=np.int64),
            ordered["orbitid"].to_numpy(dtype=np.int64),
            ordered["external_quality"].to_numpy(dtype=np.int64),
        )
    return payload


def _initialize_external_quality_worker(
    reference_by_detector: dict[
        tuple[int, int, int], tuple[np.ndarray, np.ndarray, np.ndarray]
    ],
    provenance: dict[str, Any],
    orbitid_policy: str = "strict",
) -> None:
    global _EXTERNAL_QUALITY_BY_DETECTOR
    global _EXTERNAL_QUALITY_AUTHORITY_EXCLUSIONS
    global _EXTERNAL_QUALITY_PROVENANCE
    global _ORBITID_POLICY
    exclusions = provenance.get("_authority_exclusions_by_detector")
    if not isinstance(exclusions, Mapping):
        raise ValueError(
            "external-quality provenance lacks validated authority exclusions"
        )
    _EXTERNAL_QUALITY_BY_DETECTOR = reference_by_detector
    _EXTERNAL_QUALITY_AUTHORITY_EXCLUSIONS = {
        key: np.asarray(values, dtype=np.int64) for key, values in exclusions.items()
    }
    _EXTERNAL_QUALITY_PROVENANCE = {
        key: value
        for key, value in provenance.items()
        if key != "_authority_exclusions_by_detector"
    }
    _ORBITID_POLICY = _validate_orbitid_policy(orbitid_policy)


def _apply_external_quality(lc: Any) -> dict[str, Any]:
    if (
        _EXTERNAL_QUALITY_BY_DETECTOR is None
        or _EXTERNAL_QUALITY_AUTHORITY_EXCLUSIONS is None
    ):
        raise RuntimeError("external-quality worker state was not initialized")
    lengths = {
        "time": len(lc.time),
        "cadenceno": len(lc.cadenceno),
        "orbitid": len(lc.orbitid),
        "quality": len(lc.quality),
    }
    if len(set(lengths.values())) != 1:
        raise ValueError(f"TIC {lc.tic}: compact cadence arrays disagree: {lengths}")
    cadences = np.asarray(lc.cadenceno, dtype=np.int64)
    orbits = np.asarray(lc.orbitid, dtype=np.int64).copy()
    internal_quality = np.asarray(lc.quality, dtype=np.int64)
    if len(cadences) != len(np.unique(cadences)):
        raise ValueError(f"TIC {lc.tic}: compact cadenceno values are not unique")
    detector_key = (int(lc.sector), int(lc.cam), int(lc.ccd))
    reference = _EXTERNAL_QUALITY_BY_DETECTOR.get(detector_key)
    if reference is None:
        raise ValueError(
            f"TIC {lc.tic}: no external-quality detector mapping for {detector_key}"
        )
    reference_cadence, reference_orbit, reference_quality = reference
    positions = np.searchsorted(reference_cadence, cadences)
    in_bounds = positions < len(reference_cadence)
    matched = np.zeros(len(cadences), dtype=bool)
    matched[in_bounds] = reference_cadence[positions[in_bounds]] == cadences[in_bounds]
    declared_exclusions = _EXTERNAL_QUALITY_AUTHORITY_EXCLUSIONS.get(
        detector_key,
        np.asarray([], dtype=np.int64),
    )
    authority_excluded = (~matched) & np.isin(
        cadences,
        declared_exclusions,
        assume_unique=True,
    )
    undeclared_missing = ~matched & ~authority_excluded
    if undeclared_missing.any():
        missing = cadences[undeclared_missing][:20].astype(int).tolist()
        raise ValueError(
            f"TIC {lc.tic}: external-quality coverage is missing compact cadences "
            f"for {detector_key}: {missing}"
        )
    mapped_orbits = np.zeros(len(cadences), dtype=np.int64)
    mapped_orbits[matched] = reference_orbit[positions[matched]]
    orbit_mismatch = matched & (mapped_orbits != orbits)
    corrections = sorted(
        [
            {
                "cadenceno": int(cadences[index]),
                "compact_orbitid": int(orbits[index]),
                "reference_orbitid": int(mapped_orbits[index]),
            }
            for index in np.flatnonzero(orbit_mismatch)
        ],
        key=lambda row: (
            row["cadenceno"],
            row["compact_orbitid"],
            row["reference_orbitid"],
        ),
    )
    if corrections and _ORBITID_POLICY == "strict":
        raise ValueError(
            f"TIC {lc.tic}: external-quality orbit mapping mismatch: {corrections[:20]}"
        )
    n_corrected = 0
    if _ORBITID_POLICY == "reference_by_cadence":
        reconciled_orbits = orbits.copy()
        reconciled_orbits[matched] = mapped_orbits[matched]
        lc.orbitid = reconciled_orbits
        n_corrected = len(corrections)
    correction_signature = _canonical_json_sha256(
        {
            "contract_version": ORBITID_RECONCILIATION_CONTRACT_VERSION,
            "policy": _ORBITID_POLICY,
            "tic": int(lc.tic),
            "sector": int(lc.sector),
            "camera": int(lc.cam),
            "ccd": int(lc.ccd),
            "corrections": corrections,
        }
    )
    external_quality = np.zeros(len(cadences), dtype=np.int64)
    external_quality[matched] = reference_quality[positions[matched]]
    external_quality[authority_excluded] = np.int64(
        1 << AUTHORITY_EXCLUSION_EXTERNAL_BIT
    )
    internal_bad = internal_quality != 0
    external_bad = external_quality != 0
    effective_bad = internal_bad | external_bad
    lc.quality = effective_bad.astype(np.int32)
    return {
        "n_cad_internal_bad": int(np.count_nonzero(internal_bad)),
        "n_cad_external_bad": int(np.count_nonzero(external_bad)),
        "n_cad_external_only_bad": int(np.count_nonzero(external_bad & ~internal_bad)),
        "n_cad_effective_bad": int(np.count_nonzero(effective_bad)),
        "n_cad_authority_excluded": int(np.count_nonzero(authority_excluded)),
        "n_cad_orbitid_reference_matched": int(np.count_nonzero(matched)),
        "n_cad_orbitid_mismatch": len(corrections),
        "n_cad_orbitid_corrected": n_corrected,
        "orbitid_correction_signature_sha256": correction_signature,
    }


def _target_tics(path: Path) -> list[int]:
    with h5py.File(path, "r") as h5:
        if "targets" not in h5:
            raise KeyError(f"compact export has no /targets group: {path}")
        return sorted(int(key) for key in h5["targets"].keys())


def _orbitid_summary_from_rows(
    frame: pd.DataFrame,
    *,
    orbitid_policy: str,
) -> dict[str, Any]:
    """Validate per-target reconciliation fields and fingerprint the shard."""

    policy = _validate_orbitid_policy(orbitid_policy)
    count_columns = (
        "n_cad_orbitid_reference_matched",
        "n_cad_orbitid_mismatch",
        "n_cad_orbitid_corrected",
    )
    required = {
        "tic",
        "orbitid_policy",
        "orbitid_reconciliation_contract_version",
        "orbitid_correction_signature_sha256",
        *count_columns,
    }
    missing = sorted(required - set(frame))
    if missing:
        raise ValueError(f"BLS rows lack orbit-ID reconciliation fields: {missing}")
    if set(frame["orbitid_policy"].astype(str)) != {policy}:
        raise ValueError("BLS rows disagree with the requested orbit-ID policy")
    if set(frame["orbitid_reconciliation_contract_version"].astype(str)) != {
        ORBITID_RECONCILIATION_CONTRACT_VERSION
    }:
        raise ValueError("BLS rows have the wrong orbit-ID reconciliation contract")
    targets = frame.drop_duplicates("tic", keep="first").copy()
    if (
        frame.groupby("tic", sort=False)[
            [*count_columns, "orbitid_correction_signature_sha256"]
        ]
        .nunique(dropna=False)
        .gt(1)
        .any()
        .any()
    ):
        raise ValueError("BLS orbit-ID reconciliation fields disagree within a TIC")
    for column in count_columns:
        values = pd.to_numeric(targets[column], errors="coerce")
        if (
            values.isna().any()
            or (values < 0).any()
            or (values != np.floor(values)).any()
        ):
            raise ValueError(f"BLS rows contain invalid {column}")
        targets[column] = values.astype(np.int64)
    if (
        targets["n_cad_orbitid_mismatch"] > targets["n_cad_orbitid_reference_matched"]
    ).any():
        raise ValueError("BLS orbit-ID mismatch count exceeds matched cadences")
    if policy == "strict":
        if (
            targets[["n_cad_orbitid_mismatch", "n_cad_orbitid_corrected"]].to_numpy(
                dtype=np.int64
            )
            != 0
        ).any():
            raise ValueError("strict orbit-ID policy cannot publish corrections")
    elif (
        not targets["n_cad_orbitid_corrected"]
        .eq(targets["n_cad_orbitid_mismatch"])
        .all()
    ):
        raise ValueError(
            "reference_by_cadence must correct every matched orbit-ID mismatch"
        )
    signatures = targets["orbitid_correction_signature_sha256"].astype(str)
    if not signatures.str.fullmatch(r"[0-9a-f]{64}").all():
        raise ValueError("BLS rows contain invalid orbit-ID correction signatures")
    records = [
        {
            "tic": int(row.tic),
            "n_cad_orbitid_reference_matched": int(row.n_cad_orbitid_reference_matched),
            "n_cad_orbitid_mismatch": int(row.n_cad_orbitid_mismatch),
            "n_cad_orbitid_corrected": int(row.n_cad_orbitid_corrected),
            "orbitid_correction_signature_sha256": str(
                row.orbitid_correction_signature_sha256
            ),
        }
        for row in targets.sort_values("tic", kind="stable").itertuples(index=False)
    ]
    return {
        "n_cad_orbitid_reference_matched": int(
            targets["n_cad_orbitid_reference_matched"].sum()
        ),
        "n_cad_orbitid_mismatch": int(targets["n_cad_orbitid_mismatch"].sum()),
        "n_cad_orbitid_corrected": int(targets["n_cad_orbitid_corrected"].sum()),
        "n_targets_orbitid_mismatch": int(
            targets["n_cad_orbitid_mismatch"].gt(0).sum()
        ),
        "orbitid_corrections_sha256": _canonical_json_sha256(
            {
                "contract_version": ORBITID_RECONCILIATION_CONTRACT_VERSION,
                "policy": policy,
                "targets": records,
            }
        ),
    }


def _result_rows(result: Any) -> list[dict[str, Any]]:
    base = {
        "tic": int(result.tic),
        "sector": int(result.sector),
        "cam": int(result.cam),
        "ccd": int(result.ccd),
        "tmag": float(result.tmag),
        "aperture": str(result.aperture),
        "n_cad_total": int(result.n_cad_total),
        "n_cad_quality": int(result.n_cad_quality or 0),
        "n_cad_kept": int(result.n_cad_kept),
        "n_cad_edge_trimmed": int(result.n_cad_edge_trimmed),
        "n_cad_sigma_clipped": int(result.n_cad_sigma_clipped),
        "dropout_frac": float(result.dropout_frac),
        "quality_dropout_frac": float(result.quality_dropout_frac or 0.0),
        "n_orbits": int(result.n_orbits),
        "baseline_d": float(result.baseline_d),
        "status": str(result.status),
        "bls_search_branch": "current_adp",
        "adp_only_contract_version": ADP_ONLY_CONTRACT_VERSION,
    }
    if not result.peaks:
        return [
            {
                **base,
                "peak_rank": 0,
                "period_d": np.nan,
                "t0_bjd": np.nan,
                "duration_min": np.nan,
                "depth": np.nan,
                "depth_snr": np.nan,
                "sde": np.nan,
                "log_power": np.nan,
            }
        ]
    return [{**base, **asdict(peak)} for peak in result.peaks]


def _raw_v4_light_curve(
    *,
    raw_source_h5: Path,
    tic: int,
    sector: int,
    camera: int,
    ccd: int,
    tessmag: float,
    compact_lc: Path,
    expected_n_cadences: int,
) -> tuple[HLSPLightCurve, dict[str, Any]]:
    """Build detector-consistent ADP channels from one immutable raw-v1 group."""

    key = f"{int(tic):016d}"
    with h5py.File(raw_source_h5, "r") as h5:
        if "targets" not in h5 or key not in h5["targets"]:
            raise KeyError(f"raw-v1 source lacks targets/{key}")
        group = h5["targets"][key]
        for name, expected in (
            ("sector", int(sector)),
            ("tic", int(tic)),
            ("camera", int(camera)),
            ("ccd", int(ccd)),
        ):
            if int(group.attrs.get(name, -1)) != expected:
                raise ValueError(
                    f"TIC {tic}: raw-v1 {name} mapping differs from the "
                    "frozen pool"
                )
        required = (
            "time",
            "cadenceno",
            "orbitid",
            "quality",
            "raw_flux_small",
            "raw_flux_err_small",
            "raw_flux_primary",
            "raw_flux_err_primary",
        )
        missing = [name for name in required if name not in group]
        if missing:
            raise KeyError(f"TIC {tic}: raw-v1 group lacks {missing}")
        arrays = {name: np.asarray(group[name]) for name in required}

    lengths = {name: len(values) for name, values in arrays.items()}
    if len(set(lengths.values())) != 1 or not lengths or next(iter(lengths.values())) < 1:
        raise ValueError(f"TIC {tic}: raw-v1 arrays differ in length: {lengths}")
    if next(iter(lengths.values())) != int(expected_n_cadences):
        raise ValueError(
            f"TIC {tic}: raw-v1 cadence count differs from the frozen pool: "
            f"{next(iter(lengths.values()))} != {expected_n_cadences}"
        )
    cadenceno = np.asarray(arrays["cadenceno"], dtype=np.int64)
    if len(np.unique(cadenceno)) != len(cadenceno):
        raise ValueError(f"TIC {tic}: raw-v1 cadences are not unique")
    if np.any(np.diff(cadenceno) <= 0):
        raise ValueError(f"TIC {tic}: raw-v1 cadences are not strictly sorted")
    time_bjd = np.asarray(arrays["time"], dtype=np.float64)
    if not np.isfinite(time_bjd).all():
        raise ValueError(f"TIC {tic}: raw-v1 time contains nonfinite values")
    time_btjd = time_bjd - float(BJDREFI)
    if np.any(time_btjd < 0.0) or np.any(time_btjd >= 100_000.0):
        raise ValueError(f"TIC {tic}: raw-v1 time is not absolute BJD")
    with h5py.File(compact_lc, "r") as compact_file:
        compact_path = f"targets/{int(tic):016d}"
        if compact_path not in compact_file:
            raise KeyError(f"TIC {tic}: frozen compact cadence authority is absent")
        compact_group = compact_file[compact_path]
        compact_cadence = np.asarray(compact_group["cadenceno"], dtype=np.int64)
        compact_time = np.asarray(compact_group["time"], dtype=np.float64)
    if not np.array_equal(cadenceno, compact_cadence):
        raise ValueError(
            f"TIC {tic}: raw-v1 cadence inventory differs from frozen compact"
        )
    if len(compact_time) != len(time_btjd):
        raise ValueError(f"TIC {tic}: frozen compact time length differs")
    compact_time_bjd = np.where(
        compact_time < 100_000.0,
        compact_time + float(BJDREFI),
        compact_time,
    )
    time_delta_s = np.abs(time_bjd - compact_time_bjd) * 86_400.0
    if not np.isfinite(time_delta_s).all() or float(np.max(time_delta_s)) > 2.0:
        raise ValueError(
            f"TIC {tic}: raw-v1/frozen-compact times differ by more than 2 s"
        )
    lc = HLSPLightCurve(
        tic=int(tic),
        tmag=float(tessmag),
        sector=int(sector),
        cam=int(camera),
        ccd=int(ccd),
        ra=float("nan"),
        dec=float("nan"),
        time=time_btjd,
        cadenceno=cadenceno,
        orbitid=np.asarray(arrays["orbitid"], dtype=np.int32),
        quality=np.asarray(arrays["quality"], dtype=np.int32),
        flux={},
        path=Path(f"{raw_source_h5}:targets/{key}"),
    )
    quality_counts = _apply_external_quality(lc)
    quality_counts.update(
        {
            "raw_compact_cadence_inventory_match": 1,
            "raw_compact_time_delta_max_s": float(np.max(time_delta_s)),
        }
    )
    effective_quality = np.asarray(lc.quality, dtype=np.int32)

    def rebuild_or_mark_unsearchable(
        raw_flux: np.ndarray,
        raw_error: np.ndarray,
    ) -> np.ndarray:
        flux = np.asarray(raw_flux, dtype=np.float64)
        error = np.asarray(raw_error, dtype=np.float64)
        finite_good = (
            (effective_quality == 0)
            & np.isfinite(flux)
            & np.isfinite(error)
            & (error > 0)
        )
        if not finite_good.any():
            # Preserve the observation in the BLS authority.  The downstream
            # search emits its locked rank-0 too_few_cadences status row; no
            # uncertainty or flux value is invented for this aperture.
            return np.full(len(flux), np.nan, dtype=np.float64)
        rebuilt, _ = _effective_quality_adp03q(
            time_btjd=time_btjd,
            raw_flux=flux,
            raw_error=error,
            quality=effective_quality,
        )
        return rebuilt

    small = rebuild_or_mark_unsearchable(
        arrays["raw_flux_small"], arrays["raw_flux_err_small"]
    )
    primary = rebuild_or_mark_unsearchable(
        arrays["raw_flux_primary"], arrays["raw_flux_err_primary"]
    )
    lc.flux = {
        "DET_FLUX_ADP_SML": small,
        "DET_FLUX_ADP": primary,
    }
    return lc, quality_counts


def _process_target(payload: tuple[Any, ...]) -> list[dict[str, Any]]:
    if len(payload) == 3:
        tic, compact_lc_s, cfg_payload = payload
        raw_payload = None
    else:
        tic, compact_lc_s, cfg_payload, raw_payload = payload
    compact_lc = Path(compact_lc_s)
    cfg = BLSConfig(
        apertures=tuple(cfg_payload.get("apertures", ADP_ONLY_APERTURES)),
        n_periods=int(cfg_payload["n_periods"]),
        n_peaks=int(cfg_payload["n_peaks"]),
        p_min_d=float(cfg_payload["p_min_d"]),
        p_max_cap_d=float(cfg_payload["p_max_cap_d"]),
        max_period_fraction=float(cfg_payload["max_period_fraction"]),
        durations_min=tuple(
            float(value)
            for value in cfg_payload.get("durations_min", BLSConfig.durations_min)
        ),
        period_mask_frac=float(
            cfg_payload.get("period_mask_frac", BLSConfig.period_mask_frac)
        ),
        period_bin_edges=tuple(
            float(value) for value in cfg_payload.get("period_bin_edges", ())
        ),
        max_peaks_per_period_bin=int(
            cfg_payload.get(
                "max_peaks_per_period_bin", BLSConfig.max_peaks_per_period_bin
            )
        ),
        min_cadences=int(cfg_payload.get("min_cadences", BLSConfig.min_cadences)),
        sigma_clip=float(cfg_payload["sigma_clip"]),
        orbit_edge_trim_d=float(cfg_payload["orbit_edge_trim_d"]),
    )
    config_sha256 = bls_config_sha256(cfg_payload)
    search_contract = (
        A2V1_TEACHER_BLS_SEARCH_CONTRACT
        if cfg_payload == approved_a2v1_teacher_bls_config()
        else "custom"
    )
    config_provenance = {
        "bls_n_periods": int(cfg.n_periods),
        "bls_n_peaks": int(cfg.n_peaks),
        "bls_p_min_d": float(cfg.p_min_d),
        "bls_p_max_cap_d": float(cfg.p_max_cap_d),
        "bls_max_period_fraction": float(cfg.max_period_fraction),
        "bls_sigma_clip": float(cfg.sigma_clip),
        "bls_orbit_edge_trim_d": float(cfg.orbit_edge_trim_d),
    }
    if raw_payload is None:
        lc = read_compact_lc_export(
            compact_lc, tic=tic, columns=ADP_ONLY_APERTURES
        )
    else:
        lc, quality_counts = _raw_v4_light_curve(
            raw_source_h5=Path(raw_payload["path"]),
            tic=int(tic),
            sector=int(raw_payload["sector"]),
            camera=int(raw_payload["camera"]),
            ccd=int(raw_payload["ccd"]),
            tessmag=float(raw_payload["tessmag"]),
            compact_lc=compact_lc,
            expected_n_cadences=int(raw_payload["n_cadences"]),
        )
    if lc is None:
        raise RuntimeError(
            f"TIC {tic}: compact product could not supply the locked ADP pair"
        )
    if raw_payload is None:
        quality_counts = _apply_external_quality(lc)
    if _EXTERNAL_QUALITY_PROVENANCE is None:
        raise RuntimeError("external-quality provenance was not initialized")
    rows: list[dict[str, Any]] = []
    rank_one_peaks: dict[str, dict[str, Any]] = {}
    for aperture in ADP_ONLY_APERTURES:
        result = run_bls_on_lc(lc, cfg, aperture=aperture)
        current = _result_rows(result)
        rank_one = next(
            (
                row
                for row in current
                if int(row.get("peak_rank", 0)) == 1
                and str(row.get("status", "")) == "ok"
            ),
            None,
        )
        if rank_one is not None:
            rank_one_peaks[aperture] = dict(rank_one)
        for row in current:
            row["source_product_tag"] = str(cfg_payload.get("source_product_tag", ""))
            row.update(config_provenance)
            row["bls_search_contract_version"] = search_contract
            row["bls_config_sha256"] = config_sha256
            row.update(quality_counts)
            row["external_quality_policy_contract"] = EXTERNAL_QUALITY_POLICY_CONTRACT
            row["cadence_reference_sha256"] = _EXTERNAL_QUALITY_PROVENANCE[
                "table_sha256"
            ]
            row["cadence_reference_manifest_sha256"] = _EXTERNAL_QUALITY_PROVENANCE[
                "manifest_sha256"
            ]
            row["authority_exclusion_policy_contract"] = _EXTERNAL_QUALITY_PROVENANCE[
                "authority_exclusion_policy_contract"
            ]
            row["authority_exclusion_external_bit"] = int(
                _EXTERNAL_QUALITY_PROVENANCE["authority_exclusion_external_bit"]
            )
            row["authority_exclusions_sha256"] = _EXTERNAL_QUALITY_PROVENANCE[
                "authority_exclusions_sha256"
            ]
            row["n_authority_exclusions"] = int(
                _EXTERNAL_QUALITY_PROVENANCE["n_authority_exclusions"]
            )
            row["orbitid_policy"] = _ORBITID_POLICY
            row["orbitid_reconciliation_contract_version"] = (
                ORBITID_RECONCILIATION_CONTRACT_VERSION
            )
            row["photometry_source_contract_version"] = (
                RAW_V4_INPUT_CONTRACT_VERSION if raw_payload is not None else "compact"
            )
            row["detrend_contract_version"] = (
                FULL_POOL_NATIVE_DETREND_CONTRACT_VERSION
                if raw_payload is not None
                else "precomputed_compact"
            )
            row["detrend_config_sha256"] = (
                FULL_POOL_NATIVE_DETREND_CONFIG_SHA256
                if raw_payload is not None
                else ""
            )
        rows.extend(current)
    if int(lc.sector) == 63:
        target_metadata = {
            column: float("nan") for column in S63_TARGET_METADATA_COLUMNS
        }
        anchor = rank_one_peaks.get(ADP_ONLY_APERTURES[0])
        if anchor is not None:
            measured = measure_two_aperture_candidate_metadata(
                lc,
                anchor_peak=anchor,
                own_peaks=rank_one_peaks,
                apertures=ADP_ONLY_APERTURES,
            )
            target_metadata.update(
                {
                    column: measured.get(column, float("nan"))
                    for column in S63_TARGET_METADATA_COLUMNS
                }
            )
        for row in rows:
            row.update(target_metadata)
            row["target_metadata_contract_version"] = (
                S63_TARGET_METADATA_CONTRACT_VERSION
            )
    return rows


def build_peak_table(
    *,
    sector: int = 56,
    compact_lc: Path,
    cadence_reference: Path,
    cadence_reference_manifest: Path,
    out_dir: Path,
    workers: int,
    n_periods: int,
    n_peaks: int,
    max_targets: int | None,
    progress_every: int,
    shard_index: int = 0,
    n_shards: int = 1,
    resume: bool = False,
    source_product_tag: str = "",
    target_allowlist: Path | None = None,
    orbitid_policy: str = "strict",
    expected_compact_lc_sha256: str | None = None,
    expected_target_allowlist_sha256: str | None = None,
    raw_source_h5: Path | None = None,
    raw_source_summary: Path | None = None,
    raw_export_complete: Path | None = None,
    raw_transfer_validation: Path | None = None,
    frozen_pool: Path | None = None,
    frozen_pool_summary: Path | None = None,
    execution_allowlist: Path | None = None,
) -> dict[str, Any]:
    validate_adp_only_apertures(ADP_ONLY_APERTURES)
    orbitid_policy = _validate_orbitid_policy(orbitid_policy)
    expected_compact_lc_sha256 = _validate_expected_sha256(
        expected_compact_lc_sha256,
        context="expected_compact_lc_sha256",
    )
    expected_target_allowlist_sha256 = _validate_expected_sha256(
        expected_target_allowlist_sha256,
        context="expected_target_allowlist_sha256",
    )
    compact_lc = Path(compact_lc)
    cadence_reference = Path(cadence_reference)
    cadence_reference_manifest = Path(cadence_reference_manifest)
    raw_mode = raw_source_h5 is not None
    if raw_mode:
        required_raw_inputs = {
            "raw_source_summary": raw_source_summary,
            "raw_export_complete": raw_export_complete,
            "raw_transfer_validation": raw_transfer_validation,
            "frozen_pool": frozen_pool,
            "frozen_pool_summary": frozen_pool_summary,
            "target_allowlist": target_allowlist,
        }
        missing_raw = [name for name, value in required_raw_inputs.items() if value is None]
        if missing_raw:
            raise ValueError(f"raw-v4 mode requires inputs: {missing_raw}")
        if max_targets is not None:
            raise ValueError("raw-v4 mode cannot use max_targets")
        authority = load_sector_pool_authority(
            sector=int(sector),
            pool_path=Path(frozen_pool),
            pool_summary_path=Path(frozen_pool_summary),
            allowlist_path=Path(target_allowlist),
        )
        source_rows = select_sector_shard(
            authority,
            shard_index=int(shard_index),
            n_shards=int(n_shards),
        )
        raw_release = load_production_raw_source_release(
            authority=authority,
            source_rows=source_rows,
            sector=int(sector),
            shard_index=int(shard_index),
            n_shards=int(n_shards),
            raw_source_h5=Path(raw_source_h5),
            raw_source_summary_path=Path(raw_source_summary),
            raw_export_complete_path=Path(raw_export_complete),
            raw_transfer_validation_path=Path(raw_transfer_validation),
        )
        if expected_compact_lc_sha256 != authority.compact_h5_sha256:
            raise ValueError(
                "raw-v4 compact lineage digest differs from the frozen pool"
            )
        raw_source_h5 = raw_release.raw_source.path
        compact_sha256 = authority.compact_h5_sha256
    else:
        authority = None
        source_rows = None
        raw_release = None
        compact_sha256 = _sha256(compact_lc)
    input_sha256 = {
        "compact_lc": compact_sha256,
        "cadence_reference": _sha256(cadence_reference),
        "cadence_reference_manifest": _sha256(cadence_reference_manifest),
    }
    if raw_mode and raw_release is not None:
        input_sha256.update(
            {
                "raw_source_h5": raw_release.raw_source.sha256,
                "raw_source_summary": raw_release.raw_source_summary.sha256,
                "raw_export_complete": raw_release.export_complete.sha256,
                "raw_transfer_validation": raw_release.transfer_validation.sha256,
            }
        )
    if (
        expected_compact_lc_sha256 is not None
        and input_sha256["compact_lc"] != expected_compact_lc_sha256
    ):
        raise RuntimeError(
            "compact light-curve SHA-256 disagrees with the frozen input "
            "contract"
        )
    if target_allowlist is not None:
        target_allowlist = Path(target_allowlist)
        allowlist_tics, target_selection = load_target_allowlist(target_allowlist)
        input_sha256["target_allowlist"] = target_selection["target_allowlist_sha256"]
        if max_targets is not None:
            raise ValueError(
                "max_targets cannot be combined with an exact target allowlist"
            )
        if (
            expected_target_allowlist_sha256 is not None
            and target_selection["target_allowlist_sha256"]
            != expected_target_allowlist_sha256
        ):
            raise RuntimeError(
                "target-allowlist SHA-256 disagrees with the frozen input "
                "contract"
            )
    else:
        allowlist_tics = None
        if expected_target_allowlist_sha256 is not None:
            raise ValueError(
                "expected_target_allowlist_sha256 requires target_allowlist"
            )
    if execution_allowlist is not None:
        if not raw_mode:
            raise ValueError("execution_allowlist is bounded to raw-v4 mode")
        execution_allowlist = Path(execution_allowlist)
        input_sha256["execution_allowlist"] = _sha256(execution_allowlist)
    reference, reference_provenance = load_external_quality_reference(
        table_path=cadence_reference,
        manifest_path=cadence_reference_manifest,
        sector=int(sector),
    )
    if reference_provenance["table_sha256"] != input_sha256["cadence_reference"]:
        raise RuntimeError("cadence-reference hash changed while it was validated")
    if (
        reference_provenance["manifest_sha256"]
        != input_sha256["cadence_reference_manifest"]
    ):
        raise RuntimeError("cadence-reference manifest changed while it was validated")
    reference_payload = _reference_worker_payload(reference)
    out_dir.mkdir(parents=True, exist_ok=True)
    compact_tics = (
        _target_tics(Path(raw_source_h5)) if raw_mode else _target_tics(compact_lc)
    )
    if allowlist_tics is None:
        tics = compact_tics
        if max_targets is not None:
            tics = tics[: max(0, int(max_targets))]
        target_selection = {
            "target_selection_contract_version": (TARGET_SELECTION_CONTRACT_VERSION),
            "target_allowlist": None,
            "target_allowlist_sha256": None,
            "target_allowlist_count": len(tics),
            "target_allowlist_tics_sha256": _tic_inventory_sha256(tics),
        }
    elif raw_mode:
        assert source_rows is not None
        expected_shard_tics = source_rows["tic"].astype(np.int64).tolist()
        if compact_tics != expected_shard_tics:
            raise ValueError(
                "raw-v4 target inventory differs from the authenticated raw shard"
            )
        tics = compact_tics
        n_targets_total = int(len(authority.rows))
        if execution_allowlist is not None:
            execution_tics, execution_selection = load_target_allowlist(
                Path(execution_allowlist)
            )
            outside_sector = sorted(set(execution_tics) - set(allowlist_tics))
            if outside_sector:
                raise ValueError(
                    "raw-v4 execution allowlist lies outside the sector pool; "
                    f"first={outside_sector[:10]}"
                )
            tics = sorted(set(tics) & set(execution_tics))
            if not tics:
                raise ValueError(
                    "raw-v4 execution allowlist has no targets in this raw shard"
                )
            target_selection = execution_selection
            n_targets_total = len(execution_tics)
    else:
        compact_set = set(compact_tics)
        missing_allowlist = sorted(set(allowlist_tics) - compact_set)
        if missing_allowlist:
            raise ValueError(
                "target allowlist is not exactly covered by the compact export; "
                f"missing_count={len(missing_allowlist)}, "
                f"examples={missing_allowlist[:20]}"
            )
        tics = allowlist_tics
    if n_shards < 1 or shard_index < 0 or shard_index >= n_shards:
        raise ValueError("shard_index must satisfy 0 <= shard_index < n_shards")
    if not raw_mode:
        n_targets_total = len(tics)
        tics = tics[int(shard_index) :: int(n_shards)]
    suffix = f"_{int(shard_index):03d}" if int(n_shards) > 1 else ""
    output_path = out_dir / f"real_adp_bls_peaks{suffix}.parquet"
    summary_path = out_dir / f"summary{suffix}.json"
    cfg_payload = approved_a2v1_teacher_bls_config()
    cfg_payload.update(
        {
            "n_periods": int(n_periods),
            "n_peaks": int(n_peaks),
            "source_product_tag": str(source_product_tag),
        }
    )
    cfg_sha256 = bls_config_sha256(cfg_payload)
    search_contract = (
        A2V1_TEACHER_BLS_SEARCH_CONTRACT
        if cfg_payload == approved_a2v1_teacher_bls_config()
        else "custom"
    )
    if resume and output_path.exists() and summary_path.exists():
        summary = json.loads(summary_path.read_text())
        if (
            summary.get("contract_version") == ADP_ONLY_CONTRACT_VERSION
            and summary.get("external_quality_policy_contract")
            == EXTERNAL_QUALITY_POLICY_CONTRACT
            and summary.get("compact_lc") == str(compact_lc)
            and summary.get("compact_lc_sha256") == input_sha256["compact_lc"]
            and summary.get("cadence_reference") == str(cadence_reference)
            and summary.get("cadence_reference_sha256")
            == input_sha256["cadence_reference"]
            and summary.get("cadence_reference_manifest")
            == str(cadence_reference_manifest)
            and summary.get("cadence_reference_manifest_sha256")
            == input_sha256["cadence_reference_manifest"]
            and summary.get("cadence_reference_contract_version")
            == reference_provenance["contract_version"]
            and summary.get("cadence_reference_cadence_authority")
            == reference_provenance["cadence_authority"]
            and summary.get("cadence_reference_quality_authority")
            == reference_provenance["quality_authority"]
            and summary.get("cadence_reference_source_hashes_sha256")
            == reference_provenance["source_file_sha256_declaration_sha256"]
            and summary.get("authority_exclusion_policy_contract")
            == reference_provenance["authority_exclusion_policy_contract"]
            and int(summary.get("authority_exclusion_external_bit", -1))
            == int(reference_provenance["authority_exclusion_external_bit"])
            and summary.get("authority_exclusions_sha256")
            == reference_provenance["authority_exclusions_sha256"]
            and int(summary.get("n_authority_exclusions", -1))
            == int(reference_provenance["n_authority_exclusions"])
            and int(summary.get("shard_index", -1)) == int(shard_index)
            and int(summary.get("n_shards", -1)) == int(n_shards)
            and int(summary.get("n_targets", -1)) == len(tics)
            and int(summary.get("n_targets_total", -1)) == n_targets_total
            and summary.get("source_product_tag") == str(source_product_tag)
            and summary.get("config") == cfg_payload
            and summary.get("bls_config_sha256") == cfg_sha256
            and summary.get("bls_search_contract_version") == search_contract
            and summary.get("target_selection_contract_version")
            == target_selection["target_selection_contract_version"]
            and summary.get("target_allowlist") == target_selection["target_allowlist"]
            and summary.get("target_allowlist_sha256")
            == target_selection["target_allowlist_sha256"]
            and int(summary.get("target_allowlist_count", -1))
            == int(target_selection["target_allowlist_count"])
            and summary.get("target_allowlist_tics_sha256")
            == target_selection["target_allowlist_tics_sha256"]
            and summary.get("orbitid_policy") == orbitid_policy
            and summary.get("orbitid_reconciliation_contract_version")
            == ORBITID_RECONCILIATION_CONTRACT_VERSION
            and summary.get("input_mode")
            == ("immutable_raw_v1_detector_consistent" if raw_mode else "compact")
            and summary.get("raw_source_h5_sha256")
            == input_sha256.get("raw_source_h5")
            and summary.get("execution_allowlist_sha256")
            == input_sha256.get("execution_allowlist")
            and summary.get("peak_table_sha256") == _sha256(output_path)
        ):
            print(json.dumps(summary, indent=2, sort_keys=True))
            return summary
    if raw_mode:
        assert source_rows is not None and raw_source_h5 is not None
        metadata = {
            int(row.tic): {
                "path": str(raw_source_h5),
                "sector": int(row.sector),
                "camera": int(row.camera),
                "ccd": int(row.ccd),
                "tessmag": float(row.tessmag),
                "n_cadences": int(row.n_cadences),
            }
            for row in source_rows.itertuples(index=False)
        }
        payloads = [
            (tic, str(compact_lc), cfg_payload, metadata[int(tic)]) for tic in tics
        ]
    else:
        payloads = [(tic, str(compact_lc), cfg_payload) for tic in tics]
    rows: list[dict[str, Any]] = []
    workers = max(1, int(workers))
    if workers == 1:
        _initialize_external_quality_worker(
            reference_payload,
            reference_provenance,
            orbitid_policy,
        )
        iterator = map(_process_target, payloads)
        executor = None
    else:
        executor = ProcessPoolExecutor(
            max_workers=workers,
            initializer=_initialize_external_quality_worker,
            initargs=(reference_payload, reference_provenance, orbitid_policy),
        )
        iterator = executor.map(_process_target, payloads, chunksize=1)
    try:
        for index, batch in enumerate(iterator, start=1):
            rows.extend(batch)
            if progress_every > 0 and index % int(progress_every) == 0:
                print(
                    f"[adp-real-bls] processed {index:,}/{len(payloads):,}", flush=True
                )
    finally:
        if executor is not None:
            executor.shutdown(wait=True)

    final_input_sha256 = {
        "compact_lc": (
            compact_sha256 if raw_mode else _sha256(compact_lc)
        ),
        "cadence_reference": _sha256(cadence_reference),
        "cadence_reference_manifest": _sha256(cadence_reference_manifest),
    }
    if raw_mode and raw_release is not None:
        for binding in (
            raw_release.raw_source,
            raw_release.raw_source_summary,
            raw_release.export_complete,
            raw_release.transfer_validation,
        ):
            binding.assert_unchanged()
        final_input_sha256.update(
            {
                "raw_source_h5": raw_release.raw_source.sha256,
                "raw_source_summary": raw_release.raw_source_summary.sha256,
                "raw_export_complete": raw_release.export_complete.sha256,
                "raw_transfer_validation": raw_release.transfer_validation.sha256,
            }
        )
    if target_allowlist is not None:
        final_input_sha256["target_allowlist"] = _sha256(target_allowlist)
    if execution_allowlist is not None:
        final_input_sha256["execution_allowlist"] = _sha256(execution_allowlist)
    if final_input_sha256 != input_sha256:
        raise RuntimeError(
            "BLS input source changed during the run; refusing to publish peaks"
        )

    peaks = pd.DataFrame(rows)
    if "tic" not in peaks:
        raise ValueError("BLS produced no target rows")
    observed_tics = set(
        pd.to_numeric(peaks["tic"], errors="raise").astype(np.int64).tolist()
    )
    target_metadata_summary: dict[str, Any] = {}
    if int(sector) == 63:
        missing_metadata = sorted(set(S63_TARGET_METADATA_COLUMNS) - set(peaks))
        if missing_metadata:
            raise ValueError(f"S63 BLS rows lack target metadata: {missing_metadata}")
        if set(peaks["target_metadata_contract_version"].astype(str)) != {
            S63_TARGET_METADATA_CONTRACT_VERSION
        }:
            raise ValueError("S63 BLS target metadata contract mismatch")
        target_metadata = peaks.drop_duplicates("tic", keep="first")
        finite_counts = {
            column: int(
                np.isfinite(
                    pd.to_numeric(target_metadata[column], errors="coerce").to_numpy(
                        dtype=float
                    )
                ).sum()
            )
            for column in S63_TARGET_METADATA_COLUMNS
        }
        if any(count == 0 for count in finite_counts.values()):
            raise ValueError(
                "S63 BLS target metadata contains a wholly nonfinite feature: "
                f"{finite_counts}"
            )
        target_metadata_summary = {
            "target_metadata_contract_version": (
                S63_TARGET_METADATA_CONTRACT_VERSION
            ),
            "target_metadata_columns": list(S63_TARGET_METADATA_COLUMNS),
            "target_metadata_finite_counts": finite_counts,
        }
    expected_tics = set(tics)
    if observed_tics != expected_tics:
        raise ValueError(
            "BLS target coverage does not exactly match the selected inventory; "
            f"missing={sorted(expected_tics - observed_tics)[:20]}, "
            f"unexpected={sorted(observed_tics - expected_tics)[:20]}"
        )
    output_path = write_table(peaks, output_path)
    output_sha256 = _sha256(output_path)
    status = peaks.get("status", pd.Series(dtype=str)).fillna("").astype(str)
    valid = status.eq("ok") & pd.to_numeric(peaks.get("peak_rank"), errors="coerce").gt(
        0
    )
    quality_count_columns = (
        "n_cad_internal_bad",
        "n_cad_external_bad",
        "n_cad_external_only_bad",
        "n_cad_effective_bad",
        "n_cad_authority_excluded",
    )
    quality_counts_over_unique_targets = {column: 0 for column in quality_count_columns}
    if {"tic", *quality_count_columns}.issubset(peaks.columns):
        target_quality = peaks.drop_duplicates("tic", keep="first")
        quality_counts_over_unique_targets = {
            column: int(pd.to_numeric(target_quality[column], errors="raise").sum())
            for column in quality_count_columns
        }
    raw_cadence_audit: dict[str, Any] = {}
    if raw_mode:
        required_audit = {
            "raw_compact_cadence_inventory_match",
            "raw_compact_time_delta_max_s",
        }
        missing_audit = sorted(required_audit - set(peaks.columns))
        if missing_audit:
            raise ValueError(f"raw-v4 BLS rows lack cadence audit: {missing_audit}")
        target_audit = peaks.drop_duplicates("tic", keep="first")
        inventory_match = pd.to_numeric(
            target_audit["raw_compact_cadence_inventory_match"], errors="raise"
        )
        if not inventory_match.eq(1).all():
            raise ValueError("raw-v4 cadence inventory audit did not pass")
        delta = pd.to_numeric(
            target_audit["raw_compact_time_delta_max_s"], errors="raise"
        )
        if not np.isfinite(delta.to_numpy(dtype=float)).all() or delta.gt(2.0).any():
            raise ValueError("raw-v4 cadence time audit exceeded 2 seconds")
        raw_cadence_audit = {
            "raw_compact_cadence_inventory_passed": True,
            "n_raw_compact_cadence_inventories_verified": int(len(target_audit)),
            "raw_compact_time_delta_max_s": float(delta.max()),
        }
    orbitid_summary = _orbitid_summary_from_rows(
        peaks,
        orbitid_policy=orbitid_policy,
    )
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "sector": int(sector),
        "contract_version": ADP_ONLY_CONTRACT_VERSION,
        "bls_search_contract_version": search_contract,
        "bls_config_sha256": cfg_sha256,
        "external_quality_policy_contract": EXTERNAL_QUALITY_POLICY_CONTRACT,
        "compact_lc": str(compact_lc),
        "compact_lc_sha256": input_sha256["compact_lc"],
        "compact_lc_role": (
            "cadence_time_inventory_reference_only"
            if raw_mode
            else "photometry_and_cadence_input"
        ),
        "input_mode": "immutable_raw_v1_detector_consistent" if raw_mode else "compact",
        "raw_v4_input_contract_version": (
            RAW_V4_INPUT_CONTRACT_VERSION if raw_mode else None
        ),
        "raw_source_h5": str(raw_source_h5) if raw_mode else None,
        "raw_source_h5_sha256": input_sha256.get("raw_source_h5"),
        "raw_source_summary": str(raw_source_summary) if raw_mode else None,
        "raw_source_summary_sha256": input_sha256.get("raw_source_summary"),
        "raw_export_complete": str(raw_export_complete) if raw_mode else None,
        "raw_export_complete_sha256": input_sha256.get("raw_export_complete"),
        "raw_transfer_validation": (
            str(raw_transfer_validation) if raw_mode else None
        ),
        "raw_transfer_validation_sha256": input_sha256.get(
            "raw_transfer_validation"
        ),
        "execution_allowlist": (
            str(execution_allowlist) if execution_allowlist is not None else None
        ),
        "execution_allowlist_sha256": input_sha256.get("execution_allowlist"),
        "detrend_contract_version": (
            FULL_POOL_NATIVE_DETREND_CONTRACT_VERSION if raw_mode else None
        ),
        "detrend_config_sha256": (
            FULL_POOL_NATIVE_DETREND_CONFIG_SHA256 if raw_mode else None
        ),
        "detrend_time_contract_version": (
            FULL_POOL_NATIVE_DETREND_TIME_CONTRACT_VERSION if raw_mode else None
        ),
        "cadence_reference": str(cadence_reference),
        "cadence_reference_sha256": input_sha256["cadence_reference"],
        "cadence_reference_manifest": str(cadence_reference_manifest),
        "cadence_reference_manifest_sha256": input_sha256["cadence_reference_manifest"],
        "cadence_reference_contract_version": reference_provenance["contract_version"],
        "cadence_reference_cadence_authority": reference_provenance[
            "cadence_authority"
        ],
        "cadence_reference_quality_authority": reference_provenance[
            "quality_authority"
        ],
        "cadence_reference_source_hashes_sha256": reference_provenance[
            "source_file_sha256_declaration_sha256"
        ],
        "authority_exclusion_policy_contract": reference_provenance[
            "authority_exclusion_policy_contract"
        ],
        "authority_exclusion_external_bit": int(
            reference_provenance["authority_exclusion_external_bit"]
        ),
        "authority_exclusions_sha256": reference_provenance[
            "authority_exclusions_sha256"
        ],
        "n_authority_exclusions": int(reference_provenance["n_authority_exclusions"]),
        **target_selection,
        "orbitid_policy": orbitid_policy,
        "orbitid_reconciliation_contract_version": (
            ORBITID_RECONCILIATION_CONTRACT_VERSION
        ),
        **orbitid_summary,
        "out_dir": str(out_dir),
        "apertures": list(ADP_ONLY_APERTURES),
        "n_targets": int(len(tics)),
        "n_targets_total": int(n_targets_total),
        "n_rows": int(len(peaks)),
        "n_unique_tics": int(peaks["tic"].nunique()) if "tic" in peaks else 0,
        "n_valid_peak_rows": int(valid.sum()),
        "n_periods": int(n_periods),
        "n_peaks": int(n_peaks),
        "workers": int(workers),
        "shard_index": int(shard_index),
        "n_shards": int(n_shards),
        "source_product_tag": str(source_product_tag),
        "peak_table_sha256": output_sha256,
        "config": cfg_payload,
        "status_counts": {
            str(k): int(v) for k, v in status.value_counts().sort_index().items()
        },
        "aperture_counts": {
            str(k): int(v)
            for k, v in peaks.get("aperture", pd.Series(dtype=str))
            .fillna("")
            .astype(str)
            .value_counts()
            .sort_index()
            .items()
        },
        "quality_counts_over_unique_targets": (quality_counts_over_unique_targets),
        **target_metadata_summary,
        **raw_cadence_audit,
        "outputs": {
            "peak_table": str(output_path),
            "summary": str(summary_path),
        },
    }
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True, default=json_default) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True, default=json_default))
    return summary


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sector", type=int, default=56)
    parser.add_argument("--compact-lc", type=Path, default=DEFAULT_COMPACT_LC)
    parser.add_argument("--cadence-reference", type=Path, required=True)
    parser.add_argument("--cadence-reference-manifest", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--n-periods", type=int, default=50_000)
    parser.add_argument("--n-peaks", type=int, default=10)
    parser.add_argument("--max-targets", type=int, default=None)
    parser.add_argument("--progress-every", type=int, default=100)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--n-shards", type=int, default=1)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--target-allowlist",
        type=Path,
        default=None,
        help=(
            "Optional CSV with a tic column or one-TIC-per-line text file. "
            "Every listed TIC must exist in the compact export."
        ),
    )
    parser.add_argument(
        "--orbitid-policy",
        choices=ORBITID_POLICIES,
        default="strict",
        help=(
            "strict rejects compact/reference orbit mismatches; "
            "reference_by_cadence replaces orbit IDs only where the "
            "authoritative detector-cadence reference matches."
        ),
    )
    parser.add_argument(
        "--expected-compact-lc-sha256",
        default=None,
        help="Optional frozen SHA-256 that the compact HDF5 must match.",
    )
    parser.add_argument(
        "--expected-target-allowlist-sha256",
        default=None,
        help="Optional frozen SHA-256 that the target allowlist must match.",
    )
    parser.add_argument(
        "--source-product-tag",
        default="",
        help="Set to A2v1 when the compact export came from the production A2v1 FITS tree.",
    )
    parser.add_argument("--raw-source-h5", type=Path, default=None)
    parser.add_argument("--raw-source-summary", type=Path, default=None)
    parser.add_argument("--raw-export-complete", type=Path, default=None)
    parser.add_argument("--raw-transfer-validation", type=Path, default=None)
    parser.add_argument("--frozen-pool", type=Path, default=None)
    parser.add_argument("--frozen-pool-summary", type=Path, default=None)
    parser.add_argument("--execution-allowlist", type=Path, default=None)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    build_peak_table(
        sector=args.sector,
        compact_lc=args.compact_lc,
        cadence_reference=args.cadence_reference,
        cadence_reference_manifest=args.cadence_reference_manifest,
        out_dir=args.out_dir,
        workers=args.workers,
        n_periods=args.n_periods,
        n_peaks=args.n_peaks,
        max_targets=args.max_targets,
        progress_every=args.progress_every,
        shard_index=args.shard_index,
        n_shards=args.n_shards,
        resume=args.resume,
        source_product_tag=args.source_product_tag,
        target_allowlist=args.target_allowlist,
        orbitid_policy=args.orbitid_policy,
        expected_compact_lc_sha256=args.expected_compact_lc_sha256,
        expected_target_allowlist_sha256=(
            args.expected_target_allowlist_sha256
        ),
        raw_source_h5=args.raw_source_h5,
        raw_source_summary=args.raw_source_summary,
        raw_export_complete=args.raw_export_complete,
        raw_transfer_validation=args.raw_transfer_validation,
        frozen_pool=args.frozen_pool,
        frozen_pool_summary=args.frozen_pool_summary,
        execution_allowlist=args.execution_allowlist,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
