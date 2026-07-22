"""Load and apply the authoritative A2v1 detector-cadence quality overlay.

The compact A2v1 products carry an internal quality bitmask.  Tier-1 evidence
adds a detector-specific external mask built from SPOC and QLP authorities.
This module validates that evidence pair and exposes the only accepted merge
policy for model-facing native inputs: a cadence is bad when either mask is
nonzero.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from twirl.lightcurves.a2v1_cadence_reference import (
    CADENCE_REFERENCE_COLUMNS,
    QLP_QUALITY_EXTERNAL_BIT,
    file_sha256,
)


EXTERNAL_QUALITY_POLICY_CONTRACT = (
    "a2v1_native_internal_or_authoritative_external_quality_v1"
)
EFFECTIVE_QUALITY_POLICY = "(internal_quality != 0) | (external_quality != 0)"
CADENCE_REFERENCE_CONTRACT_TEMPLATE = "s{sector}_a2v1_cadence_reference_v1"
EXPECTED_CADENCE_AUTHORITY = "qlp_cam_quat"
EXPECTED_QUALITY_AUTHORITY = "spoc_and_qlp_quality_flags"
EXPECTED_QUALITY_COMPOSITION = {
    "external_quality": "spoc_quality | (qlp_quality << 30)",
    "qlp_quality_raw_values": [0, 1],
    "qlp_quality_external_bit": QLP_QUALITY_EXTERNAL_BIT,
}


def _read_table(path: Path) -> pd.DataFrame:
    suffix = Path(path).suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path, low_memory=False)
    if suffix == ".parquet":
        return pd.read_parquet(path)
    raise ValueError("cadence reference must end in .csv or .parquet")


def _integer_array(values: Any, *, context: str) -> np.ndarray:
    numeric = pd.to_numeric(pd.Series(np.asarray(values)), errors="raise")
    array = numeric.to_numpy(dtype=np.float64)
    if array.ndim != 1:
        raise ValueError(f"{context} must be one-dimensional")
    if not np.all(np.isfinite(array)) or not np.all(array == np.floor(array)):
        raise ValueError(f"{context} must contain finite integers")
    return array.astype(np.int64)


def _canonical_json_sha256(payload: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _validate_source_inventory(
    manifest: Mapping[str, Any],
) -> tuple[
    str,
    set[tuple[int, int]],
    set[tuple[int, int, int]],
    set[tuple[int, int]],
    Mapping[str, int],
]:
    declared = manifest.get("source_file_sha256")
    if not isinstance(declared, Mapping) or not declared:
        raise ValueError("cadence-reference source_file_sha256 is missing or empty")
    normalized: dict[str, str] = {}
    for raw_path, raw_digest in declared.items():
        path = str(raw_path).strip()
        digest = str(raw_digest).strip().lower()
        if not path or not Path(path).is_absolute():
            raise ValueError("cadence-reference source paths must be absolute")
        if len(digest) != 64 or any(
            character not in "0123456789abcdef" for character in digest
        ):
            raise ValueError(f"invalid cadence-reference source hash for {path!r}")
        if path in normalized:
            raise ValueError(f"duplicate cadence-reference source path: {path}")
        normalized[path] = digest

    sources = manifest.get("sources")
    if not isinstance(sources, list) or not sources:
        raise ValueError("cadence-reference sources inventory is missing or empty")
    observed: set[str] = set()
    spoc_keys: set[tuple[int, int]] = set()
    qflag_keys: set[tuple[int, int, int]] = set()
    quat_keys: set[tuple[int, int]] = set()
    role_counts: dict[str, int] = {}
    for index, source in enumerate(sources):
        if not isinstance(source, Mapping):
            raise ValueError(f"cadence-reference source {index} is not an object")
        path = str(source.get("path", "")).strip()
        digest = str(source.get("sha256", "")).strip().lower()
        if not path or path in observed:
            raise ValueError("cadence-reference sources have an empty/duplicate path")
        observed.add(path)
        if normalized.get(path) != digest:
            raise ValueError(
                "cadence-reference sources disagree with source_file_sha256 for "
                f"{path!r}"
            )
        role = str(source.get("role", ""))
        role_counts[role] = role_counts.get(role, 0) + 1
        if role == "spoc_flag_file":
            key = (int(source.get("camera", -1)), int(source.get("ccd", -1)))
            if key in spoc_keys:
                raise ValueError("duplicate SPOC detector in source inventory")
            spoc_keys.add(key)
        elif role == "qlp_detector_qflag":
            key = (
                int(source.get("orbitid", -1)),
                int(source.get("camera", -1)),
                int(source.get("ccd", -1)),
            )
            if key in qflag_keys:
                raise ValueError("duplicate QLP qflag mapping in source inventory")
            qflag_keys.add(key)
        elif role == "qlp_cam_quat":
            key = (
                int(source.get("orbitid", -1)),
                int(source.get("camera", -1)),
            )
            if key in quat_keys:
                raise ValueError("duplicate QLP quaternion mapping in source inventory")
            quat_keys.add(key)
    if observed != set(normalized):
        raise ValueError(
            "cadence-reference sources do not exactly cover source_file_sha256"
        )
    return (
        _canonical_json_sha256(normalized),
        spoc_keys,
        qflag_keys,
        quat_keys,
        role_counts,
    )


@dataclass(frozen=True)
class QualityOverlayResult:
    """Effective binary quality plus exact external values and audit counts."""

    quality: np.ndarray
    external_quality: np.ndarray
    counts: Mapping[str, int]


@dataclass(frozen=True)
class ExternalQualityReference:
    """Validated immutable detector-cadence lookup and provenance envelope."""

    sector: int
    table_path: Path
    manifest_path: Path
    table_sha256: str
    manifest_sha256: str
    source_declaration_sha256: str
    contract_version: str
    cadence_authority: str
    quality_authority: str
    by_detector: Mapping[tuple[int, int, int], tuple[np.ndarray, np.ndarray, np.ndarray]]

    @property
    def provenance(self) -> dict[str, Any]:
        return {
            "policy_contract": EXTERNAL_QUALITY_POLICY_CONTRACT,
            "effective_quality_policy": EFFECTIVE_QUALITY_POLICY,
            "sector": int(self.sector),
            "cadence_reference_contract_version": self.contract_version,
            "cadence_reference_cadence_authority": self.cadence_authority,
            "cadence_reference_quality_authority": self.quality_authority,
            "cadence_reference_table": str(self.table_path),
            "cadence_reference_manifest": str(self.manifest_path),
            "cadence_reference_table_sha256": self.table_sha256,
            "cadence_reference_manifest_sha256": self.manifest_sha256,
            "cadence_reference_source_declaration_sha256": (
                self.source_declaration_sha256
            ),
        }

    def assert_unchanged(self) -> None:
        """Fail if either published evidence file changed after loading."""

        if file_sha256(self.table_path) != self.table_sha256:
            raise RuntimeError("cadence-reference table changed while in use")
        if file_sha256(self.manifest_path) != self.manifest_sha256:
            raise RuntimeError("cadence-reference manifest changed while in use")

    def apply(
        self,
        *,
        sector: int,
        camera: int,
        ccd: int,
        cadenceno: Any,
        orbitid: Any,
        internal_quality: Any,
        context: str = "light curve",
    ) -> QualityOverlayResult:
        """Apply the exact detector/cadence/orbit overlay to one light curve."""

        sector = int(sector)
        camera = int(camera)
        ccd = int(ccd)
        if sector != self.sector:
            raise ValueError(
                f"{context}: sector {sector} does not match reference {self.sector}"
            )
        key = (sector, camera, ccd)
        reference = self.by_detector.get(key)
        if reference is None:
            raise ValueError(f"{context}: no external-quality mapping for {key}")

        cadences = _integer_array(cadenceno, context=f"{context} cadenceno")
        orbits = _integer_array(orbitid, context=f"{context} orbitid")
        internal = _integer_array(
            internal_quality, context=f"{context} internal quality"
        )
        lengths = {len(cadences), len(orbits), len(internal)}
        if len(lengths) != 1:
            raise ValueError(f"{context}: cadence/orbit/quality lengths disagree")
        if not len(cadences):
            raise ValueError(f"{context}: cadence arrays are empty")
        if np.any(internal < 0):
            raise ValueError(f"{context}: internal quality must be nonnegative")
        if len(np.unique(cadences)) != len(cadences):
            raise ValueError(f"{context}: cadenceno values are not unique")

        reference_cadence, reference_orbit, reference_quality = reference
        positions = np.searchsorted(reference_cadence, cadences)
        in_bounds = positions < len(reference_cadence)
        matched = np.zeros(len(cadences), dtype=bool)
        matched[in_bounds] = (
            reference_cadence[positions[in_bounds]] == cadences[in_bounds]
        )
        if not matched.all():
            missing = cadences[~matched][:20].astype(int).tolist()
            raise ValueError(
                f"{context}: external-quality coverage is missing cadences for "
                f"{key}: {missing}"
            )
        mapped_orbits = reference_orbit[positions]
        mismatch = mapped_orbits != orbits
        if mismatch.any():
            examples = [
                {
                    "cadenceno": int(cadences[index]),
                    "input_orbitid": int(orbits[index]),
                    "reference_orbitid": int(mapped_orbits[index]),
                }
                for index in np.flatnonzero(mismatch)[:20]
            ]
            raise ValueError(
                f"{context}: external-quality orbit mapping mismatch: {examples}"
            )

        external = reference_quality[positions]
        internal_bad = internal != 0
        external_bad = external != 0
        effective_bad = internal_bad | external_bad
        counts = {
            "n_cad_total": int(len(cadences)),
            "n_cad_internal_bad": int(np.count_nonzero(internal_bad)),
            "n_cad_external_bad": int(np.count_nonzero(external_bad)),
            "n_cad_external_only_bad": int(
                np.count_nonzero(external_bad & ~internal_bad)
            ),
            "n_cad_effective_bad": int(np.count_nonzero(effective_bad)),
        }
        return QualityOverlayResult(
            quality=effective_bad.astype(np.int32),
            external_quality=external.astype(np.int64, copy=False),
            counts=counts,
        )


def load_external_quality_reference(
    *,
    table_path: Path,
    manifest_path: Path,
    sector: int,
    expected_orbits: Sequence[int] | None = None,
    expected_detectors: Sequence[tuple[int, int]] | None = None,
) -> ExternalQualityReference:
    """Load and validate a final A2v1 cadence-reference evidence pair."""

    table_path = Path(table_path).resolve()
    manifest_path = Path(manifest_path).resolve()
    table_sha256 = file_sha256(table_path)
    manifest_sha256 = file_sha256(manifest_path)
    with manifest_path.open("r", encoding="utf-8") as handle:
        manifest = json.load(handle)
    if not isinstance(manifest, dict):
        raise ValueError("cadence-reference manifest must contain a JSON object")

    required = {
        "contract_version",
        "sector",
        "cadence_authority",
        "quality_authority",
        "quality_composition",
        "table_sha256",
        "table_columns",
        "n_rows",
        "detectors",
        "orbits",
        "n_rows_by_detector",
        "n_nonzero_spoc_quality",
        "n_nonzero_qlp_quality",
        "n_nonzero_external_quality",
        "n_spoc_authority_files_verified",
        "n_qlp_qflag_files_verified",
        "source_file_sha256",
        "sources",
    }
    missing = sorted(required - set(manifest))
    if missing:
        raise ValueError(f"cadence-reference manifest is missing fields: {missing}")
    sector = int(sector)
    contract = CADENCE_REFERENCE_CONTRACT_TEMPLATE.format(sector=sector)
    if str(manifest["contract_version"]) != contract:
        raise ValueError("cadence-reference contract mismatch")
    if int(manifest["sector"]) != sector:
        raise ValueError("cadence-reference manifest sector mismatch")
    if str(manifest["cadence_authority"]) != EXPECTED_CADENCE_AUTHORITY:
        raise ValueError("cadence-reference cadence authority mismatch")
    if str(manifest["quality_authority"]) != EXPECTED_QUALITY_AUTHORITY:
        raise ValueError("cadence-reference quality authority mismatch")
    if manifest["quality_composition"] != EXPECTED_QUALITY_COMPOSITION:
        raise ValueError("cadence-reference quality composition mismatch")
    if str(manifest["table_sha256"]) != table_sha256:
        raise ValueError("cadence-reference table hash mismatch")
    if tuple(manifest["table_columns"]) != CADENCE_REFERENCE_COLUMNS:
        raise ValueError("cadence-reference manifest table columns mismatch")

    frame = _read_table(table_path)
    if tuple(frame.columns) != CADENCE_REFERENCE_COLUMNS:
        raise ValueError(
            "cadence-reference columns must be exactly "
            f"{CADENCE_REFERENCE_COLUMNS}; observed {tuple(frame.columns)}"
        )
    if frame.empty or len(frame) != int(manifest["n_rows"]):
        raise ValueError("cadence-reference row count disagrees with manifest")
    for column in CADENCE_REFERENCE_COLUMNS:
        frame[column] = _integer_array(
            frame[column], context=f"cadence-reference {column}"
        )
    if set(frame["sector"].astype(int)) != {sector}:
        raise ValueError("cadence-reference table sector mismatch")
    if not frame["camera"].between(1, 4).all() or not frame["ccd"].between(1, 4).all():
        raise ValueError("cadence-reference camera/ccd values must be in 1..4")
    if (frame[["orbitid", "cadenceno", "spoc_quality", "external_quality"]] < 0).any().any():
        raise ValueError("cadence-reference integer values must be nonnegative")
    if not set(frame["qlp_quality"].astype(int)).issubset({0, 1}):
        raise ValueError("cadence-reference qlp_quality must contain only 0/1")
    derived = np.bitwise_or(
        frame["spoc_quality"].to_numpy(dtype=np.int64),
        np.left_shift(
            frame["qlp_quality"].to_numpy(dtype=np.int64),
            QLP_QUALITY_EXTERNAL_BIT,
        ),
    )
    if not np.array_equal(
        derived, frame["external_quality"].to_numpy(dtype=np.int64)
    ):
        raise ValueError("cadence-reference external_quality derivation is invalid")
    key_columns = ["sector", "camera", "ccd", "cadenceno"]
    if frame.duplicated(key_columns).any():
        raise ValueError("cadence-reference has duplicate detector-cadence mappings")

    observed_detectors = sorted(
        {
            (int(camera), int(ccd))
            for camera, ccd in zip(frame["camera"], frame["ccd"], strict=True)
        }
    )
    detector_names = [f"cam{camera}_ccd{ccd}" for camera, ccd in observed_detectors]
    if sorted(str(value) for value in manifest["detectors"]) != detector_names:
        raise ValueError("cadence-reference detector inventory mismatch")
    observed_orbits = sorted(set(frame["orbitid"].astype(int)))
    if sorted(int(value) for value in manifest["orbits"]) != observed_orbits:
        raise ValueError("cadence-reference orbit inventory mismatch")
    if expected_detectors is not None and set(observed_detectors) != {
        (int(camera), int(ccd)) for camera, ccd in expected_detectors
    }:
        raise ValueError("cadence-reference does not cover the expected detectors")
    if expected_orbits is not None and set(observed_orbits) != {
        int(value) for value in expected_orbits
    }:
        raise ValueError("cadence-reference does not cover the expected orbits")

    observed_rows_by_detector = {
        f"cam{int(camera)}_ccd{int(ccd)}": int(len(group))
        for (camera, ccd), group in frame.groupby(["camera", "ccd"], sort=True)
    }
    declared_rows = {
        str(key): int(value)
        for key, value in dict(manifest["n_rows_by_detector"]).items()
    }
    if declared_rows != observed_rows_by_detector:
        raise ValueError("cadence-reference detector row counts mismatch")
    for field, column in (
        ("n_nonzero_spoc_quality", "spoc_quality"),
        ("n_nonzero_qlp_quality", "qlp_quality"),
        ("n_nonzero_external_quality", "external_quality"),
    ):
        if int(manifest[field]) != int(np.count_nonzero(frame[column])):
            raise ValueError(f"cadence-reference {field} disagrees with the table")

    (
        source_digest,
        spoc_keys,
        qflag_keys,
        quat_keys,
        role_counts,
    ) = _validate_source_inventory(manifest)
    if int(manifest["n_spoc_authority_files_verified"]) != len(spoc_keys):
        raise ValueError("cadence-reference SPOC source count mismatch")
    if int(manifest["n_qlp_qflag_files_verified"]) != len(qflag_keys):
        raise ValueError("cadence-reference QLP qflag source count mismatch")
    if spoc_keys != set(observed_detectors):
        raise ValueError("cadence-reference SPOC source inventory is incomplete")
    expected_qflag_keys = {
        (orbit, camera, ccd)
        for orbit in observed_orbits
        for camera, ccd in observed_detectors
    }
    if qflag_keys != expected_qflag_keys:
        raise ValueError("cadence-reference QLP qflag inventory is incomplete")
    expected_quat_keys = {
        (orbit, camera)
        for orbit in observed_orbits
        for camera in sorted({camera for camera, _ in observed_detectors})
    }
    if quat_keys != expected_quat_keys:
        raise ValueError("cadence-reference QLP quaternion inventory is incomplete")
    for role in ("spoc_quality_table", "spoc_quality_provenance"):
        if role_counts.get(role, 0) != 1:
            raise ValueError(f"cadence-reference must declare exactly one {role}")
    allowed_roles = {
        "qlp_cam_quat",
        "qlp_detector_qflag",
        "spoc_quality_table",
        "spoc_quality_provenance",
        "spoc_flag_file",
    }
    unknown_roles = sorted(set(role_counts) - allowed_roles)
    if unknown_roles:
        raise ValueError(f"cadence-reference has unknown source roles: {unknown_roles}")

    by_detector: dict[
        tuple[int, int, int], tuple[np.ndarray, np.ndarray, np.ndarray]
    ] = {}
    for (frame_sector, camera, ccd), group in frame.groupby(
        ["sector", "camera", "ccd"], sort=True
    ):
        ordered = group.sort_values("cadenceno", kind="stable")
        by_detector[(int(frame_sector), int(camera), int(ccd))] = (
            ordered["cadenceno"].to_numpy(dtype=np.int64),
            ordered["orbitid"].to_numpy(dtype=np.int64),
            ordered["external_quality"].to_numpy(dtype=np.int64),
        )

    # Rehash after parsing so a concurrent replacement cannot be silently used.
    if file_sha256(table_path) != table_sha256:
        raise RuntimeError("cadence-reference table changed while loading")
    if file_sha256(manifest_path) != manifest_sha256:
        raise RuntimeError("cadence-reference manifest changed while loading")
    return ExternalQualityReference(
        sector=sector,
        table_path=table_path,
        manifest_path=manifest_path,
        table_sha256=table_sha256,
        manifest_sha256=manifest_sha256,
        source_declaration_sha256=source_digest,
        contract_version=contract,
        cadence_authority=EXPECTED_CADENCE_AUTHORITY,
        quality_authority=EXPECTED_QUALITY_AUTHORITY,
        by_detector=by_detector,
    )


__all__ = [
    "CADENCE_REFERENCE_CONTRACT_TEMPLATE",
    "EFFECTIVE_QUALITY_POLICY",
    "EXPECTED_CADENCE_AUTHORITY",
    "EXPECTED_QUALITY_AUTHORITY",
    "EXTERNAL_QUALITY_POLICY_CONTRACT",
    "ExternalQualityReference",
    "QualityOverlayResult",
    "load_external_quality_reference",
]
