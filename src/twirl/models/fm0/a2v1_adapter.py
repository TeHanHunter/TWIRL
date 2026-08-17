"""Checksum-bound A2v1 HDF5 adapter for scientific FM0.1 input releases.

This adapter is deliberately separate from the strict-NPZ fixture path.  It
reads immutable TGLC/A2v1 HDF5 source files, verifies every declared digest,
joins the validated external cadence-quality authority, and emits only the
raw arrays allowed by :mod:`twirl.models.fm0.input_release`.

It never reads BLS products, labels, catalog features, or candidate tables.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import csv
import hashlib
from io import StringIO
import json
from pathlib import Path
import re
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from twirl.lightcurves.external_quality import (
    AUTHORITY_EXCLUSION_EXTERNAL_BIT,
    ExternalQualityReference,
    load_external_quality_reference,
)
from twirl.lightcurves.tglc_h5_reader import read_tglc_h5

from .registry import (
    AliasRegistry,
    FM0ContractError,
    assert_no_search_columns,
    build_observation_registry,
    publish_immutable,
    sha256_file,
)


A2V1_HDF5_ADAPTER_NAME = "a2v1_hdf5_quality_aware_v1"
A2V1_HDF5_MANIFEST_FIELDS = (
    "observation_key",
    "product_instance_id",
    "source_sha256",
    "hdf5_paths_json",
    "hdf5_sha256_json",
    "quality_table_path",
    "quality_table_sha256",
    "quality_manifest_path",
    "quality_manifest_sha256",
)
A2V1_HDF5_MANIFEST_COLUMNS = frozenset(A2V1_HDF5_MANIFEST_FIELDS)
A2V1_HDF5_SOURCE_INVENTORY_FIELDS = (
    "gaia_dr3_source_id",
    "tic_id",
    "sector",
    "a2v1_product_version",
    "product_state",
    "diagnostic_admission_receipt_path",
    "diagnostic_admission_receipt_sha256",
    "hdf5_paths_json",
    "hdf5_sha256_json",
    "quality_table_path",
    "quality_table_sha256",
    "quality_manifest_path",
    "quality_manifest_sha256",
)
A2V1_HDF5_SOURCE_INVENTORY_COLUMNS = frozenset(
    A2V1_HDF5_SOURCE_INVENTORY_FIELDS
)


@dataclass(frozen=True)
class A2V1RawObservation:
    """A source-bound raw observation ready for the six-view builder."""

    raw_arrays: Mapping[str, np.ndarray]
    source_sha256: str
    audit: Mapping[str, Any]


@dataclass
class A2V1AdapterCache:
    """Per-build cache for already validated immutable quality authorities.

    The full source manifest may contain many targets from one sector.  Parsing
    and hashing the same validated cadence table for every target would make a
    correct build needlessly expensive.  The release builder rechecks every
    cached authority after the final source is read.
    """

    references: dict[tuple[int, Path, Path], ExternalQualityReference] = field(
        default_factory=dict
    )
    components: dict[tuple[int, Path, int, int], dict[int, tuple[int, int]]] = field(
        default_factory=dict
    )

    def assert_unchanged(self) -> None:
        for reference in self.references.values():
            reference.assert_unchanged()


def _sha256_text(value: str, *, name: str) -> str:
    normalized = str(value).strip().lower()
    if not re.fullmatch(r"[0-9a-f]{64}", normalized):
        raise FM0ContractError(f"{name} must be a lowercase SHA-256 digest")
    return normalized


def _manifest_path(value: Any, *, base: Path, name: str) -> Path:
    text = str(value).strip()
    if not text:
        raise FM0ContractError(f"{name} is required")
    path = Path(text).expanduser()
    if not path.is_absolute():
        path = base / path
    path = path.resolve()
    if not path.is_file():
        raise FM0ContractError(f"{name} does not exist: {path}")
    return path


def _json_string_list(value: Any, *, name: str) -> list[str]:
    try:
        parsed = json.loads(str(value))
    except (TypeError, json.JSONDecodeError) as exc:
        raise FM0ContractError(f"{name} must be a JSON string list") from exc
    if not isinstance(parsed, list) or not parsed or not all(
        isinstance(item, str) and item.strip() for item in parsed
    ):
        raise FM0ContractError(f"{name} must be a non-empty JSON string list")
    return [item.strip() for item in parsed]


def a2v1_source_sha256(
    *,
    sector: int,
    tic_id: str,
    hdf5_sha256s: Sequence[str],
    quality_table_sha256: str,
    quality_manifest_sha256: str,
) -> str:
    """Return the portable content identity for one source-sector product.

    Paths are intentionally excluded: moving immutable source products must
    not change their identity.  The sector/TIC binding prevents the same set
    of files from being reassigned to a different registry observation.
    """

    payload = {
        "adapter": A2V1_HDF5_ADAPTER_NAME,
        "sector": int(sector),
        "tic_id": str(tic_id),
        "hdf5_sha256s": sorted(
            _sha256_text(item, name="hdf5_sha256") for item in hdf5_sha256s
        ),
        "quality_table_sha256": _sha256_text(
            quality_table_sha256, name="quality_table_sha256"
        ),
        "quality_manifest_sha256": _sha256_text(
            quality_manifest_sha256, name="quality_manifest_sha256"
        ),
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def bind_a2v1_hdf5_source_inventory(
    rows: Sequence[Mapping[str, Any]], alias_registry: AliasRegistry
) -> tuple[tuple[dict[str, Any], ...], tuple[dict[str, Any], ...]]:
    """Derive registry observations and their bound HDF5 source manifest.

    The initial inventory identifies source files by Gaia/TIC/sector but must
    not guess an ``observation_key`` or ``product_instance_id``.  Those are
    created only after the authoritative Gaia--TIC leakage registry exists.
    This function makes that two-stage binding deterministic and auditable.
    """

    observation_inputs: list[dict[str, Any]] = []
    source_rows: dict[tuple[str, str, int], Mapping[str, Any]] = {}
    for raw in rows:
        if frozenset(raw) != A2V1_HDF5_SOURCE_INVENTORY_COLUMNS:
            missing = sorted(A2V1_HDF5_SOURCE_INVENTORY_COLUMNS - frozenset(raw))
            extra = sorted(frozenset(raw) - A2V1_HDF5_SOURCE_INVENTORY_COLUMNS)
            raise FM0ContractError(
                f"A2v1 HDF5 source-inventory columns mismatch; missing={missing}, extra={extra}"
            )
        assert_no_search_columns(raw.keys(), context="A2v1 HDF5 source inventory")
        gaia = str(raw["gaia_dr3_source_id"]).strip()
        tic = str(raw["tic_id"]).strip()
        try:
            sector = int(raw["sector"])
        except (TypeError, ValueError) as exc:
            raise FM0ContractError("A2v1 source-inventory sector must be an integer") from exc
        if not gaia or not tic:
            raise FM0ContractError("A2v1 source inventory requires Gaia and TIC identifiers")
        hdf5_hashes = _json_string_list(
            raw["hdf5_sha256_json"], name="hdf5_sha256_json"
        )
        source_hash = a2v1_source_sha256(
            sector=sector,
            tic_id=tic,
            hdf5_sha256s=hdf5_hashes,
            quality_table_sha256=str(raw["quality_table_sha256"]),
            quality_manifest_sha256=str(raw["quality_manifest_sha256"]),
        )
        identity = (gaia, tic, sector)
        if identity in source_rows:
            raise FM0ContractError(
                f"duplicate A2v1 HDF5 source-inventory identity: {identity}"
            )
        source_rows[identity] = raw
        observation_inputs.append(
            {
                "gaia_dr3_source_id": gaia,
                "tic_id": tic,
                "sector": sector,
                "a2v1_product_version": str(raw["a2v1_product_version"]),
                "product_state": str(raw["product_state"]),
                "source_sha256": source_hash,
                "diagnostic_admission_receipt_path": str(
                    raw["diagnostic_admission_receipt_path"]
                ),
                "diagnostic_admission_receipt_sha256": str(
                    raw["diagnostic_admission_receipt_sha256"]
                ),
            }
        )
    observations = build_observation_registry(observation_inputs, alias_registry)
    bound_rows: list[dict[str, Any]] = []
    for observation in observations:
        identity = (
            str(observation["gaia_dr3_source_id"]),
            str(observation["tic_id"]),
            int(observation["sector"]),
        )
        source = source_rows[identity]
        bound_rows.append(
            {
                "observation_key": observation["observation_key"],
                "product_instance_id": observation["product_instance_id"],
                "source_sha256": observation["source_sha256"],
                "hdf5_paths_json": source["hdf5_paths_json"],
                "hdf5_sha256_json": source["hdf5_sha256_json"],
                "quality_table_path": source["quality_table_path"],
                "quality_table_sha256": source["quality_table_sha256"],
                "quality_manifest_path": source["quality_manifest_path"],
                "quality_manifest_sha256": source["quality_manifest_sha256"],
            }
        )
    bound_rows.sort(key=lambda row: row["observation_key"])
    return observations, tuple(bound_rows)


def write_a2v1_hdf5_bound_manifest(
    path: str | Path, rows: Sequence[Mapping[str, Any]]
) -> str:
    """Publish the registry-bound source manifest as an immutable CSV."""

    if not rows:
        raise FM0ContractError("cannot publish an empty A2v1 HDF5 manifest")
    stream = StringIO(newline="")
    writer = csv.DictWriter(
        stream, fieldnames=list(A2V1_HDF5_MANIFEST_FIELDS), lineterminator="\n"
    )
    writer.writeheader()
    for row in rows:
        if frozenset(row) != A2V1_HDF5_MANIFEST_COLUMNS:
            raise FM0ContractError("bound A2v1 HDF5 manifest row has drifted columns")
        writer.writerow({field: row[field] for field in A2V1_HDF5_MANIFEST_FIELDS})
    payload = stream.getvalue().encode("utf-8")
    publish_immutable(path, payload)
    return hashlib.sha256(payload).hexdigest()


def _quality_components(
    reference: ExternalQualityReference,
    *,
    sector: int,
    camera: int,
    ccd: int,
) -> dict[int, tuple[int, int]]:
    """Read the separately verified SPOC/QLP components for one detector."""

    try:
        frame = pd.read_csv(reference.table_path)
    except Exception as exc:  # pandas gives varied parser errors.
        raise FM0ContractError("cannot reread verified external-quality table") from exc
    selected = frame.loc[
        (frame["sector"] == int(sector))
        & (frame["camera"] == int(camera))
        & (frame["ccd"] == int(ccd)),
        ["cadenceno", "spoc_quality", "qlp_quality"],
    ]
    if selected.empty:
        raise FM0ContractError(
            f"verified quality table has no rows for S{sector} cam{camera}/ccd{ccd}"
        )
    return {
        int(row.cadenceno): (int(row.spoc_quality), int(row.qlp_quality))
        for row in selected.itertuples(index=False)
    }


def _concatenate(parts: Sequence[np.ndarray], *, dtype: Any) -> np.ndarray:
    if not parts:
        raise FM0ContractError("A2v1 source contains no HDF5 cadence arrays")
    return np.concatenate([np.asarray(part, dtype=dtype) for part in parts])


def load_a2v1_hdf5_observation(
    row: Mapping[str, Any],
    *,
    manifest_dir: str | Path,
    sector: int,
    tic_id: str,
    cache: A2V1AdapterCache | None = None,
) -> A2V1RawObservation:
    """Load one checksum-bound raw A2v1 source-sector observation.

    The caller must compare the returned ``source_sha256`` with both the
    registry and manifest row before it can publish a model release.
    """

    if frozenset(row) != A2V1_HDF5_MANIFEST_COLUMNS:
        missing = sorted(A2V1_HDF5_MANIFEST_COLUMNS - frozenset(row))
        extra = sorted(frozenset(row) - A2V1_HDF5_MANIFEST_COLUMNS)
        raise FM0ContractError(
            f"A2v1 HDF5 manifest columns mismatch; missing={missing}, extra={extra}"
        )
    assert_no_search_columns(row.keys(), context="A2v1 HDF5 manifest")
    base = Path(manifest_dir).resolve()
    path_texts = _json_string_list(row["hdf5_paths_json"], name="hdf5_paths_json")
    expected_hashes = _json_string_list(
        row["hdf5_sha256_json"], name="hdf5_sha256_json"
    )
    if len(path_texts) != len(expected_hashes):
        raise FM0ContractError("HDF5 path/hash lists have different lengths")
    hdf5_paths = [
        _manifest_path(path, base=base, name="hdf5_paths_json") for path in path_texts
    ]
    normalized_hashes = [
        _sha256_text(value, name="hdf5_sha256_json") for value in expected_hashes
    ]
    if len(set(hdf5_paths)) != len(hdf5_paths):
        raise FM0ContractError("duplicate HDF5 paths are not admissible")
    for path, expected in zip(hdf5_paths, normalized_hashes, strict=True):
        observed = sha256_file(path)
        if observed != expected:
            raise FM0ContractError(
                f"A2v1 HDF5 hash mismatch for {path}: expected {expected}, observed {observed}"
            )

    quality_table = _manifest_path(
        row["quality_table_path"], base=base, name="quality_table_path"
    )
    quality_manifest = _manifest_path(
        row["quality_manifest_path"], base=base, name="quality_manifest_path"
    )
    quality_table_hash = _sha256_text(
        row["quality_table_sha256"], name="quality_table_sha256"
    )
    quality_manifest_hash = _sha256_text(
        row["quality_manifest_sha256"], name="quality_manifest_sha256"
    )
    cache_key = (int(sector), quality_table, quality_manifest)
    reference = None if cache is None else cache.references.get(cache_key)
    if reference is None:
        if sha256_file(quality_table) != quality_table_hash:
            raise FM0ContractError("quality_table_sha256 does not match the staged table")
        if sha256_file(quality_manifest) != quality_manifest_hash:
            raise FM0ContractError("quality_manifest_sha256 does not match the staged manifest")
    elif (
        reference.table_sha256 != quality_table_hash
        or reference.manifest_sha256 != quality_manifest_hash
    ):
        raise FM0ContractError("cached external-quality authority conflicts with manifest")
    computed_source_hash = a2v1_source_sha256(
        sector=sector,
        tic_id=tic_id,
        hdf5_sha256s=normalized_hashes,
        quality_table_sha256=quality_table_hash,
        quality_manifest_sha256=quality_manifest_hash,
    )
    declared_source_hash = _sha256_text(row["source_sha256"], name="source_sha256")
    if declared_source_hash != computed_source_hash:
        raise FM0ContractError(
            "A2v1 HDF5 manifest source_sha256 does not bind its source artifacts"
        )
    if reference is None:
        try:
            reference = load_external_quality_reference(
                table_path=quality_table,
                manifest_path=quality_manifest,
                sector=int(sector),
            )
        except Exception as exc:
            raise FM0ContractError("external-quality reference failed its authority checks") from exc
        reference.assert_unchanged()
        if cache is not None:
            cache.references[cache_key] = reference

    arrays: dict[str, list[np.ndarray]] = {
        "time": [],
        "cadence": [],
        "orbit": [],
        "internal_quality": [],
        "spoc_quality": [],
        "qlp_quality": [],
        "authority_excluded": [],
        "raw_flux_1x1": [],
        "raw_flux_error_1x1": [],
        "raw_flux_3x3": [],
        "raw_flux_error_3x3": [],
    }
    seen_orbits: set[int] = set()
    quality_counts: dict[str, int] = {
        "n_cad_total": 0,
        "n_cad_internal_bad": 0,
        "n_cad_external_bad": 0,
        "n_cad_authority_excluded": 0,
    }
    for path in hdf5_paths:
        try:
            lightcurve = read_tglc_h5(path)
        except Exception as exc:
            raise FM0ContractError(f"cannot read A2v1 HDF5 source: {path}") from exc
        if lightcurve.tic != int(tic_id) or lightcurve.sector != int(sector):
            raise FM0ContractError(
                f"HDF5 {path} identifies TIC {lightcurve.tic}/S{lightcurve.sector}, "
                f"expected TIC {tic_id}/S{sector}"
            )
        if lightcurve.orbit in seen_orbits:
            raise FM0ContractError(f"duplicate orbit {lightcurve.orbit} in one source")
        seen_orbits.add(lightcurve.orbit)
        if any(aperture not in lightcurve.apertures for aperture in ("Small", "Primary")):
            raise FM0ContractError("A2v1 HDF5 lacks Small or Primary aperture")
        if lightcurve.apertures["Small"].flux_was_synthesized or lightcurve.apertures[
            "Primary"
        ].flux_was_synthesized:
            raise FM0ContractError("legacy magnitude-derived flux is not admissible for FM0.1")
        try:
            overlay = reference.apply(
                sector=sector,
                camera=lightcurve.cam,
                ccd=lightcurve.ccd,
                cadenceno=lightcurve.cadence,
                orbitid=np.full(lightcurve.cadence.size, lightcurve.orbit, dtype=np.int64),
                internal_quality=lightcurve.quality,
                context=f"FM0.1 {path.name}",
            )
        except Exception as exc:
            raise FM0ContractError(f"external-quality join failed for {path}") from exc
        component_key = (int(sector), quality_table, lightcurve.cam, lightcurve.ccd)
        components = None if cache is None else cache.components.get(component_key)
        if components is None:
            components = _quality_components(
                reference,
                sector=sector,
                camera=lightcurve.cam,
                ccd=lightcurve.ccd,
            )
            if cache is not None:
                cache.components[component_key] = components
        spoc = np.zeros(lightcurve.cadence.size, dtype=np.uint64)
        qlp = np.zeros(lightcurve.cadence.size, dtype=np.uint64)
        authority_excluded = (
            np.asarray(overlay.external_quality, dtype=np.int64)
            == np.int64(1 << AUTHORITY_EXCLUSION_EXTERNAL_BIT)
        )
        for index, cadence in enumerate(np.asarray(lightcurve.cadence, dtype=np.int64)):
            values = components.get(int(cadence))
            if values is None:
                if not authority_excluded[index]:
                    raise FM0ContractError(
                        f"quality components missing non-excluded cadence {cadence}"
                    )
                continue
            spoc[index], qlp[index] = values
        composed_external = spoc | (qlp << np.uint64(30))
        expected_external = np.asarray(overlay.external_quality, dtype=np.int64)
        expected_external = np.where(authority_excluded, 0, expected_external)
        if not np.array_equal(composed_external, expected_external.astype(np.uint64)):
            raise FM0ContractError("external quality components do not reproduce authority")
        # The HDF5 ``BJD`` array is BTJD.  It retains all physical intervals
        # while avoiding numerically ill-conditioned spline fits at absolute
        # BJD ~2.46e6; the common BJD offset is audit provenance, not a model
        # feature.
        arrays["time"].append(lightcurve.time)
        arrays["cadence"].append(lightcurve.cadence)
        arrays["orbit"].append(
            np.full(lightcurve.cadence.size, lightcurve.orbit, dtype=np.int64)
        )
        arrays["internal_quality"].append(lightcurve.quality)
        arrays["spoc_quality"].append(spoc)
        arrays["qlp_quality"].append(qlp)
        arrays["authority_excluded"].append(authority_excluded)
        arrays["raw_flux_1x1"].append(lightcurve.apertures["Small"].raw_flux)
        arrays["raw_flux_error_1x1"].append(
            lightcurve.apertures["Small"].raw_flux_err
        )
        arrays["raw_flux_3x3"].append(lightcurve.apertures["Primary"].raw_flux)
        arrays["raw_flux_error_3x3"].append(
            lightcurve.apertures["Primary"].raw_flux_err
        )
        for key in quality_counts:
            quality_counts[key] += int(overlay.counts[key])
    if cache is None:
        reference.assert_unchanged()
    raw_arrays = {
        "time": _concatenate(arrays["time"], dtype=np.float64),
        "cadence": _concatenate(arrays["cadence"], dtype=np.int64),
        "orbit": _concatenate(arrays["orbit"], dtype=np.int64),
        "internal_quality": _concatenate(arrays["internal_quality"], dtype=np.uint64),
        "spoc_quality": _concatenate(arrays["spoc_quality"], dtype=np.uint64),
        "qlp_quality": _concatenate(arrays["qlp_quality"], dtype=np.uint64),
        "authority_excluded": _concatenate(
            arrays["authority_excluded"], dtype=np.bool_
        ),
        "raw_flux_1x1": _concatenate(arrays["raw_flux_1x1"], dtype=np.float64),
        "raw_flux_error_1x1": _concatenate(
            arrays["raw_flux_error_1x1"], dtype=np.float64
        ),
        "raw_flux_3x3": _concatenate(arrays["raw_flux_3x3"], dtype=np.float64),
        "raw_flux_error_3x3": _concatenate(
            arrays["raw_flux_error_3x3"], dtype=np.float64
        ),
    }
    return A2V1RawObservation(
        raw_arrays=raw_arrays,
        source_sha256=computed_source_hash,
        audit={
            "input_adapter": A2V1_HDF5_ADAPTER_NAME,
            "scientific_training_eligible": True,
            "hdf5_paths": [str(path) for path in hdf5_paths],
            "hdf5_sha256s": normalized_hashes,
            "quality_reference": reference.provenance,
            "quality_counts": quality_counts,
        },
    )


__all__ = [
    "A2V1_HDF5_ADAPTER_NAME",
    "A2V1_HDF5_MANIFEST_COLUMNS",
    "A2V1_HDF5_MANIFEST_FIELDS",
    "A2V1_HDF5_SOURCE_INVENTORY_COLUMNS",
    "A2V1_HDF5_SOURCE_INVENTORY_FIELDS",
    "A2V1AdapterCache",
    "A2V1RawObservation",
    "a2v1_source_sha256",
    "bind_a2v1_hdf5_source_inventory",
    "load_a2v1_hdf5_observation",
    "write_a2v1_hdf5_bound_manifest",
]
