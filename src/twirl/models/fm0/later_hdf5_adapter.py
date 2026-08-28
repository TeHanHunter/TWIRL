"""Provider-neutral retained-HDF5 adapter for later TWIRL-FM sectors.

This is an additive v2 path.  It deliberately does not widen or reinterpret
the frozen FM0.1 A2v1 adapter, whose persisted quality field is correctly
named ``spoc_quality``.  Later sectors use the provider-neutral
``mission_quality`` reference so S66 records SPOC and S67+ records TICA.

Only already selected train/development rows may reach this module.  Every
HDF5 is rehashed against the immutable retained-cell output manifest before
the light curve is opened.  Sealed rows are rejected before any path is
resolved, hashed, or opened.
"""
from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath
from typing import Any

import numpy as np

from twirl.lightcurves.mission_quality_reference import MissionQualityReference
from twirl.lightcurves.tglc_h5_reader import read_tglc_h5

from .later_source_inventory import (
    SOURCE_FIELDS,
    _cell_inventory,
    _verify_source_receipt,
)
from .registry import FM0ContractError, assert_no_search_columns, sha256_file

LATER_HDF5_ADAPTER_NAME = "a2v1_retained_hdf5_mission_quality_aware_v2"
LATER_ALLOWED_SOURCE_PARTITIONS = frozenset({"poc_train", "poc_development"})
LATER_FORBIDDEN_SOURCE_PARTITION = "poc_sealed_test"
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_CELL = re.compile(
    r"^s(?P<sector>[0-9]{4})_o(?P<orbit>[0-9]+)_cam(?P<camera>[1-4])_ccd(?P<ccd>[1-4])$"
)
_EVIDENCE_FIELDS = frozenset(
    {
        "orbit",
        "camera",
        "ccd",
        "cell",
        "hdf5_sha256",
        "cell_manifest_sha256",
        "output_manifest_sha256",
        "hdf5_path",
        "output_manifest_path",
    }
)


@dataclass(frozen=True)
class LaterRawObservation:
    """One checksum-bound full sector visit ready for six-view derivation."""

    raw_arrays: Mapping[str, np.ndarray]
    source_sha256: str
    audit: Mapping[str, Any]


@dataclass(frozen=True)
class LaterCellAuthority:
    """Once-validated immutable authority for one retained production cell."""

    cell: str
    orbit: int
    camera: int
    ccd: int
    retained_root: Path
    cell_manifest_path: Path
    cell_manifest_sha256: str
    output_manifest_path: Path
    output_manifest_sha256: str
    hdf5_sha256_by_path: Mapping[Path, str]


@dataclass
class LaterAdapterCache:
    """Sector-scoped retained-cell and mission-quality authority cache.

    Cell manifests, output manifests, and completion/retention closure are
    validated exactly once while this cache is constructed.  Individual HDF5
    files are deliberately excluded from ``tracked_files``: each is instead
    rehashed immediately before its one content read.  The release builder
    calls :meth:`assert_unchanged` once after the final visit is built.
    """

    sector: int
    expected_orbits: tuple[int, int]
    reference: MissionQualityReference
    cells: dict[str, LaterCellAuthority]
    tracked_files: dict[Path, str] = field(default_factory=dict)
    tracked_markers: dict[Path, str] = field(default_factory=dict)
    final_check_performed: bool = False

    def track_file(self, path: str | Path, sha256: str) -> None:
        resolved = Path(path).expanduser().resolve()
        expected = _digest(sha256, label=f"cached authority SHA-256 for {resolved}")
        prior = self.tracked_files.setdefault(resolved, expected)
        if prior != expected:
            raise FM0ContractError(f"conflicting cached authority digest: {resolved}")

    def track_marker(self, path: str | Path, expected_text: str) -> None:
        resolved = Path(path).expanduser().resolve()
        prior = self.tracked_markers.setdefault(resolved, str(expected_text))
        if prior != str(expected_text):
            raise FM0ContractError(f"conflicting cached authority marker: {resolved}")

    def assert_unchanged(self) -> None:
        """Rehash every cached authority once, after all source reads."""

        if self.final_check_performed:
            raise FM0ContractError("later adapter cache final check ran more than once")
        try:
            self.reference.assert_unchanged()
            for path, expected in self.tracked_files.items():
                if path.is_symlink() or not path.is_file() or sha256_file(path) != expected:
                    raise FM0ContractError(f"cached later authority changed: {path}")
            for path, expected in self.tracked_markers.items():
                if (
                    path.is_symlink()
                    or not path.is_file()
                    or path.read_text(encoding="utf-8").strip() != expected
                ):
                    raise FM0ContractError(f"cached later authority marker changed: {path}")
        except FM0ContractError:
            raise
        except Exception as exc:
            raise FM0ContractError("mission-quality authority changed during build") from exc
        self.final_check_performed = True


def _digest(value: Any, *, label: str) -> str:
    result = str(value).strip().lower()
    if _SHA256.fullmatch(result) is None:
        raise FM0ContractError(f"{label} must be a lowercase SHA-256 digest")
    return result


def _json_list(value: Any, *, label: str) -> list[Any]:
    try:
        result = json.loads(str(value))
    except (TypeError, json.JSONDecodeError) as exc:
        raise FM0ContractError(f"{label} must be a JSON list") from exc
    if not isinstance(result, list) or not result:
        raise FM0ContractError(f"{label} must be a nonempty JSON list")
    return result


def _positive_identifier(value: Any, *, label: str) -> str:
    result = str(value).strip()
    if not result.isdigit() or int(result) <= 0:
        raise FM0ContractError(f"{label} must be a positive decimal integer")
    return str(int(result))


def _materialized_file(value: Any, *, label: str) -> Path:
    text = str(value).strip()
    if not text:
        raise FM0ContractError(f"{label} is required")
    declared = Path(text).expanduser()
    if not declared.is_absolute():
        raise FM0ContractError(f"{label} must be absolute")
    if declared.is_symlink() or not declared.is_file() or declared.stat().st_size <= 0:
        raise FM0ContractError(
            f"{label} must be a nonempty materialized file: {declared}"
        )
    return declared.resolve()


def _output_manifest_entries(path: Path) -> dict[str, str]:
    entries: dict[str, str] = {}
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as exc:
        raise FM0ContractError(f"cannot read retained output manifest: {path}") from exc
    for line_number, raw in enumerate(lines, start=1):
        fields = raw.split(None, 1)
        if len(fields) != 2:
            raise FM0ContractError(
                f"malformed retained output manifest row {path}:{line_number}"
            )
        digest = _digest(fields[0], label="retained output-manifest digest")
        relative = PurePosixPath(fields[1].strip())
        relative_text = relative.as_posix()
        if (
            relative.is_absolute()
            or ".." in relative.parts
            or len(relative.parts) != 2
            or relative.parts[0] not in {"epsf", "LC"}
            or relative_text in entries
        ):
            raise FM0ContractError(
                f"unsafe or duplicate retained output path: {relative_text!r}"
            )
        entries[relative_text] = digest
    if not entries:
        raise FM0ContractError(f"retained output manifest is empty: {path}")
    return entries


def build_later_adapter_cache(
    *,
    source_receipt_path: str | Path,
    expected_source_receipt_sha256: str,
    sector: int,
    expected_orbits: Sequence[int],
    expected_target_authority_sha256: str,
    reference: MissionQualityReference,
) -> LaterAdapterCache:
    """Validate every retained cell authority once for a sector build."""

    sector = int(sector)
    orbits = tuple(int(value) for value in expected_orbits)
    if sector not in range(66, 78) or orbits != (
        2 * sector + 7,
        2 * sector + 8,
    ):
        raise FM0ContractError(
            "later adapter cache is restricted to exact S66--S77 orbit pairs"
        )
    if reference.sector != sector or tuple(reference.expected_orbits) != orbits:
        raise FM0ContractError("mission-quality reference/cache identity drifted")
    target_authority_sha = _digest(
        expected_target_authority_sha256,
        label="expected target-authority SHA-256",
    )
    receipt_sha = _digest(
        expected_source_receipt_sha256,
        label="expected source-receipt SHA-256",
    )
    receipt_path = _materialized_file(
        source_receipt_path, label="ORCD retained source receipt"
    )
    # The caller-authorized digest is checked before the receipt is parsed.
    if sha256_file(receipt_path) != receipt_sha:
        raise FM0ContractError("ORCD retained source-receipt SHA-256 drifted")
    receipt, source_root, bindings = _verify_source_receipt(
        source_receipt_path=receipt_path,
        sector=sector,
        expected_orbits=orbits,
    )
    expected_cells = {
        f"s{sector:04d}_o{orbit}_cam{camera}_ccd{ccd}"
        for orbit in orbits
        for camera in range(1, 5)
        for ccd in range(1, 5)
    }
    if set(bindings) != expected_cells:
        raise FM0ContractError("retained source receipt cell inventory drifted")

    cache = LaterAdapterCache(
        sector=sector,
        expected_orbits=(orbits[0], orbits[1]),
        reference=reference,
        cells={},
    )
    cache.track_file(receipt_path, receipt_sha)
    total_hdf5 = 0
    for cell in sorted(expected_cells):
        binding = bindings[cell]
        inventory = _cell_inventory(
            source_root=source_root,
            binding=binding,
            sector=sector,
            expected_target_authority_sha256=target_authority_sha,
        )
        hdf5_rows = inventory.get("hdf5")
        if not isinstance(hdf5_rows, list) or not hdf5_rows:
            raise FM0ContractError(f"{cell} has no retained HDF5 authority")
        match = _CELL.fullmatch(cell)
        assert match is not None
        orbit = int(match.group("orbit"))
        camera = int(match.group("camera"))
        ccd = int(match.group("ccd"))
        retained_root = Path(str(binding["retained_root"])).resolve()
        manifest_hashes = {
            _digest(
                row["cell_manifest_sha256"],
                label=f"{cell} cell-manifest SHA-256",
            )
            for row in hdf5_rows
        }
        output_paths = {
            Path(str(row["output_manifest_path"])).resolve() for row in hdf5_rows
        }
        output_hashes = {
            _digest(
                row["output_manifest_sha256"],
                label=f"{cell} output-manifest SHA-256",
            )
            for row in hdf5_rows
        }
        if len(manifest_hashes) != 1 or len(output_paths) != 1 or len(output_hashes) != 1:
            raise FM0ContractError(f"{cell} has conflicting cached manifest authority")
        cell_manifest_sha = next(iter(manifest_hashes))
        output_manifest_path = next(iter(output_paths))
        output_manifest_sha = next(iter(output_hashes))
        cell_manifest_path = retained_root / "cell_manifest.json"
        hdf5_by_path: dict[Path, str] = {}
        for row in hdf5_rows:
            path = Path(str(row["hdf5_path"])).resolve()
            digest = _digest(row["hdf5_sha256"], label=f"{cell} HDF5 SHA-256")
            if path in hdf5_by_path:
                raise FM0ContractError(f"{cell} has duplicate retained HDF5 path")
            hdf5_by_path[path] = digest
        total_hdf5 += len(hdf5_by_path)
        cache.cells[cell] = LaterCellAuthority(
            cell=cell,
            orbit=orbit,
            camera=camera,
            ccd=ccd,
            retained_root=retained_root,
            cell_manifest_path=cell_manifest_path,
            cell_manifest_sha256=cell_manifest_sha,
            output_manifest_path=output_manifest_path,
            output_manifest_sha256=output_manifest_sha,
            hdf5_sha256_by_path=hdf5_by_path,
        )
        # These files were all independently parsed/hashed by the source
        # receipt and cell-inventory validators above.  Record, rather than
        # reread, their digests for the one final drift assertion.
        cache.track_file(
            binding["completion_receipt"], binding["completion_receipt_sha256"]
        )
        cache.track_file(
            binding["retention_receipt"], binding["retention_receipt_sha256"]
        )
        cache.track_file(cell_manifest_path, cell_manifest_sha)
        input_root = (
            source_root
            / "inputs"
            / f"s{sector:04d}-o{orbit}-cam{camera}-ccd{ccd}"
        )
        cache.track_file(input_root / "cell_manifest.json", cell_manifest_sha)
        cache.track_file(output_manifest_path, output_manifest_sha)
        cache.track_marker(input_root / "READY", cell_manifest_sha)
        cache.track_marker(retained_root / "READY", cell_manifest_sha)
    if total_hdf5 != int(receipt.get("n_hdf5_products_declared", -1)):
        raise FM0ContractError("cached retained HDF5 count differs from source receipt")
    return cache


def later_input_source_sha256(
    *,
    source_row_sha256: str,
    mission_quality_reference_sha256: str,
    hdf5_quality_receipt_sha256: str,
) -> str:
    """Bind retained photometry, normalized quality, and the openability gate."""

    payload = {
        "adapter": LATER_HDF5_ADAPTER_NAME,
        "source_row_sha256": _digest(
            source_row_sha256, label="source_row_sha256"
        ),
        "mission_quality_reference_sha256": _digest(
            mission_quality_reference_sha256,
            label="mission_quality_reference_sha256",
        ),
        "hdf5_quality_receipt_sha256": _digest(
            hdf5_quality_receipt_sha256,
            label="hdf5_quality_receipt_sha256",
        ),
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _concatenate(parts: Sequence[np.ndarray], *, dtype: Any) -> np.ndarray:
    if not parts:
        raise FM0ContractError("later retained source contains no HDF5 cadence arrays")
    return np.concatenate([np.asarray(value, dtype=dtype) for value in parts])


def load_later_hdf5_observation(
    row: Mapping[str, Any],
    *,
    reference: MissionQualityReference,
    hdf5_quality_receipt_sha256: str,
    cache: LaterAdapterCache | None = None,
) -> LaterRawObservation:
    """Load one full train/development visit through immutable authorities."""

    if tuple(row) != tuple(SOURCE_FIELDS):
        raise FM0ContractError("later source row columns drifted")
    assert_no_search_columns(row.keys(), context="later retained source row")

    # This check is intentionally first.  No sealed HDF5 path is parsed,
    # resolved, hashed, or opened by this adapter.
    partition = str(row.get("source_partition", "")).strip()
    if partition == LATER_FORBIDDEN_SOURCE_PARTITION:
        raise FM0ContractError("sealed-test HDF5 content must remain unopened")
    if partition not in LATER_ALLOWED_SOURCE_PARTITIONS:
        raise FM0ContractError(f"unsupported later source partition: {partition!r}")

    if cache is None:
        raise FM0ContractError("later retained HDF5 reads require a sector adapter cache")

    sector = int(row["sector"])
    if sector != reference.sector:
        raise FM0ContractError("later source sector disagrees with quality reference")
    if (
        cache.sector != sector
        or cache.expected_orbits != tuple(reference.expected_orbits)
        or cache.reference.manifest_sha256 != reference.manifest_sha256
        or cache.final_check_performed
    ):
        raise FM0ContractError("later source disagrees with its sector adapter cache")
    tic_id = _positive_identifier(row["tic_id"], label="later source TIC ID")
    if str(row.get("product_state", "")) != "ORCD_COMPLETE_DEFERRED":
        raise FM0ContractError("later source product state is not ORCD_COMPLETE_DEFERRED")
    if str(row.get("model_outcome_blind", "")).lower() != "true" or str(
        row.get("panel_admission_authorized", "")
    ).lower() != "false":
        raise FM0ContractError("later source model/panel authority flags drifted")

    path_values = _json_list(row["hdf5_paths_json"], label="hdf5_paths_json")
    hash_values = _json_list(row["hdf5_sha256_json"], label="hdf5_sha256_json")
    orbit_values = _json_list(row["orbits_json"], label="orbits_json")
    evidence_values = _json_list(row["hdf5_evidence_json"], label="hdf5_evidence_json")
    declared_count = int(row["n_hdf5_products"])
    if (
        declared_count != 2
        or len(path_values) != declared_count
        or len(hash_values) != declared_count
        or len(orbit_values) != declared_count
        or len(evidence_values) != declared_count
    ):
        raise FM0ContractError("later source must bind exactly two orbit HDF5 products")
    orbits = tuple(int(value) for value in orbit_values)
    if orbits != reference.expected_orbits:
        raise FM0ContractError("later source orbit list disagrees with quality reference")

    arrays: dict[str, list[np.ndarray]] = {
        "time": [],
        "cadence": [],
        "orbit": [],
        "internal_quality": [],
        "mission_quality": [],
        "qlp_quality": [],
        "authority_excluded": [],
        "raw_flux_1x1": [],
        "raw_flux_error_1x1": [],
        "raw_flux_3x3": [],
        "raw_flux_error_3x3": [],
    }
    observed_orbits: list[int] = []
    hdf5_paths: list[str] = []
    hdf5_hashes: list[str] = []
    output_manifest_hashes: list[str] = []
    quality_counts = {
        "n_cadences": 0,
        "n_internal_bad": 0,
        "n_external_bad": 0,
        "n_authority_excluded": 0,
    }
    for path_value, hash_value, declared_orbit, raw_evidence in zip(
        path_values, hash_values, orbit_values, evidence_values, strict=True
    ):
        if not isinstance(raw_evidence, Mapping) or set(raw_evidence) != _EVIDENCE_FIELDS:
            raise FM0ContractError("later HDF5 evidence fields drifted")
        cell = str(raw_evidence["cell"])
        match = _CELL.fullmatch(cell)
        if match is None or int(match.group("sector")) != sector:
            raise FM0ContractError("retained HDF5 evidence has an invalid cell")
        evidence_identity = (
            int(raw_evidence["orbit"]),
            int(raw_evidence["camera"]),
            int(raw_evidence["ccd"]),
        )
        cell_identity = (
            int(match.group("orbit")),
            int(match.group("camera")),
            int(match.group("ccd")),
        )
        if evidence_identity != cell_identity or evidence_identity[0] != int(declared_orbit):
            raise FM0ContractError("retained HDF5 evidence cell identity drifted")

        authority = cache.cells.get(cell)
        if authority is None or evidence_identity != (
            authority.orbit,
            authority.camera,
            authority.ccd,
        ):
            raise FM0ContractError("retained HDF5 cell is absent from cached authority")
        path = _materialized_file(path_value, label="retained HDF5")
        expected_hash = _digest(hash_value, label="retained HDF5 SHA-256")
        if str(raw_evidence["hdf5_path"]) != str(path):
            raise FM0ContractError("retained HDF5 path disagrees with source evidence")
        if (
            _digest(raw_evidence["hdf5_sha256"], label="evidence HDF5 SHA-256")
            != expected_hash
            or authority.hdf5_sha256_by_path.get(path) != expected_hash
        ):
            raise FM0ContractError("retained HDF5 digest disagrees with cached authority")
        output_manifest = Path(str(raw_evidence["output_manifest_path"])).expanduser()
        if not output_manifest.is_absolute():
            raise FM0ContractError("retained output-manifest path must be absolute")
        output_manifest = output_manifest.resolve()
        output_manifest_hash = _digest(
            raw_evidence["output_manifest_sha256"],
            label="retained output-manifest SHA-256",
        )
        cell_manifest_hash = _digest(
            raw_evidence["cell_manifest_sha256"],
            label="retained cell-manifest SHA-256",
        )
        if (
            output_manifest != authority.output_manifest_path
            or output_manifest_hash != authority.output_manifest_sha256
            or cell_manifest_hash != authority.cell_manifest_sha256
        ):
            raise FM0ContractError("retained HDF5 manifest evidence disagrees with cache")

        # Keep this hash directly adjacent to the content open.  Manifests and
        # cell receipts are sector-cached, but every individual HDF5 is always
        # rehashed at its actual read boundary.
        if sha256_file(path) != expected_hash:
            raise FM0ContractError(f"retained HDF5 drifted before read: {path}")

        try:
            lightcurve = read_tglc_h5(path)
        except Exception as exc:
            raise FM0ContractError(f"cannot read retained HDF5 source: {path}") from exc
        if (
            lightcurve.tic != int(tic_id)
            or lightcurve.sector != sector
            or (lightcurve.orbit, lightcurve.cam, lightcurve.ccd) != evidence_identity
        ):
            raise FM0ContractError(f"retained HDF5 identity drifted: {path}")
        if any(name not in lightcurve.apertures for name in ("Small", "Primary")):
            raise FM0ContractError("retained HDF5 lacks Small or Primary aperture")
        if lightcurve.apertures["Small"].flux_was_synthesized or lightcurve.apertures[
            "Primary"
        ].flux_was_synthesized:
            raise FM0ContractError("magnitude-derived flux is not admissible")
        quality = reference.lookup(
            orbit=lightcurve.orbit,
            camera=lightcurve.cam,
            ccd=lightcurve.ccd,
            cadence=lightcurve.cadence,
            context=f"later FM {path.name}",
        )
        arrays["time"].append(lightcurve.time)
        arrays["cadence"].append(lightcurve.cadence)
        arrays["orbit"].append(
            np.full(lightcurve.cadence.size, lightcurve.orbit, dtype=np.int64)
        )
        arrays["internal_quality"].append(lightcurve.quality)
        arrays["mission_quality"].append(quality.mission_quality)
        arrays["qlp_quality"].append(quality.qlp_quality)
        arrays["authority_excluded"].append(quality.authority_excluded)
        arrays["raw_flux_1x1"].append(lightcurve.apertures["Small"].raw_flux)
        arrays["raw_flux_error_1x1"].append(
            lightcurve.apertures["Small"].raw_flux_err
        )
        arrays["raw_flux_3x3"].append(lightcurve.apertures["Primary"].raw_flux)
        arrays["raw_flux_error_3x3"].append(
            lightcurve.apertures["Primary"].raw_flux_err
        )
        observed_orbits.append(lightcurve.orbit)
        hdf5_paths.append(str(path))
        hdf5_hashes.append(expected_hash)
        output_manifest_hashes.append(output_manifest_hash)
        quality_counts["n_cadences"] += int(lightcurve.cadence.size)
        quality_counts["n_internal_bad"] += int(np.count_nonzero(lightcurve.quality))
        quality_counts["n_external_bad"] += int(np.count_nonzero(quality.external_bad))
        quality_counts["n_authority_excluded"] += int(
            np.count_nonzero(quality.authority_excluded)
        )
    if tuple(observed_orbits) != reference.expected_orbits:
        raise FM0ContractError("retained HDF5 products did not preserve orbit order")

    source_sha = later_input_source_sha256(
        source_row_sha256=str(row["source_row_sha256"]),
        mission_quality_reference_sha256=reference.manifest_sha256,
        hdf5_quality_receipt_sha256=hdf5_quality_receipt_sha256,
    )
    raw_arrays = {
        "time": _concatenate(arrays["time"], dtype=np.float64),
        "cadence": _concatenate(arrays["cadence"], dtype=np.int64),
        "orbit": _concatenate(arrays["orbit"], dtype=np.int64),
        "internal_quality": _concatenate(arrays["internal_quality"], dtype=np.uint64),
        "mission_quality": _concatenate(arrays["mission_quality"], dtype=np.uint64),
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
    return LaterRawObservation(
        raw_arrays=raw_arrays,
        source_sha256=source_sha,
        audit={
            "input_adapter": LATER_HDF5_ADAPTER_NAME,
            "mission_quality_provider": reference.provider,
            "mission_quality_composition": (
                "mission_quality | (qlp_quality << 30)"
            ),
            "quality_reference": reference.provenance,
            "hdf5_quality_receipt_sha256": _digest(
                hdf5_quality_receipt_sha256,
                label="hdf5_quality_receipt_sha256",
            ),
            "hdf5_paths": hdf5_paths,
            "hdf5_sha256s": hdf5_hashes,
            "output_manifest_sha256s": output_manifest_hashes,
            "quality_counts": quality_counts,
            "full_visit_shard": True,
            "model_context_length_bound": False,
            "scientific_training_eligible": False,
            "panel_admission_authorized": False,
        },
    )


__all__ = [
    "LATER_ALLOWED_SOURCE_PARTITIONS",
    "LATER_FORBIDDEN_SOURCE_PARTITION",
    "LATER_HDF5_ADAPTER_NAME",
    "LaterAdapterCache",
    "LaterCellAuthority",
    "LaterRawObservation",
    "build_later_adapter_cache",
    "later_input_source_sha256",
    "load_later_hdf5_observation",
]
