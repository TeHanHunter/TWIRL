from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import h5py
import numpy as np
import pytest

import twirl.models.fm0.later_hdf5_adapter as later_adapter
import twirl.models.fm0.later_sector_release as later_release
from twirl.models.fm0.corpus import (
    CORPUS_SELECTION_FIELDS,
    CORPUS_SELECTION_SCHEMA_VERSION,
)
from twirl.models.fm0.hdf5_quality_admission import (
    HDF5_QUALITY_READY_STATE,
    HDF5_QUALITY_RECEIPT_SCHEMA_VERSION,
)
from twirl.models.fm0.input_release import (
    build_mission_quality_observation_release,
    build_observation_release,
    load_input_release,
)
from twirl.models.fm0.later_hdf5_adapter import (
    LATER_HDF5_ADAPTER_NAME,
    LaterAdapterCache,
    LaterCellAuthority,
    load_later_hdf5_observation,
)
from twirl.models.fm0.later_sector_release import (
    LATER_SIX_VIEW_READY_STATE,
    LATER_SIX_VIEW_RECEIPT_SCHEMA_VERSION,
    build_later_sector_release,
)
from twirl.models.fm0.later_source_inventory import (
    IDENTITY_SOURCE_CORPUS_SELECTION,
    SOURCE_FIELDS,
    SOURCE_ROW_SCHEMA_VERSION,
)
from twirl.models.fm0.later_source_inventory import (
    SUMMARY_SCHEMA_VERSION as SOURCE_INVENTORY_SUMMARY_SCHEMA_VERSION,
)
from twirl.models.fm0.registry import (
    FM0ContractError,
    deterministic_source_partition,
)

SECTOR = 67
ORBITS = (141, 142)
PRODUCER = "1" * 40
TARGET_AUTHORITY_SHA = "2" * 64
REFERENCE_MANIFEST_SHA = "3" * 64
REFERENCE_TABLE_SHA = "4" * 64
TRANSFER_MANIFEST_SHA = "5" * 64
EXCLUSION_SHA = "6" * 64


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _canonical(value: object) -> str:
    return hashlib.sha256(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
        ).encode("utf-8")
    ).hexdigest()


def _csv(path: Path, fields: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _component_for(partition: str) -> str:
    for value in range(100_000):
        component = f"leakage_{value:064x}"
        if deterministic_source_partition(component)[0] == partition:
            return component
    raise AssertionError(f"could not find component in {partition}")


def _write_hdf5(
    path: Path,
    *,
    tic: int,
    orbit: int,
    cadence_start: int,
    time_start: float,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 80
    cadence = np.arange(cadence_start, cadence_start + n, dtype=np.int64)
    time = time_start + np.arange(n, dtype=np.float64) * 200.0 / 86400.0
    with h5py.File(path, "w") as handle:
        handle.attrs["TIC ID"] = tic
        handle.attrs["Sector"] = SECTOR
        handle.attrs["Orbit"] = orbit
        handle.attrs["Camera"] = 1
        handle.attrs["CCD"] = 1
        handle.attrs["BJDoffset"] = 2457000
        handle.attrs["TessMag"] = 15.0
        lightcurve = handle.create_group("LightCurve")
        lightcurve.create_dataset("BJD", data=time)
        lightcurve.create_dataset("Cadence", data=cadence)
        lightcurve.create_dataset("QualityFlag", data=np.zeros(n, dtype=np.int64))
        lightcurve.create_dataset("X", data=np.zeros(n))
        lightcurve.create_dataset("Y", data=np.zeros(n))
        background = lightcurve.create_group("Background")
        background.create_dataset("Value", data=np.zeros(n))
        background.create_dataset("Error", data=np.ones(n))
        photometry = lightcurve.create_group("AperturePhotometry")
        for aperture, baseline in (
            ("Small", 1000.0),
            ("Primary", 2000.0),
            ("Large", 3000.0),
        ):
            group = photometry.create_group(f"{aperture}Aperture")
            flux = baseline + 4.0 * np.sin(np.arange(n) / 12.0)
            group.create_dataset("RawFlux", data=flux)
            group.create_dataset("RawFluxError", data=np.full(n, 1.5))
            group.create_dataset("RawMagnitude", data=np.full(n, 15.0))
            group.create_dataset("RawMagnitudeError", data=np.full(n, 0.01))
            group.create_dataset("X", data=np.zeros(n))
            group.create_dataset("Y", data=np.zeros(n))


class _Reference:
    def __init__(self, root: Path) -> None:
        self.sector = SECTOR
        self.provider = "tica"
        self.expected_orbits = ORBITS
        self.reference_dir = root
        self.manifest_path = root / "manifest.json"
        self.table_path = root / "cadence_quality.csv"
        self.manifest_sha256 = _sha(self.manifest_path)
        self.table_sha256 = _sha(self.table_path)
        self.transfer_manifest_sha256 = TRANSFER_MANIFEST_SHA
        self.authority_exclusions_sha256 = EXCLUSION_SHA
        self.assert_unchanged_calls = 0
        self.provenance = {
            "schema_version": "twirl_fm0_mission_quality_reference_v2",
            "mission_quality_provider": "tica",
            "manifest_sha256": self.manifest_sha256,
        }

    def assert_unchanged(self) -> None:
        self.assert_unchanged_calls += 1

    def lookup(self, *, cadence, **_kwargs):
        n = len(cadence)
        mission = np.zeros(n, dtype=np.uint64)
        qlp = np.zeros(n, dtype=np.uint64)
        excluded = np.zeros(n, dtype=bool)
        return SimpleNamespace(
            mission_quality=mission,
            qlp_quality=qlp,
            authority_excluded=excluded,
            external_bad=np.zeros(n, dtype=bool),
        )


def _adapter_cache(
    reference: _Reference, row: dict[str, object]
) -> LaterAdapterCache:
    cells: dict[str, LaterCellAuthority] = {}
    for evidence in json.loads(str(row["hdf5_evidence_json"])):
        cell = str(evidence["cell"])
        retained_root = Path(str(evidence["output_manifest_path"])).parent
        path = Path(str(evidence["hdf5_path"]))
        cells[cell] = LaterCellAuthority(
            cell=cell,
            orbit=int(evidence["orbit"]),
            camera=int(evidence["camera"]),
            ccd=int(evidence["ccd"]),
            retained_root=retained_root,
            cell_manifest_path=retained_root / "cell_manifest.json",
            cell_manifest_sha256=str(evidence["cell_manifest_sha256"]),
            output_manifest_path=Path(str(evidence["output_manifest_path"])),
            output_manifest_sha256=str(evidence["output_manifest_sha256"]),
            hdf5_sha256_by_path={path: str(evidence["hdf5_sha256"])},
        )
    return LaterAdapterCache(
        sector=SECTOR,
        expected_orbits=ORBITS,
        reference=reference,  # type: ignore[arg-type]
        cells=cells,
    )


def _source_row(
    *,
    root: Path,
    gaia: int,
    tic: int,
    partition: str,
    component: str,
    source_receipt: Path,
    identity_binding_sha: str,
    corpus_summary_sha: str,
    alias_authority_sha: str,
) -> dict[str, object]:
    evidence: list[dict[str, object]] = []
    paths: list[str] = []
    hashes: list[str] = []
    for index, orbit in enumerate(ORBITS):
        retained = root / f"{partition}-o{orbit}"
        hdf5 = retained / "LC" / f"{tic}.h5"
        _write_hdf5(
            hdf5,
            tic=tic,
            orbit=orbit,
            cadence_start=800_000 + index * 1000 + tic * 100,
            time_start=3000.0 + index * 15.0 + tic,
        )
        digest = _sha(hdf5)
        output_manifest = retained / "output_manifest.sha256"
        output_manifest.write_text(f"{digest}  LC/{tic}.h5\n", encoding="utf-8")
        evidence.append(
            {
                "orbit": orbit,
                "camera": 1,
                "ccd": 1,
                "cell": f"s{SECTOR:04d}_o{orbit}_cam1_ccd1",
                "hdf5_sha256": digest,
                "cell_manifest_sha256": "7" * 64,
                "output_manifest_sha256": _sha(output_manifest),
                "hdf5_path": str(hdf5.resolve()),
                "output_manifest_path": str(output_manifest.resolve()),
            }
        )
        paths.append(str(hdf5.resolve()))
        hashes.append(digest)
    portable = [
        {key: value for key, value in item.items() if key not in {"hdf5_path", "output_manifest_path"}}
        for item in evidence
    ]
    retained_sha = _canonical(
        {
            "schema_version": "twirl_fm0_later_retained_source_inventory_v1",
            "sector": SECTOR,
            "tic_id": str(tic),
            "target_authority_sha256": TARGET_AUTHORITY_SHA,
            "evidence": portable,
        }
    )
    row_sha = _canonical(
        {
            "retained_source_sha256": retained_sha,
            "gaia_dr3_source_id": str(gaia),
            "tic_id": str(tic),
            "leakage_component_id": component,
            "source_partition": partition,
            "identity_source": IDENTITY_SOURCE_CORPUS_SELECTION,
            "identity_binding_sha256": identity_binding_sha,
            "corpus_summary_sha256": corpus_summary_sha,
            "gaia_tic_alias_authority_sha256": alias_authority_sha,
        }
    )
    return {
        "source_row_schema_version": SOURCE_ROW_SCHEMA_VERSION,
        "sector": SECTOR,
        "gaia_dr3_source_id": str(gaia),
        "tic_id": str(tic),
        "leakage_component_id": component,
        "source_partition": partition,
        "identity_source": IDENTITY_SOURCE_CORPUS_SELECTION,
        "target_authority_sha256": TARGET_AUTHORITY_SHA,
        "corpus_selection_sha256": identity_binding_sha,
        "corpus_summary_sha256": corpus_summary_sha,
        "gaia_tic_alias_authority_sha256": alias_authority_sha,
        "identity_ambiguous": "false",
        "product_state": "ORCD_COMPLETE_DEFERRED",
        "n_hdf5_products": 2,
        "orbits_json": json.dumps(list(ORBITS), separators=(",", ":")),
        "hdf5_paths_json": json.dumps(paths, separators=(",", ":")),
        "hdf5_sha256_json": json.dumps(hashes, separators=(",", ":")),
        "hdf5_evidence_json": json.dumps(evidence, sort_keys=True, separators=(",", ":")),
        "retained_source_sha256": retained_sha,
        "source_row_sha256": row_sha,
        "source_receipt_path": str(source_receipt.resolve()),
        "source_receipt_sha256": _sha(source_receipt),
        "model_outcome_blind": "true",
        "panel_admission_authorized": "false",
    }


def _fixture(tmp_path: Path) -> tuple[Path, Path, Path, _Reference, list[dict[str, object]]]:
    source_receipt = tmp_path / "source_receipt.json"
    source_receipt.write_text("{}\n", encoding="utf-8")
    train_component = _component_for("poc_train")
    sealed_component = _component_for("poc_sealed_test")
    selection_path = tmp_path / "corpus_selection.csv"
    # The inventory rows bind this file hash, so write selection before rows.
    selection_seed = [
        (1001, 41, train_component, "poc_train"),
        (1002, 42, sealed_component, "poc_sealed_test"),
    ]
    selection_rows = [
        {
            "schema_version": CORPUS_SELECTION_SCHEMA_VERSION,
            "gaia_dr3_source_id": str(gaia),
            "tic_id": str(tic),
            "sector": SECTOR,
            "camera": 1,
            "ccd": 1,
            "leakage_component_id": component,
            "source_partition": partition,
            "orbits_json": json.dumps(list(ORBITS), separators=(",", ":")),
            "hdf5_paths_json": json.dumps(
                [
                    f"orbit-{orbit}/ffi/cam1/ccd1/LC/{tic}.h5"
                    for orbit in ORBITS
                ],
                separators=(",", ":"),
            ),
        }
        for gaia, tic, component, partition in selection_seed
    ]
    _csv(selection_path, CORPUS_SELECTION_FIELDS, selection_rows)
    selection_hash = _sha(selection_path)
    alias_authority = tmp_path / "gaia_tic_aliases.csv"
    _csv(
        alias_authority,
        ("gaia_dr3_source_id", "tic_id"),
        [
            {"gaia_dr3_source_id": str(gaia), "tic_id": str(tic)}
            for gaia, tic, _component, _partition in selection_seed
        ],
    )
    alias_authority_hash = _sha(alias_authority)
    corpus_summary = tmp_path / "corpus_summary.json"
    corpus_summary_payload = {
        "schema_version": CORPUS_SELECTION_SCHEMA_VERSION,
        "sectors": list(range(66, 78)),
        "observations_by_sector": {
            str(value): (len(selection_seed) if value == SECTOR else 1)
            for value in range(66, 78)
        },
        "outputs_sha256": {"selection.csv": selection_hash},
        "input_authorities": {
            "observation_fits": {"path": "/test/observations.fits", "sha256": TARGET_AUTHORITY_SHA},
            "gaia_tic_alias_table": {
                "path": str(alias_authority.resolve()),
                "sha256": alias_authority_hash,
            },
        },
    }
    corpus_summary.write_text(
        json.dumps(corpus_summary_payload, indent=2, sort_keys=True) + "\n"
    )
    corpus_summary_hash = _sha(corpus_summary)
    rows = [
        _source_row(
            root=tmp_path / "retained",
            gaia=gaia,
            tic=tic,
            partition=partition,
            component=component,
            source_receipt=source_receipt,
            identity_binding_sha=selection_hash,
            corpus_summary_sha=corpus_summary_hash,
            alias_authority_sha=alias_authority_hash,
        )
        for gaia, tic, component, partition in selection_seed
    ]
    inventory = tmp_path / "inventory"
    sources_path = inventory / "sources.csv"
    _csv(sources_path, SOURCE_FIELDS, rows)
    summary = {
        "summary_schema_version": SOURCE_INVENTORY_SUMMARY_SCHEMA_VERSION,
        "sector": SECTOR,
        "selection_identity_source": IDENTITY_SOURCE_CORPUS_SELECTION,
        "hdf5_content_opened": False,
        "sealed_hdf5_content_opened": False,
        "model_outcome_blind": True,
        "six_view_shards_verified": False,
        "panel_admission_authorized": False,
        "sources_sha256": _sha(sources_path),
        "n_source_rows": 2,
        "source_rows_by_partition": {
            "poc_development": 0,
            "poc_sealed_test": 1,
            "poc_train": 1,
        },
        "target_authority_sha256": TARGET_AUTHORITY_SHA,
        "identity_binding_path": str(selection_path.resolve()),
        "identity_binding_sha256": selection_hash,
        "corpus_summary_path": str(corpus_summary.resolve()),
        "corpus_summary_sha256": corpus_summary_hash,
        "gaia_tic_alias_authority_path": str(alias_authority.resolve()),
        "gaia_tic_alias_authority_sha256": alias_authority_hash,
        "gaia_tic_alias_authority_file_verified": True,
        "source_receipt_path": str(source_receipt.resolve()),
        "source_receipt_sha256": _sha(source_receipt),
        "n_hdf5_products_declared": 4,
    }
    summary_path = inventory / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    (inventory / "READY").write_text(_sha(summary_path) + "\n")

    reference_root = tmp_path / "quality_reference"
    reference_root.mkdir()
    (reference_root / "manifest.json").write_text("{}\n")
    (reference_root / "cadence_quality.csv").write_text("stub\n")
    reference = _Reference(reference_root)
    hdf5_receipt = tmp_path / "hdf5_quality_receipt.json"
    receipt = {
        "schema_version": HDF5_QUALITY_RECEIPT_SCHEMA_VERSION,
        "quality_state": HDF5_QUALITY_READY_STATE,
        "passed": True,
        "sector": SECTOR,
        "expected_orbits": list(ORBITS),
        "mission_quality_type": "tica",
        "hdf5_openability_verified": True,
        "internal_cadence_quality_verified": True,
        "external_cadence_quality_verified": True,
        "six_view_shards_verified": False,
        "panel_admission_authorized": False,
        "n_unreadable_hdf5": 0,
        "n_hdf5_opened": 4,
        "n_hdf5_products": 4,
        "source_receipt_path": str(source_receipt.resolve()),
        "source_receipt_sha256": _sha(source_receipt),
        "quality_transfer_manifest_sha256": TRANSFER_MANIFEST_SHA,
        "mission_quality_authority_exclusions_sha256": EXCLUSION_SHA,
    }
    hdf5_receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return inventory, reference_root, hdf5_receipt, reference, rows


def _release_kwargs(
    *,
    tmp_path: Path,
    inventory: Path,
    reference_root: Path,
    reference: _Reference,
    hdf5_receipt: Path,
) -> dict[str, object]:
    return {
        "sector": SECTOR,
        "source_inventory_dir": inventory,
        "expected_source_inventory_summary_sha256": _sha(
            inventory / "summary.json"
        ),
        "mission_quality_reference_dir": reference_root,
        "expected_mission_quality_reference_manifest_sha256": (
            reference.manifest_sha256
        ),
        "hdf5_quality_receipt_path": hdf5_receipt,
        "expected_hdf5_quality_receipt_sha256": _sha(hdf5_receipt),
        "output_dir": tmp_path / "release",
        "producer_git_sha": PRODUCER,
    }


def test_later_sector_release_excludes_sealed_before_hdf5_read_and_persists_tica(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    inventory, reference_root, hdf5_receipt, reference, rows = _fixture(tmp_path)
    monkeypatch.setattr(
        later_release,
        "load_mission_quality_reference",
        lambda **_kwargs: reference,
    )
    cache = _adapter_cache(reference, rows[0])
    monkeypatch.setattr(
        later_release,
        "build_later_adapter_cache",
        lambda **_kwargs: cache,
    )
    sealed = rows[1]
    for value in json.loads(str(sealed["hdf5_paths_json"])):
        Path(value).unlink()
    for value in json.loads(str(sealed["hdf5_evidence_json"])):
        Path(value["output_manifest_path"]).unlink()

    output = tmp_path / "release"
    receipt = build_later_sector_release(
        sector=SECTOR,
        source_inventory_dir=inventory,
        expected_source_inventory_summary_sha256=_sha(inventory / "summary.json"),
        mission_quality_reference_dir=reference_root,
        expected_mission_quality_reference_manifest_sha256=(
            reference.manifest_sha256
        ),
        hdf5_quality_receipt_path=hdf5_receipt,
        expected_hdf5_quality_receipt_sha256=_sha(hdf5_receipt),
        output_dir=output,
        producer_git_sha=PRODUCER,
    )

    assert receipt["schema_version"] == LATER_SIX_VIEW_RECEIPT_SCHEMA_VERSION
    assert receipt["release_state"] == LATER_SIX_VIEW_READY_STATE
    assert receipt["expected_orbits"] == [141, 142]
    assert receipt["n_observations"] == receipt["n_shards"] == 1
    assert receipt["flux_view_names"] == [
        "raw_relative_1x1",
        "raw_relative_3x3",
        "adp_1x1",
        "adp_3x3",
        "adp015_1x1",
        "adp015_3x3",
    ]
    assert receipt["mission_quality_provider"] == "tica"
    assert receipt["sealed_hdf5_content_opened"] is False
    assert receipt["sealed_shards_written"] is False
    assert receipt["model_outcome_blind"] is True
    assert receipt["scientific_training_eligible"] is False
    assert receipt["panel_admission_authorized"] is False
    assert receipt["full_visit_shards"] is True
    assert receipt["model_context_length_bound"] is False
    assert receipt["source_inventory_summary_path"] == str(
        (inventory / "summary.json").resolve()
    )
    assert receipt["source_inventory_summary_sha256"] == _sha(
        inventory / "summary.json"
    )
    assert receipt["mission_quality_reference_manifest_path"] == str(
        reference.manifest_path
    )
    assert receipt["mission_quality_reference_manifest_sha256"] == (
        reference.manifest_sha256
    )
    assert receipt["hdf5_quality_receipt_path"] == str(hdf5_receipt.resolve())
    assert receipt["hdf5_quality_receipt_sha256"] == _sha(hdf5_receipt)
    assert receipt["six_view_shards_verified"] is True
    assert (output / "READY").read_text().strip() == _sha(output / "receipt.json")

    shards = list((output / "shards").glob("*.npz"))
    assert len(shards) == 1
    assert load_input_release(shards[0]).n_cadences == 160
    manifest_text = (output / "manifest.csv").read_text()
    receipt_text = (output / "receipt.json").read_text()
    assert "mission_quality_provider" in manifest_text
    assert "tica" in manifest_text
    assert "spoc_quality" not in manifest_text + receipt_text
    assert "poc_sealed_test" not in manifest_text
    _path, _loaded, summary = later_release.validate_later_sector_release(output)
    assert summary["n_observations"] == 1
    assert reference.assert_unchanged_calls == 1
    assert all(
        path.stat().st_mode & 0o222 == 0
        for path in (output, *output.rglob("*"))
    )


def test_adapter_rejects_sealed_row_before_touching_missing_paths(tmp_path: Path) -> None:
    _inventory, _reference_root, hdf5_receipt, reference, rows = _fixture(tmp_path)
    sealed = rows[1]
    with pytest.raises(FM0ContractError, match="sealed-test HDF5 content"):
        load_later_hdf5_observation(
            sealed,
            reference=reference,
            hdf5_quality_receipt_sha256=_sha(hdf5_receipt),
            cache=None,
        )


def test_adapter_rehashes_hdf5_against_output_manifest_before_open(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _inventory, _reference_root, hdf5_receipt, reference, rows = _fixture(tmp_path)
    train = rows[0]
    first = Path(json.loads(str(train["hdf5_paths_json"]))[0])
    with first.open("ab") as handle:
        handle.write(b"drift")
    opened = False

    def forbidden_open(_path):
        nonlocal opened
        opened = True
        raise AssertionError("drifted HDF5 must not reach the reader")

    monkeypatch.setattr(later_adapter, "read_tglc_h5", forbidden_open)
    with pytest.raises(FM0ContractError, match="drifted before read"):
        load_later_hdf5_observation(
            train,
            reference=reference,
            hdf5_quality_receipt_sha256=_sha(hdf5_receipt),
            cache=_adapter_cache(reference, train),
        )
    assert opened is False


def test_provider_neutral_entry_point_preserves_frozen_numerical_core() -> None:
    n = 80
    common = {
        "time": np.linspace(3000.0, 3000.1, n),
        "cadence": np.arange(800_000, 800_000 + n),
        "orbit": np.full(n, 141),
        "internal_quality": np.zeros(n, dtype=np.uint64),
        "qlp_quality": np.zeros(n, dtype=np.uint64),
        "authority_excluded": np.zeros(n, dtype=bool),
        "raw_flux_1x1": 1000.0 + np.sin(np.arange(n) / 10),
        "raw_flux_error_1x1": np.ones(n),
        "raw_flux_3x3": 2000.0 + np.sin(np.arange(n) / 10),
        "raw_flux_error_3x3": np.ones(n),
    }
    v1 = build_observation_release(
        {**common, "spoc_quality": np.zeros(n, dtype=np.uint64)}
    )
    v2 = build_mission_quality_observation_release(
        {**common, "mission_quality": np.zeros(n, dtype=np.uint64)},
        mission_quality_provider="tica",
        input_adapter=LATER_HDF5_ADAPTER_NAME,
    )
    for name in (
        "flux",
        "flux_valid",
        "flux_error",
        "error_valid",
        "local_time_cadences",
        "delta_time_cadences",
        "time_valid",
        "segment_boundary",
        "segment_id",
        "view_present",
    ):
        assert np.array_equal(getattr(v1, name), getattr(v2, name))
    assert v1.audit["external_quality_formula"].startswith("spoc_quality")
    assert v2.audit["external_quality_formula"].startswith("mission_quality")
    assert v2.audit["mission_quality_provider"] == "tica"


def test_sector_cache_validates_each_cell_once_and_quality_once(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reference_root = tmp_path / "reference"
    reference_root.mkdir()
    (reference_root / "manifest.json").write_text("{}\n")
    (reference_root / "cadence_quality.csv").write_text("stub\n")
    reference = _Reference(reference_root)
    source_receipt = tmp_path / "source.json"
    source_receipt.write_text("{}\n")
    source_root = tmp_path / "sector"
    bindings: dict[str, dict[str, object]] = {}
    for orbit in ORBITS:
        for camera in range(1, 5):
            for ccd in range(1, 5):
                cell = f"s{SECTOR:04d}_o{orbit}_cam{camera}_ccd{ccd}"
                cell_root = source_root / "outputs" / cell
                bindings[cell] = {
                    "cell": cell,
                    "completion_receipt": str(cell_root / "complete.json"),
                    "completion_receipt_sha256": "1" * 64,
                    "retention_receipt": str(cell_root / "retention.json"),
                    "retention_receipt_sha256": "2" * 64,
                    "retained_root": str(cell_root / "retained"),
                    "n_hdf5_products_declared": 1,
                }
    receipt = {"n_hdf5_products_declared": len(bindings)}
    verify_calls = 0
    cell_calls: dict[str, int] = {}

    def fake_verify(**_kwargs):
        nonlocal verify_calls
        verify_calls += 1
        return receipt, source_root, bindings

    def fake_cell_inventory(*, binding, **_kwargs):
        cell = str(binding["cell"])
        cell_calls[cell] = cell_calls.get(cell, 0) + 1
        match = later_adapter._CELL.fullmatch(cell)
        assert match is not None
        retained = Path(str(binding["retained_root"]))
        hdf5 = retained / "LC" / f"{len(cell_calls)}.h5"
        return {
            "cell": cell,
            "hdf5": [
                {
                    "hdf5_path": str(hdf5),
                    "hdf5_sha256": "3" * 64,
                    "cell_manifest_sha256": "4" * 64,
                    "output_manifest_path": str(
                        retained / "output_manifest.sha256"
                    ),
                    "output_manifest_sha256": "5" * 64,
                }
            ],
        }

    monkeypatch.setattr(later_adapter, "_verify_source_receipt", fake_verify)
    monkeypatch.setattr(later_adapter, "_cell_inventory", fake_cell_inventory)
    cache = later_adapter.build_later_adapter_cache(
        source_receipt_path=source_receipt,
        expected_source_receipt_sha256=_sha(source_receipt),
        sector=SECTOR,
        expected_orbits=ORBITS,
        expected_target_authority_sha256=TARGET_AUTHORITY_SHA,
        reference=reference,  # type: ignore[arg-type]
    )
    assert verify_calls == 1
    assert len(cache.cells) == 32
    assert set(cell_calls.values()) == {1}
    # The synthetic paths test construction call counts; clear them so the
    # final assertion isolates the mission-quality authority call count.
    cache.tracked_files.clear()
    cache.tracked_markers.clear()
    cache.assert_unchanged()
    assert reference.assert_unchanged_calls == 1


def test_cache_final_assert_detects_authority_drift(
    tmp_path: Path,
) -> None:
    reference_root = tmp_path / "reference"
    reference_root.mkdir()
    (reference_root / "manifest.json").write_text("{}\n")
    (reference_root / "cadence_quality.csv").write_text("stub\n")
    reference = _Reference(reference_root)
    authority = tmp_path / "authority.json"
    authority.write_text("before\n")
    cache = LaterAdapterCache(
        sector=SECTOR,
        expected_orbits=ORBITS,
        reference=reference,  # type: ignore[arg-type]
        cells={},
    )
    cache.track_file(authority, _sha(authority))
    authority.write_text("after\n")
    with pytest.raises(FM0ContractError, match="cached later authority changed"):
        cache.assert_unchanged()
    assert reference.assert_unchanged_calls == 1


@pytest.mark.parametrize("sector", [65, 78])
def test_later_release_python_api_rejects_sectors_outside_66_77(
    tmp_path: Path, sector: int
) -> None:
    with pytest.raises(FM0ContractError, match="restricted to S66--S77"):
        build_later_sector_release(
            sector=sector,
            source_inventory_dir=tmp_path / "inventory",
            expected_source_inventory_summary_sha256="0" * 64,
            mission_quality_reference_dir=tmp_path / "reference",
            expected_mission_quality_reference_manifest_sha256="1" * 64,
            hdf5_quality_receipt_path=tmp_path / "hdf5.json",
            expected_hdf5_quality_receipt_sha256="2" * 64,
            output_dir=tmp_path / "release",
            producer_git_sha=PRODUCER,
        )


def test_caller_hashes_fail_before_bound_json_or_reference_reads(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    inventory, reference_root, hdf5_receipt, reference, _rows = _fixture(tmp_path)
    kwargs = _release_kwargs(
        tmp_path=tmp_path,
        inventory=inventory,
        reference_root=reference_root,
        reference=reference,
        hdf5_receipt=hdf5_receipt,
    )

    original_read_json = later_release._read_json
    monkeypatch.setattr(
        later_release,
        "_read_json",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("summary must not be parsed")
        ),
    )
    with pytest.raises(FM0ContractError, match="summary caller binding"):
        build_later_sector_release(
            **{**kwargs, "expected_source_inventory_summary_sha256": "0" * 64}
        )
    monkeypatch.setattr(later_release, "_read_json", original_read_json)

    monkeypatch.setattr(
        later_release,
        "load_mission_quality_reference",
        lambda **_kwargs: (_ for _ in ()).throw(
            AssertionError("quality reference must not be parsed")
        ),
    )
    with pytest.raises(FM0ContractError, match="reference manifest caller binding"):
        build_later_sector_release(
            **{
                **kwargs,
                "expected_mission_quality_reference_manifest_sha256": "0" * 64,
            }
        )

    monkeypatch.setattr(
        later_release,
        "load_mission_quality_reference",
        lambda **_kwargs: reference,
    )

    def guarded_read_json(path, **read_kwargs):
        if Path(path).resolve() == hdf5_receipt.resolve():
            raise AssertionError("HDF5-quality receipt must not be parsed")
        return original_read_json(path, **read_kwargs)

    monkeypatch.setattr(later_release, "_read_json", guarded_read_json)
    with pytest.raises(FM0ContractError, match="receipt caller binding"):
        build_later_sector_release(
            **{**kwargs, "expected_hdf5_quality_receipt_sha256": "0" * 64}
        )


def test_failed_build_cleans_partial_and_can_retry(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    inventory, reference_root, hdf5_receipt, reference, rows = _fixture(tmp_path)
    cache = _adapter_cache(reference, rows[0])
    monkeypatch.setattr(
        later_release,
        "load_mission_quality_reference",
        lambda **_kwargs: reference,
    )
    monkeypatch.setattr(
        later_release,
        "build_later_adapter_cache",
        lambda **_kwargs: cache,
    )
    kwargs = _release_kwargs(
        tmp_path=tmp_path,
        inventory=inventory,
        reference_root=reference_root,
        reference=reference,
        hdf5_receipt=hdf5_receipt,
    )
    original_npz = later_release.deterministic_npz_bytes
    monkeypatch.setattr(
        later_release,
        "deterministic_npz_bytes",
        lambda _release: (_ for _ in ()).throw(RuntimeError("injected failure")),
    )
    with pytest.raises(RuntimeError, match="injected failure"):
        build_later_sector_release(**kwargs)
    output = Path(str(kwargs["output_dir"]))
    assert not output.exists()
    assert not output.with_name(output.name + ".partial").exists()

    monkeypatch.setattr(later_release, "deterministic_npz_bytes", original_npz)
    receipt = build_later_sector_release(**kwargs)
    assert receipt["release_state"] == LATER_SIX_VIEW_READY_STATE
    assert output.is_dir()


def test_full_bundle_validator_rejects_missing_tampered_extra_and_wrong_provider(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    inventory, reference_root, hdf5_receipt, reference, rows = _fixture(tmp_path)
    monkeypatch.setattr(
        later_release,
        "load_mission_quality_reference",
        lambda **_kwargs: reference,
    )
    monkeypatch.setattr(
        later_release,
        "build_later_adapter_cache",
        lambda **_kwargs: _adapter_cache(reference, rows[0]),
    )
    kwargs = _release_kwargs(
        tmp_path=tmp_path,
        inventory=inventory,
        reference_root=reference_root,
        reference=reference,
        hdf5_receipt=hdf5_receipt,
    )
    build_later_sector_release(**kwargs)
    output = Path(str(kwargs["output_dir"]))
    shard = next((output / "shards").glob("*.npz"))
    shard_bytes = shard.read_bytes()

    (output / "shards").chmod(0o750)
    shard.chmod(0o600)
    shard.unlink()
    (output / "shards").chmod(0o550)
    with pytest.raises(FM0ContractError, match="shard hash drifted"):
        later_release.validate_later_sector_release(output)
    (output / "shards").chmod(0o750)
    shard.write_bytes(shard_bytes)
    shard.chmod(0o440)
    (output / "shards").chmod(0o550)

    shard.chmod(0o600)
    shard.write_bytes(shard_bytes + b"drift")
    shard.chmod(0o440)
    with pytest.raises(FM0ContractError, match="shard hash drifted"):
        later_release.validate_later_sector_release(output)
    shard.chmod(0o600)
    shard.write_bytes(shard_bytes)
    shard.chmod(0o440)

    shards_root = output / "shards"
    shards_root.chmod(0o750)
    extra_shard = shards_root / "extra.npz"
    extra_shard.write_bytes(b"extra")
    extra_shard.chmod(0o440)
    shards_root.chmod(0o550)
    with pytest.raises(FM0ContractError, match="shard closure drifted"):
        later_release.validate_later_sector_release(output)
    shards_root.chmod(0o750)
    extra_shard.unlink()
    shards_root.chmod(0o550)

    output.chmod(0o750)
    extra_root = output / "EXTRA"
    extra_root.write_text("extra\n")
    extra_root.chmod(0o440)
    output.chmod(0o550)
    with pytest.raises(FM0ContractError, match="root closure drifted"):
        later_release.validate_later_sector_release(output)
    output.chmod(0o750)
    extra_root.unlink()
    output.chmod(0o550)

    receipt_path = output / "receipt.json"
    ready_path = output / "READY"
    receipt_bytes = receipt_path.read_bytes()
    ready_bytes = ready_path.read_bytes()
    changed = json.loads(receipt_bytes)
    changed["mission_quality_provider"] = "spoc"
    changed_bytes = (json.dumps(changed, indent=2, sort_keys=True) + "\n").encode()
    receipt_path.chmod(0o600)
    ready_path.chmod(0o600)
    receipt_path.write_bytes(changed_bytes)
    ready_path.write_text(hashlib.sha256(changed_bytes).hexdigest() + "\n")
    receipt_path.chmod(0o440)
    ready_path.chmod(0o440)
    with pytest.raises(FM0ContractError, match="receipt contract drifted"):
        later_release.validate_later_sector_release(output)
    receipt_path.chmod(0o600)
    ready_path.chmod(0o600)
    receipt_path.write_bytes(receipt_bytes)
    ready_path.write_bytes(ready_bytes)
    receipt_path.chmod(0o440)
    ready_path.chmod(0o440)

    manifest = output / "manifest.csv"
    manifest.chmod(0o640)
    with pytest.raises(FM0ContractError, match="not read-only"):
        later_release.validate_later_sector_release(output)
    manifest.chmod(0o440)
