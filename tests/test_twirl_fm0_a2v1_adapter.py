from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import h5py
import numpy as np
import pandas as pd
import pytest

import twirl.models.fm0.a2v1_adapter as adapter
from twirl.models.fm0.a2v1_adapter import (
    A2V1_HDF5_ADAPTER_NAME,
    a2v1_source_sha256,
    bind_a2v1_hdf5_source_inventory,
    load_a2v1_hdf5_observation,
    write_a2v1_hdf5_bound_manifest,
)
from twirl.models.fm0.input_release import (
    load_input_release,
    write_a2v1_hdf5_input_release,
)
from twirl.models.fm0.registry import (
    FM0ContractError,
    build_alias_registry,
    build_observation_registry,
    write_registry_release,
)


def _write_hdf5(
    path: Path,
    *,
    tic: int = 42,
    sector: int = 56,
    orbit: int = 119,
    cadence_start: int = 700000,
    n_cadences: int = 80,
) -> None:
    n = n_cadences
    cadence = np.arange(cadence_start, cadence_start + n, dtype=np.int64)
    time = 2900.0 + np.arange(n, dtype=np.float64) * 200.0 / 86400.0
    with h5py.File(path, "w") as handle:
        handle.attrs["TIC ID"] = tic
        handle.attrs["Sector"] = sector
        handle.attrs["Orbit"] = orbit
        handle.attrs["Camera"] = 4
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
        for aperture, baseline in (("Small", 1000.0), ("Primary", 2000.0), ("Large", 3000.0)):
            group = photometry.create_group(f"{aperture}Aperture")
            flux = baseline + 4.0 * np.sin(np.arange(n) / 12.0)
            group.create_dataset("RawFlux", data=flux)
            group.create_dataset("RawFluxError", data=np.full(n, 1.5))
            group.create_dataset("RawMagnitude", data=np.full(n, 15.0))
            group.create_dataset("RawMagnitudeError", data=np.full(n, 0.01))
            group.create_dataset("X", data=np.zeros(n))
            group.create_dataset("Y", data=np.zeros(n))


class _Reference:
    def __init__(self, table_path: Path) -> None:
        self.table_path = table_path
        self.provenance = {"quality_table": str(table_path), "policy": "test"}

    def assert_unchanged(self) -> None:
        return None

    def apply(
        self,
        *,
        cadenceno,
        internal_quality,
        orbitid,
        orbitid_policy="strict",
        **_kwargs,
    ):
        table = pd.read_csv(self.table_path).set_index("cadenceno")
        cadence = np.asarray(cadenceno, dtype=np.int64)
        input_orbitid = np.asarray(orbitid, dtype=np.int64)
        spoc = np.asarray([table.loc[int(value), "spoc_quality"] for value in cadence], dtype=np.int64)
        qlp = np.asarray([table.loc[int(value), "qlp_quality"] for value in cadence], dtype=np.int64)
        reference_orbitid = np.asarray(
            [
                table.loc[int(value), "reference_orbitid"]
                if "reference_orbitid" in table.columns
                else input_orbitid[index]
                for index, value in enumerate(cadence)
            ],
            dtype=np.int64,
        )
        correction = input_orbitid != reference_orbitid
        if correction.any() and orbitid_policy == "strict":
            raise ValueError("strict test reference rejects an orbit mismatch")
        resolved = reference_orbitid if orbitid_policy == "reference_by_cadence" else input_orbitid
        external = spoc | (qlp << 30)
        internal = np.asarray(internal_quality, dtype=np.int64)
        return SimpleNamespace(
            external_quality=external,
            resolved_orbitid=resolved,
            orbitid_reference_correction_mask=correction,
            counts={
                "n_cad_total": int(cadence.size),
                "n_cad_internal_bad": int(np.count_nonzero(internal)),
                "n_cad_external_bad": int(np.count_nonzero(external)),
                "n_cad_authority_excluded": 0,
            },
        )


def _manifest_row(hdf5: Path, quality_table: Path, quality_manifest: Path, *, source_sha: str) -> dict[str, str]:
    return {
        "observation_key": "unused-by-adapter",
        "product_instance_id": "unused-by-adapter",
        "source_sha256": source_sha,
        "hdf5_paths_json": json.dumps([str(hdf5)]),
        "hdf5_sha256_json": json.dumps([hashlib.sha256(hdf5.read_bytes()).hexdigest()]),
        "quality_table_path": str(quality_table),
        "quality_table_sha256": hashlib.sha256(quality_table.read_bytes()).hexdigest(),
        "quality_manifest_path": str(quality_manifest),
        "quality_manifest_sha256": hashlib.sha256(quality_manifest.read_bytes()).hexdigest(),
    }


def _install_reference(monkeypatch: pytest.MonkeyPatch, table: Path) -> None:
    monkeypatch.setattr(
        adapter,
        "load_external_quality_reference",
        lambda **_kwargs: _Reference(table),
    )


def test_adapter_binds_real_hdf5_content_and_quality_components(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    hdf5 = tmp_path / "source.h5"
    _write_hdf5(hdf5)
    table = tmp_path / "quality.csv"
    pd.DataFrame(
        {
            "sector": 56,
            "camera": 4,
            "ccd": 1,
            "cadenceno": np.arange(700000, 700080),
            "spoc_quality": [0] * 79 + [4],
            "qlp_quality": [0] * 80,
        }
    ).to_csv(table, index=False)
    quality_manifest = tmp_path / "quality.json"
    quality_manifest.write_text("{}\n", encoding="utf-8")
    hdf5_hash = hashlib.sha256(hdf5.read_bytes()).hexdigest()
    table_hash = hashlib.sha256(table.read_bytes()).hexdigest()
    manifest_hash = hashlib.sha256(quality_manifest.read_bytes()).hexdigest()
    source_hash = a2v1_source_sha256(
        sector=56,
        tic_id="42",
        hdf5_sha256s=[hdf5_hash],
        quality_table_sha256=table_hash,
        quality_manifest_sha256=manifest_hash,
    )
    row = _manifest_row(hdf5, table, quality_manifest, source_sha=source_hash)
    _install_reference(monkeypatch, table)
    loaded = load_a2v1_hdf5_observation(
        row, manifest_dir=tmp_path, sector=56, tic_id="42"
    )
    assert loaded.source_sha256 == source_hash
    assert loaded.audit["input_adapter"] == A2V1_HDF5_ADAPTER_NAME
    assert loaded.raw_arrays["raw_flux_1x1"].shape == (80,)
    assert loaded.raw_arrays["spoc_quality"][-1] == 4
    assert not loaded.raw_arrays["authority_excluded"].any()


def test_adapter_reconciles_only_the_bounded_s62_orbit_boundary(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    hdf5 = tmp_path / "s62_source.h5"
    _write_hdf5(
        hdf5,
        sector=62,
        orbit=132,
        cadence_start=766048,
        n_cadences=89,
    )
    table = tmp_path / "quality.csv"
    pd.DataFrame(
        {
            "sector": 62,
            "camera": 4,
            "ccd": 1,
            "cadenceno": np.arange(766048, 766137),
            "spoc_quality": 0,
            "qlp_quality": 0,
            "reference_orbitid": 131,
        }
    ).to_csv(table, index=False)
    quality_manifest = tmp_path / "quality.json"
    quality_manifest.write_text("{}\n", encoding="utf-8")
    hdf5_hash = hashlib.sha256(hdf5.read_bytes()).hexdigest()
    source_hash = a2v1_source_sha256(
        sector=62,
        tic_id="42",
        hdf5_sha256s=[hdf5_hash],
        quality_table_sha256=hashlib.sha256(table.read_bytes()).hexdigest(),
        quality_manifest_sha256=hashlib.sha256(quality_manifest.read_bytes()).hexdigest(),
    )
    _install_reference(monkeypatch, table)

    loaded = load_a2v1_hdf5_observation(
        _manifest_row(hdf5, table, quality_manifest, source_sha=source_hash),
        manifest_dir=tmp_path,
        sector=62,
        tic_id="42",
    )

    assert np.array_equal(loaded.raw_arrays["orbit"], np.full(89, 131))
    reconciliation = loaded.audit["orbitid_reconciliation"]
    assert reconciliation["policy"] == "reference_by_cadence"
    assert reconciliation["n_cadences_reconciled"] == 89
    assert reconciliation["sources"] == [
        {
            "source_file": "s62_source.h5",
            "n_cadences": 89,
            "cadence_min": 766048,
            "cadence_max": 766136,
            "input_orbitid": 132,
            "resolved_orbitid": 131,
        }
    ]


def test_adapter_rejects_an_unbounded_s62_reference_orbit_correction(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    hdf5 = tmp_path / "s62_out_of_bounds.h5"
    _write_hdf5(
        hdf5,
        sector=62,
        orbit=132,
        cadence_start=766137,
    )
    table = tmp_path / "quality.csv"
    pd.DataFrame(
        {
            "sector": 62,
            "camera": 4,
            "ccd": 1,
            "cadenceno": np.arange(766137, 766217),
            "spoc_quality": 0,
            "qlp_quality": 0,
            "reference_orbitid": 131,
        }
    ).to_csv(table, index=False)
    quality_manifest = tmp_path / "quality.json"
    quality_manifest.write_text("{}\n", encoding="utf-8")
    hdf5_hash = hashlib.sha256(hdf5.read_bytes()).hexdigest()
    source_hash = a2v1_source_sha256(
        sector=62,
        tic_id="42",
        hdf5_sha256s=[hdf5_hash],
        quality_table_sha256=hashlib.sha256(table.read_bytes()).hexdigest(),
        quality_manifest_sha256=hashlib.sha256(quality_manifest.read_bytes()).hexdigest(),
    )
    _install_reference(monkeypatch, table)

    with pytest.raises(FM0ContractError, match="bounded 132->131"):
        load_a2v1_hdf5_observation(
            _manifest_row(hdf5, table, quality_manifest, source_sha=source_hash),
            manifest_dir=tmp_path,
            sector=62,
            tic_id="42",
        )


def test_hdf5_release_is_scientific_eligible_only_when_registry_binding_matches(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    hdf5 = tmp_path / "source.h5"
    _write_hdf5(hdf5)
    table = tmp_path / "quality.csv"
    pd.DataFrame(
        {
            "sector": 56,
            "camera": 4,
            "ccd": 1,
            "cadenceno": np.arange(700000, 700080),
            "spoc_quality": 0,
            "qlp_quality": 0,
        }
    ).to_csv(table, index=False)
    quality_manifest = tmp_path / "quality.json"
    quality_manifest.write_text("{}\n", encoding="utf-8")
    hdf5_hash = hashlib.sha256(hdf5.read_bytes()).hexdigest()
    source_hash = a2v1_source_sha256(
        sector=56,
        tic_id="42",
        hdf5_sha256s=[hdf5_hash],
        quality_table_sha256=hashlib.sha256(table.read_bytes()).hexdigest(),
        quality_manifest_sha256=hashlib.sha256(quality_manifest.read_bytes()).hexdigest(),
    )
    aliases = build_alias_registry(
        [{"gaia_dr3_source_id": "123456789012345678", "tic_id": "42"}]
    )
    observations = build_observation_registry(
        [
            {
                "gaia_dr3_source_id": "123456789012345678",
                "tic_id": "42",
                "sector": 56,
                "a2v1_product_version": "A2v1",
                "source_sha256": source_hash,
                "product_state": "A2V1_ACCEPTED",
            }
        ],
        aliases,
    )
    registry = tmp_path / "registry"
    write_registry_release(registry, aliases, observations)
    source_row = _manifest_row(hdf5, table, quality_manifest, source_sha=source_hash)
    source_row["observation_key"] = observations[0]["observation_key"]
    source_row["product_instance_id"] = observations[0]["product_instance_id"]
    manifest = tmp_path / "sources.csv"
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=sorted(source_row))
        writer.writeheader()
        writer.writerow(source_row)
    _install_reference(monkeypatch, table)
    summary = write_a2v1_hdf5_input_release(
        registry_dir=registry,
        hdf5_manifest_path=manifest,
        out_dir=tmp_path / "release",
        visit_timing_path=tmp_path / "release" / "visit_timing.csv",
    )
    assert summary["scientific_training_eligible"] is True
    assert summary["input_adapter"] == A2V1_HDF5_ADAPTER_NAME
    shard = next((tmp_path / "release" / "shards").glob("*.npz"))
    assert load_input_release(shard).flux.shape == (80, 6)
    timing_rows = list(
        csv.DictReader((tmp_path / "release" / "visit_timing.csv").open(newline=""))
    )
    assert timing_rows == [
        {
            "observation_key": observations[0]["observation_key"],
            "physical_source_id": observations[0]["physical_source_id"],
            "absolute_visit_start": "2900.0",
            "absolute_visit_end": "2900.1828703703704",
        }
    ]


def test_adapter_rejects_a_hdf5_manifest_with_an_unbound_source_hash(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    hdf5 = tmp_path / "source.h5"
    _write_hdf5(hdf5)
    table = tmp_path / "quality.csv"
    pd.DataFrame(
        {
            "sector": 56,
            "camera": 4,
            "ccd": 1,
            "cadenceno": np.arange(700000, 700080),
            "spoc_quality": 0,
            "qlp_quality": 0,
        }
    ).to_csv(table, index=False)
    quality_manifest = tmp_path / "quality.json"
    quality_manifest.write_text("{}\n", encoding="utf-8")
    _install_reference(monkeypatch, table)
    with pytest.raises(FM0ContractError, match="does not bind"):
        load_a2v1_hdf5_observation(
            _manifest_row(hdf5, table, quality_manifest, source_sha="a" * 64),
            manifest_dir=tmp_path,
            sector=56,
            tic_id="42",
        )


def test_source_inventory_is_bound_to_registry_ids_only_after_alias_grouping(
    tmp_path: Path,
) -> None:
    hdf5 = tmp_path / "source.h5"
    _write_hdf5(hdf5)
    quality_table = tmp_path / "quality.csv"
    quality_table.write_text("stub\n", encoding="utf-8")
    quality_manifest = tmp_path / "quality.json"
    quality_manifest.write_text("{}\n", encoding="utf-8")
    inventory = [
        {
            "gaia_dr3_source_id": "123456789012345678",
            "tic_id": "42",
            "sector": "56",
            "a2v1_product_version": "A2v1",
            "product_state": "A2V1_ACCEPTED",
            "diagnostic_admission_receipt_path": "",
            "diagnostic_admission_receipt_sha256": "",
            "hdf5_paths_json": json.dumps([str(hdf5)]),
            "hdf5_sha256_json": json.dumps([hashlib.sha256(hdf5.read_bytes()).hexdigest()]),
            "quality_table_path": str(quality_table),
            "quality_table_sha256": hashlib.sha256(quality_table.read_bytes()).hexdigest(),
            "quality_manifest_path": str(quality_manifest),
            "quality_manifest_sha256": hashlib.sha256(quality_manifest.read_bytes()).hexdigest(),
        }
    ]
    aliases = build_alias_registry(
        [{"gaia_dr3_source_id": "123456789012345678", "tic_id": "42"}]
    )
    observations, bound = bind_a2v1_hdf5_source_inventory(inventory, aliases)
    assert observations[0]["source_sha256"] == bound[0]["source_sha256"]
    assert bound[0]["observation_key"] == observations[0]["observation_key"]
    digest = write_a2v1_hdf5_bound_manifest(tmp_path / "bound.csv", bound)
    assert digest == hashlib.sha256((tmp_path / "bound.csv").read_bytes()).hexdigest()
