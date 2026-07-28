from __future__ import annotations

import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest

from twirl.lightcurves.a2v1_cadence_reference import (
    AUTHORITY_EXCLUSION_POLICY_CONTRACT,
)
from twirl.lightcurves.external_quality import (
    EXPECTED_CADENCE_AUTHORITY,
    EXPECTED_QUALITY_AUTHORITY,
    ExternalQualityReference,
)
from twirl.vetting import ssl_full_pool_native as native
from twirl.vetting.harmonic_inputs import read_native_light_curve_from_h5
from twirl.vetting.ssl_full_pool import (
    FULL_POOL_CONTRACT_VERSION,
    FULL_POOL_SUMMARY_SCHEMA_VERSION,
    POOL_COLUMNS,
)
from twirl.vetting.teacher_native_registry import read_table


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _identity_sha256(rows: pd.DataFrame) -> str:
    digest = hashlib.sha256()
    for row in (
        rows[["sector", "tic"]]
        .sort_values(["sector", "tic"])
        .itertuples(index=False)
    ):
        digest.update(
            json.dumps(
                {"sector": int(row.sector), "tic": int(row.tic)},
                sort_keys=True,
                separators=(",", ":"),
            ).encode("ascii")
        )
        digest.update(b"\n")
    return digest.hexdigest()


def _inventory_sha256(tics: list[int]) -> str:
    return hashlib.sha256(
        "".join(f"{tic}\n" for tic in sorted(tics)).encode("ascii")
    ).hexdigest()


def _compact_h5(
    path: Path,
    *,
    sector: int,
    tic: int,
    cadenceno: np.ndarray,
    orbitid: np.ndarray,
) -> None:
    n = len(cadenceno)
    with h5py.File(path, "w") as h5:
        group = h5.create_group(f"targets/{tic:016d}")
        group.attrs["sector"] = sector
        group.attrs["tic"] = tic
        group.attrs["camera"] = 1
        group.attrs["ccd"] = 1
        group.create_dataset("time", data=3000.0 + np.arange(n) / 720.0)
        group.create_dataset("cadenceno", data=cadenceno)
        group.create_dataset("orbitid", data=orbitid)
        group.create_dataset("quality", data=np.zeros(n, dtype=np.int32))
        group.create_dataset(
            "DET_FLUX_ADP_SML", data=np.linspace(0.98, 1.02, n)
        )
        group.create_dataset(
            "DET_FLUX_ADP", data=np.linspace(0.99, 1.01, n)
        )


def _frozen_pool_fixture(
    tmp_path: Path,
    *,
    sector: int = 62,
    tic: int = 420_062,
    cadenceno: np.ndarray | None = None,
    orbitid: np.ndarray | None = None,
) -> dict[str, Path | int | np.ndarray]:
    if cadenceno is None:
        cadenceno = np.asarray([766_047, 766_048, 766_136, 766_228])
    if orbitid is None:
        orbitid = np.asarray([131, 132, 132, 132])
    compact = tmp_path / f"s{sector}_compact.h5"
    _compact_h5(
        compact,
        sector=sector,
        tic=tic,
        cadenceno=cadenceno,
        orbitid=orbitid,
    )
    compact_hash = _sha256(compact)
    rows: list[dict[str, object]] = []
    for current_sector in range(56, 63):
        current_tic = tic if current_sector == sector else 420_000 + current_sector
        selected = current_sector == sector
        rows.append(
            {
                "pool_contract_version": FULL_POOL_CONTRACT_VERSION,
                "observation_id": (
                    f"s{current_sector:04d}-tic{current_tic:016d}"
                ),
                "sector": current_sector,
                "tic": current_tic,
                "camera": 1,
                "ccd": 1,
                "tessmag": 18.0,
                "n_cadences": len(cadenceno),
                "flux_columns_json": (
                    '["DET_FLUX_ADP_SML","DET_FLUX_ADP"]'
                ),
                "compact_h5_path": str(compact if selected else tmp_path / f"s{current_sector}.h5"),
                "compact_h5_sha256": (
                    compact_hash if selected else f"{current_sector % 10}" * 64
                ),
                "compact_h5_size_bytes": compact.stat().st_size if selected else 1,
                "compact_group_path": f"targets/{current_tic:016d}",
                "compact_manifest_path": f"/readonly/s{current_sector}.manifest.json",
                "compact_manifest_sha256": "a" * 64,
                "source_fits": f"/readonly/s{current_sector}/{current_tic}.fits",
            }
        )
    pool = pd.DataFrame(rows, columns=POOL_COLUMNS)
    pool_path = tmp_path / "teacher_ssl_full_pool_observations.csv"
    pool.to_csv(pool_path, index=False, lineterminator="\n")
    allowlist = tmp_path / f"s{sector}_tics.csv"
    allowlist.write_text(f"tic\n{tic}\n", encoding="ascii")
    compact_exports = []
    for row in rows:
        current_sector = int(row["sector"])
        compact_exports.append(
            {
                "sector": current_sector,
                "compact_h5": {
                    "path": str(row["compact_h5_path"]),
                    "size_bytes": int(row["compact_h5_size_bytes"]),
                    "sha256": str(row["compact_h5_sha256"]),
                },
            }
        )
    summary = {
        "passed": True,
        "schema_version": FULL_POOL_SUMMARY_SCHEMA_VERSION,
        "pool_contract_version": FULL_POOL_CONTRACT_VERSION,
        "sectors": list(range(56, 63)),
        "inputs": {"compact_exports": compact_exports},
        "counts": {
            "retained": {
                "n_observations": 7,
                "n_unique_tics": 7,
                "n_multisector_tics": 0,
            }
        },
        "identity_hashes": {
            "retained_observations_sha256": _identity_sha256(pool)
        },
        "leakage_audit": {
            "fixed_test_observations_retained": 0,
            "s63_reserved_observations_retained": 0,
        },
        "outputs": {
            "csv": {
                "path": str(pool_path),
                "size_bytes": pool_path.stat().st_size,
                "sha256": _sha256(pool_path),
                "n_rows": 7,
            },
            "sector_allowlists": {
                str(sector): {
                    "path": str(allowlist),
                    "size_bytes": allowlist.stat().st_size,
                    "sha256": _sha256(allowlist),
                    "n_tics": 1,
                    "tic_inventory_sha256": _inventory_sha256([tic]),
                }
            },
        },
    }
    summary_path = tmp_path / "teacher_ssl_full_pool_manifest.summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {
        "sector": sector,
        "tic": tic,
        "cadenceno": cadenceno,
        "orbitid": orbitid,
        "compact": compact,
        "pool": pool_path,
        "summary": summary_path,
        "allowlist": allowlist,
    }


def _raw_payload(
    cadenceno: np.ndarray,
    orbitid: np.ndarray,
) -> dict[str, np.ndarray]:
    n = len(cadenceno)
    return {
        "time": 2460000.0 + np.arange(n) / 720.0,
        "cadenceno": np.asarray(cadenceno),
        "orbitid": np.asarray(orbitid),
        "quality": np.zeros(n, dtype=np.int32),
        "raw_flux_small": np.linspace(10.0, 11.0, n),
        "raw_flux_err_small": np.full(n, 0.1),
        "raw_flux_primary": np.linspace(20.0, 21.0, n),
        "raw_flux_err_primary": np.full(n, 0.2),
    }


def _export_raw(
    tmp_path: Path,
    inputs: dict[str, Path | int | np.ndarray],
    monkeypatch: pytest.MonkeyPatch,
) -> Path:
    raw_root = tmp_path / "raw"
    sector = int(inputs["sector"])
    orbit_one = 119 + 2 * (sector - 56)
    orbits = (orbit_one, orbit_one + 1)
    for orbit in orbits:
        path = (
            raw_root
            / f"orbit-{orbit}"
            / "ffi/cam1/ccd1/LC"
            / f"{int(inputs['tic'])}.h5"
        )
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
    payload = _raw_payload(
        np.asarray(inputs["cadenceno"]),
        np.asarray(inputs["orbitid"]),
    )
    monkeypatch.setattr(
        native.harmonic_export,
        "merge_tglc_raw_paths",
        lambda paths: payload,
    )
    output = tmp_path / "raw_source.h5"
    result = native.export_full_pool_raw_source_shard(
        sector=int(inputs["sector"]),
        pool_path=Path(inputs["pool"]),
        pool_summary_path=Path(inputs["summary"]),
        allowlist_path=Path(inputs["allowlist"]),
        raw_root=raw_root,
        orbits=orbits,
        out_h5=output,
    )
    assert result["n_shard_observations"] == 1
    return output


def _quality_reference(
    tmp_path: Path,
    *,
    sector: int = 62,
    cadenceno: np.ndarray,
    reference_orbitid: np.ndarray,
) -> ExternalQualityReference:
    table = tmp_path / "cadence_reference.csv"
    manifest = tmp_path / "cadence_reference.json"
    table.write_text("test\n", encoding="ascii")
    manifest.write_text("{}\n", encoding="ascii")
    return ExternalQualityReference(
        sector=sector,
        table_path=table,
        manifest_path=manifest,
        table_sha256=_sha256(table),
        manifest_sha256=_sha256(manifest),
        source_declaration_sha256="b" * 64,
        contract_version=f"s{sector}_a2v1_cadence_reference_v1",
        cadence_authority=EXPECTED_CADENCE_AUTHORITY,
        quality_authority=EXPECTED_QUALITY_AUTHORITY,
        by_detector={
            (sector, 1, 1): (
                np.asarray(cadenceno, dtype=np.int64),
                np.asarray(reference_orbitid, dtype=np.int64),
                np.zeros(len(cadenceno), dtype=np.int64),
            )
        },
        authority_exclusions={},
        authority_exclusions_sha256="c" * 64,
    )


def test_sector_authority_rejects_allowlist_drift(tmp_path: Path) -> None:
    inputs = _frozen_pool_fixture(tmp_path)
    authority = native.load_sector_pool_authority(
        sector=int(inputs["sector"]),
        pool_path=Path(inputs["pool"]),
        pool_summary_path=Path(inputs["summary"]),
        allowlist_path=Path(inputs["allowlist"]),
    )
    assert authority.rows["tic"].tolist() == [int(inputs["tic"])]
    Path(inputs["allowlist"]).write_text(
        f"tic\n{int(inputs['tic']) + 1}\n", encoding="ascii"
    )
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        native.load_sector_pool_authority(
            sector=int(inputs["sector"]),
            pool_path=Path(inputs["pool"]),
            pool_summary_path=Path(inputs["summary"]),
            allowlist_path=Path(inputs["allowlist"]),
        )


def test_raw_export_is_exactly_sharded_and_bound_to_pool(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    inputs = _frozen_pool_fixture(tmp_path)
    output = _export_raw(tmp_path, inputs, monkeypatch)

    with h5py.File(output, "r") as h5:
        assert h5.attrs["real_only"] == 1
        assert h5.attrs["full_pool_raw_source_contract_version"] == (
            native.FULL_POOL_RAW_SOURCE_CONTRACT_VERSION
        )
        assert list(h5["targets"]) == [f"{int(inputs['tic']):016d}"]
        group = h5[f"targets/{int(inputs['tic']):016d}"]
        assert group.attrs["sector"] == 62
        assert group.attrs["camera"] == 1
        assert group.attrs["ccd"] == 1


def test_s62_generic_native_contract_corrects_without_frozen_997_row_binding(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    inputs = _frozen_pool_fixture(tmp_path)
    raw_source = _export_raw(tmp_path, inputs, monkeypatch)
    reference = _quality_reference(
        tmp_path,
        cadenceno=np.asarray(inputs["cadenceno"]),
        reference_orbitid=np.asarray([131, 131, 131, 132]),
    )
    monkeypatch.setattr(
        native,
        "load_external_quality_reference",
        lambda **kwargs: reference,
    )
    monkeypatch.setattr(
        native.harmonic_export,
        "_periodogram_payload",
        lambda **kwargs: {
            "bls_log_period_grid": np.asarray([-1.0, 0.0], dtype=np.float32),
            "bls_power_small": np.asarray([0.0, 1.0], dtype=np.float32),
            "bls_sde_small": np.asarray([0.0, 1.0], dtype=np.float32),
            "bls_power_primary": np.asarray([0.0, 1.0], dtype=np.float32),
            "bls_sde_primary": np.asarray([0.0, 1.0], dtype=np.float32),
        },
    )
    output = tmp_path / "native.h5"
    result = native.build_full_pool_native_shard(
        sector=62,
        pool_path=Path(inputs["pool"]),
        pool_summary_path=Path(inputs["summary"]),
        allowlist_path=Path(inputs["allowlist"]),
        raw_source_h5=raw_source,
        compact_adp_h5=Path(inputs["compact"]),
        cadence_reference_table=reference.table_path,
        cadence_reference_manifest=reference.manifest_path,
        out_h5=output,
        orbitid_policy="reference_by_cadence",
        n_periods=2,
    )

    assert result["verification"]["passed"] is True
    assert result["orbitid_reconciliation"]["n_groups_examined"] == 1
    assert result["orbitid_reconciliation"]["n_groups_corrected"] == 1
    assert result["orbitid_reconciliation"]["n_cad_corrected"] == 2
    with h5py.File(output, "r") as h5:
        assert h5.attrs["contract_version"] == (
            native.FULL_POOL_NATIVE_CONTRACT_VERSION
        )
        group_path = f"targets/{int(inputs['tic']):016d}"
        group = h5[group_path]
        assert group["orbitid"][:].tolist() == [131, 131, 131, 132]
        assert group["orbitid_pre_reference"][:].tolist() == [
            131,
            132,
            132,
            132,
        ]
        assert group[
            "orbitid_reference_correction_mask"
        ][:].tolist() == [0, 1, 1, 0]
        loaded = read_native_light_curve_from_h5(
            h5, group_path=group_path, require_errors=True
        )
        assert loaded.orbitid.tolist() == [131, 131, 131, 132]


def test_s61_strict_native_contract_preserves_orbit_ids(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    inputs = _frozen_pool_fixture(
        tmp_path,
        sector=61,
        tic=420_061,
        cadenceno=np.asarray([100, 101, 102, 103]),
        orbitid=np.asarray([129, 129, 130, 130]),
    )
    raw_source = _export_raw(tmp_path, inputs, monkeypatch)
    reference = _quality_reference(
        tmp_path,
        sector=61,
        cadenceno=np.asarray(inputs["cadenceno"]),
        reference_orbitid=np.asarray(inputs["orbitid"]),
    )
    monkeypatch.setattr(
        native,
        "load_external_quality_reference",
        lambda **kwargs: reference,
    )
    monkeypatch.setattr(
        native.harmonic_export,
        "_periodogram_payload",
        lambda **kwargs: {
            "bls_log_period_grid": np.asarray([-1.0], dtype=np.float32),
            "bls_power_small": np.asarray([1.0], dtype=np.float32),
            "bls_sde_small": np.asarray([1.0], dtype=np.float32),
            "bls_power_primary": np.asarray([1.0], dtype=np.float32),
            "bls_sde_primary": np.asarray([1.0], dtype=np.float32),
        },
    )
    output = tmp_path / "s61_native.h5"
    result = native.build_full_pool_native_shard(
        sector=61,
        pool_path=Path(inputs["pool"]),
        pool_summary_path=Path(inputs["summary"]),
        allowlist_path=Path(inputs["allowlist"]),
        raw_source_h5=raw_source,
        compact_adp_h5=Path(inputs["compact"]),
        cadence_reference_table=reference.table_path,
        cadence_reference_manifest=reference.manifest_path,
        out_h5=output,
        orbitid_policy="strict",
        n_periods=1,
    )
    assert result["verification"]["passed"] is True
    assert result["orbitid_reconciliation"]["n_cad_corrected"] == 0
    with h5py.File(output, "r") as h5:
        group = h5[f"targets/{int(inputs['tic']):016d}"]
        assert np.array_equal(
            group["orbitid"][:], group["orbitid_pre_reference"][:]
        )


def test_s62_native_build_rejects_raw_compact_pre_correction_disagreement(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    inputs = _frozen_pool_fixture(tmp_path)
    raw_source = _export_raw(tmp_path, inputs, monkeypatch)
    with h5py.File(raw_source, "a") as h5:
        h5[f"targets/{int(inputs['tic']):016d}/orbitid"][1] = 131
    reference = _quality_reference(
        tmp_path,
        cadenceno=np.asarray(inputs["cadenceno"]),
        reference_orbitid=np.asarray([131, 131, 131, 132]),
    )
    monkeypatch.setattr(
        native,
        "load_external_quality_reference",
        lambda **kwargs: reference,
    )
    with pytest.raises(RuntimeError, match="native export failed"):
        native.build_full_pool_native_shard(
            sector=62,
            pool_path=Path(inputs["pool"]),
            pool_summary_path=Path(inputs["summary"]),
            allowlist_path=Path(inputs["allowlist"]),
            raw_source_h5=raw_source,
            compact_adp_h5=Path(inputs["compact"]),
            cadence_reference_table=reference.table_path,
            cadence_reference_manifest=reference.manifest_path,
            out_h5=tmp_path / "bad_native.h5",
            orbitid_policy="reference_by_cadence",
            n_periods=2,
        )
    failures = pd.read_csv((tmp_path / "bad_native.h5").with_suffix(".failures.csv"))
    assert "raw-source and compact orbitid arrays disagree" in failures.loc[
        0, "error"
    ]


def test_registry_freeze_requires_exact_observation_coverage(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    inputs = _frozen_pool_fixture(tmp_path)
    raw_source = _export_raw(tmp_path, inputs, monkeypatch)
    reference = _quality_reference(
        tmp_path,
        cadenceno=np.asarray(inputs["cadenceno"]),
        reference_orbitid=np.asarray([131, 131, 131, 132]),
    )
    monkeypatch.setattr(
        native,
        "load_external_quality_reference",
        lambda **kwargs: reference,
    )
    monkeypatch.setattr(
        native.harmonic_export,
        "_periodogram_payload",
        lambda **kwargs: {
            "bls_log_period_grid": np.asarray([-1.0], dtype=np.float32),
            "bls_power_small": np.asarray([1.0], dtype=np.float32),
            "bls_sde_small": np.asarray([1.0], dtype=np.float32),
            "bls_power_primary": np.asarray([1.0], dtype=np.float32),
            "bls_sde_primary": np.asarray([1.0], dtype=np.float32),
        },
    )
    output = tmp_path / "native.h5"
    native.build_full_pool_native_shard(
        sector=62,
        pool_path=Path(inputs["pool"]),
        pool_summary_path=Path(inputs["summary"]),
        allowlist_path=Path(inputs["allowlist"]),
        raw_source_h5=raw_source,
        compact_adp_h5=Path(inputs["compact"]),
        cadence_reference_table=reference.table_path,
        cadence_reference_manifest=reference.manifest_path,
        out_h5=output,
        orbitid_policy="reference_by_cadence",
        n_periods=1,
    )
    with pytest.raises(ValueError, match="do not exactly cover"):
        native.write_full_pool_native_registry(
            pool_path=Path(inputs["pool"]),
            pool_summary_path=Path(inputs["summary"]),
            native_shard_paths=[output],
            source_path=tmp_path / "source.csv",
            registry_path=tmp_path / "registry.csv",
            summary_path=tmp_path / "registry.summary.json",
        )

    # A sector-truncated pool is not a valid substitute for the release.
    pool = read_table(Path(inputs["pool"]))
    selected = pool.loc[pool["sector"].eq(62)].reset_index(drop=True)
    selected_path = tmp_path / "selected.csv"
    selected.to_csv(selected_path, index=False, lineterminator="\n")
    summary = json.loads(Path(inputs["summary"]).read_text())
    summary["sectors"] = list(range(56, 63))
    summary["counts"]["retained"] = {
        "n_observations": 1,
        "n_unique_tics": 1,
        "n_multisector_tics": 0,
    }
    summary["identity_hashes"]["retained_observations_sha256"] = (
        _identity_sha256(selected)
    )
    summary["outputs"]["csv"] = {
        "path": str(selected_path),
        "size_bytes": selected_path.stat().st_size,
        "sha256": _sha256(selected_path),
        "n_rows": 1,
    }
    selected_summary = tmp_path / "selected.summary.json"
    selected_summary.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    with pytest.raises(ValueError, match="exact S56--S62"):
        native.write_full_pool_native_registry(
            pool_path=selected_path,
            pool_summary_path=selected_summary,
            native_shard_paths=[output],
            source_path=tmp_path / "source2.csv",
            registry_path=tmp_path / "registry2.csv",
            summary_path=tmp_path / "registry2.summary.json",
        )

    # Minimal contract-carrying shard fixtures are sufficient to exercise the
    # exact observation-keyed registry assembly itself.
    fake_shards: list[Path] = []
    pool_sha = _sha256(Path(inputs["pool"]))
    summary_sha = _sha256(Path(inputs["summary"]))
    for row in pool.itertuples(index=False):
        shard = tmp_path / f"s{int(row.sector)}_native.h5"
        with h5py.File(shard, "w") as h5:
            h5.attrs["contract_version"] = (
                native.FULL_POOL_NATIVE_CONTRACT_VERSION
            )
            h5.attrs["full_pool_sha256"] = pool_sha
            h5.attrs["full_pool_summary_sha256"] = summary_sha
            h5.attrs["sector"] = int(row.sector)
            group = h5.create_group(f"targets/{int(row.tic):016d}")
            group.attrs["sector"] = int(row.sector)
            group.attrs["tic"] = int(row.tic)
        fake_shards.append(shard)

    def fake_verify(path: Path, **kwargs: object) -> dict[str, object]:
        with h5py.File(path, "r") as h5:
            count = len(h5["targets"])
        return {
            "passed": True,
            "failures": [],
            "counts": {"targets": count, "injections": 0},
            "sha256": _sha256(path),
            "size_bytes": path.stat().st_size,
        }

    monkeypatch.setattr(native, "verify_full_pool_native_shard", fake_verify)
    result = native.write_full_pool_native_registry(
        pool_path=Path(inputs["pool"]),
        pool_summary_path=Path(inputs["summary"]),
        native_shard_paths=fake_shards,
        source_path=tmp_path / "source_ok.csv",
        registry_path=tmp_path / "registry_ok.csv",
        summary_path=tmp_path / "registry_ok.summary.json",
    )
    registry = pd.read_csv(tmp_path / "registry_ok.csv")
    assert len(registry) == 7
    assert registry[["sector", "tic"]].drop_duplicates().shape[0] == 7
    assert result["n_observations"] == 7


def test_fullpool_module_keeps_external_quality_contract_constants() -> None:
    assert AUTHORITY_EXCLUSION_POLICY_CONTRACT
    assert native.FULL_POOL_NATIVE_SECTORS == tuple(range(56, 63))
