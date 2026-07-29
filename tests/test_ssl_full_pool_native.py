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
from twirl.lightcurves.detrend_presets import adp03q_config
from twirl.lightcurves.flux_detrend import flux_space_detrend_result
from twirl.vetting import ssl_full_pool_native as native
from twirl.vetting.harmonic_inputs import read_native_light_curve_from_h5
from twirl.vetting.ssl_full_pool import (
    FULL_POOL_CONTRACT_VERSION,
    FULL_POOL_SUMMARY_SCHEMA_VERSION,
    POOL_COLUMNS,
)
from twirl.vetting.ssl_full_pool_eligibility import (
    ArtifactBinding,
    EligibilityAuthority,
    NATIVE_MODEL_ELIGIBILITY_CONTRACT_VERSION,
    NATIVE_MODEL_ELIGIBILITY_RELEASE_BINDING,
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
    include_adp: bool = True,
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
        if include_adp:
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
    compact_include_adp: bool = True,
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
        include_adp=compact_include_adp,
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


def _add_same_sector_companion(
    inputs: dict[str, Path | int | np.ndarray],
    *,
    tic: int,
    include_adp: bool = True,
) -> int:
    """Add a second observation while keeping the frozen fixture self-consistent."""

    sector = int(inputs["sector"])
    compact = Path(inputs["compact"])
    cadenceno = np.asarray(inputs["cadenceno"])
    orbitid = np.asarray(inputs["orbitid"])
    n = len(cadenceno)
    with h5py.File(compact, "a") as h5:
        group = h5.create_group(f"targets/{tic:016d}")
        group.attrs["sector"] = sector
        group.attrs["tic"] = tic
        group.attrs["camera"] = 1
        group.attrs["ccd"] = 1
        group.create_dataset("time", data=3000.0 + np.arange(n) / 720.0)
        group.create_dataset("cadenceno", data=cadenceno)
        group.create_dataset("orbitid", data=orbitid)
        group.create_dataset("quality", data=np.zeros(n, dtype=np.int32))
        if include_adp:
            group.create_dataset(
                "DET_FLUX_ADP_SML", data=np.full(n, -999_999.0)
            )
            group.create_dataset(
                "DET_FLUX_ADP", data=np.full(n, 999_999.0)
            )

    pool_path = Path(inputs["pool"])
    pool = read_table(pool_path)
    template = pool.loc[pool["sector"].eq(sector)].iloc[0].copy()
    template["observation_id"] = f"s{sector:04d}-tic{tic:016d}"
    template["tic"] = tic
    template["compact_group_path"] = f"targets/{tic:016d}"
    template["source_fits"] = f"/readonly/s{sector}/{tic}.fits"
    pool = pd.concat([pool, template.to_frame().T], ignore_index=True)
    compact_hash = _sha256(compact)
    selected = pool["sector"].astype(int).eq(sector)
    pool.loc[selected, "compact_h5_sha256"] = compact_hash
    pool.loc[selected, "compact_h5_size_bytes"] = compact.stat().st_size
    pool = pool[list(POOL_COLUMNS)].sort_values(
        ["sector", "tic"], kind="stable"
    ).reset_index(drop=True)
    pool.to_csv(pool_path, index=False, lineterminator="\n")

    allowlist = Path(inputs["allowlist"])
    selected_tics = (
        pool.loc[pool["sector"].astype(int).eq(sector), "tic"]
        .astype(int)
        .sort_values()
        .tolist()
    )
    allowlist.write_text(
        "tic\n" + "".join(f"{value}\n" for value in selected_tics),
        encoding="ascii",
    )

    summary_path = Path(inputs["summary"])
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["counts"]["retained"] = {
        "n_observations": int(len(pool)),
        "n_unique_tics": int(pool["tic"].astype(int).nunique()),
        "n_multisector_tics": int(
            pool.assign(tic=pool["tic"].astype(int))
            .groupby("tic")["sector"]
            .nunique()
            .gt(1)
            .sum()
        ),
    }
    summary["identity_hashes"]["retained_observations_sha256"] = (
        _identity_sha256(pool)
    )
    for item in summary["inputs"]["compact_exports"]:
        if int(item["sector"]) == sector:
            item["compact_h5"] = {
                "path": str(compact),
                "size_bytes": compact.stat().st_size,
                "sha256": compact_hash,
            }
    summary["outputs"]["csv"] = {
        "path": str(pool_path),
        "size_bytes": pool_path.stat().st_size,
        "sha256": _sha256(pool_path),
        "n_rows": int(len(pool)),
    }
    summary["outputs"]["sector_allowlists"][str(sector)] = {
        "path": str(allowlist),
        "size_bytes": allowlist.stat().st_size,
        "sha256": _sha256(allowlist),
        "n_tics": len(selected_tics),
        "tic_inventory_sha256": _inventory_sha256(selected_tics),
    }
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return tic


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


def _independent_effective_quality_adp(
    *,
    time: np.ndarray,
    raw_flux: np.ndarray,
    raw_error: np.ndarray,
    quality: np.ndarray,
) -> np.ndarray:
    result = flux_space_detrend_result(
        np.asarray(time, dtype=float),
        np.asarray(raw_flux, dtype=float),
        quality=np.asarray(quality),
        flux_err=np.asarray(raw_error, dtype=float),
        cfg=adp03q_config(),
    )
    detrended = np.asarray(result.det_flux, dtype=float)
    good = (np.asarray(quality) == 0) & np.isfinite(detrended)
    center = float(np.median(detrended[good]))
    return detrended - center + 1.0


def _export_raw(
    tmp_path: Path,
    inputs: dict[str, Path | int | np.ndarray],
    monkeypatch: pytest.MonkeyPatch,
) -> Path:
    raw_root = tmp_path / "raw"
    sector = int(inputs["sector"])
    orbit_one = 119 + 2 * (sector - 56)
    orbits = (orbit_one, orbit_one + 1)
    pool = read_table(Path(inputs["pool"]))
    sector_tics = (
        pool.loc[pool["sector"].astype(int).eq(sector), "tic"]
        .astype(int)
        .sort_values()
        .tolist()
    )
    for tic in sector_tics:
        for orbit in orbits:
            path = (
                raw_root
                / f"orbit-{orbit}"
                / "ffi/cam1/ccd1/LC"
                / f"{tic}.h5"
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
    assert result["n_shard_observations"] == len(sector_tics)
    return output


def _quality_reference(
    tmp_path: Path,
    *,
    sector: int = 62,
    cadenceno: np.ndarray,
    reference_orbitid: np.ndarray,
    external_quality: np.ndarray | None = None,
) -> ExternalQualityReference:
    table = tmp_path / "cadence_reference.csv"
    manifest = tmp_path / "cadence_reference.json"
    table.write_text("test\n", encoding="ascii")
    manifest.write_text("{}\n", encoding="ascii")
    if external_quality is None:
        external_quality = np.zeros(len(cadenceno), dtype=np.int64)
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
                np.asarray(external_quality, dtype=np.int64),
            )
        },
        authority_exclusions={},
        authority_exclusions_sha256="c" * 64,
    )


def _patch_all_eligible(
    tmp_path: Path,
    inputs: dict[str, Path | int | np.ndarray],
    monkeypatch: pytest.MonkeyPatch,
    *,
    excluded_keys: frozenset[tuple[int, int]] = frozenset(),
) -> tuple[Path, Path]:
    pool = read_table(Path(inputs["pool"]))
    full_keys = frozenset(
        zip(pool["sector"].astype(int), pool["tic"].astype(int))
    )
    if not excluded_keys.issubset(full_keys):
        raise ValueError("test exclusions lie outside fixture full pool")
    eligible_keys = full_keys - excluded_keys
    excluded_mask = [
        (int(row.sector), int(row.tic)) in excluded_keys
        for row in pool.itertuples(index=False)
    ]
    excluded_rows = pool.loc[excluded_mask].reset_index(drop=True)
    eligible_rows = pool.loc[
        np.logical_not(np.asarray(excluded_mask, dtype=bool))
    ].reset_index(drop=True)
    exclusions_path = tmp_path / "eligibility_exclusions.csv"
    summary_path = tmp_path / "eligibility.summary.json"
    exclusions_path.write_text("test\n", encoding="ascii")
    summary_path.write_text("{}\n", encoding="ascii")
    bindings = {
        "exclusions": ArtifactBinding(
            path=exclusions_path.resolve(),
            size_bytes=exclusions_path.stat().st_size,
            sha256=_sha256(exclusions_path),
        ),
        "eligibility_summary": ArtifactBinding(
            path=summary_path.resolve(),
            size_bytes=summary_path.stat().st_size,
            sha256=_sha256(summary_path),
        ),
    }
    authority = EligibilityAuthority(
        contract_version=NATIVE_MODEL_ELIGIBILITY_CONTRACT_VERSION,
        release_binding=NATIVE_MODEL_ELIGIBILITY_RELEASE_BINDING,
        anchor_aperture="DET_FLUX_ADP_SML",
        full_keys=full_keys,
        eligible_keys=eligible_keys,
        excluded_keys=excluded_keys,
        n_full=len(full_keys),
        n_eligible=len(eligible_keys),
        n_excluded=len(excluded_keys),
        full_observation_identity_sha256=_identity_sha256(pool),
        eligible_observation_identity_sha256=_identity_sha256(eligible_rows),
        excluded_observation_identity_sha256=_identity_sha256(excluded_rows),
        bindings=bindings,
        exclusions=excluded_rows[["sector", "tic"]].copy(),
        summary={},
    )
    monkeypatch.setattr(
        native,
        "load_native_model_eligibility",
        lambda *args, **kwargs: authority,
    )
    release_files = []
    for name in (
        "raw_source.summary.json",
        "raw_export.complete.json",
        "raw_transfer_validation.json",
    ):
        path = tmp_path / name
        path.write_text("{}\n", encoding="ascii")
        release_files.append(native._bind_file(path))
    monkeypatch.setattr(
        native,
        "PRODUCTION_RAW_EXPORT_COMPLETE_SHA256",
        release_files[1].sha256,
    )
    transfer_hashes = dict(native.PRODUCTION_RAW_TRANSFER_SHA256_BY_SECTOR)
    transfer_hashes[int(inputs["sector"])] = release_files[2].sha256
    monkeypatch.setattr(
        native,
        "PRODUCTION_RAW_TRANSFER_SHA256_BY_SECTOR",
        transfer_hashes,
    )

    def fake_raw_release(**kwargs: object) -> native.RawSourceReleaseAuthority:
        return native.RawSourceReleaseAuthority(
            raw_source=native._bind_file(Path(kwargs["raw_source_h5"])),
            raw_source_summary=release_files[0],
            export_complete=release_files[1],
            transfer_validation=release_files[2],
            code_revision=native.PRODUCTION_RAW_CODE_REVISION,
        )

    monkeypatch.setattr(
        native,
        "load_production_raw_source_release",
        fake_raw_release,
    )
    return exclusions_path, summary_path


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


def test_effective_quality_adp_fit_ignores_external_only_bad_cadence() -> None:
    n = 240
    time = 2459000.0 + np.arange(n) * (200.0 / 86400.0)
    clean = 1000.0 + 0.2 * np.sin(np.arange(n) / 15.0)
    contaminated = clean.copy()
    bad_index = 117
    contaminated[bad_index] = 1.0e18
    error = np.full(n, 0.1)
    final_effective_quality = np.zeros(n, dtype=np.int32)
    final_effective_quality[bad_index] = 1

    observed, observed_diagnostics = native._effective_quality_adp03q(
        time=time,
        raw_flux=contaminated,
        raw_error=error,
        quality=final_effective_quality,
    )
    reference, reference_diagnostics = native._effective_quality_adp03q(
        time=time,
        raw_flux=clean,
        raw_error=error,
        quality=final_effective_quality,
    )

    good = final_effective_quality == 0
    assert np.allclose(observed[good], reference[good], rtol=0.0, atol=1e-12)
    assert observed_diagnostics["fit_count"] == n - 1
    assert observed_diagnostics["fit_count"] == reference_diagnostics["fit_count"]


@pytest.mark.parametrize("n_periods", [1, 2])
def test_native_v3_builder_rejects_unlocked_periodogram_resolution(
    tmp_path: Path,
    n_periods: int,
) -> None:
    with pytest.raises(ValueError, match="periodogram resolution"):
        native.build_full_pool_native_shard(
            sector=61,
            pool_path=tmp_path / "pool.csv",
            pool_summary_path=tmp_path / "pool.json",
            eligibility_exclusions_path=tmp_path / "excluded.csv",
            eligibility_summary_path=tmp_path / "eligibility.json",
            allowlist_path=tmp_path / "allowlist.csv",
            raw_source_h5=tmp_path / "raw.h5",
            compact_adp_h5=tmp_path / "compact.h5",
            cadence_reference_table=tmp_path / "cadences.csv",
            cadence_reference_manifest=tmp_path / "cadences.json",
            out_h5=tmp_path / "native.h5",
            expected_git_sha="d" * 40,
            n_periods=n_periods,
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


def test_production_raw_release_uses_exact_csv_pool_authority(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Guard the CSV/Parquet split that raw-v1 provenance exposed."""

    sector = 62
    tic = 420_062
    source_rows = pd.DataFrame([{"sector": sector, "tic": tic}])
    identity = _identity_sha256(source_rows)
    pool_path = tmp_path / "teacher_ssl_full_pool_observations.csv"
    pool_path.write_text("sector,tic\n62,420062\n", encoding="ascii")
    pool_summary_path = tmp_path / "pool.summary.json"
    pool_summary_path.write_text('{"passed":true}\n', encoding="ascii")
    allowlist_path = tmp_path / "s62_tics.csv"
    allowlist_path.write_text("tic\n420062\n", encoding="ascii")
    raw_source_path = tmp_path / "raw_source_0.h5"
    raw_source_path.write_bytes(b"immutable raw-v1 fixture")
    raw_binding = native._bind_file(raw_source_path)

    pool_sha = _sha256(pool_path)
    pool_summary_sha = _sha256(pool_summary_path)
    allowlist_sha = _sha256(allowlist_path)
    monkeypatch.setattr(native, "FULL_POOL_NATIVE_SECTORS", (sector,))
    monkeypatch.setattr(native, "FULL_POOL_NATIVE_SHARDS_PER_SECTOR", 1)
    monkeypatch.setattr(
        native,
        "PRODUCTION_FROZEN_POOL_CSV_SHA256",
        pool_sha,
    )
    monkeypatch.setattr(
        native,
        "PRODUCTION_FULL_POOL_SUMMARY_SHA256",
        pool_summary_sha,
    )
    monkeypatch.setattr(
        native,
        "PRODUCTION_RAW_SECTOR_OBSERVATIONS",
        {sector: 1},
    )
    monkeypatch.setattr(
        native,
        "PRODUCTION_RAW_SECTOR_IDENTITY_SHA256",
        {sector: identity},
    )
    monkeypatch.setattr(
        native,
        "PRODUCTION_RAW_SECTOR_ALLOWLIST_SHA256",
        {sector: allowlist_sha},
    )

    raw_summary_path = tmp_path / "raw_source_0.summary.json"
    raw_summary_path.write_text(
        json.dumps(
            {
                "schema_version": (
                    native.FULL_POOL_RAW_SHARD_SUMMARY_SCHEMA_VERSION
                ),
                "stage": "pdo_raw_source_export",
                "sector": sector,
                "shard_index": 0,
                "n_shards": 1,
                "n_sector_observations": 1,
                "n_shard_observations": 1,
                "sector_observation_identity_sha256": identity,
                "shard_observation_identity_sha256": identity,
                "full_pool_sha256": pool_sha,
                "full_pool_summary_sha256": pool_summary_sha,
                "sector_allowlist_sha256": allowlist_sha,
                "out_h5_size_bytes": raw_binding.size_bytes,
                "out_h5_sha256": raw_binding.sha256,
                "real_only": True,
            },
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    transfer_path = tmp_path / "transfer_validation.json"
    transfer_path.write_text(
        json.dumps(
            {
                "passed": True,
                "schema_version": native.FULL_POOL_RAW_TRANSFER_SCHEMA_VERSION,
                "sector": sector,
                "n_observations": 1,
                "n_shards": 1,
                "full_pool_sha256": pool_sha,
                "full_pool_summary_sha256": pool_summary_sha,
                "sector_allowlist_sha256": allowlist_sha,
                "sector_observation_identity_sha256": identity,
                "shards": [
                    {
                        "shard_index": 0,
                        "sha256": raw_binding.sha256,
                        "size_bytes": raw_binding.size_bytes,
                        "n_observations": 1,
                        "shard_observation_identity_sha256": identity,
                    }
                ],
            },
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    transfer_hashes = dict(native.PRODUCTION_RAW_TRANSFER_SHA256_BY_SECTOR)
    transfer_hashes[sector] = _sha256(transfer_path)
    monkeypatch.setattr(
        native,
        "PRODUCTION_RAW_TRANSFER_SHA256_BY_SECTOR",
        transfer_hashes,
    )
    aggregate_path = tmp_path / "raw_export.complete.json"
    aggregate_path.write_text(
        json.dumps(
            {
                "passed": True,
                "schema_version": (
                    native.FULL_POOL_RAW_EXPORT_CONTROLLER_SCHEMA_VERSION
                ),
                "code_release": native.PRODUCTION_RAW_CODE_REVISION,
                "launcher_sha256": (
                    native.PRODUCTION_RAW_EXPORT_LAUNCHER_SHA256
                ),
                "n_observations": 1,
                "n_sectors": 1,
                "n_shards": 1,
                "n_shards_per_sector": 1,
                "full_pool_sha256": [pool_sha],
                "full_pool_summary_sha256": [pool_summary_sha],
                "sector_observation_identity_sha256": {
                    str(sector): identity
                },
            },
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(
        native,
        "PRODUCTION_RAW_EXPORT_COMPLETE_SHA256",
        _sha256(aggregate_path),
    )

    authority = native.SectorPoolAuthority(
        sector=sector,
        rows=source_rows,
        pool=native._bind_file(pool_path),
        pool_summary=native._bind_file(pool_summary_path),
        allowlist=native._bind_file(allowlist_path),
        compact_h5_sha256="c" * 64,
        compact_h5_size_bytes=1,
        observation_identity_sha256=identity,
    )
    release = native.load_production_raw_source_release(
        authority=authority,
        source_rows=source_rows,
        sector=sector,
        shard_index=0,
        n_shards=1,
        raw_source_h5=raw_source_path,
        raw_source_summary_path=raw_summary_path,
        raw_export_complete_path=aggregate_path,
        raw_transfer_validation_path=transfer_path,
    )
    assert release.raw_source.sha256 == raw_binding.sha256

    alternate_pool = tmp_path / "teacher_ssl_full_pool_observations.parquet"
    alternate_pool.write_bytes(pool_path.read_bytes() + b"alternate")
    wrong_authority = native.SectorPoolAuthority(
        sector=sector,
        rows=source_rows,
        pool=native._bind_file(alternate_pool),
        pool_summary=authority.pool_summary,
        allowlist=authority.allowlist,
        compact_h5_sha256=authority.compact_h5_sha256,
        compact_h5_size_bytes=authority.compact_h5_size_bytes,
        observation_identity_sha256=identity,
    )
    with pytest.raises(
        ValueError,
        match="sector authority differs from production lock",
    ):
        native.load_production_raw_source_release(
            authority=wrong_authority,
            source_rows=source_rows,
            sector=sector,
            shard_index=0,
            n_shards=1,
            raw_source_h5=raw_source_path,
            raw_source_summary_path=raw_summary_path,
            raw_export_complete_path=aggregate_path,
            raw_transfer_validation_path=transfer_path,
        )


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
            "bls_log_period_grid": np.linspace(
                -1.0, 0.0, native.FULL_POOL_NATIVE_PERIODOGRAM_N
            ).astype(np.float32),
            "bls_power_small": np.ones(
                native.FULL_POOL_NATIVE_PERIODOGRAM_N, dtype=np.float32
            ),
            "bls_sde_small": np.ones(
                native.FULL_POOL_NATIVE_PERIODOGRAM_N, dtype=np.float32
            ),
            "bls_power_primary": np.ones(
                native.FULL_POOL_NATIVE_PERIODOGRAM_N, dtype=np.float32
            ),
            "bls_sde_primary": np.ones(
                native.FULL_POOL_NATIVE_PERIODOGRAM_N, dtype=np.float32
            ),
        },
    )
    eligibility_exclusions, eligibility_summary = _patch_all_eligible(
        tmp_path, inputs, monkeypatch
    )
    output = tmp_path / "native.h5"
    result = native.build_full_pool_native_shard(
        sector=62,
        pool_path=Path(inputs["pool"]),
        pool_summary_path=Path(inputs["summary"]),
        eligibility_exclusions_path=eligibility_exclusions,
        eligibility_summary_path=eligibility_summary,
        allowlist_path=Path(inputs["allowlist"]),
        raw_source_h5=raw_source,
        compact_adp_h5=Path(inputs["compact"]),
        cadence_reference_table=reference.table_path,
        cadence_reference_manifest=reference.manifest_path,
        out_h5=output,
        expected_git_sha="d" * 40,
        orbitid_policy="reference_by_cadence",
        n_periods=native.FULL_POOL_NATIVE_PERIODOGRAM_N,
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
        compact_include_adp=False,
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
            name: np.ones(
                native.FULL_POOL_NATIVE_PERIODOGRAM_N, dtype=np.float32
            )
            for name in native.PERIODOGRAM_DATASETS
        },
    )
    eligibility_exclusions, eligibility_summary = _patch_all_eligible(
        tmp_path, inputs, monkeypatch
    )
    output = tmp_path / "s61_native.h5"
    result = native.build_full_pool_native_shard(
        sector=61,
        pool_path=Path(inputs["pool"]),
        pool_summary_path=Path(inputs["summary"]),
        eligibility_exclusions_path=eligibility_exclusions,
        eligibility_summary_path=eligibility_summary,
        allowlist_path=Path(inputs["allowlist"]),
        raw_source_h5=raw_source,
        compact_adp_h5=Path(inputs["compact"]),
        cadence_reference_table=reference.table_path,
        cadence_reference_manifest=reference.manifest_path,
        out_h5=output,
        expected_git_sha="d" * 40,
        orbitid_policy="strict",
        n_periods=native.FULL_POOL_NATIVE_PERIODOGRAM_N,
    )
    assert result["verification"]["passed"] is True
    assert result["orbitid_reconciliation"]["n_cad_corrected"] == 0
    with h5py.File(output, "r") as h5:
        group = h5[f"targets/{int(inputs['tic']):016d}"]
        assert np.array_equal(
            group["orbitid"][:], group["orbitid_pre_reference"][:]
        )
        quality = np.asarray(group["quality"])
        expected_small = _independent_effective_quality_adp(
            time=np.asarray(group["time"]),
            raw_flux=np.asarray(group["raw_flux_small"]),
            raw_error=np.asarray(group["raw_flux_err_small"]),
            quality=quality,
        )
        expected_primary = _independent_effective_quality_adp(
            time=np.asarray(group["time"]),
            raw_flux=np.asarray(group["raw_flux_primary"]),
            raw_error=np.asarray(group["raw_flux_err_primary"]),
            quality=quality,
        )
        assert np.allclose(group["det_flux_adp_sml"][:], expected_small)
        assert np.allclose(group["det_flux_adp"][:], expected_primary)
        assert group.attrs["raw_compact_internal_quality_agreement"] == 1
        assert group.attrs["compact_adp_photometry_reused"] == 0
        assert group.attrs["compact_adp_flux_reused"] == 0
        assert group.attrs["detrend_config_sha256"] == (
            native.FULL_POOL_NATIVE_DETREND_CONFIG_SHA256
        )


def test_native_v3_builder_reconstructs_adp_from_raw_with_external_bad_cadence(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    n = 240
    cadenceno = np.arange(10_000, 10_000 + n, dtype=np.int64)
    orbitid = np.concatenate(
        (
            np.full(n // 2, 129, dtype=np.int32),
            np.full(n - n // 2, 130, dtype=np.int32),
        )
    )
    inputs = _frozen_pool_fixture(
        tmp_path,
        sector=61,
        tic=420_061,
        cadenceno=cadenceno,
        orbitid=orbitid,
    )
    tic = int(inputs["tic"])
    bad_index = 117
    clean = _raw_payload(cadenceno, orbitid)
    with h5py.File(Path(inputs["compact"]), "a") as h5:
        group = h5[f"targets/{tic:016d}"]
        group["DET_FLUX_ADP_SML"][:] = -999_999.0
        group["DET_FLUX_ADP"][:] = 999_999.0
    compact = Path(inputs["compact"])
    pool_path = Path(inputs["pool"])
    pool = read_table(pool_path)
    selected = pool["sector"].astype(int).eq(61)
    pool.loc[selected, "compact_h5_sha256"] = _sha256(compact)
    pool.loc[selected, "compact_h5_size_bytes"] = compact.stat().st_size
    pool.to_csv(pool_path, index=False, lineterminator="\n")
    summary_path = Path(inputs["summary"])
    pool_summary = json.loads(summary_path.read_text(encoding="utf-8"))
    for item in pool_summary["inputs"]["compact_exports"]:
        if int(item["sector"]) == 61:
            item["compact_h5"] = {
                "path": str(compact),
                "size_bytes": compact.stat().st_size,
                "sha256": _sha256(compact),
            }
    pool_summary["outputs"]["csv"].update(
        {
            "size_bytes": pool_path.stat().st_size,
            "sha256": _sha256(pool_path),
        }
    )
    summary_path.write_text(
        json.dumps(pool_summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    raw_source = _export_raw(tmp_path, inputs, monkeypatch)
    with h5py.File(raw_source, "a") as h5:
        group = h5[f"targets/{tic:016d}"]
        group["raw_flux_small"][bad_index] = 1.0e18
        group["raw_flux_primary"][bad_index] = -1.0e18

    external_quality = np.zeros(n, dtype=np.int64)
    external_quality[bad_index] = 1
    reference = _quality_reference(
        tmp_path,
        sector=61,
        cadenceno=cadenceno,
        reference_orbitid=orbitid,
        external_quality=external_quality,
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
            name: np.ones(
                native.FULL_POOL_NATIVE_PERIODOGRAM_N, dtype=np.float32
            )
            for name in native.PERIODOGRAM_DATASETS
        },
    )
    eligibility_exclusions, eligibility_summary = _patch_all_eligible(
        tmp_path, inputs, monkeypatch
    )
    output = tmp_path / "effective_quality_native.h5"
    result = native.build_full_pool_native_shard(
        sector=61,
        pool_path=Path(inputs["pool"]),
        pool_summary_path=Path(inputs["summary"]),
        eligibility_exclusions_path=eligibility_exclusions,
        eligibility_summary_path=eligibility_summary,
        allowlist_path=Path(inputs["allowlist"]),
        raw_source_h5=raw_source,
        compact_adp_h5=Path(inputs["compact"]),
        cadence_reference_table=reference.table_path,
        cadence_reference_manifest=reference.manifest_path,
        out_h5=output,
        expected_git_sha="d" * 40,
        orbitid_policy="strict",
        n_periods=native.FULL_POOL_NATIVE_PERIODOGRAM_N,
    )

    assert result["verification"]["passed"] is True
    with h5py.File(output, "r") as h5:
        group = h5[f"targets/{tic:016d}"]
        quality = np.asarray(group["quality"])
        expected_small = _independent_effective_quality_adp(
            time=np.asarray(group["time"]),
            raw_flux=np.asarray(clean["raw_flux_small"]),
            raw_error=np.asarray(clean["raw_flux_err_small"]),
            quality=quality,
        )
        expected_primary = _independent_effective_quality_adp(
            time=np.asarray(group["time"]),
            raw_flux=np.asarray(clean["raw_flux_primary"]),
            raw_error=np.asarray(clean["raw_flux_err_primary"]),
            quality=quality,
        )
        good = quality == 0
        assert np.allclose(
            group["det_flux_adp_sml"][:][good],
            expected_small[good],
            rtol=0.0,
            atol=1.0e-12,
        )
        assert np.allclose(
            group["det_flux_adp"][:][good],
            expected_primary[good],
            rtol=0.0,
            atol=1.0e-12,
        )
        assert group.attrs["compact_adp_photometry_reused"] == 0
        assert group.attrs["compact_adp_flux_reused"] == 0


def test_excluded_observation_is_structurally_validated_but_not_published(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    inputs = _frozen_pool_fixture(
        tmp_path,
        sector=61,
        tic=420_061,
        cadenceno=np.asarray([100, 101, 102, 103]),
        orbitid=np.asarray([129, 129, 130, 130]),
    )
    excluded_tic = _add_same_sector_companion(
        inputs,
        tic=520_061,
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
            name: np.ones(
                native.FULL_POOL_NATIVE_PERIODOGRAM_N, dtype=np.float32
            )
            for name in native.PERIODOGRAM_DATASETS
        },
    )
    eligibility_exclusions, eligibility_summary = _patch_all_eligible(
        tmp_path,
        inputs,
        monkeypatch,
        excluded_keys=frozenset({(61, excluded_tic)}),
    )
    output = tmp_path / "excluded_canary.h5"
    result = native.build_full_pool_native_shard(
        sector=61,
        pool_path=Path(inputs["pool"]),
        pool_summary_path=Path(inputs["summary"]),
        eligibility_exclusions_path=eligibility_exclusions,
        eligibility_summary_path=eligibility_summary,
        allowlist_path=Path(inputs["allowlist"]),
        raw_source_h5=raw_source,
        compact_adp_h5=Path(inputs["compact"]),
        cadence_reference_table=reference.table_path,
        cadence_reference_manifest=reference.manifest_path,
        out_h5=output,
        expected_git_sha="d" * 40,
        orbitid_policy="strict",
        n_periods=native.FULL_POOL_NATIVE_PERIODOGRAM_N,
    )

    assert result["verification"]["passed"] is True
    assert result["n_source_shard_observations"] == 2
    assert result["n_shard_observations"] == 1
    assert result["n_shard_excluded_observations"] == 1
    assert result["n_shard_excluded_structurally_validated"] == 1
    with h5py.File(output, "r") as h5:
        assert sorted(h5["targets"]) == [f"{int(inputs['tic']):016d}"]
        assert h5.attrs["n_source_shard_observations"] == 2
        assert h5.attrs["n_shard_observations"] == 1
        assert h5.attrs["n_shard_excluded_observations"] == 1
        assert h5.attrs["n_shard_excluded_structurally_validated"] == 1
        assert f"{excluded_tic:016d}" not in h5["targets"]


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
    eligibility_exclusions, eligibility_summary = _patch_all_eligible(
        tmp_path, inputs, monkeypatch
    )
    with pytest.raises(RuntimeError, match="native export failed"):
        native.build_full_pool_native_shard(
            sector=62,
            pool_path=Path(inputs["pool"]),
            pool_summary_path=Path(inputs["summary"]),
            eligibility_exclusions_path=eligibility_exclusions,
            eligibility_summary_path=eligibility_summary,
            allowlist_path=Path(inputs["allowlist"]),
            raw_source_h5=raw_source,
            compact_adp_h5=Path(inputs["compact"]),
            cadence_reference_table=reference.table_path,
            cadence_reference_manifest=reference.manifest_path,
            out_h5=tmp_path / "bad_native.h5",
            expected_git_sha="d" * 40,
            orbitid_policy="reference_by_cadence",
            n_periods=native.FULL_POOL_NATIVE_PERIODOGRAM_N,
        )
    failures = pd.read_csv((tmp_path / "bad_native.h5").with_suffix(".failures.csv"))
    assert "raw-source and compact orbitid arrays disagree" in failures.loc[
        0, "error"
    ]


def test_native_build_rejects_raw_compact_quality_mismatch_without_publication(
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
    with h5py.File(raw_source, "a") as h5:
        h5[f"targets/{int(inputs['tic']):016d}/quality"][1] = 1
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
    eligibility_exclusions, eligibility_summary = _patch_all_eligible(
        tmp_path, inputs, monkeypatch
    )
    output = tmp_path / "quality_mismatch.h5"
    with pytest.raises(RuntimeError, match="native export failed"):
        native.build_full_pool_native_shard(
            sector=61,
            pool_path=Path(inputs["pool"]),
            pool_summary_path=Path(inputs["summary"]),
            eligibility_exclusions_path=eligibility_exclusions,
            eligibility_summary_path=eligibility_summary,
            allowlist_path=Path(inputs["allowlist"]),
            raw_source_h5=raw_source,
            compact_adp_h5=Path(inputs["compact"]),
            cadence_reference_table=reference.table_path,
            cadence_reference_manifest=reference.manifest_path,
            out_h5=output,
            expected_git_sha="d" * 40,
            orbitid_policy="strict",
            n_periods=native.FULL_POOL_NATIVE_PERIODOGRAM_N,
        )
    assert not output.exists()
    assert not list(tmp_path.glob(".quality_mismatch.h5.*.tmp"))
    failures = pd.read_csv(output.with_suffix(".failures.csv"))
    assert "internal quality arrays disagree" in failures.loc[0, "error"]


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
            name: np.ones(
                native.FULL_POOL_NATIVE_PERIODOGRAM_N, dtype=np.float32
            )
            for name in native.PERIODOGRAM_DATASETS
        },
    )
    eligibility_exclusions, eligibility_summary = _patch_all_eligible(
        tmp_path, inputs, monkeypatch
    )
    output = tmp_path / "native.h5"
    built = native.build_full_pool_native_shard(
        sector=62,
        pool_path=Path(inputs["pool"]),
        pool_summary_path=Path(inputs["summary"]),
        eligibility_exclusions_path=eligibility_exclusions,
        eligibility_summary_path=eligibility_summary,
        allowlist_path=Path(inputs["allowlist"]),
        raw_source_h5=raw_source,
        compact_adp_h5=Path(inputs["compact"]),
        cadence_reference_table=reference.table_path,
        cadence_reference_manifest=reference.manifest_path,
        out_h5=output,
        expected_git_sha="d" * 40,
        orbitid_policy="reference_by_cadence",
        n_periods=native.FULL_POOL_NATIVE_PERIODOGRAM_N,
    )
    output_summary = output.with_suffix(".summary.json")
    output_summary.write_text(
        json.dumps(built, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="inventory count"):
        native.write_full_pool_native_registry(
            pool_path=Path(inputs["pool"]),
            pool_summary_path=Path(inputs["summary"]),
            eligibility_exclusions_path=eligibility_exclusions,
            eligibility_summary_path=eligibility_summary,
            native_shard_paths=[output],
            native_shard_summary_paths=[output_summary],
            source_path=tmp_path / "source.csv",
            registry_path=tmp_path / "registry.csv",
            summary_path=tmp_path / "registry.summary.json",
            release_summary_path=tmp_path / "registry.release.json",
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
            eligibility_exclusions_path=eligibility_exclusions,
            eligibility_summary_path=eligibility_summary,
            native_shard_paths=[output],
            native_shard_summary_paths=[output_summary],
            source_path=tmp_path / "source2.csv",
            registry_path=tmp_path / "registry2.csv",
            summary_path=tmp_path / "registry2.summary.json",
            release_summary_path=tmp_path / "registry2.release.json",
        )

    # Minimal contract-carrying shard fixtures are sufficient to exercise the
    # exact observation-keyed registry assembly itself.
    fake_shards: list[Path] = []
    fake_summaries: list[Path] = []
    pool_sha = _sha256(Path(inputs["pool"]))
    summary_sha = _sha256(Path(inputs["summary"]))
    monkeypatch.setattr(native, "FULL_POOL_NATIVE_SHARDS_PER_SECTOR", 1)
    for row in pool.itertuples(index=False):
        shard = tmp_path / f"s{int(row.sector)}_native.h5"
        with h5py.File(shard, "w") as h5:
            h5.attrs["contract_version"] = (
                native.FULL_POOL_NATIVE_CONTRACT_VERSION
            )
            h5.attrs["release_binding"] = (
                native.FULL_POOL_NATIVE_RELEASE_BINDING
            )
            h5.attrs["builder_contract_version"] = (
                native.FULL_POOL_NATIVE_BUILDER_CONTRACT_VERSION
            )
            h5.attrs["builder_code_sha256"] = _sha256(
                Path(native.__file__)
            )
            h5.attrs["expected_git_sha"] = "d" * 40
            h5.attrs["detrend_contract_version"] = (
                native.FULL_POOL_NATIVE_DETREND_CONTRACT_VERSION
            )
            h5.attrs["detrend_config_json"] = (
                native.FULL_POOL_NATIVE_DETREND_CONFIG_JSON
            )
            h5.attrs["detrend_config_sha256"] = (
                native.FULL_POOL_NATIVE_DETREND_CONFIG_SHA256
            )
            h5.attrs["detrend_quality_source"] = (
                "final_effective_quality"
            )
            h5.attrs["raw_photometry_only"] = 1
            h5.attrs["compact_adp_photometry_reused"] = 0
            h5.attrs["compact_adp_flux_reused"] = 0
            h5.attrs["periodogram_n"] = (
                native.FULL_POOL_NATIVE_PERIODOGRAM_N
            )
            h5.attrs["full_pool_sha256"] = pool_sha
            h5.attrs["full_pool_summary_sha256"] = summary_sha
            h5.attrs["eligibility_exclusions_sha256"] = _sha256(
                eligibility_exclusions
            )
            h5.attrs["eligibility_summary_sha256"] = _sha256(
                eligibility_summary
            )
            h5.attrs["raw_source_release_code_revision"] = (
                native.PRODUCTION_RAW_CODE_REVISION
            )
            h5.attrs["raw_source_summary_sha256"] = "d" * 64
            h5.attrs["raw_export_complete_sha256"] = (
                native.PRODUCTION_RAW_EXPORT_COMPLETE_SHA256
            )
            h5.attrs["raw_transfer_validation_sha256"] = (
                native.PRODUCTION_RAW_TRANSFER_SHA256_BY_SECTOR[
                    int(row.sector)
                ]
            )
            h5.attrs["sector"] = int(row.sector)
            group = h5.create_group(f"targets/{int(row.tic):016d}")
            group.attrs["sector"] = int(row.sector)
            group.attrs["tic"] = int(row.tic)
        fake_shards.append(shard)
        shard_summary = shard.with_suffix(".summary.json")
        one_row = pd.DataFrame(
            [{"sector": int(row.sector), "tic": int(row.tic)}]
        )
        shard_summary.write_text(
            json.dumps(
                {
                    "schema_version": (
                        native.FULL_POOL_NATIVE_SUMMARY_SCHEMA_VERSION
                    ),
                    "stage": "orcd_native_shard_build",
                    "sector": int(row.sector),
                    "shard_index": 0,
                    "n_shards": 1,
                    "n_source_shard_observations": 1,
                    "n_shard_observations": 1,
                        "n_shard_excluded_observations": 0,
                        "n_shard_excluded_structurally_validated": 0,
                    "source_shard_observation_identity_sha256": (
                        _identity_sha256(one_row)
                    ),
                    "shard_observation_identity_sha256": (
                        _identity_sha256(one_row)
                    ),
                    "shard_excluded_observation_identity_sha256": (
                        _identity_sha256(one_row.iloc[0:0])
                    ),
                    "full_pool_sha256": pool_sha,
                    "full_pool_summary_sha256": summary_sha,
                    "eligibility_exclusions_sha256": _sha256(
                        eligibility_exclusions
                    ),
                    "eligibility_summary_sha256": _sha256(
                        eligibility_summary
                    ),
                    "native_model_full_identity_sha256": (
                        _identity_sha256(pool)
                    ),
                    "native_model_eligible_identity_sha256": (
                        _identity_sha256(pool)
                    ),
                    "native_model_excluded_identity_sha256": (
                        _identity_sha256(pool.iloc[0:0])
                    ),
                    "raw_source_h5_sha256": "c" * 64,
                    "raw_source_summary_sha256": "d" * 64,
                    "raw_export_complete_sha256": (
                        native.PRODUCTION_RAW_EXPORT_COMPLETE_SHA256
                    ),
                    "raw_transfer_validation_sha256": (
                        native.PRODUCTION_RAW_TRANSFER_SHA256_BY_SECTOR[
                            int(row.sector)
                        ]
                    ),
                    "raw_source_release_code_revision": (
                        native.PRODUCTION_RAW_CODE_REVISION
                    ),
                    "out_h5": str(shard),
                    "out_h5_sha256": _sha256(shard),
                    "out_h5_size_bytes": shard.stat().st_size,
                    "native_contract_version": (
                        native.FULL_POOL_NATIVE_CONTRACT_VERSION
                    ),
                    "expected_git_sha": "d" * 40,
                    "builder_contract_version": (
                        native.FULL_POOL_NATIVE_BUILDER_CONTRACT_VERSION
                    ),
                    "builder_code_sha256": _sha256(
                        Path(native.__file__)
                    ),
                    "detrend_contract_version": (
                        native.FULL_POOL_NATIVE_DETREND_CONTRACT_VERSION
                    ),
                    "detrend_config_json": (
                        native.FULL_POOL_NATIVE_DETREND_CONFIG_JSON
                    ),
                    "detrend_config_sha256": (
                        native.FULL_POOL_NATIVE_DETREND_CONFIG_SHA256
                    ),
                    "detrend_quality_source": (
                        "final_effective_quality"
                    ),
                    "raw_photometry_only": True,
                    "compact_adp_photometry_reused": False,
                    "compact_adp_flux_reused": False,
                    "periodogram_n": (
                        native.FULL_POOL_NATIVE_PERIODOGRAM_N
                    ),
                    "real_only": True,
                    "verification": {"passed": True},
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        fake_summaries.append(shard_summary)

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
        eligibility_exclusions_path=eligibility_exclusions,
        eligibility_summary_path=eligibility_summary,
        native_shard_paths=fake_shards,
        native_shard_summary_paths=fake_summaries,
        source_path=tmp_path / "source_ok.csv",
        registry_path=tmp_path / "registry_ok.csv",
        summary_path=tmp_path / "registry_ok.summary.json",
        release_summary_path=tmp_path / "registry_ok.release.json",
        expected_shards_per_sector=1,
    )
    registry = pd.read_csv(tmp_path / "registry_ok.csv")
    assert len(registry) == 7
    assert registry[["sector", "tic"]].drop_duplicates().shape[0] == 7
    assert result["n_observations"] == 7
    eligibility = native.load_native_model_eligibility(
        eligibility_exclusions,
        eligibility_summary,
    )
    loaded_registry, release = native.load_full_pool_native_registry_release(
        registry_path=tmp_path / "registry_ok.csv",
        registry_summary_path=tmp_path / "registry_ok.summary.json",
        release_summary_path=tmp_path / "registry_ok.release.json",
        eligibility=eligibility,
    )
    assert len(loaded_registry) == 7
    assert release["partition_audit"]["native_registry_equals_eligible"] is True


def test_native_v2_verifies_temporary_file_before_publication(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
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
            name: np.ones(
                native.FULL_POOL_NATIVE_PERIODOGRAM_N, dtype=np.float32
            )
            for name in native.PERIODOGRAM_DATASETS
        },
    )
    eligibility_exclusions, eligibility_summary = _patch_all_eligible(
        tmp_path, inputs, monkeypatch
    )
    monkeypatch.setattr(
        native,
        "verify_full_pool_native_shard",
        lambda *args, **kwargs: {
            "passed": False,
            "failures": ["forced prepublication verification failure"],
            "counts": {"targets": 1, "injections": 0},
        },
    )
    output = tmp_path / "must_not_publish.h5"
    with pytest.raises(RuntimeError, match="temporary.*failed verification"):
        native.build_full_pool_native_shard(
            sector=61,
            pool_path=Path(inputs["pool"]),
            pool_summary_path=Path(inputs["summary"]),
            eligibility_exclusions_path=eligibility_exclusions,
            eligibility_summary_path=eligibility_summary,
            allowlist_path=Path(inputs["allowlist"]),
            raw_source_h5=raw_source,
            compact_adp_h5=Path(inputs["compact"]),
            cadence_reference_table=reference.table_path,
            cadence_reference_manifest=reference.manifest_path,
            out_h5=output,
            expected_git_sha="d" * 40,
            orbitid_policy="strict",
            n_periods=native.FULL_POOL_NATIVE_PERIODOGRAM_N,
        )
    assert not output.exists()
    assert not list(tmp_path.glob(".must_not_publish.h5.*.tmp"))


def test_fullpool_module_keeps_external_quality_contract_constants() -> None:
    assert AUTHORITY_EXCLUSION_POLICY_CONTRACT
    assert native.FULL_POOL_NATIVE_SECTORS == tuple(range(56, 63))
