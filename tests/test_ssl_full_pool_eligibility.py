from __future__ import annotations

import json
from pathlib import Path
import shutil

import numpy as np
import pandas as pd
import pytest

import twirl.vetting.ssl_full_pool_eligibility as eligibility_module
from twirl.vetting.adp_only import ADP_ONLY_APERTURES
from twirl.vetting.ssl_full_pool import (
    FULL_POOL_CONTRACT_VERSION,
    FULL_POOL_SUMMARY_SCHEMA_VERSION,
)
from twirl.vetting.ssl_full_pool_bls import (
    GLOBAL_BLS_CONTRACT_VERSION,
    GLOBAL_BLS_SUMMARY_SCHEMA_VERSION,
)
from twirl.vetting.ssl_full_pool_eligibility import (
    ELIGIBILITY_EXCLUSION_COLUMNS,
    NATIVE_MODEL_ELIGIBILITY_CONTRACT_VERSION,
    PRODUCTION_EXCLUDED_IDENTITY_SHA256,
    derive_anchor_eligibility,
    load_native_model_eligibility,
    observation_identity_sha256,
    write_native_model_eligibility,
)
from twirl.vetting.teacher_native_registry import file_sha256


def _metadata(path: Path) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": file_sha256(path),
    }


def _write_source_authorities(
    tmp_path: Path,
) -> tuple[Path, Path, Path, Path, pd.DataFrame, pd.DataFrame]:
    pool = pd.DataFrame(
        {
            "pool_contract_version": FULL_POOL_CONTRACT_VERSION,
            "observation_id": [
                "s0056-tic0000000000000101",
                "s0056-tic0000000000000102",
                "s0057-tic0000000000000101",
            ],
            "sector": [56, 56, 57],
            "tic": [101, 102, 101],
        }
    )
    pool_path = tmp_path / "teacher_ssl_full_pool_observations.csv"
    pool.to_csv(pool_path, index=False, lineterminator="\n")
    pool_summary_path = tmp_path / "teacher_ssl_full_pool_manifest.summary.json"
    pool_identity = observation_identity_sha256(pool)
    pool_summary = {
        "passed": True,
        "schema_version": FULL_POOL_SUMMARY_SCHEMA_VERSION,
        "pool_contract_version": FULL_POOL_CONTRACT_VERSION,
        "sectors": [56, 57],
        "observation_identity_columns": ["sector", "tic"],
        "counts": {
            "retained": {
                "n_observations": 3,
                "n_unique_tics": 2,
                "n_multisector_tics": 1,
            }
        },
        "identity_hashes": {"retained_observations_sha256": pool_identity},
        "outputs": {
            "csv": {
                **_metadata(pool_path),
                "n_rows": 3,
            }
        },
    }
    pool_summary_path.write_text(
        json.dumps(pool_summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    bls = pd.DataFrame(
        [
            {
                "sector": 56,
                "tic": 101,
                "aperture": "DET_FLUX_ADP_SML",
                "peak_rank": 1,
                "status": "ok",
                "period_d": 1.2,
                "t0_bjd": 100.0,
                "duration_min": 10.0,
            },
            {
                "sector": 56,
                "tic": 101,
                "aperture": "DET_FLUX_ADP",
                "peak_rank": 1,
                "status": "ok",
                "period_d": 1.2,
                "t0_bjd": 100.0,
                "duration_min": 10.0,
            },
            {
                "sector": 56,
                "tic": 102,
                "aperture": "DET_FLUX_ADP_SML",
                "peak_rank": 0,
                "status": "too_few_cadences",
                "period_d": np.nan,
                "t0_bjd": np.nan,
                "duration_min": np.nan,
            },
            {
                "sector": 56,
                "tic": 102,
                "aperture": "DET_FLUX_ADP",
                "peak_rank": 1,
                "status": "ok",
                "period_d": 2.2,
                "t0_bjd": 101.0,
                "duration_min": 11.0,
            },
            {
                "sector": 57,
                "tic": 101,
                "aperture": "DET_FLUX_ADP_SML",
                "peak_rank": 1,
                "status": "ok",
                "period_d": 3.2,
                "t0_bjd": 102.0,
                "duration_min": 12.0,
            },
            {
                "sector": 57,
                "tic": 101,
                "aperture": "DET_FLUX_ADP",
                "peak_rank": 1,
                "status": "ok",
                "period_d": 3.2,
                "t0_bjd": 102.0,
                "duration_min": 12.0,
            },
        ]
    )
    bls_path = tmp_path / "teacher_ssl_fullpool_bls_global.parquet"
    bls.to_parquet(bls_path, index=False)
    bls_summary_path = tmp_path / "teacher_ssl_fullpool_bls_global.summary.json"
    bls_summary = {
        "passed": True,
        "schema_version": GLOBAL_BLS_SUMMARY_SCHEMA_VERSION,
        "contract_version": GLOBAL_BLS_CONTRACT_VERSION,
        "sectors": [56, 57],
        "observation_identity_columns": ["sector", "tic"],
        "frozen_pool": {
            "summary": _metadata(pool_summary_path),
            "contract_version": FULL_POOL_CONTRACT_VERSION,
            "observation_identity_sha256": pool_identity,
            "n_observations": 3,
            "n_unique_tics": 2,
        },
        "bls_contract": {"apertures": list(ADP_ONLY_APERTURES)},
        "counts": {
            "n_rows": len(bls),
            "n_observations": 3,
            "n_unique_tics": 2,
        },
        "coverage_audit": {
            "missing_frozen_observations": 0,
            "unexpected_bls_observations": 0,
            "observation_identity_sha256": pool_identity,
        },
        "sector_products": [],
        "output": {
            **_metadata(bls_path),
            "n_rows": len(bls),
            "n_observations": 3,
        },
    }
    bls_summary_path.write_text(
        json.dumps(bls_summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return (
        pool_path,
        pool_summary_path,
        bls_path,
        bls_summary_path,
        pool,
        bls,
    )


def _lock_toy_production_authority(
    monkeypatch: pytest.MonkeyPatch,
    *,
    pool_summary_path: Path,
    bls_path: Path,
    bls_summary_path: Path,
) -> None:
    bls_summary = json.loads(bls_summary_path.read_text(encoding="utf-8"))
    monkeypatch.setattr(
        eligibility_module,
        "PRODUCTION_FULL_OBSERVATIONS",
        3,
    )
    monkeypatch.setattr(
        eligibility_module,
        "PRODUCTION_ELIGIBLE_OBSERVATIONS",
        2,
    )
    monkeypatch.setattr(
        eligibility_module,
        "PRODUCTION_EXCLUDED_OBSERVATIONS",
        1,
    )
    monkeypatch.setattr(
        eligibility_module,
        "PRODUCTION_FULL_IDENTITY_SHA256",
        observation_identity_sha256({(56, 101), (56, 102), (57, 101)}),
    )
    monkeypatch.setattr(
        eligibility_module,
        "PRODUCTION_ELIGIBLE_IDENTITY_SHA256",
        observation_identity_sha256({(56, 101), (57, 101)}),
    )
    monkeypatch.setattr(
        eligibility_module,
        "PRODUCTION_EXCLUDED_IDENTITY_SHA256",
        observation_identity_sha256({(56, 102)}),
    )
    monkeypatch.setattr(
        eligibility_module,
        "PRODUCTION_FROZEN_POOL_CSV_SHA256",
        json.loads(pool_summary_path.read_text(encoding="utf-8"))["outputs"][
            "csv"
        ]["sha256"],
    )
    monkeypatch.setattr(
        eligibility_module,
        "PRODUCTION_FROZEN_POOL_SUMMARY_SHA256",
        file_sha256(pool_summary_path),
    )
    monkeypatch.setattr(
        eligibility_module,
        "PRODUCTION_GLOBAL_BLS_SHA256",
        file_sha256(bls_path),
    )
    monkeypatch.setattr(
        eligibility_module,
        "PRODUCTION_GLOBAL_BLS_SUMMARY_SHA256",
        file_sha256(bls_summary_path),
    )
    monkeypatch.setattr(
        eligibility_module,
        "PRODUCTION_EXCLUDED_BY_SECTOR",
        {56: 1},
    )
    expected_contract = bls_summary["bls_contract"]
    monkeypatch.setattr(
        eligibility_module,
        "_production_bls_contract",
        lambda: expected_contract,
    )


def test_anchor_eligibility_is_only_adp_small_rank_one_searchability(
    tmp_path: Path,
) -> None:
    *_paths, pool, bls = _write_source_authorities(tmp_path)

    decisions = derive_anchor_eligibility(bls, frozen_pool=pool)

    assert decisions[
        ["sector", "tic", "native_model_eligible", "exclusion_reason"]
    ].to_dict(orient="records") == [
        {
            "sector": 56,
            "tic": 101,
            "native_model_eligible": True,
            "exclusion_reason": "",
        },
        {
            "sector": 56,
            "tic": 102,
            "native_model_eligible": False,
            "exclusion_reason": "bls_unsearchable",
        },
        {
            "sector": 57,
            "tic": 101,
            "native_model_eligible": True,
            "exclusion_reason": "",
        },
    ]
    assert observation_identity_sha256(pool.sample(frac=1, random_state=7)) == (
        observation_identity_sha256(pool)
    )


def test_known_production_exclusion_inventory_has_locked_identity() -> None:
    keys = {
        (56, 1_201_200_937),
        (56, 1_201_317_288),
        (56, 1_973_584_484),
        (56, 2_019_898_202),
        (57, 1_400_692_573),
        (58, 1_718_236_078),
        (59, 663_454_515),
        (59, 664_713_400),
        (59, 1_201_208_167),
        (59, 1_201_208_426),
        (60, 1_271_391_891),
        (60, 1_400_673_051),
        (61, 737_545_960),
        (61, 737_624_449),
        (61, 821_867_415),
        (62, 807_785_577),
        (62, 818_817_688),
        (62, 820_202_158),
        (62, 876_333_840),
    }

    assert observation_identity_sha256(keys) == (PRODUCTION_EXCLUDED_IDENTITY_SHA256)


def test_eligibility_writer_and_loader_freeze_exact_partition(
    tmp_path: Path,
) -> None:
    pool_path, pool_summary_path, bls_path, bls_summary_path, *_ = (
        _write_source_authorities(tmp_path)
    )
    exclusions_path = tmp_path / "native_model_exclusions.csv"
    summary_path = tmp_path / "native_model_eligibility.summary.json"

    first = write_native_model_eligibility(
        pool_path=pool_path,
        pool_summary_path=pool_summary_path,
        bls_path=bls_path,
        bls_summary_path=bls_summary_path,
        exclusions_path=exclusions_path,
        summary_path=summary_path,
    )
    first_hashes = (file_sha256(exclusions_path), file_sha256(summary_path))
    second = write_native_model_eligibility(
        pool_path=pool_path,
        pool_summary_path=pool_summary_path,
        bls_path=bls_path,
        bls_summary_path=bls_summary_path,
        exclusions_path=exclusions_path,
        summary_path=summary_path,
    )

    assert first.full_keys == frozenset({(56, 101), (56, 102), (57, 101)})
    assert first.eligible_keys == frozenset({(56, 101), (57, 101)})
    assert first.excluded_keys == frozenset({(56, 102)})
    assert (first.n_full, first.n_eligible, first.n_excluded) == (3, 2, 1)
    assert first_hashes == (
        file_sha256(exclusions_path),
        file_sha256(summary_path),
    )
    assert second.excluded_keys == first.excluded_keys
    exclusions = pd.read_csv(exclusions_path)
    assert tuple(exclusions.columns) == ELIGIBILITY_EXCLUSION_COLUMNS
    assert exclusions.loc[0, "contract_version"] == (
        NATIVE_MODEL_ELIGIBILITY_CONTRACT_VERSION
    )
    assert exclusions.loc[0, "bls_status"] == "too_few_cadences"
    assert first.summary["data_usage_audit"] == {
        "labels_consumed": False,
        "injections_consumed": False,
        "raw_flux_errors_define_eligibility": False,
        "two_aperture_status_row_counts_define_eligibility": False,
    }


def test_loader_accepts_only_byte_identical_staged_authorities(
    tmp_path: Path,
) -> None:
    pool_path, pool_summary_path, bls_path, bls_summary_path, *_ = (
        _write_source_authorities(tmp_path)
    )
    exclusions_path = tmp_path / "native_model_exclusions.csv"
    summary_path = tmp_path / "native_model_eligibility.summary.json"
    write_native_model_eligibility(
        pool_path=pool_path,
        pool_summary_path=pool_summary_path,
        bls_path=bls_path,
        bls_summary_path=bls_summary_path,
        exclusions_path=exclusions_path,
        summary_path=summary_path,
    )
    staged = tmp_path / "staged"
    staged.mkdir()
    staged_paths = []
    for source in (
        pool_path,
        pool_summary_path,
        bls_path,
        bls_summary_path,
    ):
        destination = staged / source.name
        shutil.copyfile(source, destination)
        staged_paths.append(destination)

    authority = load_native_model_eligibility(
        exclusions_path,
        summary_path,
        pool_path=staged_paths[0],
        pool_summary_path=staged_paths[1],
        bls_path=staged_paths[2],
        bls_summary_path=staged_paths[3],
    )
    assert authority.n_eligible == 2

    staged_paths[0].write_text(
        staged_paths[0].read_text(encoding="utf-8") + "\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="differs from native/model"):
        load_native_model_eligibility(
            exclusions_path,
            summary_path,
            pool_path=staged_paths[0],
            pool_summary_path=staged_paths[1],
            bls_path=staged_paths[2],
            bls_summary_path=staged_paths[3],
        )


def test_shard_loader_validates_complement_without_opening_global_bls(
    tmp_path: Path,
) -> None:
    pool_path, pool_summary_path, bls_path, bls_summary_path, *_ = (
        _write_source_authorities(tmp_path)
    )
    exclusions_path = tmp_path / "native_model_exclusions.csv"
    summary_path = tmp_path / "native_model_eligibility.summary.json"
    write_native_model_eligibility(
        pool_path=pool_path,
        pool_summary_path=pool_summary_path,
        bls_path=bls_path,
        bls_summary_path=bls_summary_path,
        exclusions_path=exclusions_path,
        summary_path=summary_path,
    )
    bls_path.unlink()
    bls_summary_path.unlink()

    authority = load_native_model_eligibility(
        exclusions_path,
        summary_path,
        pool_path=pool_path,
        pool_summary_path=pool_summary_path,
        rederive_from_bls=False,
    )

    assert authority.eligible_keys == frozenset({(56, 101), (57, 101)})
    assert authority.bindings["global_bls"].sha256
    with pytest.raises(FileNotFoundError):
        load_native_model_eligibility(
            exclusions_path,
            summary_path,
            pool_path=pool_path,
            pool_summary_path=pool_summary_path,
            rederive_from_bls=True,
        )


def test_eligibility_fails_closed_on_tamper_labels_and_wrong_production_lock(
    tmp_path: Path,
) -> None:
    pool_path, pool_summary_path, bls_path, bls_summary_path, pool, bls = (
        _write_source_authorities(tmp_path)
    )
    exclusions_path = tmp_path / "native_model_exclusions.csv"
    summary_path = tmp_path / "native_model_eligibility.summary.json"
    write_native_model_eligibility(
        pool_path=pool_path,
        pool_summary_path=pool_summary_path,
        bls_path=bls_path,
        bls_summary_path=bls_summary_path,
        exclusions_path=exclusions_path,
        summary_path=summary_path,
    )
    exclusions_path.write_text(
        exclusions_path.read_text(encoding="utf-8").replace(
            "too_few_cadences",
            "ok",
        ),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="exclusions differ from summary"):
        load_native_model_eligibility(exclusions_path, summary_path)

    bls["human_label"] = "Other"
    with pytest.raises(ValueError, match="label-free"):
        derive_anchor_eligibility(bls, frozen_pool=pool)

    with pytest.raises(ValueError, match="production lock"):
        write_native_model_eligibility(
            pool_path=pool_path,
            pool_summary_path=pool_summary_path,
            bls_path=bls_path,
            bls_summary_path=bls_summary_path,
            exclusions_path=tmp_path / "locked_exclusions.csv",
            summary_path=tmp_path / "locked_summary.json",
            production_lock=True,
        )


def test_declared_production_loader_rechecks_exact_source_hashes(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    pool_path, pool_summary_path, bls_path, bls_summary_path, *_ = (
        _write_source_authorities(tmp_path)
    )
    _lock_toy_production_authority(
        monkeypatch,
        pool_summary_path=pool_summary_path,
        bls_path=bls_path,
        bls_summary_path=bls_summary_path,
    )
    exclusions_path = tmp_path / "native_model_exclusions.csv"
    summary_path = tmp_path / "native_model_eligibility.summary.json"
    write_native_model_eligibility(
        pool_path=pool_path,
        pool_summary_path=pool_summary_path,
        bls_path=bls_path,
        bls_summary_path=bls_summary_path,
        exclusions_path=exclusions_path,
        summary_path=summary_path,
        production_lock=True,
    )

    # Native shard consumers intentionally avoid reopening the global Parquet,
    # but a declared production release must still bind its exact known digest.
    monkeypatch.setattr(
        eligibility_module,
        "PRODUCTION_GLOBAL_BLS_SHA256",
        "0" * 64,
    )
    with pytest.raises(ValueError, match="source artifacts differ"):
        load_native_model_eligibility(
            exclusions_path,
            summary_path,
            pool_path=pool_path,
            pool_summary_path=pool_summary_path,
            rederive_from_bls=False,
        )


def test_production_loader_rechecks_loaded_bls_contract(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    pool_path, pool_summary_path, bls_path, bls_summary_path, *_ = (
        _write_source_authorities(tmp_path)
    )
    _lock_toy_production_authority(
        monkeypatch,
        pool_summary_path=pool_summary_path,
        bls_path=bls_path,
        bls_summary_path=bls_summary_path,
    )
    exclusions_path = tmp_path / "native_model_exclusions.csv"
    summary_path = tmp_path / "native_model_eligibility.summary.json"
    write_native_model_eligibility(
        pool_path=pool_path,
        pool_summary_path=pool_summary_path,
        bls_path=bls_path,
        bls_summary_path=bls_summary_path,
        exclusions_path=exclusions_path,
        summary_path=summary_path,
        production_lock=True,
    )

    monkeypatch.setattr(
        eligibility_module,
        "_production_bls_contract",
        lambda: {"apertures": list(ADP_ONLY_APERTURES), "tampered": True},
    )
    with pytest.raises(ValueError, match="BLS contract differs"):
        load_native_model_eligibility(
            exclusions_path,
            summary_path,
            pool_path=pool_path,
            pool_summary_path=pool_summary_path,
            bls_path=bls_path,
            bls_summary_path=bls_summary_path,
        )
