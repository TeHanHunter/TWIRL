from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest

import twirl.injections.a2v1_recovery as recovery
from twirl.injections.a2v1_recovery import (
    A2V1RecoveryConfig,
    FRESH_INJECTION_CONTRACT,
    audit_fresh_injection_shards,
    build_fresh_injection_schedule,
    compare_adp_compact_products,
    fresh_injection_contract,
    run_adp_roundtrip_parity,
    schedule_contract,
    select_parameter_spanning_rows,
    write_fresh_injection_shard,
)
from twirl.vetting.injection_teacher_recovery import (
    build_teacher_injection_holdout,
    normalize_injection_peak_candidates,
)
from twirl.vetting.harmonic_inference import prepare_inference_rows


def _config(
    sector: int = 56,
    shard_assignment: str = "contiguous_grid",
) -> A2V1RecoveryConfig:
    return A2V1RecoveryConfig(
        name="test",
        sector=sector,
        n_injections=8,
        n_shards=2,
        shard_assignment=shard_assignment,
        period_bins=2,
        radius_bins=2,
        repeats_per_cell=2,
        min_good_cadences=10,
        parity_sample_size=4,
        parity_median_abs_limit=1.0e-4,
        parity_scatter_ratio_tolerance=0.01,
    )


def _write_sources(
    tmp_path: Path,
    n_targets: int = 11,
    sector: int = 56,
) -> tuple[Path, Path]:
    raw_path = tmp_path / "raw.h5"
    adp_path = tmp_path / "adp.h5"
    n = 80
    time_rel = 2830.0 + np.arange(n) * 200.0 / 86400.0
    cadence = np.arange(1000, 1000 + n, dtype=np.int64)
    orbit = np.where(np.arange(n) < n // 2, 119, 120).astype(np.int32)
    quality = np.zeros(n, dtype=np.int32)
    with h5py.File(raw_path, "w") as raw_file, h5py.File(adp_path, "w") as adp_file:
        raw_file.attrs["contract_version"] = "s56_tglc_raw_pair_v1"
        raw_root = raw_file.create_group("targets")
        adp_root = adp_file.create_group("targets")
        for index in range(n_targets):
            tic = 10_000 + index
            key = f"{tic:016d}"
            baseline_small = 100.0 + index
            baseline_primary = 150.0 + index
            trend = 0.5 * np.sin(np.linspace(0, 3.0 * np.pi, n))
            raw_small = baseline_small + trend
            raw_primary = baseline_primary + 1.2 * trend
            if index == 0:
                raw_small[0] = -5.0
                raw_primary[0] = -7.0
            err_small = np.full(n, 2.0)
            err_primary = np.full(n, 2.5)
            det_small, _ = recovery._adp_detrend(
                time_rel, raw_small, err_small, quality
            )
            det_primary, _ = recovery._adp_detrend(
                time_rel, raw_primary, err_primary, quality
            )
            raw = raw_root.create_group(key)
            adp = adp_root.create_group(key)
            raw.attrs["tic"] = tic
            adp.attrs.update(
                {
                    "tic": tic,
                    "sector": sector,
                    "camera": 1,
                    "ccd": 1,
                    "tessmag": 16.5 + 0.35 * index,
                    "source_fits": f"source/{tic}.fits",
                }
            )
            for group, values in (
                (
                    raw,
                    {
                        "time": time_rel + 2457000.0,
                        "cadenceno": cadence,
                        "orbitid": orbit,
                        "quality": quality,
                        "raw_flux_small": raw_small,
                        "raw_flux_err_small": err_small,
                        "raw_flux_primary": raw_primary,
                        "raw_flux_err_primary": err_primary,
                    },
                ),
                (
                    adp,
                    {
                        "time": time_rel,
                        "cadenceno": cadence.astype(np.int32),
                        "orbitid": orbit.astype(np.int16),
                        "quality": quality,
                        "DET_FLUX_ADP_SML": det_small,
                        "DET_FLUX_ADP": det_primary,
                    },
                ),
            ):
                for name, value in values.items():
                    group.create_dataset(name, data=value)
    return raw_path, adp_path


def test_fresh_schedule_is_unique_host_disjoint_and_grid_complete(
    tmp_path: Path,
) -> None:
    raw_h5, adp_h5 = _write_sources(tmp_path)
    teacher = tmp_path / "teacher.csv"
    pd.DataFrame(
        {"tic": [10_000, 10_000], "human_label": ["planet_like", "planet_like"]}
    ).to_csv(teacher, index=False)
    schedule, qa, summary = build_fresh_injection_schedule(
        raw_h5=raw_h5,
        adp_h5=adp_h5,
        teacher_table=teacher,
        config=_config(),
        out_dir=tmp_path / "schedule",
    )
    assert len(schedule) == 8
    assert schedule["tic"].nunique() == 8
    assert 10_000 not in set(schedule["tic"])
    assert schedule.groupby(["grid_period_bin", "grid_radius_bin"]).size().eq(2).all()
    assert schedule.groupby("shard_index").size().to_dict() == {0: 4, 1: 4}
    assert qa.loc[qa["tic"].eq(10_000), "qa_reason"].item() == "teacher_tic"
    assert summary["n_teacher_unique_tics"] == 1


def test_sector_57_schedule_has_sector_specific_contracts_and_ids(
    tmp_path: Path,
) -> None:
    raw_h5, adp_h5 = _write_sources(tmp_path, sector=57)
    teacher = tmp_path / "teacher.csv"
    prior_sector = tmp_path / "prior_sector.parquet"
    pd.DataFrame({"tic": [99_999]}).to_csv(teacher, index=False)
    pd.DataFrame({"tic": [10_000]}).to_parquet(prior_sector, index=False)
    schedule, _, summary = build_fresh_injection_schedule(
        raw_h5=raw_h5,
        adp_h5=adp_h5,
        teacher_table=teacher,
        config=_config(57, "balanced_random"),
        out_dir=tmp_path / "schedule",
        host_overlap_audit_tables=[prior_sector],
    )
    assert schedule["sector"].eq(57).all()
    assert schedule["injection_id"].str.match(r"s57a2v1_eval_\d{6}").all()
    assert schedule["schedule_contract"].eq(schedule_contract(57)).all()
    assert summary["contract"] == schedule_contract(57)
    assert summary["injection_contract"] == fresh_injection_contract(57)
    assert summary["shard_assignment"] == "balanced_random"
    assert schedule.groupby("shard_index").size().eq(4).all()
    assert not schedule["shard_index"].is_monotonic_increasing
    assert summary["host_overlap_audits"] == [
        {
            "table": str(prior_sector.resolve()),
            "n_comparison_unique_tics": 1,
            "n_selected_host_overlap": 1,
        }
    ]


def test_schedule_excludes_prior_sector_evaluation_hosts(tmp_path: Path) -> None:
    raw_h5, adp_h5 = _write_sources(tmp_path)
    teacher = tmp_path / "teacher.csv"
    prior = tmp_path / "prior.parquet"
    pd.DataFrame({"tic": [10_000]}).to_csv(teacher, index=False)
    pd.DataFrame({"tic": [10_001, 10_002]}).to_parquet(prior, index=False)
    schedule, qa, summary = build_fresh_injection_schedule(
        raw_h5=raw_h5,
        adp_h5=adp_h5,
        teacher_table=teacher,
        config=_config(),
        out_dir=tmp_path / "schedule",
        additional_exclusion_tables=[prior],
    )
    assert not set(schedule["tic"]) & {10_000, 10_001, 10_002}
    assert set(
        qa.loc[qa["qa_reason"].eq("prior_evaluation_tic"), "tic"].astype(int)
    ) == {10_001, 10_002}
    assert summary["n_teacher_unique_tics"] == 1
    assert summary["n_prior_evaluation_unique_tics"] == 2
    assert summary["n_total_excluded_unique_tics"] == 3


def test_parameter_spanning_smoke_selection_balances_grid_marginals() -> None:
    rows = []
    index = 0
    for period_bin in range(5):
        for radius_bin in range(5):
            for slot in range(2):
                rows.append(
                    {
                        "grid_period_bin": period_bin,
                        "grid_radius_bin": radius_bin,
                        "grid_slot": slot,
                        "injection_index": index,
                        "tic": 1000 + index,
                    }
                )
                index += 1
    selected = select_parameter_spanning_rows(pd.DataFrame(rows), n_rows=10)
    assert len(selected) == 10
    assert selected["tic"].nunique() == 10
    assert selected["grid_period_bin"].value_counts().eq(2).all()
    assert selected["grid_radius_bin"].value_counts().eq(2).all()
    assert not selected.duplicated(["grid_period_bin", "grid_radius_bin"]).any()


def test_balanced_full_schedule_spans_every_grid_axis_per_shard() -> None:
    config = A2V1RecoveryConfig(sector=57, seed=570201, shard_assignment="balanced_random")
    schedule = pd.DataFrame({"injection_index": range(config.n_injections)})
    schedule["grid_cell_index"] = (
        schedule["injection_index"] // config.repeats_per_cell
    )
    schedule["grid_period_bin"] = (
        schedule["grid_cell_index"] // config.radius_bins
    )
    schedule["grid_radius_bin"] = (
        schedule["grid_cell_index"] % config.radius_bins
    )
    schedule["shard_index"] = recovery._assign_shards(schedule, config=config)
    grouped = schedule.groupby("shard_index")
    assert grouped.size().eq(config.rows_per_shard).all()
    assert grouped["grid_period_bin"].nunique().eq(config.period_bins).all()
    assert grouped["grid_radius_bin"].nunique().eq(config.radius_bins).all()


def test_adp_roundtrip_parity_uses_unmodified_raw_flux(tmp_path: Path) -> None:
    raw_h5, adp_h5 = _write_sources(tmp_path)
    teacher = tmp_path / "teacher.csv"
    pd.DataFrame({"tic": [10_010]}).to_csv(teacher, index=False)
    schedule, _, _ = build_fresh_injection_schedule(
        raw_h5=raw_h5,
        adp_h5=adp_h5,
        teacher_table=teacher,
        config=_config(),
        out_dir=tmp_path / "schedule",
    )
    metrics, summary = run_adp_roundtrip_parity(
        raw_h5=raw_h5,
        adp_h5=adp_h5,
        schedule=schedule,
        config=_config(),
        out_dir=tmp_path / "parity",
    )
    assert len(metrics) == 8
    assert summary["passed"]
    assert metrics["cadence_match"].all()
    assert metrics["median_abs_difference"].max() <= 1.0e-4


def test_rebuilt_compact_adp_requires_exact_array_parity(tmp_path: Path) -> None:
    _, reference = _write_sources(tmp_path)
    active = tmp_path / "active.h5"
    active.write_bytes(reference.read_bytes())
    comparison = compare_adp_compact_products(reference, active, progress_every=0)
    assert comparison["passed"]
    assert comparison["n_targets_compared"] == 11

    with h5py.File(active, "r+") as h5:
        dataset = h5["targets/0000000000010000/DET_FLUX_ADP_SML"]
        dataset[0] = float(dataset[0]) + 0.01
    comparison = compare_adp_compact_products(reference, active, progress_every=0)
    assert not comparison["passed"]
    assert comparison["n_mismatched_targets"] == 1


def test_injection_shard_preserves_negative_flux_and_adjusts_errors(
    tmp_path: Path, monkeypatch
) -> None:
    raw_h5, adp_h5 = _write_sources(tmp_path)
    teacher = tmp_path / "teacher.csv"
    pd.DataFrame({"tic": [99_999]}).to_csv(teacher, index=False)
    schedule, _, _ = build_fresh_injection_schedule(
        raw_h5=raw_h5,
        adp_h5=adp_h5,
        teacher_table=teacher,
        config=_config(),
        out_dir=tmp_path / "schedule",
    )

    def fake_inject(time, flux, **kwargs):
        model = np.ones(len(time), dtype=float)
        model[::9] = 0.5
        mask = model < 1.0
        baseline = float(kwargs["baseline"])
        return np.asarray(flux) + baseline * (model - 1.0), mask, model

    monkeypatch.setattr(recovery, "inject_batman_transit", fake_inject)
    shard_index = int(schedule.loc[schedule["tic"].eq(10_000), "shard_index"].item())
    out_h5 = tmp_path / f"injections_{shard_index:02d}.h5"
    summary = write_fresh_injection_shard(
        raw_h5=raw_h5,
        adp_h5=adp_h5,
        schedule=schedule,
        shard_index=shard_index,
        config=_config(),
        out_h5=out_h5,
    )
    assert summary["n_injections"] == 4
    with h5py.File(out_h5, "r") as h5:
        assert h5.attrs["contract_version"] == FRESH_INJECTION_CONTRACT
        assert len(h5["injections"]) == 4
        group = next(
            group
            for group in h5["injections"].values()
            if int(group.attrs["tic"]) == 10_000
        )
        original = np.asarray(group["RAW_FLUX_Small_original"])
        injected = np.asarray(group["RAW_FLUX_Small_injected"])
        original_error = np.asarray(group["RAW_FLUX_ERR_Small_original"])
        injected_error = np.asarray(group["RAW_FLUX_ERR_Small_injected"])
        model = np.asarray(group["transit_model"])
        assert np.any(injected < 0)
        np.testing.assert_allclose(
            injected,
            original + float(group.attrs["injection_baseline_Small"]) * (model - 1.0),
        )
        assert np.all(injected_error[model < 1.0] <= original_error[model < 1.0])
    shard_schedule = schedule.loc[schedule["shard_index"].eq(shard_index)].copy()
    audit = audit_fresh_injection_shards(
        shard_paths=[out_h5],
        schedule=shard_schedule,
        raw_h5=raw_h5,
        adp_h5=adp_h5,
        config=_config(),
    )
    assert audit["passed"]
    assert audit["n_groups"] == 4
    assert audit["n_alignment_failures"] == 0
    assert audit["n_original_copy_failures"] == 0
    assert audit["n_negative_original_points"] >= 2
    assert audit["n_negative_preservation_failures"] == 0


def test_existing_holdout_helper_excludes_all_teacher_hosts() -> None:
    manifest = pd.DataFrame(
        {
            "injection_id": ["a", "b", "c", "d"],
            "tic": [1, 2, 3, 4],
            "tessmag": [17.0, 18.0, 19.0, 20.0],
            "grid_period_bin": [0, 0, 1, 1],
            "grid_radius_bin": [0, 0, 0, 0],
        }
    )
    teacher = pd.DataFrame(
        {
            "injection_id": ["z"],
            "tic": [2],
            "is_injected_row": [False],
        }
    )
    retained, _, summary = build_teacher_injection_holdout(
        manifest,
        teacher,
        n_shards=2,
        min_cell_support=1,
    )
    assert set(retained["tic"]) == {1, 3, 4}
    assert summary["n_retained_on_any_teacher_table_host"] == 0


def test_teacher_inference_requires_explicit_injection_opt_in() -> None:
    candidates = pd.DataFrame(
        {
            "tic": [42],
            "sector": [56],
            "period_d": [1.0],
            "t0_bjd": [2459830.0],
            "duration_min": [5.0],
            "review_id": ["injection-42"],
            "source_kind": ["injected_validation_holdout"],
            "is_injected_row": [True],
            "injection_id": ["s56a2v1_eval_000042"],
        }
    )
    with pytest.raises(ValueError, match="real target groups only"):
        prepare_inference_rows(candidates)
    rows = prepare_inference_rows(candidates, allow_injections=True)
    assert rows.loc[0, "native_group_path"] == "injections/s56a2v1_eval_000042"


def test_sector_57_candidate_ids_and_source_tag_derive_from_pair(tmp_path: Path) -> None:
    pair_h5 = tmp_path / "pair.h5"
    with h5py.File(pair_h5, "w") as h5:
        root = h5.create_group("injections")
        group = root.create_group("s57a2v1_eval_000001")
        group.attrs["sector"] = 57
    rows = []
    for aperture in ("DET_FLUX_ADP_SML", "DET_FLUX_ADP"):
        rows.append(
            {
                "injection_id": "s57a2v1_eval_000001",
                "aperture": aperture,
                "status": "ok",
                "peak_rank": 1,
                "period_d": 1.0,
                "t0_bjd": 2_459_830.0,
                "duration_min": 5.0,
                "depth": 0.1,
                "depth_snr": 10.0,
                "sde": 12.0,
                "log_power": 3.0,
            }
        )
    candidates = normalize_injection_peak_candidates(
        pd.DataFrame(rows), pair_h5=pair_h5
    )
    assert candidates["review_id"].tolist() == [
        "s0057-a2v1-eval-s57a2v1_eval_000001-r1"
    ]
    assert candidates["sector"].eq(57).all()
    assert candidates["source_product_tag"].eq("S57_A2v1_fresh_eval_v1").all()
