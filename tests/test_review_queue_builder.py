from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import h5py


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_queue_builder():
    path = REPO_ROOT / "scripts" / "stage5_validation" / "build_s56_pretriage_review_queue.py"
    spec = importlib.util.spec_from_file_location("build_s56_pretriage_review_queue", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _load_recovery_sweep():
    path = REPO_ROOT / "scripts" / "stage5_validation" / "run_injection_recovery_mode_sweep.py"
    spec = importlib.util.spec_from_file_location("run_injection_recovery_mode_sweep", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    import sys

    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _load_snr_queue_builder():
    path = REPO_ROOT / "scripts" / "stage5_validation" / "build_snr_stratified_injection_review_queue.py"
    spec = importlib.util.spec_from_file_location("build_snr_stratified_injection_review_queue", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _load_bls_failure_diagnostics():
    path = REPO_ROOT / "scripts" / "stage5_validation" / "diagnose_bls_recovery_failures.py"
    spec = importlib.util.spec_from_file_location("diagnose_bls_recovery_failures", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _load_aperture_survival_summary():
    path = REPO_ROOT / "scripts" / "stage5_validation" / "summarize_predetrend_aperture_survival.py"
    spec = importlib.util.spec_from_file_location("summarize_predetrend_aperture_survival", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _load_predetrend_debug_queue_builder():
    path = REPO_ROOT / "scripts" / "stage5_validation" / "build_predetrend_debug_priority_queues.py"
    spec = importlib.util.spec_from_file_location("build_predetrend_debug_priority_queues", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _load_priority_sweep_summary():
    path = REPO_ROOT / "scripts" / "stage5_validation" / "summarize_priority_recovery_sweep.py"
    spec = importlib.util.spec_from_file_location("summarize_priority_recovery_sweep", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _load_priority_sweep_preflight():
    path = REPO_ROOT / "scripts" / "stage5_validation" / "check_s56_predetrend_priority_sweep_inputs.py"
    spec = importlib.util.spec_from_file_location("check_s56_predetrend_priority_sweep_inputs", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _load_peak_training_builder():
    path = REPO_ROOT / "scripts" / "stage5_validation" / "build_injection_peak_training_table.py"
    spec = importlib.util.spec_from_file_location("build_injection_peak_training_table", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _load_detrending_method_sweep():
    path = REPO_ROOT / "scripts" / "stage5_validation" / "sweep_predetrend_detrending_methods.py"
    spec = importlib.util.spec_from_file_location("sweep_predetrend_detrending_methods", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    import sys

    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_stratified_real_selector_tolerates_missing_aperture_agreement(tmp_path) -> None:
    module = _load_queue_builder()
    rows = []
    classes = ["planet_candidate", "sub_roche_pceb_suspect", "pceb_grid_ceiling", "systematic"]
    for idx in range(40):
        rows.append(
            {
                "tic": 100000 + idx,
                "sector": 56,
                "cam": 1 + idx % 4,
                "ccd": 1 + idx % 4,
                "tmag": 15.0 + 0.1 * idx,
                "vet_class": classes[idx % len(classes)],
                "class_rank": idx + 1,
                "blind_rank": idx + 5,
                "period_d": 0.10 + 0.15 * (idx % 20),
                "t0_bjd": 2459825.0 + idx,
                "duration_min": 5.0 + idx % 7,
                "depth": 0.01 + 0.002 * idx,
                "depth_snr": 3.0 + idx,
                "sde_max": 7.0 + idx,
                "rep_aperture": "DET_FLUX",
                "apertures_agree": "DET_FLUX",
                "centroid_status": "on_target",
                "centroid_pass": True,
                "centroid_delta_pix": 0.01,
                "centroid_z": 0.5,
                "n_in_transit": 4,
                "n_oot_band": 50,
            }
        )
    path = tmp_path / "candidates.csv"
    pd.DataFrame(rows).to_csv(path, index=False)

    selected = module.select_real_candidates_stratified(path, n_real=15, random_state=5)

    assert len(selected) == 15
    assert selected["tic"].nunique() == 15
    assert set(selected["source_kind"]) == {"real_candidate"}
    assert selected["truth_source_bucket"].astype(str).str.startswith("real_").all()


def test_select_ranker_real_candidates_keeps_peak_ephemerides(tmp_path) -> None:
    module = _load_queue_builder()
    path = tmp_path / "ranker_selected.csv"
    pd.DataFrame(
        [
            {
                "tic": 101,
                "sector": 56,
                "cam": 1,
                "ccd": 2,
                "tmag": 16.1,
                "ranker_selection_rank": 1,
                "aperture": "DET_FLUX_ADP",
                "peak_rank": 3,
                "period_d": 0.75,
                "t0_bjd": 2459825.1,
                "duration_min": 8.0,
                "depth": 0.04,
                "depth_snr": 9.0,
                "sde": 12.0,
                "ranker_p_background_peak": 0.09,
                "ranker_p_signal_peak": 0.91,
            },
            {
                "tic": 101,
                "sector": 56,
                "cam": 1,
                "ccd": 2,
                "tmag": 16.1,
                "ranker_selection_rank": 2,
                "aperture": "DET_FLUX",
                "peak_rank": 1,
                "period_d": 1.50,
                "t0_bjd": 2459825.3,
                "duration_min": 10.0,
                "depth": 0.03,
                "depth_snr": 8.0,
                "sde": 11.0,
                "ranker_p_background_peak": 0.35,
                "ranker_p_signal_peak": 0.65,
            },
            {
                "tic": 202,
                "sector": 56,
                "cam": 2,
                "ccd": 1,
                "tmag": 18.2,
                "ranker_selection_rank": 1,
                "aperture": "DET_FLUX_ADP",
                "peak_rank": 2,
                "period_d": 0.25,
                "t0_bjd": 2459826.0,
                "duration_min": 6.0,
                "depth": 0.09,
                "depth_snr": 7.0,
                "sde": 10.0,
                "ranker_p_background_peak": 0.18,
                "ranker_p_signal_peak": 0.82,
            },
            {
                "tic": 303,
                "sector": 56,
                "cam": 3,
                "ccd": 1,
                "tmag": 19.2,
                "ranker_selection_rank": 1,
                "aperture": "DET_FLUX_SML",
                "peak_rank": 2,
                "period_d": 0.35,
                "t0_bjd": 2459826.2,
                "duration_min": 7.0,
                "depth": 0.08,
                "depth_snr": 6.0,
                "sde": 9.0,
                "ranker_p_background_peak": 0.99,
                "ranker_p_signal_peak": 0.70,
            },
        ]
    ).to_csv(path, index=False)

    top_rank1 = module.select_ranker_real_candidates(path, n_real=2, random_state=56)

    assert set(top_rank1["tic"]) == {101, 202}
    assert 303 not in set(top_rank1["tic"])

    selected = module.select_ranker_real_candidates(path, n_real=4, random_state=56)

    assert len(selected) == 4
    assert set(selected["source_kind"]) == {"real_candidate"}
    assert set(selected["recovery_status"]) == {"real_ranker_selected"}
    assert selected["review_id"].nunique() == 4
    assert selected["tic"].nunique() == 3
    assert selected["review_id"].str.contains(":ranker:").all()
    assert set(selected["rep_aperture"]) == {"DET_FLUX_ADP", "DET_FLUX", "DET_FLUX_SML"}
    assert set(selected["sde_max"]) == {9.0, 10.0, 11.0, 12.0}
    assert selected["truth_source_bucket"].eq("real_ranker_selected").all()


def test_blind_review_metadata_preserves_truth_columns() -> None:
    module = _load_queue_builder()
    queue = pd.DataFrame(
        [
            {
                "tic": 1,
                "source_bucket": "injection_wd_giant_or_bd_bls_recovered",
                "source_kind": "injection_recovery",
                "vet_class": "injected_wd_giant_or_bd",
                "truth_source_bucket": "injection_wd_giant_or_bd_bls_recovered",
                "truth_source_kind": "injection_recovery",
                "truth_radius_rearth": np.float64(9.0),
            },
            {
                "tic": 2,
                "source_bucket": "real_top_sde",
                "source_kind": "real_candidate",
                "vet_class": "planet_candidate",
                "truth_source_bucket": "real_top_sde",
                "truth_source_kind": "real_candidate",
                "truth_radius_rearth": "",
            },
        ]
    )

    blinded = module._blind_review_metadata(queue)

    assert blinded["source_bucket"].tolist() == ["review_candidate", "review_candidate"]
    assert blinded.loc[0, "vet_class"] == "review_candidate"
    assert blinded.loc[1, "vet_class"] == "planet_candidate"
    assert blinded.loc[0, "truth_source_bucket"] == "injection_wd_giant_or_bd_bls_recovered"
    assert blinded.loc[0, "truth_radius_rearth"] == 9.0


def test_catalog_star_annotation_and_leo_star_dict() -> None:
    module = _load_queue_builder()
    row = {
        "source_id": 2146576589564898688,
        "teff_H": 4738.07,
        "eteff_H": 54.0,
        "logg_H": 7.87774,
        "elogg_H": 0.03,
        "mass_H": 0.506161,
        "emass_H": 0.02,
    }
    record = module._star_record_from_catalog_row(
        row,
        tic=267574918,
        source_id=row["source_id"],
        atmosphere="H",
    )
    assert record is not None
    assert record["star_source"] == "GF21_H"
    assert np.isclose(record["star_rad_rsun"], 0.01356, rtol=5.0e-3)

    queue = pd.DataFrame({"tic": [267574918, 123], "period_d": [1.0, 2.0]})
    annotated = module.annotate_star_parameters(
        queue,
        star_records={267574918: record},
    )

    assert annotated.loc[0, "star_source"] == "GF21_H"
    assert annotated.loc[1, "star_source"] == "canonical_fallback"
    assert module.star_source_counts(annotated) == {"GF21_H": 1, "canonical_fallback": 1}

    star = module._wd_star_for(annotated.loc[0])
    assert star["source_id"] == str(row["source_id"])
    assert np.isclose(star["Rs"], annotated.loc[0, "star_rad_rsun"])
    fallback = module._wd_star_for(annotated.loc[1])
    assert fallback["star_source"] == "canonical_fallback"
    assert fallback["Rs"] == module._CANONICAL_WD_STAR["Rs"]


def test_injection_transit_counts_fall_back_to_mask_dataset(tmp_path) -> None:
    module = _load_queue_builder()
    path = tmp_path / "injections.h5"
    with h5py.File(path, "w") as h5:
        group = h5.create_group("injections/inj_000000")
        group.create_dataset("in_transit", data=np.array([False, True, True, False, True]))
        group.create_dataset("quality", data=np.array([0, 0, 1, 0, 0]))

    with h5py.File(path, "r") as h5:
        group = h5["injections/inj_000000"]
        n_in, n_good = module._injection_transit_counts(group, dict(group.attrs.items()))

    assert n_in == 3
    assert n_good == 2


def test_topn_truth_match_summary_separates_exact_and_harmonic_peaks() -> None:
    module = _load_queue_builder()
    truth_t0 = 2459825.0
    peak_rows = [
        {
            "status": "ok",
            "peak_rank": 1,
            "period_d": 5.01,
            "t0_bjd": truth_t0,
            "sde": 40.0,
        },
        {
            "status": "ok",
            "peak_rank": 2,
            "period_d": 1.0005,
            "t0_bjd": truth_t0 + 1.0 / 1440.0,
            "sde": 22.0,
        },
    ]

    summary = module._topn_truth_match_summary(
        peak_rows,
        truth_period_d=1.0,
        truth_t0_bjd=truth_t0,
        truth_duration_min=5.0,
    )

    assert summary["top_harmonic_ephemeris_match"] is True
    assert summary["topn_exact_recovered"] is True
    assert summary["topn_exact_peak_rank"] == 2
    assert np.isclose(summary["topn_exact_sde"], 22.0)
    assert summary["topn_harmonic_match"] is True
    assert summary["topn_harmonic_peak_rank"] == 1
    assert np.isclose(summary["topn_harmonic_factor"], 5.0)


def test_peak_training_truth_labels_exact_harmonic_and_mismatch() -> None:
    module = _load_peak_training_builder()

    exact = module.label_peak_against_injection(
        period_d=1.0005,
        t0_bjd=2459825.001,
        duration_min=6.0,
        truth_period_d=1.0,
        truth_t0_bjd=2459825.0,
        truth_duration_min=6.0,
    )
    assert exact["match_kind"] == "exact"
    assert exact["exact_ephemeris_match"] is True
    assert exact["is_injected_signal_peak"] is True

    harmonic = module.label_peak_against_injection(
        period_d=2.001,
        t0_bjd=2459825.0,
        duration_min=6.0,
        truth_period_d=1.0,
        truth_t0_bjd=2459825.0,
        truth_duration_min=6.0,
    )
    assert harmonic["match_kind"] == "harmonic"
    assert harmonic["harmonic_ephemeris_match"] is True
    assert np.isclose(harmonic["nearest_harmonic_factor"], 2.0)

    mismatch = module.label_peak_against_injection(
        period_d=1.25,
        t0_bjd=2459825.0,
        duration_min=6.0,
        truth_period_d=1.0,
        truth_t0_bjd=2459825.0,
        truth_duration_min=6.0,
    )
    assert mismatch["match_kind"] == "mismatch"
    assert mismatch["is_injected_signal_peak"] is False


def test_peak_training_summary_reports_recall_at_k() -> None:
    module = _load_peak_training_builder()
    df = pd.DataFrame(
        {
            "injection_id": ["a", "a", "b", "b", "c"],
            "is_candidate_peak": [True, True, True, True, False],
            "peak_rank": [1, 2, 1, 2, 0],
            "is_injected_signal_peak": [False, True, False, False, False],
            "exact_ephemeris_match": [False, True, False, False, False],
            "harmonic_ephemeris_match": [False, False, False, False, False],
            "match_kind": ["mismatch", "exact", "mismatch", "mismatch", "no_peak"],
            "search_branch": ["standard", "short_pmax2", "standard", "standard", "short_pmax2"],
        }
    )

    summary = module.summarize_peak_table(df)

    assert summary["n_injections"] == 3
    assert summary["n_any_match_injections"] == 1
    assert summary["recall_at_1"]["n"] == 0
    assert summary["recall_at_2"]["n"] == 1
    assert summary["match_kind_counts"]["exact"] == 1
    assert summary["search_branch_counts"] == {"short_pmax2": 2, "standard": 3}


def test_recovery_sweep_parser_accepts_duration_overrides() -> None:
    module = _load_recovery_sweep()

    cfg = module.parse_sweep("adp_short:DET_FLUX_ADP+DET_FLUX_ADP_SML:200000:1+2+3+5")

    assert cfg.name == "adp_short"
    assert cfg.apertures == ("DET_FLUX_ADP", "DET_FLUX_ADP_SML")
    assert cfg.n_periods == 200000
    assert cfg.durations_min == (1.0, 2.0, 3.0, 5.0)


def test_recovery_sweep_loads_injection_id_file(tmp_path) -> None:
    module = _load_recovery_sweep()
    path = tmp_path / "ids.txt"
    path.write_text("\n# comment\npredet_000002\npredet_000001,extra\npredet_000002\n\n")

    ids = module.load_injection_ids(path)

    assert ids == ("predet_000002", "predet_000001")


def test_snr_stratified_queue_preserves_rare_modes(tmp_path) -> None:
    module = _load_snr_queue_builder()
    ids = [f"inj_{idx:06d}" for idx in range(6)]
    base = pd.DataFrame(
        {
            "injection_id": ids,
            "review_id": [f"inj:{item}" for item in ids],
            "tic": np.arange(100, 106),
            "recovery_status": ["bls_recovered", "bls_peak_mismatch", "bls_peak_mismatch", "bls_peak_mismatch", "bls_peak_mismatch", "bls_peak_mismatch"],
            "leo_report_name": [f"row{idx}.pdf" for idx in range(6)],
        }
    )
    topn = pd.DataFrame(
        {
            "injection_id": ids,
            "topn_recovery_status": [
                "bls_top1_recovered",
                "bls_topn_recovered",
                "bls_topn_harmonic_match",
                "bls_peak_mismatch",
                "bls_peak_mismatch",
                "bls_peak_mismatch",
            ],
        }
    )
    survival = pd.DataFrame(
        {
            "injection_id": ids,
            "aperture": ["DET_FLUX_ADP"] * len(ids),
            "delta_signal_multi_snr_mad": [30.0, 12.0, 8.0, 0.5, 4.0, 15.0],
            "depth_retention_frac": [1.0] * len(ids),
        }
    )
    base_path = tmp_path / "base.csv"
    topn_path = tmp_path / "topn.csv"
    survival_path = tmp_path / "survival.csv"
    base.to_csv(base_path, index=False)
    topn.to_csv(topn_path, index=False)
    survival.to_csv(survival_path, index=False)

    selected = module.build_queue(
        base_queue_csv=base_path,
        topn_queue_csv=topn_path,
        survival_csv=survival_path,
        out_dir=tmp_path / "out",
        aperture="DET_FLUX_ADP",
        n_rows=5,
        random_state=1,
        min_per_stratum=1,
        rare_cap=5,
        reports_mode="none",
    )

    assert len(selected) == 5
    assert {"bls_top1_recovered", "bls_topn_recovered", "bls_topn_harmonic_match"}.issubset(
        set(selected["recovery_mode"])
    )
    assert (tmp_path / "out" / "review_queue.csv").exists()
    assert (tmp_path / "out" / "summary.json").exists()


def test_bls_failure_diagnostic_summarizes_recovery_modes(tmp_path) -> None:
    module = _load_bls_failure_diagnostics()
    ids = [f"inj_{idx:06d}" for idx in range(5)]
    queue = pd.DataFrame(
        {
            "injection_id": ids,
            "recovery_status": [
                "bls_recovered",
                "bls_peak_mismatch",
                "bls_peak_mismatch",
                "bls_peak_mismatch",
                "bls_peak_mismatch",
            ],
            "topn_recovery_status": [
                "bls_top1_recovered",
                "bls_topn_recovered",
                "bls_topn_harmonic_match",
                "bls_peak_mismatch",
                "bls_peak_mismatch",
            ],
            "truth_period_d": [0.2, 0.4, 1.2, 8.0, 0.8],
            "truth_model_depth": [0.02, 0.05, 0.12, 0.2, 0.1],
            "truth_radius_rearth": [1.5, 5.0, 10.0, 18.0, 8.0],
            "tmag": [15.5, 17.5, 19.2, 20.5, 18.2],
        }
    )
    survival = pd.DataFrame(
        {
            "injection_id": ids,
            "aperture": ["DET_FLUX_ADP"] * len(ids),
            "delta_signal_multi_snr_mad": [25.0, 8.0, 4.0, 0.5, 9.0],
            "delta_signal_snr_mad": [18.0, 6.0, 3.0, 0.4, 8.5],
            "depth_retention_frac": [0.9, 0.6, 0.4, 0.1, 0.5],
            "n_good_ap_in_transit": [20, 12, 8, 5, 10],
        }
    )
    queue_path = tmp_path / "queue.csv"
    survival_path = tmp_path / "survival.csv"
    out_dir = tmp_path / "out"
    queue.to_csv(queue_path, index=False)
    survival.to_csv(survival_path, index=False)

    summary = module.diagnose(queue_path, survival_path, "DET_FLUX_ADP", out_dir)

    assert summary["recovery_mode_counts"]["bls_top1_recovered"] == 1
    assert summary["recovery_mode_counts"]["bls_topn_recovered"] == 1
    assert summary["recovery_mode_counts"]["bls_topn_harmonic_match"] == 1
    assert summary["recovery_mode_counts"]["bls_peak_mismatch"] == 2
    assert summary["failure_class_counts"]["low_snr_unmatched"] == 1
    assert summary["failure_class_counts"]["snr_qualified_bls_ranking_loss"] == 1
    assert np.isclose(summary["top1_recovered_fraction"], 0.2)
    assert np.isclose(summary["any_period_match_fraction"], 0.6)
    by_snr = pd.read_csv(out_dir / "recovery_by_snr_bin.csv")
    assert "any_period_match_frac" in by_snr.columns
    assert (out_dir / "failure_class_counts.csv").exists()
    assert (out_dir / "snr_qualified_bls_ranking_loss_cases.csv").exists()
    assert (out_dir / "bls_recovery_failure_joined.csv").exists()
    assert (out_dir / "summary.json").exists()


def test_aperture_survival_summary_finds_best_aperture_gain(tmp_path) -> None:
    module = _load_aperture_survival_summary()
    rows = []
    for inj, tmag, current_snr, small_snr, lag_snr in [
        ("inj_000000", 17.0, 8.0, 9.0, 2.0),
        ("inj_000001", 19.0, 2.0, 10.0, 1.0),
        ("inj_000002", 20.0, 0.5, 0.7, 0.2),
    ]:
        for aperture, snr in [
            ("DET_FLUX_ADP", current_snr),
            ("DET_FLUX_ADP_SML", small_snr),
            ("DET_FLUX_ADP_LAG", lag_snr),
        ]:
            rows.append(
                {
                    "injection_id": inj,
                    "tic": 100,
                    "aperture": aperture,
                    "tmag": tmag,
                    "truth_period_d": 0.5,
                    "truth_model_depth": 0.1,
                    "truth_radius_rearth": 4.0,
                    "delta_signal_multi_snr_mad": snr,
                    "depth_retention_frac": 0.5,
                    "n_good_ap_in_transit": 10,
                }
            )
    survival_path = tmp_path / "survival.csv"
    pd.DataFrame(rows).to_csv(survival_path, index=False)

    summary = module.summarize(
        survival_path,
        tmp_path / "out",
        current_aperture="DET_FLUX_ADP",
        recoverable_snr=7.0,
        best_group="adp",
    )

    assert summary["current_n_snr_ge_threshold"] == 1
    assert summary["best_group_n_snr_ge_threshold"] == 2
    assert summary["newly_recoverable_by_best_group"] == 1
    assert summary["best_single_aperture"]["aperture"] == "DET_FLUX_ADP_SML"
    assert (tmp_path / "out" / "newly_recoverable_by_best_aperture.csv").exists()
    joined = pd.read_csv(tmp_path / "out" / "best_aperture_per_injection.csv")
    assert joined["newly_recoverable_by_aperture"].sum() == 1


def test_predetrend_debug_queue_builder_merges_priority_classes(tmp_path) -> None:
    module = _load_predetrend_debug_queue_builder()
    failure = pd.DataFrame(
        {
            "injection_id": ["inj_a", "inj_b", "inj_c"],
            "review_id": ["a", "b", "c"],
            "tic": [1, 2, 3],
            "failure_class": [
                "snr_qualified_bls_ranking_loss",
                "low_snr_unmatched",
                "snr_qualified_bls_ranking_loss",
            ],
            "topn_recovery_status": ["bls_peak_mismatch", "bls_peak_mismatch", "bls_peak_mismatch"],
            "leo_report_name": ["a.pdf", "b.pdf", "c.pdf"],
        }
    )
    aperture = pd.DataFrame(
        {
            "injection_id": ["inj_a", "inj_b", "inj_c"],
            "current_snr": [8.0, 2.0, 8.5],
            "best_aperture": ["DET_FLUX_ADP", "DET_FLUX_ADP_SML", "DET_FLUX_ADP_SML"],
            "best_snr": [9.0, 10.0, 11.0],
            "best_depth_retention": [0.5, 0.6, 0.7],
            "newly_recoverable_by_aperture": [False, True, True],
            "snr_gain": [1.0, 8.0, 2.5],
        }
    )
    failure_path = tmp_path / "failure.csv"
    aperture_path = tmp_path / "aperture.csv"
    out_dir = tmp_path / "out"
    failure.to_csv(failure_path, index=False)
    aperture.to_csv(aperture_path, index=False)

    summary = module.build_priority_queues(failure_path, aperture_path, out_dir, min_best_snr=7.0)

    assert summary["n_snr_qualified_bls_ranking_loss"] == 2
    assert summary["n_newly_recoverable_by_best_aperture"] == 2
    assert summary["n_priority_unique"] == 3
    priority = pd.read_csv(out_dir / "priority_debug_queue.csv")
    assert set(priority["injection_id"]) == {"inj_a", "inj_b", "inj_c"}
    both = priority.loc[priority["injection_id"].eq("inj_c"), "priority_reason"].item()
    assert both == "search_ranking_loss+newly_recoverable_by_best_aperture"
    assert (out_dir / "priority_injection_ids.txt").exists()


def test_priority_sweep_summary_reports_delta_vs_baseline(tmp_path) -> None:
    module = _load_priority_sweep_summary()
    priority = pd.DataFrame(
        {
            "injection_id": ["inj_a", "inj_b", "inj_c"],
            "priority_reason": [
                "search_ranking_loss",
                "newly_recoverable_by_best_aperture",
                "newly_recoverable_by_best_aperture",
            ],
            "recommended_aperture": ["DET_FLUX_ADP", "DET_FLUX_ADP_SML", "DET_FLUX_SML"],
            "failure_class": [
                "snr_qualified_bls_ranking_loss",
                "low_snr_unmatched",
                "low_snr_unmatched",
            ],
        }
    )
    sweep_dir = tmp_path / "sweep"
    sweep_dir.mkdir()
    overview = pd.DataFrame(
        {
            "name": ["adp_priority_5k", "adp_sml_priority_5k"],
            "apertures": ["DET_FLUX_ADP", "DET_FLUX_ADP_SML"],
            "n_periods": [5000, 5000],
            "durations_min": ["3+5", "3+5"],
            "queue_csv": [
                str(sweep_dir / "adp_priority_5k" / "review_queue.csv"),
                str(sweep_dir / "adp_sml_priority_5k" / "review_queue.csv"),
            ],
        }
    )
    overview.to_csv(sweep_dir / "sweep_overview.csv", index=False)
    (sweep_dir / "adp_priority_5k").mkdir()
    (sweep_dir / "adp_sml_priority_5k").mkdir()
    pd.DataFrame(
        {
            "injection_id": ["inj_a", "inj_b", "inj_c"],
            "recovery_status": ["bls_peak_mismatch", "bls_peak_mismatch", "bls_recovered"],
            "topn_recovery_status": ["bls_peak_mismatch", "bls_peak_mismatch", "bls_top1_recovered"],
        }
    ).to_csv(sweep_dir / "adp_priority_5k" / "review_queue.csv", index=False)
    pd.DataFrame(
        {
            "injection_id": ["inj_a", "inj_b", "inj_c"],
            "recovery_status": ["bls_recovered", "bls_recovered", "bls_recovered"],
            "topn_recovery_status": ["bls_top1_recovered", "bls_top1_recovered", "bls_top1_recovered"],
        }
    ).to_csv(sweep_dir / "adp_sml_priority_5k" / "review_queue.csv", index=False)
    priority_path = tmp_path / "priority.csv"
    priority.to_csv(priority_path, index=False)

    summary = module.summarize_priority_sweep(
        sweep_dir,
        priority_path,
        tmp_path / "out",
        baseline_sweep="adp_priority_5k",
    )

    assert summary["baseline_found"] is True
    assert summary["best_sweep"]["name"] == "adp_sml_priority_5k"
    deltas = pd.read_csv(tmp_path / "out" / "priority_sweep_deltas_vs_baseline.csv")
    small = deltas[deltas["name"].eq("adp_sml_priority_5k")].iloc[0]
    assert small["new_matches_vs_baseline"] == 2
    assert small["net_match_gain_vs_baseline"] == 2
    assert (tmp_path / "out" / "adp_sml_priority_5k_new_matches_vs_adp_priority_5k.csv").exists()


def test_priority_sweep_preflight_validates_h5_and_csv_inputs(tmp_path) -> None:
    module = _load_priority_sweep_preflight()
    injection_h5 = tmp_path / "injections.h5"
    ids = ["predet_000001", "predet_000002"]
    apertures = ["DET_FLUX_ADP", "DET_FLUX_ADP_SML"]
    with h5py.File(injection_h5, "w") as h5:
        root = h5.create_group("injections")
        for inj_id in ids:
            group = root.create_group(inj_id)
            group.attrs["injection_id"] = inj_id
            group.attrs["apertures"] = json.dumps(apertures)
            group.create_dataset("time", data=np.arange(5, dtype=float))
            group.create_dataset("quality", data=np.zeros(5, dtype=int))
            group.create_dataset("in_transit", data=np.array([False, True, False, True, False]))
            for aperture in apertures:
                group.create_dataset(f"{aperture}_injected", data=np.ones(5))

    id_file = tmp_path / "ids.txt"
    id_file.write_text("\n".join(ids) + "\n")
    priority_csv = tmp_path / "priority.csv"
    pd.DataFrame(
        {
            "injection_id": ids,
            "priority_reason": ["search_ranking_loss", "newly_recoverable_by_best_aperture"],
            "recommended_aperture": apertures,
            "failure_class": ["snr_qualified_bls_ranking_loss", "low_snr_unmatched"],
        }
    ).to_csv(priority_csv, index=False)
    survival_rows = []
    for inj_id in ids:
        for aperture in apertures:
            survival_rows.append(
                {
                    "injection_id": inj_id,
                    "aperture": aperture,
                    "delta_signal_multi_snr_mad": 8.0,
                }
            )
    survival_csv = tmp_path / "survival.csv"
    pd.DataFrame(survival_rows).to_csv(survival_csv, index=False)

    summary = module.preflight(
        injection_h5=injection_h5,
        injection_id_file=id_file,
        priority_queue_csv=priority_csv,
        survival_csv=survival_csv,
        required_apertures=tuple(apertures),
        expected_id_count=2,
    )

    assert summary["ok"] is True
    assert summary["n_injection_ids"] == 2
    assert summary["h5"]["n_groups"] == 2
    assert summary["h5"]["available_aperture_counts"]["DET_FLUX_ADP"] == 2
    assert summary["survival_csv"]["n_priority_ids"] == 2


def test_predetrend_detrending_method_sweep_measures_signal_and_trend(tmp_path) -> None:
    module = _load_detrending_method_sweep()
    injection_h5 = tmp_path / "injections.h5"
    n = 240
    time = np.linspace(0.0, 3.0, n)
    quality = np.zeros(n, dtype=int)
    in_transit = np.abs(time - 1.5) < 0.04
    trend = 1000.0 + 80.0 * np.sin(2 * np.pi * time / 2.2)
    original = trend.copy()
    injected = original.copy()
    injected[in_transit] -= 100.0
    with h5py.File(injection_h5, "w") as h5:
        group = h5.create_group("injections/predet_000000")
        group.attrs["tic"] = 123
        group.attrs["tessmag"] = 18.0
        group.attrs["period_d"] = 1.0
        group.attrs["duration_min"] = 10.0
        group.attrs["depth"] = 0.1
        group.attrs["model_depth"] = 0.1
        group.attrs["sampled_model_depth"] = 0.1
        group.create_dataset("time", data=time)
        group.create_dataset("quality", data=quality)
        group.create_dataset("in_transit", data=in_transit)
        group.create_dataset("RAW_FLUX_Small_original", data=original)
        group.create_dataset("RAW_FLUX_Small_injected", data=injected)

    summary = module.run_sweep(
        injection_h5=injection_h5,
        out_dir=tmp_path / "out",
        methods=("constant", "poly2_gap05"),
        raw_apertures=("Small",),
        injection_ids=None,
        limit=0,
        trend_bin_d=0.5,
        snr_threshold=7.0,
    )

    assert summary["n_injection_ids"] == 1
    detail = pd.read_csv(tmp_path / "out" / "detrending_method_detail.csv")
    assert set(detail["method"]) == {"constant", "poly2_gap05"}
    assert np.isfinite(detail["depth_retention_frac"]).all()
    assert (tmp_path / "out" / "detrending_method_summary.csv").exists()
    assert (tmp_path / "out" / "summary.json").exists()


def test_leo_plot_copy_disables_pathological_fit_metrics() -> None:
    module = _load_queue_builder()

    class DummyTLC:
        pass

    tlc = DummyTLC()
    tlc.tic = 123
    tlc.planetno = 1
    tlc.time = np.linspace(2459825.0, 2459852.0, 1000)
    tlc.raw = np.ones_like(tlc.time)
    tlc.flux = np.ones_like(tlc.time)
    tlc.flux[500] = 0.8
    tlc.flux_err = np.full_like(tlc.time, 0.01)
    tlc.per = 1.0
    tlc.epo = tlc.time[500]
    tlc.dur = 10.0 / 1440.0
    tlc.dep = 0.2
    tlc.metrics = {
        "per": tlc.per,
        "epo": tlc.epo,
        "dur": tlc.dur,
        "dep": tlc.dep,
        "zpt": 1.0,
        "transit_aic": 1.0,
        "transit_per": 1.0e9,
        "transit_epo": tlc.epo,
        "transit_dur": 1.0e-9,
        "transit_RpRs": np.nan,
        "transit_aRs": np.nan,
        "transit_b": np.nan,
        "trap_aic": 1.0,
        "trap_per": 1.0e9,
        "trap_epo": tlc.epo,
        "trap_dur": 1.0e-9,
        "trap_dep": np.nan,
        "trap_qtran": np.nan,
        "trap_qin": np.nan,
        "trap_zpt": np.nan,
        "phs_sec": 0.5,
        "dep_sec": np.nan,
        "sig_sec": np.nan,
    }

    plot_tlc, notes = module._leo_plot_copy(tlc)

    assert plot_tlc is not tlc
    assert tlc.metrics["transit_aic"] == 1.0
    assert np.isnan(plot_tlc.metrics["transit_aic"])
    assert np.isnan(plot_tlc.metrics["trap_aic"])
    assert np.isnan(plot_tlc.metrics["phs_sec"])
    assert "transit_model_disabled" in notes
    assert "secondary_disabled" in notes
    assert np.any(plot_tlc.near_tran)
    assert np.any(~plot_tlc.near_tran)
