from __future__ import annotations

import json

import numpy as np
import pandas as pd

from twirl.search.peak_ranker import (
    PeakRankerConfig,
    prepare_peak_training_frame,
    score_peak_table,
    select_ranked_ephemerides,
    train_peak_ranker,
    write_peak_ranker_outputs,
)


def _synthetic_peak_table(n_injections: int = 80, n_peaks: int = 5) -> pd.DataFrame:
    rows = []
    for i in range(n_injections):
        injection_id = f"inj_{i:04d}"
        rankable = i < n_injections - 8
        true_rank = 4
        for peak_rank in range(1, n_peaks + 1):
            is_signal = rankable and peak_rank == true_rank
            rows.append(
                {
                    "injection_id": injection_id,
                    "tic": 100000 + i,
                    "tmag": 17.0 + (i % 5) * 0.3,
                    "aperture": "DET_FLUX_SML" if peak_rank % 2 else "DET_FLUX_ADP_SML",
                    "is_candidate_peak": True,
                    "peak_rank": peak_rank,
                    "period_d": 1.4 if is_signal else 0.2 + 0.1 * peak_rank,
                    "t0_bjd": 2459000.0,
                    "duration_min": 10.0 if is_signal else 25.0,
                    "qtran": 10.0 / 1440.0 / 1.4 if is_signal else 0.1,
                    "depth": 0.2 if is_signal else 0.02,
                    "depth_snr": 18.0 if is_signal else 2.0 + peak_rank * 0.2,
                    "sde": 50.0 if is_signal else 20.0 - peak_rank,
                    "log_power": 4.0 if is_signal else 1.0,
                    "n_cad_total": 11000,
                    "n_cad_kept": 10000,
                    "dropout_frac": 0.1,
                    "baseline_d": 27.0,
                    "n_orbits": 2,
                    "is_injected_signal_peak": is_signal,
                    "match_kind": "exact" if is_signal else "mismatch",
                    "truth_period_d": 1.4,
                    "truth_t0_bjd": 2459000.0,
                    "truth_radius_rearth": 10.0,
                }
            )
    return pd.DataFrame(rows)


def test_peak_ranker_reranks_signal_above_bls_rank() -> None:
    peaks = _synthetic_peak_table()
    result = train_peak_ranker(
        peaks,
        PeakRankerConfig(random_state=3, validation_fraction=0.2, test_fraction=0.2, max_iter=200),
    )

    all_summary = result.summary["splits"]["all"]
    bls_top1 = all_summary["bls_rank_recall_at_k"][1]
    model_top1 = all_summary["model_recall_at_k"][1]
    model_top5 = all_summary["model_recall_at_k"][5]

    assert result.summary["n_injections"] == 80
    assert result.summary["n_rankable_injections"] == 72
    assert result.model.metadata["n_train_background_only_injections"] > 0
    assert bls_top1["n"] == 0
    assert model_top1["rankable_fraction"] > 0.90
    assert model_top5["rankable_fraction"] == 1.0
    assert result.selected_peaks["ranker_p_signal_peak"].notna().all()


def test_prepare_peak_training_frame_marks_rankable_groups() -> None:
    frame = prepare_peak_training_frame(_synthetic_peak_table(n_injections=4, n_peaks=3))

    assert "rank_inv" in frame.columns
    assert "sde_group_rank_inv" in frame.columns
    assert frame.groupby("injection_id")["sde_group_rank_inv"].max().eq(1.0).all()
    assert set(frame["peak_label"]) == {"background_peak"}
    grouped = frame.groupby("injection_id")["rankable_injection"].first()
    assert grouped.sum() == 0  # true rank 4 is outside the 3-peak synthetic table


def test_peak_ranker_outputs_are_written(tmp_path) -> None:
    peaks = _synthetic_peak_table()
    result = train_peak_ranker(peaks, PeakRankerConfig(random_state=4, max_iter=200))
    write_peak_ranker_outputs(result, tmp_path, write_scored=True)

    assert (tmp_path / "peak_ranker_model.npz").exists()
    assert (tmp_path / "scored_peaks.csv").exists()
    assert (tmp_path / "selected_ephemerides.csv").exists()
    assert (tmp_path / "coefficients.csv").exists()
    summary = json.loads((tmp_path / "summary.json").read_text())
    assert summary["n_rankable_injections"] == 72
    coeffs = pd.read_csv(tmp_path / "coefficients.csv")
    assert np.isfinite(coeffs["coefficient"]).all()


def test_peak_ranker_scores_real_style_peaks_without_truth_columns() -> None:
    peaks = _synthetic_peak_table()
    result = train_peak_ranker(peaks, PeakRankerConfig(random_state=5, max_iter=200))

    real_peaks = (
        peaks.drop(
            columns=[
                "injection_id",
                "is_candidate_peak",
                "is_injected_signal_peak",
                "match_kind",
                "truth_period_d",
                "truth_t0_bjd",
                "truth_radius_rearth",
            ]
        )
        .rename(columns={"tic": "tic"})
        .copy()
    )
    # Simulate an older real candidate table before cadence-cleaning breakdown
    # columns were added.
    real_peaks = real_peaks.drop(
        columns=["n_cad_quality", "n_cad_edge_trimmed", "n_cad_sigma_clipped", "quality_dropout_frac"],
        errors="ignore",
    )
    real_peaks["status"] = "ok"

    scored = score_peak_table(real_peaks, result.model, group_column="tic")
    selected = select_ranked_ephemerides(scored, id_column="tic", top_n=2)

    assert "ranker_p_signal_peak" in scored.columns
    assert "n_cad_quality" in scored.columns
    assert "n_cad_edge_trimmed" in scored.columns
    assert "sde_group_rank_inv" in scored.columns
    assert selected.groupby("tic").size().max() == 2
    assert selected["ranker_selection_rank"].isin([1, 2]).all()
