from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_audit():
    path = REPO_ROOT / "scripts" / "stage5_validation" / "audit_injection_bls_failure_modes.py"
    spec = importlib.util.spec_from_file_location("audit_injection_bls_failure_modes", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _synthetic_peaks() -> pd.DataFrame:
    rows = []
    specs = [
        ("top1", 0.2, 8.0, 0.20, 20, 1),
        ("ranker", 0.4, 8.0, 0.18, 18, 3),
        ("hardmiss", 0.5, 8.0, 0.25, 25, None),
        ("lowcad", 0.6, 8.0, 0.30, 1, None),
    ]
    for injection_id, period, radius, depth, n_good, signal_rank in specs:
        for rank in (1, 2, 3):
            is_signal = signal_rank == rank
            rows.append(
                {
                    "injection_id": injection_id,
                    "tic": hash(injection_id) % 100000,
                    "tmag": 17.2,
                    "aperture": "DET_FLUX_ADP_SML",
                    "is_candidate_peak": True,
                    "peak_rank": rank,
                    "period_d": period if is_signal else 0.12 + rank * 0.02,
                    "duration_min": 8.0,
                    "depth": depth if is_signal else 0.01,
                    "depth_snr": 10.0 if is_signal else 2.0 + rank,
                    "sde": 50.0 if is_signal else 20.0 - rank,
                    "match_kind": "exact" if is_signal else "mismatch",
                    "is_injected_signal_peak": is_signal,
                    "exact_ephemeris_match": is_signal,
                    "harmonic_ephemeris_match": False,
                    "truth_period_d": period,
                    "truth_radius_rearth": radius,
                    "truth_duration_min": 8.0,
                    "truth_sampled_model_depth": depth,
                    "truth_n_good_in_transit": n_good,
                    "period_ratio": 1.0 if is_signal else (0.12 + rank * 0.02) / period,
                }
            )
    return pd.DataFrame(rows)


def test_failure_mode_audit_separates_ranker_fixable_and_hard_misses() -> None:
    module = _load_audit()
    table = module.build_failure_mode_table(
        _synthetic_peaks(),
        top_k=3,
        period_bins=(0.0, 1.0),
        radius_bins=(0.0, 12.0),
        tmag_bins=(0.0, 18.0),
        observability_quantile=0.25,
        min_reference_count=1,
    )
    modes = dict(zip(table["injection_id"], table["failure_mode"]))

    assert modes["top1"] == "top1_recovered"
    assert modes["ranker"] == "ranker_fixable"
    assert modes["hardmiss"] == "high_observability_not_in_top3"
    assert modes["lowcad"] == "low_cadence_not_in_top3"
    assert table.loc[table["injection_id"].eq("hardmiss"), "high_observability_not_in_topk"].iloc[0]


def test_failure_mode_summary_counts() -> None:
    module = _load_audit()
    table = module.build_failure_mode_table(
        _synthetic_peaks(),
        top_k=3,
        period_bins=(0.0, 1.0),
        radius_bins=(0.0, 12.0),
        tmag_bins=(0.0, 18.0),
        observability_quantile=0.25,
        min_reference_count=1,
    )
    summary = module.summarize_failure_modes(table, top_k=3)
    cells = module.build_cell_failure_table(table, top_k=3)

    assert summary["n_injections"] == 4
    assert summary["top1_n"] == 1
    assert summary["top3_n"] == 2
    assert summary["ranking_loss_n"] == 1
    assert summary["high_observability_not_in_topk_n"] == 1
    assert cells["high_observability_not_in_topk_n"].sum() == 1
