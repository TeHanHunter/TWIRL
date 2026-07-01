from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_peak_gate():
    path = REPO_ROOT / "scripts" / "stage5_validation" / "summarize_injection_peak_gate.py"
    spec = importlib.util.spec_from_file_location("summarize_injection_peak_gate", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _synthetic_peaks() -> pd.DataFrame:
    rows = []
    specs = [
        ("a", 16.5, 0.2, 10.0, 1),
        ("b", 17.5, 1.0, 8.0, 3),
        ("c", 19.5, 8.0, 2.0, None),
    ]
    for injection_id, tmag, period, radius, signal_rank in specs:
        for rank in (1, 2, 3):
            is_signal = signal_rank == rank
            rows.append(
                {
                    "injection_id": injection_id,
                    "tic": 1000 + ord(injection_id),
                    "tmag": tmag,
                    "truth_period_d": period,
                    "truth_radius_rearth": radius,
                    "truth_duration_min": 10.0,
                    "is_candidate_peak": True,
                    "peak_rank": rank,
                    "is_injected_signal_peak": is_signal,
                    "exact_ephemeris_match": is_signal,
                    "harmonic_ephemeris_match": False,
                    "match_kind": "exact" if is_signal else "mismatch",
                }
            )
    return pd.DataFrame(rows)


def test_peak_gate_separates_top1_ranking_loss_and_miss() -> None:
    module = _load_peak_gate()
    injections = module.build_injection_level_table(_synthetic_peaks(), top_k=3)

    assert len(injections) == 3
    assert injections["top1_match"].sum() == 1
    assert injections["top3_match"].sum() == 2
    assert injections["ranking_loss"].sum() == 1
    assert injections["not_in_top3"].sum() == 1


def test_peak_gate_builds_empirical_recovery_cells() -> None:
    module = _load_peak_gate()
    injections = module.build_injection_level_table(_synthetic_peaks(), top_k=3)
    injections = module.add_recovery_bins(
        injections,
        period_bins=(0.0, 0.5, 2.0, 20.0),
        radius_bins=(0.0, 4.0, 12.0),
        tmag_bins=(0.0, 17.0, 18.0, 99.0),
    )
    cells = module.build_cell_recovery_table(
        injections,
        top_k=3,
        min_cell_count=1,
        recovery_threshold=0.5,
    )
    annotated = module.annotate_injections_with_cells(injections, cells, top_k=3)
    summary = module.summarize_gate(
        injections,
        cells,
        top_k=3,
        min_topk_fraction_for_ranker=0.5,
        min_ranking_loss_fraction=0.1,
    )

    assert set(cells["recovery_gate"]) == {"above_empirical_threshold", "below_empirical_threshold"}
    assert annotated["recovery_cell_top3_fraction"].notna().all()
    assert summary["top1"]["n"] == 1
    assert summary["top3"]["n"] == 2
    assert summary["ranking_loss"]["n"] == 1
    assert summary["recommendation"] == "train_peak_ranker_before_leo_human_review"
