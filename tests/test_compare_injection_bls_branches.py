from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_compare():
    path = REPO_ROOT / "scripts" / "stage5_validation" / "compare_injection_bls_branches.py"
    spec = importlib.util.spec_from_file_location("compare_injection_bls_branches", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _peaks(branch: str, ranks_by_id: dict[str, int | None]) -> pd.DataFrame:
    rows = []
    for idx, injection_id in enumerate(("both", "rerank", "rescue", "degrade", "miss")):
        signal_rank = ranks_by_id[injection_id]
        for rank in (1, 2, 3):
            is_signal = signal_rank == rank
            rows.append(
                {
                    "search_branch": branch,
                    "injection_id": injection_id,
                    "tic": 1000 + idx,
                    "tmag": 16.5 + idx,
                    "truth_period_d": 0.2 + 0.1 * idx,
                    "truth_radius_rearth": 4.0 + idx,
                    "truth_duration_min": 8.0,
                    "is_candidate_peak": True,
                    "peak_rank": rank,
                    "is_injected_signal_peak": is_signal,
                    "exact_ephemeris_match": is_signal,
                    "harmonic_ephemeris_match": False,
                    "match_kind": "exact" if is_signal else "mismatch",
                }
            )
    return pd.DataFrame(rows)


def test_compare_branches_reports_rescues_and_degradations() -> None:
    module = _load_compare()
    baseline = _peaks(
        "standard",
        {
            "both": 1,
            "rerank": 3,
            "rescue": None,
            "degrade": 1,
            "miss": None,
        },
    )
    branch = _peaks(
        "short_pmax2",
        {
            "both": 1,
            "rerank": 1,
            "rescue": 2,
            "degrade": None,
            "miss": None,
        },
    )

    comparison, cells, summary = module.compare_branches(
        baseline,
        branch,
        top_k=3,
        baseline_label="standard",
        branch_label="short",
        baseline_filter="standard",
        branch_filter="short_pmax2",
        period_bins=(0.0, 1.0),
        radius_bins=(0.0, 20.0),
        tmag_bins=(0.0, 99.0),
    )
    by_id = comparison.set_index("injection_id")

    assert summary["standard_top1"]["n"] == 2
    assert summary["standard_top3"]["n"] == 3
    assert summary["short_top1"]["n"] == 2
    assert summary["short_top3"]["n"] == 3
    assert summary["union_top3"]["n"] == 4
    assert summary["branch_rescued_top3"]["n"] == 1
    assert summary["branch_reranked_to_top1"]["n"] == 1
    assert summary["branch_degraded_top3"]["n"] == 1
    assert by_id.loc["rescue", "branch_comparison_status"] == "branch_rescued_top3"
    assert by_id.loc["rerank", "branch_comparison_status"] == "branch_reranked_to_top1"
    assert by_id.loc["degrade", "branch_comparison_status"] == "branch_degraded_top1"
    assert by_id.loc["miss", "branch_comparison_status"] == "both_missed"
    assert not cells.empty


def test_compare_branches_defaults_to_overlap_for_subset_branch() -> None:
    module = _load_compare()
    baseline = _peaks(
        "standard",
        {
            "both": 1,
            "rerank": 3,
            "rescue": None,
            "degrade": 1,
            "miss": None,
        },
    )
    branch_full = _peaks(
        "short_pmax2",
        {
            "both": 1,
            "rerank": 1,
            "rescue": 2,
            "degrade": None,
            "miss": None,
        },
    )
    branch_subset = branch_full.loc[branch_full["injection_id"].isin(["rescue", "miss"])].copy()

    comparison, _cells, summary = module.compare_branches(
        baseline,
        branch_subset,
        top_k=3,
        baseline_label="standard",
        branch_label="short",
        baseline_filter="standard",
        branch_filter="short_pmax2",
        period_bins=(0.0, 1.0),
        radius_bins=(0.0, 20.0),
        tmag_bins=(0.0, 99.0),
    )

    assert set(comparison["injection_id"]) == {"rescue", "miss"}
    assert summary["comparison_scope"] == "intersection"
    assert summary["n_baseline_injections"] == 5
    assert summary["n_branch_injections"] == 2
    assert summary["n_overlap_injections"] == 2
    assert summary["branch_rescued_top3"]["n"] == 1
    assert summary["branch_degraded_top3"]["n"] == 0
