from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_summarizer():
    path = REPO_ROOT / "scripts" / "stage5_validation" / "summarize_ranker_selected_real_queue.py"
    spec = importlib.util.spec_from_file_location("summarize_ranker_selected_real_queue", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _selected_rows() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "tic": 101,
                "ranker_selection_rank": 1,
                "tmag": 16.5,
                "aperture": "DET_FLUX_ADP_SML",
                "peak_rank": 3,
                "period_d": 0.8,
                "t0_bjd": 2459825.1,
                "duration_min": 8.0,
                "depth": 0.08,
                "depth_snr": 8.5,
                "sde": 30.0,
                "ranker_p_background_peak": 0.05,
                "ranker_p_signal_peak": 0.95,
            },
            {
                "tic": 101,
                "ranker_selection_rank": 2,
                "tmag": 16.5,
                "aperture": "DET_FLUX_SML",
                "peak_rank": 1,
                "period_d": 1.6,
                "t0_bjd": 2459825.2,
                "duration_min": 12.0,
                "depth": 0.06,
                "depth_snr": 6.5,
                "sde": 26.0,
                "ranker_p_background_peak": 0.20,
                "ranker_p_signal_peak": 0.80,
            },
            {
                "tic": 202,
                "ranker_selection_rank": 1,
                "tmag": 18.4,
                "aperture": "DET_FLUX_ADP_SML",
                "peak_rank": 4,
                "period_d": 3.5,
                "t0_bjd": 2459825.3,
                "duration_min": 25.0,
                "depth": 0.12,
                "depth_snr": 5.5,
                "sde": 21.0,
                "ranker_p_background_peak": 0.35,
                "ranker_p_signal_peak": 0.65,
            },
        ]
    )


def _review_rows() -> pd.DataFrame:
    rows = _selected_rows().copy()
    rows["review_id"] = [
        "real:101:ranker:1",
        "real:101:ranker:2",
        "real:202:ranker:1",
    ]
    rows["source_kind"] = "real_candidate"
    rows["rep_aperture"] = rows["aperture"]
    return rows


def test_ranker_selected_summary_uses_review_queue_and_verification(tmp_path: Path) -> None:
    module = _load_summarizer()
    selected_path = tmp_path / "selected_ephemerides.csv"
    review_path = tmp_path / "review_queue.csv"
    verification_path = tmp_path / "verification.json"
    out_dir = tmp_path / "summary"
    _selected_rows().to_csv(selected_path, index=False)
    _review_rows().to_csv(review_path, index=False)
    verification_path.write_text(json.dumps({"passed": True, "failures": []}))

    summary = module.summarize_ranker_selection(
        selected_ephemerides=selected_path,
        review_queue=review_path,
        verification_json=verification_path,
        out_dir=out_dir,
    )

    assert summary["score_column"] == "ranker_p_signal_peak"
    assert summary["n_selected_rows"] == 3
    assert summary["n_selected_targets"] == 2
    assert summary["n_review_rows"] == 3
    assert summary["verification_passed"] is True
    assert summary["rows_per_target"] == {"1": 1, "2": 1}
    assert summary["selection_rank_counts"] == {"1": 2, "2": 1}
    assert (out_dir / "summary.json").exists()
    assert (out_dir / "summary.md").exists()
    assert (out_dir / "by_aperture.csv").exists()
    top = pd.read_csv(out_dir / "top_ranker_rows.csv")
    assert list(top["tic"].head(2)) == [101, 101]


def test_ranker_selected_summary_works_before_review_queue_exists(tmp_path: Path) -> None:
    module = _load_summarizer()
    selected_path = tmp_path / "selected_ephemerides.csv"
    out_dir = tmp_path / "summary"
    _selected_rows().to_csv(selected_path, index=False)

    summary = module.summarize_ranker_selection(
        selected_ephemerides=selected_path,
        review_queue=tmp_path / "missing_review_queue.csv",
        verification_json=None,
        out_dir=out_dir,
    )

    assert summary["n_selected_rows"] == 3
    assert summary["n_review_rows"] == 0
    assert summary["verification_passed"] is None
    period_bins = pd.read_csv(out_dir / "by_period_bin.csv")
    assert set(period_bins["period_bin"]) >= {"0.5-1", "1-2", "2-5"}
