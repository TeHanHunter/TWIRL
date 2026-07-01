from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_verifier():
    path = REPO_ROOT / "scripts" / "stage5_validation" / "verify_ranker_selected_ephemerides.py"
    spec = importlib.util.spec_from_file_location("verify_ranker_selected_ephemerides", path)
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
                "aperture": "DET_FLUX_ADP_SML",
                "period_d": 0.8,
                "t0_bjd": 2459825.1,
                "duration_min": 8.0,
                "ranker_p_background_peak": 0.05,
                "ranker_p_signal_peak": 0.95,
            },
            {
                "tic": 101,
                "ranker_selection_rank": 2,
                "aperture": "DET_FLUX_SML",
                "period_d": 1.6,
                "t0_bjd": 2459825.2,
                "duration_min": 12.0,
                "ranker_p_background_peak": 0.20,
                "ranker_p_signal_peak": 0.80,
            },
            {
                "tic": 202,
                "ranker_selection_rank": 1,
                "aperture": "DET_FLUX_ADP_SML",
                "period_d": 3.5,
                "t0_bjd": 2459825.3,
                "duration_min": 25.0,
                "ranker_p_background_peak": 0.35,
                "ranker_p_signal_peak": 0.65,
            },
        ]
    )


def test_ranker_selected_ephemerides_verifier_passes_valid_table(tmp_path: Path) -> None:
    module = _load_verifier()
    selected_path = tmp_path / "selected_ephemerides.csv"
    _selected_rows().to_csv(selected_path, index=False)

    result = module.verify_selected_ephemerides(
        selected_ephemerides=selected_path,
        top_n=2,
        min_rows=3,
        min_targets=2,
    )

    assert result["passed"] is True
    assert result["n_rows"] == 3
    assert result["n_targets"] == 2
    assert result["score_column"] == "ranker_p_signal_peak"
    assert result["rank_counts"] == {"1": 2, "2": 1}


def test_ranker_selected_ephemerides_verifier_rejects_bad_real_rows(tmp_path: Path) -> None:
    module = _load_verifier()
    selected = _selected_rows()
    selected.loc[0, "ranker_selection_rank"] = 3
    selected.loc[1, "period_d"] = -1.0
    selected.loc[2, "injection_id"] = "predet_000002"
    selected_path = tmp_path / "selected_ephemerides.csv"
    selected.to_csv(selected_path, index=False)

    result = module.verify_selected_ephemerides(
        selected_ephemerides=selected_path,
        top_n=2,
        min_rows=3,
        min_targets=2,
    )

    assert result["passed"] is False
    joined = "\n".join(result["failures"])
    assert "invalid ranker_selection_rank" in joined
    assert "non-positive/non-finite period_d" in joined
    assert "non-empty injection_id" in joined
