from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/stage5_validation/plot_s56_harmonic_cnn_performance.py"


def _load_plot_module():
    spec = importlib.util.spec_from_file_location("plot_s56_harmonic_cnn_performance", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_performance_plot_excludes_preserve_only_morphology_rows() -> None:
    module = _load_plot_module()
    frame = pd.DataFrame(
        {
            "morphology_target_index": [0, 1, 2, 3, -1],
            "morphology_prediction_index": [0, 1, 2, 3, 0],
            "p_planet_like": [0.8, 0.1, 0.1, 0.1, 0.8],
            "p_eclipse_contact": [0.1, 0.7, 0.1, 0.1, 0.1],
            "p_smooth_variable": [0.05, 0.1, 0.7, 0.1, 0.05],
            "p_other": [0.05, 0.1, 0.1, 0.7, 0.05],
        }
    )

    confusion = module._confusion_matrix(frame)
    calibration, _ = module._calibration_bins(frame)

    assert confusion.sum() == 4
    assert confusion.trace() == 4
    assert calibration["n"].sum() == 4
