from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from twirl.plotting.teacher_ssl_fullpool_evaluation import (
    MODEL_ORDER,
    POLICY_LABELS,
    plot_development_performance,
    plot_reference_umap,
)
from twirl.vetting.harmonic_cnn import MORPHOLOGY_CLASSES


def test_fullpool_performance_and_umap_figures_render(tmp_path: Path) -> None:
    performance_records: list[dict[str, object]] = []
    class_records: list[dict[str, object]] = []
    for policy_index, policy in enumerate(POLICY_LABELS):
        for model_index, model in enumerate(MODEL_ORDER):
            value = 0.55 + 0.05 * model_index + 0.02 * policy_index
            record: dict[str, object] = {
                "model": model,
                "label_policy": policy,
                "n_rows": 100,
                "n_tics": 95,
                "accuracy": value,
                "expected_calibration_error": 0.05,
            }
            for metric in (
                "balanced_accuracy",
                "macro_f1",
                "macro_average_precision",
            ):
                record[metric] = value
                record[f"{metric}_low"] = value - 0.04
                record[f"{metric}_high"] = value + 0.04
            record["expected_calibration_error_low"] = 0.03
            record["expected_calibration_error_high"] = 0.08
            performance_records.append(record)
            for class_index, label in enumerate(MORPHOLOGY_CLASSES):
                class_records.append(
                    {
                        "model": model,
                        "label_policy": policy,
                        "class": label,
                        "support": 10 + class_index,
                        "precision": value,
                        "recall": value,
                        "f1": value,
                        "average_precision": min(
                            0.95, value + 0.03 * class_index
                        ),
                    }
                )
    performance_path = tmp_path / "performance.csv"
    per_class_path = tmp_path / "per_class.csv"
    pd.DataFrame(performance_records).to_csv(performance_path, index=False)
    pd.DataFrame(class_records).to_csv(per_class_path, index=False)

    performance_artifacts = plot_development_performance(
        performance_path=performance_path,
        per_class_path=per_class_path,
        output_dir=tmp_path,
    )

    assert Path(performance_artifacts["png"]).is_file()
    assert Path(performance_artifacts["pdf"]).is_file()
    rng = np.random.default_rng(56)
    labels = [
        "uncertain",
        "instrumental_or_systematic",
        "stellar_variability",
        "wide_transit_like",
        "eclipsing_binary_or_pceb",
        "planet_like",
    ]
    coordinates = pd.DataFrame(
        {
            "human_label": [labels[index % len(labels)] for index in range(120)],
            "period_d": np.geomspace(0.05, 10.0, 120),
            "umap_1": rng.normal(size=120),
            "umap_2": rng.normal(size=120),
        }
    )
    umap_artifacts = plot_reference_umap(
        coordinates=coordinates,
        output_dir=tmp_path,
    )

    assert Path(umap_artifacts["png"]).is_file()
    assert Path(umap_artifacts["pdf"]).is_file()
