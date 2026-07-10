from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_summarizer():
    path = REPO_ROOT / "scripts" / "stage5_validation" / "summarize_human_vetting_labels.py"
    spec = importlib.util.spec_from_file_location("summarize_human_vetting_labels", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_summary_uses_label_csv_when_queue_has_empty_label_columns(tmp_path: Path) -> None:
    module = _load_summarizer()
    queue = pd.DataFrame(
        [
            {
                "row_id": 0,
                "tic": 1,
                "sector": 56,
                "period_d": 1.0,
                "t0_bjd": 2459825.0,
                "source_bucket": "review_candidate",
                "label": "",
                "label_source": "",
                "leo_class": "PC",
                "topn_recovery_status": "bls_top1_recovered",
                "truth_period_d": 1.0,
                "truth_model_depth": 0.1,
                "tmag": 17.2,
            },
            {
                "row_id": 1,
                "tic": 2,
                "sector": 56,
                "period_d": 2.0,
                "t0_bjd": 2459826.0,
                "source_bucket": "review_candidate",
                "label": "",
                "label_source": "",
                "leo_class": "FA",
                "topn_recovery_status": "bls_peak_mismatch",
                "truth_period_d": 2.0,
                "truth_model_depth": 0.2,
                "tmag": 19.2,
            },
        ]
    )
    queue_csv = tmp_path / "review_queue.csv"
    queue.to_csv(queue_csv, index=False)
    labels = pd.DataFrame(
        [
            {
                "row_id": 0,
                "candidate_key": "1|56|1.0|2459825.0|review_candidate",
                "label": "planet_like",
                "label_source": "human",
                "labeler": "tehan",
                "notes": "",
                "updated_utc": "2026-06-26T00:00:00+00:00",
            }
        ]
    )
    labels_csv = tmp_path / "human_labels_vetted.csv"
    labels.to_csv(labels_csv, index=False)

    out_dir = tmp_path / "summary"
    summary = module.summarize_labels(queue_csv, labels_csv, out_dir)
    saved = json.loads((out_dir / "summary.json").read_text())
    joined = pd.read_csv(out_dir / "labeled_queue_joined.csv")

    assert summary["n_labeled"] == 1
    assert saved["label_counts"] == {"planet_like": 1}
    assert joined.loc[0, "label"] == "planet_like"
    assert joined.loc[0, "queue_label"] != "planet_like"
