from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_selector():
    path = REPO_ROOT / "scripts" / "stage5_validation" / "select_next_human_labels.py"
    spec = importlib.util.spec_from_file_location("select_next_human_labels", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _training_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "row_id": 0,
                "is_labeled": True,
                "human_label": "planet_like",
                "topn_recovery_status": "bls_top1_recovered",
                "leo_class": "PC",
                "tmag": 16.1,
                "truth_period_d": 0.2,
                "truth_radius_rearth": 8.0,
                "sde_max": 30.0,
            },
            {
                "row_id": 1,
                "is_labeled": True,
                "human_label": "instrumental_or_systematic",
                "topn_recovery_status": "bls_peak_mismatch",
                "leo_class": "FA",
                "tmag": 19.4,
                "truth_period_d": 5.0,
                "truth_radius_rearth": 2.0,
                "sde_max": 8.0,
            },
            {
                "row_id": 2,
                "is_labeled": False,
                "human_label": "",
                "topn_recovery_status": "bls_top1_recovered",
                "leo_class": "PC",
                "tmag": 17.2,
                "truth_period_d": 0.3,
                "truth_radius_rearth": 10.0,
                "sde_max": 20.0,
            },
            {
                "row_id": 3,
                "is_labeled": False,
                "human_label": "",
                "topn_recovery_status": "bls_topn_harmonic_match",
                "leo_class": "FA",
                "tmag": 18.2,
                "truth_period_d": 1.5,
                "truth_radius_rearth": 5.0,
                "sde_max": 15.0,
            },
            {
                "row_id": 4,
                "is_labeled": False,
                "human_label": "",
                "topn_recovery_status": "bls_peak_mismatch",
                "leo_class": "FA",
                "tmag": 19.7,
                "truth_period_d": 8.0,
                "truth_radius_rearth": 1.0,
                "sde_max": 9.0,
            },
        ]
    )


def test_select_next_rows_prioritizes_rare_leo_and_bls_modes() -> None:
    module = _load_selector()
    selected, coverage, summary = module.select_next_rows(
        _training_table(),
        n_rows=2,
        target_per_cell=3,
    )

    assert selected.iloc[0]["row_id"] == 2
    assert "rare_leo:PC" in selected.iloc[0]["priority_reasons"]
    assert "bls_truth_mode:bls_top1_recovered" in selected.iloc[0]["priority_reasons"]
    assert 3 in selected["row_id"].tolist()
    assert summary["n_selected"] == 2
    assert (coverage["deficit"] > 0).any()


def test_select_next_rows_writes_outputs(tmp_path) -> None:
    module = _load_selector()
    selected, coverage, summary = module.select_next_rows(
        _training_table(),
        n_rows=3,
        target_per_cell=2,
    )
    out_dir = tmp_path / "priority"
    module.write_outputs(selected=selected, coverage=coverage, summary=summary, out_dir=out_dir)

    assert (out_dir / "next_label_priority.csv").exists()
    assert (out_dir / "next_label_priority_full.csv").exists()
    assert (out_dir / "next_label_row_ids.txt").read_text().strip().splitlines()[0] == "2"
    saved = json.loads((out_dir / "summary.json").read_text())
    assert saved["n_selected"] == 3
