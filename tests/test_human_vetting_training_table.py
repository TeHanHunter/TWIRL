from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_builder():
    path = REPO_ROOT / "scripts" / "stage5_validation" / "build_human_vetting_training_table.py"
    spec = importlib.util.spec_from_file_location("build_human_vetting_training_table", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    queue = pd.DataFrame(
        [
            {
                "review_id": "inj:a",
                "source_kind": "injection_recovery",
                "truth_source_kind": "injection_recovery",
                "source_bucket": "review_candidate",
                "tic": 1,
                "sector": 56,
                "period_d": 1.0,
                "t0_bjd": 2459825.0,
                "injection_id": "a",
                "topn_recovery_status": "bls_top1_recovered",
                "leo_class": "PC",
            },
            {
                "review_id": "inj:b",
                "source_kind": "injection_recovery",
                "truth_source_kind": "injection_recovery",
                "source_bucket": "review_candidate",
                "tic": 2,
                "sector": 56,
                "period_d": 2.0,
                "t0_bjd": 2459826.0,
                "injection_id": "b",
                "topn_recovery_status": "bls_peak_mismatch",
                "leo_class": "FA",
            },
            {
                "review_id": "real:3",
                "source_kind": "real_candidate",
                "truth_source_kind": "real_candidate",
                "source_bucket": "review_candidate",
                "tic": 3,
                "sector": 56,
                "period_d": 3.0,
                "t0_bjd": 2459827.0,
                "injection_id": "",
                "topn_recovery_status": "",
                "leo_class": "FA",
            },
        ]
    )
    queue_csv = tmp_path / "review_queue.csv"
    queue.to_csv(queue_csv, index=False)
    labels = pd.DataFrame(
        [
            {
                "row_id": 0,
                "candidate_key": "stale",
                "label": "uncertain",
                "label_source": "human",
                "labeler": "tehan",
                "notes": "",
                "updated_utc": "2026-06-26T00:00:00+00:00",
            },
            {
                "row_id": 0,
                "candidate_key": "1|56|1.0|2459825.0|review_candidate",
                "label": "planet_like",
                "label_source": "human",
                "labeler": "tehan",
                "notes": "latest wins",
                "updated_utc": "2026-06-26T00:01:00+00:00",
            },
            {
                "row_id": 1,
                "candidate_key": "2|56|2.0|2459826.0|review_candidate",
                "label": "instrumental_or_systematic",
                "label_source": "human",
                "labeler": "tehan",
                "notes": "",
                "updated_utc": "2026-06-26T00:02:00+00:00",
            },
            {
                "row_id": 2,
                "candidate_key": "3|56|3.0|2459827.0|review_candidate",
                "label": "skip",
                "label_source": "human",
                "labeler": "tehan",
                "notes": "",
                "updated_utc": "2026-06-26T00:03:00+00:00",
            },
        ]
    )
    labels_csv = tmp_path / "human_labels_vetted.csv"
    labels.to_csv(labels_csv, index=False)
    gate = pd.DataFrame(
        [
            {
                "injection_id": "a",
                "top1_match": True,
                "top20_match": True,
                "ranking_loss": False,
                "not_in_top20": False,
                "recovery_cell_top20_fraction": 0.9,
                "recovery_gate": "above_empirical_threshold",
            },
            {
                "injection_id": "b",
                "top1_match": False,
                "top20_match": False,
                "ranking_loss": False,
                "not_in_top20": True,
                "recovery_cell_top20_fraction": 0.1,
                "recovery_gate": "below_empirical_threshold",
            },
        ]
    )
    gate_csv = tmp_path / "gate.csv"
    gate.to_csv(gate_csv, index=False)
    return queue_csv, labels_csv, gate_csv


def test_build_human_training_table_flags_and_splits(tmp_path) -> None:
    module = _load_builder()
    queue_csv, labels_csv, gate_csv = _write_inputs(tmp_path)
    out_dir = tmp_path / "out"

    summary = module.build_training_table(
        queue_csv=queue_csv,
        labels_csv=labels_csv,
        out_dir=out_dir,
        recovery_gate_csv=gate_csv,
        validation_fraction=0.0,
        test_fraction=0.0,
        random_state=1,
    )

    table = pd.read_csv(out_dir / "human_vetting_training_table.csv")
    teacher = pd.read_csv(out_dir / "teacher_labeled_rows.csv")
    audit = pd.read_csv(out_dir / "audit_labeled_rows.csv")
    saved_summary = json.loads((out_dir / "summary.json").read_text())

    assert summary["n_labeled"] == 3
    assert saved_summary["n_teacher_rows"] == 2
    assert len(table) == 3
    assert len(teacher) == 2
    assert len(audit) == 2
    assert table.loc[0, "human_label"] == "planet_like"
    assert table.loc[0, "teacher_include"]
    assert table.loc[0, "bls_truth_match"]
    assert table.loc[1, "recovery_gate"] == "below_empirical_threshold"
    assert not table.loc[2, "teacher_include"]
    assert not table.loc[2, "audit_include"]
    assert set(teacher["training_split"]) == {"train"}


def test_latest_labels_uses_most_recent_row(tmp_path) -> None:
    module = _load_builder()
    _, labels_csv, _ = _write_inputs(tmp_path)

    latest = module.latest_labels(labels_csv)

    assert len(latest) == 3
    assert latest.loc[latest["row_id"].eq(0), "label"].item() == "planet_like"
