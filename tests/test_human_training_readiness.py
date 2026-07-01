from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_auditor():
    path = REPO_ROOT / "scripts" / "stage5_validation" / "audit_human_training_readiness.py"
    spec = importlib.util.spec_from_file_location("audit_human_training_readiness", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _row(
    idx: int,
    *,
    label: str,
    source_kind: str = "injection_recovery",
    split: str = "train",
    recovered: bool = True,
) -> dict[str, object]:
    return {
        "row_id": idx,
        "is_labeled": True,
        "teacher_include": label != "uncertain",
        "audit_include": True,
        "human_label": label,
        "teacher_target": "" if label == "uncertain" else label,
        "training_split": split if label != "uncertain" else "unlabeled_or_audit",
        "source_kind": source_kind,
        "topn_recovery_status": "bls_top1_recovered" if recovered else "bls_peak_mismatch",
        "bls_truth_match": recovered,
        "leo_class": "PC" if recovered else "FA",
        "tmag": 18.0,
        "truth_period_d": 1.0,
        "truth_radius_rearth": 6.0,
    }


def test_current_style_labels_are_smoke_only(tmp_path) -> None:
    module = _load_auditor()
    rows = []
    for i in range(25):
        rows.append(_row(i, label="planet_like", recovered=True))
    for i in range(25, 50):
        rows.append(_row(i, label="instrumental_or_systematic", recovered=False))
    for i in range(50, 60):
        rows.append(_row(i, label="uncertain", recovered=False))
    table = pd.DataFrame(rows)

    summary, deficits, markdown = module.audit_training_readiness(
        table,
        config=module.ReadinessConfig(
            min_smoke_labeled=50,
            min_smoke_positive=20,
            min_smoke_negative=20,
            min_smoke_bls_truth_match=10,
            min_smoke_bls_mismatch=20,
            min_binary_teacher_total=80,
            min_binary_teacher_per_class=30,
            min_teacher_total=120,
            min_teacher_per_class=30,
            min_train_per_class=20,
            min_validation_per_class=5,
            min_test_per_class=5,
            min_real_teacher=20,
        ),
    )

    assert summary["recommendation"] == "ready_for_injection_visibility_smoke_only"
    assert summary["gates"]["injection_visibility_smoke"]["status"] == "pass"
    assert summary["gates"]["binary_teacher_smoke"]["status"] == "fail"
    assert summary["gates"]["object_teacher_training"]["status"] == "fail"
    assert "Real false positives cannot be replaced" in deficits["note"].fillna("").str.cat(sep=" ")
    assert "ready_for_injection_visibility_smoke_only" in markdown


def test_balanced_real_table_can_pass_object_teacher_gate(tmp_path) -> None:
    module = _load_auditor()
    labels = [
        "planet_like",
        "eclipsing_binary_or_pceb",
        "stellar_variability",
        "instrumental_or_systematic",
        "centroid_contaminant",
    ]
    rows = []
    idx = 0
    for label in labels:
        for split in ("train", "train", "validation", "test"):
            rows.append(
                _row(
                    idx,
                    label=label,
                    source_kind="real_candidate",
                    split=split,
                    recovered=label == "planet_like",
                )
            )
            idx += 1
    table = pd.DataFrame(rows)

    summary, deficits, _ = module.audit_training_readiness(
        table,
        config=module.ReadinessConfig(
            min_smoke_labeled=10,
            min_smoke_positive=2,
            min_smoke_negative=2,
            min_smoke_bls_truth_match=2,
            min_smoke_bls_mismatch=2,
            min_binary_teacher_total=4,
            min_binary_teacher_per_class=2,
            min_teacher_total=20,
            min_teacher_per_class=4,
            min_train_per_class=2,
            min_validation_per_class=1,
            min_test_per_class=1,
            min_real_teacher=20,
            target_per_coverage_cell=1,
        ),
    )

    assert summary["recommendation"] == "ready_for_object_teacher_training"
    assert summary["gates"]["object_teacher_training"]["status"] == "pass"
    object_failures = deficits[
        deficits["gate"].eq("object_teacher_training") & deficits["status"].eq("fail")
    ]
    assert object_failures.empty


def test_build_readiness_audit_writes_outputs(tmp_path) -> None:
    module = _load_auditor()
    table = pd.DataFrame(
        [_row(0, label="planet_like", recovered=True), _row(1, label="instrumental_or_systematic", recovered=False)]
    )
    training_table = tmp_path / "human_vetting_training_table.csv"
    table.to_csv(training_table, index=False)
    priority = pd.DataFrame(
        [
            {"row_id": 2, "topn_recovery_status": "bls_top1_recovered", "leo_class": "PC"},
            {"row_id": 3, "topn_recovery_status": "bls_peak_mismatch", "leo_class": "FA"},
        ]
    )
    priority_table = tmp_path / "next_label_priority.csv"
    priority.to_csv(priority_table, index=False)

    summary = module.build_readiness_audit(
        training_table=training_table,
        priority_table=priority_table,
        out_dir=tmp_path / "audit",
        config=module.ReadinessConfig(
            min_smoke_labeled=2,
            min_smoke_positive=1,
            min_smoke_negative=1,
            min_smoke_bls_truth_match=1,
            min_smoke_bls_mismatch=1,
            min_binary_teacher_total=10,
            min_binary_teacher_per_class=5,
            min_teacher_total=10,
            min_teacher_per_class=5,
            min_train_per_class=2,
            min_validation_per_class=1,
            min_test_per_class=1,
            min_real_teacher=5,
        ),
    )

    out_dir = tmp_path / "audit"
    assert (out_dir / "summary.json").exists()
    assert (out_dir / "summary.md").exists()
    assert (out_dir / "label_deficits.csv").exists()
    saved = json.loads((out_dir / "summary.json").read_text())
    assert saved["priority_summary"]["n_priority_rows"] == 2
    assert summary["recommendation"] == "ready_for_injection_visibility_smoke_only"
