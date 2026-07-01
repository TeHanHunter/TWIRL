from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_audit():
    path = REPO_ROOT / "scripts" / "stage5_validation" / "audit_s56_ml_handoff_readiness.py"
    spec = importlib.util.spec_from_file_location("audit_s56_ml_handoff_readiness", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload) + "\n")


def _write_csv(path: Path, rows: int, header: str = "review_id,label") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    body = "\n".join(f"row{i},planet_like" for i in range(rows))
    path.write_text(f"{header}\n{body}\n")


def _write_common_ready_artifacts(tmp_path: Path) -> dict[str, Path]:
    queue_dir = tmp_path / "queue"
    selected_dir = tmp_path / "selected"
    tensor_summary = tmp_path / "tensor" / "summary.json"
    synthetic_train = tmp_path / "synthetic_train" / "summary.json"
    real_train = tmp_path / "real_train" / "summary.json"
    torch_info = tmp_path / "envs" / "torch-build-info.txt"

    _write_json(
        queue_dir / "verification.json",
        {
            "passed": True,
            "n_rows": 1000,
            "leo_class_counts": {"FA": 990, "FP": 5, "PC": 5},
        },
    )
    _write_json(
        tensor_summary,
        {
            "n_tensor_rows": 128,
            "shape": [128, 3, 257],
            "torch": {"cuda_available": True, "device": "cuda:0"},
        },
    )
    _write_json(
        synthetic_train,
        {
            "synthetic_label_smoke": True,
            "torch_cuda_available": True,
            "device": "cuda:0",
            "split_counts": {"train": 76, "validation": 26, "test": 26},
        },
    )
    torch_info.parent.mkdir(parents=True, exist_ok=True)
    torch_info.write_text("torch=2.11.0+cu128 cuda_available=False\n")
    return {
        "queue_dir": queue_dir,
        "orcd_selected_dir": selected_dir,
        "tensor_summary": tensor_summary,
        "synthetic_train_summary": synthetic_train,
        "real_train_summary": real_train,
        "torch_build_info": torch_info,
    }


def test_current_triage_ready_but_labels_pending(tmp_path: Path) -> None:
    module = _load_audit()
    paths = _write_common_ready_artifacts(tmp_path)

    summary = module.build_ml_handoff_readiness(
        **paths,
        out_dir=tmp_path / "out",
        min_queue_rows=1000,
        min_labels_for_audit=1,
    )

    checks = {check["name"]: check for check in summary["checks"]}
    assert checks["pdo_leo_queue_verified"]["status"] == "pass"
    assert checks["human_labels_present"]["status"] == "pending"
    assert summary["flags"]["ready_for_human_triage"]
    assert not summary["flags"]["ready_for_post_label_audit"]
    assert not summary["flags"]["ready_to_submit_real_h200_training_smoke"]
    assert summary["next_step"]["name"] == "human_labels_present"
    assert (tmp_path / "out" / "summary.json").exists()
    assert (tmp_path / "out" / "summary.md").exists()


def test_ready_to_submit_real_label_h200_smoke(tmp_path: Path) -> None:
    module = _load_audit()
    paths = _write_common_ready_artifacts(tmp_path)
    queue_dir = paths["queue_dir"]
    selected_dir = paths["orcd_selected_dir"]

    _write_csv(queue_dir / "human_labels_vetted.csv", 150)
    _write_json(
        queue_dir / "human_training_readiness" / "summary.json",
        {
            "recommendation": "ready_for_binary_teacher_smoke",
            "n_labeled": 150,
            "n_teacher_rows": 120,
            "n_real_teacher_rows": 120,
        },
    )
    _write_json(selected_dir / "selected_ephemerides_verification.json", {"passed": True})
    _write_csv(selected_dir / "selected_ephemerides.csv", 10, header="tic,period_d")

    summary = module.build_ml_handoff_readiness(
        **paths,
        out_dir=tmp_path / "out",
        min_queue_rows=1000,
        min_labels_for_audit=1,
    )

    assert summary["flags"]["ready_for_post_label_audit"]
    assert summary["flags"]["ready_to_submit_real_h200_training_smoke"]
    assert summary["flags"]["orcd_selected_outputs_ready_for_pdo_leo"]
    assert not summary["flags"]["real_h200_training_smoke_complete"]
    assert summary["next_step"]["name"] == "h200_real_label_train_smoke_complete"


def test_real_label_h200_smoke_complete(tmp_path: Path) -> None:
    module = _load_audit()
    paths = _write_common_ready_artifacts(tmp_path)
    queue_dir = paths["queue_dir"]

    _write_csv(queue_dir / "human_labels_vetted.csv", 150)
    _write_json(
        queue_dir / "human_training_readiness" / "summary.json",
        {
            "recommendation": "ready_for_object_teacher_training",
            "n_labeled": 300,
            "n_teacher_rows": 250,
            "n_real_teacher_rows": 180,
        },
    )
    _write_json(
        paths["real_train_summary"],
        {
            "synthetic_label_smoke": False,
            "torch_cuda_available": True,
            "device": "cuda:0",
            "n_training_rows": 180,
        },
    )

    summary = module.build_ml_handoff_readiness(
        **paths,
        out_dir=tmp_path / "out",
        min_queue_rows=1000,
        min_labels_for_audit=1,
    )

    assert summary["flags"]["real_h200_training_smoke_complete"]
    assert summary["next_step"]["name"] == "scale_after_smoke_review"
