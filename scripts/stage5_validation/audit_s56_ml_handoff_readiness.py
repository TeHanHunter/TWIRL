#!/usr/bin/env python3
"""Audit readiness for S56 human labels and first H200 ML handoff.

This is a local, read-only gate report. It intentionally does not SSH to PDO or
ORCD; it only inspects artifacts already present in the checkout. Use it after
syncing PDO/ORCD outputs to decide whether the next step is more human labels,
post-label audit, ORCD selected-output rendering, or the first real-label H200
training smoke.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any


def _default_repo_root() -> Path:
    script_path = Path(__file__)
    if not script_path.exists():
        return Path.cwd()
    return script_path.resolve().parents[2]


REPO_ROOT = _default_repo_root()
DEFAULT_QUEUE_DIR = REPO_ROOT / "reports/stage5_validation/s56_allhost_ranker_selected_real_leo_queue_pdo"
DEFAULT_ORCD_SELECTED_DIR = REPO_ROOT / "reports/stage5_validation/s56_ranker_selected_real_candidates_orcd"
DEFAULT_TENSOR_SUMMARY = REPO_ROOT / "reports/stage5_validation/s56_candidate_tensor_smoke_h200_torch/summary.json"
DEFAULT_SYNTHETIC_TRAIN_SUMMARY = (
    REPO_ROOT / "reports/stage5_validation/s56_candidate_tensor_train_synthetic_smoke_h200/summary.json"
)
DEFAULT_REAL_TRAIN_SUMMARY = REPO_ROOT / "reports/stage5_validation/s56_candidate_tensor_train_smoke_h200/summary.json"
DEFAULT_TORCH_BUILD_INFO = REPO_ROOT / "reports/stage5_validation/orcd_envs/twirl-s56-torch-build-info.txt"
DEFAULT_OUT_DIR = REPO_ROOT / "reports/stage5_validation/s56_ml_handoff_readiness"

TRAINING_RECOMMENDATIONS = {"ready_for_binary_teacher_smoke", "ready_for_object_teacher_training"}
ANY_LABEL_RECOMMENDATIONS = TRAINING_RECOMMENDATIONS | {"ready_for_injection_visibility_smoke_only"}


@dataclass(frozen=True)
class Check:
    name: str
    status: str
    ready: bool
    path: str
    evidence: str
    next_action: str

    def as_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "status": self.status,
            "ready": self.ready,
            "path": self.path,
            "evidence": self.evidence,
            "next_action": self.next_action,
        }


def _json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _load_json(path: Path) -> tuple[dict[str, Any], str]:
    if not path.exists():
        return {}, "missing"
    try:
        payload = json.loads(path.read_text())
    except json.JSONDecodeError as exc:
        return {"_json_error": str(exc)}, "json_error"
    if not isinstance(payload, dict):
        return {"_json_error": "top-level JSON is not an object"}, "json_error"
    return payload, "ok"


def _json_gate_status(path: Path) -> tuple[str, bool, str, dict[str, Any]]:
    payload, load_status = _load_json(path)
    if load_status == "missing":
        return "pending", False, "artifact is missing", payload
    if load_status == "json_error":
        return "fail", False, str(payload.get("_json_error", "JSON parse failed")), payload
    if "passed" in payload:
        passed = bool(payload["passed"])
        return ("pass" if passed else "fail"), passed, f"passed={passed}", payload
    failures = payload.get("failures")
    if isinstance(failures, list):
        passed = len(failures) == 0
        return ("pass" if passed else "fail"), passed, f"failures={len(failures)}", payload
    return "pass", True, "JSON exists", payload


def _csv_row_count(path: Path) -> int | None:
    if not path.exists() or path.stat().st_size == 0:
        return None
    with path.open(newline="") as handle:
        return sum(1 for _ in csv.DictReader(handle))


def _status_from_count(path: Path, count: int | None, min_rows: int) -> tuple[str, bool, str]:
    if count is None:
        return "pending", False, "artifact is missing or empty"
    ready = count >= min_rows
    return ("pass" if ready else "pending"), ready, f"rows={count}, required>={min_rows}"


def _check_queue(queue_dir: Path, min_rows: int) -> Check:
    verify = queue_dir / "verification.json"
    status, passed, evidence, payload = _json_gate_status(verify)
    n_rows = int(payload.get("n_rows", 0) or 0)
    if passed and n_rows < min_rows:
        status, passed = "pending", False
    details = [evidence, f"n_rows={n_rows}"]
    if payload.get("leo_class_counts"):
        details.append(f"leo_class_counts={payload['leo_class_counts']}")
    return Check(
        "pdo_leo_queue_verified",
        status,
        passed,
        str(verify),
        "; ".join(details),
        "Keep the PDO app running if this passes; otherwise rebuild or verify the LEO queue.",
    )


def _check_human_labels(queue_dir: Path, min_labels: int) -> Check:
    labels = queue_dir / "human_labels_vetted.csv"
    count = _csv_row_count(labels)
    status, ready, evidence = _status_from_count(labels, count, min_labels)
    return Check(
        "human_labels_present",
        status,
        ready,
        str(labels),
        evidence,
        "Use the vetting app; the PDO label-audit monitor or run_s56_allhost_real_label_audit_pdo.sh will rebuild readiness after labels appear.",
    )


def _check_human_training_readiness(queue_dir: Path) -> Check:
    summary_path = queue_dir / "human_training_readiness" / "summary.json"
    payload, load_status = _load_json(summary_path)
    if load_status == "missing":
        return Check(
            "human_training_readiness",
            "pending",
            False,
            str(summary_path),
            "post-label readiness audit has not been run",
            "After labels exist, run run_s56_allhost_real_label_audit_pdo.sh on PDO and sync the outputs.",
        )
    if load_status == "json_error":
        return Check(
            "human_training_readiness",
            "fail",
            False,
            str(summary_path),
            str(payload.get("_json_error", "JSON parse failed")),
            "Regenerate the human-training readiness audit.",
        )
    recommendation = str(payload.get("recommendation", ""))
    ready = recommendation in TRAINING_RECOMMENDATIONS
    if recommendation in ANY_LABEL_RECOMMENDATIONS and not ready:
        status = "pending"
        action = "Label more real candidates before launching the real-label H200 training smoke."
    elif ready:
        status = "pass"
        action = "This gate is ready for the bounded real-label H200 training smoke."
    else:
        status = "pending"
        action = "Continue labeling and rerun the post-label readiness audit."
    evidence = (
        f"recommendation={recommendation}; "
        f"n_labeled={payload.get('n_labeled', 0)}; "
        f"n_teacher_rows={payload.get('n_teacher_rows', 0)}; "
        f"n_real_teacher_rows={payload.get('n_real_teacher_rows', 0)}"
    )
    return Check("human_training_readiness", status, ready, str(summary_path), evidence, action)


def _check_orcd_selected(selected_dir: Path) -> Check:
    verify = selected_dir / "selected_ephemerides_verification.json"
    selected = selected_dir / "selected_ephemerides.csv"
    status, passed, evidence, _ = _json_gate_status(verify)
    count = _csv_row_count(selected)
    if passed and (count is None or count == 0):
        status, passed = "pending", False
    count_note = "rows=missing" if count is None else f"rows={count}"
    return Check(
        "orcd_selected_ephemerides_synced",
        status,
        passed,
        str(verify if verify.exists() else selected),
        f"{evidence}; {count_note}",
        "Wait for the ORCD apply dependency chain and sync selected ephemerides back locally/PDO.",
    )


def _check_torch_build_info(path: Path) -> Check:
    if not path.exists():
        return Check(
            "orcd_torch_env_documented",
            "pending",
            False,
            str(path),
            "build-info file is missing",
            "Build or sync the ORCD torch environment report.",
        )
    text = path.read_text(errors="replace")
    has_cuda_torch = "torch=" in text and "cuda_available" in text
    return Check(
        "orcd_torch_env_documented",
        "pass" if has_cuda_torch else "pending",
        has_cuda_torch,
        str(path),
        "torch build info present" if has_cuda_torch else "build info present but CUDA torch line was not recognized",
        "Use the documented twirl-s56-torch env for H200 training smokes.",
    )


def _check_tensor_summary(path: Path) -> Check:
    payload, load_status = _load_json(path)
    if load_status == "missing":
        return Check(
            "h200_tensor_smoke",
            "pending",
            False,
            str(path),
            "summary is missing",
            "Run scripts/orcd/run_s56_orcd_pilot.sh --run h200-torch-tensor-smoke.",
        )
    if load_status == "json_error":
        return Check("h200_tensor_smoke", "fail", False, str(path), str(payload.get("_json_error")), "Regenerate the tensor smoke.")
    torch = payload.get("torch", {}) if isinstance(payload.get("torch"), dict) else {}
    n_tensor = int(payload.get("n_tensor_rows", 0) or 0)
    cuda_ok = bool(torch.get("cuda_available"))
    ready = n_tensor > 0 and cuda_ok
    evidence = f"n_tensor_rows={n_tensor}; shape={payload.get('shape')}; cuda_available={cuda_ok}; device={torch.get('device')}"
    return Check(
        "h200_tensor_smoke",
        "pass" if ready else "fail",
        ready,
        str(path),
        evidence,
        "Fix tensor extraction or H200 CUDA visibility before training.",
    )


def _check_synthetic_train(path: Path) -> Check:
    payload, load_status = _load_json(path)
    if load_status == "missing":
        return Check(
            "h200_synthetic_train_smoke",
            "pending",
            False,
            str(path),
            "summary is missing",
            "Run scripts/orcd/run_s56_orcd_pilot.sh --run h200-tensor-train-synthetic-smoke.",
        )
    if load_status == "json_error":
        return Check(
            "h200_synthetic_train_smoke",
            "fail",
            False,
            str(path),
            str(payload.get("_json_error")),
            "Regenerate the synthetic H200 training smoke.",
        )
    split_counts = payload.get("split_counts", {}) if isinstance(payload.get("split_counts"), dict) else {}
    split_ok = all(int(split_counts.get(name, 0) or 0) > 0 for name in ("train", "validation", "test"))
    cuda_ok = bool(payload.get("torch_cuda_available"))
    synthetic = bool(payload.get("synthetic_label_smoke"))
    ready = split_ok and cuda_ok and synthetic
    evidence = (
        f"synthetic_label_smoke={synthetic}; cuda={cuda_ok}; device={payload.get('device')}; "
        f"split_counts={split_counts}"
    )
    return Check(
        "h200_synthetic_train_smoke",
        "pass" if ready else "fail",
        ready,
        str(path),
        evidence,
        "Fix the trainer/env plumbing before using real labels.",
    )


def _check_real_train(path: Path) -> Check:
    payload, load_status = _load_json(path)
    if load_status == "missing":
        return Check(
            "h200_real_label_train_smoke_complete",
            "pending",
            False,
            str(path),
            "real-label training smoke has not been run",
            "Run scripts/orcd/run_s56_orcd_pilot.sh --run h200-tensor-train-smoke after labels pass readiness.",
        )
    if load_status == "json_error":
        return Check(
            "h200_real_label_train_smoke_complete",
            "fail",
            False,
            str(path),
            str(payload.get("_json_error")),
            "Regenerate the real-label H200 training smoke.",
        )
    cuda_ok = bool(payload.get("torch_cuda_available"))
    synthetic = bool(payload.get("synthetic_label_smoke"))
    ready = cuda_ok and not synthetic and int(payload.get("n_training_rows", 0) or 0) > 0
    evidence = (
        f"synthetic_label_smoke={synthetic}; cuda={cuda_ok}; device={payload.get('device')}; "
        f"n_training_rows={payload.get('n_training_rows', 0)}"
    )
    return Check(
        "h200_real_label_train_smoke_complete",
        "pass" if ready else "fail",
        ready,
        str(path),
        evidence,
        "Inspect trainer outputs before scaling beyond the bounded smoke.",
    )


def build_ml_handoff_readiness(
    *,
    queue_dir: Path = DEFAULT_QUEUE_DIR,
    orcd_selected_dir: Path = DEFAULT_ORCD_SELECTED_DIR,
    tensor_summary: Path = DEFAULT_TENSOR_SUMMARY,
    synthetic_train_summary: Path = DEFAULT_SYNTHETIC_TRAIN_SUMMARY,
    real_train_summary: Path = DEFAULT_REAL_TRAIN_SUMMARY,
    torch_build_info: Path = DEFAULT_TORCH_BUILD_INFO,
    out_dir: Path = DEFAULT_OUT_DIR,
    min_queue_rows: int = 1000,
    min_labels_for_audit: int = 1,
) -> dict[str, Any]:
    checks = [
        _check_queue(queue_dir, min_queue_rows),
        _check_human_labels(queue_dir, min_labels_for_audit),
        _check_human_training_readiness(queue_dir),
        _check_orcd_selected(orcd_selected_dir),
        _check_torch_build_info(torch_build_info),
        _check_tensor_summary(tensor_summary),
        _check_synthetic_train(synthetic_train_summary),
        _check_real_train(real_train_summary),
    ]
    by_name = {check.name: check for check in checks}
    flags = {
        "ready_for_human_triage": by_name["pdo_leo_queue_verified"].ready,
        "ready_for_post_label_audit": by_name["pdo_leo_queue_verified"].ready
        and by_name["human_labels_present"].ready,
        "ready_to_submit_real_h200_training_smoke": all(
            by_name[name].ready
            for name in (
                "pdo_leo_queue_verified",
                "human_training_readiness",
                "orcd_torch_env_documented",
                "h200_tensor_smoke",
                "h200_synthetic_train_smoke",
            )
        ),
        "orcd_selected_outputs_ready_for_pdo_leo": by_name["orcd_selected_ephemerides_synced"].ready,
        "real_h200_training_smoke_complete": by_name["h200_real_label_train_smoke_complete"].ready,
    }
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "queue_dir": str(queue_dir),
        "orcd_selected_dir": str(orcd_selected_dir),
        "checks": [check.as_dict() for check in checks],
        "flags": flags,
        "next_step": _next_step(checks, flags),
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n")
    (out_dir / "summary.md").write_text(render_markdown(summary))
    print(json.dumps(summary, indent=2, sort_keys=True, default=_json_default))
    return summary


def _next_step(checks: list[Check], flags: dict[str, bool]) -> dict[str, str]:
    if not flags["ready_for_human_triage"]:
        check = next(check for check in checks if check.name == "pdo_leo_queue_verified")
        return {"name": check.name, "action": check.next_action}
    if not flags["ready_for_post_label_audit"]:
        check = next(check for check in checks if check.name == "human_labels_present")
        return {"name": check.name, "action": check.next_action}
    if not next(check for check in checks if check.name == "human_training_readiness").ready:
        check = next(check for check in checks if check.name == "human_training_readiness")
        return {"name": check.name, "action": check.next_action}
    if not flags["ready_to_submit_real_h200_training_smoke"]:
        check = next(check for check in checks if not check.ready and check.name.startswith(("orcd_torch", "h200_")))
        return {"name": check.name, "action": check.next_action}
    if not flags["real_h200_training_smoke_complete"]:
        check = next(check for check in checks if check.name == "h200_real_label_train_smoke_complete")
        return {"name": check.name, "action": check.next_action}
    return {"name": "scale_after_smoke_review", "action": "Review held-out metrics before scaling H200 training."}


def render_markdown(summary: dict[str, Any]) -> str:
    lines = [
        "# S56 ML Handoff Readiness",
        "",
        f"Created UTC: `{summary['created_utc']}`",
        "",
        "## Status",
        "",
        "| Gate | Status | Evidence |",
        "|---|---:|---|",
    ]
    for check in summary["checks"]:
        lines.append(f"| `{check['name']}` | `{check['status']}` | {check['evidence']} |")
    lines.extend(["", "## Flags", "", "| Flag | Ready |", "|---|---:|"])
    for name, ready in summary["flags"].items():
        lines.append(f"| `{name}` | `{ready}` |")
    next_step = summary["next_step"]
    lines.extend(["", "## Next Step", "", f"`{next_step['name']}`: {next_step['action']}", ""])
    return "\n".join(lines)


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queue-dir", type=Path, default=DEFAULT_QUEUE_DIR)
    parser.add_argument("--orcd-selected-dir", type=Path, default=DEFAULT_ORCD_SELECTED_DIR)
    parser.add_argument("--tensor-summary", type=Path, default=DEFAULT_TENSOR_SUMMARY)
    parser.add_argument("--synthetic-train-summary", type=Path, default=DEFAULT_SYNTHETIC_TRAIN_SUMMARY)
    parser.add_argument("--real-train-summary", type=Path, default=DEFAULT_REAL_TRAIN_SUMMARY)
    parser.add_argument("--torch-build-info", type=Path, default=DEFAULT_TORCH_BUILD_INFO)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--min-queue-rows", type=int, default=1000)
    parser.add_argument("--min-labels-for-audit", type=int, default=1)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    build_ml_handoff_readiness(
        queue_dir=args.queue_dir,
        orcd_selected_dir=args.orcd_selected_dir,
        tensor_summary=args.tensor_summary,
        synthetic_train_summary=args.synthetic_train_summary,
        real_train_summary=args.real_train_summary,
        torch_build_info=args.torch_build_info,
        out_dir=args.out_dir,
        min_queue_rows=args.min_queue_rows,
        min_labels_for_audit=args.min_labels_for_audit,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
