from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_audit():
    path = REPO_ROOT / "scripts" / "stage5_validation" / "audit_s56_peak_handoff_status.py"
    spec = importlib.util.spec_from_file_location("audit_s56_peak_handoff_status", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload) + "\n")


def test_peak_handoff_audit_reports_active_chunks(tmp_path: Path) -> None:
    module = _load_audit()
    root = tmp_path / "reports" / "s56_20k"
    chunk_root = root / "peak_training" / "chunks"
    for name, complete in (("chunk_000", True), ("chunk_001", False)):
        chunk_dir = chunk_root / name
        chunk_dir.mkdir(parents=True)
        if complete:
            (chunk_dir / "injection_bls_peaks.csv").write_text("injection_id\ninj\n")
            _write_json(chunk_dir / "injection_bls_peaks_summary.json", {"n": 1})
    manifest = root / "peak_training" / "chunk_ids" / "manifest.txt"
    manifest.parent.mkdir(parents=True)
    manifest.write_text("n_chunks=2\n")

    payload = module.audit_handoff_status(root=root, repo_root=tmp_path)

    assert payload["chunk_status"]["n_completed_chunks"] == 1
    assert payload["chunk_status"]["expected_chunks"] == 2
    assert payload["steps"][0]["name"] == "standard_injected_peak_table"
    assert not payload["steps"][0]["ready"]
    assert payload["next_step"]["name"] == "standard_injected_peak_table"
    assert not payload["ready_for_branch_comparison"]


def test_peak_handoff_default_repo_root_uses_cwd_for_stdin(tmp_path: Path) -> None:
    module = _load_audit()
    old_file = module.__file__
    old_cwd = Path.cwd()
    try:
        module.__file__ = "<stdin>"
        os.chdir(tmp_path)
        assert module._default_repo_root().resolve() == tmp_path.resolve()
    finally:
        module.__file__ = old_file
        os.chdir(old_cwd)


def test_peak_handoff_audit_detects_ready_gates(tmp_path: Path) -> None:
    module = _load_audit()
    root = tmp_path / "reports" / "s56_20k"
    repo = tmp_path
    peak_root = root / "peak_training"
    peak_root.mkdir(parents=True)
    (peak_root / "s56_20k_injection_bls_peaks_chunked.csv").write_text("injection_id\ninj\n")
    _write_json(peak_root / "s56_20k_injection_bls_peaks_chunked_verification.json", {"passed": True})
    _write_json(root / "peak_training_gate_pdo" / "summary.json", {"recall": 1})
    _write_json(root / "peak_training_gate_pdo" / "host_coverage" / "summary.json", {"n": 1})
    _write_json(root / "peak_training_gate_pdo" / "failure_modes" / "summary.json", {"n": 1})
    _write_json(root / "peak_ranker_pdo" / "peak_table_verification.json", {"passed": True})
    _write_json(root / "peak_ranker_pdo" / "summary.json", {"n": 1})
    selected_dir = repo / "reports/stage5_validation/s56_ranker_selected_real_candidates_pdo"
    _write_json(selected_dir / "real_peak_table_verification.json", {"passed": True})
    (selected_dir / "selected_ephemerides.csv").write_text("tic\n1\n")
    review_dir = repo / "reports/stage5_validation/s56_ranker_selected_real_review_queue_pdo"
    (review_dir / "review_queue.csv").parent.mkdir(parents=True)
    (review_dir / "review_queue.csv").write_text("tic\n1\n")
    _write_json(review_dir / "verification.json", {"passed": True})

    payload = module.audit_handoff_status(root=root, repo_root=repo)

    assert payload["ready_for_branch_comparison"]
    assert payload["ready_for_skip_leo_human_queue"]
    assert not payload["ready_for_leo_human_queue"]
    assert payload["next_step"]["name"] == "leo_review_queue_verified"


def test_peak_handoff_audit_supports_orcd_layout(tmp_path: Path) -> None:
    module = _load_audit()
    root = tmp_path / "reports" / "s56_20k"
    repo = tmp_path
    peak_root = root / "peak_training_orcd_full"
    peak_root.mkdir(parents=True)
    (peak_root / "s56_20k_injection_bls_peaks_orcd_full.csv").write_text("injection_id\ninj\n")
    _write_json(peak_root / "s56_20k_injection_bls_peaks_orcd_full_verification.json", {"passed": True})
    _write_json(root / "peak_training_gate_orcd" / "summary.json", {"recall": 1})
    _write_json(root / "peak_training_gate_orcd" / "host_coverage" / "summary.json", {"n": 1})
    _write_json(root / "peak_training_gate_orcd" / "failure_modes" / "summary.json", {"n": 1})
    _write_json(root / "peak_ranker_orcd" / "peak_table_verification.json", {"passed": True})
    _write_json(root / "peak_ranker_orcd" / "summary.json", {"n": 1})
    selected_dir = repo / "reports/stage5_validation/s56_ranker_selected_real_candidates_orcd"
    _write_json(selected_dir / "real_peak_table_verification.json", {"passed": True})
    (selected_dir / "selected_ephemerides.csv").parent.mkdir(parents=True, exist_ok=True)
    (selected_dir / "selected_ephemerides.csv").write_text("tic\n1\n")
    _write_json(selected_dir / "selected_ephemerides_verification.json", {"passed": True})

    payload = module.audit_handoff_status(root=root, repo_root=repo, layout="orcd")

    assert payload["layout"] == "orcd"
    assert payload["ready_for_branch_comparison"]
    assert not payload["ready_for_skip_leo_human_queue"]
    assert payload["next_step"]["name"] == "skip_leo_review_queue_verified"
    selected_step = [step for step in payload["steps"] if step["name"] == "ranker_selected_real_ephemerides"][0]
    assert selected_step["ready"]
    assert selected_step["path"].endswith("selected_ephemerides_verification.json")


def test_peak_handoff_audit_supports_allhost_pdo_layout(tmp_path: Path) -> None:
    module = _load_audit()
    root = tmp_path / "reports" / "s56_allhost"
    repo = tmp_path
    peak_root = root / "peak_training"
    peak_root.mkdir(parents=True)
    (peak_root / "s56_allhost_injection_bls_peaks.csv").write_text("injection_id\ninj\n")
    _write_json(peak_root / "s56_allhost_injection_bls_peaks_verification.json", {"passed": True})
    _write_json(root / "peak_training_gate_pdo" / "summary.json", {"recall": 1})
    _write_json(root / "peak_training_gate_pdo" / "host_coverage" / "summary.json", {"n": 1})
    _write_json(root / "peak_training_gate_pdo" / "failure_modes" / "summary.json", {"n": 1})
    _write_json(root / "peak_ranker_pdo" / "peak_table_verification.json", {"passed": True})
    _write_json(root / "peak_ranker_pdo" / "summary.json", {"n": 1})
    selected_dir = repo / "reports/stage5_validation/s56_allhost_ranker_selected_real_candidates_pdo"
    _write_json(selected_dir / "real_peak_table_verification.json", {"passed": True})
    (selected_dir / "selected_ephemerides.csv").parent.mkdir(parents=True, exist_ok=True)
    (selected_dir / "selected_ephemerides.csv").write_text("tic\n1\n")
    review_dir = repo / "reports/stage5_validation/s56_allhost_ranker_selected_real_review_queue_pdo"
    (review_dir / "review_queue.csv").parent.mkdir(parents=True)
    (review_dir / "review_queue.csv").write_text("tic\n1\n")
    _write_json(review_dir / "verification.json", {"passed": True})

    payload = module.audit_handoff_status(root=root, repo_root=repo, layout="pdo_allhost")

    assert payload["layout"] == "pdo_allhost"
    assert payload["ready_for_branch_comparison"]
    assert payload["ready_for_skip_leo_human_queue"]
    assert payload["next_step"]["name"] == "leo_review_queue_verified"


def test_peak_handoff_audit_uses_failed_selected_verification(tmp_path: Path) -> None:
    module = _load_audit()
    root = tmp_path / "reports" / "s56_20k"
    repo = tmp_path
    selected_dir = repo / "reports/stage5_validation/s56_ranker_selected_real_candidates_orcd"
    (selected_dir / "selected_ephemerides.csv").parent.mkdir(parents=True)
    (selected_dir / "selected_ephemerides.csv").write_text("tic\n1\n")
    _write_json(selected_dir / "selected_ephemerides_verification.json", {"passed": False})

    payload = module.audit_handoff_status(root=root, repo_root=repo, layout="orcd")

    selected_step = [step for step in payload["steps"] if step["name"] == "ranker_selected_real_ephemerides"][0]
    assert not selected_step["ready"]
