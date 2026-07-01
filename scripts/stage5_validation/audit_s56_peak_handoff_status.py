#!/usr/bin/env python3
"""Audit S56 injected-peak handoff readiness from existing artifacts.

This script is intentionally read-only. It checks the expected peak-table,
gate, ranker, real-candidate, and review-queue artifacts and writes a compact
JSON/Markdown status report. It is the cheap gate before launching branch
comparison or human-label scale-up.
"""
from __future__ import annotations

import argparse
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
DEFAULT_ROOT = REPO_ROOT / "reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo"
DEFAULT_OUT_DIR = REPO_ROOT / "reports/stage5_validation/s56_peak_handoff_audit"
LAYOUTS = ("pdo", "orcd", "pdo_allhost")


def _json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError as exc:
        return {"_json_error": str(exc)}


def _json_passed(path: Path) -> bool:
    payload = _load_json(path)
    if not payload:
        return False
    if "_json_error" in payload:
        return False
    if "passed" in payload:
        return bool(payload["passed"])
    failures = payload.get("failures")
    if isinstance(failures, list):
        return len(failures) == 0
    return True


def _chunk_status(chunk_root: Path) -> dict[str, Any]:
    chunk_dirs = sorted(path for path in chunk_root.glob("chunk_*") if path.is_dir()) if chunk_root.exists() else []
    completed = [
        path.name
        for path in chunk_dirs
        if (path / "injection_bls_peaks.csv").exists()
        and (path / "injection_bls_peaks_summary.json").exists()
    ]
    incomplete = [path.name for path in chunk_dirs if path.name not in set(completed)]
    manifest = chunk_root.parent / "chunk_ids" / "manifest.txt"
    expected_chunks: int | None = None
    if manifest.exists():
        for line in manifest.read_text(errors="replace").splitlines():
            key, sep, value = line.partition("=")
            if sep and key.strip() == "n_chunks":
                try:
                    expected_chunks = int(value)
                except ValueError:
                    expected_chunks = None
                break
    return {
        "chunk_root": str(chunk_root),
        "expected_chunks": expected_chunks,
        "n_chunk_dirs": len(chunk_dirs),
        "n_completed_chunks": len(completed),
        "n_incomplete_chunks": len(incomplete),
        "incomplete_chunks_preview": incomplete[:10],
    }


def _step(name: str, ready: bool, path: Path | None, note: str, action: str) -> dict[str, Any]:
    return {
        "name": name,
        "ready": bool(ready),
        "path": str(path) if path is not None else "",
        "note": note,
        "next_action_if_missing": action,
    }


def _paths_for_layout(root: Path, repo_root: Path, layout: str) -> dict[str, Path]:
    if layout not in LAYOUTS:
        raise ValueError(f"unknown layout: {layout}")
    if layout == "pdo_allhost":
        peak_root = root / "peak_training"
        return {
            "peak_root": peak_root,
            "peak_table": peak_root / "s56_allhost_injection_bls_peaks.csv",
            "peak_verify": peak_root / "s56_allhost_injection_bls_peaks_verification.json",
            "gate_dir": root / "peak_training_gate_pdo",
            "ranker_dir": root / "peak_ranker_pdo",
            "selected_dir": repo_root / "reports/stage5_validation/s56_allhost_ranker_selected_real_candidates_pdo",
            "review_dir": repo_root / "reports/stage5_validation/s56_allhost_ranker_selected_real_review_queue_pdo",
            "leo_dir": repo_root / "reports/stage5_validation/s56_allhost_ranker_selected_real_leo_queue_pdo",
        }
    if layout == "orcd":
        peak_root = root / "peak_training_orcd_full"
        return {
            "peak_root": peak_root,
            "peak_table": peak_root / "s56_20k_injection_bls_peaks_orcd_full.csv",
            "peak_verify": peak_root / "s56_20k_injection_bls_peaks_orcd_full_verification.json",
            "gate_dir": root / "peak_training_gate_orcd",
            "ranker_dir": root / "peak_ranker_orcd",
            "selected_dir": repo_root / "reports/stage5_validation/s56_ranker_selected_real_candidates_orcd",
            "review_dir": repo_root / "reports/stage5_validation/s56_ranker_selected_real_review_queue_pdo",
            "leo_dir": repo_root / "reports/stage5_validation/s56_ranker_selected_real_leo_queue_pdo",
        }
    peak_root = root / "peak_training"
    return {
        "peak_root": peak_root,
        "peak_table": peak_root / "s56_20k_injection_bls_peaks_chunked.csv",
        "peak_verify": peak_root / "s56_20k_injection_bls_peaks_chunked_verification.json",
        "gate_dir": root / "peak_training_gate_pdo",
        "ranker_dir": root / "peak_ranker_pdo",
        "selected_dir": repo_root / "reports/stage5_validation/s56_ranker_selected_real_candidates_pdo",
        "review_dir": repo_root / "reports/stage5_validation/s56_ranker_selected_real_review_queue_pdo",
        "leo_dir": repo_root / "reports/stage5_validation/s56_ranker_selected_real_leo_queue_pdo",
    }


def audit_handoff_status(
    *,
    root: Path = DEFAULT_ROOT,
    repo_root: Path = REPO_ROOT,
    layout: str = "pdo",
) -> dict[str, Any]:
    paths = _paths_for_layout(root, repo_root, layout)
    peak_root = paths["peak_root"]
    chunked_table = paths["peak_table"]
    chunked_verify = paths["peak_verify"]
    chunk_status = _chunk_status(peak_root / "chunks")

    gate_dir = paths["gate_dir"]
    host_coverage_summary = gate_dir / "host_coverage" / "summary.json"
    gate_summary = gate_dir / "summary.json"
    failure_summary = gate_dir / "failure_modes" / "summary.json"
    ranker_dir = paths["ranker_dir"]
    ranker_summary = ranker_dir / "summary.json"
    ranker_verify = ranker_dir / "peak_table_verification.json"

    selected_dir = paths["selected_dir"]
    real_verify = selected_dir / "real_peak_table_verification.json"
    selected_ephemerides = selected_dir / "selected_ephemerides.csv"
    selected_verify = selected_dir / "selected_ephemerides_verification.json"
    selected_ready = selected_ephemerides.exists()
    if selected_verify.exists():
        selected_ready = selected_ready and _json_passed(selected_verify)

    review_dir = paths["review_dir"]
    review_queue = review_dir / "review_queue.csv"
    review_verify = review_dir / "verification.json"
    leo_dir = paths["leo_dir"]
    leo_queue = leo_dir / "review_queue.csv"
    leo_verify = leo_dir / "verification.json"

    standard_ready = chunked_table.exists() and _json_passed(chunked_verify)
    steps = [
        _step(
            "standard_injected_peak_table",
            standard_ready,
            chunked_verify,
            "Merged injected BLS peak table exists and passed strict schema verification.",
            "Wait for the active chunked build to merge and self-verify.",
        ),
        _step(
            "post_bls_peak_gate",
            gate_summary.exists(),
            gate_summary,
            "Recall@K / ranker-fixable / not-in-top-N gate exists.",
            "Run summarize_injection_peak_gate.py on the verified peak table.",
        ),
        _step(
            "injection_host_coverage_audited",
            host_coverage_summary.exists(),
            host_coverage_summary,
            "Injection host coverage relative to the S56 LC export has been quantified.",
            "Run audit_s56_injection_host_coverage.py on the staged export and injection manifests.",
        ),
        _step(
            "failure_mode_audit",
            failure_summary.exists(),
            failure_summary,
            "Failure-mode audit separates ranker-fixable and search-miss injections.",
            "Run audit_injection_bls_failure_modes.py before changing BLS branches.",
        ),
        _step(
            "injected_truth_peak_ranker",
            ranker_summary.exists() and _json_passed(ranker_verify),
            ranker_summary,
            "Injected-truth peak ranker trained after peak-table verification.",
            "Run run_s56_peak_ranker_review_pdo.sh or ORCD ranker Slurm wrapper.",
        ),
        _step(
            "real_bls_peak_table_verified",
            _json_passed(real_verify),
            real_verify,
            "Real S56 multi-peak BLS table passed verifier before ranker application.",
            "Run verify_real_bls_peak_table.py on the staged real candidate table.",
        ),
        _step(
            "ranker_selected_real_ephemerides",
            selected_ready,
            selected_verify if selected_verify.exists() else selected_ephemerides,
            "Peak ranker has selected real-candidate ephemerides; verification is used when available.",
            "Run apply_injection_peak_ranker.py on the real S56 BLS table.",
        ),
        _step(
            "skip_leo_review_queue_verified",
            review_queue.exists() and _json_passed(review_verify),
            review_verify,
            "Ranker-selected real queue exists and passed pre-LEO verifier.",
            "Build and verify a skip-LEO queue before expensive LEO rendering.",
        ),
        _step(
            "leo_review_queue_verified",
            leo_queue.exists() and _json_passed(leo_verify),
            leo_verify,
            "Ranker-selected real queue has LEO reports and passed verifier.",
            "Run run_s56_ranker_selected_real_leo_pdo.sh after the skip-LEO queue passes.",
        ),
    ]

    next_step = next((step for step in steps if not step["ready"]), None)
    all_required_before_human = all(
        step["ready"]
        for step in steps
        if step["name"] in {
            "standard_injected_peak_table",
            "post_bls_peak_gate",
            "injection_host_coverage_audited",
            "failure_mode_audit",
            "injected_truth_peak_ranker",
            "real_bls_peak_table_verified",
            "ranker_selected_real_ephemerides",
            "skip_leo_review_queue_verified",
        }
    )
    payload = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "root": str(root),
        "repo_root": str(repo_root),
        "layout": layout,
        "paths": {key: str(value) for key, value in paths.items()},
        "chunk_status": chunk_status,
        "steps": steps,
        "next_step": next_step,
        "ready_for_branch_comparison": bool(
            standard_ready and gate_summary.exists() and failure_summary.exists()
        ),
        "ready_for_skip_leo_human_queue": bool(all_required_before_human),
        "ready_for_leo_human_queue": bool(leo_queue.exists() and _json_passed(leo_verify)),
    }
    return payload


def render_markdown(payload: dict[str, Any]) -> str:
    chunk = payload["chunk_status"]
    lines = [
        "# S56 Peak Handoff Audit",
        "",
        f"Created UTC: `{payload['created_utc']}`",
        "",
        "## Chunked Peak Table",
        "",
        f"- completed chunks: `{chunk['n_completed_chunks']}` / `{chunk.get('expected_chunks') or '?'}`",
        f"- chunk dirs: `{chunk['n_chunk_dirs']}`",
        f"- incomplete preview: `{', '.join(chunk['incomplete_chunks_preview']) or 'none'}`",
        "",
        "## Gates",
        "",
        "| Gate | Ready | Artifact |",
        "|---|---:|---|",
    ]
    for step in payload["steps"]:
        ready = "yes" if step["ready"] else "no"
        artifact = step["path"] or ""
        lines.append(f"| `{step['name']}` | {ready} | `{artifact}` |")
    lines.extend(
        [
            "",
            "## Next Action",
            "",
        ]
    )
    next_step = payload.get("next_step")
    if next_step:
        lines.append(f"`{next_step['name']}`: {next_step['next_action_if_missing']}")
    else:
        lines.append("All audited gates are ready.")
    lines.extend(
        [
            "",
            f"- ready for branch comparison: `{payload['ready_for_branch_comparison']}`",
            f"- ready for skip-LEO human queue: `{payload['ready_for_skip_leo_human_queue']}`",
            f"- ready for LEO human queue: `{payload['ready_for_leo_human_queue']}`",
            "",
        ]
    )
    return "\n".join(lines)


def write_outputs(payload: dict[str, Any], out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=_json_default) + "\n"
    )
    (out_dir / "summary.md").write_text(render_markdown(payload))


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--layout", choices=LAYOUTS, default="pdo")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    payload = audit_handoff_status(root=args.root, repo_root=args.repo_root, layout=args.layout)
    write_outputs(payload, args.out_dir)
    print(render_markdown(payload))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
