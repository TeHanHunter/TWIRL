#!/usr/bin/env python3
"""Write immutable provenance for a completed sector A2v1 recovery run."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from importlib import metadata
import json
from pathlib import Path
import platform
import subprocess
import sys


def _json(path: Path) -> dict:
    return json.loads(path.read_text())


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _distribution_versions() -> dict[str, str]:
    names = (
        "numpy",
        "pandas",
        "h5py",
        "astropy",
        "batman-package",
        "scipy",
        "scikit-learn",
        "matplotlib",
        "PyYAML",
    )
    versions = {}
    for name in names:
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = "not-installed"
    return versions


def _git_output(repo: Path, *args: str) -> str:
    return subprocess.check_output(
        ["git", *args], cwd=repo, text=True, stderr=subprocess.STDOUT
    ).strip()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--input-verification", type=Path, required=True)
    parser.add_argument("--schedule-summary", type=Path, required=True)
    parser.add_argument("--teacher-summary", type=Path, required=True)
    parser.add_argument("--final-summary", type=Path, required=True)
    parser.add_argument("--product-audit", type=Path, required=True)
    parser.add_argument("--checkpoints", type=Path, nargs=5, required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    args = parser.parse_args()

    import yaml

    input_verification = _json(args.input_verification)
    schedule = _json(args.schedule_summary)
    teacher = _json(args.teacher_summary)
    final = _json(args.final_summary)
    product_audit = _json(args.product_audit)
    config = yaml.safe_load(args.config.read_text())
    tracked_clean = (
        subprocess.run(["git", "diff", "--quiet"], cwd=args.repo).returncode == 0
        and subprocess.run(
            ["git", "diff", "--cached", "--quiet"], cwd=args.repo
        ).returncode
        == 0
    )
    payload = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "workflow": str(config["name"]),
        "sector": int(config["sector"]),
        "evaluation_only": True,
        "teacher_retraining_allowed": False,
        "code": {
            "repo": str(args.repo.resolve()),
            "git_sha": _git_output(args.repo, "rev-parse", "HEAD"),
            "git_describe": _git_output(args.repo, "describe", "--always", "--dirty"),
            "tracked_code_clean": tracked_clean,
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "distributions": _distribution_versions(),
            "teacher_torch_version": teacher.get("torch_version", ""),
            "teacher_torch_cuda_version": teacher.get("torch_cuda_version", ""),
            "teacher_cuda_device_name": teacher.get("cuda_device_name", ""),
        },
        "configuration": config,
        "config_sha256": _sha256(args.config),
        "input_verification": input_verification,
        "source_checksums": input_verification.get("hashes", {}),
        "checkpoint_sha256": {
            str(path.resolve()): _sha256(path) for path in args.checkpoints
        },
        "host_exclusion_and_selection": {
            key: schedule.get(key)
            for key in (
                "n_teacher_table_rows",
                "n_teacher_unique_tics",
                "n_prior_evaluation_unique_tics",
                "n_total_excluded_unique_tics",
                "n_raw_adp_intersection",
                "n_adp_without_raw",
                "n_qa_eligible_after_teacher_exclusion",
                "n_qa_eligible_after_all_exclusions",
                "n_selected",
                "n_selected_unique_tics",
                "n_cells",
                "cell_support_min",
                "cell_support_max",
                "tmag_counts",
                "qa_reason_counts",
                "host_overlap_audits",
            )
        },
        "bls": {
            "anchor_aperture": "DET_FLUX_ADP_SML",
            "supplemental_aperture": "DET_FLUX_ADP",
            "n_trial_periods": 200_000,
            "period_min_d": 0.12,
            "period_max_policy": "min(15 d, 0.45 * observed baseline)",
            "duration_min": [3, 4, 5, 6, 8, 10, 13, 16, 20, 30],
            "saved_peaks_per_aperture": 10,
            "truth_period_tolerance": 0.02,
            "truth_harmonic_factors": [0.25, 1 / 3, 0.5, 1, 2, 3, 4],
            "minimum_transit_window_overlap": 0.5,
            "publication_top_k": 5,
        },
        "teacher": teacher,
        "injection_product_audit": product_audit,
        "final_summary": final,
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
