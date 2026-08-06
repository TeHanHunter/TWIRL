#!/usr/bin/env python3
"""Apply the frozen Teacher-v3 ensemble to prospective real S63 candidates."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.vetting.teacher_v3_inference import (  # noqa: E402
    score_teacher_v3_s63_to_disk,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--expected-git-sha", required=True)
    parser.add_argument("--launch-manifest", type=Path, required=True)
    parser.add_argument("--launch-manifest-sha256", required=True)
    parser.add_argument("--preregistered-contract", type=Path, required=True)
    parser.add_argument("--selection-policy", type=Path, required=True)
    parser.add_argument("--release-summary", type=Path, required=True)
    parser.add_argument("--release-summary-sha256", required=True)
    parser.add_argument("--pretest-freeze", type=Path, required=True)
    parser.add_argument("--pretest-freeze-sha256", required=True)
    parser.add_argument("--checkpoint-manifest", type=Path, required=True)
    parser.add_argument("--checkpoint-manifest-sha256", required=True)
    parser.add_argument("--calibration", type=Path, required=True)
    parser.add_argument("--calibration-sha256", required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--candidates-sha256", required=True)
    parser.add_argument("--s63-reserved-tics", type=Path, required=True)
    parser.add_argument("--s63-reserved-tics-sha256", required=True)
    parser.add_argument("--native-h5", type=Path, required=True)
    parser.add_argument("--native-h5-sha256", required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument(
        "--allow-cpu",
        action="store_true",
        help="Allow CPU inference for controlled smoke tests; production defaults to CUDA.",
    )
    args = parser.parse_args()
    summary = score_teacher_v3_s63_to_disk(
        repo=args.repo,
        expected_git_sha=args.expected_git_sha,
        launch_manifest_path=args.launch_manifest,
        expected_launch_manifest_sha256=args.launch_manifest_sha256,
        preregistered_contract_path=args.preregistered_contract,
        selection_policy_path=args.selection_policy,
        release_summary_path=args.release_summary,
        expected_release_summary_sha256=args.release_summary_sha256,
        pretest_freeze_path=args.pretest_freeze,
        expected_pretest_freeze_sha256=args.pretest_freeze_sha256,
        checkpoint_manifest_path=args.checkpoint_manifest,
        expected_checkpoint_manifest_sha256=args.checkpoint_manifest_sha256,
        calibration_path=args.calibration,
        expected_calibration_sha256=args.calibration_sha256,
        candidates_path=args.candidates,
        expected_candidates_sha256=args.candidates_sha256,
        reserved_tics_path=args.s63_reserved_tics,
        expected_reserved_tics_sha256=args.s63_reserved_tics_sha256,
        native_h5=args.native_h5,
        expected_native_h5_sha256=args.native_h5_sha256,
        out_dir=args.out_dir,
        batch_size=args.batch_size,
        workers=args.workers,
        require_cuda=not args.allow_cpu,
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
