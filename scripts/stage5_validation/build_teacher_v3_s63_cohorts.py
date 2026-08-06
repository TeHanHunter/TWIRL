#!/usr/bin/env python3
"""Build hash-bound disjoint/repeated S63 Teacher-v3 cohorts."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.vetting.teacher_v3_prospective import (  # noqa: E402
    load_prospective_plan,
    validate_git_sha,
    write_s63_prospective_cohorts,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--expected-plan-sha256", required=True)
    parser.add_argument("--model-ready-allowlist", type=Path, required=True)
    parser.add_argument("--expected-model-ready-allowlist-sha256", required=True)
    parser.add_argument("--teacher-v3-corpus", type=Path, required=True)
    parser.add_argument("--reserved-tics", type=Path, required=True)
    parser.add_argument("--producer-git-sha", required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    producer_git_sha = validate_git_sha(
        args.producer_git_sha, context="cohort producer Git SHA"
    )

    plan, _ = load_prospective_plan(
        args.plan,
        expected_sha256=args.expected_plan_sha256,
    )
    training = plan.get("frozen_training_identity", {})
    reservation = plan.get("s63_identity_reservation", {})
    if not isinstance(training, dict) or "morphology_corpus_sha256" not in training:
        raise ValueError("prospective plan lacks the frozen morphology-corpus hash")
    if not isinstance(reservation, dict) or "reserved_tics_sha256" not in reservation:
        raise ValueError("prospective plan lacks the frozen S63 reservation hash")
    source_hashes = {
        "model_ready_allowlist_sha256": args.expected_model_ready_allowlist_sha256,
        "teacher_v3_corpus_sha256": training["morphology_corpus_sha256"],
        "reserved_tics_sha256": reservation["reserved_tics_sha256"],
    }
    summary = write_s63_prospective_cohorts(
        model_ready_path=args.model_ready_allowlist,
        teacher_v3_corpus_path=args.teacher_v3_corpus,
        reserved_tics_path=args.reserved_tics,
        out_dir=args.out_dir,
        source_hashes=source_hashes,
        producer_git_sha=producer_git_sha,
        expected_reserved_count=reservation.get("n_requested_tics"),
        expected_corpus_count=training.get("n_corpus_tics"),
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
