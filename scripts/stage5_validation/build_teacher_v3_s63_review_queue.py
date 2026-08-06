#!/usr/bin/env python3
"""Build the frozen, blinded S63 Teacher-v3 morphology review queue."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.vetting.teacher_v3_prospective import (  # noqa: E402
    COHORT_CONTRACT_VERSION,
    build_s63_prospective_review_queue,
    file_sha256,
    load_prospective_plan,
    load_selection_policy,
    read_tic_inventory,
    validate_git_sha,
    validate_queue_launch_bindings,
    validate_sha256,
    write_s63_prospective_review_queue,
)
from twirl.vetting.teacher_v3_inference import (  # noqa: E402
    verify_execution_checkout,
)


def _read_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".parquet":
        return pd.read_parquet(path)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path, low_memory=False)
    raise ValueError(f"unsupported table format: {path}")


def _bound_table(path: Path, expected_sha256: str) -> pd.DataFrame:
    expected = validate_sha256(expected_sha256, context=f"expected hash for {path}")
    before = file_sha256(path)
    if before != expected:
        raise ValueError(f"SHA-256 mismatch for {path}")
    frame = _read_table(path)
    if file_sha256(path) != before:
        raise RuntimeError(f"artifact changed while it was read: {path}")
    return frame


def _bound_json(path: Path, expected_sha256: str) -> dict[str, object]:
    expected = validate_sha256(expected_sha256, context=f"expected hash for {path}")
    before = file_sha256(path)
    if before != expected:
        raise ValueError(f"SHA-256 mismatch for {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if file_sha256(path) != before:
        raise RuntimeError(f"artifact changed while it was read: {path}")
    if not isinstance(payload, dict):
        raise ValueError(f"JSON artifact must be an object: {path}")
    return payload


def _validate_output_directories(
    *,
    repo: Path,
    public_dir: Path,
    hidden_dir: Path,
) -> tuple[Path, Path]:
    if not public_dir.is_absolute() or not hidden_dir.is_absolute():
        raise ValueError("queue output directories must be absolute paths")
    repo = Path(repo).resolve(strict=True)
    public = public_dir.expanduser().resolve(strict=False)
    hidden = hidden_dir.expanduser().resolve(strict=False)
    for name, path in (("public", public), ("hidden", hidden)):
        if path == repo or repo in path.parents:
            raise ValueError(f"{name} queue output must be outside the checkout")
        if path.exists():
            raise FileExistsError(f"{name} queue output must be a fresh path: {path}")
    if public == hidden:
        raise ValueError("public and hidden queue outputs must be distinct")
    if public in hidden.parents or hidden in public.parents:
        raise ValueError("public and hidden queue outputs must not be nested")
    return public, hidden


def _validate_candidate_score_equality(
    scores: pd.DataFrame, candidates: pd.DataFrame
) -> None:
    identity_columns = ("tic", "sector", "review_id")
    missing = sorted(
        set(identity_columns) - set(scores)
        | (set(identity_columns) - set(candidates))
    )
    if missing:
        raise KeyError(f"score/candidate identity columns are missing: {missing}")
    missing_score_columns = sorted(set(candidates.columns) - set(scores.columns))
    if missing_score_columns:
        raise KeyError(
            "score table dropped candidate-derived columns: "
            f"{missing_score_columns}"
        )
    if scores["tic"].duplicated().any() or candidates["tic"].duplicated().any():
        raise ValueError("score/candidate tables must contain one row per TIC")
    left = scores.sort_values("tic", kind="stable").reset_index(drop=True)
    right = candidates.sort_values("tic", kind="stable").reset_index(drop=True)
    if len(left) != len(right):
        raise ValueError("score/candidate row counts differ")
    for column in candidates.columns:
        candidate_values = right[column]
        score_values = left[column]
        candidate_numeric = pd.to_numeric(candidate_values, errors="coerce")
        score_numeric = pd.to_numeric(score_values, errors="coerce")
        candidate_nonempty = candidate_values.notna() & candidate_values.astype(
            str
        ).str.strip().ne("")
        score_nonempty = score_values.notna() & score_values.astype(str).str.strip().ne(
            ""
        )
        numeric_candidate_column = bool(
            candidate_nonempty.any()
            and not (candidate_nonempty & candidate_numeric.isna()).any()
        )
        if numeric_candidate_column:
            if (score_nonempty & score_numeric.isna()).any() or not np.array_equal(
                candidate_numeric.to_numpy(dtype=float),
                score_numeric.to_numpy(dtype=float),
                equal_nan=True,
            ):
                raise ValueError(
                    f"score/candidate candidate-derived field {column} differs"
                )
            continue
        candidate_text = candidate_values.fillna("<NULL>").astype(str)
        score_text = score_values.fillna("<NULL>").astype(str)
        if not candidate_text.equals(score_text):
            raise ValueError(
                f"score/candidate candidate-derived field {column} differs"
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--expected-plan-sha256", required=True)
    parser.add_argument("--selection-policy", type=Path, required=True)
    parser.add_argument("--expected-selection-policy-sha256", required=True)
    parser.add_argument("--launch-manifest", type=Path, required=True)
    parser.add_argument("--expected-launch-manifest-sha256", required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--expected-git-sha", required=True)
    parser.add_argument("--producer-git-sha", required=True)
    parser.add_argument("--scoring-summary", type=Path, required=True)
    parser.add_argument("--expected-scoring-summary-sha256", required=True)
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--expected-scores-sha256", required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--expected-candidates-sha256", required=True)
    parser.add_argument("--primary-cohort", type=Path, required=True)
    parser.add_argument("--expected-primary-cohort-sha256", required=True)
    parser.add_argument("--repeated-host-cohort", type=Path, required=True)
    parser.add_argument("--expected-repeated-host-cohort-sha256", required=True)
    parser.add_argument("--cohort-summary", type=Path, required=True)
    parser.add_argument("--expected-cohort-summary-sha256", required=True)
    parser.add_argument("--public-out-dir", type=Path, required=True)
    parser.add_argument("--hidden-out-dir", type=Path, required=True)
    args = parser.parse_args()

    plan, plan_sha256 = load_prospective_plan(
        args.plan,
        expected_sha256=args.expected_plan_sha256,
    )
    selection_policy, selection_policy_sha256 = load_selection_policy(
        args.selection_policy,
        expected_sha256=args.expected_selection_policy_sha256,
    )
    launch_manifest = _bound_json(
        args.launch_manifest,
        args.expected_launch_manifest_sha256,
    )
    git_audit = verify_execution_checkout(
        repo=args.repo,
        expected_git_sha=args.expected_git_sha,
        launch_git=launch_manifest.get("git"),
    )
    producer_git_sha = validate_git_sha(
        args.producer_git_sha, context="queue producer Git SHA"
    )
    if producer_git_sha != git_audit["observed_git_sha"]:
        raise ValueError(
            "queue producer Git SHA differs from the audited execution checkout"
        )
    launch_producer_git_sha = validate_git_sha(
        launch_manifest.get("producer_git_sha"),
        context="launch producer Git SHA",
    )
    if launch_producer_git_sha != git_audit["launch_git_sha"]:
        raise ValueError("launch producer Git SHA differs from launch git.sha")
    public_out_dir, hidden_out_dir = _validate_output_directories(
        repo=Path(git_audit["repo"]),
        public_dir=args.public_out_dir,
        hidden_dir=args.hidden_out_dir,
    )
    candidates = _bound_table(args.candidates, args.expected_candidates_sha256)
    primary = read_tic_inventory(
        args.primary_cohort,
        expected_sha256=args.expected_primary_cohort_sha256,
    )
    repeated = read_tic_inventory(
        args.repeated_host_cohort,
        expected_sha256=args.expected_repeated_host_cohort_sha256,
    )
    cohort_summary = _bound_json(
        args.cohort_summary,
        args.expected_cohort_summary_sha256,
    )
    if cohort_summary.get("schema_version") != COHORT_CONTRACT_VERSION:
        raise ValueError("cohort summary has an incompatible schema")
    checks = cohort_summary.get("partition_checks")
    if not isinstance(checks, dict) or not checks or not all(value is True for value in checks.values()):
        raise ValueError("cohort summary did not pass all partition checks")
    outputs = cohort_summary.get("outputs")
    if not isinstance(outputs, dict):
        raise ValueError("cohort summary is missing output hashes")
    expected_cohort_hashes = {
        "primary_cohort": validate_sha256(
            args.expected_primary_cohort_sha256, context="primary cohort"
        ),
        "repeated_host_cohort": validate_sha256(
            args.expected_repeated_host_cohort_sha256, context="repeated-host cohort"
        ),
    }
    for name, expected in expected_cohort_hashes.items():
        entry = outputs.get(name)
        if not isinstance(entry, dict) or entry.get("sha256") != expected:
            raise ValueError(f"cohort summary does not bind {name}")
    source_hashes = cohort_summary.get("source_hashes")
    if not isinstance(source_hashes, dict):
        raise ValueError("cohort summary is missing source hashes")
    training = plan.get("frozen_training_identity", {})
    reservation = plan.get("s63_identity_reservation", {})
    if source_hashes.get("teacher_v3_corpus_sha256") != training.get("morphology_corpus_sha256"):
        raise ValueError("cohort summary corpus hash differs from the prospective plan")
    if source_hashes.get("reserved_tics_sha256") != reservation.get("reserved_tics_sha256"):
        raise ValueError("cohort summary reservation hash differs from the prospective plan")
    hashes = {
        "prospective_plan_sha256": plan_sha256,
        "selection_policy_sha256": selection_policy_sha256,
        "launch_manifest_sha256": args.expected_launch_manifest_sha256,
        "reserved_tics_sha256": source_hashes["reserved_tics_sha256"],
        "teacher_v3_corpus_sha256": source_hashes["teacher_v3_corpus_sha256"],
        "model_ready_allowlist_sha256": source_hashes["model_ready_allowlist_sha256"],
        "primary_cohort_sha256": args.expected_primary_cohort_sha256,
        "repeated_host_cohort_sha256": args.expected_repeated_host_cohort_sha256,
        "cohort_summary_sha256": args.expected_cohort_summary_sha256,
        "candidate_table_sha256": args.expected_candidates_sha256,
        "teacher_score_table_sha256": args.expected_scores_sha256,
        "teacher_score_summary_sha256": args.expected_scoring_summary_sha256,
    }
    scoring_summary = _bound_json(
        args.scoring_summary,
        args.expected_scoring_summary_sha256,
    )
    launch_binding = validate_queue_launch_bindings(
        launch_manifest=launch_manifest,
        prospective_plan=plan,
        scoring_summary=scoring_summary,
        artifact_hashes=hashes,
        git_audit=git_audit,
    )
    scores = _bound_table(args.scores, args.expected_scores_sha256)
    _validate_candidate_score_equality(scores, candidates)
    queue, hidden, public_summary, hidden_summary = build_s63_prospective_review_queue(
        scores,
        primary_tics=primary,
        repeated_host_tics=repeated,
        artifact_hashes=hashes,
        selection_policy=selection_policy,
        launch_binding=launch_binding,
        producer_git_sha=producer_git_sha,
    )
    frozen_paths = {
        args.plan: args.expected_plan_sha256,
        args.selection_policy: args.expected_selection_policy_sha256,
        args.launch_manifest: args.expected_launch_manifest_sha256,
        args.scoring_summary: args.expected_scoring_summary_sha256,
        args.scores: args.expected_scores_sha256,
        args.candidates: args.expected_candidates_sha256,
        args.primary_cohort: args.expected_primary_cohort_sha256,
        args.repeated_host_cohort: args.expected_repeated_host_cohort_sha256,
        args.cohort_summary: args.expected_cohort_summary_sha256,
    }
    for path, expected in frozen_paths.items():
        if file_sha256(path) != validate_sha256(expected, context=f"frozen {path}"):
            raise RuntimeError(f"queue input changed before publication: {path}")
    final_git_audit = verify_execution_checkout(
        repo=args.repo,
        expected_git_sha=args.expected_git_sha,
        launch_git=launch_manifest.get("git"),
    )
    if final_git_audit != git_audit:
        raise RuntimeError("execution checkout audit changed before queue publication")
    payload = write_s63_prospective_review_queue(
        queue,
        hidden,
        public_summary,
        hidden_summary,
        public_dir=public_out_dir,
        hidden_dir=hidden_out_dir,
    )
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
