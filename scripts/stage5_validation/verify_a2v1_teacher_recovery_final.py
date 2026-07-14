#!/usr/bin/env python3
"""Verify final sector A2v1 BLS and Teacher-v1 recovery artifacts."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import pandas as pd

from twirl.injections.a2v1_recovery import load_recovery_config, schedule_contract


def _json(path: Path) -> dict:
    return json.loads(path.read_text())


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--recovery-root", type=Path, required=True)
    parser.add_argument("--input-verification", type=Path, required=True)
    parser.add_argument("--teacher-table", type=Path, required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    args = parser.parse_args()

    config = load_recovery_config(args.config)
    root = args.recovery_root
    full = root / "full"
    final = full / "final"
    schedule = pd.read_parquet(root / "schedule/injection_schedule.parquet")
    teacher = pd.read_csv(args.teacher_table, low_memory=False)
    peaks = pd.read_parquet(final / "adp_bls_top10_peaks.parquet")
    candidates = pd.read_csv(
        full / "teacher_inputs/candidates_top5.csv", low_memory=False
    )
    scores = pd.read_parquet(full / "teacher_scores.parquet")
    outcomes = pd.read_parquet(final / "injection_recovery_outcomes.parquet")
    input_verification = _json(args.input_verification)
    parity = _json(root / "schedule/adp_roundtrip_parity_summary.json")
    product_audit = _json(final / "injection_product_audit.json")
    summary = _json(final / "summary.json")
    run_manifest = _json(final / "run_manifest.json")
    failures: list[str] = []

    def require(condition: bool, message: str) -> None:
        if not condition:
            failures.append(message)

    expected = config.n_injections
    require(bool(input_verification.get("passed")), "input verification failed")
    require(
        int(input_verification.get("sector", -1)) == config.sector,
        "input verification sector mismatch",
    )
    require(bool(parity.get("passed")), "ADP no-injection parity failed")
    require(bool(product_audit.get("passed")), "full injection-product audit failed")
    require(len(schedule) == expected, "schedule row count mismatch")
    require(schedule["tic"].nunique() == expected, "schedule reuses host TICs")
    require(
        schedule["injection_id"].nunique() == expected,
        "schedule reuses injection IDs",
    )
    require(
        pd.to_numeric(schedule["sector"], errors="coerce").eq(config.sector).all(),
        "schedule sector mismatch",
    )
    require(
        schedule["schedule_contract"].astype(str).eq(
            schedule_contract(config.sector)
        ).all(),
        "schedule contract mismatch",
    )
    support = schedule.groupby(["grid_period_bin", "grid_radius_bin"]).size()
    require(
        len(support) == config.period_bins * config.radius_bins,
        "schedule grid is incomplete",
    )
    require(
        support.eq(config.repeats_per_cell).all(),
        "schedule has incorrect per-cell support",
    )
    teacher_tics = set(
        pd.to_numeric(teacher["tic"], errors="coerce").dropna().astype(int)
    )
    require(
        not (set(schedule["tic"].astype(int)) & teacher_tics),
        "schedule overlaps Teacher-v1 hosts",
    )
    require(
        bool(schedule["evaluation_only"].all())
        and bool(schedule["teacher_training_excluded"].all()),
        "schedule is not locked evaluation-only",
    )
    require(
        int(product_audit.get("n_groups", -1)) == expected,
        "product audit group count mismatch",
    )
    require(
        int(product_audit.get("n_shards", -1)) == config.n_shards,
        "product audit shard count mismatch",
    )
    expected_peak_rows = expected * 2 * 10
    require(len(peaks) == expected_peak_rows, "BLS peak row count mismatch")
    require(peaks["injection_id"].nunique() == expected, "BLS injection loss")
    require(
        set(peaks["aperture"].astype(str))
        == {"DET_FLUX_ADP_SML", "DET_FLUX_ADP"},
        "BLS table contains a non-ADP aperture",
    )
    peak_support = peaks.groupby(["injection_id", "aperture"]).size()
    require(peak_support.eq(10).all(), "BLS table does not have ten peaks/aperture")
    expected_candidate_rows = expected * 5
    require(
        len(candidates) == expected_candidate_rows,
        "Teacher candidate row count mismatch",
    )
    require(
        candidates["review_id"].is_unique,
        "Teacher candidate review IDs are duplicated",
    )
    require(
        candidates.groupby("injection_id").size().eq(5).all(),
        "Teacher candidate table does not have five candidates/injection",
    )
    require(len(scores) == expected_candidate_rows, "Teacher score row count mismatch")
    require(scores["review_id"].is_unique, "Teacher score review IDs are duplicated")
    require(
        candidates["review_id"].astype(str).tolist()
        == scores["review_id"].astype(str).tolist(),
        "Teacher scores do not preserve the candidate set and order",
    )
    require(len(outcomes) == expected, "final outcome row count mismatch")
    require(outcomes["injection_id"].is_unique, "final outcomes are duplicated")
    require(
        bool(summary.get("candidate_set_identity_verified")),
        "finalizer did not verify candidate-set identity",
    )
    require(
        int(summary.get("sector", -1)) == config.sector,
        "final summary sector mismatch",
    )
    require(
        run_manifest.get("workflow") == config.name,
        "run manifest workflow mismatch",
    )
    require(
        bool(run_manifest.get("evaluation_only"))
        and not bool(run_manifest.get("teacher_retraining_allowed")),
        "run manifest does not lock evaluation-only use",
    )
    require(
        bool(run_manifest.get("code", {}).get("tracked_code_clean")),
        "run manifest reports a dirty tracked checkout",
    )
    for field in ("producer_git_sha", "finalizer_git_sha"):
        sha = str(run_manifest.get("code", {}).get(field, "")).lower()
        require(
            len(sha) == 40 and all(value in "0123456789abcdef" for value in sha),
            f"run manifest has no valid {field}",
        )
    required_artifacts = (
        final / "bls_recovery_at_1_3_5.csv",
        final / "teacher_retention_by_tmag.csv",
        final / "period_radius_recovery_unsmoothed.csv",
        final / "teacher_v1_strict_planet_truth_matches.parquet",
        final
        / "bls_top5_recovery_map/period_radius_empirical_recovery_publication.png",
        final
        / "bls_top5_recovery_map/period_radius_empirical_recovery_publication.pdf",
        final
        / "bls_top5_recovery_map/period_radius_empirical_recovery_publication_grid.csv",
        final
        / "teacher_v1_preserve_recovery_map/period_radius_empirical_recovery_publication.png",
        final
        / "teacher_v1_preserve_recovery_map/period_radius_empirical_recovery_publication.pdf",
        final
        / "teacher_v1_preserve_recovery_map/period_radius_empirical_recovery_publication_grid.csv",
    )
    missing = [
        str(path)
        for path in required_artifacts
        if not path.is_file() or not path.stat().st_size
    ]
    require(not missing, f"missing or empty final artifacts: {missing}")

    payload = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "sector": config.sector,
        "workflow": config.name,
        "n_injections": len(schedule),
        "n_peak_rows": len(peaks),
        "n_teacher_candidate_rows": len(candidates),
        "n_teacher_score_rows": len(scores),
        "n_outcomes": len(outcomes),
        "n_required_artifacts": len(required_artifacts),
        "failures": failures,
        "passed": not failures,
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0 if payload["passed"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
