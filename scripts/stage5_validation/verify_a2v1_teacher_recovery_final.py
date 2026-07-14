#!/usr/bin/env python3
"""Verify final sector A2v1 BLS and Teacher-v1 recovery artifacts."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import pandas as pd

from twirl.injections.a2v1_recovery import (
    load_recovery_config,
    numeric_arrays_match_with_ulp_budget,
    schedule_contract,
)
from twirl.vetting.recovery50_teacher import leakage_columns


CANDIDATE_CSV_ULP_BUDGET = {
    # Candidate shards pass through two portable CSV round trips. The observed
    # decimal conversion is bounded to seven spacing units for short periods
    # and two for BJD epochs; durations are exact integer-minute grid values.
    "period_d": 8,
    "t0_bjd": 2,
    "duration_min": 0,
}


def _json(path: Path) -> dict:
    return json.loads(path.read_text())


def _bool_series(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values.fillna(False)
    return (
        values.fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
        .isin({"1", "1.0", "true", "t", "yes", "y"})
    )


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
    schedule_summary = _json(root / "schedule/schedule_summary.json")
    product_audit = _json(final / "injection_product_audit.json")
    summary = _json(final / "summary.json")
    run_manifest = _json(final / "run_manifest.json")
    teacher_summary = _json(full / "teacher_scores.summary.json")
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
    require(
        int(input_verification.get("compact_rebuild_parity", {}).get("n_mismatched_targets", -1))
        == 0,
        "rebuilt compact product has target-level mismatches",
    )
    expected_hash_names = {
        f"hlsp_s{config.sector:04d}_A2v1_tree",
        "human_vetting_training_table_adjudicated.csv",
        f"s{config.sector}_A2v1_adp_pair.h5",
        f"s{config.sector}_A2v1_tglc_raw_sources.h5",
        f"s{config.sector}_A2v1_validation_full_schema.json",
    }
    input_hashes = input_verification.get("hashes", {})
    require(
        expected_hash_names.issubset(input_hashes),
        "input verification is missing locked source hashes",
    )
    require(
        all(
            len(str(input_hashes.get(name, ""))) == 64
            and all(
                value in "0123456789abcdef"
                for value in str(input_hashes.get(name, "")).lower()
            )
            for name in expected_hash_names
        ),
        "input verification contains an invalid source SHA256",
    )
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
    shard_support = schedule.groupby("shard_index").agg(
        n_injections=("injection_id", "size"),
        n_period_bins=("grid_period_bin", "nunique"),
        n_radius_bins=("grid_radius_bin", "nunique"),
    )
    require(len(shard_support) == config.n_shards, "schedule shard count mismatch")
    require(
        shard_support["n_injections"].eq(config.rows_per_shard).all(),
        "schedule rows are not balanced across shards",
    )
    if config.shard_assignment == "balanced_random":
        require(
            shard_support["n_period_bins"].eq(config.period_bins).all()
            and shard_support["n_radius_bins"].eq(config.radius_bins).all(),
            "balanced schedule shards do not span both grid axes",
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
    host_qa = pd.read_csv(root / "schedule/host_qa.csv.gz", low_memory=False)
    eligible_bright = host_qa.loc[
        _bool_series(host_qa["eligible"])
        & pd.to_numeric(host_qa["tessmag"], errors="coerce").lt(19.0),
        "tic",
    ]
    require(
        set(pd.to_numeric(eligible_bright, errors="raise").astype(int)).issubset(
            set(schedule["tic"].astype(int))
        ),
        "schedule does not include every eligible Tmag < 19 host",
    )
    if config.sector > 56:
        require(
            bool(schedule_summary.get("host_overlap_audits")),
            "later-sector schedule has no prior-sector host-overlap audit",
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
    small_top5 = peaks.loc[
        peaks["aperture"].astype(str).eq("DET_FLUX_ADP_SML")
        & pd.to_numeric(peaks["peak_rank"], errors="coerce").between(1, 5)
    ].sort_values(["injection_id", "peak_rank"], kind="stable")
    candidate_top5 = candidates.sort_values(
        ["injection_id", "rep_peak_rank"], kind="stable"
    )
    require(
        small_top5["injection_id"].astype(str).tolist()
        == candidate_top5["injection_id"].astype(str).tolist()
        and pd.to_numeric(small_top5["peak_rank"], errors="coerce").tolist()
        == pd.to_numeric(candidate_top5["rep_peak_rank"], errors="coerce").tolist(),
        "Teacher candidate identities differ from the stored ADP-small top five",
    )
    for column, max_ulps in CANDIDATE_CSV_ULP_BUDGET.items():
        require(
            numeric_arrays_match_with_ulp_budget(
                pd.to_numeric(small_top5[column], errors="coerce"),
                pd.to_numeric(candidate_top5[column], errors="coerce"),
                max_ulps=max_ulps,
            ),
            f"Teacher candidate {column} differs from the stored ADP-small top five",
        )
    require(
        _bool_series(small_top5["is_injected_signal_peak"]).tolist()
        == _bool_series(candidate_top5["is_injected_signal_peak"]).tolist(),
        "Teacher candidate truth-match audit flags differ from the stored top five",
    )
    require(len(scores) == expected_candidate_rows, "Teacher score row count mismatch")
    require(scores["review_id"].is_unique, "Teacher score review IDs are duplicated")
    require(
        candidates["review_id"].astype(str).tolist()
        == scores["review_id"].astype(str).tolist(),
        "Teacher scores do not preserve the candidate set and order",
    )
    metadata_by_fold = teacher_summary.get("metadata_columns_by_fold", [])
    require(
        len(metadata_by_fold) == 5,
        "Teacher summary does not contain five metadata feature contracts",
    )
    metadata_columns = {
        str(column)
        for columns in metadata_by_fold
        for column in columns
    }
    require(
        not leakage_columns(metadata_columns),
        "Teacher checkpoint metadata contains truth or provenance leakage",
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
        "candidate_csv_ulp_budget": CANDIDATE_CSV_ULP_BUDGET,
        "failures": failures,
        "passed": not failures,
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0 if payload["passed"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
