#!/usr/bin/env python3
"""Verify locked/external Teacher-v2 recovery tables and figure artifacts."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--role-summary", type=Path, required=True)
    parser.add_argument("--expected-sector", type=int, required=True)
    parser.add_argument("--minimum-injections", type=int, required=True)
    parser.add_argument("--out-json", type=Path)
    args = parser.parse_args()

    root = args.root
    summary = json.loads((root / "summary.json").read_text())
    roles = json.loads(args.role_summary.read_text())
    outcomes = pd.read_parquet(root / "injection_recovery_outcomes.parquet")
    comparison = pd.read_csv(root / "workload_matched_model_comparison.csv")
    failures: list[str] = []

    def require(condition: bool, message: str) -> None:
        if not condition:
            failures.append(message)

    require(int(summary.get("sector", -1)) == args.expected_sector, "sector mismatch")
    require(len(outcomes) >= args.minimum_injections, "injection support is too small")
    require(outcomes["injection_id"].is_unique, "injection outcomes are duplicated")
    require(outcomes["tic"].nunique() == len(outcomes), "external hosts are reused")
    require(
        int(outcomes["grid_cell_id"].nunique()) == 2500,
        "period-radius support is incomplete",
    )
    require(
        bool(summary.get("candidate_set_identity_verified")),
        "candidate identity was not verified",
    )
    frozen = summary.get("frozen_selection", {})
    require(bool(frozen.get("architecture_frozen")), "architecture was not frozen")
    require(bool(frozen.get("threshold_frozen")), "threshold was not frozen")
    require(
        float(frozen.get("frozen_real_tic_workload_fraction", 1.0)) <= 0.05,
        "frozen real-TIC workload exceeds five percent",
    )
    require(
        set(comparison["model"].astype(str))
        == {"teacher_v2", "teacher_v1", "metadata_only"},
        "workload comparison does not contain all three models",
    )
    if args.expected_sector == 57:
        external = roles.get("s57_external", {})
        require(bool(external.get("passed")), "S57 role support gate failed")
        require(
            int(external.get("n_external", -1)) == len(outcomes),
            "S57 role and outcome counts differ",
        )

    required = [
        root / "bls_recovery_at_1_3_5.csv",
        root / "teacher_retention_by_tmag.csv",
        root / "period_radius_recovery_unsmoothed.csv",
        root / "teacher_v2_strict_planet_truth_matches.parquet",
    ]
    for subdir in ("bls_top5_recovery_map", "teacher_v2_compact_recovery_map"):
        required.extend(
            [
                root / subdir / "period_radius_empirical_recovery_publication.png",
                root / subdir / "period_radius_empirical_recovery_publication.pdf",
                root
                / subdir
                / "period_radius_empirical_recovery_publication_grid.csv",
            ]
        )
    missing = [str(path) for path in required if not path.is_file() or path.stat().st_size == 0]
    require(not missing, f"missing or empty artifacts: {missing}")
    result = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "root": str(root),
        "expected_sector": int(args.expected_sector),
        "n_injections": int(len(outcomes)),
        "artifact_integrity_passed": not failures,
        "scientific_acceptance_passed": bool(summary.get("acceptance_passed")),
        "scientific_acceptance": summary.get("acceptance", {}),
        "failures": failures,
    }
    out_json = args.out_json or (root / "verification.json")
    out_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
