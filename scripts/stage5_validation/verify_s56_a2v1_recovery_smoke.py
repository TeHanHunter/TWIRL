#!/usr/bin/env python3
"""Gate the full S56 A2v1 recovery run on CPU and H200 smoke artifacts."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import pandas as pd

from twirl.vetting.harmonic_cnn import MORPHOLOGY_CLASSES


def _json(path: Path) -> dict:
    return json.loads(path.read_text())


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recovery-root", type=Path, required=True)
    parser.add_argument("--input-verification", type=Path, required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    parser.add_argument("--expect-injections", type=int, default=100)
    parser.add_argument("--expect-scores", type=int, default=128)
    args = parser.parse_args()

    smoke = args.recovery_root / "smoke100"
    paths = {
        "input_verification": args.input_verification,
        "parity": args.recovery_root / "schedule/adp_roundtrip_parity_summary.json",
        "injection": smoke / "injection_shards/injections_00_summary.json",
        "peak_verification": smoke / "peak_shards/peaks_00_verification.json",
        "native": smoke / "teacher_inputs/native_00.teacher_summary.json",
        "candidates": smoke / "teacher_inputs/candidates_00.csv",
        "scores": smoke / "teacher_scores_128.parquet",
    }
    failures: list[str] = []
    for name, path in paths.items():
        if not path.exists():
            failures.append(f"missing {name}: {path}")
    payloads = {}
    for name in (
        "input_verification",
        "parity",
        "injection",
        "peak_verification",
        "native",
    ):
        if paths[name].exists():
            payloads[name] = _json(paths[name])
    if payloads.get("input_verification", {}).get("passed") is not True:
        failures.append("input transfer verification did not pass")
    if payloads.get("parity", {}).get("passed") is not True:
        failures.append("ADP round-trip parity did not pass")
    if (
        int(payloads.get("injection", {}).get("n_injections", -1))
        != args.expect_injections
    ):
        failures.append("injection smoke does not contain exactly 100 injections")
    if payloads.get("peak_verification", {}).get("passed") is not True:
        failures.append("BLS peak smoke verification did not pass")
    native_passed = (
        payloads.get("native", {}).get("native_verification", {}).get("passed")
    )
    if native_passed is not True:
        failures.append("native Teacher-v1 input verification did not pass")

    score_summary: dict[str, object] = {}
    if paths["candidates"].exists() and paths["scores"].exists():
        candidates = pd.read_csv(paths["candidates"], low_memory=False).head(
            args.expect_scores
        )
        scores = pd.read_parquet(paths["scores"])
        if len(scores) != min(args.expect_scores, len(candidates)):
            failures.append("Teacher smoke has the wrong score-row count")
        if (
            scores["review_id"].astype(str).tolist()
            != candidates["review_id"].astype(str).tolist()
        ):
            failures.append("Teacher smoke changed candidate ordering or identity")
        probability_columns = [f"p_{label}" for label in MORPHOLOGY_CLASSES]
        probability_columns.append("p_preserve")
        missing = [column for column in probability_columns if column not in scores]
        if missing:
            failures.append(f"Teacher smoke is missing probabilities: {missing}")
        else:
            probabilities = scores[probability_columns].apply(
                pd.to_numeric, errors="coerce"
            )
            if not np.isfinite(probabilities.to_numpy()).all():
                failures.append("Teacher smoke contains non-finite probabilities")
            if ((probabilities < 0) | (probabilities > 1)).any().any():
                failures.append("Teacher smoke probabilities fall outside [0, 1]")
            morphology_sum = probabilities[
                [f"p_{label}" for label in MORPHOLOGY_CLASSES]
            ].sum(axis=1)
            if not np.allclose(morphology_sum, 1.0, atol=1.0e-5):
                failures.append(
                    "Teacher smoke morphology probabilities do not sum to one"
                )
        n_classes = int(scores["predicted_morphology"].astype(str).nunique())
        if n_classes < 2:
            failures.append("Teacher smoke predictions are degenerate")
        score_summary = {
            "n_scores": int(len(scores)),
            "n_predicted_classes": n_classes,
            "predicted_class_counts": {
                str(key): int(value)
                for key, value in scores["predicted_morphology"]
                .value_counts()
                .sort_index()
                .items()
            },
        }
    result = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "recovery_root": str(args.recovery_root),
        "artifacts": {key: str(value) for key, value in paths.items()},
        "score_summary": score_summary,
        "passed": not failures,
        "failures": failures,
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["passed"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
