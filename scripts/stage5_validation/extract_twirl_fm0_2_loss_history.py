#!/usr/bin/env python3
"""Extract and validate the compact FM0.2 loss history from a milestone checkpoint."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import torch


EXPECTED_FIELDS = (
    "step",
    "total",
    "reconstruction_first",
    "reconstruction_second",
    "reconstruction_mean",
    "reconstruction",
    "vicreg_invariance",
    "vicreg_variance",
    "vicreg_covariance",
    "vicreg",
    "vicreg_weighted",
    "embedding_projection_gradient_norm",
    "learning_rate",
    "window_draws_seen",
    "masked_views_seen",
)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--checkpoint-sha256", required=True)
    parser.add_argument("--expected-step", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _close(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=2.0e-6, abs_tol=2.0e-9)


def _validated_history(checkpoint: dict[str, Any], expected_step: int) -> list[dict[str, float]]:
    progress = checkpoint.get("progress")
    if not isinstance(progress, dict) or int(progress.get("global_step", -1)) != expected_step:
        raise ValueError("checkpoint global step differs")
    raw_history = checkpoint.get("loss_history")
    if not isinstance(raw_history, list) or len(raw_history) != expected_step:
        raise ValueError("loss history length differs from the expected step")

    rows: list[dict[str, float]] = []
    for expected_row_step, raw in enumerate(raw_history, start=1):
        if not isinstance(raw, dict) or not set(EXPECTED_FIELDS).issubset(raw):
            raise ValueError(f"loss-history fields differ at step {expected_row_step}")
        row = {field: float(raw[field]) for field in EXPECTED_FIELDS}
        if int(row["step"]) != expected_row_step:
            raise ValueError("loss-history steps are not contiguous")
        if not all(math.isfinite(value) for value in row.values()):
            raise ValueError(f"nonfinite loss-history value at step {expected_row_step}")
        if not _close(row["reconstruction"], row["reconstruction_first"]):
            raise ValueError(f"optimized reconstruction identity failed at step {expected_row_step}")
        if not _close(
            row["reconstruction_mean"],
            0.5 * (row["reconstruction_first"] + row["reconstruction_second"]),
        ):
            raise ValueError(f"reconstruction mean identity failed at step {expected_row_step}")
        if not _close(row["vicreg_weighted"], 0.0002 * row["vicreg"]):
            raise ValueError(f"weighted VICReg identity failed at step {expected_row_step}")
        expected_vicreg = (
            25.0 * row["vicreg_invariance"]
            + 25.0 * row["vicreg_variance"]
            + row["vicreg_covariance"]
        )
        if not _close(row["vicreg"], expected_vicreg):
            raise ValueError(f"VICReg component identity failed at step {expected_row_step}")
        if not _close(
            row["total"], row["reconstruction_first"] + row["vicreg_weighted"]
        ):
            raise ValueError(f"total-loss identity failed at step {expected_row_step}")
        rows.append(row)
    return rows


def main() -> int:
    args = _parser().parse_args()
    observed_sha256 = _sha256(args.checkpoint)
    if observed_sha256 != args.checkpoint_sha256:
        raise ValueError("checkpoint SHA-256 differs")
    checkpoint = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    if not isinstance(checkpoint, dict):
        raise ValueError("checkpoint root is not a mapping")
    rows = _validated_history(checkpoint, args.expected_step)
    payload = {
        "schema_version": "twirl_fm0_2_loss_history_extract_v1",
        "checkpoint_path": str(args.checkpoint),
        "checkpoint_sha256": observed_sha256,
        "global_step": args.expected_step,
        "loss_history": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(args.output)
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
