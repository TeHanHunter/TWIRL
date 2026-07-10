#!/usr/bin/env python3
"""Build the fixed S56 real-label adjudication queue."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.vetting.adjudication import build_real_adjudication_queue  # noqa: E402


TEACHER_ROOT = REPO_ROOT / "reports/stage5_validation/s56_recovery50_teacher_queue_2k"
PILOT_ROOT = REPO_ROOT / "reports/stage5_validation/s56_eb_miner_adp_only/review_queue_eb_priority_pilot100"
DEPRECATED_ROOT = REPO_ROOT / "reports/stage5_validation/s56_eb_miner/review_queue_eb_priority_500"
DEFAULT_OUT = REPO_ROOT / "reports/stage5_validation/s56_label_adjudication_real343"


def _write_readme(out_dir: Path) -> None:
    text = """# S56 Real-Label Adjudication

This is a blinded, real-only consistency review. The queue contains 323 unique
sources and 20 hidden repeats. Previous labels, source cohort, repeat identity,
and model scores are not shown in the browser.

## Labels

- **Planet-like**: compact isolated transit-like event without convincing
  secondary or binary morphology.
- **Eclipse/contact**: discrete secondary, alternating depths, or convincing
  detached/contact-binary eclipse morphology.
- **Smooth variable**: continuous or sinusoidal modulation without discrete
  eclipses, including possible non-eclipsing PCEBs.
- **Systematic/artifact**: window-edge leakage, cadence aliases, detrending
  structure, or another recognizable artifact.
- **Flat/no signal**: no obvious useful signal. This is not an ambiguity label.
- **Broad isolated dip**: broad primary-like event without convincing secondary
  or continuous variability.
- **Skip** is only for broken evidence.

Eclipse/contact takes precedence over Broad isolated dip when a secondary,
odd/even difference, or contact-binary pattern is credible. Use Smooth variable
when the modulation is continuous rather than eclipse-like.

## Period Control

Leave the period at `P` unless another harmonic visibly gives the coherent
morphology. The choices are `P/4`, `P/2`, `P`, `2P`, `4P`, and `Unresolved`.
The selected factor is saved atomically with the label; notes are optional.

## Start

```bash
scripts/stage5_validation/run_s56_real_label_adjudication_app_local.sh
```

Open `http://127.0.0.1:5004/`.
"""
    (out_dir / "README.md").write_text(text)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base-training-table",
        type=Path,
        default=TEACHER_ROOT / "human_training_table/human_vetting_training_table.csv",
    )
    parser.add_argument(
        "--active-training-table",
        type=Path,
        default=TEACHER_ROOT / "human_training_table_adp_only/human_vetting_training_table.csv",
    )
    parser.add_argument("--eb-pilot-queue", type=Path, default=PILOT_ROOT / "review_queue_eb_priority_100.csv")
    parser.add_argument("--eb-pilot-labels", type=Path, default=PILOT_ROOT / "human_labels_vetted.csv")
    parser.add_argument(
        "--deprecated-eb-queue", type=Path, default=DEPRECATED_ROOT / "review_queue_eb_priority_500.csv"
    )
    parser.add_argument("--deprecated-eb-labels", type=Path, default=DEPRECATED_ROOT / "human_labels_vetted.csv")
    parser.add_argument(
        "--current-adp-bls",
        type=Path,
        default=REPO_ROOT / "reports/stage5_validation/s56_adp_real_bls_peaks/real_adp_bls_peaks.parquet",
    )
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--seed", type=int, default=56)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    summary = build_real_adjudication_queue(
        base_training_table=args.base_training_table,
        active_training_table=args.active_training_table,
        eb_pilot_queue=args.eb_pilot_queue,
        eb_pilot_labels=args.eb_pilot_labels,
        deprecated_eb_queue=args.deprecated_eb_queue,
        deprecated_eb_labels=args.deprecated_eb_labels,
        current_adp_bls=args.current_adp_bls,
        out_dir=args.out_dir,
        seed=args.seed,
    )
    _write_readme(args.out_dir)
    print(json.dumps(summary, indent=2, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
