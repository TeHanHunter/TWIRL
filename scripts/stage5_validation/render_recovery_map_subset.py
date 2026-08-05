#!/usr/bin/env python3
"""Render selected Tmag panels from an existing recovery-outcomes table."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from plot_s56_duration_aware_recovery import (  # noqa: E402
    plot_publication_period_radius_recovery_map,
)


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outcomes", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument(
        "--panels",
        choices=("bright-two", "all"),
        default="bright-two",
        help="Use the two Tmag < 18 panels or the original four-panel layout.",
    )
    parser.add_argument(
        "--recovered-col",
        default="bls_top5_recovered",
        help="Boolean outcome column used to construct the recovery surface.",
    )
    parser.add_argument(
        "--colorbar-label",
        default="Kernel-smoothed ADP BLS top-5 recovery fraction",
    )
    parser.add_argument(
        "--output-stem",
        default=None,
        help="Output filename stem; defaults to a layout-specific descriptive name.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    if args.outcomes.suffix.lower() == ".parquet":
        outcomes = pd.read_parquet(args.outcomes)
    else:
        outcomes = pd.read_csv(args.outcomes, low_memory=False)
    panel_indices = (0, 1) if args.panels == "bright-two" else (0, 1, 2, 3)
    output_stem = args.output_stem or (
        "period_radius_empirical_recovery_publication_tmag_lt18"
        if args.panels == "bright-two"
        else "period_radius_empirical_recovery_publication"
    )
    paths = plot_publication_period_radius_recovery_map(
        outcomes,
        args.out_dir,
        recovered_col=args.recovered_col,
        colorbar_label=args.colorbar_label,
        panel_indices=panel_indices,
        output_stem=output_stem,
    )
    print(json.dumps(paths, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
