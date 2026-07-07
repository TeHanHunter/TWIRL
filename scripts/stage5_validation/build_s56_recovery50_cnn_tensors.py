#!/usr/bin/env python3
"""Build candidate-observable CNN tensors for the S56 recovery50 queue."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.vetting.recovery50_cnn import (  # noqa: E402
    TensorConfig,
    build_recovery50_cnn_tensors,
)
from twirl.vetting.recovery50_teacher import DEFAULT_APERTURES, json_default  # noqa: E402

DEFAULT_ROOT = REPO_ROOT / "reports/stage5_validation/s56_recovery50_teacher_queue"
DEFAULT_TRAINING_TABLE = DEFAULT_ROOT / "human_training_table/human_vetting_training_table.csv"
DEFAULT_OUT_DIR = DEFAULT_ROOT / "cnn_tensors"
DEFAULT_COMPACT_LC = (
    REPO_ROOT
    / "data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_adp_lc_export_pdo.h5"
)
DEFAULT_HLSP_ROOT = REPO_ROOT / "data_local/stage1_lightcurves/hlsp_s0056_twirl_fs_v2_compare"
DEFAULT_INJECTION_H5 = (
    REPO_ROOT
    / "data_local/stage3_injections/s56_twirlfs_v2_injection_training/"
    "recovery50_adp_pair_subset/injected_lightcurves.h5"
)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training-table", type=Path, default=DEFAULT_TRAINING_TABLE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--compact-lc-h5", type=Path, default=DEFAULT_COMPACT_LC)
    parser.add_argument("--hlsp-root", type=Path, default=DEFAULT_HLSP_ROOT)
    parser.add_argument("--injection-h5-override", type=Path, default=DEFAULT_INJECTION_H5)
    parser.add_argument("--apertures", default=",".join(DEFAULT_APERTURES))
    parser.add_argument("--folded-points", type=int, default=512)
    parser.add_argument("--context-points", type=int, default=512)
    parser.add_argument("--event-points", type=int, default=128)
    parser.add_argument("--max-events", type=int, default=16)
    parser.add_argument("--folded-window-durations", type=float, default=4.0)
    parser.add_argument("--context-window-durations", type=float, default=12.0)
    parser.add_argument("--event-window-durations", type=float, default=4.0)
    parser.add_argument("--min-event-points", type=int, default=1)
    parser.add_argument("--max-rows", type=int, default=None)
    parser.add_argument("--progress-every", type=int, default=50)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    apertures = tuple(part.strip() for part in args.apertures.split(",") if part.strip())
    cfg = TensorConfig(
        apertures=apertures,
        folded_points=args.folded_points,
        context_points=args.context_points,
        event_points=args.event_points,
        max_events=args.max_events,
        folded_window_durations=args.folded_window_durations,
        context_window_durations=args.context_window_durations,
        event_window_durations=args.event_window_durations,
        min_event_points=args.min_event_points,
    )
    summary = build_recovery50_cnn_tensors(
        training_table=args.training_table,
        out_dir=args.out_dir,
        compact_lc_h5=args.compact_lc_h5,
        hlsp_root=args.hlsp_root,
        injection_h5_override=args.injection_h5_override,
        config=cfg,
        max_rows=args.max_rows,
        progress_every=args.progress_every,
    )
    print(json.dumps(summary, indent=2, sort_keys=True, default=json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
