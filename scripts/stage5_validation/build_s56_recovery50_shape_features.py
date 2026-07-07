#!/usr/bin/env python3
"""Build folded ADP small+primary shape features for the S56 recovery50 queue."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.vetting.recovery50_teacher import (  # noqa: E402
    DEFAULT_APERTURES,
    DEFAULT_SHAPE_BINS,
    DEFAULT_SHAPE_WINDOW_DURATIONS,
    build_folded_shape_features,
    json_default,
)


DEFAULT_QUEUE = REPO_ROOT / "reports/stage5_validation/s56_recovery50_teacher_queue/review_queue_1k.csv"
DEFAULT_OUT_DIR = DEFAULT_QUEUE.with_name("folded_shape_features")
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


def _optional_path(value: str | Path | None) -> Path | None:
    if value is None:
        return None
    text = str(value)
    if not text:
        return None
    path = Path(text)
    return path if path.exists() else None


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queue-csv", type=Path, default=DEFAULT_QUEUE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--compact-lc-h5", type=Path, default=DEFAULT_COMPACT_LC)
    parser.add_argument("--hlsp-root", type=Path, default=DEFAULT_HLSP_ROOT)
    parser.add_argument("--injection-h5-override", type=Path, default=DEFAULT_INJECTION_H5)
    parser.add_argument("--apertures", default=",".join(DEFAULT_APERTURES))
    parser.add_argument("--n-bins", type=int, default=DEFAULT_SHAPE_BINS)
    parser.add_argument("--window-durations", type=float, default=DEFAULT_SHAPE_WINDOW_DURATIONS)
    parser.add_argument("--max-rows", type=int, default=None)
    parser.add_argument("--progress-every", type=int, default=100)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    apertures = tuple(part.strip() for part in args.apertures.split(",") if part.strip())
    summary = build_folded_shape_features(
        queue_csv=args.queue_csv,
        out_dir=args.out_dir,
        compact_lc_h5=_optional_path(args.compact_lc_h5),
        hlsp_root=_optional_path(args.hlsp_root),
        injection_h5_override=_optional_path(args.injection_h5_override),
        apertures=apertures,
        n_bins=args.n_bins,
        window_durations=args.window_durations,
        max_rows=args.max_rows,
        progress_every=args.progress_every,
    )
    print(json.dumps(summary, indent=2, sort_keys=True, default=json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
