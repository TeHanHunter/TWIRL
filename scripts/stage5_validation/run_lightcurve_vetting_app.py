#!/usr/bin/env python3
"""Run a local TWIRL light-curve vetting and labeling app."""
from __future__ import annotations

import argparse
import getpass
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

DEFAULT_CANDIDATES = (
    REPO_ROOT
    / "reports/stage5_validation/s56_pretriage_review_queue_pdo/review_queue.csv"
)
DEFAULT_LABELS_OUT = (
    REPO_ROOT
    / "reports/stage5_validation/s56_pretriage_review_queue_pdo/human_labels_vetted.csv"
)
DEFAULT_HLSP_ROOT = (
    REPO_ROOT / "data_local/stage1_lightcurves/hlsp_s0056_twirl_fs_v2_compare"
)
DEFAULT_LEO_REPORT_ROOTS = (
    REPO_ROOT / "reports/stage5_validation/s56_pretriage_review_queue_pdo/vet_reports",
    REPO_ROOT / "reports/stage5_validation/s56_pretriage_review_queue/vet_reports",
    REPO_ROOT / "reports/stage5_validation/leo_vetter_s56_v2_full/vet_reports",
    REPO_ROOT / "reports/stage5_validation/leo_vetter_s56_top50_v2/vet_reports",
    REPO_ROOT / "reports/stage5_validation/leo_vetter_s56_top50/vet_reports",
)


def _build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--candidates", type=Path, default=DEFAULT_CANDIDATES)
    ap.add_argument("--labels-out", type=Path, default=DEFAULT_LABELS_OUT)
    ap.add_argument("--hlsp-root", type=Path, default=DEFAULT_HLSP_ROOT)
    ap.add_argument(
        "--leo-report-root",
        type=Path,
        action="append",
        default=None,
        help="Directory containing pre-rendered LEO-Vetter PDF reports. May be repeated.",
    )
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=5000)
    ap.add_argument("--aperture", default="DET_FLUX")
    ap.add_argument("--labeler", default=getpass.getuser())
    ap.add_argument("--shuffle-order", action="store_true",
                    help="Use a deterministic shuffled review order without rewriting the candidate CSV.")
    ap.add_argument("--shuffle-seed", type=int, default=56)
    ap.add_argument("--no-unlabeled-first", action="store_true",
                    help="Do not prioritize unlabeled rows ahead of already-labeled rows when shuffling.")
    ap.add_argument("--smoke-plot", type=Path, default=None,
                    help="Write one candidate plot PNG and exit.")
    ap.add_argument("--smoke-index", type=int, default=0)
    return ap


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    if not args.candidates.exists():
        print(f"[vet-app] missing candidate CSV: {args.candidates}", file=sys.stderr)
        return 2
    if not args.hlsp_root.exists():
        print(f"[vet-app] missing HLSP root: {args.hlsp_root}", file=sys.stderr)
        return 2

    from twirl.vetting.lightcurve_label_app import LightCurveVettingApp

    app = LightCurveVettingApp(
        candidates_path=args.candidates,
        labels_out=args.labels_out,
        hlsp_root=args.hlsp_root,
        leo_report_roots=tuple(args.leo_report_root or DEFAULT_LEO_REPORT_ROOTS),
        default_aperture=args.aperture,
        labeler=args.labeler,
        shuffle_order=args.shuffle_order,
        random_seed=args.shuffle_seed,
        unlabeled_first=not args.no_unlabeled_first,
    )

    if args.smoke_plot is not None:
        args.smoke_plot.parent.mkdir(parents=True, exist_ok=True)
        args.smoke_plot.write_bytes(app.plot_png(args.smoke_index, args.aperture))
        payload = app.candidate_payload(args.smoke_index)
        print(f"[vet-app] wrote smoke plot: {args.smoke_plot}")
        print(f"[vet-app] TIC: {payload['display'].get('tic')}")
        print(f"[vet-app] HLSP: {payload['hlsp_path']}")
        return 0

    print("[vet-app] local tunnel example:")
    print(f"  ssh -L {args.port}:localhost:{args.port} pdogpu6")
    print(f"  open http://localhost:{args.port}/")
    app.serve(host=args.host, port=args.port)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
