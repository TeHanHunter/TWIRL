#!/usr/bin/env python3
"""Apply the selected five-fold harmonic teacher to real A2v1 candidates."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.vetting.harmonic_inference import (  # noqa: E402
    SELECTED_TEACHER_PROFILE,
    score_harmonic_teacher_to_disk,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--native-h5", type=Path, required=True)
    parser.add_argument("--checkpoints", type=Path, nargs=5, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--profile", default=SELECTED_TEACHER_PROFILE)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--allow-cpu", action="store_true")
    args = parser.parse_args()
    summary = score_harmonic_teacher_to_disk(
        candidates_path=args.candidates,
        native_h5=args.native_h5,
        checkpoint_paths=args.checkpoints,
        out_dir=args.out_dir,
        profile=args.profile,
        batch_size=args.batch_size,
        workers=args.workers,
        require_cuda=not args.allow_cpu,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
