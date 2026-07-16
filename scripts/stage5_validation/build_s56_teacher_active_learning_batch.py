#!/usr/bin/env python3
"""Build one blinded S56 A2v1 teacher active-learning batch."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.vetting.teacher_active_learning import write_active_learning_batch  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--batch-index", type=int, required=True)
    parser.add_argument("--exclude", type=Path, nargs="*", default=[])
    args = parser.parse_args()
    summary = write_active_learning_batch(
        scores_path=args.scores,
        out_dir=args.out_dir,
        batch_index=args.batch_index,
        exclude_paths=args.exclude,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
