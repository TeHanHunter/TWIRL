#!/usr/bin/env python3
"""Score S56 CNN tensors with one saved ADP-only teacher checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.vetting.recovery50_cnn import score_recovery50_cnn_teacher  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tensor-npz", type=Path, required=True)
    parser.add_argument("--tensor-rows", type=Path, required=True)
    parser.add_argument("--model-path", type=Path, required=True)
    parser.add_argument("--feature-table", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--require-cuda", action="store_true")
    args = parser.parse_args(argv)
    summary = score_recovery50_cnn_teacher(
        tensor_npz=args.tensor_npz,
        tensor_rows=args.tensor_rows,
        model_path=args.model_path,
        feature_table=args.feature_table,
        out_dir=args.out_dir,
        batch_size=args.batch_size,
        require_cuda=args.require_cuda,
    )
    print(json.dumps(summary, indent=2, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
