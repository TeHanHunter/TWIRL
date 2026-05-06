#!/usr/bin/env python3
"""Thin CLI driver for the per-sector BLS first pass.

Wraps `twirl.search.sector_run.main` so it is callable as a script in addition
to `python -m twirl.search.sector_run`. All flags are forwarded.

Example
-------
    OMP_NUM_THREADS=1 \\
        python scripts/stage2_search/run_sector_bls.py --sector 56 --workers 32
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.search.sector_run import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
