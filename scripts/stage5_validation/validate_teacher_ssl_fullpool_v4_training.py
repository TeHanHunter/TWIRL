#!/usr/bin/env python3
"""Validate and freeze the detector-consistent v4 five-fold SSL run."""
from __future__ import annotations

import importlib.util
from pathlib import Path


BASE_PATH = Path(__file__).with_name(
    "validate_teacher_ssl_fullpool_v3_training.py"
)
SPEC = importlib.util.spec_from_file_location(
    "validate_teacher_ssl_fullpool_v3_training_preserved",
    BASE_PATH,
)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load the preserved v3 training validator")
BASE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)

BASE.COMPLETION_RELEASE_SCHEMA = (
    "twirl_teacher_ssl_fullpool_v4_five_fold_completion_release_v1"
)
BASE.COMPLETION_RELEASE_BINDING = (
    "teacher_ssl_fullpool_v4_detector_consistent_raw_v1_"
    "effective_quality_adp_btjd_five_fold_complete_v1"
)


if __name__ == "__main__":
    raise SystemExit(BASE.main())
