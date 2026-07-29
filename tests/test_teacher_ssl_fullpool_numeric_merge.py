from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from twirl.vetting.ssl_full_pool_eligibility import (
    observation_identity_sha256,
)


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "scripts"
    / "stage5_validation"
    / "merge_teacher_ssl_fullpool_numeric.py"
)
SPEC = importlib.util.spec_from_file_location(
    "merge_teacher_ssl_fullpool_numeric_test",
    SCRIPT,
)
assert SPEC is not None and SPEC.loader is not None
MERGE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MERGE)


def _passed_row(*, sector: int = 60, tic: int = 722_078_603) -> dict[str, object]:
    return {
        "ssl_observation_id": f"s{sector:04d}-tic{tic:019d}",
        "sector": sector,
        "tic": tic,
        "passed": True,
        "n_failures": 0,
        "failure_codes": [],
        "failures": [],
        "harmonic_max_abs": 1.0,
        "local_max_abs": 2.0,
        "periodogram_max_abs": 3.0,
    }


def test_merge_row_validator_requires_literal_true() -> None:
    row = _passed_row()
    row["passed"] = "False"
    identity = observation_identity_sha256([(60, 722_078_603)])

    with pytest.raises(ValueError, match="not an exact passed row"):
        MERGE._validate_passed_shard_rows(
            [row],
            sector=60,
            expected_count=1,
            expected_identity_sha256=identity,
            context="test-report",
        )


def test_merge_row_validator_accepts_exact_passed_evidence() -> None:
    row = _passed_row()
    identity = observation_identity_sha256([(60, 722_078_603)])

    frame = MERGE._validate_passed_shard_rows(
        [row],
        sector=60,
        expected_count=1,
        expected_identity_sha256=identity,
        context="test-report",
    )

    assert frame.to_dict(orient="records") == [row]


def test_numeric_envelope_survives_json_roundtrip() -> None:
    envelope = MERGE.FULL_POOL_NUMERIC_ENVELOPE_V1.as_dict()
    roundtripped = json.loads(json.dumps(envelope))

    assert roundtripped != envelope
    assert (
        MERGE._canonical_sha256(roundtripped)
        == MERGE.MODEL_INPUT_NUMERIC_ENVELOPE_SHA256
    )


def _small_audit(monkeypatch: pytest.MonkeyPatch) -> pd.DataFrame:
    eligible = {
        "ssl_observation_id": "s0060-tic0000000000722078603",
        "sector": 60,
        "tic": 722_078_603,
        "ssl_pool_include": True,
        "numeric_status": "passed",
        "model_input_numeric_passed": True,
        "n_failures": 0,
        "failure_codes": "[]",
        "failures_json": "[]",
        "harmonic_max_abs": 1.0,
        "local_max_abs": 2.0,
        "periodogram_max_abs": 3.0,
    }
    excluded = {
        "ssl_observation_id": "s0062-tic0000000000999999999",
        "sector": 62,
        "tic": 999_999_999,
        "ssl_pool_include": False,
        "numeric_status": "not_model_eligible",
        "model_input_numeric_passed": pd.NA,
        "n_failures": 0,
        "failure_codes": "[]",
        "failures_json": "[]",
        "harmonic_max_abs": np.nan,
        "local_max_abs": np.nan,
        "periodogram_max_abs": np.nan,
    }
    frame = pd.DataFrame.from_records(
        [eligible, excluded],
        columns=list(MERGE.AUDIT_COLUMNS),
    ).astype(MERGE.AUDIT_DTYPES)
    monkeypatch.setattr(MERGE, "PRODUCTION_FULL_OBSERVATIONS", 2)
    monkeypatch.setattr(MERGE, "PRODUCTION_ELIGIBLE_OBSERVATIONS", 1)
    monkeypatch.setattr(MERGE, "PRODUCTION_EXCLUDED_OBSERVATIONS", 1)
    monkeypatch.setattr(
        MERGE,
        "PRODUCTION_FULL_IDENTITY_SHA256",
        observation_identity_sha256([(60, 722_078_603), (62, 999_999_999)]),
    )
    monkeypatch.setattr(
        MERGE,
        "PRODUCTION_ELIGIBLE_IDENTITY_SHA256",
        observation_identity_sha256([(60, 722_078_603)]),
    )
    monkeypatch.setattr(
        MERGE,
        "PRODUCTION_EXCLUDED_IDENTITY_SHA256",
        observation_identity_sha256([(62, 999_999_999)]),
    )
    return frame


def test_merge_audit_validator_locks_roundtrip_dtypes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    frame = _small_audit(monkeypatch)
    MERGE._validate_numeric_audit_frame(frame)

    drifted = frame.copy()
    drifted["sector"] = drifted["sector"].astype("int16")
    with pytest.raises(ValueError, match="wrong exact dtypes"):
        MERGE._validate_numeric_audit_frame(drifted)
