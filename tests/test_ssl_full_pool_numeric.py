from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from twirl.vetting import ssl_full_pool_numeric as numeric_module
from twirl.vetting.ssl_full_pool_eligibility import (
    observation_identity_sha256,
)
from twirl.vetting.ssl_full_pool_numeric import (
    FLOAT32_SQUARE_SAFE_ABS_MAX,
    FULL_POOL_NUMERIC_ENVELOPE_V1,
    MODEL_INPUT_CONTRACT_VERSION,
    MODEL_INPUT_NUMERIC_AUDIT_SCHEMA,
    MODEL_INPUT_NUMERIC_AUTHORITY_NAMES,
    MODEL_INPUT_NUMERIC_ENVELOPE_SHA256,
    MODEL_INPUT_NUMERIC_RELEASE_SCHEMA,
    TEACHER_SSL_NUMERIC_ENVELOPE_V1,
    audit_collated_sample,
    audit_model_facing_sample,
    validate_numeric_gate_release,
)


def _valid_view(n_cadences: int = 4) -> tuple[np.ndarray, np.ndarray]:
    values = np.zeros((7, n_cadences), dtype=np.float32)
    values[0] = np.linspace(-0.1, 0.1, n_cadences, dtype=np.float32)
    values[1] = np.linspace(-0.2, 0.2, n_cadences, dtype=np.float32)
    values[2] = values[1] - values[0]
    values[3] = np.float32(0.02)
    values[4] = np.float32(0.03)
    values[5] = np.linspace(-0.5, 0.49, n_cadences, dtype=np.float32)
    values[6] = np.asarray([0, 1, 0, 1], dtype=np.float32)[:n_cadences]
    mask = np.ones_like(values, dtype=bool)
    mask[:5, values[6].astype(bool)] = False
    return values, mask


def _sample() -> dict[str, object]:
    harmonic_pairs = [_valid_view() for _ in range(7)]
    local_pairs = [_valid_view(3) for _ in range(14)]
    return {
        "review_id": "s0060-tic0000000722078603",
        "tic": 722_078_603,
        "period_d": np.float32(0.47312935),
        "duration_min": np.float32(3.0),
        "chronology_small": np.zeros((10, 1), dtype=np.float32),
        "chronology_small_mask": np.zeros((10, 1), dtype=bool),
        "chronology_supplemental": np.zeros((5, 1), dtype=np.float32),
        "chronology_supplemental_mask": np.zeros((5, 1), dtype=bool),
        "harmonic_values": tuple(pair[0] for pair in harmonic_pairs),
        "harmonic_mask": tuple(pair[1] for pair in harmonic_pairs),
        "local_values": tuple(pair[0] for pair in local_pairs),
        "local_mask": tuple(pair[1] for pair in local_pairs),
        "periodogram_values": np.asarray(
            [
                [-1.0, 0.0, 1.0],
                [-2.0, 0.0, 2.0],
                [-0.5, 0.0, 0.5],
                [-3.0, 0.0, 3.0],
            ],
            dtype=np.float32,
        ),
        "periodogram_mask": np.ones((4, 3), dtype=bool),
        "metadata": np.asarray([-1.0, 0.0, 2.0], dtype=np.float32),
        "morphology_target": -1,
        "preserve_target": -1,
        "harmonic_target": -1,
        "morphology_weight": np.float32(0.0),
        "preserve_weight": np.float32(0.0),
        "harmonic_weight": np.float32(0.0),
        "pretrain_target": -1,
        "compact_target": -1,
        "compact_weight": np.float32(0.0),
    }


def _failure_codes(report: dict[str, object]) -> set[str]:
    return {
        str(item["code"])
        for item in report["failures"]  # type: ignore[union-attr]
    }


def _numpy_batch(sample: dict[str, object]) -> dict[str, object]:
    return {
        "review_id": [sample["review_id"]],
        "tic": np.asarray([sample["tic"]], dtype=np.int64),
        "period_d": np.asarray([sample["period_d"]], dtype=np.float32),
        "duration_min": np.asarray([sample["duration_min"]], dtype=np.float32),
        "harmonic_values": np.stack(
            [np.stack(sample["harmonic_values"])], axis=0  # type: ignore[arg-type]
        ),
        "harmonic_mask": np.stack(
            [np.stack(sample["harmonic_mask"])], axis=0  # type: ignore[arg-type]
        ),
        "local_values": np.stack(
            [np.stack(sample["local_values"])], axis=0  # type: ignore[arg-type]
        ),
        "local_mask": np.stack(
            [np.stack(sample["local_mask"])], axis=0  # type: ignore[arg-type]
        ),
        "periodogram_values": np.stack(
            [sample["periodogram_values"]], axis=0
        ),
        "periodogram_mask": np.stack(
            [sample["periodogram_mask"]], axis=0
        ),
        "metadata": np.stack([sample["metadata"]], axis=0),
    }


def _audit(sample: dict[str, object]) -> dict[str, object]:
    return audit_collated_sample(_numpy_batch(sample))


def test_masked_extreme_is_not_a_model_facing_failure() -> None:
    sample = _sample()
    harmonic_values = list(sample["harmonic_values"])  # type: ignore[arg-type]
    harmonic_masks = list(sample["harmonic_mask"])  # type: ignore[arg-type]
    harmonic_values[0][0, 0] = np.float32(5.0e36)
    harmonic_masks[0][0, 0] = False
    sample["harmonic_values"] = tuple(harmonic_values)
    sample["harmonic_mask"] = tuple(harmonic_masks)

    report = _audit(sample)

    assert report["contract_version"] == TEACHER_SSL_NUMERIC_ENVELOPE_V1
    assert report["passed"] is True
    assert report["failures"] == []
    small = report["maxima"]["harmonic_values"]["adp_small_relative"]
    assert small["max_abs"] == pytest.approx(0.1)
    assert report["action"] == "audit_only_no_clip_no_exclusion"


def test_masked_nonfinite_payload_is_rejected_before_model_multiplication() -> None:
    sample = _sample()
    harmonic_values = list(sample["harmonic_values"])  # type: ignore[arg-type]
    harmonic_masks = list(sample["harmonic_mask"])  # type: ignore[arg-type]
    harmonic_values[0][0, 0] = np.float32(np.inf)
    harmonic_masks[0][0, 0] = False
    sample["harmonic_values"] = tuple(harmonic_values)
    sample["harmonic_mask"] = tuple(harmonic_masks)

    report = _audit(sample)

    assert report["passed"] is False
    assert "tensor_payload_nonfinite" in _failure_codes(report)


def test_masked_in_extreme_fails_without_clipping_or_exclusion() -> None:
    sample = _sample()
    harmonic_values = list(sample["harmonic_values"])  # type: ignore[arg-type]
    harmonic_values[0][0, 0] = np.float32(5.0e36)
    sample["harmonic_values"] = tuple(harmonic_values)

    report = _audit(sample)

    assert report["passed"] is False
    assert "numeric_envelope_exceeded" in _failure_codes(report)
    assert "float32_square_overflow" in _failure_codes(report)
    assert report["maxima"]["harmonic_values"]["adp_small_relative"][
        "max_abs"
    ] == pytest.approx(5.0e36)
    assert harmonic_values[0][0, 0] == np.float32(5.0e36)


def test_float32_square_boundary_is_checked_on_collated_tensor() -> None:
    sample = _sample()
    periodogram = sample["periodogram_values"]
    assert isinstance(periodogram, np.ndarray)
    periodogram[0, 0] = np.nextafter(
        np.float32(FLOAT32_SQUARE_SAFE_ABS_MAX),
        np.float32(np.inf),
    )
    report = _audit(sample)

    failures = [
        item
        for item in report["failures"]
        if item["code"] == "float32_square_overflow"
    ]
    assert failures
    assert failures[0]["tensor"] == "periodogram_values"


def test_v1_million_scale_envelope_is_fail_closed_before_float32_overflow() -> None:
    sample = _sample()
    periodogram = sample["periodogram_values"]
    assert isinstance(periodogram, np.ndarray)
    periodogram[0, 0] = np.float32(1.0e6 + 1.0)

    report = _audit(sample)

    assert "numeric_envelope_exceeded" in _failure_codes(report)
    assert "float32_square_overflow" not in _failure_codes(report)


def test_coordinate_binary_and_error_domains_are_fail_closed() -> None:
    sample = _sample()
    harmonic_values = list(sample["harmonic_values"])  # type: ignore[arg-type]
    harmonic_values[0][3, 0] = np.float32(-0.01)
    harmonic_values[0][5, 1] = np.float32(0.51)
    harmonic_values[0][6, 2] = np.float32(0.5)
    sample["harmonic_values"] = tuple(harmonic_values)

    report = _audit(sample)

    codes = _failure_codes(report)
    assert report["passed"] is False
    assert {"negative_error", "coordinate_out_of_range", "binary_domain_violation"} <= codes
    assert {
        item["channel_name"]
        for item in report["failures"]
        if "channel_name" in item
    } >= {
        "raw_error_small_scaled",
        "orbital_phase",
        "quality_nonzero",
    }


def test_masked_domain_violations_are_ignored_but_counted_as_masked() -> None:
    sample = _sample()
    local_values = list(sample["local_values"])  # type: ignore[arg-type]
    local_masks = list(sample["local_mask"])  # type: ignore[arg-type]
    local_values[0][3, 0] = np.float32(-5.0)
    local_values[0][5, 1] = np.float32(2.0)
    local_values[0][6, 2] = np.float32(0.25)
    local_masks[0][3, 0] = False
    local_masks[0][:, 1] = False
    local_masks[0][:, 2] = False
    sample["local_values"] = tuple(local_values)
    sample["local_mask"] = tuple(local_masks)

    report = _audit(sample)

    assert report["passed"] is True
    assert report["maxima"]["local_values"]["raw_error_small_scaled"][
        "n_masked_values"
    ] >= 1


def test_quality_bad_cadence_requires_all_photometry_and_errors_masked() -> None:
    sample = _sample()
    harmonic_masks = list(sample["harmonic_mask"])  # type: ignore[arg-type]
    # Cadence one has quality_nonzero=1 and retains active phase/quality
    # coordinates. Re-enabling one photometric channel violates the v1
    # effective-quality model-input contract.
    harmonic_masks[0][0, 1] = True
    sample["harmonic_mask"] = tuple(harmonic_masks)

    report = _audit(sample)

    failures = [
        item
        for item in report["failures"]
        if item["code"] == "quality_bad_photometry_unmasked"
    ]
    assert failures
    assert failures[0]["tensor"] == "harmonic_values"
    assert failures[0]["count"] == 1


def test_harmonic_coordinates_cannot_hide_quality_behind_a_false_mask() -> None:
    sample = _sample()
    harmonic_masks = list(sample["harmonic_mask"])  # type: ignore[arg-type]
    harmonic_masks[0][6, 1] = False
    harmonic_masks[0][0, 1] = True
    sample["harmonic_mask"] = tuple(harmonic_masks)

    report = _audit(sample)

    assert {
        "phase_quality_mask_mismatch",
        "photometry_without_coordinates",
    } <= _failure_codes(report)


def test_dataset_sample_wrapper_uses_production_collator() -> None:
    pytest.importorskip("torch")

    report = audit_model_facing_sample(_sample())

    assert report["passed"] is True


def _bound_artifact(path: Path) -> dict[str, object]:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    Path(str(path) + ".sha256").write_text(
        f"{digest}  {path.name}\n",
        encoding="ascii",
    )
    return {
        "path": str(path),
        "sha256": digest,
        "size_bytes": path.stat().st_size,
    }


def _write_json_artifact(
    path: Path,
    payload: object,
) -> dict[str, object]:
    path.write_text(
        json.dumps(
            payload,
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
        + "\n",
        encoding="utf-8",
    )
    return _bound_artifact(path)


def _write_numeric_release(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> tuple[Path, dict[str, object]]:
    bindings: dict[str, dict[str, object]] = {}
    for index, name in enumerate(MODEL_INPUT_NUMERIC_AUTHORITY_NAMES):
        authority = tmp_path / f"{name}.authority"
        authority.write_bytes(f"authority-{index}\n".encode("ascii"))
        bindings[name] = _bound_artifact(authority)

    eligible_rows: list[dict[str, object]] = []
    for task_id in range(112):
        sector = 56 + task_id // 16
        tic = 10_000_000 + task_id
        eligible_rows.append(
            {
                "ssl_observation_id": f"s{sector:04d}-tic{tic}",
                "sector": sector,
                "tic": tic,
                "passed": True,
                "n_failures": 0,
                "failure_codes": [],
                "failures": [],
                "harmonic_max_abs": float(task_id + 1),
                "local_max_abs": float(task_id + 2),
                "periodogram_max_abs": float(task_id + 3),
            }
        )
    excluded_rows = [
        {
            "ssl_observation_id": "s0062-tic99999999",
            "sector": 62,
            "tic": 99_999_999,
        }
    ]
    eligible_identity = observation_identity_sha256(
        [(int(row["sector"]), int(row["tic"])) for row in eligible_rows]
    )
    excluded_identity = observation_identity_sha256(
        [(int(row["sector"]), int(row["tic"])) for row in excluded_rows]
    )
    full_identity = observation_identity_sha256(
        [
            (int(row["sector"]), int(row["tic"]))
            for row in [*eligible_rows, *excluded_rows]
        ]
    )
    monkeypatch.setattr(
        numeric_module,
        "PRODUCTION_FULL_OBSERVATIONS",
        len(eligible_rows) + len(excluded_rows),
    )
    monkeypatch.setattr(
        numeric_module,
        "PRODUCTION_ELIGIBLE_OBSERVATIONS",
        len(eligible_rows),
    )
    monkeypatch.setattr(
        numeric_module,
        "PRODUCTION_EXCLUDED_OBSERVATIONS",
        len(excluded_rows),
    )
    monkeypatch.setattr(
        numeric_module,
        "PRODUCTION_FULL_IDENTITY_SHA256",
        full_identity,
    )
    monkeypatch.setattr(
        numeric_module,
        "PRODUCTION_ELIGIBLE_IDENTITY_SHA256",
        eligible_identity,
    )
    monkeypatch.setattr(
        numeric_module,
        "PRODUCTION_EXCLUDED_IDENTITY_SHA256",
        excluded_identity,
    )

    audit_rows = [
        {
            "ssl_observation_id": row["ssl_observation_id"],
            "sector": row["sector"],
            "tic": row["tic"],
            "ssl_pool_include": True,
            "numeric_status": "passed",
            "model_input_numeric_passed": True,
            "n_failures": 0,
            "failure_codes": "[]",
            "failures_json": "[]",
            "harmonic_max_abs": row["harmonic_max_abs"],
            "local_max_abs": row["local_max_abs"],
            "periodogram_max_abs": row["periodogram_max_abs"],
        }
        for row in eligible_rows
    ]
    audit_rows.extend(
        {
            "ssl_observation_id": row["ssl_observation_id"],
            "sector": row["sector"],
            "tic": row["tic"],
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
        for row in excluded_rows
    )
    audit_frame = pd.DataFrame.from_records(
        audit_rows,
        columns=list(numeric_module._MODEL_INPUT_NUMERIC_AUDIT_COLUMNS),
    )
    audit_frame["sector"] = audit_frame["sector"].astype("int64")
    audit_frame["tic"] = audit_frame["tic"].astype("int64")
    audit_frame["ssl_pool_include"] = audit_frame[
        "ssl_pool_include"
    ].astype("bool")
    audit_frame["model_input_numeric_passed"] = audit_frame[
        "model_input_numeric_passed"
    ].astype("boolean")
    audit_frame["n_failures"] = audit_frame["n_failures"].astype("int64")
    for name in (
        "harmonic_max_abs",
        "local_max_abs",
        "periodogram_max_abs",
    ):
        audit_frame[name] = audit_frame[name].astype("float64")
    audit = tmp_path / "model_input_numeric_audit.parquet"
    audit.write_bytes(numeric_module.numeric_audit_parquet_bytes(audit_frame))
    audit_binding = _bound_artifact(audit)

    shard_reports: list[dict[str, object]] = []
    for task_id, row in enumerate(eligible_rows):
        sector = 56 + task_id // 16
        shard_index = task_id % 16
        native = tmp_path / f"native_{task_id}.h5"
        native.write_bytes(f"native-{task_id}\n".encode("ascii"))
        native_binding = _bound_artifact(native)
        identity = observation_identity_sha256(
            [(int(row["sector"]), int(row["tic"]))]
        )
        report_payload = {
            "schema_version": MODEL_INPUT_NUMERIC_AUDIT_SCHEMA,
            "code_revision": "1" * 40,
            "model_input_contract_version": MODEL_INPUT_CONTRACT_VERSION,
            "envelope_contract": TEACHER_SSL_NUMERIC_ENVELOPE_V1,
            "envelope_canonical_sha256": (
                MODEL_INPUT_NUMERIC_ENVELOPE_SHA256
            ),
            "envelope": FULL_POOL_NUMERIC_ENVELOPE_V1.as_dict(),
            "sector": sector,
            "shard_index": shard_index,
            "n_shards": 16,
            "passed": True,
            "counts": {
                "scanned_observations": 1,
                "passed_observations": 1,
                "failed_observations": 0,
            },
            "observation_identity_sha256": identity,
            "native_h5": native_binding,
            "authority_bindings": bindings,
            "rows": [row],
            "action": "audit_only_no_clip_no_exclusion",
        }
        report = tmp_path / f"numeric_{task_id}.json"
        report_binding = _write_json_artifact(report, report_payload)
        shard_reports.append(
            {
                "sector": sector,
                "shard_index": shard_index,
                **report_binding,
                "native_h5": native_binding,
                "observation_identity_sha256": identity,
                "scanned_observations": 1,
            }
        )
    payload: dict[str, object] = {
        "schema_version": MODEL_INPUT_NUMERIC_RELEASE_SCHEMA,
        "passed": True,
        "model_input_contract_version": MODEL_INPUT_CONTRACT_VERSION,
        "envelope_contract": TEACHER_SSL_NUMERIC_ENVELOPE_V1,
        "envelope": FULL_POOL_NUMERIC_ENVELOPE_V1.as_dict(),
        "envelope_canonical_sha256": MODEL_INPUT_NUMERIC_ENVELOPE_SHA256,
        "counts": {
            "full_observations": len(eligible_rows) + len(excluded_rows),
            "eligible_observations": len(eligible_rows),
            "excluded_observations": len(excluded_rows),
            "scanned_observations": len(eligible_rows),
            "failed_observations": 0,
            "native_shards": 112,
        },
        "identity_hashes": {
            "full": full_identity,
            "eligible": eligible_identity,
            "excluded": excluded_identity,
        },
        "code_revision": "1" * 40,
        "authority_bindings": bindings,
        "outputs": {"numeric_audit": audit_binding},
        "shard_reports": shard_reports,
        "quality_bad_photometry_policy_verified": True,
        "float32_conversion_verified": True,
        "real_only": True,
        "labels_consumed": False,
        "injections_consumed": False,
        "action": "audit_only_no_clip_no_exclusion",
    }
    path = tmp_path / "numeric_release.json"
    _write_json_artifact(path, payload)
    return path, payload


def _two_row_numeric_audit() -> pd.DataFrame:
    frame = pd.DataFrame.from_records(
        [
            {
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
            },
            {
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
            },
        ],
        columns=list(numeric_module._MODEL_INPUT_NUMERIC_AUDIT_COLUMNS),
    )
    return frame.astype(
        dict(
            zip(
                numeric_module._MODEL_INPUT_NUMERIC_AUDIT_COLUMNS,
                numeric_module._MODEL_INPUT_NUMERIC_AUDIT_DTYPES,
                strict=True,
            )
        )
    )


def test_numeric_audit_arrow_schema_is_exact_and_nullable_pass_survives() -> None:
    frame = _two_row_numeric_audit()
    for name in numeric_module._MODEL_INPUT_NUMERIC_AUDIT_TEXT_COLUMNS:
        frame[name] = frame[name].astype("string")

    payload = numeric_module.numeric_audit_parquet_bytes(frame)
    parquet = pq.ParquetFile(pa.BufferReader(payload))

    assert parquet.schema_arrow.equals(
        numeric_module.MODEL_INPUT_NUMERIC_AUDIT_ARROW_SCHEMA,
        check_metadata=False,
    )
    observed = numeric_module.read_numeric_audit_parquet(
        pa.BufferReader(payload)
    )
    assert tuple(
        str(observed[name].dtype)
        for name in numeric_module._MODEL_INPUT_NUMERIC_AUDIT_COLUMNS
    ) == numeric_module._MODEL_INPUT_NUMERIC_AUDIT_DTYPES
    assert observed["model_input_numeric_passed"].iloc[0] == np.bool_(True)
    assert pd.isna(observed["model_input_numeric_passed"].iloc[1])


def test_numeric_audit_reader_rejects_arrow_physical_type_drift() -> None:
    payload = numeric_module.numeric_audit_parquet_bytes(
        _two_row_numeric_audit()
    )
    table = pq.ParquetFile(pa.BufferReader(payload)).read()
    sector_index = table.schema.get_field_index("sector")
    drifted = table.set_column(
        sector_index,
        pa.field("sector", pa.int32(), nullable=False),
        table.column(sector_index).cast(pa.int32()),
    )
    sink = pa.BufferOutputStream()
    pq.write_table(drifted, sink)

    with pytest.raises(ValueError, match="wrong exact Arrow schema"):
        numeric_module.read_numeric_audit_parquet(
            pa.BufferReader(sink.getvalue())
        )


def test_numeric_audit_normalization_rejects_string_pass_value() -> None:
    frame = _two_row_numeric_audit()
    frame["model_input_numeric_passed"] = frame[
        "model_input_numeric_passed"
    ].astype("object")
    frame.loc[0, "model_input_numeric_passed"] = "True"

    with pytest.raises(ValueError, match="must contain booleans or nulls"):
        numeric_module.normalize_numeric_audit_frame(frame)


def test_numeric_release_validator_rehashes_exact_authorities(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    path, payload = _write_numeric_release(tmp_path, monkeypatch)

    observed = validate_numeric_gate_release(
        path,
        expected_code_revision="1" * 40,
        expected_authority_bindings=payload["authority_bindings"],  # type: ignore[arg-type]
    )

    assert observed["passed"] is True
    assert observed["counts"]["scanned_observations"] == 112


def test_numeric_release_validator_rejects_failure_or_changed_authority(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    path, payload = _write_numeric_release(tmp_path, monkeypatch)
    payload["counts"]["failed_observations"] = 1  # type: ignore[index]
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="counts differ from production"):
        validate_numeric_gate_release(path)

    path, payload = _write_numeric_release(tmp_path, monkeypatch)
    binding = payload["authority_bindings"]["native_registry"]  # type: ignore[index]
    Path(binding["path"]).write_text("changed\n", encoding="utf-8")
    with pytest.raises(ValueError, match="file size differs|file SHA-256 differs"):
        validate_numeric_gate_release(path)


def test_numeric_release_validator_rejects_arbitrary_audit_bytes(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    path, payload = _write_numeric_release(tmp_path, monkeypatch)
    audit_binding = payload["outputs"]["numeric_audit"]  # type: ignore[index]
    audit_path = Path(audit_binding["path"])
    audit_path.write_bytes(b"not a Parquet numeric audit\n")
    audit_binding.update(_bound_artifact(audit_path))
    _write_json_artifact(path, payload)

    with pytest.raises(ValueError, match="readable Parquet"):
        validate_numeric_gate_release(path)


def test_numeric_release_validator_rejects_header_only_report(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    path, payload = _write_numeric_release(tmp_path, monkeypatch)
    report_binding = payload["shard_reports"][0]  # type: ignore[index]
    report_path = Path(report_binding["path"])
    report = json.loads(report_path.read_text(encoding="utf-8"))
    del report["rows"]
    report_binding.update(_write_json_artifact(report_path, report))
    _write_json_artifact(path, payload)

    with pytest.raises(ValueError, match="rows differ from the declared count"):
        validate_numeric_gate_release(path)


def test_numeric_release_validator_rejects_string_boolean_report_row(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    path, payload = _write_numeric_release(tmp_path, monkeypatch)
    report_binding = payload["shard_reports"][0]  # type: ignore[index]
    report_path = Path(report_binding["path"])
    report = json.loads(report_path.read_text(encoding="utf-8"))
    report["rows"][0]["passed"] = "False"
    report_binding.update(_write_json_artifact(report_path, report))
    _write_json_artifact(path, payload)

    with pytest.raises(ValueError, match="invalid strict field types"):
        validate_numeric_gate_release(path)
