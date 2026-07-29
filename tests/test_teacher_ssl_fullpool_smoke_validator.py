from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from typing import Any, Callable

import h5py
import numpy as np
import pandas as pd
import pytest

from twirl.vetting import ssl_full_pool_numeric as NUMERIC
from twirl.vetting.ssl_full_pool_eligibility import (
    PRODUCTION_ELIGIBLE_IDENTITY_SHA256,
    PRODUCTION_ELIGIBLE_OBSERVATIONS,
    PRODUCTION_EXCLUDED_IDENTITY_SHA256,
    PRODUCTION_EXCLUDED_OBSERVATIONS,
    PRODUCTION_FULL_IDENTITY_SHA256,
    PRODUCTION_FULL_OBSERVATIONS,
    observation_identity_sha256,
)
from twirl.vetting.ssl_full_pool_native import (
    FULL_POOL_NATIVE_CONTRACT_VERSION,
)
from twirl.vetting.ssl_full_pool_numeric import (
    FULL_POOL_NUMERIC_ENVELOPE_V1,
    MODEL_INPUT_CONTRACT_VERSION,
    MODEL_INPUT_NUMERIC_AUDIT_SCHEMA,
    MODEL_INPUT_NUMERIC_ENVELOPE_CONTRACT,
    MODEL_INPUT_NUMERIC_ENVELOPE_SHA256,
    MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS,
    MODEL_INPUT_NUMERIC_RELEASE_BINDING,
    MODEL_INPUT_NUMERIC_RELEASE_SCHEMA,
    numeric_native_freshness,
)
from twirl.vetting.teacher_ssl_fullpool import (
    FULLPOOL_SSL_CHECKPOINT_NAMESPACE,
    FULLPOOL_SSL_MODEL_NAMESPACE,
    FULLPOOL_SSL_RUN_CONTRACT_SCHEMA,
    FULLPOOL_SSL_RUN_ID,
    FULLPOOL_SSL_SELECTION_SCHEMA,
    FULLPOOL_SSL_SUMMARY_SCHEMA,
    FULLPOOL_SSL_TRAINING_AUTHORITY_SCHEMA,
)

torch = pytest.importorskip("torch")


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "scripts"
    / "stage5_validation"
    / "validate_teacher_ssl_fullpool_v3_smoke.py"
)
SPEC = importlib.util.spec_from_file_location(
    "validate_teacher_ssl_fullpool_v3_smoke_test",
    SCRIPT,
)
assert SPEC is not None and SPEC.loader is not None
VALIDATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VALIDATOR)


def _numeric_test_keys() -> tuple[tuple[int, int], ...]:
    keys: list[tuple[int, int]] = []
    for task_id in range(112):
        sector = 56 + task_id // 16
        shard_index = task_id % 16
        tic = (
            722_078_603
            if (sector, shard_index) == (60, 0)
            else 1_000_000_000 + task_id
        )
        keys.append((sector, tic))
    return tuple(keys)


NUMERIC_ELIGIBLE_KEYS = _numeric_test_keys()
NUMERIC_EXCLUDED_KEYS = ((62, 2_000_000_000),)
NUMERIC_FULL_KEYS = NUMERIC_ELIGIBLE_KEYS + NUMERIC_EXCLUDED_KEYS
NUMERIC_ELIGIBLE_IDENTITY_SHA256 = observation_identity_sha256(
    NUMERIC_ELIGIBLE_KEYS
)
NUMERIC_EXCLUDED_IDENTITY_SHA256 = observation_identity_sha256(
    NUMERIC_EXCLUDED_KEYS
)
NUMERIC_FULL_IDENTITY_SHA256 = observation_identity_sha256(NUMERIC_FULL_KEYS)


@pytest.fixture(autouse=True)
def _bounded_numeric_partition(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        NUMERIC,
        "PRODUCTION_FULL_OBSERVATIONS",
        len(NUMERIC_FULL_KEYS),
    )
    monkeypatch.setattr(
        NUMERIC,
        "PRODUCTION_ELIGIBLE_OBSERVATIONS",
        len(NUMERIC_ELIGIBLE_KEYS),
    )
    monkeypatch.setattr(
        NUMERIC,
        "PRODUCTION_EXCLUDED_OBSERVATIONS",
        len(NUMERIC_EXCLUDED_KEYS),
    )
    monkeypatch.setattr(
        NUMERIC,
        "PRODUCTION_FULL_IDENTITY_SHA256",
        NUMERIC_FULL_IDENTITY_SHA256,
    )
    monkeypatch.setattr(
        NUMERIC,
        "PRODUCTION_ELIGIBLE_IDENTITY_SHA256",
        NUMERIC_ELIGIBLE_IDENTITY_SHA256,
    )
    monkeypatch.setattr(
        NUMERIC,
        "PRODUCTION_EXCLUDED_IDENTITY_SHA256",
        NUMERIC_EXCLUDED_IDENTITY_SHA256,
    )


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _metadata(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": _sha256(path),
    }


def _write_json(path: Path, value: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def _write_numeric_audit(path: Path) -> None:
    records: list[dict[str, Any]] = []
    for sector, tic in NUMERIC_ELIGIBLE_KEYS:
        records.append(
            {
                "ssl_observation_id": f"s{sector:04d}-tic{tic:016d}",
                "sector": sector,
                "tic": tic,
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
        )
    for sector, tic in NUMERIC_EXCLUDED_KEYS:
        records.append(
            {
                "ssl_observation_id": f"s{sector:04d}-tic{tic:016d}",
                "sector": sector,
                "tic": tic,
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
        )
    dtypes = dict(
        zip(
            NUMERIC._MODEL_INPUT_NUMERIC_AUDIT_COLUMNS,
            NUMERIC._MODEL_INPUT_NUMERIC_AUDIT_DTYPES,
            strict=True,
        )
    )
    frame = pd.DataFrame.from_records(
        records,
        columns=list(NUMERIC._MODEL_INPUT_NUMERIC_AUDIT_COLUMNS),
    ).astype(dtypes)
    path.write_bytes(NUMERIC.numeric_audit_parquet_bytes(frame))


def _fixture(tmp_path: Path) -> dict[str, Any]:
    authority_paths: dict[str, Path] = {}
    for name in VALIDATOR.AUTHORITY_NAMES:
        if name == "numeric_gate_release":
            continue
        path = tmp_path / f"{name}.dat"
        path.write_bytes(f"{name}\n".encode("ascii"))
        authority_paths[name] = path
    authority_metadata = {
        name: _metadata(path) for name, path in authority_paths.items()
    }
    expected_revision = "1" * 40
    native_freshness = numeric_native_freshness(expected_revision)
    numeric_gate_payload = {
        "schema_version": MODEL_INPUT_NUMERIC_RELEASE_SCHEMA,
        "release_binding": MODEL_INPUT_NUMERIC_RELEASE_BINDING,
        "native_freshness": native_freshness,
        "passed": True,
        "model_input_contract_version": MODEL_INPUT_CONTRACT_VERSION,
        "envelope_contract": MODEL_INPUT_NUMERIC_ENVELOPE_CONTRACT,
        "envelope": FULL_POOL_NUMERIC_ENVELOPE_V1.as_dict(),
        "envelope_canonical_sha256": (
            MODEL_INPUT_NUMERIC_ENVELOPE_SHA256
        ),
        "counts": {
            "full_observations": len(NUMERIC_FULL_KEYS),
            "eligible_observations": len(NUMERIC_ELIGIBLE_KEYS),
            "excluded_observations": len(NUMERIC_EXCLUDED_KEYS),
            "scanned_observations": len(NUMERIC_ELIGIBLE_KEYS),
            "failed_observations": 0,
            "native_shards": 112,
        },
        "identity_hashes": {
            "full": NUMERIC_FULL_IDENTITY_SHA256,
            "eligible": NUMERIC_ELIGIBLE_IDENTITY_SHA256,
            "excluded": NUMERIC_EXCLUDED_IDENTITY_SHA256,
        },
        "code_revision": expected_revision,
        "authority_bindings": {
            "ssl_registry": authority_metadata["registry"],
            "ssl_registry_summary": authority_metadata["registry_summary"],
            "native_registry": authority_metadata["native_registry"],
            "native_registry_summary": authority_metadata[
                "native_registry_summary"
            ],
            "native_release_summary": authority_metadata[
                "native_release_summary"
            ],
        },
    }
    numeric_audit = tmp_path / "model_input_numeric_audit.parquet"
    _write_numeric_audit(numeric_audit)
    numeric_audit_digest = _sha256(numeric_audit)
    Path(str(numeric_audit) + ".sha256").write_text(
        f"{numeric_audit_digest}  {numeric_audit.name}\n",
        encoding="ascii",
    )
    shard_reports: list[dict[str, Any]] = []
    for task_id in range(112):
        sector = 56 + task_id // 16
        shard_index = task_id % 16
        expected_key = NUMERIC_ELIGIBLE_KEYS[task_id]
        assert expected_key[0] == sector
        tic = expected_key[1]
        scanned = 1
        numeric_native = tmp_path / f"numeric_native_{task_id}.h5"
        numeric_native.write_bytes(f"native-{task_id}\n".encode("ascii"))
        native_binding = _metadata(numeric_native)
        identity = observation_identity_sha256([(sector, tic)])
        row = {
            "ssl_observation_id": f"s{sector:04d}-tic{tic:016d}",
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
        report_payload = {
            "schema_version": MODEL_INPUT_NUMERIC_AUDIT_SCHEMA,
            "release_binding": MODEL_INPUT_NUMERIC_RELEASE_BINDING,
            "native_freshness": dict(
                native_freshness
            ),
            "code_revision": expected_revision,
            "model_input_contract_version": MODEL_INPUT_CONTRACT_VERSION,
            "envelope_contract": MODEL_INPUT_NUMERIC_ENVELOPE_CONTRACT,
            "envelope": FULL_POOL_NUMERIC_ENVELOPE_V1.as_dict(),
            "envelope_canonical_sha256": (
                MODEL_INPUT_NUMERIC_ENVELOPE_SHA256
            ),
            "sector": sector,
            "shard_index": shard_index,
            "n_shards": 16,
            "passed": True,
            "counts": {
                "scanned_observations": scanned,
                "passed_observations": scanned,
                "failed_observations": 0,
            },
            "observation_identity_sha256": identity,
            "rows": [row],
            "native_h5": native_binding,
            "authority_bindings": numeric_gate_payload[
                "authority_bindings"
            ],
            "action": "audit_only_no_clip_no_exclusion",
        }
        report = tmp_path / f"numeric_{task_id}.json"
        _write_json(report, report_payload)
        report_digest = _sha256(report)
        Path(str(report) + ".sha256").write_text(
            f"{report_digest}  {report.name}\n",
            encoding="ascii",
        )
        shard_reports.append(
            {
                "sector": sector,
                "shard_index": shard_index,
                **_metadata(report),
                "native_h5": native_binding,
                "observation_identity_sha256": identity,
                "scanned_observations": scanned,
            }
        )
    numeric_gate_payload.update(
        {
            "outputs": {
                "numeric_audit": {
                    "path": str(numeric_audit),
                    "size_bytes": numeric_audit.stat().st_size,
                    "sha256": numeric_audit_digest,
                }
            },
            "shard_reports": shard_reports,
            "quality_bad_photometry_policy_verified": True,
            "float32_conversion_verified": True,
            "real_only": True,
            "labels_consumed": False,
            "injections_consumed": False,
            "action": "audit_only_no_clip_no_exclusion",
        }
    )
    numeric_gate_path = tmp_path / "numeric_gate_release.json"
    _write_json(numeric_gate_path, numeric_gate_payload)
    numeric_gate_digest = _sha256(numeric_gate_path)
    Path(str(numeric_gate_path) + ".sha256").write_text(
        f"{numeric_gate_digest}  {numeric_gate_path.name}\n",
        encoding="ascii",
    )
    authority_paths["numeric_gate_release"] = numeric_gate_path
    authority_metadata["numeric_gate_release"] = _metadata(
        numeric_gate_path
    )
    numeric_gate_summary = {
        "binding": authority_metadata["numeric_gate_release"],
        "release_binding": MODEL_INPUT_NUMERIC_RELEASE_BINDING,
        "schema_version": MODEL_INPUT_NUMERIC_RELEASE_SCHEMA,
        "envelope_canonical_sha256": (
            MODEL_INPUT_NUMERIC_ENVELOPE_SHA256
        ),
        "counts": numeric_gate_payload["counts"],
        "identity_hashes": numeric_gate_payload["identity_hashes"],
        "code_revision": expected_revision,
        "passed": True,
    }
    native_path = tmp_path / "native_0.h5"
    with h5py.File(native_path, "w") as h5:
        h5.attrs["contract_version"] = FULL_POOL_NATIVE_CONTRACT_VERSION
        h5.attrs["release_binding"] = (
            MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS["release_binding"]
        )
        h5.attrs["expected_git_sha"] = expected_revision
        h5.attrs["builder_contract_version"] = (
            MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS[
                "builder_contract_version"
            ]
        )
        h5.attrs["builder_code_sha256"] = (
            MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS[
                "builder_code_sha256"
            ]
        )
        h5.attrs["detrend_contract_version"] = (
            MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS[
                "detrend_contract_version"
            ]
        )
        h5.attrs["detrend_config_sha256"] = (
            MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS[
                "detrend_config_sha256"
            ]
        )
        h5.attrs["detrend_quality_source"] = "final_effective_quality"
        h5.attrs["raw_photometry_only"] = 1
        h5.attrs["compact_adp_photometry_reused"] = 0
        h5.attrs["compact_adp_flux_reused"] = 0
        h5.attrs["periodogram_n"] = MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS[
            "periodogram_n"
        ]
    selection = {
        "selection_schema_version": FULLPOOL_SSL_SELECTION_SCHEMA,
        "held_out_fold": VALIDATOR.SMOKE_FOLD,
        "n_registry_observations": PRODUCTION_FULL_OBSERVATIONS,
        "n_eligible_observations": PRODUCTION_ELIGIBLE_OBSERVATIONS,
        "n_eligible_tics": 125_620,
        "n_held_observations": 0,
        "n_held_tics": 0,
        "n_selected_observations": VALIDATOR.SMOKE_MAX_ROWS,
        "n_selected_tics": VALIDATOR.SMOKE_MAX_ROWS,
        "max_rows": VALIDATOR.SMOKE_MAX_ROWS,
        "required_observation_ids": list(
            VALIDATOR.SMOKE_REQUIRED_OBSERVATION_IDS
        ),
        "n_required_observations": len(
            VALIDATOR.SMOKE_REQUIRED_OBSERVATION_IDS
        ),
        "required_observations_selected": True,
        "selected_rows_sha256": "2" * 64,
        "selected_tics_sha256": "3" * 64,
        "tic_disjoint": {
            "held_fold_tics": True,
            "fixed_test_tics": True,
            "reserved_prospective_tics": True,
        },
    }
    contract = {
        "schema_version": FULLPOOL_SSL_RUN_CONTRACT_SCHEMA,
        "run_id": FULLPOOL_SSL_RUN_ID,
        "model_input_contract_version": MODEL_INPUT_CONTRACT_VERSION,
        "model_input_numeric_envelope_contract": (
            MODEL_INPUT_NUMERIC_ENVELOPE_CONTRACT
        ),
        "numeric_release_binding": MODEL_INPUT_NUMERIC_RELEASE_BINDING,
        "native_freshness": dict(native_freshness),
        "numeric_gate_release": numeric_gate_summary,
        "checkpoint_namespace": FULLPOOL_SSL_CHECKPOINT_NAMESPACE,
        "fold": VALIDATOR.SMOKE_FOLD,
        "registry_path": authority_metadata["registry"]["path"],
        "registry_sha256": authority_metadata["registry"]["sha256"],
        "registry_summary_path": authority_metadata["registry_summary"]["path"],
        "registry_summary_sha256": authority_metadata["registry_summary"]["sha256"],
        "training_authority": {
            "schema_version": FULLPOOL_SSL_TRAINING_AUTHORITY_SCHEMA,
            "production_lock_passed": True,
            "partition": {
                "full": {
                    "count": PRODUCTION_FULL_OBSERVATIONS,
                    "identity_sha256": PRODUCTION_FULL_IDENTITY_SHA256,
                },
                "eligible": {
                    "count": PRODUCTION_ELIGIBLE_OBSERVATIONS,
                    "identity_sha256": PRODUCTION_ELIGIBLE_IDENTITY_SHA256,
                },
                "excluded": {
                    "count": PRODUCTION_EXCLUDED_OBSERVATIONS,
                    "identity_sha256": PRODUCTION_EXCLUDED_IDENTITY_SHA256,
                },
            },
            "native_mapping_sha256": "4" * 64,
            "model_input_contract_version": MODEL_INPUT_CONTRACT_VERSION,
            "model_input_numeric_envelope_contract": (
                MODEL_INPUT_NUMERIC_ENVELOPE_CONTRACT
            ),
            "numeric_gate_release": numeric_gate_summary,
            "native_freshness": dict(
                native_freshness
            ),
            "source_provenance_verified": True,
            "authority_bindings": authority_metadata,
        },
        "selection_audit": selection,
        "native_files": [
            {
                "native_h5_path": str(native_path.resolve()),
                "native_h5_sha256": _sha256(native_path),
                "native_h5_size_bytes": native_path.stat().st_size,
                "native_contract_version": FULL_POOL_NATIVE_CONTRACT_VERSION,
                "native_release_binding": (
                    MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS["release_binding"]
                ),
                "native_expected_git_sha": expected_revision,
                "native_builder_contract_version": (
                    MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS[
                        "builder_contract_version"
                    ]
                ),
                "native_builder_code_sha256": (
                    MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS[
                        "builder_code_sha256"
                    ]
                ),
                "native_detrend_contract_version": (
                    MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS[
                        "detrend_contract_version"
                    ]
                ),
                "native_detrend_config_sha256": (
                    MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS[
                        "detrend_config_sha256"
                    ]
                ),
                "native_detrend_quality_source": "final_effective_quality",
                "raw_photometry_only": True,
                "compact_adp_photometry_reused": False,
                "compact_adp_flux_reused": False,
                "periodogram_n": MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS[
                    "periodogram_n"
                ],
                "hash_verified_now": True,
                "root_contract_verified_now": True,
                "group_identities_verified_now": True,
                "n_selected_observations": VALIDATOR.SMOKE_MAX_ROWS,
            }
        ],
        "native_hashes_verified_now": True,
        "native_root_contracts_verified_now": True,
        "native_group_identities_verified_now": True,
        "epochs": 1,
        "batch_size": VALIDATOR.SMOKE_BATCH_SIZE,
        "workers": VALIDATOR.SMOKE_WORKERS,
        "seed": VALIDATOR.SMOKE_SEED,
        "learning_rate": VALIDATOR.SMOKE_LEARNING_RATE,
        "weight_decay": VALIDATOR.SMOKE_WEIGHT_DECAY,
        "checkpoint_every": 1,
        "require_cuda": True,
        "max_rows": VALIDATOR.SMOKE_MAX_ROWS,
        "required_observation_ids": list(
            VALIDATOR.SMOKE_REQUIRED_OBSERVATION_IDS
        ),
        "labels_loaded": False,
        "fixed_test_tensors_constructed": False,
        "prospective_sector_tensors_constructed": False,
        "embedding_export": False,
        "neighbor_probe": False,
        "code_revision": expected_revision,
    }
    fold_dir = (
        tmp_path
        / "model_runs"
        / FULLPOOL_SSL_MODEL_NAMESPACE
        / "smoke"
        / "one_epoch"
        / "encoder_pretraining"
        / "fold_2"
    )
    fold_dir.mkdir(parents=True)
    contract_path = fold_dir / "run_contract.json"
    _write_json(contract_path, contract)
    checkpoint_path = fold_dir / "encoder.pt"
    history_path = fold_dir / "history.csv"
    history_path.write_text(
        "epoch,n_observations,loss\n1,4096,0.5\n",
        encoding="ascii",
    )
    checkpoint = {
        "schema_version": VALIDATOR.FULLPOOL_SSL_CHECKPOINT_SCHEMA,
        "run_id": FULLPOOL_SSL_RUN_ID,
        "model_input_contract_version": MODEL_INPUT_CONTRACT_VERSION,
        "model_input_numeric_envelope_contract": (
            MODEL_INPUT_NUMERIC_ENVELOPE_CONTRACT
        ),
        "numeric_release_binding": MODEL_INPUT_NUMERIC_RELEASE_BINDING,
        "native_freshness": dict(native_freshness),
        "numeric_gate_release": numeric_gate_summary,
        "checkpoint_namespace": FULLPOOL_SSL_CHECKPOINT_NAMESPACE,
        "fold": VALIDATOR.SMOKE_FOLD,
        "run_contract_sha256": _sha256(contract_path),
        "selection_audit": selection,
        "epochs": 1,
        "history": [{"epoch": 1, "n_observations": 4096, "loss": 0.5}],
        "encoder_state_dict": {"weight": torch.ones(2, 2)},
        "projection_state_dict": {"weight": torch.ones(2, 2)},
    }
    torch.save(checkpoint, checkpoint_path)
    summary = {
        "schema_version": FULLPOOL_SSL_SUMMARY_SCHEMA,
        "run_id": FULLPOOL_SSL_RUN_ID,
        "model_input_contract_version": MODEL_INPUT_CONTRACT_VERSION,
        "model_input_numeric_envelope_contract": (
            MODEL_INPUT_NUMERIC_ENVELOPE_CONTRACT
        ),
        "numeric_release_binding": MODEL_INPUT_NUMERIC_RELEASE_BINDING,
        "native_freshness": dict(native_freshness),
        "numeric_gate_release": numeric_gate_summary,
        "checkpoint_namespace": FULLPOOL_SSL_CHECKPOINT_NAMESPACE,
        "fold": VALIDATOR.SMOKE_FOLD,
        "run_contract": str(contract_path.resolve()),
        "run_contract_sha256": _sha256(contract_path),
        "checkpoint": str(checkpoint_path.resolve()),
        "checkpoint_sha256": _sha256(checkpoint_path),
        "history": str(history_path.resolve()),
        "history_sha256": _sha256(history_path),
        "completed_epochs": 1,
        "global_step": VALIDATOR.SMOKE_GLOBAL_STEPS,
        "selection_audit": selection,
        "fixed_test_status": "host_excluded_tensors_not_constructed",
        "prospective_status": "host_excluded_tensors_not_constructed",
        "labels_loaded": False,
        "automatic_production_promotion": False,
    }
    summary_path = fold_dir / "summary.json"
    _write_json(summary_path, summary)
    return {
        "summary_path": summary_path,
        "contract_path": contract_path,
        "authority_paths": authority_paths,
        "expected_revision": expected_revision,
        "native_path": native_path,
    }


def _validate(
    case: dict[str, Any],
    *,
    expected_fold: int = VALIDATOR.SMOKE_FOLD,
) -> dict[str, Any]:
    paths = case["authority_paths"]
    return VALIDATOR.validate_teacher_ssl_fullpool_smoke(
        summary_path=case["summary_path"],
        expected_code_revision=case["expected_revision"],
        eligibility_exclusions_path=paths["eligibility_exclusions"],
        eligibility_summary_path=paths["eligibility_summary"],
        native_registry_path=paths["native_registry"],
        native_registry_summary_path=paths["native_registry_summary"],
        native_release_summary_path=paths["native_release_summary"],
        registry_path=paths["registry"],
        registry_summary_path=paths["registry_summary"],
        numeric_gate_release_path=paths["numeric_gate_release"],
        expected_fold=expected_fold,
    )


def _rewrite_contract(
    case: dict[str, Any],
    mutate: Callable[[dict[str, Any]], None],
    *,
    rebind_summary: bool,
) -> None:
    contract = json.loads(case["contract_path"].read_text(encoding="utf-8"))
    mutate(contract)
    _write_json(case["contract_path"], contract)
    if rebind_summary:
        summary = json.loads(case["summary_path"].read_text(encoding="utf-8"))
        summary["run_contract_sha256"] = _sha256(case["contract_path"])
        _write_json(case["summary_path"], summary)


def _rewrite_numeric_gate_release(
    case: dict[str, Any],
    mutate: Callable[[dict[str, Any]], None],
) -> None:
    path = case["authority_paths"]["numeric_gate_release"]
    release = json.loads(path.read_text(encoding="utf-8"))
    mutate(release)
    _write_json(path, release)


def test_smoke_validator_accepts_exact_production_result(tmp_path: Path) -> None:
    case = _fixture(tmp_path)

    audit = _validate(case)

    assert audit["passed"] is True
    assert audit["max_rows"] == 4096
    assert set(audit["authority_bindings"]) == set(VALIDATOR.AUTHORITY_NAMES)


def test_smoke_validator_rejects_any_fold_other_than_two(
    tmp_path: Path,
) -> None:
    case = _fixture(tmp_path)

    with pytest.raises(ValueError, match="locked smoke fold must equal 2"):
        _validate(case, expected_fold=0)


def test_smoke_validator_recomputes_run_contract_hash(tmp_path: Path) -> None:
    case = _fixture(tmp_path)
    _rewrite_contract(
        case,
        lambda contract: contract.update({"workers": 99}),
        rebind_summary=False,
    )

    with pytest.raises(ValueError, match="run-contract SHA-256"):
        _validate(case)


def test_smoke_validator_rehashes_all_authority_artifacts(
    tmp_path: Path,
) -> None:
    case = _fixture(tmp_path)
    case["authority_paths"]["native_release_summary"].write_bytes(b"changed\n")

    with pytest.raises(
        ValueError,
        match="native_release_summary",
    ):
        _validate(case)


def test_smoke_validator_rejects_failed_numeric_gate_release(
    tmp_path: Path,
) -> None:
    case = _fixture(tmp_path)
    _rewrite_numeric_gate_release(
        case,
        lambda release: release.update({"passed": False}),
    )

    with pytest.raises(ValueError, match="numeric-gate release did not pass"):
        _validate(case)


def test_smoke_validator_rejects_model_input_contract_drift(
    tmp_path: Path,
) -> None:
    case = _fixture(tmp_path)
    _rewrite_contract(
        case,
        lambda contract: contract.update(
            {"model_input_contract_version": "obsolete"}
        ),
        rebind_summary=True,
    )

    with pytest.raises(ValueError, match="wrong model-input contract"):
        _validate(case)


def test_smoke_validator_rejects_numeric_authority_hash_tamper(
    tmp_path: Path,
) -> None:
    case = _fixture(tmp_path)
    case["authority_paths"]["registry"].write_bytes(b"tampered-registry\n")

    with pytest.raises(ValueError, match="authority_bindings.ssl_registry"):
        _validate(case)


def test_smoke_validator_rejects_valid_but_stale_freshness_sha(
    tmp_path: Path,
) -> None:
    case = _fixture(tmp_path)
    summary = json.loads(case["summary_path"].read_text(encoding="utf-8"))
    summary["native_freshness"]["expected_git_sha"] = "2" * 40
    _write_json(case["summary_path"], summary)

    with pytest.raises(
        ValueError,
        match="validated numeric-release freshness",
    ):
        _validate(case)


def test_smoke_validator_rehashes_selected_native_h5(tmp_path: Path) -> None:
    case = _fixture(tmp_path)
    with case["native_path"].open("r+b") as handle:
        handle.seek(-1, 2)
        original = handle.read(1)
        handle.seek(-1, 2)
        handle.write(bytes([original[0] ^ 1]))

    with pytest.raises(ValueError, match="native file 0 hash binding changed"):
        _validate(case)


def test_smoke_validator_rejects_nonfinite_history(tmp_path: Path) -> None:
    case = _fixture(tmp_path)
    summary = json.loads(case["summary_path"].read_text(encoding="utf-8"))
    history_path = Path(summary["history"])
    history_path.write_text(
        "epoch,n_observations,loss\n1,4096,nan\n",
        encoding="ascii",
    )
    summary["history_sha256"] = _sha256(history_path)
    _write_json(case["summary_path"], summary)

    with pytest.raises(ValueError, match="history loss is non-finite"):
        _validate(case)


def test_smoke_validator_rejects_nonfinite_checkpoint_state(
    tmp_path: Path,
) -> None:
    case = _fixture(tmp_path)
    summary = json.loads(case["summary_path"].read_text(encoding="utf-8"))
    checkpoint_path = Path(summary["checkpoint"])
    checkpoint = torch.load(
        checkpoint_path,
        map_location="cpu",
        weights_only=False,
    )
    checkpoint["encoder_state_dict"]["weight"][0, 0] = float("inf")
    torch.save(checkpoint, checkpoint_path)
    summary["checkpoint_sha256"] = _sha256(checkpoint_path)
    _write_json(case["summary_path"], summary)

    with pytest.raises(ValueError, match="encoder_state_dict.weight is non-finite"):
        _validate(case)


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (
            lambda contract: contract["training_authority"].update(
                {"production_lock_passed": False}
            ),
            "production lock failed",
        ),
        (
            lambda contract: contract["training_authority"].update(
                {"source_provenance_verified": False}
            ),
            "source provenance failed",
        ),
        (
            lambda contract: contract.update(
                {"native_group_identities_verified_now": False}
            ),
            "all three native verification gates",
        ),
        (
            lambda contract: contract["training_authority"]["partition"][
                "eligible"
            ].update({"count": PRODUCTION_ELIGIBLE_OBSERVATIONS - 1}),
            "production partition eligible count",
        ),
        (
            lambda contract: contract.update({"labels_loaded": True}),
            "isolation failed",
        ),
        (
            lambda contract: contract.update({"max_rows": 4095}),
            "smoke max rows must equal 4096",
        ),
        (
            lambda contract: contract.update({"batch_size": 32}),
            "smoke batch size must equal 64",
        ),
        (
            lambda contract: contract.update({"seed": 1}),
            f"smoke seed must equal {VALIDATOR.SMOKE_SEED}",
        ),
        (
            lambda contract: contract.update(
                {"required_observation_ids": []}
            ),
            "locked diagnostic observation",
        ),
        (
            lambda contract: contract["selection_audit"].update(
                {"required_observations_selected": False}
            ),
            "omitted the locked diagnostic observation",
        ),
    ],
)
def test_smoke_validator_rejects_relaxed_gates(
    tmp_path: Path,
    mutate: Callable[[dict[str, Any]], None],
    message: str,
) -> None:
    case = _fixture(tmp_path)
    _rewrite_contract(case, mutate, rebind_summary=True)

    with pytest.raises(ValueError, match=message):
        _validate(case)
