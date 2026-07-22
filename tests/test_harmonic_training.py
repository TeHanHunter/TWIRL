from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import h5py
import numpy as np
import pandas as pd
import pytest

import twirl.vetting.harmonic_training as harmonic_training
from twirl.vetting.harmonic_training import (
    ENCODER_PRETRAINING_CACHE_SCHEMA,
    SELECTED_CHECKPOINT_MANIFEST_SCHEMA,
    TEACHER_NATIVE_V2_CHECKPOINT_NAMESPACE,
    build_teacher_input_provenance,
    classification_metrics,
    expected_calibration_error,
    injection_truth_human_audit,
    validate_encoder_pretraining_cache,
    validate_teacher_input_provenance,
    verify_deployed_teacher_checkpoint_manifest,
)
from twirl.vetting.harmonic_cnn import MODEL_VERSION
from twirl.vetting.harmonic_inputs import RAW_PAIR_CONTRACT_VERSION


def test_multiclass_metrics_report_balanced_accuracy_and_confusion() -> None:
    truth = np.asarray([0, 0, 1, 1, 2, 2, 3, 3])
    predicted = np.asarray([0, 0, 1, 0, 2, 2, 3, 1])
    probability = np.full((len(truth), 4), 0.01)
    probability[np.arange(len(truth)), predicted] = 0.97

    metrics = classification_metrics(
        truth,
        probability,
        classes=("planet", "eb", "variable", "other"),
    )

    assert np.isclose(metrics["accuracy"], 0.75)
    assert np.isclose(metrics["balanced_accuracy"], 0.75)
    assert metrics["confusion_matrix"][1][0] == 1
    assert metrics["per_class"]["planet"]["recall"] == 1.0


def test_calibration_error_is_finite_and_small_for_correct_probabilities() -> None:
    truth = np.asarray([0, 1, 0, 1])
    probability = np.asarray(
        [[0.9, 0.1], [0.1, 0.9], [0.8, 0.2], [0.2, 0.8]],
        dtype=float,
    )

    calibration = expected_calibration_error(truth, probability, n_bins=5)

    assert np.isfinite(calibration["ece"])
    assert 0.0 <= calibration["ece"] <= 0.2


def test_injection_truth_audit_keeps_human_and_truth_separate(tmp_path) -> None:
    rows = pd.DataFrame(
        {
            "review_id": ["inj:1", "inj:2", "real:3"],
            "source_kind": ["injection_recovery", "injection_recovery", "real_candidate"],
            "injection_id": ["one", "two", ""],
            "tic": [1, 2, 3],
            "human_label": ["planet_like", "uncertain", "planet_like"],
            "truth_period_d": [2.0, 1.37, np.nan],
            "period_d": [1.0, 1.0, 1.0],
            "truth_radius_rearth": [2.0, 4.0, np.nan],
            "truth_model_depth": [0.1, 0.2, np.nan],
            "tmag": [18.0, 19.0, 17.0],
            "adp_sml_sde": [20.0, 3.0, 10.0],
            "boundary_period_bin_key": ["p1", "p2", ""],
            "boundary_radius_bin_key": ["r1", "r2", ""],
            "boundary_tmag_bin_key": ["m1", "m2", ""],
        }
    )

    summary = injection_truth_human_audit(rows, out_dir=tmp_path)
    audit = pd.read_csv(tmp_path / "injection_truth_human_rows.csv")

    assert summary["n_injections"] == 2
    assert summary["n_human_visible_planet"] == 1
    assert summary["n_bls_truth_period_or_harmonic_match"] == 1
    assert audit.loc[audit["injection_id"].eq("one"), "truth_to_bls_period_factor"].iloc[0] == 2.0
    assert audit.loc[audit["injection_id"].eq("two"), "human_label"].iloc[0] == "uncertain"


def _native_v2_provenance() -> dict[str, str]:
    return {
        "checkpoint_namespace": TEACHER_NATIVE_V2_CHECKPOINT_NAMESPACE,
        "input_contract_version": RAW_PAIR_CONTRACT_VERSION,
        "native_h5_sha256": "a" * 64,
        "training_table_sha256": "b" * 64,
    }


def _native_v2_pretraining_cache() -> dict[str, object]:
    return {
        **_native_v2_provenance(),
        "pretraining_cache_schema": ENCODER_PRETRAINING_CACHE_SCHEMA,
        "pretraining_profile": "shape_plus_raw_chronology",
        "model_config": {"metadata_dim": 3},
        "train_config": {"learning_rate": 3.0e-4},
        "seed": 56,
        "pretraining_epochs": 20,
        "state_dict": {},
    }


def test_native_v2_pretraining_rejects_v1_cache() -> None:
    cache = _native_v2_pretraining_cache()
    cache["input_contract_version"] = "s56_adp_raw_pair_v1"

    with pytest.raises(RuntimeError, match="stale or unprovenanced"):
        validate_encoder_pretraining_cache(
            cache,
            expected_provenance=_native_v2_provenance(),
            model_config={"metadata_dim": 3},
            train_config={"learning_rate": 3.0e-4},
            seed=56,
            epochs=20,
        )


@pytest.mark.parametrize(
    ("field", "stale_value"),
    (
        ("native_h5_sha256", "c" * 64),
        ("training_table_sha256", "d" * 64),
    ),
)
def test_native_v2_pretraining_rejects_stale_input_hash(
    field: str,
    stale_value: str,
) -> None:
    cache = _native_v2_pretraining_cache()
    cache[field] = stale_value

    with pytest.raises(RuntimeError, match=field):
        validate_encoder_pretraining_cache(
            cache,
            expected_provenance=_native_v2_provenance(),
            model_config={"metadata_dim": 3},
            train_config={"learning_rate": 3.0e-4},
            seed=56,
            epochs=20,
        )


def test_native_v2_final_checkpoint_rejects_missing_training_binding() -> None:
    checkpoint = _native_v2_provenance()
    checkpoint.pop("training_table_sha256")

    with pytest.raises(RuntimeError, match="training_table_sha256"):
        validate_teacher_input_provenance(
            checkpoint,
            expected=_native_v2_provenance(),
            artifact="final teacher checkpoint",
        )


def test_native_v2_provenance_requires_native_table_hash_match(tmp_path) -> None:
    training_table = tmp_path / "training.csv"
    native_h5 = tmp_path / "native.h5"
    training_table.write_text("review_id,label\na,planet_like\n")
    table_sha256 = hashlib.sha256(training_table.read_bytes()).hexdigest()
    with h5py.File(native_h5, "w") as h5:
        h5.attrs["contract_version"] = RAW_PAIR_CONTRACT_VERSION
        h5.attrs["training_table_sha256"] = table_sha256

    provenance = build_teacher_input_provenance(
        training_table=training_table,
        native_h5=native_h5,
    )

    assert provenance["training_table_sha256"] == table_sha256
    assert len(provenance["native_h5_sha256"]) == 64
    training_table.write_text("review_id,label\na,eclipse_contact\n")
    with pytest.raises(RuntimeError, match="exact training table"):
        build_teacher_input_provenance(
            training_table=training_table,
            native_h5=native_h5,
        )


def test_native_v2_provenance_rejects_input_changed_while_hashing(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    training_table = tmp_path / "training.csv"
    native_h5 = tmp_path / "native.h5"
    training_table.write_text("review_id,label\na,planet_like\n")
    table_sha256 = hashlib.sha256(training_table.read_bytes()).hexdigest()
    with h5py.File(native_h5, "w") as h5:
        h5.attrs["contract_version"] = RAW_PAIR_CONTRACT_VERSION
        h5.attrs["training_table_sha256"] = table_sha256

    real_sha256 = harmonic_training._file_sha256
    calls = 0

    def changing_sha256(path: Path) -> str:
        nonlocal calls
        calls += 1
        if calls == 3:
            return "f" * 64
        return real_sha256(path)

    monkeypatch.setattr(harmonic_training, "_file_sha256", changing_sha256)
    with pytest.raises(RuntimeError, match="training table changed"):
        build_teacher_input_provenance(
            training_table=training_table,
            native_h5=native_h5,
        )


def test_deployed_manifest_rejects_changed_checkpoint(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    profile = "shape_plus_periodogram_bls"
    model_root = tmp_path / "checkpoints"
    checkpoint_root = model_root / profile
    provenance = _native_v2_provenance()
    payloads = {}
    records = []
    for fold in range(5):
        path = checkpoint_root / f"fold_{fold}" / "teacher.pt"
        path.parent.mkdir(parents=True)
        path.write_bytes(f"native-v2-checkpoint-{fold}".encode())
        payloads[str(path)] = {
            **provenance,
            "model_version": MODEL_VERSION,
            "profile": profile,
            "fold": fold,
            "encoder_pretraining_cache_schema": ENCODER_PRETRAINING_CACHE_SCHEMA,
            "encoder_pretraining_sha256": "e" * 64,
        }
        records.append(
            {
                "fold": fold,
                "path": f"{profile}/fold_{fold}/teacher.pt",
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "encoder_pretraining_sha256": "e" * 64,
            }
        )
    manifest_path = model_root / "selected_checkpoint_manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "schema_version": SELECTED_CHECKPOINT_MANIFEST_SCHEMA,
                "model_version": MODEL_VERSION,
                **provenance,
                "selected_profile": profile,
                "checkpoints": records,
            }
        )
    )
    monkeypatch.setitem(
        sys.modules,
        "torch",
        SimpleNamespace(load=lambda path, **kwargs: payloads[str(Path(path))]),
    )

    verified = verify_deployed_teacher_checkpoint_manifest(
        manifest_path=manifest_path,
        checkpoint_root=checkpoint_root,
    )
    assert verified["checkpoint_namespace"] == (
        TEACHER_NATIVE_V2_CHECKPOINT_NAMESPACE
    )

    changed = checkpoint_root / "fold_3" / "teacher.pt"
    changed.write_bytes(changed.read_bytes() + b"-stale")
    with pytest.raises(RuntimeError, match="deployed checkpoint SHA256 mismatch"):
        verify_deployed_teacher_checkpoint_manifest(
            manifest_path=manifest_path,
            checkpoint_root=checkpoint_root,
        )
