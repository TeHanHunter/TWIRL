from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
import sys

import h5py
import numpy as np
import pandas as pd
import pytest

import twirl.vetting.teacher_v3_training as teacher_v3_training
from twirl.vetting.harmonic_cnn import HarmonicTrainConfig, MODEL_VERSION
from twirl.vetting.harmonic_inputs import RAW_PAIR_CONTRACT_VERSION
from twirl.vetting.teacher_v3 import (
    TEACHER_V3_CORPUS_VERSION,
    TEACHER_V3_RUN_NAME,
)
from twirl.vetting.teacher_v3_training import (
    TEACHER_V3_BASELINE_PROFILE,
    TEACHER_V3_CHECKPOINT_MANIFEST_SCHEMA,
    TEACHER_V3_PRIMARY_LABEL_POLICY,
    TEACHER_V3_RUN_ID,
    TEACHER_V3_UNCERTAIN_MASKED_POLICY,
    build_teacher_v3_metadata_provenance,
    build_teacher_v3_native_manifest,
    calibrate_teacher_v3_profile_oof,
    import_verified_teacher_v3_metadata_baseline,
    prepare_teacher_v3_uncertain_masked_rows,
    run_teacher_v3_metadata_baseline,
    teacher_v3_tic_cluster_bootstrap,
    validate_teacher_v3_training_table,
)


def _release_rows(tmp_path: Path) -> pd.DataFrame:
    sectors = list(range(56, 63))
    tics = [1000 + value for value in sectors]
    split = ["test", *(["development"] * 6)]
    fold = [-1, 0, 1, 2, 3, 4, 0]
    return pd.DataFrame(
        {
            "review_id": [f"review:{value}" for value in sectors],
            "teacher_v3_observation_id": [
                f"observation:{value}" for value in sectors
            ],
            "teacher_run_name": [TEACHER_V3_RUN_NAME] * 7,
            "teacher_v3_corpus_version": [
                TEACHER_V3_CORPUS_VERSION
            ]
            * 7,
            "teacher_architecture_version": [MODEL_VERSION] * 7,
            "teacher_v3_training_include": [True] * 7,
            "sector": sectors,
            "tic": tics,
            "fixed_split": split,
            "cv_fold": fold,
            "native_h5_path": [
                str(tmp_path / f"s{sector}.h5") for sector in sectors
            ],
            "native_group_path": [
                f"targets/{tic:016d}" for tic in tics
            ],
            "is_injected_row": [False] * 7,
            "human_label": ["uncertain"] * 7,
        }
    )


def test_teacher_v3_training_table_validates_frozen_release(
    tmp_path: Path,
) -> None:
    rows = _release_rows(tmp_path)

    validate_teacher_v3_training_table(rows)

    rows.loc[0, "teacher_architecture_version"] = "teacher_v3_new_cnn"
    try:
        validate_teacher_v3_training_table(rows)
    except ValueError as exc:
        assert "teacher_architecture_version" in str(exc)
    else:  # pragma: no cover - defensive assertion.
        raise AssertionError("wrong architecture was accepted")


def test_metadata_provenance_needs_no_native_file(tmp_path: Path) -> None:
    table = tmp_path / "training.csv"
    table.write_text("review_id\none\n")

    provenance = build_teacher_v3_metadata_provenance(
        training_table=table,
    )

    assert provenance["native_input_mode"] == "metadata_only_no_hdf5"
    assert len(provenance["training_table_sha256"]) == 64
    assert len(provenance["native_h5_sha256"]) == 64


def test_teacher_v3_native_manifest_hashes_all_seven_files(
    tmp_path: Path,
    monkeypatch,
) -> None:
    rows = _release_rows(tmp_path)
    for row in rows.to_dict("records"):
        path = Path(row["native_h5_path"])
        with h5py.File(path, "w") as h5:
            h5.attrs["contract_version"] = RAW_PAIR_CONTRACT_VERSION
            group = h5.require_group(row["native_group_path"])
            group.attrs["sector"] = int(row["sector"])
            group.attrs["tic"] = int(row["tic"])
    registry = tmp_path / "registry.csv"
    registry.write_text("placeholder\n")
    registry_summary = tmp_path / "registry.summary.json"
    registry_summary.write_text("{}\n")
    monkeypatch.setattr(
        teacher_v3_training,
        "_validate_real_registry",
        lambda **kwargs: {"passed": True},
    )
    monkeypatch.setattr(
        teacher_v3_training,
        "verify_raw_pair_contract",
        lambda path, **kwargs: {
            "passed": True,
            "failures": [],
            "counts": {"targets": 1, "injections": 0},
            "external_quality_counts": {},
        },
    )

    manifest = build_teacher_v3_native_manifest(
        rows=rows,
        registry_path=registry,
        registry_summary_path=registry_summary,
    )

    assert manifest["run_id"] == TEACHER_V3_RUN_ID
    assert manifest["sectors"] == list(range(56, 63))
    assert len(manifest["native_files"]) == 7
    assert all(len(record["sha256"]) == 64 for record in manifest["native_files"])


def test_pooled_oof_calibration_sets_one_temperature_on_all_folds(
    tmp_path: Path,
    monkeypatch,
) -> None:
    rows: list[dict[str, object]] = []
    for fold in range(5):
        fold_dir = tmp_path / TEACHER_V3_BASELINE_PROFILE / f"fold_{fold}"
        fold_dir.mkdir(parents=True)
        predictions = pd.DataFrame(
            {
                "review_id": [f"r{fold}a", f"r{fold}b"],
                "morphology_target": [0, 1],
                "logit_planet_like": [2.0, 0.0],
                "logit_eclipse_contact": [0.0, 2.0],
                "logit_smooth_variable": [-1.0, -1.0],
                "logit_other": [-2.0, -2.0],
            }
        )
        predictions.to_csv(
            fold_dir / "validation_predictions.csv",
            index=False,
        )
        (fold_dir / "teacher.pt").write_bytes(b"pending")
        (fold_dir / "metrics.json").write_text("{}\n")
        for review_id, human_label in (
            (f"r{fold}a", "planet_like"),
            (f"r{fold}b", "eclipsing_binary_or_pceb"),
        ):
            rows.append(
                {
                    "review_id": review_id,
                    "fixed_split": "development",
                    "is_injected_row": False,
                    "human_label": human_label,
                }
            )
    saved: list[dict[str, object]] = []

    def fake_load(*args, **kwargs):
        return {
            "temperature_calibration_scope": (
                "pending_pooled_oof_development"
            )
        }

    def fake_save(payload, path):
        saved.append(dict(payload))
        Path(path).write_bytes(b"calibrated")

    monkeypatch.setitem(
        sys.modules,
        "torch",
        SimpleNamespace(load=fake_load, save=fake_save),
    )
    monkeypatch.setattr(
        teacher_v3_training,
        "fit_temperature",
        lambda logits, truth: 2.5,
    )
    monkeypatch.setattr(
        teacher_v3_training,
        "validate_teacher_input_provenance",
        lambda *args, **kwargs: None,
    )

    result = calibrate_teacher_v3_profile_oof(
        rows=pd.DataFrame(rows),
        out_dir=tmp_path,
        profile=TEACHER_V3_BASELINE_PROFILE,
        input_provenance={
            "checkpoint_namespace": "test",
            "input_contract_version": "test",
            "native_h5_sha256": "a" * 64,
            "training_table_sha256": "b" * 64,
        },
    )

    assert result["temperature"] == 2.5
    assert len(saved) == 5
    assert {item["temperature"] for item in saved} == {2.5}
    assert {
        item["temperature_calibration_scope"] for item in saved
    } == {"pooled_oof_development"}


def test_metadata_baseline_skips_native_and_encoder_pretraining(
    tmp_path: Path,
    monkeypatch,
) -> None:
    source = _release_rows(tmp_path)
    training_table = tmp_path / "training.csv"
    source.to_csv(training_table, index=False)
    calls: list[dict[str, object]] = []
    monkeypatch.setattr(
        teacher_v3_training,
        "prepare_harmonic_training_rows",
        lambda rows, seed: rows.copy(),
    )
    monkeypatch.setattr(
        teacher_v3_training,
        "injection_truth_human_audit",
        lambda *args, **kwargs: {"n_injections": 0},
    )
    monkeypatch.setattr(
        teacher_v3_training,
        "_train_one_fold",
        lambda **kwargs: calls.append(kwargs) or {},
    )
    metrics = {
        "all": {
            "macro_f1": 0.5,
            "balanced_accuracy": 0.5,
            "per_class": {},
        },
        "real": {
            "per_class": {
                "planet_like": {"recall": 0.5},
                "eclipse_contact": {"recall": 0.5},
                "smooth_variable": {"recall": 0.5},
                "other": {"recall": 0.5},
            }
        },
    }
    monkeypatch.setattr(
        teacher_v3_training,
        "calibrate_teacher_v3_profile_oof",
        lambda **kwargs: {
            "profile": TEACHER_V3_BASELINE_PROFILE,
            "temperature": 1.5,
            "metrics": {
                "morphology_by_source": metrics,
                "calibration": {"all": {"ece": 0.1}},
            },
        },
    )
    monkeypatch.setattr(
        teacher_v3_training,
        "_evaluate_fixed_profile",
        lambda **kwargs: {"morphology_by_source": metrics},
    )

    def fake_manifest(*, out_dir, **kwargs):
        path = Path(out_dir) / "selected_checkpoint_manifest.json"
        path.write_text(json.dumps({"ok": True}) + "\n")
        return path

    monkeypatch.setattr(
        teacher_v3_training,
        "_write_selected_checkpoint_manifest",
        fake_manifest,
    )

    summary = run_teacher_v3_metadata_baseline(
        training_table=training_table,
        out_dir=tmp_path / "out",
        workers=0,
        require_cuda=False,
    )

    assert summary["native_inputs_used"] is False
    assert summary["encoder_pretraining_used"] is False
    assert len(calls) == 5
    assert all(call["native_h5"] is None for call in calls)
    assert all(call["skip_encoder_pretraining"] is True for call in calls)
    assert all(call["pretrain_epochs"] == 0 for call in calls)


def test_tic_cluster_bootstrap_is_deterministic_and_tracks_support() -> None:
    rows = pd.DataFrame(
        {
            "review_id": [f"r{index}" for index in range(9)],
            "tic": [1, 1, 2, 3, 4, 5, 6, 7, 99],
            "is_injected_row": [False] * 8 + [True],
            "fixed_split": ["test"] * 9,
        }
    )
    truth = np.asarray([0, 0, 1, 1, 2, 2, 3, 3, 0], dtype=int)
    probability = np.asarray(
        [
            [0.80, 0.10, 0.05, 0.05],
            [0.60, 0.25, 0.10, 0.05],
            [0.20, 0.65, 0.10, 0.05],
            [0.35, 0.55, 0.05, 0.05],
            [0.05, 0.10, 0.75, 0.10],
            [0.05, 0.10, 0.55, 0.30],
            [0.05, 0.05, 0.10, 0.80],
            [0.10, 0.10, 0.35, 0.45],
            [0.90, 0.05, 0.03, 0.02],
        ],
        dtype=float,
    )

    first = teacher_v3_tic_cluster_bootstrap(
        rows=rows,
        truth=truth,
        probability=probability,
        profile=TEACHER_V3_BASELINE_PROFILE,
        label_policy=TEACHER_V3_PRIMARY_LABEL_POLICY,
        draws=250,
        seed=1234,
    )
    permutation = np.asarray([8, 6, 4, 2, 0, 7, 5, 3, 1])
    second = teacher_v3_tic_cluster_bootstrap(
        rows=rows.iloc[permutation].reset_index(drop=True),
        truth=truth[permutation],
        probability=probability[permutation],
        profile=TEACHER_V3_BASELINE_PROFILE,
        label_policy=TEACHER_V3_PRIMARY_LABEL_POLICY,
        draws=250,
        seed=1234,
    )

    pd.testing.assert_frame_equal(first[0], second[0])
    pd.testing.assert_frame_equal(first[1], second[1])
    assert first[2] == second[2]
    assert first[2]["n_rows"] == 8
    assert first[2]["n_tics"] == 7
    planet = first[0].set_index("metric").loc["planet_recall"]
    assert planet["support_rows"] == 2
    assert planet["support_tics"] == 1
    assert planet["support_warning"] == "small_real_test_support"
    assert planet["draws"] == 250
    assert 0 < planet["finite_draws"] <= 250


def test_uncertain_masked_sensitivity_removes_rows_before_retraining() -> None:
    rows = pd.DataFrame(
        {
            "review_id": [
                "test_keep",
                "test_mask",
                *[f"fold_{fold}" for fold in range(5)],
                "development_mask",
            ],
            "human_label": [
                "planet_like",
                "uncertain",
                *(["eclipsing_binary_or_pceb"] * 5),
                "uncertain",
            ],
            "fixed_split": [
                "test",
                "test",
                *(["development"] * 5),
                "development",
            ],
            "cv_fold": [-1, -1, 0, 1, 2, 3, 4, 0],
            "tic": np.arange(100, 108),
        }
    )

    retained, audit = prepare_teacher_v3_uncertain_masked_rows(rows)

    assert len(rows) == 8
    assert len(retained) == 6
    assert not retained["human_label"].eq("uncertain").any()
    assert set(retained.loc[
        retained["fixed_split"].eq("development"),
        "cv_fold",
    ]) == set(range(5))
    assert audit["label_policy"] == TEACHER_V3_UNCERTAIN_MASKED_POLICY
    assert audit["operation"] == (
        "rows_removed_before_all_model_fitting_and_evaluation"
    )
    assert audit["n_masked_rows"] == 2
    assert audit["n_masked_fixed_test_rows"] == 1
    assert len(audit["masked_review_ids_sha256"]) == 64


def test_import_metadata_baseline_requires_sealed_hash_bound_run(
    tmp_path: Path,
    monkeypatch,
) -> None:
    train_config = HarmonicTrainConfig(seed=560062)
    monkeypatch.setitem(
        sys.modules,
        "torch",
        SimpleNamespace(
            load=lambda path, **kwargs: json.loads(Path(path).read_text())
        ),
    )
    training_table = tmp_path / "training.csv"
    training_table.write_text("review_id\none\n")
    provenance = build_teacher_v3_metadata_provenance(
        training_table=training_table,
    )
    source = tmp_path / "metadata_source"
    profile_dir = source / TEACHER_V3_BASELINE_PROFILE
    profile_dir.mkdir(parents=True)
    calibration_path = profile_dir / "pooled_oof_calibration.json"
    calibration_payload = {
        "schema_version": teacher_v3_training.TEACHER_V3_CALIBRATION_SCHEMA,
        "run_id": TEACHER_V3_RUN_ID,
        "profile": TEACHER_V3_BASELINE_PROFILE,
        "scope": "concatenated_five_fold_development_oof_logits",
        "n_oof_rows": 1,
        "oof_logits_sha256": "c" * 64,
        "temperature": 1.2,
        **provenance,
    }
    calibration_path.write_text(
        json.dumps(calibration_payload, sort_keys=True) + "\n"
    )
    calibration_sha256 = teacher_v3_training._file_sha256(calibration_path)
    checkpoint_records = []
    for fold in range(5):
        fold_dir = profile_dir / f"fold_{fold}"
        fold_dir.mkdir()
        checkpoint_path = fold_dir / "teacher.pt"
        checkpoint_path.write_text(
            json.dumps(
                {
                **provenance,
                "run_id": TEACHER_V3_RUN_ID,
                "release_name": TEACHER_V3_RUN_NAME,
                "model_version": MODEL_VERSION,
                "profile": TEACHER_V3_BASELINE_PROFILE,
                "fold": fold,
                "temperature_calibration_scope": (
                    "pooled_oof_development"
                ),
                "pooled_oof_calibration_sha256": calibration_sha256,
                "temperature": 1.2,
                "encoder_pretraining_skipped": True,
                "train_config": train_config.__dict__,
                },
                sort_keys=True,
            )
            + "\n"
        )
        checkpoint_records.append(
            {
                "fold": fold,
                "path": (
                    f"{TEACHER_V3_BASELINE_PROFILE}/"
                    f"fold_{fold}/teacher.pt"
                ),
                "sha256": teacher_v3_training._file_sha256(
                    checkpoint_path
                ),
                "pooled_oof_calibration_sha256": calibration_sha256,
            }
        )
    manifest_path = source / "metadata_baseline_checkpoint_manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "schema_version": TEACHER_V3_CHECKPOINT_MANIFEST_SCHEMA,
                "run_id": TEACHER_V3_RUN_ID,
                "release_name": TEACHER_V3_RUN_NAME,
                "model_version": MODEL_VERSION,
                "selected_profile": TEACHER_V3_BASELINE_PROFILE,
                "selection_policy": (
                    "fixed_metadata_baseline; "
                    "fixed_test_sealed_until_primary_freeze"
                ),
                **provenance,
                "checkpoints": checkpoint_records,
            },
            sort_keys=True,
        )
        + "\n"
    )
    freeze_path = source / "pretest_model_freeze.json"
    freeze_path.write_text(
        json.dumps(
            {
                "schema_version": "twirl_teacher_v3_model_freeze_v1",
                "run_id": TEACHER_V3_RUN_ID,
                "model_version": MODEL_VERSION,
                "profile": TEACHER_V3_BASELINE_PROFILE,
                "fixed_role": "metadata_baseline",
                "native_inputs_used": False,
                "encoder_pretraining_used": False,
                "temperature_calibration": (
                    "one_temperature_from_concatenated_five_fold_"
                    "development_oof_logits"
                ),
                "test_rows_used_for_selection_or_calibration": False,
                **provenance,
            },
            sort_keys=True,
        )
        + "\n"
    )
    summary_path = source / "summary.json"
    summary_payload = {
        "schema_version": (
            "twirl_teacher_v3_metadata_baseline_training_summary_v1"
        ),
        "run_id": TEACHER_V3_RUN_ID,
        "release_name": TEACHER_V3_RUN_NAME,
        "model_version": MODEL_VERSION,
        "profiles": [TEACHER_V3_BASELINE_PROFILE],
        "fixed_role": "metadata_baseline",
        "native_inputs_used": False,
        "encoder_pretraining_used": False,
        "train_config": train_config.__dict__,
        "fixed_test_status": "sealed_pending_primary_profile_freeze",
        "test_metrics": {},
        "calibration": {
            "profile": TEACHER_V3_BASELINE_PROFILE,
            "temperature": 1.2,
            "oof_logits_sha256": "c" * 64,
            "calibration_path": str(calibration_path),
            "calibration_sha256": calibration_sha256,
        },
        "development_fixed_profile_comparison": [
            {
                "profile": TEACHER_V3_BASELINE_PROFILE,
                "fixed_role": "metadata_baseline",
            }
        ],
        "pretest_model_freeze": str(freeze_path),
        "pretest_model_freeze_sha256": (
            teacher_v3_training._file_sha256(freeze_path)
        ),
        "metadata_baseline_checkpoint_manifest": str(manifest_path),
        "metadata_baseline_checkpoint_manifest_sha256": (
            teacher_v3_training._file_sha256(manifest_path)
        ),
        **provenance,
    }
    summary_path.write_text(
        json.dumps(summary_payload, sort_keys=True) + "\n"
    )

    calibration, observed_provenance, audit = (
        import_verified_teacher_v3_metadata_baseline(
            source_dir=source,
            out_dir=tmp_path / "combined",
            training_table=training_table,
            expected_train_config=train_config,
        )
    )

    assert calibration["temperature"] == 1.2
    assert observed_provenance == provenance
    assert audit["n_verified_checkpoints"] == 5
    assert (
        tmp_path
        / "combined"
        / TEACHER_V3_BASELINE_PROFILE
        / "fold_4"
        / "teacher.pt"
    ).is_file()

    summary_payload["test_metrics"] = {"opened": True}
    summary_path.write_text(
        json.dumps(summary_payload, sort_keys=True) + "\n"
    )
    with pytest.raises(RuntimeError, match="fixed-test metrics"):
        import_verified_teacher_v3_metadata_baseline(
            source_dir=source,
            out_dir=tmp_path / "rejected",
            training_table=training_table,
            expected_train_config=train_config,
        )

    summary_payload["test_metrics"] = {}
    summary_payload["metadata_baseline_checkpoint_manifest_sha256"] = (
        "0" * 64
    )
    summary_path.write_text(
        json.dumps(summary_payload, sort_keys=True) + "\n"
    )
    with pytest.raises(RuntimeError, match="checkpoint-manifest hash"):
        import_verified_teacher_v3_metadata_baseline(
            source_dir=source,
            out_dir=tmp_path / "bad_manifest_binding",
            training_table=training_table,
            expected_train_config=train_config,
        )

    summary_payload["metadata_baseline_checkpoint_manifest_sha256"] = (
        teacher_v3_training._file_sha256(manifest_path)
    )
    summary_path.write_text(
        json.dumps(summary_payload, sort_keys=True) + "\n"
    )
    checkpoint_zero = (
        profile_dir / "fold_0" / "teacher.pt"
    )
    checkpoint_payload = json.loads(checkpoint_zero.read_text())
    checkpoint_payload["pooled_oof_calibration_sha256"] = "0" * 64
    checkpoint_zero.write_text(
        json.dumps(checkpoint_payload, sort_keys=True) + "\n"
    )
    manifest_payload = json.loads(manifest_path.read_text())
    manifest_payload["checkpoints"][0]["sha256"] = (
        teacher_v3_training._file_sha256(checkpoint_zero)
    )
    manifest_path.write_text(
        json.dumps(manifest_payload, sort_keys=True) + "\n"
    )
    summary_payload["metadata_baseline_checkpoint_manifest_sha256"] = (
        teacher_v3_training._file_sha256(manifest_path)
    )
    summary_path.write_text(
        json.dumps(summary_payload, sort_keys=True) + "\n"
    )
    with pytest.raises(RuntimeError, match="internal calibration hash"):
        import_verified_teacher_v3_metadata_baseline(
            source_dir=source,
            out_dir=tmp_path / "bad_internal_calibration",
            training_table=training_table,
            expected_train_config=train_config,
        )

    checkpoint_payload["pooled_oof_calibration_sha256"] = calibration_sha256
    checkpoint_zero.write_text(
        json.dumps(checkpoint_payload, sort_keys=True) + "\n"
    )
    manifest_payload["checkpoints"][0]["sha256"] = (
        teacher_v3_training._file_sha256(checkpoint_zero)
    )
    manifest_path.write_text(
        json.dumps(manifest_payload, sort_keys=True) + "\n"
    )
    summary_payload["metadata_baseline_checkpoint_manifest_sha256"] = (
        teacher_v3_training._file_sha256(manifest_path)
    )
    summary_path.write_text(
        json.dumps(summary_payload, sort_keys=True) + "\n"
    )
    real_copytree = teacher_v3_training.shutil.copytree

    def copy_then_mutate_source(src, dst, *args, **kwargs):
        result = real_copytree(src, dst, *args, **kwargs)
        calibration_path.write_text(
            calibration_path.read_text() + " ",
        )
        return result

    monkeypatch.setattr(
        teacher_v3_training.shutil,
        "copytree",
        copy_then_mutate_source,
    )
    with pytest.raises(
        RuntimeError,
        match="changed while it was copied|source profile mutated",
    ):
        import_verified_teacher_v3_metadata_baseline(
            source_dir=source,
            out_dir=tmp_path / "mutated_during_copy",
            training_table=training_table,
            expected_train_config=train_config,
        )


def test_teacher_v3_output_directories_are_clean_run_only(
    tmp_path: Path,
) -> None:
    dirty_baseline = tmp_path / "dirty_baseline"
    dirty_baseline.mkdir()
    (dirty_baseline / "partial.pt").write_bytes(b"partial")
    with pytest.raises(FileExistsError, match="fresh empty output"):
        run_teacher_v3_metadata_baseline(
            training_table=tmp_path / "not_read.csv",
            out_dir=dirty_baseline,
            workers=0,
            require_cuda=False,
        )

    dirty_full = tmp_path / "dirty_full"
    dirty_full.mkdir()
    (dirty_full / "fixed_test_predictions.csv").write_text("opened\n")
    with pytest.raises(FileExistsError, match="fresh empty output"):
        teacher_v3_training.run_teacher_v3_training(
            training_table=tmp_path / "not_read.csv",
            native_registry=tmp_path / "missing_registry.csv",
            native_registry_summary=tmp_path / "missing_summary.json",
            out_dir=dirty_full,
            require_cuda=False,
        )


def test_teacher_v3_initial_table_read_detects_mutation(
    tmp_path: Path,
    monkeypatch,
) -> None:
    table = tmp_path / "training.csv"
    table.write_text("review_id\none\n")
    original_read_csv = teacher_v3_training.pd.read_csv

    def read_then_mutate(path, *args, **kwargs):
        frame = original_read_csv(path, *args, **kwargs)
        Path(path).write_text("review_id\ntwo\n")
        return frame

    monkeypatch.setattr(
        teacher_v3_training.pd,
        "read_csv",
        read_then_mutate,
    )
    with pytest.raises(RuntimeError, match="changed while it was initially read"):
        teacher_v3_training._read_training_table_with_stable_hash(
            table,
            artifact="test table",
        )


def test_teacher_v3_frozen_manifest_rejects_checkpoint_mutation(
    tmp_path: Path,
    monkeypatch,
) -> None:
    monkeypatch.setitem(
        sys.modules,
        "torch",
        SimpleNamespace(
            load=lambda path, **kwargs: json.loads(Path(path).read_text())
        ),
    )
    training_table = tmp_path / "training.csv"
    training_table.write_text("review_id\none\n")
    provenance = build_teacher_v3_metadata_provenance(
        training_table=training_table,
    )
    profile_dir = tmp_path / TEACHER_V3_BASELINE_PROFILE
    profile_dir.mkdir()
    calibration_path = profile_dir / "pooled_oof_calibration.json"
    calibration_path.write_text(
        json.dumps(
            {
                "schema_version": (
                    teacher_v3_training.TEACHER_V3_CALIBRATION_SCHEMA
                ),
                "run_id": TEACHER_V3_RUN_ID,
                "profile": TEACHER_V3_BASELINE_PROFILE,
                "scope": (
                    "concatenated_five_fold_development_oof_logits"
                ),
                "temperature": 1.25,
                **provenance,
            },
            sort_keys=True,
        )
        + "\n"
    )
    calibration_sha256 = teacher_v3_training._file_sha256(
        calibration_path
    )
    records: list[dict[str, object]] = []
    for fold in range(5):
        fold_dir = profile_dir / f"fold_{fold}"
        fold_dir.mkdir()
        checkpoint_path = fold_dir / "teacher.pt"
        checkpoint_path.write_text(
            json.dumps(
                {
                    **provenance,
                    "run_id": TEACHER_V3_RUN_ID,
                    "release_name": TEACHER_V3_RUN_NAME,
                    "model_version": MODEL_VERSION,
                    "profile": TEACHER_V3_BASELINE_PROFILE,
                    "fold": fold,
                    "temperature": 1.25,
                    "temperature_calibration_scope": (
                        "pooled_oof_development"
                    ),
                    "pooled_oof_calibration_sha256": calibration_sha256,
                },
                sort_keys=True,
            )
            + "\n"
        )
        records.append(
            {
                "fold": fold,
                "path": (
                    f"{TEACHER_V3_BASELINE_PROFILE}/"
                    f"fold_{fold}/teacher.pt"
                ),
                "sha256": teacher_v3_training._file_sha256(
                    checkpoint_path
                ),
                "pooled_oof_calibration_sha256": calibration_sha256,
            }
        )
    manifest_path = tmp_path / "metadata_baseline_checkpoint_manifest.json"
    selection_policy = (
        "fixed_metadata_baseline; fixed_test_sealed_until_primary_freeze"
    )
    manifest_path.write_text(
        json.dumps(
            {
                "schema_version": TEACHER_V3_CHECKPOINT_MANIFEST_SCHEMA,
                "run_id": TEACHER_V3_RUN_ID,
                "release_name": TEACHER_V3_RUN_NAME,
                "model_version": MODEL_VERSION,
                "selected_profile": TEACHER_V3_BASELINE_PROFILE,
                "selection_policy": selection_policy,
                **provenance,
                "checkpoints": records,
            },
            sort_keys=True,
        )
        + "\n"
    )
    manifest_sha256 = teacher_v3_training._file_sha256(manifest_path)

    audit = teacher_v3_training.verify_teacher_v3_checkpoint_manifest(
        manifest_path=manifest_path,
        expected_manifest_sha256=manifest_sha256,
        expected_profile=TEACHER_V3_BASELINE_PROFILE,
        expected_selection_policy=selection_policy,
        expected_provenance=provenance,
        calibration_path=calibration_path,
        expected_calibration_sha256=calibration_sha256,
    )
    assert audit["n_verified_checkpoints"] == 5

    checkpoint_zero = profile_dir / "fold_0" / "teacher.pt"
    checkpoint_zero.write_text(checkpoint_zero.read_text() + " ")
    with pytest.raises(RuntimeError, match="checkpoint SHA-256 mismatch"):
        teacher_v3_training.verify_teacher_v3_checkpoint_manifest(
            manifest_path=manifest_path,
            expected_manifest_sha256=manifest_sha256,
            expected_profile=TEACHER_V3_BASELINE_PROFILE,
            expected_selection_policy=selection_policy,
            expected_provenance=provenance,
            calibration_path=calibration_path,
            expected_calibration_sha256=calibration_sha256,
        )
