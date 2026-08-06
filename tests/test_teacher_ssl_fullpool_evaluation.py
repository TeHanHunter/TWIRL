from __future__ import annotations

from pathlib import Path
import json
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import twirl.vetting.teacher_ssl_fullpool_evaluation as evaluation
from twirl.vetting.harmonic_cnn import MORPHOLOGY_CLASSES


def _development_rows() -> pd.DataFrame:
    n_rows = evaluation.EXPECTED_REAL_DEVELOPMENT_ROWS
    n_tics = evaluation.EXPECTED_REAL_DEVELOPMENT_TICS
    tic = 100_000 + np.arange(n_rows, dtype=np.int64)
    tic[n_tics:] = tic[: n_rows - n_tics]
    sector = 56 + np.arange(n_rows, dtype=np.int16) % 7
    sector[: n_rows - n_tics] = 56
    sector[n_tics:] = 57
    sector[n_tics : n_tics + 2] = sector[:2]
    period = np.ones(n_rows, dtype=float)
    epoch = np.full(n_rows, 2_459_000.0, dtype=float)
    duration = np.full(n_rows, 10.0, dtype=float)
    period[n_tics : n_tics + 2] = (2.0, 3.0)
    epoch[n_tics : n_tics + 2] += (0.1, 0.2)
    duration[n_tics : n_tics + 2] = (5.0, 6.0)
    human_label = np.asarray(
        ["uncertain"] * evaluation.EXPECTED_REAL_DEVELOPMENT_UNCERTAIN_ROWS
        + ["planet_like", "eclipsing_binary_or_pceb", "stellar_variability"]
        * 796,
        dtype=object,
    )[:n_rows]
    fold_by_tic = {
        int(value): int(value % 5) for value in np.unique(tic)
    }
    return pd.DataFrame(
        {
            "review_id": [f"r{index}" for index in range(n_rows)],
            "sector": sector,
            "tic": tic,
            "fixed_split": "development",
            "cv_fold": [fold_by_tic[int(value)] for value in tic],
            "human_label": human_label,
            "is_injected_row": False,
            "period_d": period,
            "t0_bjd": epoch,
            "duration_min": duration,
            "native_h5_path": "/old/native.h5",
            "native_group_path": "/old/group",
        }
    )


def _registry(rows: pd.DataFrame) -> pd.DataFrame:
    registry = pd.DataFrame(
        {
            "sector": rows["sector"],
            "tic": rows["tic"],
            "ssl_observation_id": [f"ssl-{value}" for value in rows["review_id"]],
            "period_d": 0.25,
            "t0_bjd": 2_460_000.0,
            "duration_min": 4.0,
            "bls_status": "ok",
            "native_h5_path": "/new/native.h5",
            "native_group_path": [f"/new/{value}" for value in rows["review_id"]],
            "native_h5_sha256": "a" * 64,
            "native_contract_version": "native-v-test",
            "ssl_pool_include": True,
            "fixed_test_member": False,
            "reserved_prospective_member": False,
            "ssl_held_out_fold": rows["cv_fold"],
            "fold_assignment_source": "frozen_development_split",
        }
    )
    return registry.drop_duplicates(["sector", "tic"], keep="first").reset_index(
        drop=True
    )


def _patch_all_labels_matched(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        evaluation,
        "EXPECTED_MATCHED_DEVELOPMENT_ROWS",
        evaluation.EXPECTED_REAL_DEVELOPMENT_ROWS,
    )
    monkeypatch.setattr(
        evaluation,
        "EXPECTED_MATCHED_DEVELOPMENT_TICS",
        evaluation.EXPECTED_REAL_DEVELOPMENT_TICS,
    )
    monkeypatch.setattr(
        evaluation,
        "EXPECTED_MATCHED_DEVELOPMENT_UNCERTAIN_ROWS",
        evaluation.EXPECTED_REAL_DEVELOPMENT_UNCERTAIN_ROWS,
    )


def test_probability_metrics_reports_multiclass_map() -> None:
    truth = np.arange(len(MORPHOLOGY_CLASSES), dtype=int)
    probability = np.eye(len(MORPHOLOGY_CLASSES), dtype=float)

    metrics = evaluation.probability_metrics(truth, probability)

    assert metrics["accuracy"] == pytest.approx(1.0)
    assert metrics["macro_f1"] == pytest.approx(1.0)
    assert metrics["macro_average_precision"] == pytest.approx(1.0)
    assert set(metrics["per_class_average_precision"]) == set(
        MORPHOLOGY_CLASSES
    )


def test_teacher_v3_baseline_binds_legacy_predictions_to_authorized_tics(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    root = tmp_path / "teacher_v3"
    summary_path = root / "summary.json"
    summary_path.parent.mkdir(parents=True)
    summary_path.write_text("{}\n", encoding="utf-8")
    monkeypatch.setattr(
        evaluation,
        "TEACHER_V3_EXPECTED_SUMMARY_SHA256",
        evaluation._file_sha256(summary_path),
    )
    review_ids = [f"review-{fold}" for fold in range(5)]
    truth = {
        review_id: fold % len(MORPHOLOGY_CLASSES)
        for fold, review_id in enumerate(review_ids)
    }
    tics = {review_id: 900_000 + fold for fold, review_id in enumerate(review_ids)}
    for fold, review_id in enumerate(review_ids):
        fold_dir = (
            root
            / evaluation.FULLPOOL_SSL_PROFILE
            / f"fold_{fold}"
        )
        fold_dir.mkdir(parents=True)
        target = truth[review_id]
        probability = np.full(len(MORPHOLOGY_CLASSES), 0.1, dtype=float)
        probability[target] = 0.7
        payload = {
            "review_id": [review_id],
            "morphology_target": [target],
            **{
                f"p_{label}": [probability[index]]
                for index, label in enumerate(MORPHOLOGY_CLASSES)
            },
        }
        pd.DataFrame(payload).to_csv(
            fold_dir / "validation_predictions.csv",
            index=False,
        )

    baseline, summary = evaluation.load_teacher_v3_baseline_oof(
        baseline_root=root,
        policy="uncertain_as_other",
        expected_review_ids=review_ids,
        expected_truth=truth,
        expected_tic=tics,
    )

    assert baseline.set_index("review_id")["tic"].to_dict() == tics
    assert summary["n_tics"] == 5
    assert summary["metrics"]["accuracy"] == pytest.approx(1.0)


def test_development_label_binding_uses_fullpool_inputs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _patch_all_labels_matched(monkeypatch)
    rows = _development_rows()
    registry = _registry(rows)
    training = tmp_path / "training.csv"
    registry_path = tmp_path / "registry.parquet"
    registry_summary = tmp_path / "registry.summary.json"
    training.write_text("training\n", encoding="utf-8")
    registry_path.write_text("registry\n", encoding="utf-8")
    registry_summary.write_text("{}\n", encoding="utf-8")
    monkeypatch.setattr(
        evaluation,
        "_read_training_table_with_stable_hash",
        lambda *args, **kwargs: (
            pd.DataFrame({"placeholder": [1]}),
            evaluation._file_sha256(training),
            {"stable_during_initial_read": True},
        ),
    )
    monkeypatch.setattr(
        evaluation,
        "validate_teacher_v3_training_table",
        lambda source: None,
    )
    monkeypatch.setattr(
        evaluation,
        "_active_release_rows",
        lambda source, seed: rows.copy(),
    )
    monkeypatch.setattr(
        evaluation,
        "load_fullpool_ssl_registry",
        lambda **kwargs: (
            registry.copy(),
            {"summary_schema_version": "test"},
        ),
    )

    bound, audit = evaluation.load_fullpool_development_labels(
        training_table_path=training,
        registry_path=registry_path,
        registry_summary_path=registry_summary,
        seed=56,
    )

    assert len(bound) == evaluation.EXPECTED_REAL_DEVELOPMENT_ROWS
    assert bound["tic"].nunique() == evaluation.EXPECTED_REAL_DEVELOPMENT_TICS
    assert np.array_equal(
        np.sort(bound["period_d"].unique()),
        np.asarray([1.0, 2.0, 3.0]),
    )
    assert bound["fullpool_period_d"].eq(0.25).all()
    assert bound["native_h5_path"].eq("/new/native.h5").all()
    assert np.array_equal(
        bound["cv_fold"].to_numpy(dtype=int),
        bound["fullpool_ssl_held_out_fold"].to_numpy(dtype=int),
    )
    assert audit["fixed_test_tensors_constructed"] is False
    assert audit["sector_63_rows_present"] is False
    assert audit["injected_rows_present"] is False
    assert audit["evaluation_unit"] == "candidate_review_record"
    assert audit["n_duplicate_target_sector_label_rows"] == 4
    assert audit["n_duplicate_target_sector_keys"] == 2


def test_development_label_binding_rejects_fold_disagreement(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _patch_all_labels_matched(monkeypatch)
    rows = _development_rows()
    registry = _registry(rows)
    registry.loc[0, "ssl_held_out_fold"] = (
        int(registry.loc[0, "ssl_held_out_fold"]) + 1
    ) % 5
    training = tmp_path / "training.csv"
    registry_path = tmp_path / "registry.parquet"
    registry_summary = tmp_path / "registry.summary.json"
    for path in (training, registry_path, registry_summary):
        path.write_text("test\n", encoding="utf-8")
    monkeypatch.setattr(
        evaluation,
        "_read_training_table_with_stable_hash",
        lambda *args, **kwargs: (
            pd.DataFrame({"placeholder": [1]}),
            evaluation._file_sha256(training),
            {},
        ),
    )
    monkeypatch.setattr(
        evaluation,
        "validate_teacher_v3_training_table",
        lambda source: None,
    )
    monkeypatch.setattr(
        evaluation,
        "_active_release_rows",
        lambda source, seed: rows.copy(),
    )
    monkeypatch.setattr(
        evaluation,
        "load_fullpool_ssl_registry",
        lambda **kwargs: (registry.copy(), {}),
    )

    with pytest.raises(RuntimeError, match="disagree"):
        evaluation.load_fullpool_development_labels(
            training_table_path=training,
            registry_path=registry_path,
            registry_summary_path=registry_summary,
            seed=56,
        )


def test_development_label_binding_proves_s63_whole_tic_exclusions(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rows = _development_rows()
    excluded_indices = [1000, 4000]
    rows.loc[excluded_indices, "sector"] = np.asarray([61, 62], dtype=np.int16)
    excluded = rows.loc[excluded_indices].copy()
    registry = _registry(rows).loc[
        lambda value: ~value.set_index(["sector", "tic"]).index.isin(
            excluded.set_index(["sector", "tic"]).index
        )
    ].reset_index(drop=True)

    reserved_path = tmp_path / "s63_reserved_tics.txt"
    reserved_path.write_text(
        "".join(f"{int(value)}\n" for value in excluded["tic"]),
        encoding="utf-8",
    )
    compact_exports = []
    for sector in (61, 62):
        tic = int(excluded.loc[excluded["sector"].eq(sector), "tic"].iloc[0])
        manifest_path = tmp_path / f"s{sector}.manifest.json"
        manifest_path.write_text(
            json.dumps({"records": [{"sector": sector, "tic": tic}]}) + "\n",
            encoding="utf-8",
        )
        compact_exports.append(
            {
                "sector": sector,
                "n_exported_targets": 1,
                "manifest": {
                    "path": str(manifest_path),
                    "sha256": evaluation._file_sha256(manifest_path),
                    "size_bytes": manifest_path.stat().st_size,
                },
            }
        )
    fullpool_summary_path = tmp_path / "fullpool.summary.json"
    fullpool_summary_path.write_text(
        json.dumps(
            {
                "counts": {
                    "input": {"n_observations": len(rows)},
                    "retained": {"n_observations": len(rows) - 2},
                    "excluded_input": {"n_observations": 2},
                },
                "exclusion_scope": (
                    "whole host: exclude every S56--S62 observation of any TIC "
                    "in the frozen fixed-test registry or prospective S63 "
                    "reservation"
                ),
                "inputs": {
                    "s63_reserved_tics": {
                        "path": str(reserved_path),
                        "sha256": evaluation._file_sha256(reserved_path),
                        "size_bytes": reserved_path.stat().st_size,
                    },
                    "compact_exports": compact_exports,
                },
            }
        )
        + "\n",
        encoding="utf-8",
    )
    registry_summary_payload = {
        "summary_schema_version": "test",
        "source_provenance": {
            "frozen_pool_summary": {
                "path": str(fullpool_summary_path),
                "sha256": evaluation._file_sha256(fullpool_summary_path),
                "size_bytes": fullpool_summary_path.stat().st_size,
            },
            "reserved_hosts": {
                "path": str(reserved_path),
                "sha256": evaluation._file_sha256(reserved_path),
                "size_bytes": reserved_path.stat().st_size,
            },
        },
    }
    training = tmp_path / "training.csv"
    registry_path = tmp_path / "registry.parquet"
    registry_summary = tmp_path / "registry.summary.json"
    training.write_text("training\n", encoding="utf-8")
    registry_path.write_text("registry\n", encoding="utf-8")
    registry_summary.write_text("{}\n", encoding="utf-8")
    monkeypatch.setattr(evaluation, "EXPECTED_S63_RESERVED_LABEL_ROWS", 2)
    monkeypatch.setattr(evaluation, "EXPECTED_S63_RESERVED_LABEL_TICS", 2)
    monkeypatch.setattr(
        evaluation,
        "EXPECTED_S63_RESERVED_LABEL_SECTORS",
        {61: 1, 62: 1},
    )
    monkeypatch.setattr(evaluation, "EXPECTED_S63_RESERVED_TICS", 2)
    monkeypatch.setattr(evaluation, "EXPECTED_FULLPOOL_RAW_ROWS", len(rows))
    monkeypatch.setattr(
        evaluation, "EXPECTED_FULLPOOL_RETAINED_ROWS", len(rows) - 2
    )
    monkeypatch.setattr(evaluation, "EXPECTED_FULLPOOL_EXCLUDED_ROWS", 2)
    monkeypatch.setattr(
        evaluation, "EXPECTED_MATCHED_DEVELOPMENT_ROWS", len(rows) - 2
    )
    monkeypatch.setattr(
        evaluation,
        "EXPECTED_MATCHED_DEVELOPMENT_TICS",
        evaluation.EXPECTED_REAL_DEVELOPMENT_TICS - 2,
    )
    monkeypatch.setattr(
        evaluation,
        "EXPECTED_MATCHED_DEVELOPMENT_UNCERTAIN_ROWS",
        evaluation.EXPECTED_REAL_DEVELOPMENT_UNCERTAIN_ROWS - 1,
    )
    monkeypatch.setattr(
        evaluation,
        "_read_training_table_with_stable_hash",
        lambda *args, **kwargs: (
            pd.DataFrame({"placeholder": [1]}),
            evaluation._file_sha256(training),
            {},
        ),
    )
    monkeypatch.setattr(
        evaluation, "validate_teacher_v3_training_table", lambda source: None
    )
    monkeypatch.setattr(
        evaluation, "_active_release_rows", lambda source, seed: rows.copy()
    )
    monkeypatch.setattr(
        evaluation,
        "load_fullpool_ssl_registry",
        lambda **kwargs: (registry.copy(), registry_summary_payload),
    )

    bound, audit = evaluation.load_fullpool_development_labels(
        training_table_path=training,
        registry_path=registry_path,
        registry_summary_path=registry_summary,
        seed=56,
    )

    assert len(bound) == len(rows) - 2
    assert not set(excluded["review_id"]).intersection(bound["review_id"])
    authority = audit["excluded_label_authority"]
    assert authority["reason"] == "prospective_s63_whole_tic_holdout"
    assert authority["n_rows"] == 2
    assert authority["raw_compact_inventory_present"] is True
    assert set(authority["raw_compact_manifests"]) == {"61", "62"}


def test_resume_recomputes_and_validates_interrupted_embedding_exports(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rows = pd.DataFrame(
        {
            "review_id": [f"review-{fold}" for fold in range(5)],
            "tic": np.arange(100, 105, dtype=np.int64),
            "cv_fold": np.arange(5, dtype=np.int8),
        }
    )
    encoders = [
        {"fold": fold, "encoder_state_dict": {"test_fold": fold}}
        for fold in range(5)
    ]

    def fake_extract(**kwargs: object) -> np.ndarray:
        state = kwargs["encoder_state"]
        assert isinstance(state, dict)
        return np.full(
            (len(rows), 128),
            float(state["test_fold"]),
            dtype=np.float32,
        )

    monkeypatch.setattr(evaluation, "_extract_encoder_embeddings", fake_extract)
    output = tmp_path / "representation"
    output.mkdir()
    reference_path = output / "reference_fold_0_embeddings.npz"
    held_path = output / "fold_safe_held_embeddings.npz"
    with reference_path.open("wb") as stream:
        np.savez_compressed(
            stream,
            review_id=rows["review_id"].to_numpy(dtype=object),
            embedding=np.zeros((len(rows), 128), dtype=np.float32),
        )
    with held_path.open("wb") as stream:
        np.savez_compressed(
            stream,
            review_id=rows["review_id"].to_numpy(dtype=object),
            embedding=np.repeat(
                np.arange(5, dtype=np.float32)[:, None], 128, axis=1
            ),
        )
    summary_path = output / "embedding_export.summary.json"
    summary_path.write_text(
        json.dumps(
            {
                "reference_embeddings": str(reference_path),
                "reference_embeddings_sha256": evaluation._file_sha256(
                    reference_path
                ),
                "fold_safe_embeddings": str(held_path),
                "fold_safe_embeddings_sha256": evaluation._file_sha256(
                    held_path
                ),
            }
        )
        + "\n",
        encoding="utf-8",
    )

    embeddings, audit = evaluation.validate_and_recompute_existing_embeddings(
        rows=rows,
        encoders=encoders,
        output_dir=output,
        batch_size=2,
        workers=0,
        seed=56,
        require_cuda=False,
    )

    assert len(embeddings) == 5
    assert audit["resume_reused"] is True
    assert "exactly_equal" in audit["resume_validation"]


def test_resume_validates_calibrated_finetuned_policy(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fake_torch = SimpleNamespace()
    monkeypatch.setitem(sys.modules, "torch", fake_torch)

    rows = pd.DataFrame(
        {
            "review_id": [f"review-{fold}" for fold in range(5)],
            "tic": np.arange(900, 905, dtype=np.int64),
            "cv_fold": np.arange(5, dtype=np.int8),
            "morphology_target_index": np.arange(5, dtype=int) % 4,
        }
    )
    fine_root = tmp_path / "fine_tuned"
    profile = fine_root / evaluation.FULLPOOL_SSL_PROFILE
    profile.mkdir(parents=True)
    probability = np.full((5, 4), 0.1, dtype=float)
    probability[np.arange(5), rows["morphology_target_index"]] = 0.7
    oof = pd.DataFrame(
        {
            "review_id": rows["review_id"],
            "tic": rows["tic"],
            "cv_fold": rows["cv_fold"],
            "morphology_target": rows["morphology_target_index"],
            **{
                f"p_{label}": probability[:, index]
                for index, label in enumerate(MORPHOLOGY_CLASSES)
            },
        }
    )
    oof_path = profile / "development_oof_predictions.csv"
    oof.to_csv(oof_path, index=False)
    embedding_path = profile / "development_oof_embeddings.npz"
    with embedding_path.open("wb") as stream:
        np.savez_compressed(
            stream,
            review_id=rows["review_id"].to_numpy(dtype=str),
            tic=rows["tic"].to_numpy(dtype=np.int64),
            embedding=np.zeros((5, 128), dtype=np.float32),
        )
    interrupted_revision = "a" * 40
    current_provenance = {
        "evaluation_code_revision": "b" * 40,
        "training_table_sha256": "1" * 64,
        "fullpool_registry_sha256": "2" * 64,
        "fullpool_registry_summary_sha256": "3" * 64,
        "fullpool_completion_release_sha256": "4" * 64,
        "fullpool_training_code_revision": "c" * 40,
    }
    calibration_path = profile / "pooled_oof_calibration.json"
    calibration = {
        "schema_version": evaluation.FULLPOOL_EVALUATION_OOF_SCHEMA,
        "scope": "five_fold_real_development_oof_fine_tuned",
        "label_policy": "uncertain_as_other",
        "n_rows": 5,
        "n_real_rows": 5,
        "n_tics": 5,
        "evaluation_code_revision": interrupted_revision,
        "oof_predictions": str(oof_path),
        "oof_predictions_sha256": evaluation._file_sha256(oof_path),
        **{
            name: value
            for name, value in current_provenance.items()
            if name != "evaluation_code_revision"
        },
    }
    calibration_path.write_text(
        json.dumps(calibration) + "\n",
        encoding="utf-8",
    )
    calibration_sha256 = evaluation._file_sha256(calibration_path)
    encoders = []
    for fold in range(5):
        fold_dir = profile / f"fold_{fold}"
        fold_dir.mkdir()
        checkpoint_path = fold_dir / "teacher.pt"
        checkpoint_path.write_text(f"checkpoint-{fold}\n", encoding="utf-8")
        encoders.append({"fold": fold, "checkpoint_sha256": f"ssl-{fold}"})

    def fake_torch_load(path: Path, **kwargs: object) -> dict[str, object]:
        fold = int(Path(path).parent.name.removeprefix("fold_"))
        return {
            "schema_version": evaluation.FULLPOOL_EVALUATION_CHECKPOINT_SCHEMA,
            "run_id": evaluation.FULLPOOL_EVALUATION_RUN_ID,
            "model_facing_name": evaluation.FULLPOOL_SSL_MODEL_FACING_NAME,
            "checkpoint_namespace": evaluation.FULLPOOL_EVALUATION_NAMESPACE,
            "label_policy": "uncertain_as_other",
            "profile": evaluation.FULLPOOL_SSL_PROFILE,
            "fold": fold,
            "ssl_encoder_checkpoint_sha256": f"ssl-{fold}",
            "temperature_calibration_scope": "pooled_development_oof",
            "pooled_oof_calibration_sha256": calibration_sha256,
            "evaluation_code_revision": interrupted_revision,
        }

    monkeypatch.setattr(fake_torch, "load", fake_torch_load, raising=False)

    loaded_oof, summary = evaluation.load_existing_finetuned_policy(
        rows=rows,
        policy="uncertain_as_other",
        fine_tune_root=fine_root,
        encoders=encoders,
        expected_input_provenance=current_provenance,
        expected_evaluation_code_revision=interrupted_revision,
    )

    assert len(loaded_oof) == 5
    assert summary["resume_reused"] is True
    assert len(summary["checkpoint_sha256"]) == 5


def test_orcd_evaluation_asset_uses_one_h200_and_sealed_development() -> None:
    root = Path(__file__).resolve().parents[1]
    asset = (
        root
        / "scripts"
        / "orcd"
        / "slurm_teacher_ssl_fullpool_v4_evaluate_h200.sbatch"
    ).read_text(encoding="utf-8")

    assert "#SBATCH --gres=gpu:h200:1" in asset
    assert "--completion-release" in asset
    assert "--baseline-teacher-v3-root" in asset
    assert "--bootstrap-draws" in asset
    assert "TWIRL_EVALUATION_RESUME" in asset
    assert "--resume" in asset
    assert "fixed_test_predictions" not in asset
    assert "sector_63" not in asset.lower()
