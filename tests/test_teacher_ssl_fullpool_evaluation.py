from __future__ import annotations

from pathlib import Path

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


def test_development_label_binding_uses_fullpool_inputs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
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
    assert "fixed_test_predictions" not in asset
    assert "sector_63" not in asset.lower()
