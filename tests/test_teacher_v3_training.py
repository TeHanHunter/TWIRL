from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
import sys

import h5py
import numpy as np
import pandas as pd

import twirl.vetting.teacher_v3_training as teacher_v3_training
from twirl.vetting.harmonic_cnn import MODEL_VERSION
from twirl.vetting.harmonic_inputs import RAW_PAIR_CONTRACT_VERSION
from twirl.vetting.teacher_v3 import (
    TEACHER_V3_CORPUS_VERSION,
    TEACHER_V3_RUN_NAME,
)
from twirl.vetting.teacher_v3_training import (
    TEACHER_V3_BASELINE_PROFILE,
    TEACHER_V3_RUN_ID,
    build_teacher_v3_metadata_provenance,
    build_teacher_v3_native_manifest,
    calibrate_teacher_v3_profile_oof,
    run_teacher_v3_metadata_baseline,
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
