from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import twirl.vetting.teacher_ssl_training as ssl_training
from twirl.vetting.harmonic_cnn import MORPHOLOGY_CLASSES


def _prepared_rows() -> pd.DataFrame:
    records: list[dict[str, object]] = [
        {
            "review_id": "fixed_test",
            "sector": 56,
            "tic": 9000,
            "fixed_split": "test",
            "cv_fold": -1,
            "human_label": "uncertain",
            "is_injected_row": False,
        }
    ]
    for fold in range(5):
        for kind, label, injected in (
            ("uncertain", "uncertain", False),
            ("planet", "planet_like", False),
            ("injection", "planet_like", True),
        ):
            records.append(
                {
                    "review_id": f"{kind}:{fold}",
                    "sector": 56 + fold,
                    "tic": 1000 + 10 * fold + int(injected),
                    "fixed_split": "development",
                    "cv_fold": fold,
                    "human_label": label,
                    "is_injected_row": injected,
                }
            )
    return pd.DataFrame(records)


def test_teacher_ssl_output_must_be_fresh(tmp_path: Path) -> None:
    output = ssl_training.require_fresh_teacher_ssl_output_dir(
        tmp_path / "fresh"
    )
    assert output.is_dir()
    (output / "prior.txt").write_text("prior\n", encoding="utf-8")

    with pytest.raises(FileExistsError, match="fresh empty"):
        ssl_training.require_fresh_teacher_ssl_output_dir(output)


def test_json_serialization_replaces_nested_nonfinite_metrics() -> None:
    payload = {
        "python_nan": float("nan"),
        "numpy_values": np.asarray([1.0, np.inf, -np.inf], dtype=np.float64),
        "nested": {"finite": np.float64(2.5)},
    }

    encoded = ssl_training._json_bytes(payload)
    decoded = json.loads(encoded)

    assert decoded == {
        "python_nan": None,
        "numpy_values": [1.0, None, None],
        "nested": {"finite": 2.5},
    }
    json.dumps(decoded, allow_nan=False)


def test_development_guard_rejects_every_forbidden_ssl_row() -> None:
    clean = _prepared_rows().query(
        "fixed_split == 'development' and not is_injected_row and cv_fold != 2"
    )
    ssl_training._assert_development_only(
        clean,
        artifact="test",
        held_fold=2,
        real_only=True,
    )

    for mutation, message in (
        ({"fixed_split": "test"}, "fixed-test"),
        ({"cv_fold": 2}, "held development fold"),
        ({"is_injected_row": True}, "injected"),
        ({"sector": 63}, "Sector 63"),
    ):
        bad = clean.iloc[[0]].copy()
        for name, value in mutation.items():
            bad[name] = value
        with pytest.raises(RuntimeError, match=message):
            ssl_training._assert_development_only(
                bad,
                artifact="test",
                held_fold=2,
                real_only=True,
            )


def test_pilot_orchestration_never_passes_test_or_injections_to_ssl(
    tmp_path: Path,
    monkeypatch,
) -> None:
    training_table = tmp_path / "training.csv"
    training_table.write_text("frozen\n", encoding="utf-8")
    registry = tmp_path / "registry.csv"
    registry.write_text("registry\n", encoding="utf-8")
    registry_summary = tmp_path / "registry.summary.json"
    registry_summary.write_text("{}\n", encoding="utf-8")
    table_sha = ssl_training._file_sha256(training_table)
    manifest = {"schema_version": "test_manifest"}
    manifest_sha = hashlib.sha256(
        ssl_training._json_bytes(manifest)
    ).hexdigest()
    rows = _prepared_rows()
    captured_ssl: list[pd.DataFrame] = []
    captured_fine: list[pd.DataFrame] = []

    monkeypatch.setattr(
        ssl_training,
        "_read_training_table_with_stable_hash",
        lambda *args, **kwargs: (
            pd.DataFrame({"placeholder": [1]}),
            table_sha,
            {"stable_during_initial_read": True},
        ),
    )
    monkeypatch.setattr(
        ssl_training,
        "validate_teacher_v3_training_table",
        lambda source: None,
    )
    monkeypatch.setattr(
        ssl_training,
        "_active_release_rows",
        lambda source, seed: rows.copy(),
    )
    monkeypatch.setattr(
        ssl_training,
        "build_teacher_v3_native_manifest",
        lambda **kwargs: manifest,
    )
    monkeypatch.setattr(
        ssl_training,
        "build_teacher_v3_input_provenance",
        lambda **kwargs: {
            "checkpoint_namespace": "teacher_v3",
            "input_contract_version": "native",
            "native_h5_sha256": manifest_sha,
            "native_manifest_sha256": manifest_sha,
            "training_table_sha256": table_sha,
            "native_registry_sha256": "a" * 64,
            "native_registry_summary_sha256": "b" * 64,
        },
    )
    monkeypatch.setattr(
        ssl_training,
        "_split_registry_sha256",
        lambda **kwargs: ("c" * 64, "test"),
    )
    monkeypatch.setattr(
        ssl_training,
        "_code_revision",
        lambda: "d" * 40,
    )
    monkeypatch.setattr(
        ssl_training,
        "prepare_teacher_v3_uncertain_masked_rows",
        lambda frame: (
            frame.loc[~frame["human_label"].eq("uncertain")]
            .reset_index(drop=True),
            {"n_masked_rows": int(frame["human_label"].eq("uncertain").sum())},
        ),
    )

    def fake_ssl(**kwargs):
        selected, _ = ssl_training.select_ssl_fold_rows(
            kwargs["rows"],
            held_out_fold=kwargs["fold"],
        )
        captured_ssl.append(selected)
        return {
            "fold": kwargs["fold"],
            "checkpoint": f"encoder_{kwargs['fold']}.pt",
            "checkpoint_sha256": "e" * 64,
            "cache_identity": {"identity_sha256": "f" * 64},
            "encoder_state_dict": {},
        }

    def fake_fine(**kwargs):
        captured_fine.append(kwargs["rows"].copy())
        return {"fold": kwargs["fold"], "checkpoint": "fine.pt"}

    monkeypatch.setattr(ssl_training, "_run_ssl_fold", fake_ssl)
    monkeypatch.setattr(ssl_training, "_run_finetune_fold", fake_fine)
    monkeypatch.setattr(
        ssl_training,
        "_pool_development_oof",
        lambda **kwargs: (
            pd.DataFrame(
                {
                    "review_id": kwargs["rows"]["review_id"],
                    "tic": kwargs["rows"]["tic"],
                }
            ),
            {"scope": "mock_oof"},
        ),
    )
    monkeypatch.setattr(
        ssl_training,
        "_assert_teacher_v3_inputs_unchanged",
        lambda **kwargs: None,
    )

    summary = ssl_training.run_teacher_ssl_pilot(
        training_table=training_table,
        native_registry=registry,
        native_registry_summary=registry_summary,
        out_dir=tmp_path / "out",
        ssl_epochs=1,
        fine_tune_epochs=1,
        workers=0,
        require_cuda=False,
        baseline_teacher_v3_root=None,
    )

    assert len(captured_ssl) == len(captured_fine) == 5
    for frame in captured_ssl:
        assert frame["fixed_split"].eq("development").all()
        assert not frame["is_injected_row"].any()
        assert frame["human_label"].eq("uncertain").any()
    for frame in captured_fine:
        assert frame["fixed_split"].eq("development").all()
        assert not frame["human_label"].eq("uncertain").any()
        assert frame["is_injected_row"].any()
    assert summary["fixed_test_status"] == (
        "containing_files_integrity_validated_"
        "tensors_not_constructed_not_evaluated"
    )
    assert summary["fixed_test_metrics"] == {}
    assert summary["training_table_sha256"] == table_sha
    json.dumps(summary, allow_nan=False)


def _oof_fixture() -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    predictions: list[dict[str, object]] = []
    for fold in range(5):
        target = fold % len(MORPHOLOGY_CLASSES)
        probability = np.full(len(MORPHOLOGY_CLASSES), 0.05)
        probability[target] = 0.85
        rows.append(
            {
                "review_id": f"r{fold}",
                "tic": 200 + fold,
                "sector": 56 + fold,
                "fixed_split": "development",
                "cv_fold": fold,
                "human_label": MORPHOLOGY_CLASSES[target],
                "is_injected_row": False,
            }
        )
        predictions.append(
            {
                "review_id": f"r{fold}",
                "tic": 200 + fold,
                "morphology_target": target,
                **{
                    f"p_{label}": probability[index]
                    for index, label in enumerate(MORPHOLOGY_CLASSES)
                },
            }
        )
    return pd.DataFrame(rows), pd.DataFrame(predictions)


def test_teacher_v3_comparison_requires_exact_common_oof_and_summary_hash(
    tmp_path: Path,
) -> None:
    rows, current = _oof_fixture()
    root = tmp_path / "baseline"
    root.mkdir()
    summary_path = root / "summary.json"
    summary_path.write_text('{"frozen": true}\n', encoding="utf-8")
    summary_sha = ssl_training._file_sha256(summary_path)
    profile = (
        root
        / "sensitivities"
        / "uncertain_masked"
        / ssl_training.TEACHER_SSL_PROFILE
    )
    for fold in range(5):
        fold_dir = profile / f"fold_{fold}"
        fold_dir.mkdir(parents=True)
        current.loc[[fold]].drop(columns=["tic"]).to_csv(
            fold_dir / "validation_predictions.csv",
            index=False,
        )

    result = ssl_training._teacher_v3_baseline_comparison(
        current_oof=current,
        rows=rows,
        baseline_root=root,
        expected_summary_sha256=summary_sha,
        input_provenance={"training_table_sha256": "a" * 64},
        out_dir=tmp_path,
    )

    assert result["n_rows"] == 5
    assert result["baseline_summary_sha256_verified"] is True
    assert set(result["metrics"]) == {
        "teacher_v3_uncertain_masked",
        "teacher_v4_ssl",
    }

    bad_path = profile / "fold_4" / "validation_predictions.csv"
    bad = pd.read_csv(bad_path)
    bad["review_id"] = "different"
    bad.to_csv(bad_path, index=False)
    with pytest.raises(RuntimeError, match="not exactly equal"):
        ssl_training._teacher_v3_baseline_comparison(
            current_oof=current,
            rows=rows,
            baseline_root=root,
            expected_summary_sha256=summary_sha,
            input_provenance={"training_table_sha256": "a" * 64},
            out_dir=tmp_path / "bad",
        )


def test_pre_finetune_neighbor_probe_uses_held_queries(
    tmp_path: Path,
) -> None:
    rows, _ = _oof_fixture()
    # Duplicate the tiny labeled reference so top-k behavior is meaningful.
    rows = pd.concat([rows] * 3, ignore_index=True)
    rows["review_id"] = [
        f"{value}:{index}"
        for index, value in enumerate(rows["review_id"].astype(str))
    ]
    rows["tic"] = np.arange(500, 500 + len(rows))
    rng = np.random.default_rng(56)
    embedding = rng.normal(size=(len(rows), 128)).astype(np.float32)

    audit = ssl_training._write_ssl_neighbor_probe(
        rows=rows,
        embeddings=embedding,
        held_fold=2,
        path=tmp_path / "neighbors.csv",
        neighbors=3,
    )

    table = pd.read_csv(audit["path"])
    assert audit["n_query_rows"] == int(rows["cv_fold"].eq(2).sum())
    assert audit["neighbors"] == 3
    assert table["cv_fold"].eq(2).all()
