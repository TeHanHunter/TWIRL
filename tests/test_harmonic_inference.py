from __future__ import annotations

import hashlib
from pathlib import Path
import sys
from types import SimpleNamespace
from typing import Any

import pandas as pd
import pytest

from twirl.vetting import harmonic_inference


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _inputs(tmp_path: Path) -> tuple[Path, Path, Path, list[Path]]:
    candidates = tmp_path / "candidates.csv"
    pd.DataFrame(
        {
            "review_id": ["real:1", "real:2"],
            "tic": [1, 2],
            "period_d": [1.0, 2.0],
            "t0_bjd": [2459825.0, 2459826.0],
            "duration_min": [10.0, 12.0],
            "native_input_include": [True, True],
        }
    ).to_csv(candidates, index=False)
    candidate_summary = tmp_path / "candidate_summary.json"
    candidate_summary.write_text("{}\n")
    native = tmp_path / "native.h5"
    native.write_bytes(b"native-v2")
    checkpoints = []
    for fold in range(5):
        path = tmp_path / f"fold_{fold}.pt"
        path.write_bytes(f"checkpoint-{fold}".encode())
        checkpoints.append(path)
    return candidates, candidate_summary, native, checkpoints


def _fake_scores(candidates: pd.DataFrame) -> pd.DataFrame:
    scored = candidates.copy()
    scored["p_planet_like"] = [0.9, 0.4]
    scored["p_preserve"] = [0.8, 0.7]
    scored["std_p_planet_like"] = [0.01, 0.02]
    scored["sde_max"] = [12.0, 9.0]
    return scored


def _checkpoint_payload(
    fold: int,
    *,
    native_h5_sha256: str = "a" * 64,
) -> dict[str, Any]:
    return {
        "model_version": harmonic_inference.MODEL_VERSION,
        "checkpoint_namespace": (
            harmonic_inference.TEACHER_NATIVE_V2_CHECKPOINT_NAMESPACE
        ),
        "input_contract_version": harmonic_inference.RAW_PAIR_CONTRACT_VERSION,
        "native_h5_sha256": native_h5_sha256,
        "training_table_sha256": "b" * 64,
        "encoder_pretraining_cache_schema": (
            harmonic_inference.ENCODER_PRETRAINING_CACHE_SCHEMA
        ),
        "encoder_pretraining_sha256": "c" * 64,
        "profile": harmonic_inference.SELECTED_TEACHER_PROFILE,
        "fold": fold,
        "model_config": {},
        "model_state_dict": {},
        "metadata_normalization": {"columns": [], "center": [], "scale": []},
    }


def test_load_ensemble_rejects_checkpoint_without_native_v2_namespace(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    checkpoint = _checkpoint_payload(0)
    checkpoint.pop("checkpoint_namespace")
    monkeypatch.setitem(
        sys.modules,
        "torch",
        SimpleNamespace(load=lambda *args, **kwargs: checkpoint),
    )

    with pytest.raises(ValueError, match="checkpoint_namespace"):
        harmonic_inference._load_ensemble(
            [tmp_path / f"fold_{fold}.pt" for fold in range(5)],
            profile=harmonic_inference.SELECTED_TEACHER_PROFILE,
            device="cpu",
        )


def test_load_ensemble_rejects_inconsistent_native_training_hashes(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    checkpoints = {
        str(tmp_path / f"fold_{fold}.pt"): _checkpoint_payload(
            fold,
            native_h5_sha256=("d" * 64 if fold == 1 else "a" * 64),
        )
        for fold in range(5)
    }

    class DummyModel:
        def to(self, device: Any) -> "DummyModel":
            return self

        def load_state_dict(self, state: dict[str, Any], strict: bool) -> None:
            return None

        def eval(self) -> None:
            return None

    monkeypatch.setitem(
        sys.modules,
        "torch",
        SimpleNamespace(
            load=lambda path, **kwargs: checkpoints[str(Path(path))],
        ),
    )
    monkeypatch.setattr(
        harmonic_inference,
        "build_harmonic_cnn",
        lambda *args, **kwargs: DummyModel(),
    )

    with pytest.raises(ValueError, match="inconsistent ensemble training provenance"):
        harmonic_inference._load_ensemble(
            [tmp_path / f"fold_{fold}.pt" for fold in range(5)],
            profile=harmonic_inference.SELECTED_TEACHER_PROFILE,
            device="cpu",
        )


def test_score_artifacts_are_hash_bound_and_cache_detects_input_change(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    candidates, candidate_summary, native, checkpoints = _inputs(tmp_path)
    monkeypatch.setattr(
        harmonic_inference,
        "verify_native_candidate_binding",
        lambda *args, **kwargs: {"passed": True, "failures": []},
    )

    def fake_ensemble(
        *,
        candidates: pd.DataFrame,
        checkpoint_paths: tuple[Path, ...],
        **kwargs: Any,
    ) -> tuple[pd.DataFrame, dict[str, Any]]:
        return _fake_scores(candidates), {
            "checkpoint_sha256": {
                str(path): _sha256(path) for path in checkpoint_paths
            },
            "n_rows": len(candidates),
        }

    monkeypatch.setattr(
        harmonic_inference, "score_harmonic_teacher_ensemble", fake_ensemble
    )
    out_dir = tmp_path / "scores"
    summary = harmonic_inference.score_harmonic_teacher_to_disk(
        candidates_path=candidates,
        candidate_summary_path=candidate_summary,
        native_h5=native,
        checkpoint_paths=checkpoints,
        out_dir=out_dir,
        require_cuda=False,
    )

    assert summary["artifact_contract_version"] == (
        harmonic_inference.TEACHER_SCORE_ARTIFACT_CONTRACT
    )
    assert summary["strict_provenance_passed"] is True
    for record in (
        summary["input_artifacts"]["candidate_table"],
        summary["input_artifacts"]["candidate_summary"],
        summary["input_artifacts"]["native_h5"],
        *summary["input_artifacts"]["checkpoints"],
    ):
        assert record["sha256_before"] == record["sha256_after"]
    for name in ("scores", "planet_enrichment_ranking"):
        assert summary["output_sha256"][name] == _sha256(
            Path(summary["outputs"][name])
        )

    verification = harmonic_inference.verify_teacher_score_cache(
        summary_path=out_dir / "summary.json",
        candidates_path=candidates,
        candidate_summary_path=candidate_summary,
        native_h5=native,
        checkpoint_paths=checkpoints,
    )
    assert verification["passed"], verification["failures"]

    candidates.write_text(candidates.read_text() + "\n")
    verification = harmonic_inference.verify_teacher_score_cache(
        summary_path=out_dir / "summary.json",
        candidates_path=candidates,
        candidate_summary_path=candidate_summary,
        native_h5=native,
        checkpoint_paths=checkpoints,
    )
    assert not verification["passed"]
    assert "candidate_table SHA256 mismatch" in verification["failures"]


def test_score_publish_rejects_native_change_during_inference(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    candidates, candidate_summary, native, checkpoints = _inputs(tmp_path)
    monkeypatch.setattr(
        harmonic_inference,
        "verify_native_candidate_binding",
        lambda *args, **kwargs: {"passed": True, "failures": []},
    )

    def mutating_ensemble(
        *,
        candidates: pd.DataFrame,
        native_h5: Path,
        checkpoint_paths: tuple[Path, ...],
        **kwargs: Any,
    ) -> tuple[pd.DataFrame, dict[str, Any]]:
        native_h5.write_bytes(native_h5.read_bytes() + b"-changed")
        return _fake_scores(candidates), {
            "checkpoint_sha256": {
                str(path): _sha256(path) for path in checkpoint_paths
            }
        }

    monkeypatch.setattr(
        harmonic_inference, "score_harmonic_teacher_ensemble", mutating_ensemble
    )
    out_dir = tmp_path / "scores"
    with pytest.raises(RuntimeError, match="inputs changed"):
        harmonic_inference.score_harmonic_teacher_to_disk(
            candidates_path=candidates,
            candidate_summary_path=candidate_summary,
            native_h5=native,
            checkpoint_paths=checkpoints,
            out_dir=out_dir,
            require_cuda=False,
        )
    assert not (out_dir / "summary.json").exists()
    assert not list(out_dir.glob("teacher_v1_*"))
