from __future__ import annotations

import importlib.util
import hashlib
import json
from pathlib import Path
import sys

import h5py
import pandas as pd
import pytest


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts/stage5_validation/score_s56_a2v1_teacher_recovery.py"
)


def _load_script():
    spec = importlib.util.spec_from_file_location("recovery_teacher_scoring", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _inputs(tmp_path: Path) -> tuple[Path, Path, Path, tuple[Path, ...]]:
    candidates = tmp_path / "candidates.csv"
    pd.DataFrame(
        {
            "review_id": ["injection:a", "injection:b"],
            "injection_id": ["a", "b"],
            "tic": [1, 2],
            "period_d": [1.0, 2.0],
            "t0_bjd": [2459000.0, 2459001.0],
            "duration_min": [10.0, 12.0],
        }
    ).to_csv(candidates, index=False)
    native_h5 = tmp_path / "native.h5"
    with h5py.File(native_h5, "w") as h5:
        h5.attrs["training_table"] = str(candidates.resolve())
        h5.attrs["training_table_sha256"] = hashlib.sha256(
            candidates.read_bytes()
        ).hexdigest()
    model_root = tmp_path / "model"
    checkpoint_root = model_root / "shape_plus_periodogram_bls"
    checkpoints = tuple(
        checkpoint_root / f"fold_{fold}" / "teacher.pt" for fold in range(5)
    )
    for fold, path in enumerate(checkpoints):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(f"checkpoint-{fold}".encode())
    manifest = model_root / "selected_checkpoint_manifest.json"
    manifest.write_text("{}\n")
    return candidates, native_h5, manifest, checkpoints


def _patch_model_calls(module, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        module,
        "verify_deployed_teacher_checkpoint_manifest",
        lambda **_: {
            "selected_profile": "shape_plus_periodogram_bls",
            "checkpoint_namespace": "s56_harmonic_cnn_v1_native_v2",
            "input_contract_version": "s56_adp_raw_pair_v2",
        },
    )

    def score(*, candidates, checkpoint_paths, **_):
        scored = candidates.copy()
        scored["p_planet_like"] = 0.5
        return scored, {
            "checkpoint_sha256": {
                str(path): module._sha256(path) for path in checkpoint_paths
            }
        }

    monkeypatch.setattr(module, "score_harmonic_teacher_ensemble", score)


def test_recovery_scoring_publishes_manifest_bound_outputs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    module = _load_script()
    candidates, native_h5, manifest, checkpoints = _inputs(tmp_path)
    _patch_model_calls(module, monkeypatch)
    out_scores = tmp_path / "scores.parquet"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(SCRIPT),
            "--candidates",
            str(candidates),
            "--native-h5",
            str(native_h5),
            "--checkpoints",
            *(str(path) for path in checkpoints),
            "--checkpoint-manifest",
            str(manifest),
            "--out-scores",
            str(out_scores),
            "--allow-cpu",
        ],
    )

    assert module.main() == 0
    summary_path = out_scores.with_suffix(".summary.json")
    summary = json.loads(summary_path.read_text())
    assert summary["strict_provenance_passed"] is True
    assert summary["artifact_contract_version"] == (
        module.RECOVERY_SCORE_ARTIFACT_CONTRACT
    )
    assert summary["checkpoint_manifest"]["sha256"] == module._sha256(manifest)
    assert summary["out_scores_sha256"] == module._sha256(out_scores)
    assert len(pd.read_parquet(out_scores)) == 2


def test_recovery_scoring_rejects_input_changed_during_inference(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    module = _load_script()
    candidates, native_h5, manifest, checkpoints = _inputs(tmp_path)
    _patch_model_calls(module, monkeypatch)
    stable_score = module.score_harmonic_teacher_ensemble

    def mutating_score(**kwargs):
        result = stable_score(**kwargs)
        native_h5.write_bytes(b"changed-native-v2")
        return result

    monkeypatch.setattr(module, "score_harmonic_teacher_ensemble", mutating_score)
    out_scores = tmp_path / "scores.parquet"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(SCRIPT),
            "--candidates",
            str(candidates),
            "--native-h5",
            str(native_h5),
            "--checkpoints",
            *(str(path) for path in checkpoints),
            "--checkpoint-manifest",
            str(manifest),
            "--out-scores",
            str(out_scores),
            "--allow-cpu",
        ],
    )

    with pytest.raises(RuntimeError, match="inputs changed"):
        module.main()
    assert not out_scores.exists()
    assert not out_scores.with_suffix(".summary.json").exists()


def test_recovery_scoring_rejects_stale_native_candidate_binding(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    module = _load_script()
    candidates, native_h5, manifest, checkpoints = _inputs(tmp_path)
    _patch_model_calls(module, monkeypatch)
    candidates.write_text(candidates.read_text() + "\n")
    out_scores = tmp_path / "scores.parquet"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(SCRIPT),
            "--candidates",
            str(candidates),
            "--native-h5",
            str(native_h5),
            "--checkpoints",
            *(str(path) for path in checkpoints),
            "--checkpoint-manifest",
            str(manifest),
            "--out-scores",
            str(out_scores),
            "--allow-cpu",
        ],
    )

    with pytest.raises(ValueError, match="training_table_sha256"):
        module.main()
    assert not out_scores.exists()
