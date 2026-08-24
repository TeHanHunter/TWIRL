from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from twirl.plotting.fm0_development_comparison import (
    EXPECTED_LOGGED_STEPS,
    make_fm0_development_comparison,
    parse_fm0_training_log,
)


def _write_training_log(path: Path, *, loss_offset: float) -> None:
    lines = ["unrelated Slurm preamble"]
    for step in EXPECTED_LOGGED_STEPS:
        loss = loss_offset + 0.02 / (1.0 + step / 400.0)
        learning_rate = 3.0e-4 * min(1.0, step / 1_000.0)
        lines.append(
            f"[fm0-train] step={step}/20000 loss={loss:.9g} "
            f"lr={learning_rate:.9g} elapsed_s={step * 0.5:.1f}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _interval(mean: float, *, lower: float, upper: float) -> dict[str, object]:
    return {
        "n_source_clusters": 80,
        "mean": mean,
        "confidence_level": 0.95,
        "bootstrap_replicates": 1_000,
        "lower": lower,
        "upper": upper,
    }


def _encoder(*, rank: float, separation: float, retrieval: float) -> dict[str, object]:
    return {
        "embedding_health": {
            "n_embeddings": 256,
            "embedding_dim": 256,
            "per_dimension_variance_min": 0.01,
            "per_dimension_variance_median": 0.05,
            "per_dimension_variance_max": 0.2,
            "zero_or_constant_dimensions": 0,
            "near_duplicate_dimension_pairs": 0,
            "effective_rank": rank,
            "leading_principal_component_share": 0.12,
        },
        "safe_mask_pair_separation": {
            "n_pairs": 256,
            "paired_cosine_mean": 0.8,
            "unrelated_cosine_mean": 0.8 - separation,
            "paired_minus_unrelated_mean": separation,
            "paired_minus_unrelated_median": separation,
            "paired_minus_unrelated_source_clustered_95_interval": _interval(
                separation,
                lower=separation - 0.02,
                upper=separation + 0.02,
            ),
        },
        "same_component_cross_visit_retrieval": {
            "status": "available",
            "n_visit_embeddings": 300,
            "n_repeated_component_queries": 100,
            "top1_same_component_retrieval": retrieval,
            "top1_source_clustered_95_interval": _interval(
                retrieval,
                lower=retrieval - 0.05,
                upper=retrieval + 0.05,
            ),
        },
    }


def _write_representation(
    path: Path,
    *,
    variant: str,
    architecture: str,
    development_loss: float,
    component_hash: str = "components-matched",
    observation_hash: str = "observations-matched",
) -> None:
    payload = {
        "schema_version": "twirl_fm0_1_representation_health_v1",
        "passed": True,
        "evaluator_git_sha": "9cb8c96000000000000000000000000000000000",
        "run": {
            "run_dir": f"/immutable/{variant}",
            "variant": variant,
            "architecture": architecture,
            "run_git_sha": "ece8619fd72fbb7fc0d993cf6422a3f58fb880cb",
            "global_step": 20_000,
            "checkpoint_sha256": f"checkpoint-{variant}",
        },
        "evaluation_population": {
            "source_partition": "poc_development",
            "selected_leakage_components": 256,
            "selected_observation_visits": 300,
            "selected_leakage_components_sha256": component_hash,
            "selected_observation_keys_sha256": observation_hash,
        },
        "masked_reconstruction": {
            "masked_valid_target_count": 12_000,
            "masked_huber_mean": development_loss,
            "zero_prediction_masked_huber_mean": 0.01,
            "model_to_zero_baseline_ratio": development_loss / 0.01,
        },
        "trained_encoder": _encoder(
            rank=41.0 if architecture == "tcn" else 43.0,
            separation=0.16 if architecture == "tcn" else 0.18,
            retrieval=0.68 if architecture == "tcn" else 0.70,
        ),
        "random_encoder_control": _encoder(
            rank=32.0 if architecture == "tcn" else 34.0,
            separation=0.04,
            retrieval=0.32,
        ),
    }
    path.write_text(json.dumps(payload) + "\n", encoding="utf-8")


def _write_summary(
    path: Path,
    *,
    variant: str,
    architecture: str,
    parameter_count: int,
) -> None:
    path.write_text(
        json.dumps(
            {
                "variant": variant,
                "architecture": architecture,
                "global_step": 20_000,
                "parameter_count": parameter_count,
                "elapsed_seconds_this_invocation": 10_000.125,
            }
        )
        + "\n",
        encoding="utf-8",
    )


def test_fm0_development_comparison_renders_matched_report(tmp_path: Path) -> None:
    log_1 = tmp_path / "fm0_1_1.out"
    log_2 = tmp_path / "fm0_1_2.out"
    representation_1 = tmp_path / "fm0_1_1.representation.json"
    representation_2 = tmp_path / "fm0_1_2.representation.json"
    summary_1 = tmp_path / "fm0_1_1.summary.json"
    summary_2 = tmp_path / "fm0_1_2.summary.json"
    _write_training_log(log_1, loss_offset=0.00460)
    _write_training_log(log_2, loss_offset=0.00464)
    _write_representation(
        representation_1,
        variant="TWIRL-FM0.1.1",
        architecture="tcn",
        development_loss=0.0044,
    )
    _write_representation(
        representation_2,
        variant="TWIRL-FM0.1.2",
        architecture="conformer",
        development_loss=0.0043,
    )
    _write_summary(
        summary_1,
        variant="TWIRL-FM0.1.1",
        architecture="tcn",
        parameter_count=8_825_602,
    )
    _write_summary(
        summary_2,
        variant="TWIRL-FM0.1.2",
        architecture="conformer",
        parameter_count=9_345_282,
    )

    result = make_fm0_development_comparison(
        fm0_1_1_log_path=log_1,
        fm0_1_2_log_path=log_2,
        fm0_1_1_representation_path=representation_1,
        fm0_1_2_representation_path=representation_2,
        fm0_1_1_summary_path=summary_1,
        fm0_1_2_summary_path=summary_2,
        output_dir=tmp_path / "report",
    )

    for key in ("comparison_metrics_csv", "training_curve_bins_csv", "png", "pdf"):
        assert Path(result[key]).is_file()
        assert Path(result[key]).stat().st_size > 0
    metrics = pd.read_csv(result["comparison_metrics_csv"])
    assert metrics["variant"].tolist() == ["TWIRL-FM0.1.1", "TWIRL-FM0.1.2"]
    assert metrics["parameter_count"].tolist() == [8_825_602, 9_345_282]
    assert metrics["model_to_zero_baseline_ratio"].tolist() == pytest.approx(
        [0.44, 0.43]
    )
    assert metrics["training_final_1000_step_median_loss"].notna().all()
    assert metrics["training_effective_windows_per_second"].tolist() == pytest.approx(
        [128.0, 128.0]
    )
    assert metrics["trained_near_duplicate_dimension_pairs"].tolist() == [0, 0]
    assert metrics["all_current_representation_hard_gates_pass"].tolist() == [
        True,
        True,
    ]
    curves = pd.read_csv(result["training_curve_bins_csv"])
    assert len(curves) == 400
    assert set(curves["n_logged_steps"]) == {10, 11}
    provenance = json.loads(Path(result["provenance_json"]).read_text())
    binding = provenance["matched_evaluation_binding"]
    assert binding["source_partition"] == "poc_development"
    assert binding["selected_leakage_components_sha256"] == "components-matched"
    assert "not a cross-sector test" in provenance["claim_limit"]
    assert (
        provenance["training_curve_contract"][
            "minimum_training_loss_is_model_selection_metric"
        ]
        is False
    )


def test_training_log_requires_exact_real_progress_schedule(tmp_path: Path) -> None:
    path = tmp_path / "incomplete.out"
    _write_training_log(path, loss_offset=0.004)
    lines = path.read_text(encoding="utf-8").splitlines()
    path.write_text("\n".join(line for line in lines if "step=100/" not in line) + "\n")

    with pytest.raises(ValueError, match="exact step-1 plus 10-step schedule"):
        parse_fm0_training_log(path)


def test_comparison_rejects_mismatched_evaluation_population(tmp_path: Path) -> None:
    log_1 = tmp_path / "fm0_1_1.out"
    log_2 = tmp_path / "fm0_1_2.out"
    representation_1 = tmp_path / "fm0_1_1.representation.json"
    representation_2 = tmp_path / "fm0_1_2.representation.json"
    _write_training_log(log_1, loss_offset=0.00460)
    _write_training_log(log_2, loss_offset=0.00464)
    _write_representation(
        representation_1,
        variant="TWIRL-FM0.1.1",
        architecture="tcn",
        development_loss=0.0044,
    )
    _write_representation(
        representation_2,
        variant="TWIRL-FM0.1.2",
        architecture="conformer",
        development_loss=0.0043,
        observation_hash="different-observations",
    )

    with pytest.raises(
        ValueError,
        match="populations differ for selected_observation_keys_sha256",
    ):
        make_fm0_development_comparison(
            fm0_1_1_log_path=log_1,
            fm0_1_2_log_path=log_2,
            fm0_1_1_representation_path=representation_1,
            fm0_1_2_representation_path=representation_2,
            output_dir=tmp_path / "report",
        )
