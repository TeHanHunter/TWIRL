from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import yaml

from twirl.models.fm0 import event_transfer_canary as canary

CONFIG = (
    Path(__file__).parents[1]
    / "configs/models/twirl_fm0_2_s66_s77_event_transfer_canary_v1.yaml"
)


def _sample(length: int = 2_048) -> dict[str, np.ndarray]:
    return {
        "flux": np.zeros((2, length), dtype=np.float32),
        "flux_valid": np.ones((2, length), dtype=bool),
        "flux_error": np.full((2, length), 0.01, dtype=np.float32),
        "error_valid": np.ones((2, length), dtype=bool),
        "local_time_cadences": np.arange(length, dtype=np.float32),
        "delta_time_cadences": np.ones(length, dtype=np.float32),
        "time_valid": np.ones(length, dtype=bool),
        "segment_boundary": np.zeros(length, dtype=bool),
        "view_present": np.ones(2, dtype=bool),
        "temporal_mask": np.zeros(length, dtype=bool),
        "reconstruction_mask": np.zeros((2, length), dtype=bool),
    }


def test_frozen_config_enforces_no_cadence_merging(tmp_path: Path) -> None:
    config, source, digest = canary.load_event_transfer_config(CONFIG)
    assert source == CONFIG.resolve()
    assert len(digest) == 64
    assert config["temporal_resolution"]["one_encoder_token_per_input_cadence"] is True
    changed = yaml.safe_load(CONFIG.read_text())
    changed["temporal_resolution"]["patching_or_cadence_averaging"] = True
    path = tmp_path / "changed.yaml"
    path.write_text(yaml.safe_dump(changed))
    with pytest.raises(ValueError, match="patching_or_cadence_averaging"):
        canary.load_event_transfer_config(path)


def test_sector_blocks_are_disjoint_and_complete() -> None:
    assert [canary.sector_block(value) for value in range(66, 72)] == ["train"] * 6
    assert [canary.sector_block(value) for value in range(72, 75)] == ["validation"] * 3
    assert [canary.sector_block(value) for value in range(75, 78)] == [
        "locked_development_test"
    ] * 3
    with pytest.raises(ValueError):
        canary.sector_block(65)


def _row(component: str, observation: str, sector: int, cohort: str) -> dict[str, str]:
    return {
        "leakage_component_id": component,
        "observation_key": observation,
        "sector": str(sector),
        "temporal_cohort": cohort,
        "source_partition": "poc_development",
        # FM0.2.1 requires ADP1x1/3x3 at six-view indices 2 and 3.  Raw and
        # ADP015 availability must not substitute for those views.
        "view_present_json": "[0,0,1,1,0,0]",
    }


def test_schedule_quarantines_cross_block_components_and_balances_cells() -> None:
    rows: list[dict[str, str]] = []
    for block_index, sector in enumerate((66, 72, 75)):
        for cohort in canary.TEMPORAL_COHORTS:
            for index in range(12):
                component = f"{cohort}-{block_index}-{index}"
                rows.append(_row(component, f"obs-{component}", sector, cohort))
    rows.extend(
        (
            _row("cross", "cross-train", 66, "new"),
            _row("cross", "cross-test", 75, "new"),
        )
    )
    loaded: list[str] = []

    def loader(row: dict[str, str]) -> dict[str, object]:
        loaded.append(row["observation_key"])
        return {"sample": _sample(), "binding": {}}

    schedule, audit = canary.freeze_component_schedule(
        rows, visit_loader=loader, target_per_cohort_block=12
    )
    assert len(schedule) == 3 * 2 * 12
    assert audit["quarantined_cross_block_components"] == 1
    assert not any(item["component_id"] == "cross" for item in schedule)
    assert "cross-train" not in loaded and "cross-test" not in loaded
    for block in canary.SPLIT_SECTORS:
        for cohort in canary.TEMPORAL_COHORTS:
            cells = [
                (item["duration_cadences"], item["fractional_depth"])
                for item in schedule
                if item["block"] == block and item["cohort"] == cohort
            ]
            assert set(cells) == set(canary.EVENT_CELLS)

    with pytest.raises(ValueError, match="exact target is 12"):
        canary.freeze_component_schedule(
            rows[:-3], visit_loader=loader, target_per_cohort_block=12
        )


def test_direct_crops_and_jittered_injection_preserve_every_cadence() -> None:
    schedule = [
        {
            "visit": {"sample": _sample()},
            "duration_cadences": 3,
            "fractional_depth": 0.1,
            "jitter_cadences": 7,
            "component_id": "component-a",
            "cohort": "new",
        }
    ]
    samples, labels, components, cohorts = canary.build_paired_samples(
        schedule, context_length=128
    )
    assert [sample["flux"].shape[-1] for sample in samples] == [128, 128]
    assert labels.tolist() == [0, 1]
    assert components == ("component-a", "component-a")
    assert cohorts == ("new", "new")
    changed = np.flatnonzero(samples[1]["flux"][0] != samples[0]["flux"][0])
    assert changed.tolist() == [70, 71, 72]
    for key in samples[0]:
        if key != "flux":
            np.testing.assert_array_equal(samples[0][key], samples[1][key])


def test_linear_max_probe_finds_single_cadence_signal_without_pooling() -> None:
    rng = np.random.default_rng(3)
    components = 48
    length = 32
    clean = rng.normal(0.0, 0.15, size=(components, length, 2))
    injected = clean.copy()
    positions = rng.integers(0, length, size=components)
    injected[np.arange(components), positions, 0] += 3.0
    tokens = np.empty((2 * components, length, 2))
    tokens[0::2] = clean
    tokens[1::2] = injected
    labels = np.tile([0, 1], components)
    valid = np.ones((2 * components, length), dtype=bool)
    fit = canary.fit_shared_linear_max_probe(tokens, valid, labels, epochs=300, seed=8)
    scores = canary.score_shared_linear_max_probe(fit, tokens, valid)
    ids = tuple(f"c{index // 2}" for index in range(2 * components))
    metrics = canary.probe_metrics(
        labels,
        scores,
        ids,
        threshold=canary.select_balanced_accuracy_threshold(labels, scores),
    )
    assert metrics["paired_component_ranking_accuracy"] > 0.95
    assert metrics["roc_auc"] > 0.95


def test_raw_and_quality_controls_keep_one_row_per_cadence() -> None:
    sample = _sample(128)
    raw, valid = canary.raw_cadence_features([sample], quality_only=False)
    quality, quality_valid = canary.raw_cadence_features([sample], quality_only=True)
    assert raw.shape == (1, 128, 12)
    assert quality.shape == (1, 128, 8)
    assert valid.shape == quality_valid.shape == (1, 128)
    np.testing.assert_array_equal(raw[0, :, 8], np.arange(128))
    np.testing.assert_array_equal(quality[0, :, 4], np.arange(128))


def test_conformer_stride_greater_than_one_is_rejected() -> None:
    output = {"h_cadence": np.zeros((2, 128, 16))}
    model = SimpleNamespace(
        config=SimpleNamespace(architecture="conformer", patch_stride=4, d_model=16)
    )
    with pytest.raises(ValueError, match="patch_stride=1"):
        canary.assert_cadence_preserving_model(model, context_length=128, output=output)


def test_paired_bootstrap_resamples_whole_pairs() -> None:
    labels = np.tile([0, 1], 6)
    scores = np.tile([0.0, 1.0], 6)
    components = tuple(f"c{index // 2}" for index in range(12))
    result = canary.paired_component_bootstrap(
        labels, scores, components, threshold=0.5, replicates=50, seed=2
    )
    assert result["paired_component_ranking_accuracy"]["estimate"] == 1.0
    assert result["roc_auc"]["paired_component_bootstrap_95_interval"] == [1.0, 1.0]


def test_validation_frozen_operating_point_and_extended_metrics() -> None:
    labels = np.tile([0, 1], 20)
    scores = np.tile([0.0, 2.0], 20)
    components = tuple(f"c{index // 2}" for index in range(40))
    fpr_threshold = canary.select_fpr_threshold(labels, scores, maximum_fpr=0.05)
    prediction = scores >= fpr_threshold
    assert np.mean(prediction[labels == 0]) <= 0.05
    metrics = canary.probe_metrics(
        labels,
        scores,
        components,
        threshold=1.0,
        fpr_threshold=fpr_threshold,
    )
    assert metrics["tpr_at_validation_frozen_5_percent_fpr"] == 1.0
    assert metrics["fpr_at_validation_frozen_5_percent_fpr"] == 0.0
    assert metrics["macro_f1"] == 1.0
    assert 0.0 <= metrics["brier_score"] <= 1.0
    assert 0.0 <= metrics["expected_calibration_error"] <= 1.0


def test_rank_metrics_handle_tied_scores_exactly() -> None:
    labels = np.asarray([0, 0, 1, 1])
    scores = np.asarray([0.0, 1.0, 1.0, 2.0])
    assert canary._rank_auc(labels, scores) == 0.875
    assert 0.0 <= canary._average_precision(labels, scores) <= 1.0


def test_readiness_uses_frozen_auc_deltas_raw_noninferiority_and_tpr() -> None:
    labels = np.tile([0, 1], 24)
    components = tuple(f"c{index // 2}" for index in range(48))
    cohorts = np.repeat(["repeated", "new"], 24)
    primary_scores = np.tile([0.0, 2.0], 24)
    random_control = np.zeros_like(primary_scores)
    scores = {
        "step2000_h_cadence_128": primary_scores,
        "step0_h_cadence_128": random_control,
        "step2000_h_cadence_2048": random_control,
        "raw_adp_validity_error_cadence_128": primary_scores,
    }
    overall = canary.paired_component_bootstrap(
        labels,
        primary_scores,
        components,
        threshold=1.0,
        fpr_threshold=1.0,
        replicates=50,
        seed=9,
    )
    by_cohort = {}
    for cohort in canary.TEMPORAL_COHORTS:
        mask = cohorts == cohort
        by_cohort[cohort] = canary.paired_component_bootstrap(
            labels[mask],
            primary_scores[mask],
            np.asarray(components, dtype=object)[mask],
            threshold=1.0,
            fpr_threshold=1.0,
            replicates=50,
            seed=10,
        )
    feature_results = {
        "step2000_h_cadence_128": {
            "locked_development_test": overall,
            "locked_development_test_by_cohort": by_cohort,
        }
    }
    config, _, _ = canary.load_event_transfer_config(CONFIG)
    readiness = canary.summarize_readiness(
        feature_results,
        scores,
        labels,
        components,
        readiness_config=config["readiness"],
        bootstrap_replicates=50,
        bootstrap_seed=11,
    )
    assert readiness["ready_for_next_real_training"] is True
    assert readiness["useful_representation_claim_supported"] is False
