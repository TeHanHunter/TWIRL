from __future__ import annotations

import hashlib
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from twirl.models.fm0 import matched_canary_evaluation as evaluation
from twirl.models.fm0 import matched_canary_payload_plan as payload_plan
from twirl.models.fm0.input_release import (
    ObservationRelease,
    deterministic_npz_bytes,
    validate_observation_release,
)
from twirl.models.fm0.matched_canary_plan import SPLIT_ORDER, TEMPORAL_COHORTS
from twirl.models.fm0.registry import FM0ContractError


def _release() -> ObservationRelease:
    n_cadences = 128
    flux = np.zeros((n_cadences, 6), dtype=np.float32)
    flux[:, 2] = np.linspace(-0.02, 0.02, n_cadences, dtype=np.float32)
    flux[:, 3] = np.linspace(0.01, -0.01, n_cadences, dtype=np.float32)
    boundary = np.zeros(n_cadences, dtype=bool)
    boundary[0] = True
    release = ObservationRelease(
        flux=flux,
        flux_valid=np.ones((n_cadences, 6), dtype=bool),
        flux_error=np.full((n_cadences, 2), 0.01, dtype=np.float32),
        error_valid=np.ones((n_cadences, 2), dtype=bool),
        local_time_cadences=100.0 + np.arange(n_cadences, dtype=np.float32),
        delta_time_cadences=np.ones(n_cadences, dtype=np.float32),
        time_valid=np.ones(n_cadences, dtype=bool),
        segment_boundary=boundary,
        segment_id=np.zeros(n_cadences, dtype=np.int64),
        view_present=np.ones(6, dtype=bool),
    )
    validate_observation_release(release)
    return release


def _metadata(index: int = 0) -> dict[str, str]:
    split = SPLIT_ORDER[(index // 480) % len(SPLIT_ORDER)]
    cohort = TEMPORAL_COHORTS[(index // 240) % len(TEMPORAL_COHORTS)]
    duration = (1, 3, 9)[index % 3]
    depth = ("0.01", "0.03", "0.1", "0.3")[index % 4]
    digest = hashlib.sha256(f"evaluation-{index}".encode()).hexdigest()
    return {
        "split": split,
        "cohort": cohort,
        "sector": "77" if split == "fresh_s77_test" else "66",
        "observation_key": f"observation_{digest}",
        "leakage_component_id": f"leakage_{digest}",
        "event_pair_id": f"event_pair_{digest}",
        "duration_cadences": str(duration),
        "fractional_depth": depth,
    }


def _pair(index: int = 0) -> evaluation.FrozenPair:
    clean = evaluation._normalize_crop_time(
        payload_plan.crop_arrays(_release(), start=0)
    )
    injected = payload_plan.inject_centered_event(
        clean,
        duration_cadences=3,
        fractional_depth=0.03,
    )
    return evaluation.FrozenPair(
        metadata=_metadata(index), clean=clean, injected=injected
    )


def _payload_row(tmp_path: Path) -> dict[str, str]:
    release = _release()
    observation = f"observation_{'5' * 64}"
    shard_dir = tmp_path / "shards"
    shard_dir.mkdir()
    shard_path = shard_dir / f"{observation}.npz"
    shard_bytes = deterministic_npz_bytes(release)
    shard_path.write_bytes(shard_bytes)
    identity = {
        "schema_version": "twirl_fm0_3_matched_canary_plan_v1",
        "event_pair_id": f"event_pair_{'1' * 64}",
        "split": "probe_train",
        "cohort": "new",
        "selection_rank_one_based": "1",
        "source_authority": "temporal_panel",
        "source_id": "later_s0066",
        "source_partition": "poc_development",
        "source_role": "development_evaluation",
        "sector": "66",
        "observation_key": observation,
        "leakage_component_id": f"leakage_{'2' * 64}",
        "duration_cadences": "3",
        "fractional_depth": "0.03",
        "event_center_index_zero_based": "64",
        "context_cadences": "128",
        "nominal_cadence_seconds": "200",
    }
    binding = payload_plan.SourceBinding(
        release_root=tmp_path,
        relative_path=f"shards/{observation}.npz",
        shard_sha256=hashlib.sha256(shard_bytes).hexdigest(),
        n_cadences=release.n_cadences,
        n_segments=release.n_segments,
        view_present=tuple(bool(value) for value in release.view_present),
    )
    return payload_plan._screen_source(
        identity,
        binding,
        identity_plan_receipt_sha256="3" * 64,
        identity_schedule_sha256="4" * 64,
        require_read_only=False,
    )


def test_bound_crop_is_reopened_hashed_and_injected_without_resampling(
    tmp_path: Path,
) -> None:
    row = _payload_row(tmp_path)
    pair = evaluation._load_one_frozen_pair(row, require_read_only=False)
    assert pair.clean["flux"].shape == (128, 6)
    assert pair.clean["local_time_cadences"][0] == 0.0
    assert not np.array_equal(pair.clean["flux"][:, 2:4], pair.injected["flux"][:, 2:4])
    assert np.array_equal(
        pair.clean["flux"][:, (0, 1, 4, 5)], pair.injected["flux"][:, (0, 1, 4, 5)]
    )

    drifted = {**row, "crop_payload_sha256": "0" * 64}
    with pytest.raises(FM0ContractError, match="crop payload binding drifted"):
        evaluation._load_one_frozen_pair(drifted, require_read_only=False)
    drifted = {**row, "source_shard_sha256": "0" * 64}
    with pytest.raises(FM0ContractError, match="source shard SHA-256 drifted"):
        evaluation._load_one_frozen_pair(drifted, require_read_only=False)


def test_exact_center_controls_keep_native_cadence_and_quality_is_paired() -> None:
    pair = _pair()
    raw = evaluation.exact_center_control_features((pair,), quality_only=False)
    quality = evaluation.exact_center_control_features((pair,), quality_only=True)
    assert raw.shape == (2, len(payload_plan.RAW_FEATURE_CHANNELS))
    assert quality.shape == (2, len(payload_plan.QUALITY_FEATURE_CHANNELS))
    assert np.array_equal(quality[0], quality[1])
    assert np.all(raw[1, :2] < raw[0, :2])
    assert np.array_equal(raw[0, 2:], raw[1, 2:])

    samples = evaluation._model_samples((pair,), variant="TWIRL-FM0.3.1")
    assert len(samples) == 2
    assert samples[0]["flux"].shape == (2, 128)
    assert not np.any(samples[0]["temporal_mask"])
    assert not np.any(samples[0]["reconstruction_mask"])


def test_fixed_probe_and_component_bootstrap_are_deterministic() -> None:
    features = np.tile(np.asarray([[-1.0], [1.0]]), (24, 1))
    labels = np.tile(np.asarray([0, 1], dtype=np.int8), 24)
    first = evaluation.fit_frozen_linear_probe(features, labels)
    second = evaluation.fit_frozen_linear_probe(features, labels)
    np.testing.assert_array_equal(first.weight, second.weight)
    np.testing.assert_array_equal(first.center, np.asarray([0.0]))
    scores = evaluation.score_frozen_linear_probe(first, features)
    components = tuple(f"component_{index // 2}" for index in range(48))
    metric = evaluation.paired_component_bootstrap(labels, scores, components)
    assert metric["sample_roc_auc"]["estimate"] == 1.0
    assert metric["sample_roc_auc"]["component_pair_bootstrap_95_interval"] == [
        1.0,
        1.0,
    ]

    chance_scores = np.repeat(np.arange(24, dtype=np.float64), 2)
    chance = evaluation.paired_component_bootstrap(labels, chance_scores, components)
    assert chance["sample_roc_auc"]["estimate"] == 0.5
    assert chance["paired_clean_injected_ranking_accuracy"]["estimate"] == 0.5
    delta = evaluation.paired_auc_delta_bootstrap(
        labels, scores, chance_scores, components
    )
    assert delta["estimate"] == 0.5
    assert delta["component_pair_bootstrap_95_interval"] == [0.5, 0.5]


def test_temporal_panel_authority_is_exact_and_recordable() -> None:
    authority = {
        "root": "/frozen/panel",
        "receipt_sha256": evaluation.TEMPORAL_PANEL_RECEIPT_SHA256,
        "panel_sha256": "1" * 64,
        "sector_bindings_sha256": "2" * 64,
    }
    receipt = {"identity_plan": {"input_authorities": {"temporal_panel": authority}}}
    assert evaluation._temporal_panel_authority(receipt) == authority
    unauthorized = {
        "identity_plan": {
            "input_authorities": {
                "temporal_panel": {**authority, "receipt_sha256": "3" * 64}
            }
        }
    }
    with pytest.raises(FM0ContractError, match="unauthorized"):
        evaluation._temporal_panel_authority(unauthorized)


def test_controls_only_bundle_is_checksum_bound_and_fail_closed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    payload_root = tmp_path / "payload"
    payload_root.mkdir()
    payload_schedule = payload_root / "schedule.csv"
    payload_schedule.write_text("fixture\n", encoding="utf-8")
    temporal = {
        "root": str(tmp_path / "panel"),
        "receipt_sha256": evaluation.TEMPORAL_PANEL_RECEIPT_SHA256,
        "panel_sha256": "1" * 64,
        "sector_bindings_sha256": "2" * 64,
    }
    fake_payload = SimpleNamespace(
        root=payload_root,
        schedule_path=payload_schedule,
        receipt_sha256="5" * 64,
        schedule_sha256="6" * 64,
        receipt={
            "producer_git_sha": "4" * 40,
            "identity_plan": {"input_authorities": {"temporal_panel": temporal}},
            "payload_bindings": {
                "source_shard_bindings_sha256": "7" * 64,
                "crop_payload_bindings_sha256": "8" * 64,
            },
        },
    )
    schedule = [_metadata(index) for index in range(1_440)]
    base = _pair()
    pairs = tuple(
        evaluation.FrozenPair(
            metadata=schedule[index],
            clean=base.clean,
            injected=base.injected,
        )
        for index in range(1_440)
    )
    monkeypatch.setattr(
        evaluation,
        "validate_matched_canary_payload_plan",
        lambda *args, **kwargs: fake_payload,
    )
    monkeypatch.setattr(evaluation, "_read_schedule", lambda path: schedule)
    monkeypatch.setattr(
        evaluation,
        "load_frozen_pairs",
        lambda rows, require_read_only: pairs,
    )
    monkeypatch.setattr(
        evaluation,
        "_metric_breakdowns",
        lambda *args, **kwargs: {"fixture_metrics": True},
    )
    monkeypatch.setattr(
        evaluation,
        "_control_gates",
        lambda results: {"passed": True, "criteria": {"fixture": True}},
    )
    output = tmp_path / "controls"
    result = evaluation.evaluate_matched_canary(
        output,
        payload_plan_root=payload_root,
        payload_plan_receipt_sha256="5" * 64,
        producer_git_sha="4" * 40,
        require_read_only=False,
    )
    assert result.receipt["feature_arms"] == list(evaluation.CONTROL_ARMS)
    assert result.receipt["n_score_rows"] == 5_760
    assert result.receipt["payload_plan"]["temporal_panel_authority"] == temporal
    repeated = evaluation.evaluate_matched_canary(
        output,
        payload_plan_root=payload_root,
        payload_plan_receipt_sha256="5" * 64,
        producer_git_sha="4" * 40,
        require_read_only=False,
    )
    assert repeated.receipt_sha256 == result.receipt_sha256

    result.scores_path.chmod(0o644)
    with result.scores_path.open("a", encoding="utf-8") as handle:
        handle.write("tampered\n")
    with pytest.raises(FM0ContractError, match="artifact hashes drifted"):
        evaluation.validate_matched_canary_evaluation(
            output,
            expected_receipt_sha256=result.receipt_sha256,
            require_read_only=False,
        )
    output.chmod(0o755)
    for path in output.iterdir():
        path.chmod(0o644)


def test_contract_has_no_search_averaging_sealed_access_or_tuning() -> None:
    contract = evaluation.evaluation_algorithm_contract()
    assert (
        contract["source_loading"]["required_temporal_panel_receipt_sha256"]
        == evaluation.TEMPORAL_PANEL_RECEIPT_SHA256
    )
    assert contract["pair_construction"]["cadence_count"] == 128
    assert (
        contract["pair_construction"][
            "cadence_averaging_merging_resampling_or_downsampling"
        ]
        is False
    )
    assert contract["feature_arms"]["model_feature"] == "h_cadence_token_64_only"
    assert contract["probe"]["epochs"] == 400
    assert contract["probe"]["validation_or_test_tuning"] is False
    assert all(contract["forbidden"].values())


def test_controls_cpu_wrapper_does_not_request_a_gpu() -> None:
    root = Path(__file__).resolve().parents[1]
    wrapper = (
        root
        / "scripts"
        / "orcd"
        / "slurm_twirl_fm0_3_matched_canary_controls_cpu.sbatch"
    ).read_text(encoding="utf-8")
    assert "--controls-only" in wrapper
    assert "--device cpu" in wrapper
    assert "#SBATCH --gres=" not in wrapper
    assert evaluation.TEMPORAL_PANEL_RECEIPT_SHA256 in wrapper
