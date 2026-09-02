from __future__ import annotations

import csv
import hashlib
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from twirl.models.fm0 import matched_canary_payload_plan as payload_plan
from twirl.models.fm0.input_release import (
    ObservationRelease,
    deterministic_npz_bytes,
    validate_observation_release,
)
from twirl.models.fm0.matched_canary_plan import SCHEDULE_FIELDS
from twirl.models.fm0.registry import FM0ContractError


def _release(segment_lengths: tuple[int, ...] = (100, 220)) -> ObservationRelease:
    n_cadences = sum(segment_lengths)
    segment_id = np.concatenate(
        [np.full(length, index, dtype=np.int64) for index, length in enumerate(segment_lengths)]
    )
    boundary = np.zeros(n_cadences, dtype=bool)
    boundary[np.cumsum((0, *segment_lengths[:-1]))] = True
    flux = np.zeros((n_cadences, 6), dtype=np.float32)
    flux[:, 2] = np.linspace(-0.02, 0.02, n_cadences, dtype=np.float32)
    flux[:, 3] = np.linspace(0.01, -0.01, n_cadences, dtype=np.float32)
    release = ObservationRelease(
        flux=flux,
        flux_valid=np.ones((n_cadences, 6), dtype=bool),
        flux_error=np.full((n_cadences, 2), 0.01, dtype=np.float32),
        error_valid=np.ones((n_cadences, 2), dtype=bool),
        local_time_cadences=np.arange(n_cadences, dtype=np.float32),
        delta_time_cadences=np.ones(n_cadences, dtype=np.float32),
        time_valid=np.ones(n_cadences, dtype=bool),
        segment_boundary=boundary,
        segment_id=segment_id,
        view_present=np.ones(6, dtype=bool),
    )
    validate_observation_release(release)
    return release


def _identity(observation_key: str) -> dict[str, str]:
    row = {
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
        "observation_key": observation_key,
        "leakage_component_id": f"leakage_{'2' * 64}",
        "duration_cadences": "9",
        "fractional_depth": "0.03",
        "event_center_index_zero_based": "64",
        "context_cadences": "128",
        "nominal_cadence_seconds": "200",
    }
    assert tuple(row) == SCHEDULE_FIELDS
    return row


def test_direct_128_crop_is_deterministic_and_segment_coordinates_are_explicit() -> None:
    release = _release()
    kwargs = {
        "identity_plan_receipt_sha256": "3" * 64,
        "identity_schedule_sha256": "4" * 64,
        "event_pair_id": f"event_pair_{'1' * 64}",
        "observation_key": f"observation_{'5' * 64}",
    }
    first = payload_plan.choose_crop(release, **kwargs)
    second = payload_plan.choose_crop(release, **kwargs)
    moved_receipt = payload_plan.choose_crop(
        release,
        **{**kwargs, "identity_plan_receipt_sha256": "9" * 64},
    )
    assert first == second
    assert first == moved_receipt
    assert first.segment_id == 1
    assert first.segment_start == 100
    assert first.crop_stop - first.crop_start == 128
    assert first.crop_start_offset == first.crop_start - 100
    assert first.crop_stop_offset == first.crop_stop - 100
    assert first.crop_start_offset != first.crop_start
    assert first.joint_valid == 128


def test_validity_threshold_and_complete_center_support_fail_closed() -> None:
    exact = _release((128,))
    exact.flux_valid[:25, 2] = False
    assert payload_plan.eligible_crop_starts(exact) == (0,)

    below = _release((128,))
    below.flux_valid[:26, 2] = False
    assert payload_plan.eligible_crop_starts(below) == ()

    broken_support = _release((128,))
    broken_support.time_valid[64] = False
    assert payload_plan.eligible_crop_starts(broken_support) == ()
    with pytest.raises(FM0ContractError, match="no eligible"):
        payload_plan.choose_crop(
            broken_support,
            identity_plan_receipt_sha256="3" * 64,
            identity_schedule_sha256="4" * 64,
            event_pair_id=f"event_pair_{'1' * 64}",
            observation_key=f"observation_{'5' * 64}",
        )


def test_selected_shard_is_hashed_and_screened_without_padding(tmp_path: Path) -> None:
    release = _release()
    observation = f"observation_{'5' * 64}"
    shard_dir = tmp_path / "shards"
    shard_dir.mkdir()
    shard_path = shard_dir / f"{observation}.npz"
    shard_bytes = deterministic_npz_bytes(release)
    shard_path.write_bytes(shard_bytes)
    binding = payload_plan.SourceBinding(
        release_root=tmp_path,
        relative_path=f"shards/{observation}.npz",
        shard_sha256=hashlib.sha256(shard_bytes).hexdigest(),
        n_cadences=release.n_cadences,
        n_segments=release.n_segments,
        view_present=tuple(bool(value) for value in release.view_present),
    )
    row = payload_plan._screen_source(
        _identity(observation),
        binding,
        identity_plan_receipt_sha256="3" * 64,
        identity_schedule_sha256="4" * 64,
        require_read_only=False,
    )
    assert int(row["crop_stop_index_exclusive"]) - int(
        row["crop_start_index_zero_based"]
    ) == 128
    assert int(row["n_joint_valid"]) == 128
    assert row["source_shard_sha256"] == binding.shard_sha256
    assert len(row["crop_payload_sha256"]) == 64

    bad = payload_plan.SourceBinding(
        release_root=binding.release_root,
        relative_path=binding.relative_path,
        shard_sha256="0" * 64,
        n_cadences=binding.n_cadences,
        n_segments=binding.n_segments,
        view_present=binding.view_present,
    )
    with pytest.raises(FM0ContractError, match="SHA-256 drifted"):
        payload_plan._screen_source(
            _identity(observation),
            bad,
            identity_plan_receipt_sha256="3" * 64,
            identity_schedule_sha256="4" * 64,
            require_read_only=False,
        )


@pytest.mark.parametrize("duration", (1, 3, 9))
def test_event_formula_changes_only_two_adp_flux_views(duration: int) -> None:
    release = _release((128,))
    clean = payload_plan.crop_arrays(release, start=0)
    original = {name: value.copy() for name, value in clean.items()}
    injected = payload_plan.inject_centered_event(
        clean, duration_cadences=duration, fractional_depth=0.03
    )
    support_start = 64 - duration // 2
    support_stop = support_start + duration
    profile = payload_plan.centered_trapezoid(duration, 0.03)
    for view in (2, 3):
        expected = (1.0 + original["flux"][support_start:support_stop, view]) * (
            1.0 - profile
        ) - 1.0
        np.testing.assert_allclose(
            injected["flux"][support_start:support_stop, view], expected
        )
    assert np.array_equal(clean["flux"], original["flux"])
    assert np.array_equal(injected["flux"][:, (0, 1, 4, 5)], clean["flux"][:, (0, 1, 4, 5)])
    for name in payload_plan.SHARD_ARRAY_KEYS - {"flux"}:
        assert np.array_equal(injected[name], clean[name])


def test_science_contract_freezes_probe_features_metrics_and_boundaries() -> None:
    mechanics = payload_plan.science_mechanics()
    assert mechanics["bls_or_periodic_search"] is False
    assert mechanics["features"]["raw_adp_validity_error_exact_center"] == list(
        payload_plan.RAW_FEATURE_CHANNELS
    )
    assert mechanics["features"]["quality_only_exact_center"] == list(
        payload_plan.QUALITY_FEATURE_CHANNELS
    )
    assert mechanics["probe"] == {
        "family": "linear_logit_binary_classifier",
        "one_separate_fitted_weight_vector_per_feature_arm": True,
        "identical_fixed_algorithm_for_every_arm": True,
        "input": "exact_center_token_or_exact_center_raw_control_only",
        "standardization": "dimensionwise_mean_and_std_fit_on_probe_train_only",
        "fit_split": "probe_train",
        "validation_role": "diagnostic_only_no_tuning_or_selection",
        "fresh_s77_test_used_for_fit_tuning_threshold_or_selection": False,
        "optimizer": "adam",
        "epochs": 400,
        "learning_rate": 0.02,
        "l2_weight": 0.001,
        "initialization_seed": 560306,
        "hyperparameter_tuning": False,
    }
    bootstrap = mechanics["metrics"]["bootstrap"]
    assert mechanics["metrics"]["primary"] == "sample_roc_auc"
    assert bootstrap["method"] == "component_pair_clustered_bootstrap"
    assert bootstrap["replicates"] == 1000
    assert bootstrap["seed"] == 560305
    assert mechanics["inference"]["temporal_mask_all_128_cadences"] == "zero_false"


def test_freezer_publishes_and_revalidates_exact_immutable_bundle(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    identity_root = tmp_path / "identity"
    identity_root.mkdir()
    identity_schedule = identity_root / "schedule.csv"
    identity_rows: list[dict[str, str]] = []
    index = 0
    for split in payload_plan.SPLIT_ORDER:
        for cohort in payload_plan.TEMPORAL_COHORTS:
            for rank in range(1, payload_plan.TARGET_COMPONENTS_PER_COHORT + 1):
                index += 1
                label = f"{split}-{cohort}-{rank}"
                observation = f"observation_{hashlib.sha256(('o-' + label).encode()).hexdigest()}"
                component = f"leakage_{hashlib.sha256(('c-' + label).encode()).hexdigest()}"
                identity_rows.append(
                    {
                        "schema_version": "twirl_fm0_3_matched_canary_plan_v1",
                        "event_pair_id": f"event_pair_{hashlib.sha256(('p-' + label).encode()).hexdigest()}",
                        "split": split,
                        "cohort": cohort,
                        "selection_rank_one_based": str(rank),
                        "source_authority": (
                            "composite_release"
                            if split == "fresh_s77_test"
                            else "temporal_panel"
                        ),
                        "source_id": "later_s0077",
                        "source_partition": (
                            "poc_train"
                            if split == "fresh_s77_test"
                            else "poc_development"
                        ),
                        "source_role": (
                            "temporal_holdout"
                            if split == "fresh_s77_test"
                            else "development_evaluation"
                        ),
                        "sector": "77" if split == "fresh_s77_test" else "66",
                        "observation_key": observation,
                        "leakage_component_id": component,
                        "duration_cadences": ("1", "3", "9")[(rank - 1) % 3],
                        "fractional_depth": ("0.01", "0.03", "0.1", "0.3")[
                            (rank - 1) % 4
                        ],
                        "event_center_index_zero_based": "64",
                        "context_cadences": "128",
                        "nominal_cadence_seconds": "200",
                    }
                )
    with identity_schedule.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(SCHEDULE_FIELDS))
        writer.writeheader()
        writer.writerows(identity_rows)
    identity_schedule_sha = hashlib.sha256(identity_schedule.read_bytes()).hexdigest()
    fake_plan = SimpleNamespace(
        root=identity_root,
        schedule_path=identity_schedule,
        receipt_sha256="3" * 64,
        schedule_sha256=identity_schedule_sha,
        receipt={
            "producer_git_sha": "4" * 40,
            "input_authorities": {"synthetic": True},
        },
    )
    monkeypatch.setattr(
        payload_plan,
        "validate_matched_canary_evaluation_plan",
        lambda *args, **kwargs: fake_plan,
    )
    monkeypatch.setattr(
        payload_plan,
        "_resolve_sources",
        lambda rows, **kwargs: {row["observation_key"]: None for row in rows},
    )

    def fake_screen(identity: dict[str, str], binding: object, **kwargs: object):
        identity_row = {field: identity[field] for field in SCHEDULE_FIELDS}
        row = {
            "payload_schema_version": payload_plan.PAYLOAD_PLAN_SCHEMA_VERSION,
            **identity_row,
        }
        observation = row["observation_key"]
        row.update(
            {
                "identity_row_sha256": payload_plan._canonical_sha256(
                    identity_row
                ),
                "source_release_root": str(tmp_path / "release"),
                "source_relative_path": f"shards/{observation}.npz",
                "source_shard_sha256": hashlib.sha256(
                    f"source-{observation}".encode()
                ).hexdigest(),
                "source_n_cadences": "128",
                "source_n_segments": "1",
                "crop_segment_id": "0",
                "segment_start_index_zero_based": "0",
                "segment_stop_index_exclusive": "128",
                "crop_start_index_zero_based": "0",
                "crop_stop_index_exclusive": "128",
                "crop_start_offset_within_segment": "0",
                "crop_stop_offset_within_segment_exclusive": "128",
                "n_time_valid": "128",
                "n_adp_1x1_valid": "128",
                "n_adp_3x3_valid": "128",
                "n_joint_valid": "128",
                "joint_valid_fraction": "1.00000000",
                "center_support_start_index_zero_based": "60",
                "center_support_stop_index_exclusive": "69",
                "n_center_support_joint_valid": "9",
                "n_eligible_crops": "1",
                "selected_eligible_crop_rank_zero_based": "0",
                "crop_payload_sha256": hashlib.sha256(
                    f"crop-{observation}".encode()
                ).hexdigest(),
            }
        )
        return row

    monkeypatch.setattr(payload_plan, "_screen_source", fake_screen)
    output = tmp_path / "payload"
    result = payload_plan.freeze_matched_canary_payload_plan(
        output,
        identity_plan_root=identity_root,
        identity_plan_receipt_sha256="3" * 64,
        producer_git_sha="4" * 40,
        require_read_only=False,
    )
    assert result.receipt["schedule"]["n_rows"] == 1440
    assert result.receipt["identity_plan"]["receipt_sha256"] == "3" * 64
    assert result.receipt["identity_plan"]["schedule_sha256"] == identity_schedule_sha
    assert result.receipt["limits"]["bound_source_payloads_opened"] is True
    assert result.receipt["limits"]["evaluation_executed"] is False
    repeated = payload_plan.freeze_matched_canary_payload_plan(
        output,
        identity_plan_root=identity_root,
        identity_plan_receipt_sha256="3" * 64,
        producer_git_sha="4" * 40,
        require_read_only=False,
    )
    assert repeated.receipt_sha256 == result.receipt_sha256

    output.chmod(0o755)
    for path in output.iterdir():
        path.chmod(0o644)
