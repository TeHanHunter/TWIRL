from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

from twirl.models.fm0 import matched_canary_plan as plan
from twirl.models.fm0.composite_release import (
    HOLDOUT_ROLE,
    CompositeObservation,
)
from twirl.models.fm0.registry import FM0ContractError

ROOT = Path(__file__).resolve().parents[1]
DRIVER = (
    ROOT / "scripts/stage5_validation/"
    "freeze_twirl_fm0_3_matched_canary_evaluation_plan.py"
)
ORCD_WRAPPER = (
    ROOT / "scripts/orcd/slurm_twirl_fm0_3_matched_canary_evaluation_plan_cpu.sbatch"
)
RECEIPT_SHA = "1" * 64
SOURCE_SHA = "2" * 64
ROLE_SHA = "3" * 64
PANEL_SHA = "4" * 64
PRODUCER_SHA = "5" * 40


def _identity(prefix: str, value: str) -> str:
    return f"{prefix}_{hashlib.sha256(value.encode()).hexdigest()}"


def _panel_row(*, split: str, cohort: str, index: int) -> dict[str, str]:
    sector = 66 if split == "probe_train" else 72
    label = f"panel-{split}-{cohort}-{index}"
    return {
        "observation_key": _identity("observation", label),
        "leakage_component_id": _identity("leakage", label),
        "temporal_cohort": cohort,
        "sector": str(sector),
        "source_partition": "poc_development",
        "view_present_json": "[1,1,1,1,0,0]",
    }


def _composite_observation(
    *, label: str, component: str | None = None, views: str = "[1,1,1,1,0,0]"
) -> CompositeObservation:
    observation = _identity("observation", label)
    return CompositeObservation(
        source_id="later_s0077",
        sector=77,
        observation_key=observation,
        component=component or _identity("leakage", label),
        relative_path=f"shards/{observation}.npz",
        shard_sha256=hashlib.sha256(f"shard-{label}".encode()).hexdigest(),
        n_cadences=128,
        view_present_json=views,
    )


def _install_authorities(
    monkeypatch: pytest.MonkeyPatch,
    *,
    omit_one_validation_new: bool = False,
) -> tuple[list[dict[str, str]], set[str], set[str]]:
    panel_rows: list[dict[str, str]] = []
    for split in ("probe_train", "validation"):
        for cohort in plan.TEMPORAL_COHORTS:
            count = plan.TARGET_COMPONENTS_PER_COHORT
            if split == "validation" and cohort == "new" and omit_one_validation_new:
                count -= 1
            panel_rows.extend(
                _panel_row(split=split, cohort=cohort, index=index)
                for index in range(count)
            )

    cross_component = _identity("leakage", "cross-block")
    for split, sector in (("probe_train", 67), ("validation", 73)):
        panel_rows.append(
            {
                **_panel_row(split=split, cohort="repeated", index=99_999),
                "observation_key": _identity("observation", f"cross-{split}"),
                "leakage_component_id": cross_component,
                "sector": str(sector),
            }
        )
    old_panel_component = _identity("leakage", "old-panel-s77")
    panel_rows.append(
        {
            **_panel_row(split="validation", cohort="new", index=88_888),
            "observation_key": _identity("observation", "old-panel-s77"),
            "leakage_component_id": old_panel_component,
            "sector": "77",
        }
    )

    def fake_load_temporal_panel(root: str | Path, *, receipt_sha256: str):
        assert receipt_sha256 == PANEL_SHA
        return panel_rows, {
            "root": str(Path(root).resolve()),
            "receipt_sha256": PANEL_SHA,
            "panel_sha256": "6" * 64,
            "sector_bindings_sha256": "7" * 64,
        }

    repeated: list[CompositeObservation] = []
    new: list[CompositeObservation] = []
    for cohort, target in (
        ("repeated", repeated),
        ("new", new),
    ):
        target.extend(
            _composite_observation(label=f"test-{cohort}-{index}")
            for index in range(plan.TARGET_COMPONENTS_PER_COHORT)
        )
    overlap = _composite_observation(
        label="test-overlaps-panel", component=old_panel_component
    )
    missing = _composite_observation(label="test-missing-view", views="[1,1,1,0,0,0]")
    holdout = [*repeated, *new, overlap, missing]
    excluded = [
        _composite_observation(label=f"earlier-{index}", component=item.component)
        for index, item in enumerate(repeated)
    ]
    excluded.append(
        _composite_observation(label="earlier-overlap", component=overlap.component)
    )

    def fake_load_composite_identity_authority(
        root: str | Path,
        *,
        receipt_sha256: str,
        source_bindings_sha256: str,
        role_index_sha256: str,
        require_read_only: bool,
    ):
        assert receipt_sha256 == RECEIPT_SHA
        assert source_bindings_sha256 == SOURCE_SHA
        assert role_index_sha256 == ROLE_SHA
        assert require_read_only is True
        return plan._CompositeIdentityAuthority(
            receipt={"producer_git_sha": "8" * 40},
            holdout=tuple(holdout),
            repeated_components=frozenset(item.component for item in excluded),
        )

    monkeypatch.setattr(plan, "load_temporal_panel", fake_load_temporal_panel)
    monkeypatch.setattr(
        plan,
        "_load_composite_identity_authority",
        fake_load_composite_identity_authority,
    )
    return (
        panel_rows,
        {item.component for item in repeated},
        {item.component for item in new},
    )


def _freeze_args(tmp_path: Path, output_name: str = "plan") -> dict[str, object]:
    return {
        "output_dir": tmp_path / output_name,
        "temporal_panel_root": tmp_path / "panel-authority",
        "temporal_panel_receipt_sha256": PANEL_SHA,
        "composite_root": tmp_path / "composite-authority",
        "composite_receipt_sha256": RECEIPT_SHA,
        "composite_source_bindings_sha256": SOURCE_SHA,
        "composite_role_index_sha256": ROLE_SHA,
        "producer_git_sha": PRODUCER_SHA,
    }


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def _write_csv(path: Path, fields: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def test_composite_identity_loader_never_requires_a_shard_directory(
    tmp_path: Path,
) -> None:
    composite_root = tmp_path / "composite"
    s77_root = tmp_path / "s77-release"
    repeated_component = _identity("leakage", "loader-repeated")
    new_component = _identity("leakage", "loader-new")
    repeated_observation = _identity("observation", "loader-repeated")
    new_observation = _identity("observation", "loader-new")

    def manifest_row(observation: str, component: str) -> dict[str, str]:
        values = {field: "fixture" for field in plan.LATER_SIX_VIEW_MANIFEST_FIELDS}
        values.update(
            {
                "manifest_schema_version": (
                    plan.LATER_SIX_VIEW_MANIFEST_SCHEMA_VERSION
                ),
                "observation_key": observation,
                "sector": "77",
                "leakage_component_id": component,
                "source_partition": "poc_train",
                "relative_path": f"shards/{observation}.npz",
                "sha256": hashlib.sha256(f"shard-{observation}".encode()).hexdigest(),
                "n_cadences": "128",
                "view_present_json": "[1,1,1,1,0,0]",
            }
        )
        return values

    manifest = s77_root / "manifest.csv"
    _write_csv(
        manifest,
        plan.LATER_SIX_VIEW_MANIFEST_FIELDS,
        [
            manifest_row(repeated_observation, repeated_component),
            manifest_row(new_observation, new_component),
        ],
    )
    manifest_sha = hashlib.sha256(manifest.read_bytes()).hexdigest()
    binding_rows: list[dict[str, str]] = []
    source_ids = [
        "legacy_s56_s64",
        *(f"later_s{sector:04d}" for sector in range(66, 78)),
    ]
    for source_id in source_ids:
        is_s77 = source_id == "later_s0077"
        sector = "77" if is_s77 else "66"
        binding_rows.append(
            {
                "schema_version": plan.SOURCE_BINDING_SCHEMA_VERSION,
                "source_id": source_id,
                "source_kind": "later" if source_id != "legacy_s56_s64" else "legacy",
                "sector_min": "56" if source_id == "legacy_s56_s64" else sector,
                "sector_max": "64" if source_id == "legacy_s56_s64" else sector,
                "release_root": str(s77_root if is_s77 else tmp_path / source_id),
                "manifest_sha256": manifest_sha if is_s77 else "9" * 64,
                "receipt_path": str(tmp_path / source_id / "receipt.json"),
                "receipt_sha256": "a" * 64,
                "n_rows": "2" if is_s77 else "1",
                "n_cadences": "256" if is_s77 else "1",
                "rows_by_partition_json": (
                    '{"poc_development":0,"poc_sealed_test":0,"poc_train":2}'
                    if is_s77
                    else '{"poc_development":0,"poc_sealed_test":0,"poc_train":1}'
                ),
            }
        )
    source_bindings = composite_root / "source_bindings.csv"
    _write_csv(source_bindings, plan.SOURCE_BINDING_FIELDS, binding_rows)
    source_sha = hashlib.sha256(source_bindings.read_bytes()).hexdigest()

    role_rows = [
        {
            "schema_version": plan.ROLE_INDEX_SCHEMA_VERSION,
            "source_id": "later_s0077",
            "observation_key": repeated_observation,
            "leakage_component_id": repeated_component,
            "role": HOLDOUT_ROLE,
        },
        {
            "schema_version": plan.ROLE_INDEX_SCHEMA_VERSION,
            "source_id": "later_s0077",
            "observation_key": new_observation,
            "leakage_component_id": new_component,
            "role": HOLDOUT_ROLE,
        },
        {
            "schema_version": plan.ROLE_INDEX_SCHEMA_VERSION,
            "source_id": "later_s0066",
            "observation_key": _identity("observation", "loader-earlier"),
            "leakage_component_id": repeated_component,
            "role": plan.EXCLUDED_OVERLAP_ROLE,
        },
    ]
    role_index = composite_root / "role_index.csv"
    _write_csv(role_index, plan.ROLE_INDEX_FIELDS, role_rows)
    role_sha = hashlib.sha256(role_index.read_bytes()).hexdigest()
    receipt = {
        "schema_version": plan.COMPOSITE_RELEASE_SCHEMA_VERSION,
        "release_state": plan.COMPOSITE_RELEASE_STATE,
        "passed": True,
        "producer_git_sha": "b" * 40,
        "sources": {"source_bindings_sha256": source_sha},
        "selection": {
            "temporal_holdout_sector": 77,
            "role_index_sha256": role_sha,
            "role_counts": {
                plan.ROLE_ORDER[0]: 0,
                HOLDOUT_ROLE: 2,
                plan.EXCLUDED_OVERLAP_ROLE: 1,
            },
        },
        "limits": {
            "identity_only": True,
            "source_shards_opened": False,
            "sealed_rows_selected": 0,
        },
    }
    receipt_path = composite_root / "receipt.json"
    receipt_path.write_text(json.dumps(receipt, sort_keys=True) + "\n")
    receipt_sha = hashlib.sha256(receipt_path.read_bytes()).hexdigest()
    ready = composite_root / "READY"
    ready.write_text(receipt_sha + "\n")
    for path in (manifest, source_bindings, role_index, receipt_path, ready):
        path.chmod(0o444)
    s77_root.chmod(0o555)
    composite_root.chmod(0o555)

    authority = plan._load_composite_identity_authority(
        composite_root,
        receipt_sha256=receipt_sha,
        source_bindings_sha256=source_sha,
        role_index_sha256=role_sha,
        require_read_only=True,
    )
    assert len(authority.holdout) == 2
    assert authority.repeated_components == {repeated_component}
    assert not (s77_root / "shards").exists()
    assert all(not (tmp_path / source_id).exists() for source_id in source_ids[:-1])


def test_freeze_builds_one_balanced_shared_component_disjoint_schedule(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    panel_rows, repeated_test, new_test = _install_authorities(monkeypatch)
    result = plan.freeze_matched_canary_evaluation_plan(**_freeze_args(tmp_path))
    rows = _rows(result.schedule_path)

    expected_rows = (
        len(plan.SPLIT_ORDER)
        * len(plan.TEMPORAL_COHORTS)
        * plan.TARGET_COMPONENTS_PER_COHORT
    )
    assert len(rows) == expected_rows == 1_440
    assert {row["split"] for row in rows} == set(plan.SPLIT_ORDER)
    assert {row["cohort"] for row in rows} == set(plan.TEMPORAL_COHORTS)
    assert {row["event_center_index_zero_based"] for row in rows} == {"64"}
    assert {row["context_cadences"] for row in rows} == {"128"}
    assert {row["nominal_cadence_seconds"] for row in rows} == {"200"}
    assert {row["duration_cadences"] for row in rows} == {"1", "3", "9"}
    assert {row["fractional_depth"] for row in rows} == {
        "0.01",
        "0.03",
        "0.1",
        "0.3",
    }

    split_components = {
        split: {row["leakage_component_id"] for row in rows if row["split"] == split}
        for split in plan.SPLIT_ORDER
    }
    assert not (split_components["probe_train"] & split_components["validation"])
    assert not (split_components["probe_train"] & split_components["fresh_s77_test"])
    assert not (split_components["validation"] & split_components["fresh_s77_test"])
    panel_components = {row["leakage_component_id"] for row in panel_rows}
    assert not (split_components["fresh_s77_test"] & panel_components)

    test_rows = [row for row in rows if row["split"] == "fresh_s77_test"]
    assert {row["sector"] for row in test_rows} == {"77"}
    assert {row["source_role"] for row in test_rows} == {HOLDOUT_ROLE}
    assert {row["source_partition"] for row in test_rows} == {"poc_train"}
    assert {
        row["leakage_component_id"] for row in test_rows if row["cohort"] == "repeated"
    } == repeated_test
    assert {
        row["leakage_component_id"] for row in test_rows if row["cohort"] == "new"
    } == new_test
    assert not any(
        row["sector"] in {"75", "76", "77"}
        for row in rows
        if row["source_authority"] == "temporal_panel"
    )

    for split in plan.SPLIT_ORDER:
        for cohort in plan.TEMPORAL_COHORTS:
            subset = [
                row for row in rows if row["split"] == split and row["cohort"] == cohort
            ]
            cells = {
                (row["duration_cadences"], row["fractional_depth"]) for row in subset
            }
            assert len(cells) == 12
            assert {
                sum(
                    (row["duration_cadences"], row["fractional_depth"]) == cell
                    for row in subset
                )
                for cell in cells
            } == {20}

    receipt = result.receipt
    assert receipt["architecture_comparison"]["variants"] == list(plan.ARCHITECTURES)
    assert (
        receipt["architecture_comparison"]["shared_schedule_sha256"]
        == result.schedule_sha256
    )
    assert (
        receipt["architecture_comparison"]["identical_schedule_across_architectures"]
        is True
    )
    assert "architecture" not in plan.SCHEDULE_FIELDS
    assert receipt["limits"] == {
        "identity_only": True,
        "s77_shard_payloads_opened": False,
        "any_shard_payloads_opened": False,
        "events_injected": False,
        "checkpoints_loaded": False,
        "probe_fitted": False,
        "metrics_computed": False,
        "sealed_test_opened": False,
        "model_training_authorized": False,
    }
    analysis = receipt["analysis_contract"]
    assert analysis["probe"] == {
        "family": "full_batch_adam_linear_logit",
        "fit_split": "probe_train",
        "fit_sectors": [66, 67, 68, 69, 70, 71],
        "validation_split": "validation",
        "validation_sectors": [72, 73, 74],
        "validation_role": "diagnostic_only_no_tuning_or_selection",
        "test_split": "fresh_s77_test",
        "test_sector": 77,
        "standardization_fit_on_probe_train_only": True,
        "test_used_for_fit_tuning_or_threshold": False,
        "one_separate_fitted_weight_vector_per_feature_arm": True,
        "optimizer": "adam",
        "epochs": 400,
        "learning_rate": 0.02,
        "l2_weight": 0.001,
        "initialization_seed": 560306,
        "hyperparameter_tuning": False,
    }
    assert analysis["center_readout"] == {
        "candidate_center_index_zero_based": 64,
        "encoder_tokens_required_before_readout": 128,
        "encoder_token_stride": 1,
        "score_definition": "linear_logit_from_h_cadence_token_64_only",
        "all_128_tokens_exist_before_center_token_selection": True,
        "token_or_window_pooling_before_score": False,
        "off_center_tokens_used_as_probe_inputs": False,
    }
    assert analysis["arms"] == {
        "controls": [
            "raw_adp_validity_error_exact_center",
            "quality_only_exact_center",
        ],
        "TWIRL-FM0.3.1": [
            "step0_h_cadence_token64",
            "step2000_h_cadence_token64",
        ],
        "TWIRL-FM0.3.2": [
            "step0_h_cadence_token64",
            "step2000_h_cadence_token64",
        ],
        "same_schedule_fit_and_metric_code_for_every_arm": True,
    }
    assert not any(
        analysis["feature_boundaries"][key]
        for key in (
            "synthetic_support_used_as_input",
            "event_duration_used_as_input",
            "event_depth_used_as_input",
            "event_center_numeric_value_used_as_input",
            "bls_or_period_features_used_as_input",
            "search_or_candidate_score_used_as_input",
        )
    )
    assert analysis["metrics"]["primary"] == "sample_roc_auc"
    assert analysis["metrics"]["diagnostic"] == "paired_clean_injected_ranking_accuracy"
    assert analysis["metrics"]["confidence_interval"] == {
        "level": 0.95,
        "method": "deterministic_cluster_bootstrap",
        "replicates": 1000,
        "seed": 560305,
        "resampling_unit": "whole_leakage_component_with_clean_injected_pair",
        "same_resample_indices_across_arms_and_paired_deltas": True,
    }
    gates = analysis["primary_test_gates"]
    assert gates["raw_control"]["minimum_overall_sample_roc_auc_lower_95"] == 0.90
    assert gates["quality_only_control"] == {
        "maximum_absolute_overall_sample_roc_auc_from_chance": 0.03,
        "overall_sample_roc_auc_interval_must_contain": 0.50,
    }
    fm_gates = gates["each_architecture_step2000"]
    assert fm_gates["minimum_overall_sample_roc_auc_lower_95"] == 0.75
    assert fm_gates["minimum_each_cohort_sample_roc_auc_estimate"] == 0.70
    assert (
        fm_gates["minimum_each_cohort_sample_roc_auc_lower_95_strictly_above"] == 0.50
    )
    assert fm_gates["paired_step2000_minus_own_step0_auc"] == {
        "minimum_estimate": 0.02,
        "lower_95_strictly_above": 0.0,
    }
    assert fm_gates["blocking_depth_aggregates"] == {
        "fractional_depths": [0.1, 0.3],
        "minimum_each_depth_sample_roc_auc_lower_95": 0.80,
        "aggregate_over_durations": [1, 3, 9],
    }
    assert fm_gates["nonblocking_depth_aggregates"] == {
        "fractional_depths": [0.01, 0.03],
        "report_with_intervals": True,
    }
    assert analysis["decision_rule"] == {
        "one_seed": 560067,
        "raw_and_quality_control_gates_must_pass_before_fm_interpretation": True,
        "one_seed_pass_authorizes_only": (
            "longer_matched_two_seed_architecture_comparison"
        ),
        "one_seed_pass_authorizes_architecture_promotion": False,
        "one_seed_pass_authorizes_foundation_model_claim": False,
        "architecture_comparison_interpretable_only_if_both_step2000_arms_"
        "pass_all_matched_gates": True,
        "otherwise_architecture_result": "inconclusive",
    }
    assert result.ready_path.read_text().strip() == result.receipt_sha256
    assert {path.name for path in result.root.iterdir()} == {
        "schedule.csv",
        "receipt.json",
        "READY",
    }
    assert not result.root.stat().st_mode & 0o222
    assert all(not path.stat().st_mode & 0o222 for path in result.root.iterdir())


def test_freeze_is_byte_deterministic_and_idempotent(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _install_authorities(monkeypatch)
    first = plan.freeze_matched_canary_evaluation_plan(**_freeze_args(tmp_path, "a"))
    second = plan.freeze_matched_canary_evaluation_plan(**_freeze_args(tmp_path, "b"))
    again = plan.freeze_matched_canary_evaluation_plan(**_freeze_args(tmp_path, "a"))
    assert first.schedule_path.read_bytes() == second.schedule_path.read_bytes()
    assert first.receipt_path.read_bytes() == second.receipt_path.read_bytes()
    assert first.receipt_sha256 == second.receipt_sha256 == again.receipt_sha256


def test_freeze_fails_closed_on_split_shortfall(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _install_authorities(monkeypatch, omit_one_validation_new=True)
    args = _freeze_args(tmp_path)
    with pytest.raises(
        FM0ContractError, match="validation/new has only 239 eligible components"
    ):
        plan.freeze_matched_canary_evaluation_plan(**args)
    assert not Path(args["output_dir"]).exists()


def test_driver_exposes_only_identity_authorities_and_output() -> None:
    spec = importlib.util.spec_from_file_location("matched_canary_plan_driver", DRIVER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    destinations = {action.dest for action in module.build_parser()._actions}
    assert {
        "temporal_panel_root",
        "temporal_panel_receipt_sha256",
        "composite_root",
        "composite_receipt_sha256",
        "composite_source_bindings_sha256",
        "composite_role_index_sha256",
        "producer_git_sha",
        "output",
    } <= destinations
    assert (
        not {
            "checkpoint",
            "architecture",
            "shard",
            "bls",
            "period",
        }
        & destinations
    )


def test_fixed_contract_preserves_native_cadences_and_balanced_cells() -> None:
    assert plan.TARGET_COMPONENTS_PER_COHORT == 240
    assert plan.CONTEXT_CADENCES == 128
    assert plan.NOMINAL_CADENCE_SECONDS == 200
    assert plan.EVENT_CENTER_INDEX_ZERO_BASED == 64
    assert plan.EVENT_DURATIONS == (1, 3, 9)
    assert plan.EVENT_DEPTHS_TEXT == ("0.01", "0.03", "0.1", "0.3")
    assert len(plan.EVENT_CELLS) == 12
    assert plan.TARGET_COMPONENTS_PER_COHORT % len(plan.EVENT_CELLS) == 0
    assert "relative_path" not in plan.SCHEDULE_FIELDS
    assert "shard_sha256" not in plan.SCHEDULE_FIELDS


def test_cpu_wrapper_binds_authorities_and_only_freezes_the_plan() -> None:
    script = ORCD_WRAPPER.read_text(encoding="utf-8")
    assert "#SBATCH -p pg_mki_aryeh" in script
    assert "#SBATCH --exclude=node4900" in script
    assert "#SBATCH --gres" not in script
    assert "nvidia-smi" not in script
    assert "78c370e10c556472c5997c20cfe95207a0b334bafe7f024bf7ba4fc7ec4de624" in script
    assert "cc5cc3bce4c24e74bef1fbf084f407855233de7893183e6bc3486e284a2f44d9" in script
    assert "8cbcac99409ab89fe2dee0c36687f50122e78376206495073698489ca5424f2b" in script
    assert "abe9c616523f2486bf1b7be69dfcfda6193d534b40b3f4afc49d8ebf3e40ce5a" in script
    assert "TWIRL_FM0_TEMPORAL_PANEL_DIR:?" in script
    assert "TWIRL_FM0_MATCHED_CANARY_EVALUATION_PLAN_DIR:?" in script
    assert 'git -C "${REPO}" rev-parse HEAD' in script
    assert "status --porcelain=v1 --untracked-files=all" in script
    assert (
        script.count(
            "scripts/stage5_validation/"
            "freeze_twirl_fm0_3_matched_canary_evaluation_plan.py"
        )
        == 1
    )
    assert "validate_matched_canary_evaluation_plan" in script
    assert 'result.receipt["schedule"]["n_rows"] == 1440' in script
    assert "counts[split][cohort] == 240" in script
    assert "train_twirl_fm0.py" not in script
    assert "import torch" not in script
