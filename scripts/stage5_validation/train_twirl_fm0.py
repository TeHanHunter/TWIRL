#!/usr/bin/env python3
"""Run one checksum-bound synthetic or real TWIRL-FM0 training invocation."""

from __future__ import annotations

import argparse
import json
import os
import sys
from dataclasses import asdict
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from twirl.models.fm0.composite_release import (  # noqa: E402
    COMPOSITE_RELEASE_SCHEMA_VERSION,
    FM0CompositeDatasetConfig,
    FM0CompositeReleaseDataset,
)
from twirl.models.fm0.composite_release import (  # noqa: E402
    TRAIN_ROLE as FM0_3_TRAIN_ROLE,
)
from twirl.models.fm0.dataset import (  # noqa: E402
    FM0ReleaseDataset,
    FM0ReleaseDatasetConfig,
    SyntheticFM0Config,
    SyntheticFM0Dataset,
)
from twirl.models.fm0.matched_canary_payload_plan import (  # noqa: E402
    validate_matched_canary_payload_plan,
)
from twirl.models.fm0.model import (  # noqa: E402
    architecture_for_variant,
    build_fm0_model,
    count_trainable_parameters,
)
from twirl.models.fm0.training import (  # noqa: E402
    CADENCE_VICREG_OBJECTIVE_IDENTITY,
    FM0_3_MASK_SPAN_RANGE,
    FM0_3_MASK_TARGET_FRACTION,
    FM0_3_WINDOW_LENGTH,
    FM0OptimizationConfig,
    run_real_training,
    run_synthetic_training,
    seed_everything,
)
from twirl.models.fm0.validation import (  # noqa: E402
    FM0_3_TEMPORAL_PANEL_RECEIPT_SHA256,
    REAL_RUN_CONTRACT_SCHEMA_VERSION,
    REAL_RUN_SUMMARY_SCHEMA_VERSION,
    RUN_CONTRACT_SCHEMA_VERSION,
    RUN_SUMMARY_SCHEMA_VERSION,
    read_json,
    require_clean_git_revision,
    sha256_file,
    validate_frozen_authorities,
    validate_fm0_3_prestart_smoke,
    validate_real_run_release,
    validate_run_release,
    write_json_with_sha256,
    write_sha256_sidecar,
)

FM0_3_CONFIG_SCHEMA_VERSION = "twirl_fm0_3_cadence_objective_config_v1"
FM0_2_CONFIG_SCHEMA_VERSION = "twirl_fm0_2_objective_canary_config_v1"
RUNTIME_ONLY_STOP_CONFIG_SCHEMAS = frozenset(
    {FM0_2_CONFIG_SCHEMA_VERSION, FM0_3_CONFIG_SCHEMA_VERSION}
)
FM0_3_AUTHORIZED_VARIANTS = ("TWIRL-FM0.3.1", "TWIRL-FM0.3.2")
FM0_3_CONFIG_STATUS = "frozen_matched_canary_authorized_not_started"
FM0_3_AUTHORIZED_MAX_STOP = 2_000
FM0_3_SYNTHETIC_SMOKE_STEP = 8
FM0_3_REAL_RESTART_STEP = 64
FM0_3_MILESTONE_STEPS = (0, 64, 2_000)
FM0_3_AUTHORIZED_STOPS = (64, 2_000)
FM0_3_VARIANT_ARCHITECTURES = {
    "TWIRL-FM0.3.1": "tcn",
    "TWIRL-FM0.3.2": "conformer",
}


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument(
        "--design", type=Path, default=ROOT / "doc" / "foundation_model_design.md"
    )
    parser.add_argument("--freeze-receipt", type=Path, required=True)
    parser.add_argument("--variant", required=True)
    parser.add_argument("--architecture", choices=("tcn", "conformer"))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--synthetic-smoke", action="store_true")
    parser.add_argument("--input-release", type=Path)
    parser.add_argument("--input-release-receipt", type=Path)
    parser.add_argument("--input-reuse-receipt", type=Path)
    parser.add_argument("--evaluation-plan", type=Path)
    parser.add_argument("--evaluation-plan-receipt-sha256")
    parser.add_argument("--prestart-smoke-run", type=Path)
    parser.add_argument("--prestart-smoke-summary-sha256")
    parser.add_argument("--expected-git-sha")
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument(
        "--stop-after-step",
        "--target-step",
        "--max-steps",
        dest="target_step",
        type=int,
        default=1,
        help=(
            "runtime-only invocation stop; FM0.2/FM0.3 keep the frozen "
            "20,000-step optimizer horizon in the invariant run contract"
        ),
    )
    parser.add_argument("--micro-batch-windows", type=int, default=2)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--precision", choices=("fp32", "bf16"), default="fp32")
    parser.add_argument("--resume-checkpoint", type=Path)
    return parser


def _load_config(path: Path) -> dict[str, object]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("FM0 config must be a YAML mapping")
    if payload.get("schema_version") not in {
        "twirl_fm0_1_poc_config_v1",
        FM0_2_CONFIG_SCHEMA_VERSION,
        FM0_3_CONFIG_SCHEMA_VERSION,
    }:
        raise ValueError("FM0 config schema mismatch")
    return payload


def _objective_settings(
    config: dict[str, object],
    variant_payload: dict[str, object],
) -> tuple[str, bool, bool, str | None]:
    """Resolve only an objective explicitly authorized by its config schema."""

    schema_version = config.get("schema_version")
    objective_name = str(variant_payload.get("objective", ""))
    objective_identity = None
    if schema_version == FM0_3_CONFIG_SCHEMA_VERSION:
        if objective_name != CADENCE_VICREG_OBJECTIVE_IDENTITY:
            raise ValueError(
                "FM0.3 requires the frozen cadence-level objective identity"
            )
        objective_identity = CADENCE_VICREG_OBJECTIVE_IDENTITY
        use_vicreg = True
    elif objective_name == "masked_reconstruction":
        use_vicreg = False
    elif objective_name in {
        "masked_reconstruction_plus_same_window_vicreg",
        "dual_mask_reconstruction_plus_same_window_vicreg",
    }:
        use_vicreg = True
    else:
        raise ValueError(f"unsupported FM0 objective: {objective_name!r}")

    if schema_version == "twirl_fm0_1_poc_config_v1":
        reconstruct_second_view = use_vicreg
    else:
        objective_payload = config.get("objective")
        if not isinstance(objective_payload, dict):
            raise ValueError("FM0 reconstruction objective is malformed")
        if schema_version == FM0_3_CONFIG_SCHEMA_VERSION and (
            objective_payload.get("name") != objective_name
        ):
            raise ValueError("FM0.3 objective declarations disagree")
        reconstruction_payload = objective_payload.get("reconstruction", {})
        if not isinstance(reconstruction_payload, dict):
            raise ValueError("FM0 reconstruction objective is malformed")
        optimized_masks = reconstruction_payload.get("optimized_mask_views")
        if optimized_masks == ["first"]:
            reconstruct_second_view = False
        elif optimized_masks == ["first", "second"]:
            reconstruct_second_view = True
        else:
            raise ValueError("FM0 optimized mask-view declaration is invalid")
    if objective_identity is not None and reconstruct_second_view:
        raise ValueError("FM0.3 cadence objective optimizes only the first mask")
    return (
        objective_name,
        use_vicreg,
        reconstruct_second_view,
        objective_identity,
    )


def _dataset_geometry(schema_version: object) -> dict[str, object]:
    """Return the schema-bound model window and masking geometry."""

    if schema_version == FM0_3_CONFIG_SCHEMA_VERSION:
        return {
            "window_length": FM0_3_WINDOW_LENGTH,
            "mask_target_fraction": FM0_3_MASK_TARGET_FRACTION,
            "mask_span_range": FM0_3_MASK_SPAN_RANGE,
        }
    return {}


def _require_lowercase_sha256(value: object, *, label: str) -> str:
    text = str(value)
    if len(text) != 64 or any(
        character not in "0123456789abcdef" for character in text
    ):
        raise ValueError(f"{label} must be a lowercase SHA-256")
    return text


def _is_exact_int(value: object, expected: int) -> bool:
    return type(value) is int and value == expected


def _fm0_3_evaluation_freeze_contract() -> dict[str, object]:
    """Return the exact pre-checkpoint evaluation-freeze contract."""

    return {
        "required_before_real_checkpoint": True,
        "shared_between_variants": True,
        "temporal_panel_receipt_sha256": (
            FM0_3_TEMPORAL_PANEL_RECEIPT_SHA256
        ),
        "n_identity_crops": 1_440,
        "crop_length_cadences": 128,
        "crop_must_be_unpadded_and_within_one_segment": True,
        "minimum_joint_time_and_both_adp_valid_cadences": 103,
        "complete_joint_valid_support_indices_inclusive": [60, 68],
        "ineligible_preselected_identity_policy": "fail_closed_no_replacement",
        "splits": {
            "probe_train_sectors": [66, 67, 68, 69, 70, 71],
            "validation_sectors": [72, 73, 74],
            "fresh_test_sector": 77,
            "cohorts": ["repeated", "new"],
            "component_pairs_per_split_and_cohort": 240,
        },
        "injection": {
            "center_index_zero_based": 64,
            "durations_cadences": [1, 3, 9],
            "fractional_depths": [0.01, 0.03, 0.1, 0.3],
            "profile": "symmetric_trapezoid_sampled_at_cadence_centers",
            "ingress_egress_fraction_of_duration": 1.0 / 3.0,
            "flux_views_only": ["adp_1x1", "adp_3x3"],
            "formula": "injected=(1+flux)*(1-profile)-1",
        },
        "inference": {
            "temporal_and_reconstruction_masks_all_false": True,
            "readout": "h_cadence_token_64_only",
            "token_or_window_pooling": False,
            "truth_or_search_features_allowed": False,
        },
        "probe": {
            "family": "full_batch_adam_linear_logit",
            "train_only_dimensionwise_standardization": True,
            "epochs": 400,
            "learning_rate": 0.02,
            "l2_weight": 0.001,
            "seed": 560306,
            "hyperparameter_tuning": False,
        },
        "metric": {
            "primary": "sample_roc_auc",
            "bootstrap_method": "component_pair_clustered",
            "bootstrap_replicates": 1_000,
            "bootstrap_seed": 560305,
            "same_resamples_across_arms_and_paired_deltas": True,
        },
        "bls_or_periodic_search": False,
    }


def _validate_fm0_3_execution_contract(
    config: dict[str, object],
    *,
    variant: str,
    target_step: int,
    synthetic_smoke: bool,
    precision: str,
    device: str,
    optimizer_horizon: int,
) -> dict[str, object]:
    """Validate the complete, bounded FM0.3 invocation authorization."""

    authorization = config.get("authorization")
    if config.get("status") != FM0_3_CONFIG_STATUS or authorization != {
        "scientific_contract_frozen": True,
        "gpu_training_authorized": True,
        "gpu_submission_requires_prestart_smoke": True,
        "payload_screened_evaluation_freeze_required": True,
        "authorized_variants": list(FM0_3_AUTHORIZED_VARIANTS),
        "authorized_max_stop_after_step": FM0_3_AUTHORIZED_MAX_STOP,
        "sealed_test_access_authorized": False,
        "production_model_claim": False,
        "foundation_model_claim": False,
    }:
        raise ValueError(
            "FM0.3 config does not authorize only the frozen matched canary"
        )
    if variant not in FM0_3_AUTHORIZED_VARIANTS:
        raise ValueError("only the matched FM0.3.1/0.3.2 variants are authorized")

    variants = config.get("variants")
    if not isinstance(variants, dict) or set(variants) != set(
        FM0_3_VARIANT_ARCHITECTURES
    ):
        raise ValueError("FM0.3 config does not define only the matched variants")
    for variant_name, architecture in FM0_3_VARIANT_ARCHITECTURES.items():
        payload = variants.get(variant_name)
        if (
            not isinstance(payload, dict)
            or payload.get("architecture") != architecture
            or payload.get("objective") != CADENCE_VICREG_OBJECTIVE_IDENTITY
            or payload.get("reconstruction_input")
            != "contextual_h_cadence_without_stem_skip"
        ):
            raise ValueError("FM0.3 matched-variant architecture contract differs")

    if config.get("evaluation_freeze_contract") != (
        _fm0_3_evaluation_freeze_contract()
    ):
        raise ValueError("FM0.3 evaluation-freeze contract differs")

    data_contract = config.get("data_contract")
    if (
        not isinstance(data_contract, dict)
        or data_contract.get("gradient_role") != FM0_3_TRAIN_ROLE
        or not _is_exact_int(
            data_contract.get("temporal_holdout_gradient_rows_allowed"), 0
        )
        or not _is_exact_int(data_contract.get("sealed_rows_allowed"), 0)
        or data_contract.get("host_first_sampling") is not True
        or data_contract.get("flux_views") != ["adp_1x1", "adp_3x3"]
    ):
        raise ValueError("FM0.3 config does not authorize only the training role")

    cadence_contract = config.get("cadence_contract")
    if (
        not isinstance(cadence_contract, dict)
        or not _is_exact_int(cadence_contract.get("nominal_cadence_seconds"), 200)
        or not _is_exact_int(
            cadence_contract.get("context_length_cadences"), FM0_3_WINDOW_LENGTH
        )
        or cadence_contract.get("one_encoder_token_per_cadence") is not True
        or not _is_exact_int(cadence_contract.get("patch_stride"), 1)
        or cadence_contract.get("cadence_averaging") is not False
        or cadence_contract.get("cadence_merging") is not False
        or cadence_contract.get("temporal_downsampling") is not False
        or cadence_contract.get("temporal_pooling") is not False
        or cadence_contract.get("optimized_representation") != "h_cadence"
        or cadence_contract.get("pooled_window_objective") is not False
        or cadence_contract.get("mask_target_fraction") != FM0_3_MASK_TARGET_FRACTION
        or cadence_contract.get("mask_span_range_cadences")
        != list(FM0_3_MASK_SPAN_RANGE)
    ):
        raise ValueError("FM0.3 native-cadence geometry contract differs")

    if (
        isinstance(target_step, bool)
        or not isinstance(target_step, int)
        or target_step <= 0
        or target_step > FM0_3_AUTHORIZED_MAX_STOP
    ):
        raise ValueError("FM0.3 invocation exceeds the authorized canary stop")
    if optimizer_horizon < FM0_3_AUTHORIZED_MAX_STOP:
        raise ValueError("FM0.3 optimizer horizon cannot reach the authorized canary")

    canary = config.get("canary")
    if canary != {
        "fp32_synthetic_smoke_steps": FM0_3_SYNTHETIC_SMOKE_STEP,
        "bf16_real_throughput_and_restart_step": FM0_3_REAL_RESTART_STEP,
        "immutable_milestone_steps": list(FM0_3_MILESTONE_STEPS),
        "authorized_stop_after_steps": list(FM0_3_AUTHORIZED_STOPS),
        "one_h200_per_variant": True,
    }:
        raise ValueError(
            "FM0.3 config lacks the exact matched-canary execution contract"
        )
    if synthetic_smoke:
        if target_step != FM0_3_SYNTHETIC_SMOKE_STEP:
            raise ValueError("FM0.3 synthetic smoke must stop at step 8")
        if precision != "fp32":
            raise ValueError("FM0.3 synthetic smoke must use FP32")
    else:
        if target_step not in FM0_3_AUTHORIZED_STOPS:
            raise ValueError("FM0.3 real-data stop is not preauthorized")
        if precision != "bf16":
            raise ValueError("FM0.3 real-data matched canary must use BF16")
    if str(device) not in {"cuda", "cuda:0"}:
        raise ValueError("FM0.3 training invocations require one CUDA device")
    return canary


def _validate_fm0_3_resume_invocation(
    *,
    target_step: int,
    synthetic_smoke: bool,
    output_dir: Path,
    resume_checkpoint: Path | None,
) -> int | None:
    """Enforce the frozen fresh-64 then exact-64-to-2000 sequence."""

    if synthetic_smoke or target_step == FM0_3_REAL_RESTART_STEP:
        if resume_checkpoint is not None:
            raise ValueError("FM0.3 step 8 and step 64 must start fresh")
        return None
    if target_step != FM0_3_AUTHORIZED_MAX_STOP:
        raise ValueError("FM0.3 resume sequence received an unauthorized stop")
    if resume_checkpoint is None:
        raise ValueError("FM0.3 step 2000 requires the step-64 checkpoint")
    output = output_dir.resolve(strict=True)
    expected = (output / "checkpoint_step_00000064.pt").resolve(strict=True)
    observed = resume_checkpoint.resolve(strict=True)
    if observed != expected:
        raise ValueError(
            "FM0.3 step 2000 must resume this run's exact step-64 checkpoint"
        )
    return FM0_3_REAL_RESTART_STEP


def _require_fm0_3_expected_git_sha(value: str | None) -> str:
    """Require the exact clean code identity for every FM0.3 invocation."""

    if value is None:
        raise ValueError("FM0.3 requires --expected-git-sha")
    if len(value) != 40 or any(character not in "0123456789abcdef" for character in value):
        raise ValueError("FM0.3 expected Git SHA must be 40 lowercase hex characters")
    return value


def _build_fm0_3_composite_dataset(
    *,
    config: dict[str, object],
    release_root: Path,
    receipt_path: Path,
    receipt_sha256: str,
    variant: str,
    seed: int,
    windows_per_epoch: int,
) -> tuple[FM0CompositeReleaseDataset, dict[str, object]]:
    """Construct the exact config-bound training role of the FM0.3 release."""

    composite = config.get("composite_input")
    if not isinstance(composite, dict):
        raise ValueError("FM0.3 config lacks its composite-input binding")
    expected_receipt = _require_lowercase_sha256(
        composite.get("receipt_sha256"), label="FM0.3 composite receipt hash"
    )
    source_hash = _require_lowercase_sha256(
        composite.get("source_bindings_sha256"),
        label="FM0.3 source-bindings hash",
    )
    role_hash = _require_lowercase_sha256(
        composite.get("role_index_sha256"), label="FM0.3 role-index hash"
    )
    root = release_root.resolve(strict=True)
    bound_receipt_path = (root / "receipt.json").resolve(strict=True)
    if receipt_path.resolve(strict=True) != bound_receipt_path:
        raise ValueError("FM0.3 receipt must be the requested composite receipt")
    if receipt_sha256 != expected_receipt:
        raise ValueError("FM0.3 config does not bind the composite receipt")
    receipt = read_json(bound_receipt_path)
    limits = receipt.get("limits")
    if (
        receipt.get("schema_version") != COMPOSITE_RELEASE_SCHEMA_VERSION
        or receipt.get("passed") is not True
        or not isinstance(limits, dict)
        or limits.get("scientific_training_eligible") is not True
        or limits.get("sealed_rows_selected") != 0
        or limits.get("sealed_test_access_authorized") is not False
    ):
        raise ValueError("FM0.3 composite receipt is not training-eligible")

    source_path = (root / "source_bindings.csv").resolve(strict=True)
    role_path = (root / "role_index.csv").resolve(strict=True)
    if sha256_file(source_path) != source_hash:
        raise ValueError("FM0.3 source-bindings file differs from the frozen config")
    if sha256_file(role_path) != role_hash:
        raise ValueError("FM0.3 role-index file differs from the frozen config")
    dataset = FM0CompositeReleaseDataset(
        FM0CompositeDatasetConfig(
            composite_root=str(root),
            receipt_sha256=expected_receipt,
            source_bindings_sha256=source_hash,
            role_index_sha256=role_hash,
            variant=variant,
            role=FM0_3_TRAIN_ROLE,
            seed=seed,
            windows_per_epoch=windows_per_epoch,
            window_length=FM0_3_WINDOW_LENGTH,
            mask_target_fraction=FM0_3_MASK_TARGET_FRACTION,
            mask_span_range=FM0_3_MASK_SPAN_RANGE,
        )
    )
    binding = {
        "kind": "fm0_3_composite_release",
        "release_root": str(root),
        "receipt_path": str(bound_receipt_path),
        "receipt_sha256": expected_receipt,
        "source_bindings_path": str(source_path),
        "source_bindings_sha256": source_hash,
        "role_index_path": str(role_path),
        "role_index_sha256": role_hash,
        "role": FM0_3_TRAIN_ROLE,
        "n_sources": dataset.contract["n_sources"],
        "n_observations": dataset.contract["n_observations"],
        "n_excluded_missing_required_views": dataset.contract[
            "n_excluded_missing_required_views"
        ],
    }
    return dataset, binding


def _build_fm0_3_evaluation_plan_binding(
    *,
    config: dict[str, object],
    plan_root: Path,
    receipt_sha256: str,
    expected_git_sha: str,
) -> dict[str, object]:
    """Validate and bind the shared payload-screened evaluation schedule."""

    expected_receipt = _require_lowercase_sha256(
        receipt_sha256, label="FM0.3 evaluation-plan receipt hash"
    )
    result = validate_matched_canary_payload_plan(
        plan_root,
        expected_receipt_sha256=expected_receipt,
        require_read_only=True,
    )
    receipt = result.receipt
    if receipt.get("producer_git_sha") != expected_git_sha:
        raise ValueError("FM0.3 evaluation plan was not frozen by this exact code")

    identity = receipt.get("identity_plan")
    payload_bindings = receipt.get("payload_bindings")
    if not isinstance(identity, dict) or not isinstance(payload_bindings, dict):
        raise ValueError("FM0.3 evaluation-plan binding is malformed")
    input_authorities = identity.get("input_authorities")
    temporal_authority = (
        input_authorities.get("temporal_panel")
        if isinstance(input_authorities, dict)
        else None
    )
    composite_authority = (
        input_authorities.get("composite_release")
        if isinstance(input_authorities, dict)
        else None
    )
    evaluation_contract = config.get("evaluation_freeze_contract")
    composite_config = config.get("composite_input")
    if (
        not isinstance(temporal_authority, dict)
        or not isinstance(composite_authority, dict)
        or not isinstance(evaluation_contract, dict)
        or not isinstance(composite_config, dict)
    ):
        raise ValueError("FM0.3 evaluation plan lacks a frozen input authority")
    if temporal_authority.get("receipt_sha256") != evaluation_contract.get(
        "temporal_panel_receipt_sha256"
    ):
        raise ValueError(
            "FM0.3 evaluation-plan temporal-panel receipt differs from training"
        )
    for name in (
        "receipt_sha256",
        "source_bindings_sha256",
        "role_index_sha256",
    ):
        if composite_authority.get(name) != composite_config.get(name):
            raise ValueError(
                f"FM0.3 evaluation-plan composite {name} differs from training"
            )

    return {
        "kind": "fm0_3_payload_screened_evaluation_plan",
        "root": str(result.root),
        "receipt_path": str(result.receipt_path),
        "receipt_sha256": result.receipt_sha256,
        "schedule_path": str(result.schedule_path),
        "schedule_sha256": result.schedule_sha256,
        "producer_git_sha": expected_git_sha,
        "identity_plan_receipt_sha256": identity["receipt_sha256"],
        "identity_plan_schedule_sha256": identity["schedule_sha256"],
        "identity_plan_producer_git_sha": identity["producer_git_sha"],
        "temporal_panel_receipt_sha256": temporal_authority["receipt_sha256"],
        "temporal_panel_sha256": temporal_authority["panel_sha256"],
        "temporal_panel_sector_bindings_sha256": temporal_authority[
            "sector_bindings_sha256"
        ],
        "source_shard_bindings_sha256": payload_bindings[
            "source_shard_bindings_sha256"
        ],
        "crop_payload_bindings_sha256": payload_bindings[
            "crop_payload_bindings_sha256"
        ],
        "n_crops": payload_bindings["n_crops_frozen"],
    }


def _build_fm0_3_prestart_smoke_binding(
    *,
    is_fm0_3: bool,
    synthetic_smoke: bool,
    smoke_run: Path | None,
    smoke_summary_sha256: str | None,
    campaign_id: object,
    variant: str,
    architecture: str,
    expected_git_sha: str | None,
    authorities: dict[str, object],
) -> dict[str, object] | None:
    """Require the exact same-variant smoke before any real FM0.3 run."""

    supplied = (smoke_run, smoke_summary_sha256)
    if (smoke_run is None) != (smoke_summary_sha256 is None):
        raise ValueError(
            "prestart smoke root and summary SHA-256 must be supplied together"
        )
    if not is_fm0_3:
        if any(value is not None for value in supplied):
            raise ValueError("only FM0.3 accepts a prestart-smoke binding")
        return None
    if synthetic_smoke:
        if any(value is not None for value in supplied):
            raise ValueError("FM0.3 synthetic smoke cannot depend on another smoke")
        return None
    if smoke_run is None:
        raise ValueError("FM0.3 real training requires its same-variant 8-step smoke")
    if expected_git_sha is None:  # pragma: no cover - checked by caller
        raise RuntimeError("FM0.3 exact Git identity was not resolved")
    return validate_fm0_3_prestart_smoke(
        smoke_run,
        expected_summary_sha256=str(smoke_summary_sha256),
        expected_campaign_id=campaign_id,
        expected_variant=variant,
        expected_architecture=architecture,
        expected_git_sha=expected_git_sha,
        expected_authorities=authorities,
    )


def _bind_training_stop_contract(
    run_contract: dict[str, object],
    *,
    schema_version: object,
    optimizer_horizon: int,
    invocation_target_step: int,
    canary_payload: dict[str, object] | None = None,
) -> None:
    """Bind invariant training science while leaving staged stops resumable."""

    if schema_version in {FM0_2_CONFIG_SCHEMA_VERSION, FM0_3_CONFIG_SCHEMA_VERSION}:
        if canary_payload is None:
            raise ValueError("FM0 canary stop contract lacks its declaration")
        run_contract.update(
            {
                "training_horizon_step": optimizer_horizon,
                "immutable_milestone_steps": [
                    int(value) for value in canary_payload["immutable_milestone_steps"]
                ],
                "authorized_stop_after_steps": [
                    int(value)
                    for value in canary_payload["authorized_stop_after_steps"]
                ],
                "synthetic_smoke_step": int(
                    canary_payload["fp32_synthetic_smoke_steps"]
                ),
            }
        )
        if schema_version == FM0_3_CONFIG_SCHEMA_VERSION:
            run_contract[
                "stop_after_step_is_execution_state_not_scientific_contract"
            ] = True
    else:
        # Preserve the byte-level FM0.1 run/checkpoint contract.
        run_contract["target_step"] = invocation_target_step


def main() -> int:
    args = _parser().parse_args()
    if args.synthetic_smoke == (args.input_release is not None):
        raise SystemExit("choose exactly one of --synthetic-smoke or --input-release")
    authorities = validate_frozen_authorities(
        design_path=args.design,
        config_path=args.config,
        freeze_receipt_path=args.freeze_receipt,
    )
    config = _load_config(args.config)
    schema_version = config.get("schema_version")
    is_fm0_2 = schema_version == FM0_2_CONFIG_SCHEMA_VERSION
    is_fm0_3 = schema_version == FM0_3_CONFIG_SCHEMA_VERSION
    runtime_only_stop = schema_version in RUNTIME_ONLY_STOP_CONFIG_SCHEMAS
    if is_fm0_2:
        authorization = config.get("authorization")
        if (
            config.get("status") != "frozen_gated_canary_authorized_not_started"
            or not isinstance(authorization, dict)
            or authorization.get("scientific_contract_frozen") is not True
            or authorization.get("gpu_training_authorized") is not True
            or authorization.get("gpu_submission_requires_all_prestart_gates")
            is not True
            or authorization.get("sealed_test_access_authorized") is not False
            or authorization.get("production_model_claim") is not False
            or authorization.get("foundation_model_claim") is not False
        ):
            raise ValueError("FM0.2 config does not authorize only the gated canary")
        if args.variant != authorization.get("authorized_variant"):
            raise ValueError("only the frozen FM0.2.1 canary variant is authorized")
        maximum_stop = authorization.get("authorized_max_stop_after_step")
        if (
            isinstance(maximum_stop, bool)
            or not isinstance(maximum_stop, int)
            or args.target_step > maximum_stop
        ):
            raise ValueError("FM0.2 invocation exceeds the authorized canary stop")
    variants = config.get("variants", {})
    if not isinstance(variants, dict) or args.variant not in variants:
        raise ValueError("requested variant is absent from the frozen config")
    variant_payload = variants[args.variant]
    if not isinstance(variant_payload, dict):
        raise ValueError("requested FM0 variant is malformed")
    (
        objective_name,
        use_vicreg,
        reconstruct_second_view,
        objective_identity,
    ) = _objective_settings(config, variant_payload)
    dataset_geometry = _dataset_geometry(schema_version)
    optimization_payload = config.get("optimization", {})
    if not isinstance(optimization_payload, dict):
        raise ValueError("frozen optimization config is malformed")
    frozen_seeds = tuple(int(value) for value in optimization_payload.get("seeds", []))
    if args.seed not in frozen_seeds:
        raise ValueError(f"seed must be one of the frozen values {frozen_seeds}")
    resolved_architecture = architecture_for_variant(
        args.variant, development_winner=args.architecture
    )
    optimization = FM0OptimizationConfig(
        learning_rate=float(optimization_payload["learning_rate"]),
        weight_decay=float(optimization_payload["weight_decay"]),
        warmup_steps=int(optimization_payload["warmup_steps"]),
        max_optimizer_steps=int(optimization_payload["max_optimizer_steps"]),
        effective_batch_windows=int(optimization_payload["effective_batch_windows"]),
        huber_delta=float(
            config["objective"]["reconstruction"]["huber_delta_fractional_flux"]
        ),
        vicreg_total_weight=float(config["objective"]["vicreg"]["total_weight"]),
        vicreg_invariance_weight=float(
            config["objective"]["vicreg"]["invariance_weight"]
        ),
        vicreg_variance_weight=float(config["objective"]["vicreg"]["variance_weight"]),
        vicreg_covariance_weight=float(
            config["objective"]["vicreg"]["covariance_weight"]
        ),
    )
    canary_payload: dict[str, object] = {}
    expected_resume_step: int | None = None
    fm0_3_expected_git_sha: str | None = None
    if is_fm0_2:
        raw_canary = config.get("canary")
        if not isinstance(raw_canary, dict):
            raise ValueError("FM0.2 config lacks its canary execution contract")
        canary_payload = raw_canary
        if int(raw_canary.get("run_contract_target_step", -1)) != (
            optimization.max_optimizer_steps
        ):
            raise ValueError("FM0.2 optimizer horizon differs from its run contract")
        if args.synthetic_smoke:
            if args.target_step != int(raw_canary["fp32_synthetic_smoke_steps"]):
                raise ValueError("FM0.2 synthetic smoke must stop at the frozen step")
            if args.precision != "fp32":
                raise ValueError("FM0.2 synthetic smoke must use FP32")
        else:
            allowed_stops = tuple(
                int(value) for value in raw_canary["authorized_stop_after_steps"]
            )
            if args.target_step not in allowed_stops:
                raise ValueError("FM0.2 real-data stop is not preauthorized")
            if args.precision != "bf16":
                raise ValueError("FM0.2 real-data canary must use BF16")
        if not str(args.device).startswith("cuda"):
            raise ValueError("FM0.2 training invocations require one CUDA device")
    elif is_fm0_3:
        canary_payload = _validate_fm0_3_execution_contract(
            config,
            variant=args.variant,
            target_step=args.target_step,
            synthetic_smoke=args.synthetic_smoke,
            precision=args.precision,
            device=args.device,
            optimizer_horizon=optimization.max_optimizer_steps,
        )
        expected_resume_step = _validate_fm0_3_resume_invocation(
            target_step=args.target_step,
            synthetic_smoke=args.synthetic_smoke,
            output_dir=args.output_dir,
            resume_checkpoint=args.resume_checkpoint,
        )
        fm0_3_expected_git_sha = _require_fm0_3_expected_git_sha(
            args.expected_git_sha
        )

    git_sha = None
    if args.expected_git_sha:
        git_sha = require_clean_git_revision(ROOT, args.expected_git_sha)

    prestart_smoke_binding = _build_fm0_3_prestart_smoke_binding(
        is_fm0_3=is_fm0_3,
        synthetic_smoke=args.synthetic_smoke,
        smoke_run=args.prestart_smoke_run,
        smoke_summary_sha256=args.prestart_smoke_summary_sha256,
        campaign_id=config.get("campaign_id"),
        variant=args.variant,
        architecture=resolved_architecture,
        expected_git_sha=fm0_3_expected_git_sha,
        authorities=authorities,
    )

    evaluation_plan_binding = None
    evaluation_plan_arguments = (
        args.evaluation_plan,
        args.evaluation_plan_receipt_sha256,
    )
    if (args.evaluation_plan is None) != (
        args.evaluation_plan_receipt_sha256 is None
    ):
        raise ValueError(
            "evaluation plan root and receipt SHA-256 must be supplied together"
        )
    if not is_fm0_3 and any(value is not None for value in evaluation_plan_arguments):
        raise ValueError("only FM0.3 accepts a matched-canary evaluation plan")
    if is_fm0_3 and args.synthetic_smoke:
        if any(value is not None for value in evaluation_plan_arguments):
            raise ValueError("FM0.3 synthetic smoke must not consume evaluation data")
    elif is_fm0_3:
        if args.evaluation_plan is None:
            raise ValueError(
                "FM0.3 real training requires the frozen payload evaluation plan"
            )
        if fm0_3_expected_git_sha is None:  # pragma: no cover - checked above
            raise RuntimeError("FM0.3 exact Git identity was not resolved")
        evaluation_plan_binding = _build_fm0_3_evaluation_plan_binding(
            config=config,
            plan_root=args.evaluation_plan,
            receipt_sha256=str(args.evaluation_plan_receipt_sha256),
            expected_git_sha=fm0_3_expected_git_sha,
        )
    input_receipt = None
    receipt_payload = None
    receipt_argument = args.input_release_receipt
    if receipt_argument is None and os.environ.get("TWIRL_FM0_INPUT_RECEIPT"):
        receipt_argument = Path(os.environ["TWIRL_FM0_INPUT_RECEIPT"])
    if receipt_argument is not None:
        receipt_path = receipt_argument.resolve(strict=True)
        expected_receipt_hash = os.environ.get("TWIRL_FM0_INPUT_RECEIPT_SHA256")
        observed_receipt_hash = sha256_file(receipt_path)
        if expected_receipt_hash and observed_receipt_hash != expected_receipt_hash:
            raise ValueError(
                "input-release receipt differs from expected environment hash"
            )
        input_receipt = {
            "path": str(receipt_path),
            "sha256": observed_receipt_hash,
        }
        receipt_payload = read_json(receipt_path)
    if not args.synthetic_smoke and input_receipt is None:
        raise ValueError("real-data training requires an input-release receipt")

    input_reuse_receipt = None
    input_reuse_payload = None
    reuse_argument = args.input_reuse_receipt
    if reuse_argument is None and os.environ.get("TWIRL_FM0_INPUT_REUSE_RECEIPT"):
        reuse_argument = Path(os.environ["TWIRL_FM0_INPUT_REUSE_RECEIPT"])
    if reuse_argument is not None:
        reuse_path = reuse_argument.resolve(strict=True)
        expected_reuse_hash = os.environ.get("TWIRL_FM0_INPUT_REUSE_RECEIPT_SHA256")
        observed_reuse_hash = sha256_file(reuse_path)
        if expected_reuse_hash and observed_reuse_hash != expected_reuse_hash:
            raise ValueError(
                "input-reuse receipt differs from expected environment hash"
            )
        input_reuse_receipt = {
            "path": str(reuse_path),
            "sha256": observed_reuse_hash,
        }
        input_reuse_payload = read_json(reuse_path)
    if is_fm0_2 and not args.synthetic_smoke and input_reuse_receipt is None:
        raise ValueError("FM0.2 real-data training requires its reuse receipt")
    if is_fm0_3 and input_reuse_receipt is not None:
        raise ValueError("FM0.3 composite training does not accept a reuse receipt")

    release_dataset = None
    release_binding = None
    if not args.synthetic_smoke and is_fm0_3:
        if input_receipt is None:  # pragma: no cover - checked above
            raise RuntimeError("FM0.3 composite receipt was not loaded")
        release_dataset, release_binding = _build_fm0_3_composite_dataset(
            config=config,
            release_root=args.input_release,
            receipt_path=Path(input_receipt["path"]),
            receipt_sha256=input_receipt["sha256"],
            variant=args.variant,
            seed=args.seed,
            windows_per_epoch=(
                int(optimization.max_optimizer_steps)
                * int(optimization.effective_batch_windows)
            ),
        )
    if not args.synthetic_smoke and not is_fm0_3:
        if receipt_payload is None or receipt_payload.get("schema_version") != (
            "twirl_fm0_1_input_release_receipt_v1"
        ):
            raise ValueError("real-data input receipt has the wrong schema")
        if (
            receipt_payload.get("passed") is not True
            or receipt_payload.get("scientific_training_eligible") is not True
        ):
            raise ValueError("real-data input receipt does not authorize training")
        if config.get("schema_version") == "twirl_fm0_1_poc_config_v1":
            if receipt_payload.get("science_bindings") != authorities:
                raise ValueError(
                    "input receipt differs from the frozen science authorities"
                )
        else:
            reuse = config.get("input_release_reuse")
            if not isinstance(reuse, dict):
                raise ValueError("FM0.2 config lacks its input-release reuse binding")
            if reuse.get("fm0_1_input_receipt_sha256") != input_receipt["sha256"]:
                raise ValueError(
                    "FM0.2 config does not bind the upstream input receipt"
                )
            if (
                input_reuse_receipt is None
                or input_reuse_payload is None
                or input_reuse_payload.get("schema_version")
                != "twirl_fm0_2_input_reuse_receipt_v1"
                or input_reuse_payload.get("passed") is not True
                or input_reuse_payload.get("scientific_training_eligible") is not True
                or input_reuse_payload.get("sealed_test", {}).get("access_count") != 0
                or input_reuse_payload.get("upstream_fm0_1_input_receipt", {}).get(
                    "sha256"
                )
                != input_receipt["sha256"]
            ):
                raise ValueError(
                    "FM0.2 input-reuse receipt does not authorize training"
                )
            supporting = authorities.get("supporting_authorities")
            bound_reuse = (
                supporting.get("input_reuse_receipt")
                if isinstance(supporting, dict)
                else None
            )
            if (
                not isinstance(bound_reuse, dict)
                or bound_reuse.get("sha256") != input_reuse_receipt["sha256"]
            ):
                raise ValueError("FM0.2 freeze does not bind the input-reuse receipt")
        release_root = args.input_release.resolve(strict=True)
        manifest_path = (release_root / "manifest.csv").resolve(strict=True)
        receipt_release = receipt_payload.get("release")
        if not isinstance(receipt_release, dict):
            raise ValueError("input receipt release binding is malformed")
        manifest_sha = sha256_file(manifest_path)
        if config.get("schema_version") == "twirl_fm0_2_objective_canary_config_v1":
            reuse = config["input_release_reuse"]
            if reuse.get("manifest_sha256") != manifest_sha:
                raise ValueError("FM0.2 config does not bind the reused manifest")
            if (
                input_reuse_payload is None
                or input_reuse_payload.get("release", {}).get("id")
                != reuse.get("release_id")
                or input_reuse_payload.get("release", {}).get("manifest_sha256")
                != manifest_sha
                or input_reuse_payload.get("release", {}).get("sectors")
                != reuse.get("sectors")
            ):
                raise ValueError("FM0.2 reuse receipt differs from the frozen release")
        if (
            Path(str(receipt_release.get("manifest_path", ""))).resolve(strict=True)
            != manifest_path
            or receipt_release.get("manifest_sha256") != manifest_sha
        ):
            raise ValueError("input receipt does not bind the requested release")
        release_dataset = FM0ReleaseDataset(
            FM0ReleaseDatasetConfig(
                release_root=str(release_root),
                manifest_sha256=manifest_sha,
                variant=args.variant,
                seed=args.seed,
                source_partition="poc_train",
                windows_per_epoch=int(optimization.max_optimizer_steps)
                * int(optimization.effective_batch_windows),
                **dataset_geometry,
            )
        )
        release_binding = {
            "release_root": str(release_root),
            "manifest_path": str(manifest_path),
            "manifest_sha256": manifest_sha,
            "receipt_path": input_receipt["path"],
            "receipt_sha256": input_receipt["sha256"],
            "n_sources": release_dataset.contract["n_sources"],
            "n_observations": release_dataset.contract["n_observations"],
            "source_partition": "poc_train",
            "certifies_full_campaign": bool(
                receipt_payload.get("certifies_full_campaign", False)
            ),
        }
    output = args.output_dir.resolve()
    if args.resume_checkpoint is None:
        if output.exists():
            raise FileExistsError(f"refusing non-fresh output directory: {output}")
        output.mkdir(parents=True)
    elif not output.is_dir():
        raise FileNotFoundError("resume output directory does not exist")

    run_contract = {
        "schema_version": (
            RUN_CONTRACT_SCHEMA_VERSION
            if args.synthetic_smoke
            else REAL_RUN_CONTRACT_SCHEMA_VERSION
        ),
        "campaign_id": config["campaign_id"],
        "variant": args.variant,
        "architecture": resolved_architecture,
        "objective": objective_name,
        "use_vicreg": use_vicreg,
        "reconstruct_second_view": reconstruct_second_view,
        "seed": args.seed,
        "synthetic_smoke": bool(args.synthetic_smoke),
        "real_data_consumed": not args.synthetic_smoke,
        "precision": args.precision,
        "device_request": args.device,
        "micro_batch_windows": args.micro_batch_windows,
        "authorities": authorities,
        "input_release_receipt_precondition": (
            input_receipt if args.synthetic_smoke else None
        ),
        "input_release": release_binding,
        "input_release_reuse": (
            input_reuse_receipt if is_fm0_2 and not args.synthetic_smoke else None
        ),
        "expected_git_sha": git_sha,
        "optimization": asdict(optimization),
    }
    if is_fm0_3:
        run_contract.update(
            {
                "mask_target_fraction": FM0_3_MASK_TARGET_FRACTION,
                "mask_span_range": list(FM0_3_MASK_SPAN_RANGE),
                "evaluation_plan": evaluation_plan_binding,
                "prestart_smoke": prestart_smoke_binding,
            }
        )
    _bind_training_stop_contract(
        run_contract,
        schema_version=schema_version,
        optimizer_horizon=optimization.max_optimizer_steps,
        invocation_target_step=args.target_step,
        canary_payload=canary_payload if (is_fm0_2 or is_fm0_3) else None,
    )
    contract_path = output / "run_contract.json"
    if contract_path.exists():
        if read_json(contract_path) != run_contract:
            raise ValueError("resume run contract differs from existing contract")
        contract_sha = sha256_file(contract_path)
    else:
        contract_sha = write_json_with_sha256(contract_path, run_contract)

    seed_everything(args.seed)
    model = build_fm0_model(
        args.variant,
        development_winner=args.architecture,
        enforce_parameter_budget=True,
    )
    parameter_count = count_trainable_parameters(model)
    if args.synthetic_smoke:
        dataset = SyntheticFM0Dataset(
            SyntheticFM0Config(
                variant=args.variant,
                seed=args.seed,
                **dataset_geometry,
            )
        )
        result = run_synthetic_training(
            model=model,
            dataset=dataset,
            output_dir=output,
            run_contract=run_contract,
            optimization=optimization,
            target_step=args.target_step,
            micro_batch_windows=args.micro_batch_windows,
            device=args.device,
            precision=args.precision,
            use_vicreg=use_vicreg,
            reconstruct_second_view=reconstruct_second_view,
            objective_identity=objective_identity,
            resume_checkpoint=args.resume_checkpoint,
            expected_resume_step=expected_resume_step,
        )
    else:
        if release_dataset is None:  # pragma: no cover - defensive
            raise RuntimeError("real release dataset was not constructed")
        result = run_real_training(
            model=model,
            dataset=release_dataset,
            output_dir=output,
            run_contract=run_contract,
            optimization=optimization,
            target_step=args.target_step,
            micro_batch_windows=args.micro_batch_windows,
            device=args.device,
            precision=args.precision,
            use_vicreg=use_vicreg,
            reconstruct_second_view=reconstruct_second_view,
            objective_identity=objective_identity,
            resume_checkpoint=args.resume_checkpoint,
            expected_resume_step=expected_resume_step,
        )
    checkpoint_path = output / "checkpoint.pt"
    checkpoint_sha = write_sha256_sidecar(checkpoint_path)
    summary = {
        "schema_version": (
            RUN_SUMMARY_SCHEMA_VERSION
            if args.synthetic_smoke
            else REAL_RUN_SUMMARY_SCHEMA_VERSION
        ),
        "passed": True,
        "synthetic_only": bool(args.synthetic_smoke),
        "real_data_consumed": not args.synthetic_smoke,
        "scientific_result": False,
        "variant": args.variant,
        "architecture": model.config.architecture,
        "objective": objective_name,
        "parameter_count": parameter_count,
        "global_step": result["global_step"],
        "requested_stop_after_step": args.target_step if runtime_only_stop else None,
        "final_metrics": result["final_metrics"],
        "precision": result["precision"],
        "device": result["device"],
        "elapsed_seconds_this_invocation": result["elapsed_seconds_this_invocation"],
        "run_contract_sha256": contract_sha,
        "checkpoint_sha256": checkpoint_sha,
        "immutable_milestone_checkpoints": result.get(
            "immutable_milestone_checkpoints", []
        ),
    }
    write_json_with_sha256(output / "summary.json", summary)
    validation = (
        validate_run_release(output)
        if args.synthetic_smoke
        else validate_real_run_release(output)
    )
    print(json.dumps(validation, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
