from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from twirl.models.fm0.training import (
    CADENCE_VICREG_OBJECTIVE_IDENTITY,
    FM0_3_MASK_SPAN_RANGE,
    FM0_3_MASK_TARGET_FRACTION,
    FM0_3_WINDOW_LENGTH,
)

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "stage5_validation" / "train_twirl_fm0.py"
SPEC = importlib.util.spec_from_file_location("train_twirl_fm0_cli_test", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
CLI = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CLI)


def _fm03_config() -> dict[str, object]:
    return {
        "schema_version": CLI.FM0_3_CONFIG_SCHEMA_VERSION,
        "status": CLI.FM0_3_CONFIG_STATUS,
        "authorization": {
            "scientific_contract_frozen": True,
            "gpu_training_authorized": True,
            "gpu_submission_requires_prestart_smoke": True,
            "payload_screened_evaluation_freeze_required": True,
            "authorized_variants": list(CLI.FM0_3_AUTHORIZED_VARIANTS),
            "authorized_max_stop_after_step": CLI.FM0_3_AUTHORIZED_MAX_STOP,
            "sealed_test_access_authorized": False,
            "production_model_claim": False,
            "foundation_model_claim": False,
        },
        "canary": {
            "fp32_synthetic_smoke_steps": CLI.FM0_3_SYNTHETIC_SMOKE_STEP,
            "bf16_real_throughput_and_restart_step": CLI.FM0_3_REAL_RESTART_STEP,
            "immutable_milestone_steps": list(CLI.FM0_3_MILESTONE_STEPS),
            "authorized_stop_after_steps": list(CLI.FM0_3_AUTHORIZED_STOPS),
            "one_h200_per_variant": True,
        },
        "data_contract": {
            "gradient_role": CLI.FM0_3_TRAIN_ROLE,
            "temporal_holdout_gradient_rows_allowed": 0,
            "sealed_rows_allowed": 0,
            "host_first_sampling": True,
            "flux_views": ["adp_1x1", "adp_3x3"],
        },
        "cadence_contract": {
            "nominal_cadence_seconds": 200,
            "context_length_cadences": FM0_3_WINDOW_LENGTH,
            "one_encoder_token_per_cadence": True,
            "patch_stride": 1,
            "cadence_averaging": False,
            "cadence_merging": False,
            "temporal_downsampling": False,
            "temporal_pooling": False,
            "optimized_representation": "h_cadence",
            "pooled_window_objective": False,
            "mask_target_fraction": FM0_3_MASK_TARGET_FRACTION,
            "mask_span_range_cadences": list(FM0_3_MASK_SPAN_RANGE),
        },
        "variants": {
            name: {
                "architecture": architecture,
                "objective": CADENCE_VICREG_OBJECTIVE_IDENTITY,
                "reconstruction_input": (
                    "contextual_h_cadence_without_stem_skip"
                ),
            }
            for name, architecture in CLI.FM0_3_VARIANT_ARCHITECTURES.items()
        },
        "evaluation_freeze_contract": CLI._fm0_3_evaluation_freeze_contract(),
        "objective": {
            "name": CADENCE_VICREG_OBJECTIVE_IDENTITY,
            "reconstruction": {"optimized_mask_views": ["first"]},
        },
    }


def test_cli_accepts_only_the_explicit_fm03_training_schema(tmp_path: Path) -> None:
    path = tmp_path / "fm03.yaml"
    path.write_text(yaml.safe_dump(_fm03_config()), encoding="utf-8")
    assert CLI._load_config(path)["schema_version"] == (
        "twirl_fm0_3_cadence_objective_config_v1"
    )

    with pytest.raises(ValueError, match="config schema mismatch"):
        CLI._load_config(
            ROOT / "configs" / "models" / "twirl_fm0_3_stride1_mechanics_v1.yaml"
        )


def test_cli_resolves_fm03_cadence_objective_without_window_vicreg() -> None:
    settings = CLI._objective_settings(
        _fm03_config(),
        {"objective": CADENCE_VICREG_OBJECTIVE_IDENTITY},
    )
    assert settings == (
        CADENCE_VICREG_OBJECTIVE_IDENTITY,
        True,
        False,
        CADENCE_VICREG_OBJECTIVE_IDENTITY,
    )


@pytest.mark.parametrize(
    ("variant_objective", "optimized_masks", "message"),
    [
        (
            "masked_reconstruction_plus_same_window_vicreg",
            ["first"],
            "cadence-level objective identity",
        ),
        (
            CADENCE_VICREG_OBJECTIVE_IDENTITY,
            ["first", "second"],
            "optimizes only the first mask",
        ),
    ],
)
def test_cli_rejects_fm03_objective_contract_drift(
    variant_objective: str,
    optimized_masks: list[str],
    message: str,
) -> None:
    config = _fm03_config()
    config["objective"]["name"] = variant_objective
    config["objective"]["reconstruction"]["optimized_mask_views"] = optimized_masks
    with pytest.raises(ValueError, match=message):
        CLI._objective_settings(config, {"objective": variant_objective})


def test_cli_freezes_fm03_native_cadence_dataset_geometry() -> None:
    assert CLI._dataset_geometry(CLI.FM0_3_CONFIG_SCHEMA_VERSION) == {
        "window_length": FM0_3_WINDOW_LENGTH,
        "mask_target_fraction": FM0_3_MASK_TARGET_FRACTION,
        "mask_span_range": FM0_3_MASK_SPAN_RANGE,
    }
    assert CLI._dataset_geometry("twirl_fm0_2_objective_canary_config_v1") == {}


def test_cli_keeps_fm03_invocation_stop_out_of_invariant_contract() -> None:
    first: dict[str, object] = {}
    later: dict[str, object] = {}
    for contract, target in ((first, 1), (later, 500)):
        CLI._bind_training_stop_contract(
            contract,
            schema_version=CLI.FM0_3_CONFIG_SCHEMA_VERSION,
            optimizer_horizon=20_000,
            invocation_target_step=target,
            canary_payload=_fm03_config()["canary"],
        )
    assert (
        first
        == later
        == {
            "training_horizon_step": 20_000,
            "immutable_milestone_steps": [0, 64, 2000],
            "authorized_stop_after_steps": [64, 2000],
            "synthetic_smoke_step": 8,
            "stop_after_step_is_execution_state_not_scientific_contract": True,
        }
    )
    assert "target_step" not in first

    legacy: dict[str, object] = {}
    CLI._bind_training_stop_contract(
        legacy,
        schema_version="twirl_fm0_1_poc_config_v1",
        optimizer_horizon=20_000,
        invocation_target_step=7,
    )
    assert legacy == {"target_step": 7}


def test_cli_preserves_fm02_invariant_canary_stop_contract() -> None:
    contract: dict[str, object] = {}
    CLI._bind_training_stop_contract(
        contract,
        schema_version=CLI.FM0_2_CONFIG_SCHEMA_VERSION,
        optimizer_horizon=20_000,
        invocation_target_step=500,
        canary_payload={
            "immutable_milestone_steps": [0, 64, 500, 1000, 2000],
            "authorized_stop_after_steps": [64, 500, 1000, 2000],
            "fp32_synthetic_smoke_steps": 8,
        },
    )
    assert contract == {
        "training_horizon_step": 20_000,
        "immutable_milestone_steps": [0, 64, 500, 1000, 2000],
        "authorized_stop_after_steps": [64, 500, 1000, 2000],
        "synthetic_smoke_step": 8,
    }
    assert "target_step" not in contract


@pytest.mark.parametrize(
    ("synthetic", "target_step", "precision"),
    [(True, 8, "fp32"), (False, 64, "bf16"), (False, 2000, "bf16")],
)
def test_cli_accepts_only_frozen_fm03_cuda_canary_invocations(
    synthetic: bool,
    target_step: int,
    precision: str,
) -> None:
    canary = CLI._validate_fm0_3_execution_contract(
        _fm03_config(),
        variant="TWIRL-FM0.3.2",
        target_step=target_step,
        synthetic_smoke=synthetic,
        precision=precision,
        device="cuda:0",
        optimizer_horizon=20_000,
    )
    assert canary["one_h200_per_variant"] is True


def test_cli_enforces_fm03_fresh_smoke_and_step64(tmp_path: Path) -> None:
    checkpoint = tmp_path / "checkpoint.pt"
    checkpoint.write_bytes(b"not-used")
    for synthetic, target_step in ((True, 8), (False, 64)):
        assert (
            CLI._validate_fm0_3_resume_invocation(
                target_step=target_step,
                synthetic_smoke=synthetic,
                output_dir=tmp_path / "fresh",
                resume_checkpoint=None,
            )
            is None
        )
        with pytest.raises(ValueError, match="must start fresh"):
            CLI._validate_fm0_3_resume_invocation(
                target_step=target_step,
                synthetic_smoke=synthetic,
                output_dir=tmp_path / "fresh",
                resume_checkpoint=checkpoint,
            )


def test_cli_requires_exact_fm03_step64_path_for_step2000(tmp_path: Path) -> None:
    output = tmp_path / "run"
    output.mkdir()
    expected = output / "checkpoint_step_00000064.pt"
    expected.write_bytes(b"checkpoint")
    assert CLI._validate_fm0_3_resume_invocation(
        target_step=2_000,
        synthetic_smoke=False,
        output_dir=output,
        resume_checkpoint=expected,
    ) == 64

    with pytest.raises(ValueError, match="requires the step-64"):
        CLI._validate_fm0_3_resume_invocation(
            target_step=2_000,
            synthetic_smoke=False,
            output_dir=output,
            resume_checkpoint=None,
        )
    other = tmp_path / "checkpoint_step_00000064.pt"
    other.write_bytes(b"wrong run")
    with pytest.raises(ValueError, match="this run's exact step-64"):
        CLI._validate_fm0_3_resume_invocation(
            target_step=2_000,
            synthetic_smoke=False,
            output_dir=output,
            resume_checkpoint=other,
        )


def test_cli_requires_exact_git_sha_for_fm03() -> None:
    expected = "a" * 40
    assert CLI._require_fm0_3_expected_git_sha(expected) == expected
    for invalid in (None, "a" * 39, "A" * 40, "g" * 40):
        with pytest.raises(ValueError, match="FM0.3"):
            CLI._require_fm0_3_expected_git_sha(invalid)


def test_cli_requires_same_variant_smoke_for_real_fm03(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    common = {
        "is_fm0_3": True,
        "synthetic_smoke": False,
        "campaign_id": "twirl_fm0_3_matched_canary_test",
        "variant": "TWIRL-FM0.3.2",
        "architecture": "conformer",
        "expected_git_sha": "a" * 40,
        "authorities": {
            "design_sha256": "b" * 64,
            "config_sha256": "c" * 64,
            "freeze_receipt_sha256": "d" * 64,
        },
    }
    with pytest.raises(ValueError, match="same-variant 8-step smoke"):
        CLI._build_fm0_3_prestart_smoke_binding(
            smoke_run=None,
            smoke_summary_sha256=None,
            **common,
        )

    expected = {"kind": "fm0_3_same_variant_prestart_smoke"}
    calls: list[dict[str, object]] = []

    def fake_validate(root: Path, **kwargs: object) -> dict[str, object]:
        calls.append({"root": root, **kwargs})
        return expected

    monkeypatch.setattr(CLI, "validate_fm0_3_prestart_smoke", fake_validate)
    smoke = tmp_path / "smoke"
    result = CLI._build_fm0_3_prestart_smoke_binding(
        smoke_run=smoke,
        smoke_summary_sha256="e" * 64,
        **common,
    )

    assert result == expected
    assert calls == [
        {
            "root": smoke,
            "expected_summary_sha256": "e" * 64,
            "expected_campaign_id": "twirl_fm0_3_matched_canary_test",
            "expected_variant": "TWIRL-FM0.3.2",
            "expected_architecture": "conformer",
            "expected_git_sha": "a" * 40,
            "expected_authorities": common["authorities"],
        }
    ]


def test_cli_smoke_does_not_accept_a_prestart_smoke(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="cannot depend on another smoke"):
        CLI._build_fm0_3_prestart_smoke_binding(
            is_fm0_3=True,
            synthetic_smoke=True,
            smoke_run=tmp_path / "smoke",
            smoke_summary_sha256="e" * 64,
            campaign_id="twirl_fm0_3_matched_canary_test",
            variant="TWIRL-FM0.3.1",
            architecture="tcn",
            expected_git_sha="a" * 40,
            authorities={},
        )


@pytest.mark.parametrize(
    ("drift", "message"),
    [
        ("status", "frozen matched canary"),
        ("authorization", "frozen matched canary"),
        ("variant", "matched FM0.3.1/0.3.2"),
        ("architecture", "architecture contract"),
    ],
)
def test_cli_rejects_fm03_authority_or_variant_drift(
    drift: str,
    message: str,
) -> None:
    config = _fm03_config()
    variant = "TWIRL-FM0.3.1"
    if drift == "status":
        config["status"] = "draft"
    elif drift == "authorization":
        config["authorization"]["gpu_training_authorized"] = False
    elif drift == "variant":
        variant = "TWIRL-FM0.3.3"
    else:
        config["variants"][variant]["architecture"] = "conformer"
    with pytest.raises(ValueError, match=message):
        CLI._validate_fm0_3_execution_contract(
            config,
            variant=variant,
            target_step=64,
            synthetic_smoke=False,
            precision="bf16",
            device="cuda:0",
            optimizer_horizon=20_000,
        )


@pytest.mark.parametrize(
    ("synthetic", "target_step", "precision", "device", "message"),
    [
        (True, 7, "fp32", "cuda", "stop at step 8"),
        (True, 8, "bf16", "cuda", "must use FP32"),
        (False, 65, "bf16", "cuda", "not preauthorized"),
        (False, 64, "fp32", "cuda", "must use BF16"),
        (False, 64, "bf16", "cpu", "require one CUDA"),
    ],
)
def test_cli_rejects_fm03_invocation_drift(
    synthetic: bool,
    target_step: int,
    precision: str,
    device: str,
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        CLI._validate_fm0_3_execution_contract(
            _fm03_config(),
            variant="TWIRL-FM0.3.1",
            target_step=target_step,
            synthetic_smoke=synthetic,
            precision=precision,
            device=device,
            optimizer_horizon=20_000,
        )


@pytest.mark.parametrize(
    ("section", "field", "value", "message"),
    [
        ("data_contract", "gradient_role", "temporal_holdout", "training role"),
        ("data_contract", "sealed_rows_allowed", 1, "training role"),
        ("cadence_contract", "context_length_cadences", 2_048, "geometry"),
        ("cadence_contract", "cadence_averaging", True, "geometry"),
    ],
)
def test_cli_rejects_fm03_role_or_cadence_contract_drift(
    section: str,
    field: str,
    value: object,
    message: str,
) -> None:
    config = _fm03_config()
    config[section][field] = value
    with pytest.raises(ValueError, match=message):
        CLI._validate_fm0_3_execution_contract(
            config,
            variant="TWIRL-FM0.3.1",
            target_step=64,
            synthetic_smoke=False,
            precision="bf16",
            device="cuda:0",
            optimizer_horizon=20_000,
        )


def test_cli_builds_only_the_config_bound_fm03_train_role(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    root = tmp_path / "composite"
    root.mkdir()
    source = root / "source_bindings.csv"
    role = root / "role_index.csv"
    source.write_text("source\n", encoding="utf-8")
    role.write_text("role\n", encoding="utf-8")
    receipt = root / "receipt.json"
    receipt.write_text(
        json.dumps(
            {
                "schema_version": CLI.COMPOSITE_RELEASE_SCHEMA_VERSION,
                "passed": True,
                "limits": {
                    "scientific_training_eligible": True,
                    "sealed_rows_selected": 0,
                    "sealed_test_access_authorized": False,
                },
            }
        )
        + "\n",
        encoding="utf-8",
    )

    def digest(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    config = _fm03_config()
    config["composite_input"] = {
        "receipt_sha256": digest(receipt),
        "source_bindings_sha256": digest(source),
        "role_index_sha256": digest(role),
    }
    captured: dict[str, object] = {}

    class FakeDataset:
        def __init__(self, dataset_config: object) -> None:
            captured["config"] = dataset_config
            self.contract = {
                "n_sources": 3,
                "n_observations": 5,
                "n_excluded_missing_required_views": 1,
            }

    monkeypatch.setattr(CLI, "FM0CompositeReleaseDataset", FakeDataset)
    dataset, binding = CLI._build_fm0_3_composite_dataset(
        config=config,
        release_root=root,
        receipt_path=receipt,
        receipt_sha256=digest(receipt),
        variant="TWIRL-FM0.3.1",
        seed=560067,
        windows_per_epoch=1_280_000,
    )

    dataset_config = captured["config"]
    assert dataset.contract["n_observations"] == 5
    assert dataset_config.role == CLI.FM0_3_TRAIN_ROLE
    assert dataset_config.window_length == FM0_3_WINDOW_LENGTH
    assert dataset_config.mask_target_fraction == FM0_3_MASK_TARGET_FRACTION
    assert dataset_config.mask_span_range == FM0_3_MASK_SPAN_RANGE
    assert binding["kind"] == "fm0_3_composite_release"
    assert binding["role"] == CLI.FM0_3_TRAIN_ROLE
    assert binding["receipt_sha256"] == config["composite_input"]["receipt_sha256"]

    config["composite_input"]["role_index_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="role-index file differs"):
        CLI._build_fm0_3_composite_dataset(
            config=config,
            release_root=root,
            receipt_path=receipt,
            receipt_sha256=digest(receipt),
            variant="TWIRL-FM0.3.1",
            seed=560067,
            windows_per_epoch=1_280_000,
        )


def test_cli_binds_exact_shared_fm03_payload_evaluation_plan(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    root = tmp_path / "payload-plan"
    receipt_path = root / "receipt.json"
    schedule_path = root / "schedule.csv"
    config = _fm03_config()
    config["composite_input"] = {
        "receipt_sha256": "a" * 64,
        "source_bindings_sha256": "b" * 64,
        "role_index_sha256": "c" * 64,
    }
    temporal_authority = {
        "receipt_sha256": CLI.FM0_3_TEMPORAL_PANEL_RECEIPT_SHA256,
        "panel_sha256": "8" * 64,
        "sector_bindings_sha256": "9" * 64,
    }
    receipt = {
        "producer_git_sha": "d" * 40,
        "identity_plan": {
            "receipt_sha256": "e" * 64,
            "schedule_sha256": "f" * 64,
            "producer_git_sha": "d" * 40,
            "input_authorities": {
                "temporal_panel": temporal_authority,
                "composite_release": dict(config["composite_input"]),
            },
        },
        "payload_bindings": {
            "source_shard_bindings_sha256": "1" * 64,
            "crop_payload_bindings_sha256": "2" * 64,
            "n_crops_frozen": 1_440,
        },
    }
    result = SimpleNamespace(
        root=root.resolve(),
        receipt_path=receipt_path.resolve(),
        receipt_sha256="3" * 64,
        schedule_path=schedule_path.resolve(),
        schedule_sha256="4" * 64,
        receipt=receipt,
    )
    calls: list[dict[str, object]] = []

    def fake_validate(plan_root: Path, **kwargs: object) -> object:
        calls.append({"root": plan_root, **kwargs})
        return result

    monkeypatch.setattr(CLI, "validate_matched_canary_payload_plan", fake_validate)
    binding = CLI._build_fm0_3_evaluation_plan_binding(
        config=config,
        plan_root=root,
        receipt_sha256="3" * 64,
        expected_git_sha="d" * 40,
    )

    assert calls == [
        {
            "root": root,
            "expected_receipt_sha256": "3" * 64,
            "require_read_only": True,
        }
    ]
    assert binding["kind"] == "fm0_3_payload_screened_evaluation_plan"
    assert binding["receipt_sha256"] == "3" * 64
    assert binding["schedule_sha256"] == "4" * 64
    assert binding["identity_plan_producer_git_sha"] == "d" * 40
    assert binding["temporal_panel_receipt_sha256"] == (
        CLI.FM0_3_TEMPORAL_PANEL_RECEIPT_SHA256
    )
    assert binding["temporal_panel_sha256"] == "8" * 64
    assert binding["temporal_panel_sector_bindings_sha256"] == "9" * 64
    assert binding["n_crops"] == 1_440

    receipt["producer_git_sha"] = "9" * 40
    with pytest.raises(ValueError, match="exact code"):
        CLI._build_fm0_3_evaluation_plan_binding(
            config=config,
            plan_root=root,
            receipt_sha256="3" * 64,
            expected_git_sha="d" * 40,
        )

    receipt["producer_git_sha"] = "d" * 40
    temporal_authority["receipt_sha256"] = "7" * 64
    with pytest.raises(ValueError, match="temporal-panel receipt differs"):
        CLI._build_fm0_3_evaluation_plan_binding(
            config=config,
            plan_root=root,
            receipt_sha256="3" * 64,
            expected_git_sha="d" * 40,
        )
