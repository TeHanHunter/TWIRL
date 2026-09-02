from __future__ import annotations

import hashlib
import importlib.util
import inspect
import json
from pathlib import Path

import numpy as np
import pytest
import yaml

from twirl.models.fm0 import event_transfer_candidate_band as band
from twirl.models.fm0.event_transfer_canary import ProbeFit
from twirl.models.fm0.event_transfer_mechanics import (
    EVENT_TRANSFER_MECHANICS_CAMPAIGN_ID,
    EVENT_TRANSFER_MECHANICS_RESULT_SCHEMA_VERSION,
)

ROOT = Path(__file__).parents[1]
CONFIG = (
    ROOT / "configs/models/twirl_fm0_2_s66_s77_event_transfer_candidate_band_v4.yaml"
)
DRIVER = (
    ROOT
    / "scripts/stage5_validation/evaluate_twirl_fm0_event_transfer_candidate_band.py"
)


def _v3_result() -> dict[str, object]:
    return {
        "schema_version": EVENT_TRANSFER_MECHANICS_RESULT_SCHEMA_VERSION,
        "campaign_id": EVENT_TRANSFER_MECHANICS_CAMPAIGN_ID,
        "evaluation_completed": True,
        "passed": True,
        "probe_mechanics": {"passed": True},
        "cadence_preservation": {
            "nominal_cadence_seconds": 200,
            "one_logit_per_native_cadence": True,
            "patching_or_cadence_averaging": False,
            "representation_pooling": False,
            "held_out_reduction": (
                "exact_injected_support_center_vs_exact_paired_clean_cadence"
            ),
            "off_support_used_for_fit_or_evaluation": False,
            "support_used_as_input_feature": False,
        },
        "boundaries": {
            "fm_checkpoint_loaded": False,
            "fm_encoder_trained_or_evaluated": False,
            "fm_metrics_interpreted": False,
            "bls_or_period_features_used": False,
            "sealed_test_opened": False,
            "formal_model_gate": False,
            "development_only": True,
        },
    }


def test_frozen_config_has_only_native_128_candidate_band_arms() -> None:
    config, source, digest = band.load_event_transfer_band_config(CONFIG)
    assert source == CONFIG.resolve()
    assert len(digest) == 64
    temporal = config["temporal_resolution"]
    assert temporal["context_cadences"] == 128
    assert temporal["fixed_band_start_index_zero_based"] == 36
    assert temporal["fixed_band_stop_index_inclusive_zero_based"] == 92
    assert temporal["fixed_band_cadences"] == 57
    assert temporal["one_encoder_token_per_native_cadence"] is True
    assert temporal["score_all_128_tokens_before_reduction"] is True
    assert temporal["temporal_downsampling"] is False
    assert temporal["patching_or_cadence_averaging"] is False
    assert temporal["representation_pooling"] is False
    assert tuple(config["features"]["arms"]) == band.FEATURE_SPECS
    assert all(
        "2048" not in arm and "z_cadence" not in arm for arm in band.FEATURE_SPECS
    )
    assert config["purpose"]["architecture_selection_claim"] is False


def test_config_rejects_support_router_and_changed_band(tmp_path: Path) -> None:
    changed = yaml.safe_load(CONFIG.read_text())
    changed["probe"]["held_out_support_used_as_router"] = True
    path = tmp_path / "support-router.yaml"
    path.write_text(yaml.safe_dump(changed))
    with pytest.raises(ValueError, match="held_out_support_used_as_router"):
        band.load_event_transfer_band_config(path)

    changed = yaml.safe_load(CONFIG.read_text())
    changed["temporal_resolution"]["fixed_band_start_index_zero_based"] = 37
    path = tmp_path / "changed-band.yaml"
    path.write_text(yaml.safe_dump(changed))
    with pytest.raises(ValueError, match="cadence grid"):
        band.load_event_transfer_band_config(path)


def test_fixed_band_is_literal_and_has_no_support_argument() -> None:
    mask = band.fixed_candidate_band_mask()
    assert mask.shape == (128,)
    np.testing.assert_array_equal(np.flatnonzero(mask), np.arange(36, 93))
    assert (
        "support" not in inspect.signature(band.reduce_fixed_candidate_band).parameters
    )
    assert (
        "support"
        not in inspect.signature(band.validate_fixed_band_scoring_validity).parameters
    )

    logits = np.zeros((2, 128), dtype=np.float64)
    valid = np.ones((2, 128), dtype=bool)
    logits[:, 50] = [2.0, 3.0]
    logits[:, 0] = 100.0
    first = band.reduce_fixed_candidate_band(logits, valid)
    # Two possible per-sample truth supports cannot route the reducer because
    # no truth array enters it; both are deliberately disjoint from cadence 50.
    support_a = np.zeros_like(valid)
    support_b = np.zeros_like(valid)
    support_a[:, 40] = True
    support_b[:, 88] = True
    assert not np.array_equal(support_a, support_b)
    second = band.reduce_fixed_candidate_band(logits.copy(), valid.copy())
    np.testing.assert_array_equal(first, [2.0, 3.0])
    np.testing.assert_array_equal(first, second)
    band.validate_fixed_band_scoring_validity(valid)


def test_fit_validity_uses_only_fixed_band_and_requires_support_inside() -> None:
    valid = np.ones((4, 128), dtype=bool)
    targets = np.zeros_like(valid)
    targets[1, 40] = True
    targets[3, 88] = True
    fit_valid = band.fixed_band_training_validity(valid, targets)
    assert np.all(fit_valid[:, 36:93])
    assert not np.any(fit_valid[:, :36])
    assert not np.any(fit_valid[:, 93:])

    escaped = targets.copy()
    escaped[1, 10] = True
    with pytest.raises(ValueError, match="escaped"):
        band.fixed_band_training_validity(valid, escaped)

    missing = valid.copy()
    missing[0, 50] = False
    missing[1, 50] = False
    with pytest.raises(ValueError, match="must be valid"):
        band.fixed_band_training_validity(missing, targets)


def test_shared_linear_scorer_emits_all_tokens_before_fixed_max() -> None:
    fit = ProbeFit(
        weight=np.asarray([2.0, -1.0], dtype=np.float32),
        bias=0.5,
        center=np.zeros(2, dtype=np.float32),
        scale=np.ones(2, dtype=np.float32),
        objective_history=(),
    )
    tokens = np.zeros((2, 128, 2), dtype=np.float32)
    valid = np.ones((2, 128), dtype=bool)
    tokens[0, 50] = [3.0, 1.0]
    tokens[1, 90] = [2.0, -2.0]
    tokens[:, 0] = [100.0, 0.0]
    scores, logits = band.score_fixed_candidate_band(fit, tokens, valid)
    assert logits.shape == (2, 128)
    assert logits[0, 50] == pytest.approx(5.5)
    assert logits[1, 90] == pytest.approx(6.5)
    np.testing.assert_allclose(scores, [5.5, 6.5])
    # A token mutation changes exactly that token logit. No cadence
    # representation is pooled or averaged before the fixed scalar maximum.
    changed = tokens.copy()
    changed[0, 75, 0] = 4.0
    changed_scores, changed_logits = band.score_fixed_candidate_band(
        fit, changed, valid
    )
    difference = np.flatnonzero(changed_logits[0] != logits[0])
    np.testing.assert_array_equal(difference, [75])
    assert changed_scores[0] == pytest.approx(8.5)


def test_v3_prerequisite_is_byte_bound_and_must_pass(tmp_path: Path) -> None:
    result_path = tmp_path / "result.json"
    payload = _v3_result()
    result_path.write_text(json.dumps(payload, sort_keys=True) + "\n")
    digest = hashlib.sha256(result_path.read_bytes()).hexdigest()
    verified = band.validate_mechanics_prerequisite(
        result_path,
        expected_sha256=digest,
    )
    assert verified["passed"] is True
    assert verified["sha256"] == digest

    with pytest.raises(ValueError, match="SHA-256 differs"):
        band.validate_mechanics_prerequisite(
            result_path,
            expected_sha256="0" * 64,
        )
    payload["passed"] = False
    result_path.write_text(json.dumps(payload, sort_keys=True) + "\n")
    digest = hashlib.sha256(result_path.read_bytes()).hexdigest()
    with pytest.raises(ValueError, match="did not pass"):
        band.validate_mechanics_prerequisite(
            result_path,
            expected_sha256=digest,
        )


def test_driver_requires_mechanics_result_and_has_no_architecture_choice() -> None:
    spec = importlib.util.spec_from_file_location(
        "event_transfer_candidate_band_driver",
        DRIVER,
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    parser = module.build_parser()
    destinations = {action.dest for action in parser._actions}
    assert "mechanics_result" in destinations
    assert "step0_checkpoint" in destinations
    assert "step2000_checkpoint" in destinations
    assert "architecture" not in destinations
    assert "support" not in destinations
