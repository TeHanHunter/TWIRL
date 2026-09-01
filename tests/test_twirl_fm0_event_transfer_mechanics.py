from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest
import yaml

from twirl.models.fm0 import event_transfer_canary as canary
from twirl.models.fm0 import event_transfer_mechanics as mechanics

ROOT = Path(__file__).parents[1]
CONFIG = ROOT / "configs/models/twirl_fm0_2_s66_s77_event_transfer_mechanics_v3.yaml"
DRIVER = (
    ROOT / "scripts/stage5_validation/evaluate_twirl_fm0_event_transfer_mechanics.py"
)


def _paired_cube(
    *,
    components: int = 48,
    length: int = 32,
    feature_count: int = 2,
    quality_only: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, tuple[str, ...]]:
    rng = np.random.default_rng(12)
    clean = rng.normal(0.0, 0.05, size=(components, length, feature_count))
    injected = clean.copy()
    positions = rng.integers(4, length - 4, size=components)
    if not quality_only:
        injected[np.arange(components), positions, :2] -= 1.0
    tokens = np.empty((2 * components, length, feature_count), dtype=np.float32)
    tokens[0::2] = clean
    tokens[1::2] = injected
    valid = np.ones((2 * components, length), dtype=bool)
    targets = np.zeros_like(valid)
    targets[2 * np.arange(components) + 1, positions] = True
    identities = tuple(f"component-{index // 2}" for index in range(2 * components))
    return tokens, valid, targets, identities


def test_v3_config_freezes_mechanics_only_no_merging_contract(tmp_path: Path) -> None:
    config, source, digest = mechanics.load_event_transfer_mechanics_config(CONFIG)
    assert source == CONFIG.resolve()
    assert len(digest) == 64
    assert config["temporal_resolution"]["one_logit_per_native_cadence"] is True
    assert config["temporal_resolution"]["representation_pooling"] is False
    assert config["authorization"]["frozen_checkpoint_inference"] is False
    changed = yaml.safe_load(CONFIG.read_text())
    changed["temporal_resolution"]["patching_or_cadence_averaging"] = True
    path = tmp_path / "changed.yaml"
    path.write_text(yaml.safe_dump(changed))
    with pytest.raises(ValueError, match="patching_or_cadence_averaging"):
        mechanics.load_event_transfer_mechanics_config(path)


def test_center_probe_emits_one_logit_per_cadence_and_ignores_off_support() -> None:
    tokens, valid, targets, identities = _paired_cube()
    # An identical, much deeper off-support dip defeats the legacy all-window
    # maximum but must be irrelevant to the v3 center-only fit and evaluation.
    tokens[:, 0, :2] = -100.0
    fit = mechanics.fit_shared_linear_center_probe(
        tokens,
        valid,
        targets,
        identities,
        epochs=300,
        seed=9,
    )
    legacy_scores = canary.score_shared_linear_max_probe(fit, tokens, valid)
    np.testing.assert_array_equal(legacy_scores[0::2], legacy_scores[1::2])
    logits = mechanics.score_shared_linear_cadence_logits(fit, tokens, valid)
    assert logits.shape == valid.shape
    clean, injected, selection = mechanics.oracle_support_center_pair_scores(
        logits,
        valid,
        targets,
        identities,
    )
    result = mechanics.paired_center_ranking_bootstrap(
        clean,
        injected,
        selection.components,
        replicates=100,
        seed=3,
    )
    assert result["estimate"] > 0.99
    assert np.all(fit.weight[:2] < 0.0)

    mutated = tokens.copy()
    rng = np.random.default_rng(91)
    off_support = ~targets[1::2]
    for pair_index in range(len(selection.components)):
        changed = rng.normal(
            0.0, 1.0e6, size=(np.count_nonzero(off_support[pair_index]), 2)
        )
        mutated[2 * pair_index, off_support[pair_index], :2] = changed
        mutated[2 * pair_index + 1, off_support[pair_index], :2] = -changed
    mutated_fit = mechanics.fit_shared_linear_center_probe(
        mutated,
        valid,
        targets,
        identities,
        epochs=300,
        seed=9,
    )
    np.testing.assert_array_equal(mutated_fit.weight, fit.weight)
    assert mutated_fit.bias == fit.bias
    mutated_logits = mechanics.score_shared_linear_cadence_logits(
        mutated_fit,
        mutated,
        valid,
    )
    mutated_clean, mutated_injected, _ = mechanics.oracle_support_center_pair_scores(
        mutated_logits,
        valid,
        targets,
        identities,
    )
    np.testing.assert_array_equal(mutated_clean, clean)
    np.testing.assert_array_equal(mutated_injected, injected)


def test_center_probe_ignores_noncentral_cadences_inside_duration_nine_support() -> (
    None
):
    tokens, valid, targets, identities = _paired_cube(components=12, length=32)
    centers = np.asarray(
        [int(np.flatnonzero(targets[row])[0]) for row in range(1, 24, 2)]
    )
    targets[:] = False
    for pair_index, center in enumerate(centers):
        targets[2 * pair_index + 1, center - 4 : center + 5] = True
    fit = mechanics.fit_shared_linear_center_probe(
        tokens,
        valid,
        targets,
        identities,
        epochs=100,
        seed=6,
    )
    changed = tokens.copy()
    for pair_index, center in enumerate(centers):
        support_edges = np.r_[center - 4 : center, center + 1 : center + 5]
        changed[2 * pair_index, support_edges] = 1.0e6
        changed[2 * pair_index + 1, support_edges] = -1.0e6
    changed_fit = mechanics.fit_shared_linear_center_probe(
        changed,
        valid,
        targets,
        identities,
        epochs=100,
        seed=6,
    )
    np.testing.assert_array_equal(changed_fit.center, fit.center)
    np.testing.assert_array_equal(changed_fit.scale, fit.scale)
    np.testing.assert_array_equal(changed_fit.weight, fit.weight)
    assert changed_fit.bias == fit.bias


def test_identical_quality_pairs_return_exact_chance() -> None:
    tokens, valid, targets, identities = _paired_cube(
        feature_count=8,
        quality_only=True,
    )
    fit = mechanics.fit_shared_linear_center_probe(
        tokens,
        valid,
        targets,
        identities,
        epochs=100,
        seed=4,
    )
    logits = mechanics.score_shared_linear_cadence_logits(fit, tokens, valid)
    clean, injected, selection = mechanics.oracle_support_center_pair_scores(
        logits,
        valid,
        targets,
        identities,
    )
    result = mechanics.paired_center_ranking_bootstrap(
        clean,
        injected,
        selection.components,
        replicates=100,
        seed=7,
    )
    assert result["estimate"] == 0.5
    assert result["component_pair_bootstrap_95_interval"] == [0.5, 0.5]
    assert result["tied_pairs"] == len(selection.components)


@pytest.mark.parametrize(
    ("mutation", "message"),
    (
        ("clean_support", "clean row contains"),
        ("invalid_support", "invalid injected cadence"),
        ("noncontiguous", "contiguous with duration"),
        ("unsupported_duration", "contiguous with duration"),
        ("component_mismatch", "identity or adjacency"),
        ("validity_mismatch", "validity masks differ"),
    ),
)
def test_center_pair_rejects_invalid_support_or_pair_alignment(
    mutation: str,
    message: str,
) -> None:
    _tokens, valid, targets, identities = _paired_cube(components=2, length=16)
    valid = valid.copy()
    targets = targets.copy()
    identities = list(identities)
    center = int(np.flatnonzero(targets[1])[0])
    if mutation == "clean_support":
        targets[0, center] = True
    elif mutation == "invalid_support":
        valid[1, center] = False
    elif mutation == "noncontiguous":
        targets[1, center + 2] = True
        targets[1, center + 3] = True
    elif mutation == "unsupported_duration":
        targets[1, center - 2 : center + 3] = True
    elif mutation == "component_mismatch":
        identities[1] = "wrong-component"
    elif mutation == "validity_mismatch":
        valid[0, 0] = False
    with pytest.raises(ValueError, match=message):
        mechanics.paired_support_center_indices(valid, targets, identities)


def test_center_pair_metric_weights_components_not_support_duration() -> None:
    valid = np.ones((4, 16), dtype=bool)
    targets = np.zeros_like(valid)
    targets[1, 7] = True
    targets[3, 4:13] = True
    identities = ("one", "one", "nine", "nine")
    selection = mechanics.paired_support_center_indices(valid, targets, identities)
    assert selection.center_cadences.tolist() == [7, 8]
    result = mechanics.paired_center_ranking_bootstrap(
        np.asarray([0.0, 1.0]),
        np.asarray([1.0, 1.0]),
        selection.components,
        replicates=100,
        seed=5,
    )
    # One winning component and one tied component receive one vote each.
    assert result["estimate"] == 0.75
    with pytest.raises(ValueError, match="replicates must be positive"):
        mechanics.paired_center_ranking_bootstrap(
            np.asarray([0.0, 1.0]),
            np.asarray([1.0, 1.0]),
            selection.components,
            replicates=0,
        )


def test_fast_driver_has_no_checkpoint_or_fm_arguments() -> None:
    spec = importlib.util.spec_from_file_location(
        "event_transfer_mechanics_driver", DRIVER
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    destinations = {action.dest for action in module.build_parser()._actions}
    assert "temporal_panel_dir" in destinations
    assert "step0_checkpoint" not in destinations
    assert "step2000_checkpoint" not in destinations
    assert "run_dir" not in destinations
