from __future__ import annotations

import numpy as np
import pandas as pd

from twirl.vetting.adjudication_audit import HARMONIC_CNN_TARGET_POLICY
from twirl.vetting.harmonic_dataset import (
    _cap_injected_task_weights,
    attach_fold_training_weights,
    build_injection_pretraining_rows,
    build_metadata_matrix,
    candidate_bls_ephemeris,
    native_group_path,
    prepare_harmonic_training_rows,
)


def test_candidate_inputs_ignore_human_corrected_periods() -> None:
    row = {
        "period_d": 2.0,
        "t0_bjd": 2459825.25,
        "duration_min": 12.0,
        "model_period_d": 6.0,
        "effective_period_d": 6.0,
        "adjudicated_period_d": 6.0,
        "effective_period_factor": 3.0,
    }

    assert candidate_bls_ephemeris(row) == (2.0, 2459825.25, 12.0)


def test_native_group_uses_tic_or_injection_id_only() -> None:
    assert native_group_path({"tic": 42, "source_kind": "real_candidate"}) == "targets/0000000000000042"
    assert native_group_path(
        {"tic": 42, "source_kind": "injection_recovery", "injection_id": "predet_000042"}
    ) == "injections/predet_000042"


def test_training_rows_keep_broad_for_preserve_and_not_morphology() -> None:
    rows = []
    labels = [
        ("planet_like", "planet_like", True, "preserve", True, "p", True, False),
        ("eclipsing_binary_or_pceb", "eclipse_contact", True, "preserve", True, "2p", True, False),
        ("stellar_variability", "smooth_variable", True, "preserve", True, "", False, False),
        ("instrumental_or_systematic", "other", True, "reject", True, "", False, False),
        ("wide_transit_like", "", False, "preserve", True, "p_over_2", True, True),
    ]
    for index in range(25):
        (
            human,
            morphology,
            morph_include,
            preserve,
            preserve_include,
            harmonic,
            harmonic_include,
            broad,
        ) = labels[index % 5]
        rows.append(
            {
                "review_id": f"real:{index}",
                "tic": 1000 + index,
                "source_kind": "real_candidate",
                "is_injected_row": False,
                "human_label": human,
                "period_d": 1.0,
                "t0_bjd": 2459825.0,
                "duration_min": 10.0,
                "morphology_target_v1": morphology,
                "morphology_include_v1": morph_include,
                "preserve_target_v1": preserve,
                "preserve_include_v1": preserve_include,
                "harmonic_target_v1": harmonic,
                "harmonic_include_v1": harmonic_include,
                "broad_preserve_only": broad,
                "model_target_policy_version": HARMONIC_CNN_TARGET_POLICY,
            }
        )

    prepared = prepare_harmonic_training_rows(pd.DataFrame(rows))
    broad = prepared[prepared["human_label"].eq("wide_transit_like")]

    assert broad["morphology_target_index"].eq(-1).all()
    assert broad["preserve_target_index"].ge(0).all()
    assert broad["harmonic_target_index"].ge(0).all()
    assert set(prepared["fixed_split"]) == {"development", "test"}


def test_metadata_normalization_uses_fit_rows_and_excludes_targets() -> None:
    rows = pd.DataFrame(
        {
            "period_d": [1.0, 2.0, 100.0],
            "sde_max": [4.0, 6.0, 1000.0],
            "truth_period_d": [7.0, 8.0, 9.0],
            "effective_period_factor": [1.0, 2.0, 3.0],
        }
    )
    matrix, normalization = build_metadata_matrix(
        rows,
        fit_mask=np.asarray([True, True, False]),
    )

    assert normalization.columns == ("period_d", "sde_max")
    assert np.allclose(normalization.center, (1.5, 5.0))
    assert matrix.shape == (3, 2)
    assert matrix[2, 0] > 50


def test_pretraining_uses_visible_injections_pairs_and_equal_real_negatives() -> None:
    rows = []
    for index in range(6):
        rows.append(
            {
                "review_id": f"inj:{index}",
                "tic": index,
                "is_injected_row": True,
                "human_label": "planet_like",
                "morphology_target_v1": "planet_like",
                "native_group_path": f"injections/inj_{index}",
            }
        )
    for index in range(20):
        rows.append(
            {
                "review_id": f"real:{index}",
                "tic": 100 + index,
                "is_injected_row": False,
                "human_label": "instrumental_or_systematic" if index % 2 else "uncertain",
                "morphology_target_v1": "other",
                "native_group_path": f"targets/{100 + index:016d}",
            }
        )

    pretrain = build_injection_pretraining_rows(pd.DataFrame(rows))

    assert len(pretrain) == 18
    assert pretrain["pretrain_target"].value_counts().to_dict() == {0: 12, 1: 6}
    assert pretrain["input_variant"].value_counts().to_dict() == {
        "injected": 6,
        "paired_original": 6,
        "observed": 6,
    }


def test_auxiliary_weights_cap_injected_nonplanet_contribution() -> None:
    rows = pd.DataFrame(
        {
            "preserve_target_index": [0] * 10 + [0] * 50 + [0] * 30,
            "is_injected_row": [False] * 10 + [True] * 50 + [True] * 30,
            "human_label": (
                ["instrumental_or_systematic"] * 10
                + ["uncertain"] * 50
                + ["planet_like"] * 30
            ),
        }
    )
    weight = _cap_injected_task_weights(
        rows,
        target_column="preserve_target_index",
        weights=np.ones(len(rows)),
    )
    real = ~rows["is_injected_row"].to_numpy()
    injected_nonplanet = rows["is_injected_row"].to_numpy() & rows["human_label"].ne("planet_like").to_numpy()
    injected_planet = rows["is_injected_row"].to_numpy() & rows["human_label"].eq("planet_like").to_numpy()

    assert weight[injected_nonplanet].sum() <= 0.25 * weight[real].sum() + 1.0e-5
    assert weight[injected_planet].sum() <= weight[real].sum() + 1.0e-5


def test_fold_weights_do_not_use_or_weight_held_out_labels() -> None:
    rows = pd.DataFrame(
        {
            "morphology_target_v1": ["planet_like", "other", "planet_like"],
            "morphology_include_v1": [True, True, True],
            "morphology_target_index": [0, 3, 0],
            "preserve_target_index": [1, 0, 1],
            "harmonic_target_index": [3, -1, 3],
            "is_injected_row": [False, False, True],
            "human_label": ["planet_like", "uncertain", "planet_like"],
        }
    )

    weighted = attach_fold_training_weights(
        rows,
        fit_mask=np.asarray([True, True, False]),
    )

    assert weighted.loc[2, ["morphology_weight", "preserve_weight", "harmonic_weight"]].eq(0).all()
    assert weighted.loc[:1, "morphology_weight"].gt(0).all()
