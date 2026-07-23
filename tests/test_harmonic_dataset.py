from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from twirl.lightcurves.a2v1_cadence_reference import (
    AUTHORITY_EXCLUSION_EXTERNAL_BIT,
    AUTHORITY_EXCLUSION_POLICY_CONTRACT,
)
from twirl.lightcurves.external_quality import (
    EFFECTIVE_QUALITY_POLICY,
    EXTERNAL_QUALITY_POLICY_CONTRACT,
)
from twirl.vetting.adjudication_audit import HARMONIC_CNN_TARGET_POLICY
from twirl.vetting.harmonic_dataset import (
    HarmonicNativeDataset,
    _cap_injected_task_weights,
    attach_fold_training_weights,
    build_injection_pretraining_rows,
    build_metadata_matrix,
    candidate_bls_ephemeris,
    native_group_path,
    prepare_harmonic_training_rows,
)
from twirl.vetting.harmonic_inputs import (
    CHANNEL_CONTRACT,
    NATIVE_DATASETS,
    RAW_PAIR_CONTRACT_VERSION,
)


def _write_native_target(path: Path) -> None:
    n_cadences = 48
    with h5py.File(path, "w") as h5:
        h5.attrs["contract_version"] = RAW_PAIR_CONTRACT_VERSION
        h5.attrs["time_system"] = "BJD"
        for name, channels in CHANNEL_CONTRACT.items():
            h5.attrs[name] = json.dumps(channels)
        group = h5.create_group("targets/0000000000000001")
        payload = {
            "time": 2459825.0 + np.arange(n_cadences) * 200.0 / 86400.0,
            "cadenceno": np.arange(n_cadences),
            "orbitid": np.ones(n_cadences, dtype=np.int32),
            "quality": np.zeros(n_cadences, dtype=np.int32),
            "raw_flux_small": np.full(n_cadences, 100.0),
            "raw_flux_err_small": np.full(n_cadences, 2.0),
            "raw_flux_primary": np.full(n_cadences, 300.0),
            "raw_flux_err_primary": np.full(n_cadences, 3.0),
            "det_flux_adp_sml": np.ones(n_cadences),
            "det_flux_adp": np.ones(n_cadences),
        }
        for name in NATIVE_DATASETS:
            group.create_dataset(name, data=payload[name])
        h5.attrs["external_quality_policy_contract"] = (
            EXTERNAL_QUALITY_POLICY_CONTRACT
        )
        h5.attrs["effective_quality_policy"] = EFFECTIVE_QUALITY_POLICY
        h5.attrs["cadence_reference_contract_version"] = (
            "s56_a2v1_cadence_reference_v1"
        )
        h5.attrs["cadence_reference_cadence_authority"] = "qlp_cam_quat"
        h5.attrs["cadence_reference_quality_authority"] = (
            "spoc_and_qlp_quality_flags"
        )
        h5.attrs["cadence_reference_table"] = "/authority/reference.csv"
        h5.attrs["cadence_reference_manifest"] = "/authority/reference.json"
        h5.attrs["cadence_reference_table_sha256"] = "1" * 64
        h5.attrs["cadence_reference_manifest_sha256"] = "2" * 64
        h5.attrs["cadence_reference_source_declaration_sha256"] = "3" * 64
        h5.attrs["authority_exclusion_policy_contract"] = (
            AUTHORITY_EXCLUSION_POLICY_CONTRACT
        )
        h5.attrs["authority_exclusion_external_bit"] = (
            AUTHORITY_EXCLUSION_EXTERNAL_BIT
        )
        h5.attrs["authority_exclusions_sha256"] = "4" * 64
        h5.attrs["n_authority_exclusions"] = 0
        quality_counts = {
            "n_cad_total": n_cadences,
            "n_cad_internal_bad": 0,
            "n_cad_external_bad": 0,
            "n_cad_external_only_bad": 0,
            "n_cad_authority_excluded": 0,
            "n_cad_effective_bad": 0,
        }
        group.attrs["quality_policy_contract"] = EXTERNAL_QUALITY_POLICY_CONTRACT
        for name, value in quality_counts.items():
            group.attrs[name] = value
            h5.attrs[f"quality_overlay_{name}"] = value


def test_native_dataset_reuses_one_hdf5_handle(tmp_path: Path) -> None:
    path = tmp_path / "native.h5"
    _write_native_target(path)
    rows = pd.DataFrame(
        {
            "review_id": ["one", "two"],
            "tic": [1, 1],
            "native_group_path": ["targets/0000000000000001"] * 2,
            "period_d": [0.2, 0.3],
            "t0_bjd": [2459825.01, 2459825.01],
            "duration_min": [10.0, 10.0],
            "morphology_target_index": [3, 3],
            "preserve_target_index": [0, 0],
            "harmonic_target_index": [-1, -1],
            "morphology_weight": [1.0, 1.0],
            "preserve_weight": [1.0, 1.0],
            "harmonic_weight": [0.0, 0.0],
        }
    )
    dataset = HarmonicNativeDataset(
        rows,
        native_h5=path,
        metadata=np.empty((2, 0), dtype=np.float32),
        cache_size=0,
        profile="seven_harmonic_shape",
    )

    dataset[0]
    first_handle = dataset._native_handles[str(path)]
    dataset[1]

    assert dataset._native_handles[str(path)] is first_handle
    assert first_handle.id.valid
    dataset.close()
    assert not first_handle.id.valid


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
