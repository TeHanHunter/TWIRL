from __future__ import annotations

import json
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest

import twirl.vetting.teacher_v2_training as teacher_v2_training
from twirl.vetting.harmonic_cnn import HarmonicModelConfig
from twirl.vetting.teacher_v2 import (
    TEACHER_V2_ROLE_POLICY,
    assert_teacher_v2_feature_columns,
    assign_s56_injection_roles,
    attach_prior_human_splits,
    build_global_tic_split_registry,
    build_franklin_a2v1_rereview_queue,
    build_s56_injection_training_rows,
    build_s57_external_role_manifest,
    build_s56_real_workload_candidates,
    freeze_real_tic_workload_threshold,
    injection_recall_at_fold_thresholds,
    join_franklin_labels_to_queue,
    mark_native_input_availability,
    normalize_franklin_labels,
    normalize_real_adp_candidates,
    transfer_human_labels_to_a2v1_candidates,
)
from twirl.vetting.teacher_v2_cnn import (
    TeacherV2TrainConfig,
    build_teacher_v2_cnn,
    teacher_v2_multitask_loss,
)
from twirl.vetting.teacher_v2_training import (
    _validate_training_rows,
    attach_teacher_v2_fold_weights,
    build_teacher_v2_metadata_matrix,
    compact_metrics,
    write_teacher_v2_native_verification_cache,
)
from twirl.vetting.teacher_v2_inference import prepare_teacher_v2_inference_rows
from twirl.vetting.teacher_v2_evaluation import evaluate_locked_human_holdout
from twirl.vetting.teacher_v2_selection import freeze_teacher_v2_selection
from twirl.vetting.teacher_v2_recovery import (
    aggregate_compact_recovery,
    bls_topk_recovery_table,
    period_radius_tmag_support,
)


def _full_schedule() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for cell in range(2500):
        for slot in range(8):
            index = 8 * cell + slot
            rows.append(
                {
                    "injection_id": f"s56eval_{index:06d}",
                    "tic": 1_000_000 + index,
                    "grid_cell_id": f"p{cell // 50:02d}_r{cell % 50:02d}",
                    "grid_slot": slot,
                }
            )
    return pd.DataFrame(rows)


def test_s56_injection_roles_are_exactly_six_two_per_cell() -> None:
    roles = assign_s56_injection_roles(_full_schedule())
    counts = pd.crosstab(roles["grid_cell_id"], roles["teacher_v2_partition"])
    assert len(roles) == 20_000
    assert counts["development"].eq(6).all()
    assert counts["locked_holdout"].eq(2).all()
    assert roles["tic"].nunique() == 20_000


def test_injection_rows_use_only_recovered_hosts_and_keep_paired_negatives() -> None:
    schedule = _full_schedule()
    injection_a = schedule.loc[0, "injection_id"]
    injection_b = schedule.loc[1, "injection_id"]
    candidates = pd.DataFrame(
        [
            {
                "injection_id": injection_a,
                "review_id": f"a-{rank}",
                "tic": schedule.loc[0, "tic"],
                "period_d": 1.0 + rank / 10,
                "t0_bjd": 2_459_000.0,
                "duration_min": 10.0,
                "rep_peak_rank": rank,
                "is_injected_signal_peak": match,
                "nearest_harmonic_factor": factor,
                "is_injected_row": True,
            }
            for rank, match, factor in ((1, True, 1.0), (2, False, 1.0), (3, True, 0.5))
        ]
        + [
            {
                "injection_id": injection_b,
                "review_id": f"b-{rank}",
                "tic": schedule.loc[1, "tic"],
                "period_d": 2.0 + rank / 10,
                "t0_bjd": 2_459_000.0,
                "duration_min": 10.0,
                "rep_peak_rank": rank,
                "is_injected_signal_peak": False,
                "nearest_harmonic_factor": 1.0,
                "is_injected_row": True,
            }
            for rank in (1, 2)
        ]
    )
    rows, manifest = build_s56_injection_training_rows(
        schedule,
        candidates,
        native_h5=Path("native.h5"),
    )
    assert set(rows["injection_id"]) == {injection_a}
    assert rows["teacher_v2_role"].value_counts().to_dict() == {
        "compact_positive": 2,
        "paired_pre_injection": 2,
        "same_lc_unmatched_peak": 1,
    }
    assert rows.loc[rows["teacher_v2_role"].eq("compact_positive"), "compact_target_index"].eq(1).all()
    assert rows.loc[rows["teacher_v2_role"].ne("compact_positive"), "compact_target_index"].eq(0).all()
    assert rows.loc[rows["teacher_v2_role"].eq("compact_positive"), "harmonic_target_index"].tolist() == [3, 2]
    recovered = manifest.set_index("injection_id")["bls_top5_recovered"]
    assert bool(recovered[injection_a])
    assert not bool(recovered[injection_b])


def test_franklin_identity_ignores_legacy_row_id() -> None:
    labels = pd.DataFrame(
        {
            "row_id": [4, 4],
            "candidate_key": ["real:a", "real:b"],
            "tic": [10, 11],
            "sector": [56, 56],
            "label": ["uncertain", "stellar_variability"],
        }
    )
    normalized = normalize_franklin_labels(labels)
    assert normalized["legacy_row_id"].tolist() == [4, 4]
    assert normalized["franklin_identity"].is_unique
    assert not normalized["teacher_v2_training_include"].any()


def test_franklin_candidate_key_join_and_blinded_rereview_queue() -> None:
    queue = pd.DataFrame(
        {
            "review_id": ["real:10", "real:11"],
            "tic": ["10", "11"],
            "sector": ["56", "56"],
            "period_d": ["1.25", "2.5"],
            "t0_bjd": ["2459000.25", "2459001.5"],
            "duration_min": ["10", "12"],
        }
    )
    labels = pd.DataFrame(
        {
            "row_id": [99, 99],
            "candidate_key": [
                "real:10|10|56|1.25|2459000.25",
                "real:11|11|56|2.5|2459001.5",
            ],
            "tic": [10, 11],
            "sector": [56, 56],
            "label": ["stellar_variability", "uncertain"],
            "label_source": ["human", "human"],
            "labeler": ["franklin", "franklin"],
            "notes": ["", ""],
            "updated_utc": ["2026-01-01", "2026-01-01"],
        }
    )
    joined = join_franklin_labels_to_queue(labels, queue)
    assert joined["franklin_queue_review_id"].tolist() == ["real:10", "real:11"]
    assert joined["legacy_row_id"].tolist() == [99, 99]
    candidates = pd.DataFrame(
        {
            "tic": [10, 11],
            "sector": [56, 56],
            "review_id": ["a2-10", "a2-11"],
            "rep_peak_rank": [1, 1],
            "period_d": [1.2, 2.4],
            "t0_bjd": [2_459_000.2, 2_459_001.4],
            "duration_min": [10.0, 12.0],
            "source_kind": ["real_candidate", "real_candidate"],
        }
    )
    public, private, summary = build_franklin_a2v1_rereview_queue(
        joined,
        candidates,
        n_controls=1,
        expected_signal_rows=1,
    )
    assert len(public) == 2
    assert len(private) == 2
    assert summary["n_signal_rows"] == 1
    assert "label" not in public
    assert not any("franklin" in column.lower() for column in public)


def test_s57_manifest_excludes_every_registered_s56_or_franklin_tic() -> None:
    injection_roles = pd.DataFrame(
        {"tic": [10, 11], "teacher_v2_partition": ["development", "locked_holdout"]}
    )
    human = pd.DataFrame({"tic": [12]})
    franklin = pd.DataFrame({"tic": [13]})
    s57 = pd.DataFrame(
        {
            "tic": [10, 12, 13, 14, 15],
            "grid_cell_id": ["a", "a", "b", "a", "b"],
        }
    )
    registry = build_global_tic_split_registry(
        s56_injection_roles=injection_roles,
        s56_human_rows=human,
        franklin_rows=franklin,
        s57_schedule=s57,
    )
    selected, summary = build_s57_external_role_manifest(
        s57, registry, minimum_injections=2
    )
    assert selected["tic"].tolist() == [14, 15]
    assert summary["passed"]
    assert summary["cell_support_min"] == 1


def test_real_workload_pool_excludes_every_registered_host() -> None:
    candidates = pd.DataFrame(
        {
            "tic": [1, 1, 2, 3, 4],
            "review_id": ["1a", "1b", "2a", "3a", "4a"],
        }
    )
    registry = pd.DataFrame({"tic": [1, 3]})
    candidates = mark_native_input_availability(
        candidates,
        available_tics={1, 2, 3},
    )
    selected, summary = build_s56_real_workload_candidates(
        candidates, registry, minimum_tics=1
    )
    assert set(selected["tic"]) == {2}
    assert summary["n_workload_tics"] == 1
    assert summary["n_native_unavailable_tics"] == 1
    assert summary["passed"]


def test_native_availability_preserves_upstream_exclusions() -> None:
    candidates = pd.DataFrame(
        {
            "tic": [1, 2, 3],
            "native_input_include": [True, True, False],
        }
    )
    marked = mark_native_input_availability(candidates, available_tics={1, 3})
    assert marked["native_input_include"].tolist() == [True, False, False]
    assert marked["native_input_status"].tolist() == [
        "available",
        "missing_raw_source",
        "excluded_upstream",
    ]


def test_tmag_support_audit_flags_sparse_bright_cells_without_host_reuse() -> None:
    manifest = pd.DataFrame(
        {
            "injection_id": ["a", "b", "c", "d"],
            "tic": [1, 2, 3, 4],
            "tmag": [16.5, 17.5, 19.5, 19.7],
            "grid_cell_id": ["cell0", "cell1", "cell2", "cell3"],
        }
    )

    support = period_radius_tmag_support(manifest).set_index("tmag_bin")

    assert support.loc["Tmag < 17", "n_unique_hosts"] == 1
    assert support.loc["Tmag < 17", "cell_coverage_fraction"] == 0.25
    assert bool(support.loc["Tmag < 17", "bright_enrichment_recommended"])
    assert support.loc["Tmag >= 19", "host_reuse_factor"] == 1.0
    assert not bool(
        support.loc["Tmag >= 19", "bright_enrichment_recommended"]
    )


def test_workload_threshold_never_exceeds_budget_when_scores_tie() -> None:
    scores = pd.DataFrame(
        {
            "tic": np.arange(100),
            "p_compact_transit": [0.9] * 6 + [0.1] * 94,
        }
    )
    frozen = freeze_real_tic_workload_threshold(scores, max_fraction=0.05)
    assert frozen.max_pass_tics == 5
    assert frozen.n_pass_tics == 0
    assert frozen.threshold > 0.9


def test_oof_injection_recall_uses_each_rows_fold_threshold() -> None:
    scores = pd.DataFrame(
        {
            "injection_id": ["a", "b"],
            "teacher_v2_partition": ["development", "development"],
            "teacher_v2_role": ["compact_positive", "compact_positive"],
            "is_injected_signal_peak": [True, True],
            "p_compact_transit": [0.8, 0.8],
            "cv_fold": [0, 1],
        }
    )
    result = injection_recall_at_fold_thresholds(
        scores,
        thresholds={0: 0.7, 1: 0.9, 2: 0.5, 3: 0.5, 4: 0.5},
        partition="development",
    )
    assert result["n_injections"] == 2
    assert result["n_recovered"] == 1
    assert result["recovery_fraction"] == pytest.approx(0.5)


def test_teacher_v2_leakage_guard_rejects_truth_and_roles() -> None:
    assert_teacher_v2_feature_columns(["period_d", "adp_sml_sde"])
    with pytest.raises(ValueError, match="leakage"):
        assert_teacher_v2_feature_columns(["period_d", "truth_radius_rearth"])
    with pytest.raises(ValueError, match="leakage"):
        assert_teacher_v2_feature_columns(["period_d", "teacher_v2_role"])


def test_real_adp_normalization_and_human_transfer_use_current_ephemeris() -> None:
    peaks = pd.DataFrame(
        [
            {
                "tic": 42,
                "sector": 56,
                "aperture": aperture,
                "status": "ok",
                "peak_rank": rank,
                "period_d": period,
                "t0_bjd": t0,
                "duration_min": 12.0,
                "depth": 0.1,
                "depth_snr": 15.0,
                "sde": 20.0,
                "log_power": 3.0,
            }
            for aperture, rank, period, t0 in (
                ("DET_FLUX_ADP_SML", 1, 1.0, 2_459_000.0),
                ("DET_FLUX_ADP_SML", 2, 2.0, 2_459_000.0),
                ("DET_FLUX_ADP", 1, 1.01, 2_459_000.0),
            )
        ]
    )
    candidates = normalize_real_adp_candidates(peaks, small_peaks_per_tic=2)
    human = pd.DataFrame(
        {
            "source_uid": ["old:42"],
            "review_id": ["old-review"],
            "tic": [42],
            "period_d": [1.0],
            "effective_period_d": [2.0],
            "t0_bjd": [2_459_000.0],
            "duration_min": [12.0],
            "human_label": ["planet_like"],
            "morphology_target_v1": ["planet_like"],
            "morphology_include_v1": [True],
            "preserve_target_v1": ["preserve"],
            "preserve_include_v1": [True],
            "harmonic_target_v1": ["2p"],
            "harmonic_include_v1": [True],
            "broad_preserve_only": [False],
            "model_target_policy_version": ["s56_harmonic_cnn_v1"],
            "is_injected_row": [False],
        }
    )
    transferred, audit = transfer_human_labels_to_a2v1_candidates(human, candidates)
    assert len(transferred) == 1
    assert transferred.loc[0, "period_d"] == pytest.approx(1.0)
    assert audit.loc[0, "a2v1_period_factor"] == pytest.approx(2.0)
    assert transferred.loc[0, "source_review_id"] == "old-review"
    assert bool(audit.loc[0, "a2v1_transfer_ok"])


def test_teacher_v2_model_and_loss_include_compact_head() -> None:
    torch = pytest.importorskip("torch")
    config = HarmonicModelConfig(metadata_dim=4, embedding_dim=8)
    model = build_teacher_v2_cnn(config, profile="metadata_only")
    batch = {
        "harmonic_values": torch.zeros(3, 7, 7, 4),
        "harmonic_mask": torch.ones(3, 7, 7, 4, dtype=torch.bool),
        "metadata": torch.zeros(3, 4),
        "morphology_target": torch.tensor([-1, -1, -1]),
        "preserve_target": torch.tensor([-1, -1, -1]),
        "harmonic_target": torch.tensor([-1, -1, -1]),
        "morphology_weight": torch.zeros(3),
        "preserve_weight": torch.zeros(3),
        "harmonic_weight": torch.zeros(3),
        "compact_target": torch.tensor([1, 0, 0]),
        "compact_weight": torch.tensor([1.0, 0.5, 0.5]),
    }
    output = model(batch)
    loss, parts = teacher_v2_multitask_loss(
        output, batch, config=TeacherV2TrainConfig()
    )
    assert output["compact_logits"].shape == (3, 2)
    assert torch.isfinite(loss)
    assert float(parts["compact"]) > 0


def test_teacher_v2_fold_weights_balance_sources_and_mask_paired_bls() -> None:
    rows = pd.DataFrame(
        {
            "morphology_target_index": [0, 3, -1, -1, -1, -1],
            "preserve_target_index": [1, 0, -1, -1, -1, -1],
            "harmonic_target_index": [3, -1, 3, -1, -1, -1],
            "compact_target_index": [-1, -1, 1, 0, 0, 0],
            "teacher_v2_role": [
                "real_human_morphology",
                "real_human_morphology",
                "compact_positive",
                "same_lc_unmatched_peak",
                "paired_pre_injection",
                "paired_pre_injection",
            ],
            "human_label": [
                "planet_like",
                "uncertain",
                "",
                "",
                "",
                "",
            ],
        }
    )
    fit = np.ones(len(rows), dtype=bool)
    shape = attach_teacher_v2_fold_weights(
        rows, fit_mask=fit, profile="seven_harmonic_shape"
    )
    bls = attach_teacher_v2_fold_weights(
        rows, fit_mask=fit, profile="shape_plus_periodogram_bls"
    )
    paired = rows["teacher_v2_role"].eq("paired_pre_injection")
    assert shape.loc[paired, "compact_weight"].gt(0).all()
    assert bls.loc[paired, "compact_weight"].eq(0).all()
    assert bls.loc[rows["compact_target_index"].eq(1), "compact_weight"].gt(0).all()


def test_teacher_v2_metadata_is_train_normalized_and_paired_neutral() -> None:
    rows = pd.DataFrame(
        {
            "period_d": [1.0, 2.0, 3.0],
            "sde_max": [10.0, 20.0, 30.0],
            "teacher_v2_role": [
                "compact_positive",
                "same_lc_unmatched_peak",
                "paired_pre_injection",
            ],
        }
    )
    values, normalization = build_teacher_v2_metadata_matrix(
        rows, fit_mask=np.array([True, True, False])
    )
    assert normalization.columns == ("period_d", "sde_max")
    assert np.all(values[2] == 0)


def test_compact_metrics_handles_tied_scores() -> None:
    metrics = compact_metrics(
        np.array([1, 0, 1, 0]), np.array([0.8, 0.8, 0.2, 0.2])
    )
    assert metrics["n"] == 4
    assert metrics["average_precision"] == pytest.approx(5.0 / 6.0)
    assert metrics["roc_auc"] == pytest.approx(0.5)


def _teacher_v2_validation_rows(native_h5: Path) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "review_id": ["development-row"],
            "tic": [42],
            "teacher_v2_role": ["compact_positive"],
            "teacher_v2_role_policy": [TEACHER_V2_ROLE_POLICY],
            "fixed_split": ["development"],
            "cv_fold": [0],
            "native_h5_path": [str(native_h5)],
            "native_group_path": ["targets/0000000000000042"],
            "morphology_target_index": [-1],
            "preserve_target_index": [-1],
            "harmonic_target_index": [3],
            "compact_target_index": [1],
        }
    )


def _teacher_v2_fingerprint_h5(path: Path) -> None:
    with h5py.File(path, "w") as h5:
        h5.attrs["contract_version"] = teacher_v2_training.RAW_PAIR_CONTRACT_VERSION
        h5.attrs["time_system"] = "BJD"
        for name, channels in teacher_v2_training.CHANNEL_CONTRACT.items():
            h5.attrs[name] = json.dumps(channels)


def test_teacher_v2_native_verification_cache_reuses_exact_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    native_h5 = tmp_path / "native.h5"
    cache_path = tmp_path / "verification.json"
    _teacher_v2_fingerprint_h5(native_h5)
    full_verification = {
        "passed": True,
        "failures": [],
        "counts": {"targets": 1, "injections": 0},
    }
    write_teacher_v2_native_verification_cache(
        cache_path, {str(native_h5): full_verification}
    )

    def unexpected_scan(*args: object, **kwargs: object) -> dict[str, object]:
        raise AssertionError("exactly fingerprinted input should not be rescanned")

    monkeypatch.setattr(
        teacher_v2_training, "verify_raw_pair_contract", unexpected_scan
    )
    result = _validate_training_rows(
        _teacher_v2_validation_rows(native_h5),
        native_verification_cache=cache_path,
    )

    assert result["native_verification_cache"]["cache_hits"] == [
        str(native_h5.resolve())
    ]
    assert result["native_verification_cache"]["scanned"] == []


def test_teacher_v2_native_verification_cache_invalidates_changed_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    native_h5 = tmp_path / "native.h5"
    cache_path = tmp_path / "verification.json"
    _teacher_v2_fingerprint_h5(native_h5)
    full_verification = {
        "passed": True,
        "failures": [],
        "counts": {"targets": 1, "injections": 0},
    }
    write_teacher_v2_native_verification_cache(
        cache_path, {str(native_h5): full_verification}
    )
    stat = native_h5.stat()
    os.utime(native_h5, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1))
    scans: list[Path] = []

    def record_scan(path: Path, **kwargs: object) -> dict[str, object]:
        scans.append(Path(path))
        return full_verification

    monkeypatch.setattr(
        teacher_v2_training, "verify_raw_pair_contract", record_scan
    )
    result = _validate_training_rows(
        _teacher_v2_validation_rows(native_h5),
        native_verification_cache=cache_path,
    )

    assert scans == [native_h5]
    assert result["native_verification_cache"]["cache_hits"] == []
    assert result["native_verification_cache"]["scanned"] == [
        str(native_h5.resolve())
    ]


def test_teacher_v2_inference_rows_have_no_active_targets() -> None:
    candidates = pd.DataFrame(
        {
            "tic": [42],
            "period_d": [1.2],
            "t0_bjd": [2_459_000.0],
            "duration_min": [10.0],
            "review_id": ["candidate-42"],
            "source_kind": ["real_candidate"],
            "is_injected_row": [False],
        }
    )
    rows = prepare_teacher_v2_inference_rows(candidates, native_h5=Path("native.h5"))
    assert rows.loc[0, "native_group_path"] == "targets/0000000000000042"
    assert rows.loc[0, "compact_target_index"] == -1
    assert rows.loc[0, "native_h5_path"] == "native.h5"


def test_teacher_v2_selection_freezes_threshold_before_holdout(tmp_path: Path) -> None:
    profile = "seven_harmonic_shape"
    predictions = pd.DataFrame(
        {
            "review_id": ["positive", "negative", "human"],
            "tic": [1, 2, 3],
            "teacher_v2_role": [
                "compact_positive",
                "same_lc_unmatched_peak",
                "real_human_morphology",
            ],
            "teacher_v2_partition": ["development"] * 3,
            "cv_fold": [0, 1, 2],
            "injection_id": ["inj-1", "inj-2", ""],
            "is_injected_signal_peak": [True, False, False],
            "compact_target": [1, 0, -1],
            "p_compact_transit": [0.95, 0.05, 0.2],
            "morphology_target": [-1, -1, 3],
            "p_planet_like": [0.1, 0.1, 0.05],
            "p_eclipse_contact": [0.1, 0.1, 0.05],
            "p_smooth_variable": [0.1, 0.1, 0.05],
            "p_other": [0.7, 0.7, 0.85],
        }
    )
    prediction_path = tmp_path / "oof.parquet"
    predictions.to_parquet(prediction_path, index=False)
    real_scores = pd.DataFrame(
        {"tic": np.arange(20), "p_compact_transit": np.linspace(0.0, 0.9, 20)}
    )
    for fold in range(5):
        real_scores[f"member_{fold}_p_compact_transit"] = np.linspace(
            0.0, 0.9, 20
        )
    real_path = tmp_path / "real.parquet"
    real_scores.to_parquet(real_path, index=False)
    checkpoints = []
    for fold in range(5):
        path = tmp_path / f"fold-{fold}.pt"
        path.write_bytes(f"checkpoint-{fold}".encode())
        checkpoints.append(path)
    ranking, frozen = freeze_teacher_v2_selection(
        profile_prediction_paths={profile: [prediction_path]},
        profile_real_score_paths={profile: real_path},
        profile_checkpoint_paths={profile: checkpoints},
    )
    assert ranking.loc[0, "profile"] == profile
    assert frozen["architecture_frozen"]
    assert not frozen["s56_holdout_opened"]
    assert not frozen["s57_opened"]


def test_teacher_v2_recovery_uses_only_truth_matched_candidate_scores() -> None:
    manifest = pd.DataFrame({"injection_id": ["a", "b", "c"]})
    candidates = pd.DataFrame(
        {
            "injection_id": ["a", "a", "b", "c"],
            "rep_peak_rank": [1, 2, 3, 5],
            "is_injected_signal_peak": [True, False, True, False],
            "p_compact_transit": [0.8, 0.99, 0.4, 0.99],
        }
    )
    outcomes, summary = aggregate_compact_recovery(
        manifest, candidates, threshold=0.5
    )
    by_id = outcomes.set_index("injection_id")
    assert bool(by_id.loc["a", "teacher_v2_compact_recovered"])
    assert not bool(by_id.loc["b", "teacher_v2_compact_recovered"])
    assert not bool(by_id.loc["c", "teacher_v2_compact_recovered"])
    assert summary["n_bls_top5_recovered"] == 2
    assert summary["n_compact_recovered"] == 1
    topk = bls_topk_recovery_table(manifest, candidates)
    assert topk.set_index("top_k").loc[1, "n_recovered"] == 1
    assert topk.set_index("top_k").loc[3, "n_recovered"] == 2


def test_locked_human_holdout_is_same_row_and_enforces_recall() -> None:
    labels = ["planet_like", "eclipse_contact", "smooth_variable", "other"]
    human = pd.DataFrame(
        {
            "review_id": [f"row-{index}" for index in range(4)],
            "fixed_split": ["test"] * 4,
            "teacher_v2_role": ["real_human_morphology"] * 4,
            "morphology_target_index": np.arange(4),
            "human_label": [
                "planet_like",
                "eclipsing_binary_or_pceb",
                "stellar_variability",
                "uncertain",
            ],
        }
    )

    def scores() -> pd.DataFrame:
        frame = pd.DataFrame({"review_id": human["review_id"]})
        for label_index, label in enumerate(labels):
            frame[f"p_{label}"] = [
                0.97 if row_index == label_index else 0.01
                for row_index in range(4)
            ]
        return frame

    joined, summary = evaluate_locked_human_holdout(human, scores(), scores())
    assert len(joined) == 4
    assert summary["acceptance"]["passed"]
    assert summary["required_teacher_v2_recall"] == {
        "eclipse_contact": 1.0,
        "smooth_variable": 1.0,
        "other": 1.0,
    }


def test_locked_human_holdout_rejects_missing_model_row() -> None:
    human = pd.DataFrame(
        {
            "review_id": ["row-1"],
            "fixed_split": ["test"],
            "teacher_v2_role": ["real_human_morphology"],
            "morphology_target_index": [3],
            "human_label": ["uncertain"],
        }
    )
    empty = pd.DataFrame(
        columns=["review_id", *[f"p_{label}" for label in ("planet_like", "eclipse_contact", "smooth_variable", "other")]]
    )
    with pytest.raises(ValueError, match="missing 1 locked holdout"):
        evaluate_locked_human_holdout(human, empty, empty)


def test_prior_human_splits_allow_unmatched_skip_provenance() -> None:
    human = pd.DataFrame(
        {
            "review_id": ["new-active", "new-skip"],
            "source_review_id": ["old-active", "old-skip"],
            "tic": [1, 2],
            "human_label": ["planet_like", "skip"],
        }
    )
    prior = pd.DataFrame(
        {
            "review_id": ["old-active"],
            "fixed_split": ["development"],
            "cv_fold": [3],
        }
    )

    attached = attach_prior_human_splits(human, prior)

    assert attached.loc[0, "fixed_split"] == "development"
    assert attached.loc[0, "cv_fold"] == 3
    assert pd.isna(attached.loc[1, "fixed_split"])
    assert pd.isna(attached.loc[1, "cv_fold"])


def test_prior_human_splits_reject_unmatched_trainable_row() -> None:
    human = pd.DataFrame(
        {
            "review_id": ["new-active"],
            "source_review_id": ["missing-old-active"],
            "tic": [1],
            "human_label": ["planet_like"],
        }
    )
    prior = pd.DataFrame(
        {
            "review_id": ["different-row"],
            "fixed_split": ["development"],
            "cv_fold": [0],
        }
    )

    with pytest.raises(RuntimeError, match="trainable human transfers"):
        attach_prior_human_splits(human, prior)
