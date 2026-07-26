from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from twirl.vetting.adjudication_audit import HARMONIC_CNN_TARGET_POLICY
from twirl.vetting.teacher_v3 import (
    TEACHER_V3_CORPUS_VERSION,
    TEACHER_V3_RUN_NAME,
    build_teacher_v3_corpus,
)
from twirl.vetting.teacher_v3_preparation import (
    prepare_teacher_v3_native_tables,
)
from twirl.vetting.teacher_split_registry import build_tic_split_registry


def _source_row(
    *,
    sector: int | str,
    tic: int,
    source_uid: str,
    label: str = "uncertain",
) -> dict[str, object]:
    morphology = {
        "planet_like": "planet_like",
        "eclipsing_binary_or_pceb": "eclipse_contact",
        "stellar_variability": "smooth_variable",
        "instrumental_or_systematic": "other",
        "uncertain": "other",
    }.get(label, "")
    preserve = (
        "preserve"
        if label
        in {
            "planet_like",
            "eclipsing_binary_or_pceb",
            "stellar_variability",
            "wide_transit_like",
        }
        else ("reject" if label in {"uncertain", "instrumental_or_systematic"} else "")
    )
    return {
        "review_id": f"review-{source_uid}",
        "source_uid": source_uid,
        "source_kind": "real_candidate",
        "source_bucket": "review_candidate",
        "tic": tic,
        "sector": sector,
        "period_d": 1.0 + 0.01 * tic,
        "t0_bjd": 2_459_000.0 + tic,
        "duration_min": 10.0,
        "human_label": label,
        "human_label_source": "human",
        "human_labeler": "tehan",
        "human_notes": "",
        "human_updated_utc": "2026-07-24T00:00:00Z",
        "morphology_target_v1": morphology,
        "morphology_include_v1": bool(morphology),
        "preserve_target_v1": preserve,
        "preserve_include_v1": bool(preserve),
        "harmonic_target_v1": "",
        "harmonic_include_v1": False,
        "broad_preserve_only": label == "wide_transit_like",
        "model_target_policy_version": HARMONIC_CNN_TARGET_POLICY,
        "is_injected_row": False,
        "injection_id": "",
    }


def _decision(
    source: dict[str, object],
    *,
    row_id: int,
    label: str,
) -> dict[str, object]:
    targets = _source_row(
        sector=source["sector"],
        tic=int(source["tic"]),
        source_uid=str(source["source_uid"]),
        label=label,
    )
    return {
        **source,
        "row_id": row_id,
        "selected_source_uid": source["source_uid"],
        "observation_candidate_key": (
            f"{source['tic']}|{source['sector']}|{source['period_d']}|"
            f"{source['t0_bjd']}|review_candidate"
        ),
        "morphology_review_status": "explicit_override",
        "rereview_period_factor": "1",
        "rereview_period_status": "review_period",
        "rereview_period_is_audit_only": True,
        "final_human_label": label,
        "final_labeler": "tehan",
        "final_notes": "",
        "final_updated_utc": "2026-07-24T20:11:02Z",
        "human_label_source": "human_adjudication",
        "morphology_target_v1": targets["morphology_target_v1"],
        "morphology_include_v1": targets["morphology_include_v1"],
        "preserve_target_v1": targets["preserve_target_v1"],
        "preserve_include_v1": targets["preserve_include_v1"],
        "harmonic_target_v1": "",
        "harmonic_include_v1": False,
        "harmonic_supervision_verified": False,
        "broad_preserve_only": targets["broad_preserve_only"],
        "model_target_policy_version": HARMONIC_CNN_TARGET_POLICY,
    }


def test_teacher_v3_corpus_applies_matches_and_adds_only_s56() -> None:
    s56 = _source_row(
        sector="",
        tic=5601,
        source_uid="s56:existing",
    )
    s57 = _source_row(
        sector=57,
        tic=5701,
        source_uid="s57:existing",
    )
    s58 = _source_row(
        sector=58,
        tic=5801,
        source_uid="s58:existing",
    )
    s59 = _source_row(
        sector=59,
        tic=5901,
        source_uid="s59:existing",
    )
    s60 = _source_row(
        sector=60,
        tic=6001,
        source_uid="s60:existing",
    )
    s61 = _source_row(
        sector=61,
        tic=6101,
        source_uid="s61:existing",
    )
    s62 = _source_row(
        sector=62,
        tic=6201,
        source_uid="s62:existing",
    )
    additional = _source_row(
        sector=56,
        tic=5699,
        source_uid="s56:additional",
    )
    decisions = pd.DataFrame(
        [
            _decision(s57, row_id=0, label="planet_like"),
            _decision(
                additional,
                row_id=1,
                label="eclipsing_binary_or_pceb",
            ),
        ]
    )

    corpus, summary = build_teacher_v3_corpus(
        s56_training=pd.DataFrame([s56]),
        s57_s59_labels=pd.DataFrame([s57, s58, s59]),
        s60_s62_labels=pd.DataFrame([s60, s61, s62]),
        signal_freeze=decisions,
    )

    assert len(corpus) == 8
    assert set(corpus["sector"]) == set(range(56, 63))
    assert corpus["teacher_v3_observation_id"].is_unique
    assert corpus["observation_candidate_key"].ne("").all()
    assert corpus["signal_rereview_applied"].astype(bool).sum() == 2
    assert (
        corpus.set_index("source_uid").loc["s57:existing", "human_label"]
        == "planet_like"
    )
    assert (
        corpus.set_index("source_uid").loc[
            "s56:additional", "human_label"
        ]
        == "eclipsing_binary_or_pceb"
    )
    assert summary["n_signal_review_matched_existing"] == 1
    assert summary["n_signal_review_additional_s56"] == 1
    assert summary["teacher_run_name"] == TEACHER_V3_RUN_NAME
    assert summary["corpus_contract_version"] == TEACHER_V3_CORPUS_VERSION

    registry, _ = build_tic_split_registry(
        corpus,
        seed=560062,
        stratum_column="split_stratum",
        identity_columns=(
            "sector",
            "tic",
            "teacher_v3_observation_id",
        ),
    )
    (
        prepared,
        per_sector,
        native_source,
        injection_mapping,
        preparation_summary,
    ) = prepare_teacher_v3_native_tables(
        corpus=corpus,
        split_registry=registry,
        native_root=Path("/orcd/native"),
    )
    assert len(prepared) == len(corpus)
    assert set(per_sector) == set(range(56, 63))
    assert len(native_source) == len(corpus)
    assert injection_mapping.empty
    assert prepared["native_group_path"].str.startswith("targets/").all()
    assert preparation_summary["group_leakage_count"] == 0


def test_teacher_v3_rejects_unmatched_non_s56_decision() -> None:
    sources = {
        sector: _source_row(
            sector=sector,
            tic=sector * 100,
            source_uid=f"s{sector}:existing",
        )
        for sector in range(56, 63)
    }
    unmatched = _source_row(
        sector=62,
        tic=6299,
        source_uid="s62:missing",
    )
    with pytest.raises(ValueError, match="only additional S56"):
        build_teacher_v3_corpus(
            s56_training=pd.DataFrame([sources[56]]),
            s57_s59_labels=pd.DataFrame(
                [sources[57], sources[58], sources[59]]
            ),
            s60_s62_labels=pd.DataFrame(
                [sources[60], sources[61], sources[62]]
            ),
            signal_freeze=pd.DataFrame(
                [_decision(unmatched, row_id=0, label="planet_like")]
            ),
        )
