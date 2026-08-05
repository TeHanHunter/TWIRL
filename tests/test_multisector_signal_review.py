from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd
import pytest

from twirl.vetting.label_io import candidate_key
from twirl.vetting.multisector_signal_review import (
    ACCEPTED_UNCHANGED_REVIEW_STATUS,
    EXPLICIT_OVERRIDE_REVIEW_STATUS,
    build_signal_rereview_queue,
    finalize_signal_rereview,
    materialize_signal_rereview_acceptance,
    normalize_accepted_franklin_signals,
    normalize_browser_signal_rows,
    normalize_s56_adjudicated_signals,
    tehan_app_candidate_key,
)


ROOT = Path(__file__).resolve().parents[1]


def _load_app():
    path = ROOT / "scripts/stage5_validation/franklin_vetting_app.py"
    spec = importlib.util.spec_from_file_location("franklin_vetting_app_signal_test", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _load_finalizer():
    path = (
        ROOT
        / "scripts/stage5_validation/finalize_s56_s62_signal_rereview.py"
    )
    spec = importlib.util.spec_from_file_location(
        "finalize_s56_s62_signal_rereview_test",
        path,
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _s56() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "review_id": ["real:10", "real:20"],
            "tic": [10, 20],
            "sector": [56, 56],
            "period_d": [1.5, 2.5],
            "t0_bjd": [2_459_000.1, 2_459_000.2],
            "duration_min": [10.0, 12.0],
            "human_label": ["planet_like", "eclipsing_binary_or_pceb"],
            "human_labeler": ["tehan", "tehan"],
            "human_notes": ["", ""],
            "human_updated_utc": ["2026-07-01", "2026-07-01"],
            "human_label_adjudicated": [
                "planet_like",
                "eclipsing_binary_or_pceb",
            ],
            "adjudication_final": [True, True],
            "adjudicated_labeler": ["tehan", "tehan"],
            "adjudicated_notes": ["best at 3P", ""],
            "adjudicated_updated_utc": ["2026-07-02", "2026-07-02"],
            "adjudicated_period_status": ["resolved", "resolved"],
            "effective_period_factor": [3.0, 1.0],
            "harmonic_target_v1": ["3p", "p"],
            "harmonic_include_v1": [True, True],
            "is_injected_row": [False, False],
            "source_uid": ["s56:10", "s56:20"],
            "twirl_vet_sheet_name": ["real_10.png", "real_20.png"],
        }
    )


def _revisit() -> tuple[pd.DataFrame, pd.DataFrame]:
    queue = pd.DataFrame(
        {
            "row_id": [0],
            "review_id": ["s0056-revisit-10"],
            "tic": [10],
            "sector": [56],
            "period_d": [1.5],
            "t0_bjd": [2_459_000.1],
            "duration_min": [10.0],
            "twirl_vet_sheet_name": ["s0056-revisit-10.png"],
            "source_bucket": [""],
        }
    )
    queue["candidate_key"] = queue.apply(candidate_key, axis=1)
    labels = pd.DataFrame(
        {
            "row_id": [0],
            "candidate_key": queue["candidate_key"],
            "tic": [10],
            "sector": [56],
            "label": ["planet_like"],
            "labeler": ["tehan"],
            "notes": ["later morphology pass"],
            "period_factor": ["1"],
            "period_status": ["review_period"],
            "updated_utc": ["2026-07-03"],
        }
    )
    return queue, labels


def _franklin() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "row_id": [0],
            "review_id": ["s0057-franklin-10"],
            "tic": [10],
            "sector": [57],
            "period_d": [1.6],
            "t0_bjd": [2_459_030.1],
            "duration_min": [11.0],
            "human_label": ["eclipsing_binary_or_pceb"],
            "human_labeler": ["franklin"],
            "human_notes": ["possible 2P"],
            "human_updated_utc": ["2026-07-21"],
            "reported_period_factor": ["2"],
            "reported_period_status": ["resolved"],
            "source_batch_id": ["s57_s59"],
            "source_uid": ["s57_s59:10"],
            "pipeline_candidate_key": ["10|57|1.6|2459030.1|"],
            "factor_review_status": ["not_explicitly_reviewed"],
            "harmonic_supervision_verified": [False],
            "harmonic_include_v1": [False],
            "twirl_vet_sheet_name": ["s0057-franklin-10.png"],
        }
    )


def _current_queue() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    revisit_queue, revisit_labels = _revisit()
    sources = [
        normalize_s56_adjudicated_signals(_s56()),
        normalize_browser_signal_rows(
            revisit_queue,
            revisit_labels,
            source_batch_id="s56_revisit",
        ),
        normalize_accepted_franklin_signals(_franklin()),
    ]
    queue, provenance, assets, _ = build_signal_rereview_queue(sources)
    return queue, provenance, assets


def test_queue_collapses_exact_duplicate_but_keeps_sector_observations() -> None:
    queue, provenance, assets = _current_queue()
    assert len(provenance) == 4
    assert len(queue) == 3
    assert len(assets) == 3
    assert queue["candidate_key"].ne("").all()
    assert queue["candidate_key"].equals(queue["observation_candidate_key"])
    assert queue.groupby(["sector", "tic"]).size().to_dict() == {
        (56, 10): 1,
        (56, 20): 1,
        (57, 10): 1,
    }
    s56_tic10 = queue.loc[(queue["sector"] == 56) & (queue["tic"] == 10)].iloc[0]
    assert s56_tic10["prior_labeler"] == "tehan"
    assert s56_tic10["initial_period_factor"] == "3"
    assert bool(s56_tic10["original_harmonic_supervision_verified"])
    assert s56_tic10["original_harmonic_target"] == "3p"
    s57 = queue.loc[queue["sector"] == 57].iloc[0]
    assert s57["initial_period_factor"] == "2"
    assert not bool(s57["original_harmonic_supervision_verified"])
    label_order = queue["prior_label"].tolist()
    assert label_order == [
        "planet_like",
        "eclipsing_binary_or_pceb",
        "eclipsing_binary_or_pceb",
    ]


def test_app_prefills_without_marking_reviewed_and_searches_multiple_roots(
    tmp_path: Path,
) -> None:
    queue, _, _ = _current_queue()
    queue_path = tmp_path / "queue.csv"
    queue.to_csv(queue_path, index=False, float_format="%.15g")
    labels_path = tmp_path / "labels.csv"
    pd.DataFrame(
        columns=[
            "row_id",
            "candidate_key",
            "tic",
            "sector",
            "label",
            "label_source",
            "labeler",
            "notes",
            "period_factor",
            "period_status",
            "updated_utc",
        ]
    ).to_csv(labels_path, index=False)
    root_a = tmp_path / "a"
    root_b = tmp_path / "b"
    root_a.mkdir()
    root_b.mkdir()
    names = queue["twirl_vet_sheet_name"].tolist()
    (root_a / names[0]).write_bytes(b"\x89PNG\r\n\x1a\n")
    (root_b / names[1]).write_bytes(b"\x89PNG\r\n\x1a\n")
    (root_b / names[2]).write_bytes(b"\x89PNG\r\n\x1a\n")

    app = _load_app()
    store = app.Store(queue_path, labels_path, [root_a, root_b], "tehan")
    assert store.summary()["reviewed"] == 0
    assert store.summary()["pending_review"] == 3
    assert store.payload(0)["label"] in {
        "planet_like",
        "eclipsing_binary_or_pceb",
    }
    assert not store.payload(0)["reviewed"]
    assert all(store.sheet_path(index) is not None for index in range(3))
    row = store.row(0)
    record = store.save_label(
        0,
        row["label"],
        "tehan",
        row["notes"],
        period_factor=row["period_factor"],
    )
    assert record["period_factor"] in {
        "0.25",
        "0.3333333333333333",
        "0.5",
        "1",
        "2",
        "3",
        "4",
        "unresolved",
    }
    assert store.summary()["reviewed"] == 1

    exact_store = app.Store(
        queue_path,
        labels_path,
        [root_a, root_b],
        "tehan",
        allow_sheet_fallback=False,
    )
    assert exact_store.summary()["allow_sheet_fallback"] is False
    assert all(exact_store.sheet_path(index) is not None for index in range(3))

    exact_name = names[0]
    (root_a / exact_name).unlink()
    wildcard_name = (
        f"fallback_{queue.iloc[0]['tic']}_twirl_twoap_unrelated.png"
    )
    (root_a / wildcard_name).write_bytes(b"\x89PNG\r\n\x1a\n")
    assert store.sheet_path(0) is not None
    assert exact_store.sheet_path(0) is None


def test_final_review_preserves_only_previously_verified_harmonics() -> None:
    queue, _, _ = _current_queue()
    string_queue = queue.astype(object).where(pd.notna(queue), "")
    labels = []
    for _, row in string_queue.iterrows():
        labels.append(
            {
                "row_id": row["row_id"],
                "candidate_key": tehan_app_candidate_key(row),
                "tic": row["tic"],
                "sector": row["sector"],
                "label": row["initial_label"],
                "label_source": "human",
                "labeler": "tehan",
                "notes": "",
                "period_factor": row["initial_period_factor"],
                "period_status": "resolved",
                "updated_utc": "2026-07-24T00:00:00+00:00",
            }
        )
    final = finalize_signal_rereview(
        string_queue,
        pd.DataFrame(labels),
        adjudicator="tehan",
        accepted_utc="2026-07-24T00:01:00+00:00",
    )
    s56 = final.loc[final["sector"].astype(int).eq(56)]
    s57 = final.loc[final["sector"].astype(int).eq(57)]
    assert s56["harmonic_include_v1"].astype(bool).all()
    assert set(s56["harmonic_target_v1"]) == {"3p", "p"}
    assert not s57["harmonic_include_v1"].astype(bool).any()
    assert s57["harmonic_target_v1"].eq("").all()

    incomplete = pd.DataFrame(labels).iloc[:-1]
    with pytest.raises(ValueError, match="incomplete"):
        finalize_signal_rereview(
            string_queue,
            incomplete,
            adjudicator="tehan",
            accepted_utc="2026-07-24T00:01:00+00:00",
        )


def test_declared_full_pass_materializes_only_unchanged_rows() -> None:
    queue, _, _ = _current_queue()
    string_queue = queue.astype(object).where(pd.notna(queue), "")
    changed_row = string_queue.iloc[0]
    saved = pd.DataFrame(
        [
            {
                "row_id": changed_row["row_id"],
                "candidate_key": tehan_app_candidate_key(changed_row),
                "tic": changed_row["tic"],
                "sector": changed_row["sector"],
                "label": "wide_transit_like",
                "label_source": "human",
                "labeler": "tehan",
                "notes": "explicit change",
                "period_factor": "2",
                "period_status": "refolded",
                "updated_utc": "2026-07-24T19:00:00+00:00",
            }
        ]
    )
    accepted_utc = "2026-07-24T20:10:00+00:00"
    materialized = materialize_signal_rereview_acceptance(
        string_queue,
        saved,
        adjudicator="tehan",
        accepted_utc=accepted_utc,
    )

    assert len(materialized) == len(string_queue)
    assert materialized["row_id"].tolist() == string_queue["row_id"].map(str).tolist()
    assert materialized["review_acceptance_mode"].value_counts().to_dict() == {
        ACCEPTED_UNCHANGED_REVIEW_STATUS: len(string_queue) - 1,
        EXPLICIT_OVERRIDE_REVIEW_STATUS: 1,
    }
    explicit = materialized.loc[
        materialized["review_acceptance_mode"].eq(
            EXPLICIT_OVERRIDE_REVIEW_STATUS
        )
    ].iloc[0]
    assert explicit["label"] == "wide_transit_like"
    assert explicit["notes"] == "explicit change"
    assert explicit["updated_utc"] == "2026-07-24T19:00:00+00:00"
    unchanged = materialized.loc[
        materialized["review_acceptance_mode"].eq(
            ACCEPTED_UNCHANGED_REVIEW_STATUS
        )
    ]
    expected = string_queue.set_index(string_queue["row_id"].map(str))
    assert unchanged["label"].tolist() == [
        expected.loc[row_id, "initial_label"]
        for row_id in unchanged["row_id"]
    ]
    assert unchanged["updated_utc"].eq(accepted_utc).all()

    final = finalize_signal_rereview(
        string_queue,
        materialized,
        adjudicator="tehan",
        accepted_utc=accepted_utc,
    )
    assert final["morphology_review_status"].value_counts().to_dict() == {
        ACCEPTED_UNCHANGED_REVIEW_STATUS: len(string_queue) - 1,
        EXPLICIT_OVERRIDE_REVIEW_STATUS: 1,
    }
    assert set(
        final.loc[
            final["morphology_review_status"].eq(
                ACCEPTED_UNCHANGED_REVIEW_STATUS
            ),
            "human_label_source",
        ]
    ) == {"human_acceptance_of_prior_review"}

    tampered = saved.copy()
    tampered.loc[0, "candidate_key"] = "wrong"
    with pytest.raises(ValueError, match="exact frozen queue"):
        materialize_signal_rereview_acceptance(
            string_queue,
            tampered,
            adjudicator="tehan",
            accepted_utc=accepted_utc,
        )


def test_disk_finalizer_binds_exact_selected_source_provenance(
    tmp_path: Path,
) -> None:
    queue, provenance, _ = _current_queue()
    string_queue = queue.astype(object).where(pd.notna(queue), "")
    labels = pd.DataFrame(
        [
            {
                "row_id": row["row_id"],
                "candidate_key": tehan_app_candidate_key(row),
                "tic": row["tic"],
                "sector": row["sector"],
                "label": row["initial_label"],
                "label_source": "human",
                "labeler": "tehan",
                "notes": "",
                "period_factor": row["initial_period_factor"],
                "period_status": "resolved",
                "updated_utc": "2026-07-24T00:00:00+00:00",
            }
            for _, row in string_queue.iterrows()
        ]
    )
    queue_path = tmp_path / "queue.csv"
    labels_path = tmp_path / "labels.csv"
    provenance_path = tmp_path / "provenance.csv"
    string_queue.to_csv(queue_path, index=False)
    labels.to_csv(labels_path, index=False)
    provenance.to_csv(provenance_path, index=False)

    finalizer = _load_finalizer()
    summary = finalizer.finalize(
        queue_path=queue_path,
        labels_path=labels_path,
        provenance_path=provenance_path,
        out_dir=tmp_path / "accepted",
        adjudicator="tehan",
        accepted_utc="2026-07-24T00:01:00+00:00",
    )
    assert summary["n_rows"] == len(queue)
    assert summary["acceptance"]["review_status_counts"] == {
        EXPLICIT_OVERRIDE_REVIEW_STATUS: len(queue)
    }
    assert (
        tmp_path / "accepted" / "accepted_label_decisions.csv"
    ).is_file()

    tampered = provenance.copy()
    selected = tampered["selected_for_queue"].astype(bool)
    selected_indices = tampered.index[selected].tolist()
    assert len(selected_indices) >= 2
    first, second = selected_indices[:2]
    tampered.loc[
        [first, second],
        "observation_candidate_key",
    ] = tampered.loc[
        [second, first],
        "observation_candidate_key",
    ].to_numpy()
    tampered_path = tmp_path / "tampered_provenance.csv"
    tampered.to_csv(tampered_path, index=False)
    with pytest.raises(ValueError, match="does not match the exact queue"):
        finalizer.finalize(
            queue_path=queue_path,
            labels_path=labels_path,
            provenance_path=tampered_path,
            out_dir=tmp_path / "tampered",
            adjudicator="tehan",
            accepted_utc="2026-07-24T00:01:00+00:00",
        )


def test_disk_finalizer_requires_declaration_for_unchanged_priors(
    tmp_path: Path,
) -> None:
    queue, provenance, _ = _current_queue()
    string_queue = queue.astype(object).where(pd.notna(queue), "")
    row = string_queue.iloc[0]
    saved = pd.DataFrame(
        [
            {
                "row_id": row["row_id"],
                "candidate_key": tehan_app_candidate_key(row),
                "tic": row["tic"],
                "sector": row["sector"],
                "label": "wide_transit_like",
                "label_source": "human",
                "labeler": "tehan",
                "notes": "",
                "period_factor": "1",
                "period_status": "review_period",
                "updated_utc": "2026-07-24T19:00:00+00:00",
            }
        ]
    )
    queue_path = tmp_path / "queue.csv"
    labels_path = tmp_path / "labels.csv"
    provenance_path = tmp_path / "provenance.csv"
    string_queue.to_csv(queue_path, index=False)
    saved.to_csv(labels_path, index=False)
    provenance.to_csv(provenance_path, index=False)
    finalizer = _load_finalizer()

    with pytest.raises(ValueError, match="acceptance declaration"):
        finalizer.finalize(
            queue_path=queue_path,
            labels_path=labels_path,
            provenance_path=provenance_path,
            out_dir=tmp_path / "missing_declaration",
            adjudicator="tehan",
            accepted_utc="2026-07-24T20:10:00+00:00",
            accept_unchanged_priors=True,
        )

    summary = finalizer.finalize(
        queue_path=queue_path,
        labels_path=labels_path,
        provenance_path=provenance_path,
        out_dir=tmp_path / "accepted_declared_pass",
        adjudicator="tehan",
        accepted_utc="2026-07-24T20:10:00+00:00",
        accept_unchanged_priors=True,
        acceptance_declaration="Reviewed all rows; unchanged priors accepted.",
    )
    assert summary["acceptance"]["review_status_counts"] == {
        ACCEPTED_UNCHANGED_REVIEW_STATUS: len(queue) - 1,
        EXPLICIT_OVERRIDE_REVIEW_STATUS: 1,
    }
