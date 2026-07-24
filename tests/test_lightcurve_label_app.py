from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd
import pytest

from twirl.vetting.lightcurve_label_app import (
    CandidateStore,
    LightCurveVettingApp,
    PERIOD_FACTOR_OPTIONS,
    _index_html,
    _normalize_period_factor,
    find_leo_report,
    find_leo_report_for_row,
    find_hlsp_path,
    find_twirl_vet_sheet_for_row,
    leo_class_from_report,
    tic_shard_path,
)

ROOT = Path(__file__).resolve().parents[1]


def _load_app_runner():
    path = ROOT / "scripts/stage5_validation/run_lightcurve_vetting_app.py"
    spec = importlib.util.spec_from_file_location("run_lightcurve_vetting_app_test", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_tic_shard_path_and_find_hlsp_path(tmp_path) -> None:
    root = tmp_path / "hlsp"
    tic = 267574918
    shard = tic_shard_path(root, tic)
    shard.mkdir(parents=True)
    path = shard / "hlsp_twirlfs_tess_ffi_s0056-0000000267574918_tess_v01_llc.fits"
    path.write_text("placeholder")

    assert shard == root / "0000" / "0002" / "6757" / "4918"
    assert find_hlsp_path(root, tic, 56) == path


def test_candidate_store_saves_and_reloads_labels(tmp_path) -> None:
    candidates = tmp_path / "candidates.csv"
    labels = tmp_path / "labels.csv"
    pd.DataFrame(
        [
            {
                "tic": 267574918,
                "sector": 56,
                "source_bucket": "wd1856_benchmark",
                "period_d": 1.4079,
                "t0_bjd": 2459825.47,
            }
        ]
    ).to_csv(candidates, index=False)

    store = CandidateStore(candidates_path=candidates, labels_out=labels)
    record = store.save_label(
        row_id=0,
        label="planet_like",
        labeler="tester",
        notes="clear transit",
        period_factor="0.5",
    )

    assert record["label"] == "planet_like"
    assert record["period_factor"] == "0.5"
    assert record["period_status"] == "refolded"
    assert labels.exists()

    reloaded = CandidateStore(candidates_path=candidates, labels_out=labels)
    row = reloaded.row(0)
    assert row["label"] == "planet_like"
    assert row["labeler"] == "tester"
    assert row["notes"] == "clear transit"
    assert row["period_factor"] == "0.5"
    assert row["period_status"] == "refolded"
    assert bool(row["reviewed"])


def test_candidate_store_prefills_without_marking_reviewed(tmp_path) -> None:
    candidates = tmp_path / "candidates.csv"
    labels = tmp_path / "labels.csv"
    pd.DataFrame(
        [
            {
                "tic": 267574918,
                "sector": 56,
                "period_d": 1.4079,
                "t0_bjd": 2459825.47,
                "initial_label": "planet_like",
                "initial_notes": "prior morphology",
                "initial_period_factor": "3",
                "initial_period_status": "refolded",
            }
        ]
    ).to_csv(candidates, index=False)

    store = CandidateStore(candidates_path=candidates, labels_out=labels)
    row = store.row(0)

    assert row["label"] == "planet_like"
    assert row["notes"] == "prior morphology"
    assert row["period_factor"] == "3"
    assert row["period_status"] == "refolded"
    assert not bool(row["reviewed"])
    assert store.summary()["reviewed"] == 0
    assert store.summary()["pending_review"] == 1
    assert store.summary()["labeled"] == 0
    assert store.summary()["unlabeled"] == 1
    assert store.summary()["label_counts"] == {}

    store.save_label(
        row_id=0,
        label="planet_like",
        labeler="tehan",
        notes="confirmed",
        period_factor="3",
    )

    assert store.summary()["reviewed"] == 1
    assert store.summary()["pending_review"] == 0
    assert store.summary()["labeled"] == 1
    assert store.summary()["unlabeled"] == 0
    assert store.summary()["label_counts"] == {"planet_like": 1}
    reloaded = CandidateStore(candidates_path=candidates, labels_out=labels)
    assert bool(reloaded.row(0)["reviewed"])
    assert reloaded.row(0)["notes"] == "confirmed"


def test_blank_note_save_remains_pending_review(tmp_path) -> None:
    candidates = tmp_path / "candidates.csv"
    labels = tmp_path / "labels.csv"
    pd.DataFrame(
        [
            {
                "tic": 267574918,
                "sector": 56,
                "period_d": 1.4079,
                "t0_bjd": 2459825.47,
                "initial_label": "planet_like",
            }
        ]
    ).to_csv(candidates, index=False)
    store = CandidateStore(candidates_path=candidates, labels_out=labels)

    store.save_label(
        row_id=0,
        label="",
        labeler="tehan",
        notes="note only",
    )

    assert not bool(store.row(0)["reviewed"])
    assert store.summary()["reviewed"] == 0
    assert store.summary()["pending_review"] == 1
    reloaded = CandidateStore(candidates_path=candidates, labels_out=labels)
    assert not bool(reloaded.row(0)["reviewed"])
    assert reloaded.row(0)["notes"] == "note only"


def test_candidate_store_hard_fails_candidate_key_mismatch(tmp_path) -> None:
    candidates = tmp_path / "candidates.csv"
    labels = tmp_path / "labels.csv"
    pd.DataFrame(
        [{"tic": 10, "sector": 56, "period_d": 1.0, "t0_bjd": 2.0, "source_bucket": "real"}]
    ).to_csv(candidates, index=False)
    pd.DataFrame(
        [
            {
                "row_id": 0,
                "candidate_key": "wrong",
                "label": "planet_like",
                "updated_utc": "2026-07-10T00:00:00Z",
            }
        ]
    ).to_csv(labels, index=False)

    with pytest.raises(ValueError, match="candidate_key mismatch"):
        CandidateStore(candidates_path=candidates, labels_out=labels)


def test_candidate_store_replaces_upstream_row_id(tmp_path) -> None:
    candidates = tmp_path / "scored_candidates.csv"
    labels = tmp_path / "labels.csv"
    pd.DataFrame(
        [
            {"row_id": 374, "tic": 10, "period_d": 1.0, "t0_bjd": 2.0},
            {"row_id": 812, "tic": 20, "period_d": 2.0, "t0_bjd": 3.0},
        ]
    ).to_csv(candidates, index=False)

    store = CandidateStore(candidates_path=candidates, labels_out=labels)

    assert store.frame["row_id"].tolist() == [0, 1]
    assert store.frame["tic"].tolist() == [10, 20]


def test_candidate_store_shuffles_review_order_without_changing_row_ids(tmp_path) -> None:
    candidates = tmp_path / "candidates.csv"
    labels = tmp_path / "labels.csv"
    pd.DataFrame(
        [
            {"tic": 1, "sector": 56, "source_bucket": "benchmark", "period_d": 1.0, "t0_bjd": 1.0},
            {"tic": 2, "sector": 56, "source_bucket": "candidate", "period_d": 2.0, "t0_bjd": 2.0},
            {"tic": 3, "sector": 56, "source_bucket": "injection", "period_d": 3.0, "t0_bjd": 3.0},
        ]
    ).to_csv(candidates, index=False)
    pd.DataFrame(
        [
            {
                "row_id": 0,
                "candidate_key": "1|56|1.0|1.0|benchmark",
                "label": "planet_like",
                "label_source": "human",
                "labeler": "tester",
                "notes": "",
                "updated_utc": "2026-06-18T00:00:00+00:00",
            }
        ]
    ).to_csv(labels, index=False)

    store = CandidateStore(
        candidates_path=candidates,
        labels_out=labels,
        shuffle_order=True,
        random_seed=4,
        unlabeled_first=True,
    )
    first = store.row(0)

    assert first["row_id"] != 0
    assert first["label"] == ""
    record = store.save_label(
        row_id=int(first["row_id"]),
        label="uncertain",
        labeler="tester",
        notes="shuffled triage",
    )
    assert record["tic"] == int(first["tic"])

    reloaded = CandidateStore(candidates_path=candidates, labels_out=labels)
    labeled = reloaded.frame.set_index("row_id").loc[int(first["row_id"])]
    assert labeled["label"] == "uncertain"


def test_candidate_store_prioritizes_unreviewed_prefilled_rows(tmp_path) -> None:
    candidates = tmp_path / "candidates.csv"
    labels = tmp_path / "labels.csv"
    rows = [
        {
            "tic": tic,
            "sector": 56,
            "period_d": float(tic),
            "t0_bjd": float(tic),
            "initial_label": "planet_like",
        }
        for tic in (1, 2, 3)
    ]
    pd.DataFrame(rows).to_csv(candidates, index=False)
    pd.DataFrame(
        [
            {
                "row_id": 0,
                "candidate_key": "1|56|1|1|",
                "label": "planet_like",
                "label_source": "human",
                "labeler": "tester",
                "notes": "",
                "updated_utc": "2026-07-24T00:00:00+00:00",
                "period_factor": "1",
                "period_status": "review_period",
            }
        ]
    ).to_csv(labels, index=False)

    store = CandidateStore(
        candidates_path=candidates,
        labels_out=labels,
        shuffle_order=True,
        random_seed=4,
        unlabeled_first=True,
    )

    assert not bool(store.row(0)["reviewed"])
    assert int(store.row(0)["row_id"]) != 0
    assert bool(store.row(store.count - 1)["reviewed"])


def test_candidate_payload_hides_source_until_labeled(tmp_path) -> None:
    candidates = tmp_path / "candidates.csv"
    labels = tmp_path / "labels.csv"
    hlsp_root = tmp_path / "hlsp"
    hlsp_root.mkdir()
    pd.DataFrame(
        [
            {
                "tic": 267574918,
                "sector": 56,
                "source_bucket": "wd1856_benchmark",
                "vet_class": "injected_wd1856_like",
                "period_d": 1.4079,
                "t0_bjd": 2459825.47,
                "label": "",
                "label_source": "human",
            }
        ]
    ).to_csv(candidates, index=False)

    app = LightCurveVettingApp(
        candidates_path=candidates,
        labels_out=labels,
        hlsp_root=hlsp_root,
        labeler="tester",
    )
    payload = app.candidate_payload(0)

    assert payload["label"] == ""
    assert payload["label_source"] == ""
    assert "source_bucket" not in payload["display"]
    assert payload["display"]["vet_class"] == "review_candidate"


def test_leo_report_lookup_and_payload(tmp_path) -> None:
    candidates = tmp_path / "candidates.csv"
    labels = tmp_path / "labels.csv"
    hlsp_root = tmp_path / "hlsp"
    leo_root = tmp_path / "leo" / "vet_reports"
    hlsp_root.mkdir()
    leo_root.mkdir(parents=True)
    report = leo_root / "PC_rank11_tic0267574918_T16.34_P001.4080d.pdf"
    report.write_bytes(b"%PDF-1.4\n")
    pd.DataFrame(
        [
            {
                "tic": 267574918,
                "sector": 56,
                "source_bucket": "wd1856_benchmark",
                "period_d": 1.4079,
                "t0_bjd": 2459825.47,
            }
        ]
    ).to_csv(candidates, index=False)

    assert find_leo_report((leo_root,), 267574918) == report
    assert leo_class_from_report(report) == "PC"

    app = LightCurveVettingApp(
        candidates_path=candidates,
        labels_out=labels,
        hlsp_root=hlsp_root,
        leo_report_roots=(leo_root,),
        labeler="tester",
    )
    payload = app.candidate_payload(0)

    assert payload["leo_report_path"] == str(report)
    assert payload["leo_report_name"] == report.name
    assert payload["leo_class"] == "PC"


def test_candidate_payload_prefers_row_specific_leo_report(tmp_path) -> None:
    candidates = tmp_path / "candidates.csv"
    labels = tmp_path / "labels.csv"
    hlsp_root = tmp_path / "hlsp"
    leo_root = tmp_path / "leo" / "vet_reports"
    hlsp_root.mkdir()
    leo_root.mkdir(parents=True)
    stale = leo_root / "FA_row0001_tic0267574918_DET_FLUX_P001.4079d.pdf"
    wanted = leo_root / "PC_row0099_tic0267574918_DET_FLUX_SML_P001.4079d.pdf"
    stale.write_bytes(b"%PDF-1.4\n")
    wanted.write_bytes(b"%PDF-1.4\n")
    pd.DataFrame(
        [
            {
                "tic": 267574918,
                "sector": 56,
                "period_d": 1.4079,
                "t0_bjd": 2459825.47,
                "leo_report_name": wanted.name,
            }
        ]
    ).to_csv(candidates, index=False)

    app = LightCurveVettingApp(
        candidates_path=candidates,
        labels_out=labels,
        hlsp_root=hlsp_root,
        leo_report_roots=(leo_root,),
        labeler="tester",
    )
    row = app.store.row(0)
    payload = app.candidate_payload(0)

    assert find_leo_report_for_row((leo_root,), row) == wanted
    assert payload["leo_report_path"] == str(wanted)
    assert payload["leo_class"] == "PC"


def test_candidate_payload_marks_fallback_leo_report(tmp_path) -> None:
    candidates = tmp_path / "candidates.csv"
    labels = tmp_path / "labels.csv"
    hlsp_root = tmp_path / "hlsp"
    leo_root = tmp_path / "leo" / "vet_reports"
    hlsp_root.mkdir()
    leo_root.mkdir(parents=True)
    report = leo_root / "FA_row0012_real_candidate_tic2041296453_DET_FLUX_LAG_P008.4710d.pdf"
    report.write_bytes(b"%PDF-1.4\n")
    pd.DataFrame(
        [
            {
                "review_id": "real:2041296453",
                "tic": 2041296453,
                "sector": 56,
                "source_bucket": "planet_centroid_pass",
                "period_d": 8.4709,
                "t0_bjd": 2459828.9449,
            }
        ]
    ).to_csv(candidates, index=False)
    pd.DataFrame(
        [
            {
                "review_id": "real:2041296453",
                "tic": 2041296453,
                "plot_error": "ValueError: cannot convert float NaN to integer",
                "error": "",
                "leo_class": "FA",
                "leo_report_name": report.name,
                "leo_report_path": str(report),
            }
        ]
    ).to_csv(tmp_path / "leo_metrics.csv", index=False)

    app = LightCurveVettingApp(
        candidates_path=candidates,
        labels_out=labels,
        hlsp_root=hlsp_root,
        leo_report_roots=(leo_root,),
        labeler="tester",
    )
    payload = app.candidate_payload(0)

    assert payload["leo_report_kind"] == "fallback_plot"
    assert payload["leo_plot_error"].startswith("ValueError:")


def test_period_factor_controls_include_thirds() -> None:
    assert ("0.3333333333333333", "P/3") in PERIOD_FACTOR_OPTIONS
    assert ("3", "3P") in PERIOD_FACTOR_OPTIONS
    assert _normalize_period_factor("0.333333")[0] == "0.3333333333333333"
    assert _normalize_period_factor("3.0") == ("3", "refolded")


def test_exact_twirl_vet_sheet_mode_disables_all_fallbacks(tmp_path) -> None:
    candidates = tmp_path / "candidates.csv"
    labels = tmp_path / "labels.csv"
    hlsp_root = tmp_path / "hlsp"
    vet_root = tmp_path / "vet_sheets"
    hlsp_root.mkdir()
    vet_root.mkdir()
    fallback = vet_root / "real_42_twirl_twoap_twirl_fs_v2_adp015q.png"
    fallback.write_bytes(b"\x89PNG\r\n\x1a\nfallback")
    pd.DataFrame(
        [
            {
                "review_id": "real:42",
                "tic": 42,
                "sector": 56,
                "period_d": 1.0,
                "t0_bjd": 2.0,
                "twirl_vet_sheet_name": "exact_current.png",
            }
        ]
    ).to_csv(candidates, index=False)
    row = pd.read_csv(candidates).iloc[0]

    assert find_twirl_vet_sheet_for_row((vet_root,), row) == fallback
    assert (
        find_twirl_vet_sheet_for_row(
            (vet_root,),
            row,
            exact_name_only=True,
        )
        is None
    )

    app = LightCurveVettingApp(
        candidates_path=candidates,
        labels_out=labels,
        hlsp_root=hlsp_root,
        twirl_vet_roots=(vet_root,),
        labeler="tehan",
        exact_twirl_vet_sheets=True,
    )
    assert app.candidate_payload(0)["twirl_vet_sheet_path"] is None

    exact = vet_root / "exact_current.png"
    exact.write_bytes(b"\x89PNG\r\n\x1a\nexact")
    payload = app.candidate_payload(0)
    assert payload["twirl_vet_sheet_path"] == str(exact)
    assert payload["twirl_vet_sheet_name"] == exact.name

    runner = _load_app_runner()
    args = runner._build_arg_parser().parse_args(["--exact-twirl-vet-sheets"])
    assert args.exact_twirl_vet_sheets


def test_vetting_app_has_quick_label_autosave_controls() -> None:
    html = _index_html()

    assert 'data-shortcut="1"' in html
    assert 'data-label="planet_like"' in html
    assert 'data-label="centroid_contaminant"' not in html
    assert "async function labelAndNext" in html
    assert 'document.addEventListener("keydown"' in html
    assert "function preloadNextCandidates" in html
    assert "function preloadPdf" in html
    assert "`<div>leo_report</div>" not in html
    assert "Save notes" in html
