from __future__ import annotations

import pandas as pd

from twirl.vetting.lightcurve_label_app import (
    CandidateStore,
    LightCurveVettingApp,
    _index_html,
    find_leo_report,
    find_hlsp_path,
    leo_class_from_report,
    tic_shard_path,
)


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
    )

    assert record["label"] == "planet_like"
    assert labels.exists()

    reloaded = CandidateStore(candidates_path=candidates, labels_out=labels)
    row = reloaded.row(0)
    assert row["label"] == "planet_like"
    assert row["labeler"] == "tester"
    assert row["notes"] == "clear transit"


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


def test_vetting_app_has_quick_label_autosave_controls() -> None:
    html = _index_html()

    assert 'data-shortcut="1"' in html
    assert 'data-label="planet_like"' in html
    assert 'data-label="centroid_contaminant"' not in html
    assert "async function labelAndNext" in html
    assert 'document.addEventListener("keydown"' in html
    assert "function preloadNextCandidates" in html
    assert "function preloadPdf" in html
    assert "Save notes" in html
