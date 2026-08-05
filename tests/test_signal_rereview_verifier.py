from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest
from PIL import Image

from twirl.vetting.label_io import candidate_key


ROOT = Path(__file__).resolve().parents[1]


def _load_verifier():
    path = (
        ROOT
        / "scripts/stage5_validation/verify_s56_s62_signal_rereview.py"
    )
    spec = importlib.util.spec_from_file_location(
        "verify_s56_s62_signal_rereview_test",
        path,
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_bundle(bundle: Path) -> None:
    rows = [
        {
            "row_id": 0,
            "review_id": "planet-row",
            "selected_source_uid": "planet-source",
            "sector": 56,
            "tic": 101,
            "period_d": 0.2,
            "t0_bjd": 2459000.1,
            "duration_min": 15.0,
            "prior_label": "planet_like",
            "twirl_vet_sheet_name": "planet.png",
        },
        {
            "row_id": 1,
            "review_id": "eb-row",
            "selected_source_uid": "eb-source",
            "sector": 56,
            "tic": 202,
            "period_d": 0.4,
            "t0_bjd": 2459000.2,
            "duration_min": 30.0,
            "prior_label": "eclipsing_binary_or_pceb",
            "twirl_vet_sheet_name": "eb.png",
        },
    ]
    queue = pd.DataFrame(rows)
    queue["candidate_key"] = queue.apply(candidate_key, axis=1)
    queue["observation_candidate_key"] = queue["candidate_key"]
    queue.to_csv(bundle / "review_queue_planet_eb.csv", index=False)
    queue[
        [
            "row_id",
            "review_id",
            "sector",
            "tic",
            "period_d",
            "t0_bjd",
            "duration_min",
            "twirl_vet_sheet_name",
        ]
    ].to_csv(bundle / "required_sheet_manifest.csv", index=False)
    pd.DataFrame(
        columns=[
            "row_id",
            "candidate_key",
            "tic",
            "sector",
            "label",
            "labeler",
        ]
    ).to_csv(
        bundle / "human_labels_final.csv",
        index=False,
    )

    sheets = bundle / "sheets"
    provenance = bundle / "provenance"
    sheets.mkdir()
    provenance.mkdir()
    for name in ("planet.png", "eb.png"):
        Image.new("RGB", (2, 2), "white").save(sheets / name)

    metrics = queue.rename(
        columns={
            "period_d": "anchor_period_d",
            "t0_bjd": "anchor_t0_bjd",
            "duration_min": "anchor_duration_min",
        }
    )[
        [
            "review_id",
            "tic",
            "sector",
            "twirl_vet_sheet_name",
            "anchor_period_d",
            "anchor_t0_bjd",
            "anchor_duration_min",
        ]
    ].copy()
    metrics["twirl_vet_status"] = "ok"
    metrics["anchor_source"] = "row_metadata"
    metrics["vet_branch"] = "current_adp"
    metrics["vet_sheet_version"] = "S56-ADP-HV2"
    metrics.to_csv(provenance / "s0056_metrics.csv", index=False)
    (provenance / "s0056_summary.json").write_text(
        json.dumps(
            {
                "n_rows": 2,
                "status_counts": {"ok": 2},
                "branch_name": "current_adp",
                "apertures": ["DET_FLUX_ADP_SML", "DET_FLUX_ADP"],
                "use_row_ephemeris": True,
                "write_pdf": False,
                "harmonic_factors": [0.25, 0.5, 1.0, 2.0, 4.0],
                "n_periods": 20_000,
                "n_peaks": 10,
                "lc_export_h5": "/pdo/test/s56_current.h5",
            }
        )
    )
    (bundle / "sheet_set.json").write_text(
        json.dumps(
            {
                "sheet_set_version": "test-v1",
                "sheet_root": "sheets",
                "render_provenance_root": "provenance",
                "require_render_provenance": True,
                "branch_name": "current_adp",
                "renderer_version": "S56-ADP-HV2",
                "apertures": ["DET_FLUX_ADP_SML", "DET_FLUX_ADP"],
                "use_row_ephemeris": True,
                "write_pdf": False,
                "harmonic_factors": [0.25, 0.5, 1.0, 2.0, 4.0],
                "n_periods": 20_000,
                "n_peaks": 10,
                "lc_export_h5_by_sector": {
                    "56": "/pdo/test/s56_current.h5"
                },
            }
        )
    )


def test_verifier_accepts_exact_planet_first_render_bundle(tmp_path) -> None:
    _write_bundle(tmp_path)
    verifier = _load_verifier()

    summary = verifier.verify(tmp_path)

    assert summary["assets_complete"]
    assert summary["n_render_provenance_verified"] == 2
    assert summary["n_present_sheets"] == 2
    assert not summary["review_complete"]


def test_verifier_does_not_count_blank_saved_label_as_reviewed(
    tmp_path,
) -> None:
    _write_bundle(tmp_path)
    queue = pd.read_csv(tmp_path / "review_queue_planet_eb.csv")
    pd.DataFrame(
        [
            {
                "row_id": row.row_id,
                "candidate_key": row.candidate_key,
                "tic": row.tic,
                "sector": row.sector,
                "label": "planet_like" if row.row_id == 0 else "",
                "labeler": "tehan",
            }
            for row in queue.itertuples()
        ]
    ).to_csv(tmp_path / "human_labels_final.csv", index=False)
    verifier = _load_verifier()

    summary = verifier.verify(tmp_path)

    assert summary["n_saved_rows"] == 2
    assert summary["n_reviewed_rows"] == 1
    assert summary["n_pending_rows"] == 1
    assert not summary["review_complete"]


def test_verifier_rejects_corrupt_png_and_recipe_drift(tmp_path) -> None:
    _write_bundle(tmp_path)
    (tmp_path / "sheets/planet.png").write_bytes(
        b"\x89PNG\r\n\x1a\ntruncated"
    )
    summary_path = tmp_path / "provenance/s0056_summary.json"
    render_summary = json.loads(summary_path.read_text())
    render_summary["n_periods"] = 100
    summary_path.write_text(json.dumps(render_summary))
    verifier = _load_verifier()

    summary = verifier.verify(tmp_path)

    assert not summary["assets_complete"]
    assert summary["invalid_png_names"] == ["planet.png"]
    assert summary["n_render_provenance_verified"] == 0
    assert summary["render_provenance_errors"] == [
        "S56: renderer summary recipe mismatch"
    ]


def test_verifier_rejects_eb_before_planet(tmp_path) -> None:
    _write_bundle(tmp_path)
    queue_path = tmp_path / "review_queue_planet_eb.csv"
    queue = pd.read_csv(queue_path).iloc[::-1]
    queue.to_csv(queue_path, index=False)
    verifier = _load_verifier()

    with pytest.raises(ValueError, match="Planet-like before EB"):
        verifier.verify(tmp_path)
