from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_verifier():
    path = REPO_ROOT / "scripts" / "stage5_validation" / "verify_ranker_selected_real_queue.py"
    spec = importlib.util.spec_from_file_location("verify_ranker_selected_real_queue", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _ranker_row(tic: int, rank: int = 1) -> dict[str, object]:
    return {
        "review_id": f"real:{tic}:ranker:{rank}",
        "source_kind": "real_candidate",
        "source_bucket": "real_ranker_selected",
        "tic": tic,
        "sector": 56,
        "period_d": 1.25 + 0.1 * rank,
        "t0_bjd": 2459825.1,
        "duration_min": 8.0,
        "rep_aperture": "DET_FLUX_ADP_SML",
        "sde_max": 12.5,
        "recovery_status": "real_ranker_selected",
        "ranker_selection_rank": rank,
        "injection_id": "",
        "truth_period_d": "",
        "truth_t0_bjd": "",
        "truth_radius_rearth": "",
        "truth_source_kind": "real_candidate",
        "leo_class": "",
        "leo_report_name": "",
    }


def test_verify_ranker_queue_passes_skip_leo_preflight(tmp_path: Path) -> None:
    module = _load_verifier()
    queue = pd.DataFrame([_ranker_row(101, 1), _ranker_row(101, 2), _ranker_row(202, 1)])
    queue_path = tmp_path / "review_queue.csv"
    queue.to_csv(queue_path, index=False)
    summary = {
        "real_selection": "ranker_selected",
        "n_review_rows": 3,
        "n_real": 3,
        "n_injections": 0,
        "n_leo_reports_attempted": 0,
        "n_leo_errors": 0,
        "n_leo_plot_errors": 0,
    }
    summary_path = tmp_path / "summary.json"
    summary_path.write_text(json.dumps(summary))

    result = module.verify_ranker_queue(
        queue_path=queue_path,
        summary_path=summary_path,
        min_rows=3,
        expect_real=3,
        expect_injected=0,
        expect_skip_leo=True,
    )

    assert result["passed"]
    assert result["n_rows"] == 3
    assert result["n_unique_tic"] == 2


def test_verify_ranker_queue_fails_on_injection_and_bad_ephemeris(tmp_path: Path) -> None:
    module = _load_verifier()
    rows = [_ranker_row(101, 1)]
    bad = _ranker_row(202, 1)
    bad.update(
        {
            "review_id": "inj:bad",
            "source_kind": "injection_recovery",
            "period_d": float("nan"),
            "duration_min": -1,
            "injection_id": "predet_001",
            "truth_period_d": 2.0,
        }
    )
    rows.append(bad)
    queue_path = tmp_path / "review_queue.csv"
    pd.DataFrame(rows).to_csv(queue_path, index=False)

    result = module.verify_ranker_queue(
        queue_path=queue_path,
        min_rows=2,
        expect_real=1,
        expect_injected=0,
        expect_skip_leo=True,
    )

    assert not result["passed"]
    assert any("injection rows" in failure for failure in result["failures"])
    assert any("non-finite period_d" in failure for failure in result["failures"])
    assert any("duration_min" in failure for failure in result["failures"])


def test_verify_ranker_queue_requires_report_files_when_requested(tmp_path: Path) -> None:
    module = _load_verifier()
    row = _ranker_row(101, 1)
    row["leo_report_name"] = "PC_row00001.pdf"
    queue_path = tmp_path / "review_queue.csv"
    pd.DataFrame([row]).to_csv(queue_path, index=False)
    reports_dir = tmp_path / "vet_reports"
    reports_dir.mkdir()
    (reports_dir / "PC_row00001.pdf").write_bytes(b"%PDF-1.4\n")

    result = module.verify_ranker_queue(
        queue_path=queue_path,
        reports_dir=reports_dir,
        min_rows=1,
        expect_real=1,
        expect_injected=0,
        require_reports=True,
    )

    assert result["passed"]


def test_verify_ranker_queue_rejects_mixed_source_bucket(tmp_path: Path) -> None:
    module = _load_verifier()
    good = _ranker_row(101, 1)
    bad = _ranker_row(202, 1)
    bad["source_bucket"] = "real_random"
    queue_path = tmp_path / "review_queue.csv"
    pd.DataFrame([good, bad]).to_csv(queue_path, index=False)

    result = module.verify_ranker_queue(
        queue_path=queue_path,
        min_rows=2,
        expect_real=2,
        expect_injected=0,
        expect_skip_leo=True,
    )

    assert not result["passed"]
    assert any("source_bucket=real_ranker_selected" in failure for failure in result["failures"])
