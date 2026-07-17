from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from twirl.vetting.franklin_multisector import (
    build_franklin_multisector_batch,
    prepare_single_teacher_scores,
    verify_franklin_multisector_batch,
)
from twirl.vetting.teacher_v2_active_learning import EnrichmentQuotas


ROOT = Path(__file__).resolve().parents[1]


def _load_script(name: str, relative: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _teacher_scores(
    *,
    sector: int,
    tic_start: int,
    n_tics: int = 80,
) -> pd.DataFrame:
    rng = np.random.default_rng(5600 + sector)
    tic = np.repeat(np.arange(tic_start, tic_start + n_tics), 2)
    peak_rank = np.tile([1, 2], n_tics)
    n_rows = len(tic)
    probabilities = rng.dirichlet(np.ones(4), size=n_rows)
    frame = pd.DataFrame(
        {
            "review_id": [
                f"s{sector}-tic{value}-rank{rank}"
                for value, rank in zip(tic, peak_rank)
            ],
            "tic": tic,
            "sector": sector,
            "period_d": 0.2 + rng.random(n_rows) * 8.0,
            "t0_bjd": 2_459_000.0 + rng.random(n_rows),
            "duration_min": 3.0 + rng.random(n_rows) * 20.0,
            "tmag": 15.0 + rng.random(n_rows) * 7.0,
            "sde_max": 3.0 + rng.random(n_rows) * 30.0,
            "rep_peak_rank": peak_rank,
            "cam": 1,
            "ccd": 2,
            "source_kind": "real_candidate",
            "p_preserve": rng.random(n_rows),
            "std_p_preserve": rng.random(n_rows) * 0.05,
            "model_profile": "shape_plus_periodogram_bls",
            "model_version": "s56_harmonic_cnn_v1",
        }
    )
    for index, label in enumerate(
        ("planet_like", "eclipse_contact", "smooth_variable", "other")
    ):
        frame[f"p_{label}"] = probabilities[:, index]
        frame[f"std_p_{label}"] = rng.random(n_rows) * 0.05
    return frame


def test_single_teacher_adapter_derives_compact_score_and_rejects_injections() -> None:
    scores = _teacher_scores(sector=56, tic_start=10_000, n_tics=5)
    adapted = prepare_single_teacher_scores(scores)
    expected = np.sqrt(
        scores["p_planet_like"].to_numpy() * scores["p_preserve"].to_numpy()
    )
    assert np.allclose(adapted["p_compact_transit"], expected)
    assert adapted["morphology_entropy"].between(0.0, 1.0).all()

    injected = scores.copy()
    injected.loc[0, "source_kind"] = "injected_signal"
    with pytest.raises(ValueError, match="injected"):
        prepare_single_teacher_scores(injected)


def test_multisector_batch_uses_rank1_unequal_quotas_and_global_tic_dedup() -> None:
    scores = {
        56: _teacher_scores(sector=56, tic_start=10_000),
        57: _teacher_scores(sector=57, tic_start=10_060),
        58: _teacher_scores(sector=58, tic_start=20_000),
        59: _teacher_scores(sector=59, tic_start=20_060),
    }
    small = EnrichmentQuotas(4, 3, 2, 1, 1)
    large = EnrichmentQuotas(8, 6, 4, 2, 2)
    quotas = {56: small, 57: small, 58: large, 59: large}
    excluded = pd.DataFrame({"tic": [10_000, 20_000]})
    queue, overlap, hidden, summary = build_franklin_multisector_batch(
        scores,
        sector_quotas=quotas,
        excluded_tables=[excluded],
        double_review_count=4,
    )
    assert len(queue) == 66
    assert queue["tic"].nunique() == 66
    assert queue["sector"].value_counts().to_dict() == {
        58: 22,
        59: 22,
        56: 11,
        57: 11,
    }
    assert queue["rep_peak_rank"].eq(1).all()
    assert set(queue["tic"]).isdisjoint({10_000, 20_000})
    assert len(overlap) == 4
    assert summary["rank_policy"].endswith("rank 1 only")
    assert not any(
        column.startswith(("p_", "std_p_", "selection_", "model_"))
        for column in queue
    )
    assert verify_franklin_multisector_batch(
        queue,
        overlap,
        hidden,
        sector_quotas=quotas,
        expected_overlap_count=4,
        excluded_tics={10_000, 20_000},
    )["passed"]


def test_handoff_builder_and_app_roundtrip_period_factor(tmp_path: Path) -> None:
    scores = {56: _teacher_scores(sector=56, tic_start=30_000, n_tics=20)}
    quotas = {56: EnrichmentQuotas(2, 2, 1, 1, 1)}
    queue, _, _, _ = build_franklin_multisector_batch(
        scores,
        sector_quotas=quotas,
    )
    source_queue = tmp_path / "queue.csv"
    queue.to_csv(source_queue, index=False)
    source_sheets = tmp_path / "source_sheets"
    source_sheets.mkdir()
    for name in queue["twirl_vet_sheet_name"]:
        (source_sheets / name).write_bytes(b"\x89PNG\r\n\x1a\n")
    reference_root = tmp_path / "references"
    reference_root.mkdir()
    (reference_root / "example.png").write_bytes(b"\x89PNG\r\n\x1a\n")
    reference_csv = tmp_path / "reference_examples.csv"
    pd.DataFrame({"label": ["planet_like"], "example_sheet": ["example.png"]}).to_csv(
        reference_csv, index=False
    )

    handoff_module = _load_script(
        "build_franklin_multisector_handoff",
        "scripts/stage5_validation/build_franklin_multisector_handoff.py",
    )
    out_dir = tmp_path / "handoff"
    summary = handoff_module.build_handoff(
        queue_path=source_queue,
        sheet_root=source_sheets,
        out_dir=out_dir,
        app_source=ROOT / "scripts/stage5_validation/franklin_vetting_app.py",
        expected_sector_counts={56: 7},
        reference_root=reference_root,
        reference_csv=reference_csv,
        port=5003,
    )
    assert summary["n_rows"] == 7
    assert summary["scores_in_package"] is False

    app_module = _load_script(
        "franklin_vetting_app",
        "scripts/stage5_validation/franklin_vetting_app.py",
    )
    labels = out_dir / "franklin_labels_vetted.csv"
    store = app_module.Store(
        out_dir / "franklin_review_queue_7_real.csv",
        labels,
        out_dir / "vet_sheets",
        "franklin",
    )
    record = store.save_label(
        0,
        "eclipsing_binary_or_pceb",
        "franklin",
        "best at twice the displayed period",
        period_factor="2",
    )
    assert record["period_factor"] == "2"
    assert record["period_status"] == "resolved"
    reloaded = app_module.Store(
        out_dir / "franklin_review_queue_7_real.csv",
        labels,
        out_dir / "vet_sheets",
        "franklin",
    )
    assert reloaded.payload(0)["period_factor"] == "2"

    broken = pd.read_csv(labels)
    broken.loc[0, "candidate_key"] = "wrong"
    broken.to_csv(labels, index=False)
    with pytest.raises(ValueError, match="candidate_key mismatch"):
        app_module.Store(
            out_dir / "franklin_review_queue_7_real.csv",
            labels,
            out_dir / "vet_sheets",
            "franklin",
        )
