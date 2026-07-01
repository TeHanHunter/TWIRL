from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_verifier():
    path = REPO_ROOT / "scripts" / "stage5_validation" / "verify_injection_peak_training_table.py"
    spec = importlib.util.spec_from_file_location("verify_injection_peak_training_table", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _peak_rows(include_cadence: bool = True) -> pd.DataFrame:
    rows = []
    for inj_idx in range(4):
        injection_id = f"inj_{inj_idx:04d}"
        for aperture in ("DET_FLUX_ADP_SML", "DET_FLUX_SML"):
            for peak_rank in range(1, 4):
                is_signal = peak_rank == 2
                row = {
                    "injection_id": injection_id,
                    "tic": 100000 + inj_idx,
                    "aperture": aperture,
                    "status": "ok",
                    "is_candidate_peak": True,
                    "peak_rank": peak_rank,
                    "is_injected_signal_peak": is_signal,
                    "exact_ephemeris_match": is_signal,
                    "harmonic_ephemeris_match": False,
                    "match_kind": "exact" if is_signal else "mismatch",
                    "truth_period_d": 1.0,
                    "truth_t0_bjd": 2459825.0,
                    "truth_duration_min": 8.0,
                    "period_d": 1.0 if is_signal else 0.5 + peak_rank,
                    "t0_bjd": 2459825.0,
                    "duration_min": 8.0,
                    "sde": 20.0 - peak_rank,
                }
                if include_cadence:
                    row.update(
                        {
                            "n_cad_total": 11000,
                            "n_cad_quality": 10500,
                            "n_cad_kept": 10400,
                            "n_cad_edge_trimmed": 50,
                            "n_cad_sigma_clipped": 50,
                            "dropout_frac": 0.05,
                            "quality_dropout_frac": 0.04,
                        }
                    )
                row.update(
                    {
                        "transit_window_match": is_signal,
                        "transit_window_phase_period_d": 1.0,
                        "transit_window_delta_min": 0.0 if is_signal else 120.0,
                        "transit_window_overlap_min": 8.0 if is_signal else 0.0,
                        "transit_window_overlap_fraction": 1.0 if is_signal else 0.0,
                    }
                )
                rows.append(row)
    return pd.DataFrame(rows)


def test_injection_peak_table_verifier_accepts_cadence_diagnostics(tmp_path: Path) -> None:
    module = _load_verifier()
    path = tmp_path / "peaks.csv"
    _peak_rows(include_cadence=True).to_csv(path, index=False)

    result = module.verify_injection_peak_table(
        peak_table=path,
        min_injections=4,
        min_candidate_rows=20,
        min_apertures=2,
        min_positive_peak_ranks=3,
        require_cadence_diagnostics=True,
    )

    assert result["passed"]
    assert result["n_injections"] == 4
    assert result["n_signal_peak_rows"] == 8
    assert result["n_apertures"] == 2
    assert "n_cad_quality" in result["cadence_diagnostic_columns"]
    assert "transit_window_overlap_fraction" in result["transit_window_overlap_columns"]


def test_injection_peak_table_verifier_rejects_missing_cadence_diagnostics(tmp_path: Path) -> None:
    module = _load_verifier()
    path = tmp_path / "peaks.csv"
    _peak_rows(include_cadence=False).to_csv(path, index=False)

    result = module.verify_injection_peak_table(
        peak_table=path,
        min_injections=4,
        min_candidate_rows=20,
        min_apertures=2,
        min_positive_peak_ranks=3,
        require_cadence_diagnostics=True,
    )

    assert not result["passed"]
    assert any("missing cadence diagnostic columns" in failure for failure in result["failures"])


def test_injection_peak_table_verifier_can_skip_signal_requirement_for_smoke(tmp_path: Path) -> None:
    module = _load_verifier()
    path = tmp_path / "peaks.csv"
    rows = _peak_rows(include_cadence=True)
    rows["is_injected_signal_peak"] = False
    rows["exact_ephemeris_match"] = False
    rows["match_kind"] = "mismatch"
    rows.to_csv(path, index=False)

    result = module.verify_injection_peak_table(
        peak_table=path,
        min_injections=4,
        min_candidate_rows=20,
        min_apertures=2,
        min_positive_peak_ranks=3,
        require_cadence_diagnostics=True,
        require_signal_peaks=False,
    )

    assert result["passed"]
    assert result["n_signal_peak_rows"] == 0


def test_injection_peak_table_verifier_rejects_signal_without_window_overlap(tmp_path: Path) -> None:
    module = _load_verifier()
    path = tmp_path / "peaks.csv"
    rows = _peak_rows(include_cadence=True)
    signal = rows["is_injected_signal_peak"].astype(bool)
    rows.loc[signal, "transit_window_match"] = False
    rows.to_csv(path, index=False)

    result = module.verify_injection_peak_table(
        peak_table=path,
        min_injections=4,
        min_candidate_rows=20,
        min_apertures=2,
        min_positive_peak_ranks=3,
        require_cadence_diagnostics=True,
    )

    assert not result["passed"]
    assert any("signal peak rows lack transit-window overlap" in failure for failure in result["failures"])
