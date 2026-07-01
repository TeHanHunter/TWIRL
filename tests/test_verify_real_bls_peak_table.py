from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_verifier():
    path = REPO_ROOT / "scripts" / "stage5_validation" / "verify_real_bls_peak_table.py"
    spec = importlib.util.spec_from_file_location("verify_real_bls_peak_table", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _real_peak_rows(
    n_tics: int = 4,
    n_apertures: int = 2,
    n_peaks: int = 5,
    include_provenance: bool = False,
) -> pd.DataFrame:
    apertures = ["DET_FLUX_SML", "DET_FLUX"][:n_apertures]
    rows = []
    for tic_idx in range(n_tics):
        for aperture in apertures:
            for peak_rank in range(1, n_peaks + 1):
                row = {
                    "tic": 1000 + tic_idx,
                    "sector": 56,
                    "cam": 1,
                    "ccd": 1,
                    "tmag": 17.0,
                    "aperture": aperture,
                    "peak_rank": peak_rank,
                    "period_d": 0.2 + 0.1 * peak_rank,
                    "t0_bjd": 2459825.0 + 0.01 * peak_rank,
                    "duration_min": 5.0 + peak_rank,
                    "depth": 0.01,
                    "depth_snr": 5.0,
                    "sde": 12.0 - peak_rank * 0.1,
                    "status": "ok",
                }
                if include_provenance:
                    row.update(
                        {
                            "bls_search_branch": "short_pmax2_quota",
                            "bls_n_periods": 100000,
                            "bls_n_peaks": 20,
                            "bls_p_min_d": 0.12,
                            "bls_p_max_cap_d": 2.0,
                            "bls_max_period_fraction": 0.45,
                            "bls_period_mask_frac": 0.005,
                            "bls_period_bin_edges": "[0.12, 0.25, 0.5, 1.0, 2.0]",
                            "bls_max_peaks_per_period_bin": 5,
                            "bls_sigma_clip": 5.0,
                            "bls_orbit_edge_trim_d": 0.0,
                            "bls_durations_min": "[3.0, 4.0, 5.0]",
                        }
                    )
                rows.append(row)
    return pd.DataFrame(rows)


def test_real_bls_peak_table_verifier_accepts_multi_peak_table(tmp_path: Path) -> None:
    module = _load_verifier()
    path = tmp_path / "candidates.csv"
    _real_peak_rows().to_csv(path, index=False)

    result = module.verify_real_peak_table(
        peak_table=path,
        min_rows=40,
        min_tics=4,
        min_peak_ranks=5,
        min_apertures=2,
        min_rows_per_tic_max=10,
    )

    assert result["passed"]
    assert result["n_unique_tic"] == 4
    assert result["n_positive_peak_ranks"] == 5
    assert result["n_apertures"] == 2
    assert result["max_rows_per_tic"] == 10
    assert result["bls_provenance_columns"] == []


def test_real_bls_peak_table_verifier_rejects_one_peak_table(tmp_path: Path) -> None:
    module = _load_verifier()
    path = tmp_path / "candidates.csv"
    _real_peak_rows(n_tics=4, n_apertures=1, n_peaks=1).to_csv(path, index=False)

    result = module.verify_real_peak_table(
        peak_table=path,
        min_rows=4,
        min_tics=4,
        min_peak_ranks=5,
        min_apertures=2,
        min_rows_per_tic_max=10,
    )

    assert not result["passed"]
    assert any("positive peak-rank count" in failure for failure in result["failures"])
    assert any("aperture count" in failure for failure in result["failures"])
    assert any("max rows per TIC" in failure for failure in result["failures"])


def test_real_bls_peak_table_verifier_reports_branch_provenance(tmp_path: Path) -> None:
    module = _load_verifier()
    path = tmp_path / "candidates.csv"
    _real_peak_rows(include_provenance=True).to_csv(path, index=False)

    result = module.verify_real_peak_table(
        peak_table=path,
        min_rows=40,
        min_tics=4,
        min_peak_ranks=5,
        min_apertures=2,
        min_rows_per_tic_max=10,
        require_bls_provenance=True,
    )

    assert result["passed"]
    assert result["bls_search_branch_counts"] == {"short_pmax2_quota": 40}
    assert "bls_n_periods" in result["bls_provenance_numeric_summaries"]


def test_real_bls_peak_table_verifier_can_require_branch_provenance(tmp_path: Path) -> None:
    module = _load_verifier()
    path = tmp_path / "candidates.csv"
    _real_peak_rows().to_csv(path, index=False)

    result = module.verify_real_peak_table(
        peak_table=path,
        min_rows=40,
        min_tics=4,
        min_peak_ranks=5,
        min_apertures=2,
        min_rows_per_tic_max=10,
        require_bls_provenance=True,
    )

    assert not result["passed"]
    assert any("missing BLS provenance columns" in failure for failure in result["failures"])
