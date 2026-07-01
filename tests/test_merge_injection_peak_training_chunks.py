from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_merger():
    path = REPO_ROOT / "scripts" / "stage5_validation" / "merge_injection_peak_training_chunks.py"
    spec = importlib.util.spec_from_file_location("merge_injection_peak_training_chunks", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _chunk_frame(branch: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "search_branch": [branch, branch],
            "injection_id": ["inj_0001", "inj_0001"],
            "is_candidate_peak": [True, True],
            "peak_rank": [1, 2],
            "is_injected_signal_peak": [False, True],
            "exact_ephemeris_match": [False, True],
            "harmonic_ephemeris_match": [False, False],
            "match_kind": ["mismatch", "exact"],
        }
    )


def test_merge_peak_chunks_allows_duplicate_injections_when_branch_is_key(tmp_path: Path) -> None:
    module = _load_merger()
    chunk_root = tmp_path / "chunks"
    for idx, branch in enumerate(("standard", "short_pmax2")):
        chunk_dir = chunk_root / f"chunk_{idx:03d}"
        chunk_dir.mkdir(parents=True)
        _chunk_frame(branch).to_csv(chunk_dir / "injection_bls_peaks.csv", index=False)

    summary = module.merge_peak_chunks(
        chunk_root=chunk_root,
        out_table=tmp_path / "merged.csv",
        table_name="injection_bls_peaks.csv",
        expect_injections=1,
        id_columns=("search_branch", "injection_id"),
    )
    merged = pd.read_csv(tmp_path / "merged.csv")

    assert summary["n_unique_ids"] == 2
    assert summary["n_injections"] == 1
    assert set(merged["search_branch"]) == {"standard", "short_pmax2"}


def test_merge_peak_chunks_rejects_duplicate_injections_with_default_key(tmp_path: Path) -> None:
    module = _load_merger()
    chunk_root = tmp_path / "chunks"
    for idx, branch in enumerate(("standard", "short_pmax2")):
        chunk_dir = chunk_root / f"chunk_{idx:03d}"
        chunk_dir.mkdir(parents=True)
        _chunk_frame(branch).to_csv(chunk_dir / "injection_bls_peaks.csv", index=False)

    try:
        module.merge_peak_chunks(
            chunk_root=chunk_root,
            out_table=tmp_path / "merged.csv",
            table_name="injection_bls_peaks.csv",
            expect_injections=0,
        )
    except ValueError as exc:
        assert "duplicated IDs" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("expected duplicate injection IDs to fail with default merge key")
