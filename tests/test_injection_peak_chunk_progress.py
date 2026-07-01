from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_module():
    path = REPO_ROOT / "scripts" / "stage5_validation" / "summarize_injection_peak_chunk_progress.py"
    spec = importlib.util.spec_from_file_location("summarize_injection_peak_chunk_progress", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_chunk(root: Path, name: str, rows: list[dict[str, object]]) -> None:
    chunk = root / name
    chunk.mkdir(parents=True)
    pd.DataFrame(rows).to_csv(chunk / "injection_bls_peaks.csv", index=False)
    (chunk / "injection_bls_peaks_summary.json").write_text(json.dumps({"rows": len(rows)}) + "\n")


def _peak_row(
    injection_id: str,
    *,
    tmag: float,
    peak_rank: int,
    is_signal: bool,
    match_kind: str,
) -> dict[str, object]:
    return {
        "injection_id": injection_id,
        "tmag": tmag,
        "is_candidate_peak": True,
        "is_injected_signal_peak": is_signal,
        "peak_rank": peak_rank,
        "match_kind": match_kind,
    }


def test_chunk_progress_summarizes_partial_recall_and_tmag_bins(tmp_path: Path) -> None:
    module = _load_module()
    chunk_root = tmp_path / "peak_training" / "chunks"
    manifest_dir = tmp_path / "peak_training" / "chunk_ids"
    manifest_dir.mkdir(parents=True)
    (manifest_dir / "manifest.txt").write_text("n_chunks=3\n")

    _write_chunk(
        chunk_root,
        "chunk_000",
        [
            _peak_row("inj_a", tmag=16.5, peak_rank=1, is_signal=True, match_kind="exact"),
            _peak_row("inj_a", tmag=16.5, peak_rank=2, is_signal=False, match_kind="mismatch"),
            _peak_row("inj_b", tmag=18.5, peak_rank=3, is_signal=True, match_kind="harmonic"),
        ],
    )
    _write_chunk(
        chunk_root,
        "chunk_001",
        [
            _peak_row("inj_c", tmag=19.5, peak_rank=1, is_signal=False, match_kind="mismatch"),
            _peak_row("inj_c", tmag=19.5, peak_rank=2, is_signal=False, match_kind="mismatch"),
        ],
    )
    (chunk_root / "chunk_002").mkdir(parents=True)

    summary = module.summarize_chunk_progress(chunk_root)

    assert summary["expected_chunks"] == 3
    assert summary["n_completed_chunks"] == 2
    assert summary["n_incomplete_chunks"] == 1
    assert summary["n_injections"] == 3
    assert summary["recall"]["top1"]["n"] == 1
    assert summary["recall"]["top5"]["n"] == 2
    assert summary["recall"]["top20"]["n"] == 2
    assert summary["by_tmag"]["<17"]["top20"] == 1
    assert summary["by_tmag"]["18-19"]["top20"] == 1
    assert summary["by_tmag"][">19"]["top20"] == 0
    assert summary["match_kind_counts"]["mismatch"] == 3


def test_chunk_progress_default_repo_root_uses_cwd_for_stdin(tmp_path: Path) -> None:
    module = _load_module()
    old_file = module.__file__
    old_cwd = Path.cwd()
    try:
        module.__file__ = "<stdin>"
        os.chdir(tmp_path)
        assert module._default_repo_root().resolve() == tmp_path.resolve()
    finally:
        module.__file__ = old_file
        os.chdir(old_cwd)


def test_chunk_progress_writes_json_and_markdown(tmp_path: Path) -> None:
    module = _load_module()
    chunk_root = tmp_path / "chunks"
    out_dir = tmp_path / "out"
    _write_chunk(
        chunk_root,
        "chunk_000",
        [_peak_row("inj_a", tmag=16.5, peak_rank=1, is_signal=True, match_kind="exact")],
    )

    code = module.main(["--chunk-root", str(chunk_root), "--out-dir", str(out_dir)])

    assert code == 0
    assert (out_dir / "summary.json").exists()
    markdown = (out_dir / "summary.md").read_text()
    assert "top1" in markdown
    assert "1/1" in markdown
