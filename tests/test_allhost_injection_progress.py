from __future__ import annotations

import importlib.util
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_module(name: str, relative: str):
    path = REPO_ROOT / relative
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_summarize_allhost_progress_counts_completed_and_active(tmp_path: Path) -> None:
    module = _load_module(
        "summarize_s56_allhost_injection_progress",
        "scripts/stage5_validation/summarize_s56_allhost_injection_progress.py",
    )
    out_root = tmp_path / "allhost"
    shard_dir = out_root / "target_shards"
    chunk_root = out_root / "chunks"
    shard_dir.mkdir(parents=True)
    (shard_dir / "manifest.json").write_text(
        json.dumps(
            {
                "n_targets": 3,
                "n_paths": 6,
                "n_shards": 3,
                "shards": [{"shard": "chunk_000"}, {"shard": "chunk_001"}, {"shard": "chunk_002"}],
            }
        )
    )
    chunk0 = chunk_root / "chunk_000"
    chunk1 = chunk_root / "chunk_001"
    chunk0.mkdir(parents=True)
    chunk1.mkdir(parents=True)
    (chunk0 / "summary.json").write_text(
        json.dumps(
            {
                "n_injections": 2,
                "n_unique_injected_tics": 2,
                "skipped": {"read_failed": 1, "model_failed": 0},
            }
        )
    )
    (chunk0 / "injected_lightcurves.h5").write_bytes(b"not-real-hdf5")

    summary = module.summarize_progress(out_root)

    assert summary["n_shards"] == 3
    assert summary["n_started_shards"] == 2
    assert summary["n_completed_shards"] == 1
    assert summary["n_completed_h5_shards"] == 1
    assert summary["n_active_or_incomplete_shards"] == 1
    assert summary["n_pending_shards"] == 1
    assert summary["n_injections_completed"] == 2
    assert summary["n_unique_injected_tics_completed"] == 2
    assert summary["skipped_totals_completed"]["read_failed"] == 1
    assert summary["active_or_incomplete_shards"] == ["chunk_001"]
