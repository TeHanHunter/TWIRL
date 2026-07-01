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


def test_target_shards_group_orbit_paths_by_tic(tmp_path: Path) -> None:
    module = _load_module(
        "build_s56_predetrend_target_shards",
        "scripts/stage5_validation/build_s56_predetrend_target_shards.py",
    )
    paths = (
        tmp_path / "orbit119" / "cam4" / "ccd1" / "LC" / "100.h5",
        tmp_path / "orbit120" / "cam4" / "ccd1" / "LC" / "100.h5",
        tmp_path / "orbit119" / "cam4" / "ccd1" / "LC" / "200.h5",
    )

    targets = module.discover_raw_h5_targets(orbit_roots=(), target_h5_paths=paths)
    manifest = module.write_target_shards(targets=targets, out_dir=tmp_path / "shards", shard_size=1)

    assert manifest["n_targets"] == 2
    assert manifest["n_paths"] == 3
    assert manifest["n_shards"] == 2
    first = json.loads((tmp_path / "shards" / "chunk_000.meta.json").read_text())
    assert first["n_tics"] == 1
    assert first["n_paths"] == 2
    assert first["start_index"] == 0
    assert (tmp_path / "shards" / "chunk_000.txt").read_text().count("100.h5") == 2


def test_merge_shard_metadata_combines_manifests(tmp_path: Path) -> None:
    module = _load_module(
        "merge_s56_predetrend_injection_shards",
        "scripts/stage5_validation/merge_s56_predetrend_injection_shards.py",
    )
    chunk_root = tmp_path / "chunks"
    for idx, tic in enumerate((100, 200)):
        chunk = chunk_root / f"chunk_{idx:03d}"
        chunk.mkdir(parents=True)
        (chunk / "injection_manifest.csv").write_text(
            "injection_id,tic,split\n"
            f"predet_{idx:06d},{tic},train\n"
        )
        (chunk / "injection_labels.csv").write_text(
            "injection_id,tic,label\n"
            f"predet_{idx:06d},{tic},planet_like\n"
        )
        (chunk / "summary.json").write_text(
            json.dumps(
                {
                    "n_source_targets": 1,
                    "n_injections": 1,
                    "n_unique_injected_tics": 1,
                    "all_source_targets_injected": True,
                    "skipped": {"read_failed": 0},
                }
            )
        )

    summary = module.merge_shard_metadata(chunk_root=chunk_root, out_dir=tmp_path / "merged")

    assert summary["n_completed_shards"] == 2
    assert summary["n_source_targets"] == 2
    assert summary["n_injections"] == 2
    assert summary["n_unique_injected_tics"] == 2
    assert summary["all_source_targets_injected"] is True
    merged = (tmp_path / "merged" / "injection_manifest.csv").read_text()
    assert "source_chunk" in merged
    assert "source_h5" in merged
