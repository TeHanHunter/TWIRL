from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = (
    REPO_ROOT
    / "scripts"
    / "stage5_validation"
    / "merge_adp_real_bls_shards.py"
)


def _load_module():
    spec = importlib.util.spec_from_file_location(
        "merge_adp_real_bls_shards_test", SCRIPT_PATH
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_shards(tmp_path: Path) -> tuple[object, Path]:
    module = _load_module()
    shard_dir = tmp_path / "shards"
    shard_dir.mkdir()
    apertures = list(module.ADP_ONLY_APERTURES)
    for shard_index, tic in enumerate((1001, 1002)):
        table_path = shard_dir / f"real_adp_bls_peaks_{shard_index:03d}.parquet"
        frame = pd.DataFrame(
            {
                "tic": [tic, tic],
                "aperture": apertures,
                "peak_rank": [1, 1],
                "status": ["ok", "ok"],
                "period_d": [0.5 + shard_index, 0.5 + shard_index],
                "n_cad_total": [100, 100],
                "n_cad_internal_bad": [1, 1],
                "n_cad_external_bad": [2, 2],
                "n_cad_external_only_bad": [2, 2],
                "n_cad_effective_bad": [3, 3],
                "n_cad_authority_excluded": [shard_index, shard_index],
            }
        )
        frame.to_parquet(table_path, index=False)
        summary = {
            "sector": 56,
            "contract_version": "s56_adp_only_v1",
            "bls_search_contract_version": "s56_a2v1_real_bls_v1",
            "bls_config_sha256": "a" * 64,
            "external_quality_policy_contract": "s56_external_quality_v1",
            "compact_lc": "/immutable/compact.h5",
            "compact_lc_sha256": "b" * 64,
            "cadence_reference": "/immutable/cadence.csv",
            "cadence_reference_sha256": "c" * 64,
            "cadence_reference_manifest": "/immutable/cadence.json",
            "cadence_reference_manifest_sha256": "d" * 64,
            "cadence_reference_contract_version": "s56_a2v1_cadence_reference_v1",
            "cadence_reference_cadence_authority": "qlp_cam_quat",
            "cadence_reference_quality_authority": "spoc_and_qlp_quality_flags",
            "cadence_reference_source_hashes_sha256": "e" * 64,
            "authority_exclusion_policy_contract": (
                "a2v1_quat_absent_spoc_cadence_exclusions_v1"
            ),
            "authority_exclusion_external_bit": 62,
            "authority_exclusions_sha256": "f" * 64,
            "n_authority_exclusions": 4,
            "apertures": apertures,
            "n_targets_total": 2,
            "n_periods": 50000,
            "n_peaks": 10,
            "source_product_tag": "A2v1",
            "config": {
                "p_min_d": 0.12,
                "p_max_cap_d": 15.0,
                "max_period_fraction": 0.45,
            },
            "shard_index": shard_index,
            "n_shards": 2,
            "n_targets": 1,
            "n_rows": len(frame),
            "peak_table_sha256": _sha256(table_path),
            "outputs": {"peak_table": str(table_path.resolve())},
        }
        (shard_dir / f"summary_{shard_index:03d}.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    return module, shard_dir


def test_merge_revalidates_provenance_and_emits_validator_summary(
    tmp_path: Path,
) -> None:
    module, shard_dir = _write_shards(tmp_path)
    out_path = tmp_path / "merged" / "real_adp_bls_peaks.parquet"

    summary = module.merge_shards(
        shard_dir=shard_dir,
        out_path=out_path,
        n_shards=2,
    )

    merged = pd.read_parquet(out_path)
    assert len(merged) == 4
    assert set(merged["tic"].astype(int)) == {1001, 1002}
    assert summary["n_shards"] == 1
    assert summary["shard_index"] == 0
    assert summary["n_source_shards"] == 2
    assert summary["n_targets"] == summary["n_targets_total"] == 2
    assert summary["peak_table_sha256"] == _sha256(out_path)
    assert summary["authority_exclusion_external_bit"] == 62
    assert summary["n_authority_exclusions"] == 4
    assert summary["quality_counts_over_unique_targets"] == {
        "n_cad_internal_bad": 2,
        "n_cad_external_bad": 4,
        "n_cad_external_only_bad": 4,
        "n_cad_effective_bad": 6,
        "n_cad_authority_excluded": 1,
    }
    assert summary["outputs"]["peak_table"] == str(out_path)
    assert len(summary["source_shards"]) == 2
    persisted = json.loads(
        out_path.with_suffix(".summary.json").read_text(encoding="utf-8")
    )
    assert persisted["peak_table_sha256"] == _sha256(out_path)
    assert persisted["passed"] is True


def test_merge_rejects_shard_provenance_disagreement(tmp_path: Path) -> None:
    module, shard_dir = _write_shards(tmp_path)
    summary_path = shard_dir / "summary_001.json"
    payload = json.loads(summary_path.read_text(encoding="utf-8"))
    payload["compact_lc_sha256"] = "f" * 64
    summary_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    with pytest.raises(ValueError, match="disagrees on provenance field"):
        module.merge_shards(
            shard_dir=shard_dir,
            out_path=tmp_path / "merged.parquet",
            n_shards=2,
        )


def test_merge_does_not_publish_before_quality_validation(tmp_path: Path) -> None:
    module, shard_dir = _write_shards(tmp_path)
    shard_path = shard_dir / "real_adp_bls_peaks_001.parquet"
    frame = pd.read_parquet(shard_path)
    frame.loc[0, "n_cad_effective_bad"] = 99
    frame.to_parquet(shard_path, index=False)
    summary_path = shard_dir / "summary_001.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["peak_table_sha256"] = _sha256(shard_path)
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    out_path = tmp_path / "merged.parquet"
    accepted = pd.DataFrame({"tic": [9999], "status": ["accepted"]})
    accepted.to_parquet(out_path, index=False)
    accepted_sha256 = _sha256(out_path)

    with pytest.raises(ValueError, match="counts disagree across target rows"):
        module.merge_shards(
            shard_dir=shard_dir,
            out_path=out_path,
            n_shards=2,
        )

    assert _sha256(out_path) == accepted_sha256


def test_merge_rejects_shard_table_hash_mismatch(tmp_path: Path) -> None:
    module, shard_dir = _write_shards(tmp_path)
    summary_path = shard_dir / "summary_001.json"
    payload = json.loads(summary_path.read_text(encoding="utf-8"))
    payload["peak_table_sha256"] = "0" * 64
    summary_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    with pytest.raises(ValueError, match="peak-table SHA256 mismatch"):
        module.merge_shards(
            shard_dir=shard_dir,
            out_path=tmp_path / "merged.parquet",
            n_shards=2,
        )
