from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def _text(name: str) -> str:
    return (ROOT / "scripts" / "orcd" / name).read_text(encoding="utf-8")


def test_v4_bls_rebuilds_uniformly_from_all_112_raw_v1_shards() -> None:
    text = _text("slurm_teacher_ssl_fullpool_v4_bls_cpu.sbatch")

    assert "#SBATCH --array=0-111%28" in text
    assert "N_SHARDS=16" in text
    assert "raw_source_${SHARD_INDEX}.h5" in text
    assert "--raw-source-h5" in text
    assert "--raw-source-summary" in text
    assert "--raw-export-complete" in text
    assert "--raw-transfer-validation" in text
    assert "--frozen-pool" in text
    assert "--n-periods 50000" in text
    assert "--n-peaks 10" in text
    assert "A2v1-fullpool-v1" in text
    assert "TWIRL_SSL_FULLPOOL_EXECUTION_ALLOWLIST_ROOT" in text
    assert "--execution-allowlist" in text


def test_v4_bls_merges_exactly_16_shards_per_sector() -> None:
    text = _text("slurm_teacher_ssl_fullpool_v4_bls_merge_cpu.sbatch")

    assert "#SBATCH --array=0-6%7" in text
    assert "N_SHARDS=16" in text
    assert "--n-shards \"${N_SHARDS}\"" in text


def test_v4_jobs_bind_the_exact_deployed_commit() -> None:
    for name in (
        "slurm_teacher_ssl_fullpool_v4_bls_cpu.sbatch",
        "slurm_teacher_ssl_fullpool_v4_bls_merge_cpu.sbatch",
        "slurm_teacher_ssl_fullpool_v4_bls_global_cpu.sbatch",
    ):
        text = _text(name)
        assert "TWIRL_EXPECTED_GIT_SHA" in text
        assert 'git -C "${REPO}" rev-parse HEAD' in text


def test_v4_canary_is_the_exact_34_affected_observation_set() -> None:
    root = ROOT / "configs" / "teacher_ssl_fullpool_v4_canary"
    by_sector = {
        int(path.stem[1:3]): pd.read_csv(path)["tic"].astype(int).tolist()
        for path in sorted(root.glob("s*_affected.csv"))
    }

    assert set(by_sector) == {56, 58, 61, 62}
    assert {sector: len(tics) for sector, tics in by_sector.items()} == {
        56: 7,
        58: 3,
        61: 16,
        62: 8,
    }
    assert len({(sector, tic) for sector, tics in by_sector.items() for tic in tics}) == 34
