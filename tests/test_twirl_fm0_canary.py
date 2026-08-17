from __future__ import annotations

from pathlib import Path

from twirl.models.fm0.canary import (
    WD1856_GAIA_DR3_SOURCE_ID,
    WD1856_TIC_ID,
    select_s56_canary,
    write_s56_canary_selection,
)
from twirl.models.fm0.registry import build_alias_registry


def _rows() -> tuple[list[dict[str, object]], list[dict[str, str]]]:
    observations: list[dict[str, object]] = []
    aliases: list[dict[str, str]] = []
    counter = 1000
    for camera in range(1, 5):
        for ccd in range(1, 5):
            for _ in range(3):
                gaia, tic = str(counter), str(counter + 100000)
                aliases.append({"gaia_dr3_source_id": gaia, "tic_id": tic})
                for orbit in (119, 120):
                    observations.append(
                        {
                            "source_id": gaia,
                            "tic_id": tic,
                            "sector": 56,
                            "orbit": orbit,
                            "camera": camera,
                            "ccd": ccd,
                            "edge_warn": False,
                        }
                    )
                counter += 1
    aliases.append(
        {"gaia_dr3_source_id": WD1856_GAIA_DR3_SOURCE_ID, "tic_id": WD1856_TIC_ID}
    )
    for orbit in (119, 120):
        observations.append(
            {
                "source_id": WD1856_GAIA_DR3_SOURCE_ID,
                "tic_id": WD1856_TIC_ID,
                "sector": 56,
                "orbit": orbit,
                "camera": 4,
                "ccd": 1,
                "edge_warn": False,
            }
        )
    return observations, aliases


def test_canary_is_balanced_deterministic_and_includes_benchmark(tmp_path: Path) -> None:
    observations, aliases = _rows()
    observations.append(
        {
            "source_id": "999999",
            "tic_id": "--",
            "sector": 56,
            "orbit": 119,
            "camera": 1,
            "ccd": 1,
            "edge_warn": False,
        }
    )
    registry = build_alias_registry(aliases)
    first, first_aliases = select_s56_canary(
        observations, registry, pdo_root="/pdo/a2v1", per_detector=2
    )
    second, second_aliases = select_s56_canary(
        reversed(observations), registry, pdo_root="/pdo/a2v1", per_detector=2
    )
    assert first == second
    assert first_aliases == second_aliases
    assert len(first) in {32, 33}  # Benchmark may already be in the detector sample.
    assert sum(bool(row["is_benchmark"]) for row in first) == 1
    assert all(row["orbits_json"] == "[119,120]" for row in first)
    assert all("orbit-119" in row["hdf5_paths_json"] for row in first)
    assert all("orbit-120" in row["hdf5_paths_json"] for row in first)

    summary = write_s56_canary_selection(tmp_path, first, first_aliases)
    assert summary["n_sources"] == len(first)
    assert summary["n_benchmark_sources"] == 1
    assert write_s56_canary_selection(tmp_path, first, first_aliases) == summary
