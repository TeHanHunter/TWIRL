from __future__ import annotations

import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest

from twirl.vetting.teacher_v3_injection_pairs import (
    CANONICAL_REQUIRED_DATASETS,
    PAIR_REQUIRED_DATASETS,
    TEACHER_V3_INJECTION_PAIR_MERGE_CONTRACT,
    promote_teacher_v3_injection_pairs,
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_pair_group(
    h5: h5py.File,
    injection_id: str,
    *,
    tic: int,
    n: int = 8,
) -> None:
    group = h5.create_group(f"injections/{injection_id}")
    group.attrs["injection_id"] = injection_id
    group.attrs["tic"] = tic
    group.attrs["sector"] = 56
    group.attrs["camera"] = 1
    group.attrs["ccd"] = 2
    group.attrs["cadence_s"] = 200.0
    group.attrs["injection_baseline_Small"] = 100.0
    group.attrs["injection_baseline_Primary"] = 200.0
    group.attrs["source_injection_h5"] = "/obsolete/raw_injections.h5"
    group.attrs["custom_attr"] = "preserve-me"
    values = {
        "time": 2825.0 + np.arange(n) * 200.0 / 86400.0,
        "cadenceno": np.arange(1000, 1000 + n, dtype=np.int64),
        "orbitid": np.full(n, 119, dtype=np.int32),
        "quality": np.zeros(n, dtype=np.int32),
        "transit_model": np.linspace(1.0, 0.9, n),
        "DET_FLUX_ADP_SML_injected": np.linspace(1.0, 0.9, n),
        "DET_FLUX_ADP_injected": np.linspace(1.0, 0.92, n),
        "DET_FLUX_ADP_SML_original": np.ones(n),
        "DET_FLUX_ADP_original": np.ones(n),
    }
    assert set(values) == set(PAIR_REQUIRED_DATASETS)
    for name, value in values.items():
        group.create_dataset(name, data=value, compression="gzip")
    group.create_dataset("extra_dataset", data=np.arange(n))


def _write_canonical_group(
    h5: h5py.File,
    injection_id: str,
    *,
    tic: int,
    n: int = 8,
) -> None:
    group = h5.create_group(f"injections/{injection_id}")
    group.attrs["injection_id"] = injection_id
    group.attrs["tic"] = tic
    values = {
        "RAW_FLUX_Small_injected": np.full(n, 99.0),
        "RAW_FLUX_Small_original": np.full(n, 100.0),
        "RAW_FLUX_Primary_injected": np.full(n, 198.0),
        "RAW_FLUX_Primary_original": np.full(n, 200.0),
    }
    assert set(values) == set(CANONICAL_REQUIRED_DATASETS)
    for name, value in values.items():
        group.create_dataset(name, data=value)
    group.create_dataset(
        "time", data=2825.0 + np.arange(n) * 200.0 / 86400.0
    )
    group.create_dataset(
        "cadenceno", data=np.arange(1000, 1000 + n, dtype=np.int64)
    )
    group.create_dataset("orbitid", data=np.full(n, 119, dtype=np.int32))
    group.create_dataset("transit_model", data=np.linspace(1.0, 0.9, n))


def _write_fixture(tmp_path: Path) -> dict[str, Path]:
    repo = tmp_path / "repo"
    pairs = repo / "pairs"
    pairs.mkdir(parents=True)
    pair_a = pairs / "pair_a.h5"
    pair_b = pairs / "pair_b.h5"
    canonical = tmp_path / "canonical.h5"
    with h5py.File(pair_a, "w") as h5:
        h5.attrs["source_injection_h5"] = "/obsolete/raw_a.h5"
        _write_pair_group(h5, "inj-a", tic=101)
        _write_pair_group(h5, "not-selected", tic=999)
    with h5py.File(pair_b, "w") as h5:
        h5.attrs["source_injection_h5"] = "/obsolete/raw_b.h5"
        _write_pair_group(h5, "inj-b", tic=202)
    with h5py.File(canonical, "w") as h5:
        _write_canonical_group(h5, "inj-a", tic=101)
        _write_canonical_group(h5, "inj-b", tic=202)
        _write_canonical_group(h5, "not-selected", tic=999)
    table = repo / "selection.csv"
    pd.DataFrame(
        [
            {
                "native_input_include": True,
                "injection_id": "inj-a",
                "tic": 101,
                "sector": 56,
                "source_h5": "pairs/pair_a.h5",
                "h5_group": "/injections/inj-a",
            },
            {
                "native_input_include": True,
                "injection_id": "inj-b",
                "tic": 202,
                "sector": 56,
                "source_h5": "pairs/pair_b.h5",
                "h5_group": "/injections/inj-b",
            },
            {
                "native_input_include": False,
                "injection_id": "not-selected",
                "tic": 999,
                "sector": 56,
                "source_h5": "pairs/pair_a.h5",
                "h5_group": "/injections/not-selected",
            },
            {
                "native_input_include": True,
                "injection_id": "",
                "tic": 303,
                "sector": 56,
                "source_h5": "",
                "h5_group": "",
            },
        ]
    ).to_csv(table, index=False)
    return {
        "repo": repo,
        "pair_a": pair_a,
        "pair_b": pair_b,
        "canonical": canonical,
        "table": table,
        "out": tmp_path / "promoted.h5",
    }


def _promote(paths: dict[str, Path]) -> dict[str, object]:
    return promote_teacher_v3_injection_pairs(
        selection_table=paths["table"],
        source_pair_h5s=[paths["pair_a"], paths["pair_b"]],
        canonical_injection_h5=paths["canonical"],
        out_h5=paths["out"],
        repo_root=paths["repo"],
        expected_count=2,
        expected_source_counts=(1, 1),
    )


def test_promote_exact_selected_groups_and_rebind_canonical_path(
    tmp_path: Path,
) -> None:
    paths = _write_fixture(tmp_path)

    summary = _promote(paths)

    assert summary["n_injections"] == 2
    assert summary["out_h5_sha256"] == _sha256(paths["out"])
    with h5py.File(paths["out"], "r") as h5:
        assert h5.attrs["contract_version"] == (
            TEACHER_V3_INJECTION_PAIR_MERGE_CONTRACT
        )
        assert sorted(h5["injections"].keys()) == ["inj-a", "inj-b"]
        assert h5.attrs["source_injection_h5"] == str(
            paths["canonical"].resolve()
        )
        assert h5.attrs["source_injection_h5_sha256"] == _sha256(
            paths["canonical"]
        )
        source_hashes = json.loads(h5.attrs["source_pair_h5_sha256"])
        assert source_hashes[str(paths["pair_a"].resolve())] == _sha256(
            paths["pair_a"]
        )
        assert source_hashes[str(paths["pair_b"].resolve())] == _sha256(
            paths["pair_b"]
        )
        group = h5["injections/inj-a"]
        assert group.attrs["custom_attr"] == "preserve-me"
        assert group.attrs["source_injection_h5_original"] == (
            "/obsolete/raw_injections.h5"
        )
        assert group.attrs["source_injection_h5"] == str(
            paths["canonical"].resolve()
        )
        assert group.attrs["promoted_source_pair_h5_sha256"] == _sha256(
            paths["pair_a"]
        )
        assert "extra_dataset" in group
        assert group["time"].compression == "gzip"


@pytest.mark.parametrize("failure", ["missing", "duplicate"])
def test_promote_rejects_missing_or_duplicate_source_group(
    tmp_path: Path,
    failure: str,
) -> None:
    paths = _write_fixture(tmp_path)
    if failure == "missing":
        with h5py.File(paths["pair_a"], "a") as h5:
            del h5["injections/inj-a"]
    else:
        with h5py.File(paths["pair_b"], "a") as h5:
            _write_pair_group(h5, "inj-a", tic=101)

    with pytest.raises(ValueError, match=("missing" if failure == "missing" else "appears")):
        _promote(paths)
    assert not paths["out"].exists()


def test_promote_rejects_duplicate_selection_id(tmp_path: Path) -> None:
    paths = _write_fixture(tmp_path)
    table = pd.read_csv(paths["table"])
    table = pd.concat([table, table.iloc[[0]]], ignore_index=True)
    table.to_csv(paths["table"], index=False)

    with pytest.raises(ValueError, match="duplicate active injection IDs"):
        promote_teacher_v3_injection_pairs(
            selection_table=paths["table"],
            source_pair_h5s=[paths["pair_a"], paths["pair_b"]],
            canonical_injection_h5=paths["canonical"],
            out_h5=paths["out"],
            repo_root=paths["repo"],
            expected_count=3,
            expected_source_counts=(2, 1),
        )
    assert not paths["out"].exists()


@pytest.mark.parametrize(
    ("mutation", "match"),
    [
        ("pair_tic", "pair TIC"),
        ("canonical_tic", "canonical TIC"),
        ("pair_missing_dataset", "missing required datasets"),
        ("canonical_length", "dataset lengths disagree"),
    ],
)
def test_promote_rejects_mismatched_tics_and_datasets(
    tmp_path: Path,
    mutation: str,
    match: str,
) -> None:
    paths = _write_fixture(tmp_path)
    if mutation == "pair_tic":
        with h5py.File(paths["pair_a"], "a") as h5:
            h5["injections/inj-a"].attrs["tic"] = 999
    elif mutation == "canonical_tic":
        with h5py.File(paths["canonical"], "a") as h5:
            h5["injections/inj-a"].attrs["tic"] = 999
    elif mutation == "pair_missing_dataset":
        with h5py.File(paths["pair_a"], "a") as h5:
            del h5["injections/inj-a/DET_FLUX_ADP_original"]
    else:
        with h5py.File(paths["canonical"], "a") as h5:
            group = h5["injections/inj-a"]
            del group["RAW_FLUX_Primary_injected"]
            group.create_dataset("RAW_FLUX_Primary_injected", data=np.ones(7))

    with pytest.raises(ValueError, match=match):
        _promote(paths)
    assert not paths["out"].exists()
