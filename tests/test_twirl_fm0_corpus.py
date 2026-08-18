from __future__ import annotations

import csv
import json

from twirl.models.fm0.corpus import (
    expected_orbits,
    select_corpus_observations,
    write_corpus_selection,
)
from twirl.models.fm0.registry import build_alias_registry


def _row(source, tic, sector, orbit, camera=1, ccd=1, edge=False):
    return {
        "source_id": str(source),
        "tic_id": str(tic),
        "sector": sector,
        "orbit": orbit,
        "camera": camera,
        "ccd": ccd,
        "edge_warn": edge,
    }


def test_corpus_selection_keeps_one_complete_product_and_full_alias_closure(tmp_path):
    aliases = build_alias_registry(
        [
            {"gaia_dr3_source_id": "100", "tic_id": "10"},
            {"gaia_dr3_source_id": "100", "tic_id": "11"},
            {"gaia_dr3_source_id": "200", "tic_id": "20"},
        ]
    )
    s56 = expected_orbits(56)
    s57 = expected_orbits(57)
    observations = [
        *[_row(100, 10, 56, orbit, camera=2, ccd=3) for orbit in s56],
        _row(100, 11, 56, s56[0], camera=4, ccd=4),
        *[_row(200, 20, 57, orbit) for orbit in s57],
        _row(200, 20, 56, s56[0]),
        _row(200, 20, 56, s56[1], edge=True),
    ]
    selection, selected_aliases, audit = select_corpus_observations(
        observations, aliases, sectors=(56, 57)
    )
    assert [(row["gaia_dr3_source_id"], row["sector"]) for row in selection] == [
        ("100", 56),
        ("200", 57),
    ]
    assert json.loads(selection[0]["hdf5_paths_json"]) == [
        "orbit-119/ffi/cam2/ccd3/LC/10.h5",
        "orbit-120/ffi/cam2/ccd3/LC/10.h5",
    ]
    assert {(row["gaia_dr3_source_id"], row["tic_id"]) for row in selected_aliases} == {
        ("100", "10"),
        ("100", "11"),
        ("200", "20"),
    }
    assert audit[0]["reason"] == "no_complete_two_orbit_product"

    summary = write_corpus_selection(
        tmp_path,
        selection,
        selected_aliases,
        audit,
        sectors=(56, 57),
        input_authorities={"test": {"sha256": "a" * 64}},
    )
    assert summary["n_observations"] == 2
    assert summary["n_repeated_sources"] == 0
    assert summary["outputs_sha256"]["selection.csv"]
    with (tmp_path / "selection.csv").open(newline="") as handle:
        assert len(list(csv.DictReader(handle))) == 2


def test_corpus_selection_quarantines_multiple_complete_products():
    aliases = build_alias_registry(
        [
            {"gaia_dr3_source_id": "100", "tic_id": "10"},
            {"gaia_dr3_source_id": "100", "tic_id": "11"},
        ]
    )
    orbits = expected_orbits(56)
    rows = [
        *[_row(100, 10, 56, orbit, camera=1, ccd=1) for orbit in orbits],
        *[_row(100, 11, 56, orbit, camera=2, ccd=2) for orbit in orbits],
    ]
    try:
        select_corpus_observations(rows, aliases, sectors=(56,))
    except ValueError as exc:
        # With every visit quarantined, the builder must fail rather than
        # publish an apparently valid empty corpus.
        assert "no complete observations" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("ambiguous corpus unexpectedly passed")
