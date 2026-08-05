from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

from twirl.vetting.teacher_split_registry import (
    TIC_SPLIT_POLICY_VERSION,
    attach_tic_split_registry,
    build_tic_split_registry,
    validate_tic_split_assignments,
    validate_tic_split_registry,
    write_tic_split_registry,
)


def _observation_corpus(n_tics: int = 60) -> pd.DataFrame:
    strata = (
        "planet_like",
        "eclipse_contact",
        "smooth_variable",
        "other",
        "broad_preserve",
    )
    rows: list[dict[str, object]] = []
    for index in range(n_tics):
        tic = 100_000 + index
        for sector in (56, 57):
            rows.append(
                {
                    "sector": sector,
                    "tic": tic,
                    "candidate_key": f"s{sector:04d}-tic{tic}",
                    "split_stratum": strata[index % len(strata)],
                }
            )
    return pd.DataFrame(rows)


def test_registry_is_deterministic_order_independent_and_tic_grouped() -> None:
    corpus = _observation_corpus()
    first, first_summary = build_tic_split_registry(corpus, seed=560062)
    shuffled, shuffled_summary = build_tic_split_registry(
        corpus.sample(frac=1.0, random_state=91).reset_index(drop=True),
        seed=560062,
    )

    pd.testing.assert_frame_equal(first, shuffled)
    assert (
        first_summary["observation_identity_content_sha256"]
        == shuffled_summary["observation_identity_content_sha256"]
    )
    assert first_summary["target_test_fraction"] == 0.20
    assert first_summary["test_observation_fraction"] == pytest.approx(0.20)
    assert set(first["fixed_split"]) == {"development", "test"}
    assert set(first.loc[first["fixed_split"].eq("development"), "cv_fold"]) == set(
        range(5)
    )
    assert first.loc[first["fixed_split"].eq("test"), "cv_fold"].eq(-1).all()
    assert first["tic"].is_unique
    assert first["split_policy_version"].eq(TIC_SPLIT_POLICY_VERSION).all()

    attached = attach_tic_split_registry(corpus, first)
    assert len(attached) == len(corpus)
    assert (
        attached.groupby("tic")[["fixed_split", "cv_fold"]].nunique().to_numpy()
        == 1
    ).all()

    with pytest.raises(ValueError, match="observation counts disagree"):
        attach_tic_split_registry(corpus.iloc[1:].reset_index(drop=True), first)

    changed_stratum = corpus.copy()
    changed_stratum.loc[0, "split_stratum"] = "other"
    with pytest.raises(ValueError, match="identities/strata bound"):
        attach_tic_split_registry(changed_stratum, first)

    swapped_strata = corpus.copy()
    first_tic = swapped_strata["tic"].iloc[0]
    tic_rows = swapped_strata.index[swapped_strata["tic"].eq(first_tic)]
    swapped_strata.loc[tic_rows, "split_stratum"] = [
        "planet_like",
        "other",
    ]
    registry_for_swapped, _ = build_tic_split_registry(
        swapped_strata, seed=560062
    )
    swapped_strata.loc[tic_rows, "split_stratum"] = [
        "other",
        "planet_like",
    ]
    with pytest.raises(ValueError, match="identities/strata bound"):
        attach_tic_split_registry(swapped_strata, registry_for_swapped)

    changed_identity = corpus.copy()
    changed_identity.loc[0, "candidate_key"] = "different-but-unique-key"
    with pytest.raises(ValueError, match="identities/strata bound"):
        attach_tic_split_registry(changed_identity, first)


def test_registry_locks_twenty_percent_of_tics_despite_group_imbalance() -> None:
    rows = [
        {
            "sector": 56 + index,
            "tic": 100_000,
            "candidate_key": f"large-group-{index}",
            "split_stratum": "planet_like",
        }
        for index in range(1_000)
    ]
    rows.extend(
        {
            "sector": 56,
            "tic": 100_000 + index,
            "candidate_key": f"small-group-{index}",
            "split_stratum": "planet_like",
        }
        for index in range(1, 21)
    )

    registry, summary = build_tic_split_registry(
        pd.DataFrame(rows), seed=0
    )

    assert len(registry) == 21
    assert registry["fixed_split"].eq("test").sum() == 4
    assert summary["target_test_tic_count"] == 4
    assert summary["test_fraction_unit"] == "unique_tics"
    assert summary["test_tic_fraction"] == pytest.approx(4 / 21)
    assert set(
        registry.loc[
            registry["fixed_split"].eq("development"), "cv_fold"
        ]
    ) == set(range(5))


def test_immutable_writer_is_byte_idempotent_and_hash_bound(tmp_path: Path) -> None:
    corpus_path = tmp_path / "corpus.csv"
    registry_path = tmp_path / "tic_split_registry.csv"
    summary_path = tmp_path / "tic_split_registry.summary.json"
    _observation_corpus().to_csv(corpus_path, index=False)

    first = write_tic_split_registry(
        corpus_path=corpus_path,
        registry_path=registry_path,
        summary_path=summary_path,
        seed=560062,
    )
    registry_bytes = registry_path.read_bytes()
    summary_bytes = summary_path.read_bytes()
    second = write_tic_split_registry(
        corpus_path=corpus_path,
        registry_path=registry_path,
        summary_path=summary_path,
        seed=560062,
    )

    assert first == second
    assert registry_path.read_bytes() == registry_bytes
    assert summary_path.read_bytes() == summary_bytes
    assert first["input_corpus"]["sha256"] == hashlib.sha256(
        corpus_path.read_bytes()
    ).hexdigest()
    assert first["output_registry"]["sha256"] == hashlib.sha256(
        registry_bytes
    ).hexdigest()
    assert json.loads(summary_bytes) == first

    with pytest.raises(FileExistsError, match="immutable output"):
        write_tic_split_registry(
            corpus_path=corpus_path,
            registry_path=registry_path,
            summary_path=summary_path,
            seed=123,
        )


@pytest.mark.parametrize(
    ("mutator", "message"),
    [
        (
            lambda frame: pd.concat([frame, frame.iloc[[0]]], ignore_index=True),
            "duplicate identities",
        ),
        (
            lambda frame: frame.assign(
                candidate_key=[
                    "reused-key",
                    "reused-key",
                    *frame["candidate_key"].iloc[2:].tolist(),
                ]
            ),
            "maps to multiple sector/TIC",
        ),
        (
            lambda frame: frame.assign(
                candidate_key=[
                    "",
                    *frame["candidate_key"].iloc[1:].tolist(),
                ]
            ),
            "blank identities",
        ),
        (
            lambda frame: frame.assign(
                split_stratum=[
                    "",
                    *frame["split_stratum"].iloc[1:].tolist(),
                ]
            ),
            "blank strata",
        ),
        (
            lambda frame: frame.assign(
                tic=[0, *frame["tic"].iloc[1:].tolist()]
            ),
            "non-positive",
        ),
    ],
)
def test_registry_rejects_bad_observation_identities(mutator, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        build_tic_split_registry(mutator(_observation_corpus()), seed=560062)


def test_split_validation_rejects_group_leakage_and_incomplete_columns() -> None:
    leaking = pd.DataFrame(
        {
            "tic": [1, 1, 2, 3, 4, 5, 6],
            "fixed_split": [
                "test",
                "development",
                "development",
                "development",
                "development",
                "development",
                "development",
            ],
            "cv_fold": [-1, 0, 0, 1, 2, 3, 4],
        }
    )
    with pytest.raises(ValueError, match="group leakage"):
        validate_tic_split_assignments(leaking)

    with pytest.raises(KeyError, match="cv_fold"):
        validate_tic_split_assignments(
            pd.DataFrame({"tic": [1], "fixed_split": ["test"]})
        )

    wrapped = pd.DataFrame(
        {
            "tic": [1, 2, 3, 4, 5, 6],
            "fixed_split": ["test"] + ["development"] * 5,
            "cv_fold": [65_535, 65_536, 1, 2, 3, 4],
        }
    )
    with pytest.raises(ValueError, match="test TICs must have cv_fold=-1"):
        validate_tic_split_assignments(wrapped)


def test_full_registry_validation_rejects_corrupt_provenance_metadata() -> None:
    registry, _ = build_tic_split_registry(_observation_corpus(), seed=560062)

    corrupt = registry.copy()
    corrupt.loc[0, "observation_record_sha256"] = "not-a-digest"
    with pytest.raises(ValueError, match="invalid digests"):
        validate_tic_split_registry(corrupt)

    corrupt = registry.copy()
    corrupt.loc[0, "split_policy_version"] = "unversioned"
    with pytest.raises(ValueError, match="split_policy_version"):
        validate_tic_split_registry(corrupt)


def test_registry_rejects_conflicting_preassigned_seed_contract() -> None:
    corpus = _observation_corpus()
    registry, _ = build_tic_split_registry(corpus, seed=560062)
    attached = attach_tic_split_registry(corpus, registry)

    with pytest.raises(ValueError, match="deterministic declared seed"):
        build_tic_split_registry(attached, seed=560063)
