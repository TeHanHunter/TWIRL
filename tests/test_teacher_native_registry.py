from __future__ import annotations

from pathlib import Path

import h5py
import pandas as pd
import pytest

from twirl.vetting import harmonic_dataset
from twirl.vetting.adjudication_audit import HARMONIC_CNN_TARGET_POLICY
from twirl.vetting.harmonic_inputs import RAW_PAIR_CONTRACT_VERSION
from twirl.vetting.teacher_native_registry import (
    NATIVE_INPUT_REGISTRY_CONTRACT_VERSION,
    attach_native_input_registry,
    build_native_input_registry,
    file_sha256,
    validate_native_input_registry,
    validate_native_input_registry_path,
    write_native_input_registry,
    write_native_registry_attachment,
)


def _write_native(
    path: Path,
    groups: list[str],
    *,
    contract: str = RAW_PAIR_CONTRACT_VERSION,
) -> None:
    with h5py.File(path, "w") as h5:
        h5.attrs["contract_version"] = contract
        for group_path in groups:
            group = h5.create_group(group_path)
            parts = group_path.split("/")
            group.attrs["sector"] = (
                int(parts[1]) if parts[0] == "sectors" else 56
            )
            group.attrs["tic"] = int(parts[-1])


def _mapping_source(tmp_path: Path) -> tuple[pd.DataFrame, Path]:
    native_h5 = tmp_path / "native.h5"
    groups = [
        "sectors/0056/targets/0000000000000042",
        "sectors/0057/targets/0000000000000042",
    ]
    _write_native(native_h5, groups)
    source = pd.DataFrame(
        {
            "sector": [56, 56, 57],
            "tic": [42, 42, 42],
            "candidate_key": ["a", "b", "c"],
            "native_h5_path": [native_h5.name] * 3,
            "native_group_path": [groups[0], groups[0], groups[1]],
        }
    )
    return source, native_h5


def test_build_registry_keys_repeated_tics_by_sector_and_hashes_files(
    tmp_path: Path,
) -> None:
    source, native_h5 = _mapping_source(tmp_path)

    registry = build_native_input_registry(source, path_base=tmp_path)

    assert registry[["sector", "tic"]].to_dict(orient="records") == [
        {"sector": 56, "tic": 42},
        {"sector": 57, "tic": 42},
    ]
    assert registry["native_group_path"].nunique() == 2
    assert registry["native_h5_path"].eq(str(native_h5.resolve())).all()
    assert registry["native_h5_sha256"].eq(file_sha256(native_h5)).all()
    assert registry["native_contract_version"].eq(
        RAW_PAIR_CONTRACT_VERSION
    ).all()
    assert registry["native_registry_contract_version"].eq(
        NATIVE_INPUT_REGISTRY_CONTRACT_VERSION
    ).all()


def test_build_registry_rejects_blank_and_conflicting_source_mappings(
    tmp_path: Path,
) -> None:
    source, _ = _mapping_source(tmp_path)
    blank = source.copy()
    blank.loc[0, "native_group_path"] = " "
    with pytest.raises(ValueError, match="native_group_path.*blank"):
        build_native_input_registry(blank, path_base=tmp_path)

    conflict = source.copy()
    conflict.loc[1, "native_group_path"] = (
        "sectors/0057/targets/0000000000000042"
    )
    with pytest.raises(ValueError, match="conflicting native mappings"):
        build_native_input_registry(conflict, path_base=tmp_path)


def test_registry_rejects_duplicate_keys_and_tic_only_storage_collision(
    tmp_path: Path,
) -> None:
    source, _ = _mapping_source(tmp_path)
    registry = build_native_input_registry(source, path_base=tmp_path)

    duplicate = pd.concat([registry, registry.iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match=r"duplicate \(sector, tic\)"):
        validate_native_input_registry(duplicate)

    collision = registry.copy()
    for column in (
        "native_h5_path",
        "native_group_path",
        "native_h5_sha256",
        "native_contract_version",
    ):
        collision.loc[1, column] = collision.loc[0, column]
    with pytest.raises(ValueError, match="TIC-only/native-storage collision"):
        validate_native_input_registry(collision)


def test_registry_rejects_missing_files_groups_and_changed_hash(
    tmp_path: Path,
) -> None:
    source, native_h5 = _mapping_source(tmp_path)
    registry = build_native_input_registry(source, path_base=tmp_path)

    missing_file = registry.copy()
    missing_file.loc[0, "native_h5_path"] = str(tmp_path / "missing.h5")
    with pytest.raises(FileNotFoundError, match="missing native HDF5"):
        validate_native_input_registry(missing_file)

    missing_group = registry.copy()
    missing_group.loc[0, "native_group_path"] = "sectors/0056/targets/missing"
    with pytest.raises(KeyError, match="missing native HDF5 group"):
        validate_native_input_registry(missing_group)

    with h5py.File(native_h5, "a") as h5:
        h5.attrs["post_registry_change"] = True
    with pytest.raises(ValueError, match="SHA-256 changed"):
        validate_native_input_registry(registry)


def test_registry_rejects_group_identity_swaps(tmp_path: Path) -> None:
    source, _ = _mapping_source(tmp_path)
    source.loc[source["sector"].eq(56), "native_group_path"] = (
        "sectors/0057/targets/0000000000000042"
    )
    source.loc[source["sector"].eq(57), "native_group_path"] = (
        "sectors/0056/targets/0000000000000042"
    )

    with pytest.raises(ValueError, match=r"declares \(sector, tic\)"):
        build_native_input_registry(source, path_base=tmp_path)


def test_registry_rejects_wrong_or_blank_file_contract(tmp_path: Path) -> None:
    wrong = tmp_path / "wrong.h5"
    _write_native(wrong, ["targets/0000000000000042"], contract="wrong_v1")
    source = pd.DataFrame(
        {
            "sector": [56],
            "tic": [42],
            "native_h5_path": [wrong],
            "native_group_path": ["targets/0000000000000042"],
        }
    )
    with pytest.raises(ValueError, match="expected"):
        build_native_input_registry(source, path_base=tmp_path)

    blank = tmp_path / "blank.h5"
    with h5py.File(blank, "w") as h5:
        group = h5.create_group("targets/0000000000000042")
        group.attrs["sector"] = 56
        group.attrs["tic"] = 42
    source["native_h5_path"] = blank
    with pytest.raises(ValueError, match="blank contract_version"):
        build_native_input_registry(
            source,
            path_base=tmp_path,
            expected_contract_version=None,
        )


def test_attach_uses_sector_and_tic_and_rejects_tic_only_or_conflicts(
    tmp_path: Path,
) -> None:
    source, native_h5 = _mapping_source(tmp_path)
    registry = build_native_input_registry(source, path_base=tmp_path)
    corpus = pd.DataFrame(
        {
            "sector": [57, 56, 56],
            "tic": [42, 42, 42],
            "candidate_key": ["c", "a", "b"],
        }
    )

    attached = attach_native_input_registry(corpus, registry)

    expected_groups = {
        56: "sectors/0056/targets/0000000000000042",
        57: "sectors/0057/targets/0000000000000042",
    }
    assert attached["native_group_path"].tolist() == [
        expected_groups[int(sector)] for sector in corpus["sector"]
    ]
    assert attached["native_h5_path"].eq(str(native_h5.resolve())).all()
    assert attached["candidate_key"].tolist() == ["c", "a", "b"]

    with pytest.raises(KeyError, match="by TIC alone"):
        attach_native_input_registry(corpus.drop(columns="sector"), registry)

    conflicting = corpus.copy()
    conflicting["native_group_path"] = "targets/0000000000000042"
    with pytest.raises(ValueError, match="conflicts"):
        attach_native_input_registry(conflicting, registry)

    missing = corpus.iloc[[0]].copy()
    missing["sector"] = 58
    with pytest.raises(ValueError, match="no exact"):
        attach_native_input_registry(missing, registry)


def test_registry_and_attachment_summaries_bind_all_tables(tmp_path: Path) -> None:
    source, _ = _mapping_source(tmp_path)
    source_path = tmp_path / "source.csv"
    registry_path = tmp_path / "registry.csv"
    registry_summary_path = tmp_path / "registry.json"
    source.to_csv(source_path, index=False)

    summary = write_native_input_registry(
        source_path=source_path,
        registry_path=registry_path,
        summary_path=registry_summary_path,
    )
    audit = validate_native_input_registry_path(
        registry_path=registry_path,
        summary_path=registry_summary_path,
    )

    assert audit["passed"] is True
    assert summary["native_registry_sha256"] == file_sha256(registry_path)
    assert summary["source_table_sha256"] == file_sha256(source_path)
    assert summary["n_observations"] == 2
    assert summary["n_multisector_tics"] == 1
    registry_bytes = registry_path.read_bytes()
    summary_bytes = registry_summary_path.read_bytes()
    repeated = write_native_input_registry(
        source_path=source_path,
        registry_path=registry_path,
        summary_path=registry_summary_path,
    )
    assert repeated == summary
    assert registry_path.read_bytes() == registry_bytes
    assert registry_summary_path.read_bytes() == summary_bytes

    corpus_path = tmp_path / "corpus.parquet"
    output_path = tmp_path / "attached.parquet"
    attachment_summary_path = tmp_path / "attached.json"
    source[["sector", "tic", "candidate_key"]].to_parquet(
        corpus_path,
        index=False,
    )
    attachment = write_native_registry_attachment(
        corpus_path=corpus_path,
        registry_path=registry_path,
        registry_summary_path=registry_summary_path,
        output_path=output_path,
        summary_path=attachment_summary_path,
    )
    attached = pd.read_parquet(output_path)

    assert attachment["corpus_table_sha256"] == file_sha256(corpus_path)
    assert attachment["output_table_sha256"] == file_sha256(output_path)
    assert attached["native_registry_sha256"].eq(
        file_sha256(registry_path)
    ).all()
    assert len(attached) == len(source)
    output_bytes = output_path.read_bytes()
    attachment_summary_bytes = attachment_summary_path.read_bytes()
    repeated_attachment = write_native_registry_attachment(
        corpus_path=corpus_path,
        registry_path=registry_path,
        registry_summary_path=registry_summary_path,
        output_path=output_path,
        summary_path=attachment_summary_path,
    )
    assert repeated_attachment == attachment
    assert output_path.read_bytes() == output_bytes
    assert attachment_summary_path.read_bytes() == attachment_summary_bytes

    registry_path.write_text(registry_path.read_text() + "\n")
    with pytest.raises(ValueError, match="summary SHA-256"):
        validate_native_input_registry_path(
            registry_path=registry_path,
            summary_path=registry_summary_path,
        )


def test_prepare_training_rows_preserves_explicit_groups_and_s56_fallback(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rows = pd.DataFrame(
        {
            "review_id": ["explicit", "legacy"],
            "sector": [57, 56],
            "tic": [42, 43],
            "source_kind": ["real_candidate", "real_candidate"],
            "native_h5_path": ["/native/s57.h5", ""],
            "native_group_path": [
                "sectors/0057/targets/0000000000000042",
                "",
            ],
            "morphology_target_v1": ["other", "other"],
            "morphology_include_v1": [True, True],
            "preserve_target_v1": ["reject", "reject"],
            "preserve_include_v1": [True, True],
            "harmonic_target_v1": ["", ""],
            "harmonic_include_v1": [False, False],
            "broad_preserve_only": [False, False],
            "model_target_policy_version": [HARMONIC_CNN_TARGET_POLICY] * 2,
        }
    )

    def fake_splits(work: pd.DataFrame, *, seed: int) -> pd.DataFrame:
        assert seed == 56
        return pd.DataFrame(
            {
                "fixed_split": ["development"] * len(work),
                "cv_fold": [0] * len(work),
            }
        )

    monkeypatch.setattr(
        harmonic_dataset,
        "build_grouped_test_and_cv_folds",
        fake_splits,
    )

    prepared = harmonic_dataset.prepare_harmonic_training_rows(rows)

    assert prepared["native_group_path"].tolist() == [
        "sectors/0057/targets/0000000000000042",
        "targets/0000000000000043",
    ]
    assert prepared["native_h5_path"].tolist() == ["/native/s57.h5", ""]
