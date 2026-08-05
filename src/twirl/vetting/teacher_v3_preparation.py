"""Training-table and native-export preparation for Teacher v3."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import pandas as pd

from twirl.vetting.harmonic_inputs import (
    RAW_PAIR_CONTRACT_VERSION,
    native_group_path,
)
from twirl.vetting.teacher_split_registry import (
    attach_tic_split_registry,
    validate_tic_split_registry,
)
from twirl.vetting.teacher_v3 import (
    TEACHER_V3_CORPUS_VERSION,
    TEACHER_V3_RUN_NAME,
    file_sha256,
)


TEACHER_V3_NATIVE_PREPARATION_VERSION = (
    "teacher_v3_s56_s62_native_preparation_v1"
)


def _truth(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values
    return (
        values.fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
        .isin({"1", "1.0", "true", "t", "yes", "y"})
    )


def _publish(path: Path, payload: bytes) -> str:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.read_bytes() != payload:
        raise FileExistsError(
            f"refusing to replace immutable Teacher v3 preparation: {path}"
        )
    if not path.exists():
        path.write_bytes(payload)
    return hashlib.sha256(payload).hexdigest()


def _csv_payload(frame: pd.DataFrame) -> bytes:
    return frame.to_csv(
        index=False,
        lineterminator="\n",
        float_format="%.15g",
    ).encode("utf-8")


def _verify_bound_summaries(
    *,
    corpus_path: Path,
    corpus_summary_path: Path,
    split_registry_path: Path,
    split_summary_path: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    corpus_summary = json.loads(Path(corpus_summary_path).read_text())
    split_summary = json.loads(Path(split_summary_path).read_text())
    corpus_hash = file_sha256(corpus_path)
    split_hash = file_sha256(split_registry_path)
    if (
        corpus_summary.get("corpus_contract_version")
        != TEACHER_V3_CORPUS_VERSION
        or corpus_summary.get("teacher_run_name") != TEACHER_V3_RUN_NAME
    ):
        raise ValueError("corpus summary is not the frozen Teacher v3 contract")
    if (
        corpus_summary.get("output_corpus", {}).get("sha256")
        != corpus_hash
    ):
        raise ValueError("corpus summary SHA-256 does not match the corpus")
    if (
        split_summary.get("input_corpus", {}).get("sha256") != corpus_hash
        or split_summary.get("output_registry", {}).get("sha256") != split_hash
    ):
        raise ValueError(
            "split summary is not bound to the supplied corpus and registry"
        )
    return corpus_summary, split_summary


def prepare_teacher_v3_native_tables(
    *,
    corpus: pd.DataFrame,
    split_registry: pd.DataFrame,
    native_root: Path,
) -> tuple[
    pd.DataFrame,
    dict[int, pd.DataFrame],
    pd.DataFrame,
    pd.DataFrame,
    dict[str, Any],
]:
    """Attach immutable splits and assign each active row a native location."""

    required = {
        "sector",
        "tic",
        "source_uid",
        "is_injected_row",
        "injection_id",
        "teacher_v3_training_include",
        "teacher_v3_corpus_version",
        "teacher_run_name",
    }
    missing = sorted(required - set(corpus.columns))
    if missing:
        raise KeyError(f"Teacher v3 corpus is missing columns: {missing}")
    if not corpus["teacher_v3_corpus_version"].astype(str).eq(
        TEACHER_V3_CORPUS_VERSION
    ).all():
        raise ValueError("corpus rows have the wrong Teacher v3 contract")
    if not corpus["teacher_run_name"].astype(str).eq(
        TEACHER_V3_RUN_NAME
    ).all():
        raise ValueError("corpus rows have the wrong Teacher v3 run name")

    checked_registry = validate_tic_split_registry(split_registry)
    normalized_corpus = corpus.copy()
    if {
        "candidate_key",
        "observation_candidate_key",
    }.issubset(normalized_corpus.columns):
        candidate = (
            normalized_corpus["candidate_key"]
            .fillna("")
            .astype(str)
            .str.strip()
        )
        observation = (
            normalized_corpus["observation_candidate_key"]
            .fillna("")
            .astype(str)
            .str.strip()
        )
        normalized_corpus["candidate_key"] = candidate.where(
            candidate.ne(""),
            observation,
        )
    attached = attach_tic_split_registry(
        normalized_corpus,
        checked_registry,
    )
    active = _truth(attached["teacher_v3_training_include"])
    injected = _truth(attached["is_injected_row"])
    native_root = Path(native_root)
    sector_paths = {
        sector: str(
            native_root / f"s{sector:02d}_teacher_v3_native_v2.h5"
        )
        for sector in range(56, 63)
    }
    attached["native_h5_path"] = ""
    attached["native_group_path"] = ""
    attached["native_contract_version_expected"] = (
        RAW_PAIR_CONTRACT_VERSION
    )
    attached["native_input_include"] = active
    for sector, path in sector_paths.items():
        in_sector = attached["sector"].astype(int).eq(sector)
        attached.loc[in_sector & active, "native_h5_path"] = path

    for index, row in attached.loc[active].iterrows():
        attached.at[index, "native_group_path"] = native_group_path(row)
    missing_native_identity = active & (
        attached["native_h5_path"].eq("")
        | attached["native_group_path"].eq("")
    )
    if missing_native_identity.any():
        raise ValueError(
            "active Teacher v3 rows are missing native HDF5 identities"
        )
    active_injected = active & injected
    if attached.loc[active_injected, "injection_id"].fillna("").astype(
        str
    ).str.strip().eq("").any():
        raise ValueError("active injected rows are missing injection_id")
    if not attached.loc[
        active_injected, "native_group_path"
    ].str.startswith("injections/").all():
        raise ValueError("injected rows have non-injection native groups")
    active_real = active & ~injected
    if not attached.loc[
        active_real, "native_group_path"
    ].str.startswith("targets/").all():
        raise ValueError("real rows have non-target native groups")

    per_sector = {
        sector: attached.loc[
            attached["sector"].astype(int).eq(sector)
        ].copy()
        for sector in range(56, 63)
    }
    real_registry_source = (
        attached.loc[
            active_real,
            ["sector", "tic", "native_h5_path", "native_group_path"],
        ]
        .drop_duplicates()
        .sort_values(["sector", "tic"], kind="stable")
        .reset_index(drop=True)
    )
    if real_registry_source.duplicated(["sector", "tic"]).any():
        raise ValueError(
            "real native-registry source has conflicting observation mappings"
        )
    injection_mapping = (
        attached.loc[
            active_injected,
            [
                "source_uid",
                "review_id",
                "injection_id",
                "sector",
                "tic",
                "native_h5_path",
                "native_group_path",
            ],
        ]
        .sort_values(["sector", "source_uid"], kind="stable")
        .reset_index(drop=True)
    )
    summary: dict[str, Any] = {
        "preparation_contract_version": (
            TEACHER_V3_NATIVE_PREPARATION_VERSION
        ),
        "teacher_run_name": TEACHER_V3_RUN_NAME,
        "native_contract_version_expected": RAW_PAIR_CONTRACT_VERSION,
        "native_root": str(native_root),
        "sector_native_h5_paths": {
            str(key): value for key, value in sector_paths.items()
        },
        "n_rows": int(len(attached)),
        "n_active_rows": int(active.sum()),
        "n_active_real_rows": int(active_real.sum()),
        "n_active_injected_rows": int(active_injected.sum()),
        "n_real_registry_observations": int(len(real_registry_source)),
        "sector_row_counts": {
            str(sector): int(len(frame))
            for sector, frame in per_sector.items()
        },
        "sector_active_counts": {
            str(sector): int(
                _truth(frame["native_input_include"]).sum()
            )
            for sector, frame in per_sector.items()
        },
        "fixed_split_counts": {
            str(key): int(value)
            for key, value in attached["fixed_split"]
            .value_counts()
            .sort_index()
            .items()
        },
        "group_leakage_count": 0,
    }
    return (
        attached,
        per_sector,
        real_registry_source,
        injection_mapping,
        summary,
    )


def write_teacher_v3_native_preparation(
    *,
    corpus_path: Path,
    corpus_summary_path: Path,
    split_registry_path: Path,
    split_summary_path: Path,
    native_root: Path,
    out_dir: Path,
) -> dict[str, Any]:
    """Publish split-bound global/per-sector tables and native mappings."""

    corpus_path = Path(corpus_path)
    corpus_summary_path = Path(corpus_summary_path)
    split_registry_path = Path(split_registry_path)
    split_summary_path = Path(split_summary_path)
    corpus_summary, split_summary = _verify_bound_summaries(
        corpus_path=corpus_path,
        corpus_summary_path=corpus_summary_path,
        split_registry_path=split_registry_path,
        split_summary_path=split_summary_path,
    )
    input_paths = {
        "corpus": corpus_path,
        "corpus_summary": corpus_summary_path,
        "split_registry": split_registry_path,
        "split_summary": split_summary_path,
    }
    input_hashes = {
        name: file_sha256(path) for name, path in input_paths.items()
    }
    corpus = pd.read_csv(
        corpus_path,
        dtype=str,
        keep_default_na=False,
        low_memory=False,
    )
    split_registry = pd.read_csv(split_registry_path, low_memory=False)
    (
        attached,
        per_sector,
        real_registry_source,
        injection_mapping,
        summary,
    ) = prepare_teacher_v3_native_tables(
        corpus=corpus,
        split_registry=split_registry,
        native_root=native_root,
    )
    if {
        name: file_sha256(path) for name, path in input_paths.items()
    } != input_hashes:
        raise RuntimeError(
            "Teacher v3 corpus or split inputs changed during preparation"
        )

    out_dir = Path(out_dir)
    outputs: dict[str, dict[str, Any]] = {}

    def write_frame(name: str, path: Path, frame: pd.DataFrame) -> None:
        payload = _csv_payload(frame)
        outputs[name] = {
            "path": str(path.resolve()),
            "sha256": _publish(path, payload),
            "n_rows": int(len(frame)),
        }

    write_frame(
        "training_table_with_splits",
        out_dir / "training_table_with_splits.csv",
        attached,
    )
    for sector, frame in per_sector.items():
        write_frame(
            f"s{sector}_native_export_table",
            out_dir / f"s{sector}_native_export_table.csv",
            frame,
        )
    write_frame(
        "real_native_registry_source",
        out_dir / "real_native_registry_source.csv",
        real_registry_source,
    )
    write_frame(
        "injection_native_mapping",
        out_dir / "injection_native_mapping.csv",
        injection_mapping,
    )
    final_summary = {
        **summary,
        "inputs": {
            name: {
                "path": str(path.resolve()),
                "sha256": input_hashes[name],
            }
            for name, path in input_paths.items()
        },
        "corpus_summary_contract": corpus_summary[
            "corpus_contract_version"
        ],
        "split_policy_version": split_summary["split_policy_version"],
        "outputs": outputs,
    }
    summary_path = out_dir / "native_preparation.summary.json"
    summary_payload = (
        json.dumps(
            final_summary,
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")
    _publish(summary_path, summary_payload)
    return final_summary


__all__ = [
    "TEACHER_V3_NATIVE_PREPARATION_VERSION",
    "prepare_teacher_v3_native_tables",
    "write_teacher_v3_native_preparation",
]
