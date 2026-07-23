"""Observation-keyed native-input registries for Teacher-v1.

The legacy S56 teacher used one HDF5 file and addressed real targets by TIC.
That is ambiguous once the same TIC contributes observations from more than
one sector.  This module freezes the complete native storage identity as
``(sector, tic) -> (native_h5_path, native_group_path)`` and binds every HDF5
file by its SHA-256 digest and root contract version.
"""
from __future__ import annotations

import hashlib
import io
import json
import os
from pathlib import Path
import tempfile
from typing import Any, Mapping

import numpy as np
import pandas as pd

from twirl.vetting.harmonic_inputs import RAW_PAIR_CONTRACT_VERSION


NATIVE_INPUT_REGISTRY_CONTRACT_VERSION = (
    "twirl_teacher_v1_observation_native_registry_v1"
)
NATIVE_INPUT_REGISTRY_SUMMARY_VERSION = (
    "twirl_teacher_v1_observation_native_registry_summary_v1"
)
NATIVE_INPUT_ATTACHMENT_SUMMARY_VERSION = (
    "twirl_teacher_v1_native_registry_attachment_summary_v1"
)
REGISTRY_KEY_COLUMNS: tuple[str, str] = ("sector", "tic")
REGISTRY_STORAGE_COLUMNS: tuple[str, str] = (
    "native_h5_path",
    "native_group_path",
)
REGISTRY_FILE_METADATA_COLUMNS: tuple[str, str] = (
    "native_h5_sha256",
    "native_contract_version",
)
REGISTRY_COLUMNS: tuple[str, ...] = (
    "native_registry_contract_version",
    *REGISTRY_KEY_COLUMNS,
    *REGISTRY_STORAGE_COLUMNS,
    *REGISTRY_FILE_METADATA_COLUMNS,
)
MAX_EXACT_INTEGER = 2**53 - 1


def file_sha256(path: Path) -> str:
    """Return the SHA-256 digest of one file."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_table(path: Path) -> pd.DataFrame:
    """Read one CSV or Parquet table."""

    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path, low_memory=False)
    if suffix in {".parquet", ".pq"}:
        return pd.read_parquet(path)
    raise ValueError(f"unsupported table format: {path}")


def _table_payload(table: pd.DataFrame, path: Path) -> bytes:
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return table.to_csv(index=False, lineterminator="\n").encode("utf-8")
    if suffix in {".parquet", ".pq"}:
        buffer = io.BytesIO()
        table.to_parquet(buffer, index=False)
        return buffer.getvalue()
    raise ValueError(f"unsupported table format: {path}")


def _publish_immutable(path: Path, payload: bytes) -> None:
    """Publish bytes once; identical retries succeed and differences fail."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != payload:
            raise FileExistsError(
                "refusing to replace immutable output with different bytes: "
                f"{path}"
            )
        return
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.link(temporary, path)
        except FileExistsError:
            if path.read_bytes() != payload:
                raise
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _write_table(table: pd.DataFrame, path: Path) -> None:
    """Byte-idempotently publish one immutable CSV or Parquet table."""

    _publish_immutable(path, _table_payload(table, path))


def _write_json(payload: Mapping[str, Any], path: Path) -> None:
    """Byte-idempotently publish one immutable JSON object."""

    serialized = (
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    _publish_immutable(path, serialized)


def _positive_integer(values: pd.Series, *, name: str) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    array = numeric.to_numpy(dtype=float)
    invalid = (
        ~np.isfinite(array)
        | (array <= 0)
        | (array > MAX_EXACT_INTEGER)
        | ~np.isclose(array, np.rint(array), rtol=0.0, atol=0.0)
    )
    if invalid.any():
        examples = values.loc[invalid].head(5).tolist()
        raise ValueError(
            f"{name} contains {int(invalid.sum())} blank, non-integral, "
            f"non-positive, or out-of-range values; first={examples}"
        )
    return pd.Series(
        np.rint(array).astype(np.int64),
        index=values.index,
        name=name,
    )


def _group_positive_integer(group: Any, *, name: str, context: str) -> int:
    if name not in group.attrs:
        raise ValueError(f"{context} is missing required {name!r} attribute")
    raw = np.asarray(group.attrs[name])
    if raw.size != 1:
        raise ValueError(f"{context} has a non-scalar {name!r} attribute")
    try:
        value = float(raw.reshape(-1)[0])
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{context} has an invalid {name!r} attribute") from exc
    if (
        not np.isfinite(value)
        or value <= 0
        or value > MAX_EXACT_INTEGER
        or not np.isclose(value, np.rint(value), rtol=0.0, atol=0.0)
    ):
        raise ValueError(f"{context} has an invalid {name!r} attribute")
    return int(np.rint(value))


def _nonblank_text(values: pd.Series, *, name: str) -> pd.Series:
    normalized = values.fillna("").astype(str).str.strip()
    blank = normalized.eq("")
    if blank.any():
        raise ValueError(f"{name} contains {int(blank.sum())} blank values")
    return normalized


def _absolute_native_paths(values: pd.Series, *, path_base: Path) -> pd.Series:
    normalized = _nonblank_text(values, name="native_h5_path")
    base = Path(path_base).expanduser().resolve()

    def resolve_one(value: str) -> str:
        candidate = Path(value).expanduser()
        if not candidate.is_absolute():
            candidate = base / candidate
        return str(candidate.resolve())

    return normalized.map(resolve_one)


def _normalize_source_mappings(
    rows: pd.DataFrame,
    *,
    path_base: Path,
) -> pd.DataFrame:
    required = {*REGISTRY_KEY_COLUMNS, *REGISTRY_STORAGE_COLUMNS}
    missing = sorted(required - set(rows.columns))
    if missing:
        raise KeyError(f"native mapping table is missing columns: {missing}")
    if rows.empty:
        raise ValueError("native mapping table is empty")
    work = rows.copy()
    work["sector"] = _positive_integer(work["sector"], name="sector")
    work["tic"] = _positive_integer(work["tic"], name="tic")
    work["native_h5_path"] = _absolute_native_paths(
        work["native_h5_path"],
        path_base=path_base,
    )
    work["native_group_path"] = _nonblank_text(
        work["native_group_path"],
        name="native_group_path",
    )
    return work


def _collapse_source_observations(rows: pd.DataFrame) -> pd.DataFrame:
    """Collapse repeated candidate rows only when their native mapping agrees."""

    records: list[pd.Series] = []
    for key, group in rows.groupby(list(REGISTRY_KEY_COLUMNS), sort=True):
        mappings = group.loc[:, list(REGISTRY_STORAGE_COLUMNS)].drop_duplicates()
        if len(mappings) != 1:
            examples = mappings.head(5).to_dict(orient="records")
            raise ValueError(
                "conflicting native mappings for observation "
                f"(sector={int(key[0])}, tic={int(key[1])}); first={examples}"
            )
        records.append(group.iloc[0])
    return pd.DataFrame(records).reset_index(drop=True)


def _validate_registry_shape(
    registry: pd.DataFrame,
    *,
    path_base: Path,
) -> pd.DataFrame:
    missing = sorted(set(REGISTRY_COLUMNS) - set(registry.columns))
    if missing:
        raise KeyError(f"native registry is missing columns: {missing}")
    if registry.empty:
        raise ValueError("native registry is empty")
    work = registry.loc[:, list(REGISTRY_COLUMNS)].copy()
    work["sector"] = _positive_integer(work["sector"], name="sector")
    work["tic"] = _positive_integer(work["tic"], name="tic")
    work["native_h5_path"] = _absolute_native_paths(
        work["native_h5_path"],
        path_base=path_base,
    )
    work["native_group_path"] = _nonblank_text(
        work["native_group_path"],
        name="native_group_path",
    )
    work["native_h5_sha256"] = _nonblank_text(
        work["native_h5_sha256"],
        name="native_h5_sha256",
    ).str.lower()
    valid_digest = work["native_h5_sha256"].str.fullmatch(r"[0-9a-f]{64}")
    if not valid_digest.all():
        raise ValueError(
            "native_h5_sha256 contains "
            f"{int((~valid_digest).sum())} invalid SHA-256 digests"
        )
    work["native_contract_version"] = _nonblank_text(
        work["native_contract_version"],
        name="native_contract_version",
    )
    work["native_registry_contract_version"] = _nonblank_text(
        work["native_registry_contract_version"],
        name="native_registry_contract_version",
    )
    wrong_registry_contract = work["native_registry_contract_version"].ne(
        NATIVE_INPUT_REGISTRY_CONTRACT_VERSION
    )
    if wrong_registry_contract.any():
        observed = sorted(
            work.loc[
                wrong_registry_contract, "native_registry_contract_version"
            ].unique()
        )
        raise ValueError(
            "native registry has unexpected contract versions: "
            f"{observed}"
        )

    duplicate_key = work.duplicated(list(REGISTRY_KEY_COLUMNS), keep=False)
    if duplicate_key.any():
        examples = (
            work.loc[duplicate_key, list(REGISTRY_KEY_COLUMNS)]
            .head(5)
            .to_dict(orient="records")
        )
        raise ValueError(
            "native registry contains duplicate (sector, tic) observations; "
            f"first={examples}"
        )

    duplicate_storage = work.duplicated(
        list(REGISTRY_STORAGE_COLUMNS),
        keep=False,
    )
    if duplicate_storage.any():
        examples = (
            work.loc[
                duplicate_storage,
                [*REGISTRY_KEY_COLUMNS, *REGISTRY_STORAGE_COLUMNS],
            ]
            .head(5)
            .to_dict(orient="records")
        )
        raise ValueError(
            "native registry has a TIC-only/native-storage collision: one "
            "(native_h5_path, native_group_path) maps to multiple "
            f"(sector, tic) observations; first={examples}"
        )

    for path, group in work.groupby("native_h5_path", sort=False):
        if group["native_h5_sha256"].nunique() != 1:
            raise ValueError(
                f"native file {path} has conflicting SHA-256 declarations"
            )
        if group["native_contract_version"].nunique() != 1:
            raise ValueError(
                f"native file {path} has conflicting contract declarations"
            )
    return work.sort_values(list(REGISTRY_KEY_COLUMNS)).reset_index(drop=True)


def _inspect_native_files(
    registry: pd.DataFrame,
    *,
    expected_contract_version: str | None,
    compare_declared_metadata: bool,
) -> dict[str, dict[str, Any]]:
    import h5py

    inspections: dict[str, dict[str, Any]] = {}
    for path_text, rows in registry.groupby("native_h5_path", sort=True):
        path = Path(path_text)
        if not path.is_file():
            raise FileNotFoundError(f"missing native HDF5 file: {path}")
        digest = file_sha256(path)
        with h5py.File(path, "r") as h5:
            contract = str(h5.attrs.get("contract_version", "")).strip()
            if not contract:
                raise ValueError(
                    f"native HDF5 file has blank contract_version: {path}"
                )
            if (
                expected_contract_version is not None
                and contract != expected_contract_version
            ):
                raise ValueError(
                    f"native HDF5 file {path} has contract_version={contract!r}; "
                    f"expected {expected_contract_version!r}"
                )
            for registry_row in rows.to_dict("records"):
                group_path = str(registry_row["native_group_path"])
                if group_path not in h5:
                    raise KeyError(
                        f"missing native HDF5 group {group_path!r} in {path}"
                    )
                if not isinstance(h5[group_path], h5py.Group):
                    raise ValueError(
                        f"native HDF5 path {group_path!r} in {path} is not a group"
                    )
                group = h5[group_path]
                context = f"native HDF5 group {group_path!r} in {path}"
                observed_sector = _group_positive_integer(
                    group,
                    name="sector",
                    context=context,
                )
                observed_tic = _group_positive_integer(
                    group,
                    name="tic",
                    context=context,
                )
                expected_sector = int(registry_row["sector"])
                expected_tic = int(registry_row["tic"])
                if (observed_sector, observed_tic) != (
                    expected_sector,
                    expected_tic,
                ):
                    raise ValueError(
                        f"{context} declares (sector, tic)="
                        f"({observed_sector}, {observed_tic}), expected "
                        f"({expected_sector}, {expected_tic})"
                    )
        if compare_declared_metadata:
            declared_digest = str(rows["native_h5_sha256"].iloc[0]).lower()
            declared_contract = str(rows["native_contract_version"].iloc[0])
            if digest != declared_digest:
                raise ValueError(
                    f"native HDF5 SHA-256 changed for {path}: "
                    f"{digest} != {declared_digest}"
                )
            if contract != declared_contract:
                raise ValueError(
                    f"native HDF5 contract changed for {path}: "
                    f"{contract!r} != {declared_contract!r}"
                )
        inspections[str(path)] = {
            "native_h5_path": str(path),
            "native_h5_sha256": digest,
            "native_contract_version": contract,
            "native_h5_size_bytes": int(path.stat().st_size),
            "n_observations": int(len(rows)),
        }
    return inspections


def build_native_input_registry(
    source_rows: pd.DataFrame,
    *,
    path_base: Path = Path("."),
    expected_contract_version: str | None = RAW_PAIR_CONTRACT_VERSION,
) -> pd.DataFrame:
    """Build and fully validate one observation-keyed native registry.

    A source corpus may contain several candidate rows for the same
    ``(sector, tic)``.  Such rows collapse to one registry entry only when
    their explicit HDF5 path and group agree exactly.
    """

    source = _normalize_source_mappings(source_rows, path_base=path_base)
    observations = _collapse_source_observations(source)
    inspections = _inspect_native_files(
        observations,
        expected_contract_version=expected_contract_version,
        compare_declared_metadata=False,
    )
    observations["native_registry_contract_version"] = (
        NATIVE_INPUT_REGISTRY_CONTRACT_VERSION
    )
    observations["native_h5_sha256"] = observations["native_h5_path"].map(
        lambda value: inspections[str(value)]["native_h5_sha256"]
    )
    observations["native_contract_version"] = observations[
        "native_h5_path"
    ].map(lambda value: inspections[str(value)]["native_contract_version"])
    registry = observations.loc[:, list(REGISTRY_COLUMNS)]
    return validate_native_input_registry(
        registry,
        path_base=path_base,
        expected_contract_version=expected_contract_version,
    )


def validate_native_input_registry(
    registry: pd.DataFrame,
    *,
    path_base: Path = Path("."),
    expected_contract_version: str | None = RAW_PAIR_CONTRACT_VERSION,
) -> pd.DataFrame:
    """Validate identities, storage uniqueness, files, hashes, and groups."""

    normalized = _validate_registry_shape(registry, path_base=path_base)
    _inspect_native_files(
        normalized,
        expected_contract_version=expected_contract_version,
        compare_declared_metadata=True,
    )
    return normalized


def _native_file_records(registry: pd.DataFrame) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for path, rows in registry.groupby("native_h5_path", sort=True):
        records.append(
            {
                "native_h5_path": str(path),
                "native_h5_sha256": str(rows["native_h5_sha256"].iloc[0]),
                "native_contract_version": str(
                    rows["native_contract_version"].iloc[0]
                ),
                "native_h5_size_bytes": int(Path(path).stat().st_size),
                "n_observations": int(len(rows)),
            }
        )
    return records


def _registry_counts(registry: pd.DataFrame) -> dict[str, int]:
    tic_sector_counts = registry.groupby("tic", sort=False)["sector"].nunique()
    return {
        "n_observations": int(len(registry)),
        "n_unique_tics": int(registry["tic"].nunique()),
        "n_multisector_tics": int(tic_sector_counts.gt(1).sum()),
        "n_sectors": int(registry["sector"].nunique()),
        "n_native_h5_files": int(registry["native_h5_path"].nunique()),
    }


def write_native_input_registry(
    *,
    source_path: Path,
    registry_path: Path,
    summary_path: Path,
    path_base: Path | None = None,
    expected_contract_version: str | None = RAW_PAIR_CONTRACT_VERSION,
) -> dict[str, Any]:
    """Build a registry and write a hash-bound JSON summary."""

    source_path = Path(source_path).resolve()
    registry_path = Path(registry_path)
    summary_path = Path(summary_path)
    source = read_table(source_path)
    registry = build_native_input_registry(
        source,
        path_base=path_base or source_path.parent,
        expected_contract_version=expected_contract_version,
    )
    _write_table(registry, registry_path)
    summary: dict[str, Any] = {
        "summary_contract_version": NATIVE_INPUT_REGISTRY_SUMMARY_VERSION,
        "native_registry_contract_version": (
            NATIVE_INPUT_REGISTRY_CONTRACT_VERSION
        ),
        "key_columns": list(REGISTRY_KEY_COLUMNS),
        "storage_columns": list(REGISTRY_STORAGE_COLUMNS),
        "source_table": str(source_path),
        "source_table_sha256": file_sha256(source_path),
        "n_source_rows": int(len(source)),
        "native_registry": str(registry_path.resolve()),
        "native_registry_sha256": file_sha256(registry_path),
        "expected_native_contract_version": (
            expected_contract_version if expected_contract_version is not None else ""
        ),
        **_registry_counts(registry),
        "native_files": _native_file_records(registry),
    }
    _write_json(summary, summary_path)
    return summary


def validate_native_input_registry_path(
    *,
    registry_path: Path,
    summary_path: Path | None = None,
    expected_contract_version: str | None = RAW_PAIR_CONTRACT_VERSION,
) -> dict[str, Any]:
    """Validate a registry file and, when supplied, its frozen summary."""

    registry_path = Path(registry_path).resolve()
    registry = validate_native_input_registry(
        read_table(registry_path),
        path_base=registry_path.parent,
        expected_contract_version=expected_contract_version,
    )
    digest = file_sha256(registry_path)
    audit: dict[str, Any] = {
        "passed": True,
        "native_registry": str(registry_path),
        "native_registry_sha256": digest,
        "native_registry_contract_version": (
            NATIVE_INPUT_REGISTRY_CONTRACT_VERSION
        ),
        "expected_native_contract_version": (
            expected_contract_version if expected_contract_version is not None else ""
        ),
        **_registry_counts(registry),
        "native_files": _native_file_records(registry),
    }
    if summary_path is not None:
        summary_path = Path(summary_path).resolve()
        try:
            summary = json.loads(summary_path.read_text())
        except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise ValueError(
                f"invalid native registry summary {summary_path}: {exc}"
            ) from exc
        if not isinstance(summary, dict):
            raise ValueError("native registry summary must be a JSON object")
        if (
            summary.get("summary_contract_version")
            != NATIVE_INPUT_REGISTRY_SUMMARY_VERSION
        ):
            raise ValueError("native registry summary has the wrong contract")
        if summary.get("native_registry_sha256") != digest:
            raise ValueError(
                "native registry summary SHA-256 does not match the registry"
            )
        if (
            summary.get("native_registry_contract_version")
            != NATIVE_INPUT_REGISTRY_CONTRACT_VERSION
        ):
            raise ValueError(
                "native registry summary has the wrong registry contract"
            )
        for name, value in _registry_counts(registry).items():
            if int(summary.get(name, -1)) != value:
                raise ValueError(
                    f"native registry summary {name} does not match the registry"
                )
        if summary.get("native_files") != _native_file_records(registry):
            raise ValueError(
                "native registry summary file metadata does not match the registry"
            )
        audit["summary"] = str(summary_path)
        audit["summary_sha256"] = file_sha256(summary_path)
    return audit


def attach_native_input_registry(
    corpus: pd.DataFrame,
    registry: pd.DataFrame,
    *,
    corpus_path_base: Path = Path("."),
    registry_path_base: Path = Path("."),
    expected_contract_version: str | None = RAW_PAIR_CONTRACT_VERSION,
) -> pd.DataFrame:
    """Attach native identities to a corpus using both sector and TIC."""

    required = set(REGISTRY_KEY_COLUMNS)
    missing = sorted(required - set(corpus.columns))
    if missing:
        raise KeyError(
            "cannot attach a native registry by TIC alone; corpus is missing "
            f"observation key columns: {missing}"
        )
    if corpus.empty:
        raise ValueError("native attachment corpus is empty")
    normalized_registry = validate_native_input_registry(
        registry,
        path_base=registry_path_base,
        expected_contract_version=expected_contract_version,
    )
    work = corpus.copy()
    work["sector"] = _positive_integer(work["sector"], name="sector")
    work["tic"] = _positive_integer(work["tic"], name="tic")
    work["_native_registry_row_order"] = np.arange(len(work), dtype=np.int64)

    attachment_columns = [
        "native_registry_contract_version",
        *REGISTRY_STORAGE_COLUMNS,
        *REGISTRY_FILE_METADATA_COLUMNS,
    ]
    mapping = normalized_registry.loc[
        :, [*REGISTRY_KEY_COLUMNS, *attachment_columns]
    ].rename(columns={name: f"_registry_{name}" for name in attachment_columns})
    joined = work.merge(
        mapping,
        on=list(REGISTRY_KEY_COLUMNS),
        how="left",
        validate="many_to_one",
        sort=False,
        indicator="_native_registry_merge",
    )
    missing_mapping = joined["_native_registry_merge"].ne("both")
    if missing_mapping.any():
        examples = (
            joined.loc[missing_mapping, list(REGISTRY_KEY_COLUMNS)]
            .drop_duplicates()
            .head(5)
            .to_dict(orient="records")
        )
        raise ValueError(
            "native registry has no exact (sector, tic) mapping for "
            f"{int(missing_mapping.sum())} corpus rows; first={examples}"
        )

    for column in attachment_columns:
        registry_column = f"_registry_{column}"
        if column in joined:
            existing = joined[column]
            nonblank = existing.notna() & existing.astype(str).str.strip().ne("")
            if nonblank.any():
                if column == "native_h5_path":
                    existing_normalized = _absolute_native_paths(
                        existing.loc[nonblank],
                        path_base=corpus_path_base,
                    )
                else:
                    existing_normalized = (
                        existing.loc[nonblank].astype(str).str.strip()
                    )
                    if column == "native_h5_sha256":
                        existing_normalized = existing_normalized.str.lower()
                expected = joined.loc[nonblank, registry_column].astype(str)
                mismatch = existing_normalized.ne(expected)
                if mismatch.any():
                    rows = mismatch[mismatch].index[:5]
                    examples = joined.loc[
                        rows,
                        [*REGISTRY_KEY_COLUMNS, column, registry_column],
                    ].to_dict(orient="records")
                    raise ValueError(
                        f"corpus {column} conflicts with the observation-keyed "
                        f"native registry; first={examples}"
                    )
            joined[column] = joined[registry_column]
        else:
            joined[column] = joined[registry_column]
        joined = joined.drop(columns=registry_column)
    joined = (
        joined.sort_values("_native_registry_row_order")
        .drop(columns=["_native_registry_row_order", "_native_registry_merge"])
        .reset_index(drop=True)
    )
    return joined


def write_native_registry_attachment(
    *,
    corpus_path: Path,
    registry_path: Path,
    output_path: Path,
    summary_path: Path,
    registry_summary_path: Path | None = None,
    expected_contract_version: str | None = RAW_PAIR_CONTRACT_VERSION,
) -> dict[str, Any]:
    """Attach a validated registry and write a hash-bound output summary."""

    corpus_path = Path(corpus_path).resolve()
    registry_path = Path(registry_path).resolve()
    output_path = Path(output_path)
    registry_audit = validate_native_input_registry_path(
        registry_path=registry_path,
        summary_path=registry_summary_path,
        expected_contract_version=expected_contract_version,
    )
    registry = read_table(registry_path)
    corpus = read_table(corpus_path)
    attached = attach_native_input_registry(
        corpus,
        registry,
        corpus_path_base=corpus_path.parent,
        registry_path_base=registry_path.parent,
        expected_contract_version=expected_contract_version,
    )
    attached["native_registry_sha256"] = registry_audit[
        "native_registry_sha256"
    ]
    _write_table(attached, output_path)
    summary: dict[str, Any] = {
        "summary_contract_version": NATIVE_INPUT_ATTACHMENT_SUMMARY_VERSION,
        "native_registry_contract_version": (
            NATIVE_INPUT_REGISTRY_CONTRACT_VERSION
        ),
        "corpus_table": str(corpus_path),
        "corpus_table_sha256": file_sha256(corpus_path),
        "native_registry": str(registry_path),
        "native_registry_sha256": registry_audit["native_registry_sha256"],
        "native_registry_summary": (
            str(Path(registry_summary_path).resolve())
            if registry_summary_path is not None
            else ""
        ),
        "output_table": str(output_path.resolve()),
        "output_table_sha256": file_sha256(output_path),
        "n_input_rows": int(len(corpus)),
        "n_output_rows": int(len(attached)),
        "n_observations": int(
            attached.loc[:, list(REGISTRY_KEY_COLUMNS)].drop_duplicates().shape[0]
        ),
        "n_unique_tics": int(attached["tic"].nunique()),
        "n_native_h5_files": int(attached["native_h5_path"].nunique()),
        "native_files": registry_audit["native_files"],
    }
    _write_json(summary, summary_path)
    return summary


__all__ = [
    "NATIVE_INPUT_ATTACHMENT_SUMMARY_VERSION",
    "NATIVE_INPUT_REGISTRY_CONTRACT_VERSION",
    "NATIVE_INPUT_REGISTRY_SUMMARY_VERSION",
    "REGISTRY_COLUMNS",
    "REGISTRY_FILE_METADATA_COLUMNS",
    "REGISTRY_KEY_COLUMNS",
    "REGISTRY_STORAGE_COLUMNS",
    "attach_native_input_registry",
    "build_native_input_registry",
    "file_sha256",
    "read_table",
    "validate_native_input_registry",
    "validate_native_input_registry_path",
    "write_native_input_registry",
    "write_native_registry_attachment",
]
