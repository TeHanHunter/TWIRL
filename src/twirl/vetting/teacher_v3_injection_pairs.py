"""Promote the frozen Teacher-v3 injection pairs into one exact HDF5.

The frozen S56 table points at two ADP-pair products.  The native-input
builder accepts one explicit pair file, so this module copies only the active
injection groups named by the table into one provenance-bound product.  Raw
injected flux remains in a separate canonical injection HDF5; the promoted
root and every copied group are rebound to that one explicit path.
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd


TEACHER_V3_INJECTION_PAIR_MERGE_CONTRACT = (
    "teacher_v3_s56_injection_pair_merge_v1"
)
TEACHER_V3_S56_EXPECTED_INJECTIONS = 456
TEACHER_V3_S56_EXPECTED_SOURCE_COUNTS: tuple[int, int] = (300, 156)

PAIR_REQUIRED_ATTRS: tuple[str, ...] = ("tic", "sector", "camera", "ccd")
PAIR_REQUIRED_DATASETS: tuple[str, ...] = (
    "time",
    "cadenceno",
    "orbitid",
    "quality",
    "transit_model",
    "DET_FLUX_ADP_SML_injected",
    "DET_FLUX_ADP_injected",
    "DET_FLUX_ADP_SML_original",
    "DET_FLUX_ADP_original",
)
CANONICAL_REQUIRED_DATASETS: tuple[str, ...] = (
    "RAW_FLUX_Small_injected",
    "RAW_FLUX_Small_original",
    "RAW_FLUX_Primary_injected",
    "RAW_FLUX_Primary_original",
)
_TRUE_VALUES = {"1", "true", "t", "yes", "y"}
_FALSE_VALUES = {"", "0", "false", "f", "no", "n"}


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _read_table(path: Path) -> pd.DataFrame:
    suffix = Path(path).suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path, low_memory=False)
    if suffix == ".parquet":
        return pd.read_parquet(path)
    raise ValueError(f"unsupported selection-table format: {path}")


def _strict_bool(value: Any, *, context: str) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return False
    text = str(value).strip().lower()
    if text in _TRUE_VALUES:
        return True
    if text in _FALSE_VALUES:
        return False
    raise ValueError(f"{context}: invalid boolean value {value!r}")


def _strict_int(value: Any, *, context: str) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{context}: boolean is not a valid integer")
    try:
        numeric = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{context}: invalid integer {value!r}") from exc
    if not np.isfinite(numeric) or not numeric.is_integer():
        raise ValueError(f"{context}: invalid integer {value!r}")
    return int(numeric)


def _resolve_path(value: Any, *, repo_root: Path) -> Path:
    path = Path(str(value))
    if not path.is_absolute():
        path = Path(repo_root) / path
    return path.resolve()


def _normalized_selection_digest(records: Sequence[dict[str, Any]]) -> str:
    payload = json.dumps(
        list(records),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _validate_one_dimensional_datasets(
    group: Any,
    *,
    required: Sequence[str],
    context: str,
    expected_length: int | None = None,
) -> int:
    import h5py

    missing = [name for name in required if name not in group]
    if missing:
        raise ValueError(f"{context}: missing required datasets {missing}")
    lengths: dict[str, int] = {}
    for name in required:
        dataset = group[name]
        if not isinstance(dataset, h5py.Dataset):
            raise ValueError(f"{context}: {name} is not a dataset")
        if dataset.ndim != 1:
            raise ValueError(
                f"{context}: dataset {name} must be one-dimensional; "
                f"got shape {dataset.shape}"
            )
        lengths[name] = int(dataset.shape[0])
    if len(set(lengths.values())) != 1:
        raise ValueError(f"{context}: required dataset lengths disagree: {lengths}")
    length = next(iter(lengths.values()))
    if length < 1:
        raise ValueError(f"{context}: required datasets are empty")
    if expected_length is not None and length != int(expected_length):
        raise ValueError(
            f"{context}: dataset length {length} does not match pair length "
            f"{int(expected_length)}"
        )
    return length


def _validate_pair_group(
    *,
    pair: Any,
    canonical: Any,
    injection_id: str,
    table_tic: int,
    expected_sector: int,
) -> None:
    context = f"injection {injection_id}"
    missing_attrs = [name for name in PAIR_REQUIRED_ATTRS if name not in pair.attrs]
    if missing_attrs:
        raise ValueError(f"{context}: pair group is missing attrs {missing_attrs}")
    pair_tic = _strict_int(pair.attrs["tic"], context=f"{context} pair tic")
    if pair_tic != int(table_tic):
        raise ValueError(
            f"{context}: pair TIC {pair_tic} does not match table TIC "
            f"{int(table_tic)}"
        )
    sector = _strict_int(pair.attrs["sector"], context=f"{context} sector")
    camera = _strict_int(pair.attrs["camera"], context=f"{context} camera")
    ccd = _strict_int(pair.attrs["ccd"], context=f"{context} ccd")
    if sector != int(expected_sector):
        raise ValueError(
            f"{context}: pair sector {sector} does not match {int(expected_sector)}"
        )
    if camera not in range(1, 5) or ccd not in range(1, 5):
        raise ValueError(
            f"{context}: invalid detector mapping camera={camera}, ccd={ccd}"
        )
    if "injection_id" in pair.attrs and str(pair.attrs["injection_id"]) != injection_id:
        raise ValueError(
            f"{context}: pair injection_id attr is {pair.attrs['injection_id']!r}"
        )
    for aperture in ("Small", "Primary"):
        candidates = (
            f"injection_baseline_{aperture}",
            f"raw_baseline_{aperture}",
        )
        present = [name for name in candidates if name in pair.attrs]
        if not present:
            raise ValueError(
                f"{context}: pair group lacks a {aperture} injection baseline"
            )
        baseline = float(pair.attrs[present[0]])
        if not np.isfinite(baseline) or baseline < 0:
            raise ValueError(
                f"{context}: {present[0]} must be finite and nonnegative"
            )

    pair_length = _validate_one_dimensional_datasets(
        pair,
        required=PAIR_REQUIRED_DATASETS,
        context=f"{context} pair group",
    )
    if "tic" not in canonical.attrs:
        raise ValueError(f"{context}: canonical group is missing tic attr")
    canonical_tic = _strict_int(
        canonical.attrs["tic"], context=f"{context} canonical tic"
    )
    if canonical_tic != int(table_tic):
        raise ValueError(
            f"{context}: canonical TIC {canonical_tic} does not match table TIC "
            f"{int(table_tic)}"
        )
    if (
        "injection_id" in canonical.attrs
        and str(canonical.attrs["injection_id"]) != injection_id
    ):
        raise ValueError(
            f"{context}: canonical injection_id attr is "
            f"{canonical.attrs['injection_id']!r}"
        )
    _validate_one_dimensional_datasets(
        canonical,
        required=CANONICAL_REQUIRED_DATASETS,
        context=f"{context} canonical group",
        expected_length=pair_length,
    )

    for name in ("cadenceno", "orbitid"):
        if name in canonical and not np.array_equal(
            np.asarray(pair[name]), np.asarray(canonical[name])
        ):
            raise ValueError(
                f"{context}: pair and canonical {name} datasets do not match"
            )
    for name in ("time", "transit_model"):
        if name in canonical and not np.allclose(
            np.asarray(pair[name], dtype=float),
            np.asarray(canonical[name], dtype=float),
            rtol=0.0,
            atol=1.0e-8,
            equal_nan=True,
        ):
            raise ValueError(
                f"{context}: pair and canonical {name} datasets do not match"
            )


def promote_teacher_v3_injection_pairs(
    *,
    selection_table: Path,
    source_pair_h5s: Sequence[Path],
    canonical_injection_h5: Path,
    out_h5: Path,
    repo_root: Path,
    expected_count: int = TEACHER_V3_S56_EXPECTED_INJECTIONS,
    expected_source_counts: Sequence[int] = (
        TEACHER_V3_S56_EXPECTED_SOURCE_COUNTS
    ),
    expected_sector: int = 56,
    overwrite: bool = False,
) -> dict[str, Any]:
    """Copy exactly the selected S56 injection-pair groups into one HDF5.

    Every selected ID must appear in exactly one supplied pair file, at the
    exact source named by the frozen table.  The separate canonical HDF5 must
    contain the corresponding raw-flux group.  Inputs are hashed before and
    after the copy, and the destination is published atomically.
    """

    import h5py

    selection_table = Path(selection_table).resolve()
    repo_root = Path(repo_root).resolve()
    out_h5 = Path(out_h5).resolve()
    source_paths = [Path(path).resolve() for path in source_pair_h5s]
    canonical_path = Path(canonical_injection_h5).resolve()
    if len(source_paths) != 2 or len(set(source_paths)) != 2:
        raise ValueError("exactly two distinct source pair HDF5 files are required")
    if len(expected_source_counts) != len(source_paths):
        raise ValueError(
            "expected_source_counts must have one value per source pair HDF5"
        )
    if int(expected_count) < 1:
        raise ValueError("expected_count must be positive")
    if sum(int(value) for value in expected_source_counts) != int(expected_count):
        raise ValueError("expected_source_counts must sum to expected_count")
    for path in (selection_table, *source_paths, canonical_path):
        if not path.is_file():
            raise FileNotFoundError(f"required input does not exist: {path}")
    if out_h5.exists() and not overwrite:
        raise FileExistsError(f"output exists; pass overwrite=True: {out_h5}")
    if canonical_path == out_h5 or canonical_path in source_paths:
        raise ValueError(
            "canonical injection HDF5 must be distinct from the pair inputs and output"
        )

    table_sha256 = _file_sha256(selection_table)
    source_sha256 = {str(path): _file_sha256(path) for path in source_paths}
    canonical_sha256 = _file_sha256(canonical_path)
    table = _read_table(selection_table)
    required_columns = {
        "native_input_include",
        "injection_id",
        "tic",
        "source_h5",
        "h5_group",
    }
    missing_columns = sorted(required_columns - set(table.columns))
    if missing_columns:
        raise ValueError(
            f"selection table is missing required columns {missing_columns}"
        )
    include = pd.Series(
        [
            _strict_bool(
                value,
                context=f"selection row {index} native_input_include",
            )
            for index, value in table["native_input_include"].items()
        ],
        index=table.index,
        dtype=bool,
    )
    injection_text = table["injection_id"].fillna("").astype(str).str.strip()
    group_text = table["h5_group"].fillna("").astype(str).str.strip()
    injection_like = injection_text.ne("") | group_text.str.startswith(
        "/injections/"
    )
    if "source_kind" in table:
        injection_like |= (
            table["source_kind"]
            .fillna("")
            .astype(str)
            .str.strip()
            .str.lower()
            .str.startswith("injection")
        )
    if "is_injected_row" in table:
        injection_like |= pd.Series(
            [
                _strict_bool(value, context=f"selection row {index} is_injected_row")
                for index, value in table["is_injected_row"].items()
            ],
            index=table.index,
            dtype=bool,
        )
    missing_active_ids = table.index[
        include & injection_like & injection_text.eq("")
    ].tolist()
    if missing_active_ids:
        raise ValueError(
            "active injection rows are missing injection_id values: "
            f"{missing_active_ids[:10]}"
        )
    selected = table.loc[include & injection_like].copy()
    selected["injection_id"] = injection_text.loc[selected.index]
    if len(selected) != int(expected_count):
        raise ValueError(
            f"selected {len(selected)} active injection rows; "
            f"expected {int(expected_count)}"
        )
    duplicate_ids = sorted(
        selected.loc[
            selected["injection_id"].duplicated(keep=False), "injection_id"
        ].unique()
    )
    if duplicate_ids:
        raise ValueError(
            "selection table contains duplicate active injection IDs: "
            f"{duplicate_ids[:10]}"
        )

    source_set = set(source_paths)
    normalized_rows: list[dict[str, Any]] = []
    for index, row in selected.iterrows():
        injection_id = str(row["injection_id"])
        if "/" in injection_id:
            raise ValueError(
                f"selection row {index}: injection_id may not contain '/': "
                f"{injection_id!r}"
            )
        tic = _strict_int(row["tic"], context=f"injection {injection_id} table tic")
        if tic <= 0:
            raise ValueError(f"injection {injection_id}: TIC must be positive")
        source_path = _resolve_path(row["source_h5"], repo_root=repo_root)
        if source_path not in source_set:
            raise ValueError(
                f"injection {injection_id}: table source {source_path} is not "
                "one of the two explicit pair files"
            )
        group_path = str(row["h5_group"]).strip()
        expected_group_path = f"/injections/{injection_id}"
        if group_path != expected_group_path:
            raise ValueError(
                f"injection {injection_id}: h5_group {group_path!r} does not "
                f"match {expected_group_path!r}"
            )
        if "sector" in selected.columns:
            table_sector = _strict_int(
                row["sector"], context=f"injection {injection_id} table sector"
            )
            if table_sector != int(expected_sector):
                raise ValueError(
                    f"injection {injection_id}: table sector {table_sector} "
                    f"does not match {int(expected_sector)}"
                )
        normalized_rows.append(
            {
                "injection_id": injection_id,
                "tic": tic,
                "source_h5": str(source_path),
                "h5_group": expected_group_path,
            }
        )
    normalized_rows.sort(key=lambda row: row["injection_id"])
    selected_rows_sha256 = _normalized_selection_digest(normalized_rows)

    observed_source_counts = {
        str(path): sum(row["source_h5"] == str(path) for row in normalized_rows)
        for path in source_paths
    }
    for path, expected in zip(source_paths, expected_source_counts, strict=True):
        observed = observed_source_counts[str(path)]
        if observed != int(expected):
            raise ValueError(
                f"source {path} contributes {observed} selected groups; "
                f"expected {int(expected)}"
            )

    out_h5.parent.mkdir(parents=True, exist_ok=True)
    temporary = out_h5.with_suffix(out_h5.suffix + ".tmp")
    if temporary.exists():
        raise FileExistsError(f"temporary output already exists: {temporary}")
    created_utc = datetime.now(timezone.utc).isoformat()
    try:
        with (
            h5py.File(source_paths[0], "r") as source0,
            h5py.File(source_paths[1], "r") as source1,
            h5py.File(canonical_path, "r") as canonical,
            h5py.File(temporary, "w") as output,
        ):
            sources = {
                source_paths[0]: source0,
                source_paths[1]: source1,
            }
            for path, handle in sources.items():
                if "injections" not in handle:
                    raise ValueError(f"source pair HDF5 lacks /injections: {path}")
            if "injections" not in canonical:
                raise ValueError(
                    f"canonical injection HDF5 lacks /injections: {canonical_path}"
                )

            output.attrs["contract_version"] = (
                TEACHER_V3_INJECTION_PAIR_MERGE_CONTRACT
            )
            output.attrs["created_utc"] = created_utc
            output.attrs["sector"] = int(expected_sector)
            output.attrs["selection_table"] = str(selection_table)
            output.attrs["selection_table_sha256"] = table_sha256
            output.attrs["selected_rows_sha256"] = selected_rows_sha256
            output.attrs["n_injections"] = int(len(normalized_rows))
            output.attrs["source_pair_h5s"] = json.dumps(
                [str(path) for path in source_paths]
            )
            output.attrs["source_pair_h5_sha256"] = json.dumps(
                source_sha256, sort_keys=True
            )
            output.attrs["source_pair_group_counts"] = json.dumps(
                observed_source_counts, sort_keys=True
            )
            output.attrs["source_injection_h5"] = str(canonical_path)
            output.attrs["source_injection_h5_sha256"] = canonical_sha256
            injection_root = output.create_group("injections")

            for row in normalized_rows:
                injection_id = row["injection_id"]
                group_path = row["h5_group"]
                expected_source = Path(row["source_h5"])
                locations = [
                    path for path, handle in sources.items() if group_path in handle
                ]
                if not locations:
                    raise ValueError(
                        f"injection {injection_id}: group is missing from both "
                        "source pair HDF5 files"
                    )
                if len(locations) != 1:
                    raise ValueError(
                        f"injection {injection_id}: group appears in "
                        f"{len(locations)} source pair HDF5 files"
                    )
                if locations[0] != expected_source:
                    raise ValueError(
                        f"injection {injection_id}: group is in {locations[0]} "
                        f"but the table binds {expected_source}"
                    )
                if group_path not in canonical:
                    raise ValueError(
                        f"injection {injection_id}: canonical HDF5 is missing "
                        f"{group_path}"
                    )
                pair = sources[expected_source][group_path]
                canonical_group = canonical[group_path]
                _validate_pair_group(
                    pair=pair,
                    canonical=canonical_group,
                    injection_id=injection_id,
                    table_tic=int(row["tic"]),
                    expected_sector=int(expected_sector),
                )
                sources[expected_source].copy(
                    pair,
                    injection_root,
                    name=injection_id,
                )
                destination = injection_root[injection_id]
                original_source = destination.attrs.get(
                    "source_injection_h5", ""
                )
                destination.attrs["source_injection_h5_original"] = str(
                    original_source
                )
                destination.attrs["source_injection_h5"] = str(canonical_path)
                destination.attrs["source_injection_h5_sha256"] = (
                    canonical_sha256
                )
                destination.attrs["promoted_source_pair_h5"] = str(
                    expected_source
                )
                destination.attrs["promoted_source_pair_h5_sha256"] = (
                    source_sha256[str(expected_source)]
                )
                destination.attrs["promoted_source_group"] = group_path
                destination.attrs["selection_table_sha256"] = table_sha256
                destination.attrs["selected_rows_sha256"] = selected_rows_sha256

            observed_ids = sorted(output["injections"].keys())
            expected_ids = sorted(
                row["injection_id"] for row in normalized_rows
            )
            if observed_ids != expected_ids:
                raise RuntimeError(
                    "promoted injection groups do not exactly match the selection"
                )
            if str(output.attrs["source_injection_h5"]) != str(canonical_path):
                raise RuntimeError("promoted root canonical path changed")
            for injection_id in expected_ids:
                if (
                    str(
                        output[f"injections/{injection_id}"].attrs[
                            "source_injection_h5"
                        ]
                    )
                    != str(canonical_path)
                ):
                    raise RuntimeError(
                        f"injection {injection_id}: promoted canonical path changed"
                    )

        observed_inputs = {
            "selection_table": _file_sha256(selection_table),
            "canonical_injection_h5": _file_sha256(canonical_path),
            **{
                f"source_pair_h5:{path}": _file_sha256(path)
                for path in source_paths
            },
        }
        expected_inputs = {
            "selection_table": table_sha256,
            "canonical_injection_h5": canonical_sha256,
            **{
                f"source_pair_h5:{path}": source_sha256[str(path)]
                for path in source_paths
            },
        }
        if observed_inputs != expected_inputs:
            raise RuntimeError(
                "one or more selection/pair/canonical inputs changed during merge"
            )
        temporary.replace(out_h5)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise

    output_sha256 = _file_sha256(out_h5)
    return {
        "contract_version": TEACHER_V3_INJECTION_PAIR_MERGE_CONTRACT,
        "created_utc": created_utc,
        "selection_table": str(selection_table),
        "selection_table_sha256": table_sha256,
        "selected_rows_sha256": selected_rows_sha256,
        "source_pair_h5_sha256": source_sha256,
        "source_pair_group_counts": observed_source_counts,
        "canonical_injection_h5": str(canonical_path),
        "canonical_injection_h5_sha256": canonical_sha256,
        "n_injections": int(len(normalized_rows)),
        "out_h5": str(out_h5),
        "out_h5_sha256": output_sha256,
    }


__all__ = [
    "CANONICAL_REQUIRED_DATASETS",
    "PAIR_REQUIRED_ATTRS",
    "PAIR_REQUIRED_DATASETS",
    "TEACHER_V3_INJECTION_PAIR_MERGE_CONTRACT",
    "TEACHER_V3_S56_EXPECTED_INJECTIONS",
    "TEACHER_V3_S56_EXPECTED_SOURCE_COUNTS",
    "promote_teacher_v3_injection_pairs",
]
