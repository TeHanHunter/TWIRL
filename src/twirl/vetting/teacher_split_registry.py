"""Immutable TIC-grouped split registry for Teacher-v1 training."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import tempfile
from typing import Any, Sequence

import numpy as np
import pandas as pd


TIC_SPLIT_REGISTRY_SCHEMA_VERSION = "twirl_teacher_v1_tic_split_registry_v1"
TIC_SPLIT_POLICY_VERSION = "tic_grouped_stratified_test20_dev5_v1"
DEFAULT_OBSERVATION_IDENTITY_COLUMNS: tuple[str, ...] = (
    "sector",
    "tic",
    "candidate_key",
)
N_DEVELOPMENT_FOLDS = 5
TEST_FRACTION = 0.20
MAX_EXACT_INTEGER = 2**53 - 1
REGISTRY_COLUMNS: tuple[str, ...] = (
    "tic",
    "fixed_split",
    "cv_fold",
    "identity_columns_json",
    "n_observations",
    "observation_record_sha256",
    "split_seed",
    "split_stratum_column",
    "split_policy_version",
)


def _positive_integers(values: pd.Series, *, name: str) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    invalid = (
        ~np.isfinite(numeric)
        | (numeric <= 0)
        | (numeric > MAX_EXACT_INTEGER)
        | ~np.isclose(numeric, np.rint(numeric), rtol=0.0, atol=0.0)
    )
    if invalid.any():
        examples = values.loc[invalid].head(5).tolist()
        raise ValueError(
            f"{name} contains blank, non-integral, non-positive, or "
            f"out-of-range values; first={examples}"
        )
    return pd.Series(
        np.rint(numeric).astype(np.int64),
        index=values.index,
        name=name,
    )


def _normalize_observations(
    observations: pd.DataFrame,
    *,
    stratum_column: str,
    identity_columns: Sequence[str],
) -> pd.DataFrame:
    identity_columns = tuple(identity_columns)
    if not identity_columns or len(set(identity_columns)) != len(identity_columns):
        raise ValueError("identity_columns must be nonempty and unique")
    if "tic" not in identity_columns:
        raise ValueError("identity_columns must contain 'tic'")
    if not str(stratum_column).strip() or stratum_column in identity_columns:
        raise ValueError(
            "stratum_column must be nonblank and separate from identity_columns"
        )
    missing = sorted(
        {*identity_columns, stratum_column} - set(observations.columns)
    )
    if missing:
        raise KeyError(f"observation corpus is missing columns: {missing}")
    if observations.empty:
        raise ValueError("observation corpus is empty")

    work = observations.copy()
    work["tic"] = _positive_integers(work["tic"], name="tic")
    if "sector" in identity_columns:
        work["sector"] = _positive_integers(work["sector"], name="sector")
    for column in identity_columns:
        if column in {"tic", "sector"}:
            continue
        work[column] = work[column].fillna("").astype(str).str.strip()
        if work[column].eq("").any():
            raise ValueError(f"{column} contains blank identities")
    work[stratum_column] = (
        work[stratum_column].fillna("").astype(str).str.strip()
    )
    if work[stratum_column].eq("").any():
        raise ValueError(f"{stratum_column} contains blank strata")

    duplicate = work.duplicated(list(identity_columns), keep=False)
    if duplicate.any():
        examples = (
            work.loc[duplicate, list(identity_columns)]
            .head(5)
            .to_dict(orient="records")
        )
        raise ValueError(
            "observation corpus contains duplicate identities; "
            f"first={examples}"
        )
    if {"sector", "tic", "candidate_key"}.issubset(work.columns):
        key_counts = work.groupby("candidate_key", sort=False).agg(
            n_sectors=("sector", "nunique"),
            n_tics=("tic", "nunique"),
        )
        inconsistent = key_counts[
            key_counts["n_sectors"].gt(1) | key_counts["n_tics"].gt(1)
        ]
        if not inconsistent.empty:
            raise ValueError(
                "candidate_key maps to multiple sector/TIC identities; "
                f"first={inconsistent.head(5).index.tolist()}"
            )
    return work.reset_index(drop=True)


def validate_tic_split_assignments(
    rows: pd.DataFrame,
    *,
    require_unique_tics: bool = False,
    require_complete_partitions: bool = True,
) -> pd.DataFrame:
    """Normalize assignments and reject invalid values or TIC leakage."""

    missing = sorted({"tic", "fixed_split", "cv_fold"} - set(rows.columns))
    if missing:
        raise KeyError(f"split assignments are missing columns: {missing}")
    if rows.empty:
        raise ValueError("split assignments are empty")
    work = rows.loc[:, ["tic", "fixed_split", "cv_fold"]].copy()
    work["tic"] = _positive_integers(work["tic"], name="tic")
    if require_unique_tics and work["tic"].duplicated().any():
        raise ValueError("TIC split registry contains duplicate TIC keys")

    work["fixed_split"] = (
        work["fixed_split"].fillna("").astype(str).str.strip().str.lower()
    )
    invalid_split = ~work["fixed_split"].isin({"development", "test"})
    if invalid_split.any():
        values = sorted(work.loc[invalid_split, "fixed_split"].unique())
        raise ValueError(f"invalid fixed_split values: {values}")

    folds = pd.to_numeric(work["cv_fold"], errors="coerce").to_numpy(
        dtype=float
    )
    invalid_fold = ~np.isfinite(folds) | ~np.isclose(
        folds, np.rint(folds), rtol=0.0, atol=0.0
    )
    if invalid_fold.any():
        raise ValueError("cv_fold contains invalid integers")
    test = work["fixed_split"].eq("test").to_numpy()
    development = ~test
    if np.any(test & ~np.isclose(folds, -1.0, rtol=0.0, atol=0.0)):
        raise ValueError("test TICs must have cv_fold=-1")
    if np.any(
        development
        & ((folds < 0) | (folds >= N_DEVELOPMENT_FOLDS))
    ):
        raise ValueError("development TICs must have cv_fold in [0,4]")
    work["cv_fold"] = np.rint(folds).astype(np.int16)

    counts = work.groupby("tic", sort=False).agg(
        n_splits=("fixed_split", "nunique"),
        n_folds=("cv_fold", "nunique"),
    )
    leaking = counts["n_splits"].ne(1) | counts["n_folds"].ne(1)
    if leaking.any():
        raise ValueError(
            "TIC group leakage across fixed_split/cv_fold assignments; "
            f"first={counts.index[leaking].tolist()[:5]}"
        )
    if require_complete_partitions:
        if set(work["fixed_split"]) != {"development", "test"}:
            raise ValueError(
                "split assignments must contain development and test TICs"
            )
        observed_folds = set(
            work.loc[work["fixed_split"].eq("development"), "cv_fold"].astype(
                int
            )
        )
        if observed_folds != set(range(N_DEVELOPMENT_FOLDS)):
            raise ValueError(
                "development assignments must contain all five CV folds; "
                f"observed={sorted(observed_folds)}"
            )
    return (
        work.drop_duplicates("tic")
        .sort_values("tic", kind="stable")
        .reset_index(drop=True)
    )


def _target_test_tic_count(n_tics: int) -> int:
    if n_tics < N_DEVELOPMENT_FOLDS + 1:
        raise ValueError(
            "at least six unique TICs are required for a test set and five "
            "nonempty development folds"
        )
    nearest = max(1, int(np.floor(TEST_FRACTION * n_tics + 0.5)))
    return min(nearest, n_tics - N_DEVELOPMENT_FOLDS)


def _balanced_capacities(total: int, n_folds: int) -> tuple[int, ...]:
    base, remainder = divmod(total, n_folds)
    return tuple(base + int(fold < remainder) for fold in range(n_folds))


def _assign_groups(
    rows: pd.DataFrame,
    *,
    stratum_column: str,
    capacities: Sequence[int],
    seed: int,
) -> pd.Series:
    """Greedily balance strata under exact whole-TIC fold capacities."""

    capacities = np.asarray(capacities, dtype=int)
    groups = rows["tic"].astype(str)
    table = pd.crosstab(groups, rows[stratum_column]).sort_index()
    if capacities.ndim != 1 or len(capacities) < 2:
        raise ValueError("at least two fold capacities are required")
    if np.any(capacities < 0) or capacities.sum() != len(table):
        raise ValueError("fold capacities must cover every unique TIC")

    totals = table.sum(axis=0).to_numpy(dtype=float)
    fractions = capacities / float(len(table))
    target_labels = fractions[:, None] * totals[None, :]
    target_rows = fractions * len(rows)
    rarity = (
        table.to_numpy(dtype=float) / np.maximum(totals, 1.0)
    ).max(axis=1)
    rng = np.random.default_rng(seed)
    order = np.lexsort(
        (
            rng.random(len(table)),
            -table.sum(axis=1).to_numpy(dtype=float),
            -rarity,
        )
    )
    label_counts = np.zeros_like(target_labels)
    row_counts = np.zeros(len(capacities), dtype=float)
    group_counts = np.zeros(len(capacities), dtype=int)
    assignment: dict[str, int] = {}
    for position in order:
        vector = table.iloc[position].to_numpy(dtype=float)
        size = float(vector.sum())
        choices = np.flatnonzero(group_counts < capacities)
        costs: list[float] = []
        for fold in choices:
            candidate_labels = label_counts.copy()
            candidate_rows = row_counts.copy()
            candidate_labels[fold] += vector
            candidate_rows[fold] += size
            label_cost = np.sum(
                np.square(candidate_labels - target_labels)
                / np.maximum(target_labels, 1.0)
            )
            row_cost = np.sum(
                np.square(candidate_rows - target_rows)
                / np.maximum(target_rows, 1.0)
            )
            costs.append(float(label_cost + 0.1 * row_cost))
        best = np.min(costs)
        finalists = choices[np.isclose(costs, best)]
        fold = int(rng.choice(finalists))
        key = str(table.index[position])
        assignment[key] = fold
        label_counts[fold] += vector
        row_counts[fold] += size
        group_counts[fold] += 1
    return groups.map(assignment).astype(np.int16)


def _frame_sha256(frame: pd.DataFrame, columns: Sequence[str]) -> str:
    canonical = frame.loc[:, list(columns)].sort_values(
        list(columns), kind="stable"
    )
    digest = hashlib.sha256()
    for record in canonical.to_dict(orient="records"):
        digest.update(
            json.dumps(
                record,
                sort_keys=True,
                separators=(",", ":"),
                ensure_ascii=False,
            ).encode("utf-8")
        )
        digest.update(b"\n")
    return digest.hexdigest()


def _tic_observation_metadata(
    observations: pd.DataFrame,
    *,
    identity_columns: Sequence[str],
    stratum_column: str,
) -> pd.DataFrame:
    records = []
    hash_columns = (*identity_columns, stratum_column)
    for tic, group in observations.groupby("tic", sort=True):
        records.append(
            {
                "tic": int(tic),
                "n_observations": int(len(group)),
                "observation_record_sha256": _frame_sha256(
                    group, hash_columns
                ),
            }
        )
    return pd.DataFrame.from_records(records)


def validate_tic_split_registry(registry: pd.DataFrame) -> pd.DataFrame:
    """Validate the immutable one-row-per-TIC registry contract."""

    missing = sorted(set(REGISTRY_COLUMNS) - set(registry.columns))
    if missing:
        raise KeyError(f"TIC split registry is missing columns: {missing}")
    work = registry.loc[:, REGISTRY_COLUMNS].copy()
    assignments = validate_tic_split_assignments(
        work,
        require_unique_tics=True,
        require_complete_partitions=True,
    )
    expected_test = _target_test_tic_count(len(assignments))
    observed_test = int(assignments["fixed_split"].eq("test").sum())
    if observed_test != expected_test:
        raise ValueError(
            "registry must hold the nearest feasible 20% of unique TICs in "
            f"test; expected={expected_test}, observed={observed_test}"
        )

    work["tic"] = _positive_integers(work["tic"], name="tic")
    work["n_observations"] = _positive_integers(
        work["n_observations"], name="n_observations"
    )
    seeds = pd.to_numeric(work["split_seed"], errors="coerce").to_numpy(
        dtype=float
    )
    invalid_seed = (
        ~np.isfinite(seeds)
        | (seeds < 0)
        | (seeds > MAX_EXACT_INTEGER)
        | ~np.isclose(seeds, np.rint(seeds), rtol=0.0, atol=0.0)
    )
    if invalid_seed.any() or len(set(seeds)) != 1:
        raise ValueError("split_seed must be one non-negative integer")
    work["split_seed"] = np.rint(seeds).astype(np.int64)

    for column in (
        "identity_columns_json",
        "split_stratum_column",
        "split_policy_version",
    ):
        work[column] = work[column].fillna("").astype(str).str.strip()
        if work[column].eq("").any() or work[column].nunique() != 1:
            raise ValueError(f"{column} must contain one nonblank value")
    if not work["split_policy_version"].eq(TIC_SPLIT_POLICY_VERSION).all():
        raise ValueError(
            f"split_policy_version must equal {TIC_SPLIT_POLICY_VERSION!r}"
        )
    try:
        identity_columns = json.loads(work["identity_columns_json"].iloc[0])
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        raise ValueError("identity_columns_json is invalid") from exc
    if (
        not isinstance(identity_columns, list)
        or not identity_columns
        or "tic" not in identity_columns
        or len(identity_columns) != len(set(identity_columns))
        or any(not isinstance(value, str) or not value for value in identity_columns)
    ):
        raise ValueError("identity_columns_json is invalid")
    work["identity_columns_json"] = json.dumps(
        identity_columns, separators=(",", ":")
    )
    work["observation_record_sha256"] = (
        work["observation_record_sha256"]
        .fillna("")
        .astype(str)
        .str.lower()
    )
    if (
        ~work["observation_record_sha256"].str.fullmatch(r"[0-9a-f]{64}")
    ).any():
        raise ValueError("observation_record_sha256 contains invalid digests")

    work = work.merge(
        assignments,
        on="tic",
        validate="one_to_one",
        suffixes=("", "_normalized"),
    )
    work["fixed_split"] = work.pop("fixed_split_normalized")
    work["cv_fold"] = work.pop("cv_fold_normalized")
    return (
        work.loc[:, REGISTRY_COLUMNS]
        .sort_values("tic", kind="stable")
        .reset_index(drop=True)
    )


def build_tic_split_registry(
    observations: pd.DataFrame,
    *,
    seed: int,
    stratum_column: str = "split_stratum",
    identity_columns: Sequence[str] = DEFAULT_OBSERVATION_IDENTITY_COLUMNS,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Build a deterministic 20%-test/five-fold registry keyed by TIC."""

    try:
        seed = int(str(seed).strip())
    except (TypeError, ValueError) as exc:
        raise ValueError("seed must be a non-negative integer") from exc
    if seed < 0 or seed > MAX_EXACT_INTEGER:
        raise ValueError("seed must be a non-negative integer")
    identity_columns = tuple(identity_columns)
    work = _normalize_observations(
        observations,
        stratum_column=stratum_column,
        identity_columns=identity_columns,
    )

    n_tics = int(work["tic"].nunique())
    n_test_tics = _target_test_tic_count(n_tics)
    outer = _assign_groups(
        work,
        stratum_column=stratum_column,
        capacities=(
            n_test_tics,
            *_balanced_capacities(n_tics - n_test_tics, 4),
        ),
        seed=seed,
    )
    test = outer.eq(0)
    generated = pd.DataFrame(
        {
            "tic": work["tic"],
            "fixed_split": np.where(test, "test", "development"),
            "cv_fold": np.full(len(work), -1, dtype=np.int16),
        }
    )
    development = work.loc[~test]
    dev_folds = _assign_groups(
        development,
        stratum_column=stratum_column,
        capacities=_balanced_capacities(
            int(development["tic"].nunique()), N_DEVELOPMENT_FOLDS
        ),
        seed=seed + 1,
    )
    generated.loc[development.index, "cv_fold"] = dev_folds.to_numpy(
        dtype=np.int16
    )
    assignments = validate_tic_split_assignments(generated)

    has_fixed = "fixed_split" in work.columns
    has_fold = "cv_fold" in work.columns
    if has_fixed != has_fold:
        raise ValueError(
            "preassigned split columns are incomplete; fixed_split and cv_fold "
            "must be supplied together"
        )
    if has_fixed:
        declared = validate_tic_split_assignments(work)
        comparison = assignments.merge(
            declared,
            on="tic",
            validate="one_to_one",
            suffixes=("_generated", "_declared"),
        )
        mismatch = (
            comparison["fixed_split_generated"].ne(
                comparison["fixed_split_declared"]
            )
            | comparison["cv_fold_generated"].ne(
                comparison["cv_fold_declared"]
            )
        )
        if mismatch.any():
            raise ValueError(
                "preassigned splits disagree with the deterministic declared "
                "seed/stratum contract"
            )

    registry = assignments.merge(
        _tic_observation_metadata(
            work,
            identity_columns=identity_columns,
            stratum_column=stratum_column,
        ),
        on="tic",
        validate="one_to_one",
    )
    registry["identity_columns_json"] = json.dumps(
        list(identity_columns), separators=(",", ":")
    )
    registry["split_seed"] = seed
    registry["split_stratum_column"] = stratum_column
    registry["split_policy_version"] = TIC_SPLIT_POLICY_VERSION
    registry = validate_tic_split_registry(registry)

    observed = work.loc[:, ["tic"]].merge(
        registry.loc[:, ["tic", "fixed_split", "cv_fold"]],
        on="tic",
        validate="many_to_one",
    )
    development_registry = registry["fixed_split"].eq("development")
    test_observations = int(observed["fixed_split"].eq("test").sum())
    summary: dict[str, Any] = {
        "passed": True,
        "schema_version": TIC_SPLIT_REGISTRY_SCHEMA_VERSION,
        "split_policy_version": TIC_SPLIT_POLICY_VERSION,
        "seed": seed,
        "group_column": "tic",
        "stratum_column": stratum_column,
        "identity_columns": list(identity_columns),
        "target_test_fraction": TEST_FRACTION,
        "test_fraction_unit": "unique_tics",
        "target_test_tic_count": n_test_tics,
        "n_development_folds": N_DEVELOPMENT_FOLDS,
        "n_observations": int(len(work)),
        "n_unique_tics": int(len(registry)),
        "n_test_observations": test_observations,
        "n_test_tics": n_test_tics,
        "test_observation_fraction": test_observations / float(len(work)),
        "test_tic_fraction": n_test_tics / float(len(registry)),
        "development_fold_observation_counts": {
            str(fold): int(
                (
                    observed["fixed_split"].eq("development")
                    & observed["cv_fold"].eq(fold)
                ).sum()
            )
            for fold in range(N_DEVELOPMENT_FOLDS)
        },
        "development_fold_tic_counts": {
            str(fold): int(
                (
                    development_registry & registry["cv_fold"].eq(fold)
                ).sum()
            )
            for fold in range(N_DEVELOPMENT_FOLDS)
        },
        "stratum_observation_counts": {
            str(key): int(value)
            for key, value in work[stratum_column]
            .value_counts()
            .sort_index()
            .items()
        },
        "observation_identity_content_sha256": _frame_sha256(
            work, (*identity_columns, stratum_column)
        ),
        "registry_content_sha256": _frame_sha256(
            registry, REGISTRY_COLUMNS
        ),
        "group_leakage_count": 0,
    }
    return registry, summary


def attach_tic_split_registry(
    observations: pd.DataFrame,
    registry: pd.DataFrame,
) -> pd.DataFrame:
    """Attach a registry only when it exactly matches the bound corpus."""

    checked = validate_tic_split_registry(registry)
    identity_columns = tuple(
        json.loads(checked["identity_columns_json"].iloc[0])
    )
    stratum_column = str(checked["split_stratum_column"].iloc[0])
    work = _normalize_observations(
        observations,
        stratum_column=stratum_column,
        identity_columns=identity_columns,
    )
    observed_metadata = _tic_observation_metadata(
        work,
        identity_columns=identity_columns,
        stratum_column=stratum_column,
    ).sort_values("tic", kind="stable").reset_index(drop=True)
    expected_metadata = checked.loc[
        :, ["tic", "n_observations", "observation_record_sha256"]
    ]
    if not observed_metadata["tic"].equals(expected_metadata["tic"]):
        raise ValueError(
            "observation TIC coverage disagrees with the immutable registry"
        )
    if not observed_metadata["n_observations"].equals(
        expected_metadata["n_observations"]
    ):
        raise ValueError(
            "observation counts disagree with the immutable TIC registry"
        )
    if not observed_metadata["observation_record_sha256"].equals(
        expected_metadata["observation_record_sha256"]
    ):
        raise ValueError(
            "observations disagree with the identities/strata bound to the "
            "immutable TIC registry"
        )

    assignments = checked.loc[:, ["tic", "fixed_split", "cv_fold"]]
    has_fixed = "fixed_split" in work.columns
    has_fold = "cv_fold" in work.columns
    if has_fixed != has_fold:
        raise ValueError(
            "preassigned split columns are incomplete; fixed_split and cv_fold "
            "must be supplied together"
        )
    if has_fixed:
        declared = validate_tic_split_assignments(work)
        try:
            pd.testing.assert_frame_equal(
                declared,
                assignments,
                check_dtype=False,
            )
        except AssertionError as exc:
            raise ValueError(
                "existing observation splits disagree with the immutable TIC "
                "registry"
            ) from exc
        work = work.drop(columns=["fixed_split", "cv_fold"])
    attached = work.merge(
        assignments,
        on="tic",
        how="left",
        validate="many_to_one",
    )
    validate_tic_split_assignments(attached)
    return attached


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _publish_immutable(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != payload:
            raise FileExistsError(
                f"refusing to replace immutable output with different bytes: {path}"
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


def write_tic_split_registry(
    *,
    corpus_path: Path,
    registry_path: Path,
    summary_path: Path,
    seed: int,
    stratum_column: str = "split_stratum",
    identity_columns: Sequence[str] = DEFAULT_OBSERVATION_IDENTITY_COLUMNS,
) -> dict[str, Any]:
    """Build and byte-idempotently publish one CSV registry/JSON summary."""

    corpus_path = Path(corpus_path)
    registry_path = Path(registry_path)
    summary_path = Path(summary_path)
    if registry_path.suffix.lower() != ".csv":
        raise ValueError("registry_path must use CSV")
    if summary_path.suffix.lower() != ".json":
        raise ValueError("summary_path must use JSON")
    if len(
        {
            corpus_path.resolve(),
            registry_path.resolve(),
            summary_path.resolve(),
        }
    ) != 3:
        raise ValueError("corpus, registry, and summary paths must be distinct")

    input_hash = _file_sha256(corpus_path)
    suffix = corpus_path.suffix.lower()
    if suffix == ".csv":
        observations = pd.read_csv(corpus_path, low_memory=False)
    elif suffix in {".parquet", ".pq"}:
        observations = pd.read_parquet(corpus_path)
    else:
        raise ValueError(f"unsupported observation table: {corpus_path}")
    registry, summary = build_tic_split_registry(
        observations,
        seed=seed,
        stratum_column=stratum_column,
        identity_columns=identity_columns,
    )
    if _file_sha256(corpus_path) != input_hash:
        raise RuntimeError("observation corpus changed while splits were built")

    registry_payload = registry.to_csv(
        index=False, lineterminator="\n"
    ).encode("utf-8")
    registry_hash = hashlib.sha256(registry_payload).hexdigest()
    summary = {
        **summary,
        "input_corpus": {
            "path": str(corpus_path.resolve()),
            "sha256": input_hash,
        },
        "output_registry": {
            "path": str(registry_path.resolve()),
            "sha256": registry_hash,
        },
    }
    summary_payload = (
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")

    for path, payload in (
        (registry_path, registry_payload),
        (summary_path, summary_payload),
    ):
        if path.exists() and path.read_bytes() != payload:
            raise FileExistsError(
                f"refusing to replace immutable output with different bytes: {path}"
            )
    _publish_immutable(registry_path, registry_payload)
    if _file_sha256(registry_path) != registry_hash:
        raise RuntimeError("published split registry hash mismatch")
    _publish_immutable(summary_path, summary_payload)
    return summary


__all__ = [
    "DEFAULT_OBSERVATION_IDENTITY_COLUMNS",
    "TIC_SPLIT_POLICY_VERSION",
    "attach_tic_split_registry",
    "build_tic_split_registry",
    "validate_tic_split_assignments",
    "validate_tic_split_registry",
    "write_tic_split_registry",
]
