"""Shared browser-label table normalization and validation."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd


BASE_LABEL_COLUMNS: tuple[str, ...] = (
    "row_id",
    "candidate_key",
    "tic",
    "sector",
    "label",
    "label_source",
    "labeler",
    "notes",
    "updated_utc",
)

ADJUDICATION_LABEL_COLUMNS: tuple[str, ...] = (
    *BASE_LABEL_COLUMNS,
    "period_factor",
    "period_status",
)


def clean_value(value: Any) -> Any:
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except TypeError:
        pass
    if isinstance(value, np.generic):
        return value.item()
    return value


def _canonical_numeric(value: Any, *, integer: bool = False) -> str:
    cleaned = clean_value(value)
    if cleaned is None or str(cleaned).strip() == "":
        return ""
    try:
        numeric = float(cleaned)
    except (TypeError, ValueError):
        return str(cleaned).strip()
    if not np.isfinite(numeric):
        return ""
    if integer and np.isclose(numeric, round(numeric), rtol=0.0, atol=1.0e-9):
        return str(int(round(numeric)))
    return format(numeric, ".15g")


def candidate_key(row: pd.Series) -> str:
    """Return the cross-pandas stable key written by the browser label app."""

    return "|".join(
        (
            _canonical_numeric(row.get("tic", ""), integer=True),
            _canonical_numeric(row.get("sector", ""), integer=True),
            _canonical_numeric(row.get("period_d", "")),
            _canonical_numeric(row.get("t0_bjd", "")),
            str(clean_value(row.get("source_bucket", "")) or "").strip(),
        )
    )


def canonicalize_candidate_key(value: Any) -> str:
    """Normalize legacy float spellings without weakening candidate identity."""

    text = str(clean_value(value) or "")
    parts = text.split("|")
    if len(parts) != 5:
        return text
    return "|".join(
        (
            _canonical_numeric(parts[0], integer=True),
            _canonical_numeric(parts[1], integer=True),
            _canonical_numeric(parts[2]),
            _canonical_numeric(parts[3]),
            parts[4].strip(),
        )
    )


def normalize_review_queue(frame: pd.DataFrame) -> pd.DataFrame:
    """Apply the exact display-order identity contract used by the browser."""

    out = frame.reset_index(drop=True).copy()
    out = out.drop(columns=["row_id"], errors="ignore")
    out.insert(0, "row_id", np.arange(len(out), dtype=int))
    out["candidate_key"] = out.apply(candidate_key, axis=1)
    if out["row_id"].duplicated().any():
        raise ValueError("normalized review queue contains duplicate row_id values")
    return out


def latest_label_records(
    labels_csv: Path,
    *,
    columns: Iterable[str] = ADJUDICATION_LABEL_COLUMNS,
) -> pd.DataFrame:
    """Return the latest browser label per normalized display row."""

    columns = tuple(columns)
    path = Path(labels_csv)
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=columns)
    labels = pd.read_csv(path)
    for col in columns:
        if col not in labels:
            labels[col] = ""
    if labels.empty:
        return labels.loc[:, columns]
    labels["row_id"] = pd.to_numeric(labels["row_id"], errors="coerce")
    labels = labels.dropna(subset=["row_id"]).copy()
    labels["row_id"] = labels["row_id"].astype(int)
    labels["_updated_sort"] = pd.to_datetime(labels["updated_utc"], errors="coerce")
    labels = labels.sort_values(["row_id", "_updated_sort"], na_position="first", kind="stable")
    return labels.drop_duplicates("row_id", keep="last").loc[:, columns].reset_index(drop=True)


def validate_label_records(queue: pd.DataFrame, labels: pd.DataFrame) -> None:
    """Reject labels that cannot be tied exactly to the normalized queue."""

    if labels.empty:
        return
    if labels["row_id"].duplicated().any():
        raise ValueError("latest label records contain duplicate row_id values")
    queue_by_id = queue.set_index("row_id")
    unknown = sorted(set(labels["row_id"].astype(int)) - set(queue_by_id.index.astype(int)))
    if unknown:
        raise ValueError(f"labels reference unknown normalized row_id values: {unknown[:10]}")
    expected = labels["row_id"].map(queue_by_id["candidate_key"]).fillna("").map(canonicalize_candidate_key)
    observed = labels["candidate_key"].fillna("").map(canonicalize_candidate_key)
    mismatch = ~observed.eq(expected)
    if mismatch.any():
        bad = labels.loc[mismatch, ["row_id", "candidate_key"]].copy()
        bad["expected_candidate_key"] = expected.loc[mismatch]
        raise ValueError(
            "candidate_key mismatch for normalized browser labels: "
            + bad.head(5).to_dict("records").__repr__()
        )


def normalized_queue_and_labels(
    queue: pd.DataFrame,
    labels_csv: Path,
    *,
    columns: Iterable[str] = ADJUDICATION_LABEL_COLUMNS,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    normalized = normalize_review_queue(queue)
    labels = latest_label_records(Path(labels_csv), columns=columns)
    validate_label_records(normalized, labels)
    return normalized, labels


__all__ = [
    "ADJUDICATION_LABEL_COLUMNS",
    "BASE_LABEL_COLUMNS",
    "candidate_key",
    "canonicalize_candidate_key",
    "clean_value",
    "latest_label_records",
    "normalize_review_queue",
    "normalized_queue_and_labels",
    "validate_label_records",
]
