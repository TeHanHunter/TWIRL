"""Blinded, rank-1 multisector enrichment queues for Franklin review."""
from __future__ import annotations

from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from twirl.vetting.label_io import candidate_key
from twirl.vetting.teacher_v2_active_learning import (
    EnrichmentQuotas,
    PUBLIC_COLUMNS,
    build_existing_teacher_enrichment_batch,
)


FRANKLIN_MULTISECTOR_POLICY_VERSION = "s57_s59_franklin_rank1_teacher_v1"
FRANKLIN_LABEL_RETURN_POLICY_VERSION = "franklin_morphology_accept_harmonic_mask_v1"

_IDENTITY_COLUMNS: tuple[str, ...] = (
    "review_id",
    "tic",
    "sector",
    "period_d",
    "t0_bjd",
    "duration_min",
)
_MORPHOLOGY_CLASSES: tuple[str, ...] = (
    "planet_like",
    "eclipse_contact",
    "smooth_variable",
    "other",
)
_FRANKLIN_LABELS: frozenset[str] = frozenset(
    {
        "planet_like",
        "wide_transit_like",
        "eclipsing_binary_or_pceb",
        "stellar_variability",
        "instrumental_or_systematic",
        "uncertain",
        "skip",
    }
)


def _string_value(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except TypeError:
        pass
    return str(value)


def standalone_app_candidate_key(row: Mapping[str, Any] | pd.Series) -> str:
    """Return the identity written by the standalone Franklin browser app.

    This intentionally differs from :func:`twirl.vetting.label_io.candidate_key`.
    Exact string preservation matters, so production callers should read both
    queue and label CSV files with ``dtype=str`` and ``keep_default_na=False``.
    """

    return "|".join(
        _string_value(row.get(column, ""))
        for column in ("review_id", "tic", "sector", "period_d", "t0_bjd")
    )


def normalize_franklin_label_return(
    queue: pd.DataFrame,
    labels: pd.DataFrame,
    *,
    source_batch_id: str,
    morphology_adjudicator: str,
    morphology_accepted_utc: str,
    expected_sector_counts: Mapping[int, int] | None = None,
    native_h5_by_sector: Mapping[int, Path | str] | None = None,
    expected_labeler: str = "franklin",
    expected_label_source: str = "human",
) -> pd.DataFrame:
    """Strictly join one completed Franklin return to its exact frozen queue.

    Franklin's labels are accepted at the sector-observation morphology level.
    Reported period factors/statuses are retained as raw audit metadata but are
    deliberately excluded from harmonic supervision unless a later, explicit
    factor-only review produces a separate verified decision.
    """

    required_queue = {
        "row_id",
        "review_id",
        "tic",
        "sector",
        "period_d",
        "t0_bjd",
        "duration_min",
        "rep_peak_rank",
        "source_kind",
        "candidate_key",
    }
    required_labels = {
        "row_id",
        "candidate_key",
        "tic",
        "sector",
        "label",
        "label_source",
        "labeler",
        "notes",
        "period_factor",
        "period_status",
        "updated_utc",
    }
    missing_queue = sorted(required_queue - set(queue.columns))
    missing_labels = sorted(required_labels - set(labels.columns))
    if missing_queue:
        raise KeyError(f"Franklin queue is missing columns: {missing_queue}")
    if missing_labels:
        raise KeyError(f"Franklin labels are missing columns: {missing_labels}")
    if not source_batch_id.strip():
        raise ValueError("source_batch_id must be nonblank")
    if not morphology_adjudicator.strip():
        raise ValueError("morphology_adjudicator must be nonblank")
    if not morphology_accepted_utc.strip():
        raise ValueError("morphology_accepted_utc must be nonblank")
    if not expected_labeler.strip() or not expected_label_source.strip():
        raise ValueError("expected labeler/source must be nonblank")

    public = queue.copy()
    returned = labels.copy()
    for name, frame in (("queue", public), ("labels", returned)):
        row_ids = frame["row_id"].map(_string_value)
        if row_ids.eq("").any() or row_ids.duplicated().any():
            raise ValueError(f"Franklin {name} contains blank or duplicate row_id values")
        frame["row_id"] = row_ids
    if set(public["row_id"]) != set(returned["row_id"]):
        missing = sorted(set(public["row_id"]) - set(returned["row_id"]))[:5]
        extra = sorted(set(returned["row_id"]) - set(public["row_id"]))[:5]
        raise ValueError(
            "Franklin row coverage differs between queue and labels: "
            f"missing={missing}, extra={extra}"
        )

    public["standalone_app_candidate_key"] = public.apply(
        standalone_app_candidate_key, axis=1
    )
    if public["standalone_app_candidate_key"].duplicated().any():
        raise ValueError("Franklin queue contains duplicate standalone app keys")
    returned_key = returned["candidate_key"].map(_string_value)
    if returned_key.eq("").any() or returned_key.duplicated().any():
        raise ValueError("Franklin labels contain blank or duplicate candidate_key values")
    returned["standalone_app_candidate_key"] = returned_key

    public = public.rename(columns={"candidate_key": "pipeline_candidate_key"})
    returned = returned.rename(
        columns={
            "tic": "returned_tic",
            "sector": "returned_sector",
            "label": "human_label",
            "label_source": "human_label_source",
            "labeler": "human_labeler",
            "notes": "human_notes",
            "period_factor": "reported_period_factor",
            "period_status": "reported_period_status",
            "updated_utc": "human_updated_utc",
        }
    )
    keep = [
        "row_id",
        "standalone_app_candidate_key",
        "returned_tic",
        "returned_sector",
        "human_label",
        "human_label_source",
        "human_labeler",
        "human_notes",
        "reported_period_factor",
        "reported_period_status",
        "human_updated_utc",
    ]
    joined = public.merge(
        returned.loc[:, keep],
        on="row_id",
        how="left",
        validate="one_to_one",
        suffixes=("", "_returned"),
    )
    key_mismatch = joined["standalone_app_candidate_key"].ne(
        joined["standalone_app_candidate_key_returned"]
    )
    tic_mismatch = joined["tic"].map(_string_value).ne(
        joined["returned_tic"].map(_string_value)
    )
    sector_mismatch = joined["sector"].map(_string_value).ne(
        joined["returned_sector"].map(_string_value)
    )
    mismatch = key_mismatch | tic_mismatch | sector_mismatch
    if mismatch.any():
        sample = joined.loc[
            mismatch,
            [
                "row_id",
                "standalone_app_candidate_key",
                "standalone_app_candidate_key_returned",
            ],
        ].head(3)
        raise ValueError(
            "Franklin labels do not match the exact frozen queue; "
            f"first={sample.to_dict(orient='records')}"
        )
    joined = joined.drop(
        columns=[
            "standalone_app_candidate_key_returned",
            "returned_tic",
            "returned_sector",
        ]
    )

    labels_seen = set(joined["human_label"].fillna("").astype(str))
    invalid_labels = sorted(labels_seen - _FRANKLIN_LABELS)
    if invalid_labels or "" in labels_seen:
        raise ValueError(f"Franklin return contains invalid or blank labels: {invalid_labels}")
    observed_labelers = set(joined["human_labeler"].fillna("").astype(str))
    if observed_labelers != {expected_labeler}:
        raise ValueError(
            "Franklin labeler differs: "
            f"observed={sorted(observed_labelers)}, expected={expected_labeler!r}"
        )
    observed_label_sources = set(
        joined["human_label_source"].fillna("").astype(str)
    )
    if observed_label_sources != {expected_label_source}:
        raise ValueError(
            "Franklin label source differs: "
            f"observed={sorted(observed_label_sources)}, "
            f"expected={expected_label_source!r}"
        )
    if not pd.to_numeric(joined["rep_peak_rank"], errors="coerce").eq(1).all():
        raise ValueError("Franklin return contains a non-rank-1 review ephemeris")
    if set(joined["source_kind"].fillna("").astype(str)) != {"real_candidate"}:
        raise ValueError("Franklin return must be real-candidate only")

    observed_sector_counts = {
        int(key): int(value)
        for key, value in pd.to_numeric(
            joined["sector"], errors="raise"
        ).value_counts().sort_index().items()
    }
    if expected_sector_counts is not None:
        expected = {int(key): int(value) for key, value in expected_sector_counts.items()}
        if observed_sector_counts != expected:
            raise ValueError(
                "Franklin sector counts differ: "
                f"observed={observed_sector_counts}, expected={expected}"
            )

    joined["source_uid"] = (
        source_batch_id + ":" + joined["standalone_app_candidate_key"].astype(str)
    )
    joined["source_batch_id"] = source_batch_id
    joined["label_unit"] = "sector_observation"
    joined["is_injected_row"] = False
    joined["morphology_review_status"] = "accepted_batch_level"
    joined["morphology_adjudicator"] = morphology_adjudicator
    joined["morphology_accepted_utc"] = morphology_accepted_utc
    joined["factor_review_status"] = "not_explicitly_reviewed"
    joined["harmonic_supervision_verified"] = False
    joined["label_return_policy_version"] = FRANKLIN_LABEL_RETURN_POLICY_VERSION
    native_paths = {int(key): str(value) for key, value in (native_h5_by_sector or {}).items()}
    joined["native_h5_path"] = pd.to_numeric(
        joined["sector"], errors="raise"
    ).map(native_paths).fillna("")

    from twirl.vetting.adjudication_audit import add_harmonic_cnn_targets

    joined = add_harmonic_cnn_targets(joined)
    if joined["harmonic_include_v1"].astype(bool).any():
        raise AssertionError("unverified Franklin factors entered harmonic supervision")
    if joined["harmonic_target_v1"].fillna("").astype(str).ne("").any():
        raise AssertionError("unverified Franklin factors produced harmonic targets")
    return joined.reset_index(drop=True)


def _read_table(path: Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path, low_memory=False)
    if path.suffix.lower() == ".parquet":
        return pd.read_parquet(path)
    raise ValueError(f"unsupported table format: {path}")


def _write_table(frame: pd.DataFrame, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix.lower() == ".csv":
        frame.to_csv(path, index=False, float_format="%.15g")
    elif path.suffix.lower() == ".parquet":
        frame.to_parquet(path, compression="zstd", index=False)
    else:
        raise ValueError(f"unsupported table format: {path}")
    return path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _sum_quotas(values: Sequence[EnrichmentQuotas]) -> EnrichmentQuotas:
    fields = tuple(asdict(EnrichmentQuotas()))
    return EnrichmentQuotas(
        **{
            field: int(sum(getattr(value, field) for value in values))
            for field in fields
        }
    )


def prepare_single_teacher_scores(scores: pd.DataFrame) -> pd.DataFrame:
    """Adapt one Teacher-v1 score table to the shared enrichment selector.

    Teacher v1 has no separate compact head. Its compact score is therefore a
    conservative geometric mean of the Planet-like morphology and preserve
    probabilities. This derived score is selection provenance, never a label.
    """

    required = {
        *_IDENTITY_COLUMNS,
        "tmag",
        "sde_max",
        "rep_peak_rank",
        "p_preserve",
        "std_p_preserve",
        *(f"p_{label}" for label in _MORPHOLOGY_CLASSES),
        *(f"std_p_{label}" for label in _MORPHOLOGY_CLASSES),
    }
    missing = sorted(required - set(scores.columns))
    if missing:
        raise KeyError(f"single-teacher score table is missing columns: {missing}")
    work = scores.copy()
    if work["review_id"].fillna("").astype(str).duplicated().any():
        raise ValueError("single-teacher score table contains duplicate review_id values")
    if "source_kind" in work:
        injected = work["source_kind"].fillna("").astype(str).str.contains(
            "inject", case=False
        )
        if injected.any():
            raise ValueError("Franklin enrichment score table contains injected rows")
    if "is_injected_row" in work:
        injected = work["is_injected_row"].fillna(False).astype(bool)
        if injected.any():
            raise ValueError("Franklin enrichment score table contains injected rows")
    work["source_kind"] = "real_candidate"
    work["is_injected_row"] = False

    probability_columns = [f"p_{label}" for label in _MORPHOLOGY_CLASSES]
    probabilities = (
        work.loc[:, probability_columns]
        .apply(pd.to_numeric, errors="coerce")
        .fillna(0.0)
        .clip(0.0, 1.0)
    )
    preserve = (
        pd.to_numeric(work["p_preserve"], errors="coerce")
        .fillna(0.0)
        .clip(0.0, 1.0)
    )
    work["p_compact_transit"] = np.sqrt(
        probabilities["p_planet_like"].to_numpy(float) * preserve.to_numpy(float)
    )
    planet_std = (
        pd.to_numeric(work["std_p_planet_like"], errors="coerce")
        .fillna(0.0)
        .clip(lower=0.0)
    )
    preserve_std = (
        pd.to_numeric(work["std_p_preserve"], errors="coerce")
        .fillna(0.0)
        .clip(lower=0.0)
    )
    work["std_p_compact_transit"] = np.maximum(planet_std, preserve_std)
    work["morphology_entropy"] = (
        -(probabilities * np.log(probabilities.clip(lower=1.0e-12))).sum(axis=1)
        / np.log(len(probability_columns))
    )
    sorted_probability = np.sort(probabilities.to_numpy(float), axis=1)
    work["morphology_margin"] = (
        sorted_probability[:, -1] - sorted_probability[:, -2]
    )
    work["compact_profile_disagreement"] = 0.0
    work["compact_model_profile"] = work.get(
        "model_profile", pd.Series("shape_plus_periodogram_bls", index=work.index)
    )
    work["morphology_model_profile"] = work["compact_model_profile"]
    work["compact_model_version"] = work.get(
        "model_version", pd.Series("s56_harmonic_cnn_v1", index=work.index)
    )
    work["morphology_model_version"] = work["compact_model_version"]
    return work


def verify_franklin_multisector_batch(
    queue: pd.DataFrame,
    overlap: pd.DataFrame,
    hidden: pd.DataFrame,
    *,
    sector_quotas: Mapping[int, EnrichmentQuotas],
    expected_overlap_count: int,
    excluded_tics: set[int] | None = None,
    sheet_dir: Path | None = None,
) -> dict[str, Any]:
    total = int(sum(value.total for value in sector_quotas.values()))
    if len(queue) != total or len(hidden) != total:
        raise ValueError(
            f"Franklin queue size mismatch: queue={len(queue)}, hidden={len(hidden)}, "
            f"expected={total}"
        )
    if queue["tic"].nunique() != total:
        raise ValueError("Franklin queue must contain one globally unique TIC per row")
    if queue["row_id"].tolist() != list(range(total)):
        raise ValueError("Franklin queue row_id is not normalized display order")
    observed_sector_counts = {
        int(key): int(value)
        for key, value in queue["sector"].value_counts().sort_index().items()
    }
    expected_sector_counts = {
        int(sector): int(quota.total) for sector, quota in sorted(sector_quotas.items())
    }
    if observed_sector_counts != expected_sector_counts:
        raise ValueError(
            f"sector counts differ: observed={observed_sector_counts}, "
            f"expected={expected_sector_counts}"
        )
    if not pd.to_numeric(queue["rep_peak_rank"], errors="coerce").eq(1).all():
        raise ValueError("Franklin queue contains a non-rank-1 review ephemeris")
    if set(queue["source_kind"].fillna("").astype(str)) != {"real_candidate"}:
        raise ValueError("Franklin queue is not real-only")
    if queue["twirl_vet_sheet_pdf_name"].fillna("").astype(str).ne("").any():
        raise ValueError("Franklin queue requests PDF vet sheets")
    expected_key = queue.apply(candidate_key, axis=1)
    if not queue["candidate_key"].astype(str).equals(expected_key.astype(str)):
        raise ValueError("Franklin queue candidate_key mismatch")

    hidden_tokens = (
        "p_",
        "std_p_",
        "member_",
        "selection_",
        "model_",
        "source_candidate",
        "active_learning_",
    )
    exposed = [column for column in queue if column.startswith(hidden_tokens)]
    if exposed:
        raise ValueError(f"public Franklin queue exposes provenance: {exposed}")
    if set(queue["candidate_key"]) != set(hidden["candidate_key"]):
        raise ValueError("hidden provenance does not match the public Franklin queue")

    expected_bucket_counts = asdict(_sum_quotas(list(sector_quotas.values())))
    expected_bucket_counts = {
        key: int(value) for key, value in expected_bucket_counts.items() if value > 0
    }
    observed_bucket_counts = {
        str(key): int(value)
        for key, value in hidden["selection_bucket"].value_counts().items()
    }
    if observed_bucket_counts != expected_bucket_counts:
        raise ValueError(
            f"selection bucket counts differ: observed={observed_bucket_counts}, "
            f"expected={expected_bucket_counts}"
        )

    expected_overlap = min(max(0, int(expected_overlap_count)), total)
    if len(overlap) != expected_overlap:
        raise ValueError(
            f"double-review count differs: observed={len(overlap)}, "
            f"expected={expected_overlap}"
        )
    if not set(overlap.get("candidate_key", pd.Series(dtype=str))).issubset(
        set(queue["candidate_key"])
    ):
        raise ValueError("double-review rows are absent from the Franklin queue")
    excluded_tics = excluded_tics or set()
    overlap_tics = set(pd.to_numeric(queue["tic"], errors="raise").astype(int)) & excluded_tics
    if overlap_tics:
        raise ValueError(
            f"Franklin queue overlaps {len(overlap_tics)} excluded TICs; "
            f"first={sorted(overlap_tics)[:5]}"
        )

    sheet_status = "not_checked"
    if sheet_dir is not None:
        root = Path(sheet_dir)
        missing = [
            str(name)
            for name in queue["twirl_vet_sheet_name"].fillna("").astype(str)
            if not name or not (root / name).is_file()
        ]
        if missing:
            raise FileNotFoundError(
                f"Franklin queue is missing {len(missing)} PNG sheets; first={missing[:3]}"
            )
        if list(root.glob("*.pdf")):
            raise ValueError("Franklin sheet directory contains PDFs")
        sheet_status = "passed"
    return {
        "passed": True,
        "n_rows": total,
        "n_unique_tics": int(queue["tic"].nunique()),
        "sector_counts": observed_sector_counts,
        "selection_bucket_counts": observed_bucket_counts,
        "n_double_review": int(len(overlap)),
        "sheet_status": sheet_status,
    }


def build_franklin_multisector_batch(
    score_tables: Mapping[int, pd.DataFrame],
    *,
    sector_quotas: Mapping[int, EnrichmentQuotas],
    excluded_tables: Sequence[pd.DataFrame] = (),
    batch_index: int = 0,
    double_review_count: int = 0,
    seed: int = 560717,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Build one blinded queue with unequal sector quotas and rank-1 folds."""

    sectors = tuple(sorted(int(value) for value in score_tables))
    if sectors != tuple(sorted(int(value) for value in sector_quotas)):
        raise ValueError("score-table sectors and quota sectors differ")
    if not sectors:
        raise ValueError("at least one sector score table is required")

    excluded_tics: set[int] = set()
    for frame in excluded_tables:
        if "tic" in frame:
            excluded_tics.update(
                pd.to_numeric(frame["tic"], errors="coerce")
                .dropna()
                .astype(np.int64)
                .tolist()
            )
    cumulative_exclusions = list(excluded_tables)
    selected: list[pd.DataFrame] = []
    sector_summaries: dict[str, Any] = {}
    for sector in sectors:
        prepared = prepare_single_teacher_scores(score_tables[sector])
        _, _, hidden, summary = build_existing_teacher_enrichment_batch(
            prepared,
            sector=sector,
            batch_index=int(batch_index),
            excluded_tables=cumulative_exclusions,
            quotas=sector_quotas[sector],
            seed=int(seed),
            required_peak_rank=1,
            double_review_count=0,
        )
        selected.append(hidden)
        cumulative_exclusions.append(hidden.loc[:, ["tic"]].copy())
        summary["policy_version"] = FRANKLIN_MULTISECTOR_POLICY_VERSION
        summary["compact_ranker"] = (
            "Teacher-v1 sqrt(p_planet_like * p_preserve)"
        )
        summary["morphology_ranker"] = (
            "Teacher-v1 shape_plus_periodogram_bls ensemble"
        )
        sector_summaries[str(sector)] = summary

    hidden = pd.concat(selected, ignore_index=True, sort=False)
    hidden = hidden.sample(
        frac=1.0,
        random_state=int(seed) + int(batch_index) + sum(sectors),
    ).reset_index(drop=True)
    hidden["row_id"] = np.arange(len(hidden), dtype=int)
    hidden["candidate_key"] = hidden.apply(candidate_key, axis=1)
    hidden["active_learning_policy_version"] = FRANKLIN_MULTISECTOR_POLICY_VERSION
    hidden["double_review"] = False
    n_overlap = min(max(0, int(double_review_count)), len(hidden))
    overlap_indices = hidden.sample(
        n=n_overlap,
        random_state=int(seed) + int(batch_index) + 10_000 + sum(sectors),
    ).index
    hidden.loc[overlap_indices, "double_review"] = True

    queue = hidden.loc[
        :, [column for column in PUBLIC_COLUMNS if column in hidden]
    ].copy()
    overlap = queue.loc[overlap_indices].reset_index(drop=True).copy()
    overlap["row_id"] = np.arange(len(overlap), dtype=int)
    verification = verify_franklin_multisector_batch(
        queue,
        overlap,
        hidden,
        sector_quotas=sector_quotas,
        expected_overlap_count=n_overlap,
        excluded_tics=excluded_tics,
    )
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "policy_version": FRANKLIN_MULTISECTOR_POLICY_VERSION,
        "batch_index": int(batch_index),
        "seed": int(seed),
        "rank_policy": "ADP-small representative BLS peak rank 1 only",
        "compact_score_policy": "sqrt(p_planet_like * p_preserve)",
        "n_initial_excluded_tics": int(len(excluded_tics)),
        "sector_quotas": {
            str(sector): asdict(sector_quotas[sector]) for sector in sectors
        },
        "sector_summaries": sector_summaries,
        "scores_hidden_from_browser": True,
        "training_performed": False,
        "verification": verification,
    }
    return queue, overlap, hidden, summary


def write_franklin_multisector_batch(
    *,
    sector_score_paths: Mapping[int, Path],
    sector_quotas: Mapping[int, EnrichmentQuotas],
    out_dir: Path,
    exclude_paths: Sequence[Path] = (),
    batch_index: int = 0,
    double_review_count: int = 0,
    seed: int = 560717,
) -> dict[str, Any]:
    scores = {
        int(sector): _read_table(Path(path))
        for sector, path in sector_score_paths.items()
    }
    exclusions = [
        _read_table(Path(path)) for path in exclude_paths if Path(path).is_file()
    ]
    queue, overlap, hidden, summary = build_franklin_multisector_batch(
        scores,
        sector_quotas=sector_quotas,
        excluded_tables=exclusions,
        batch_index=batch_index,
        double_review_count=double_review_count,
        seed=seed,
    )
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    count_token = (
        f"{len(queue) // 1000}k" if len(queue) % 1000 == 0 else str(len(queue))
    )
    queue_path = _write_table(
        queue, out_dir / f"franklin_review_queue_{count_token}_real.csv"
    )
    overlap_path = _write_table(
        overlap, out_dir / f"double_review_queue_{len(overlap)}.csv"
    )
    hidden_path = _write_table(
        hidden, out_dir / "hidden_selection_provenance.parquet"
    )
    sector_paths: dict[str, str] = {}
    for sector in sorted(sector_score_paths):
        sector_queue = queue.loc[queue["sector"].eq(int(sector))].copy()
        path = _write_table(
            sector_queue,
            out_dir
            / f"sector_{int(sector):04d}_review_queue_{len(sector_queue)}.csv",
        )
        sector_paths[str(sector)] = str(path)

    summary["score_provenance"] = {
        str(sector): {
            "path": str(Path(path)),
            "sha256": _sha256(Path(path)),
        }
        for sector, path in sorted(sector_score_paths.items())
    }
    summary["exclusion_provenance"] = [
        {"path": str(path), "sha256": _sha256(Path(path))}
        for path in exclude_paths
        if Path(path).is_file()
    ]
    summary["outputs"] = {
        "review_queue": str(queue_path),
        "double_review_queue": str(overlap_path),
        "hidden_selection_provenance": str(hidden_path),
        "sector_queues": sector_paths,
        "summary": str(out_dir / "summary.json"),
    }
    (out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    return summary


__all__ = [
    "FRANKLIN_LABEL_RETURN_POLICY_VERSION",
    "FRANKLIN_MULTISECTOR_POLICY_VERSION",
    "build_franklin_multisector_batch",
    "normalize_franklin_label_return",
    "prepare_single_teacher_scores",
    "standalone_app_candidate_key",
    "verify_franklin_multisector_batch",
    "write_franklin_multisector_batch",
]
