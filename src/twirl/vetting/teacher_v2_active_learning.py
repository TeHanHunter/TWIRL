"""Blinded active-learning queues from existing Teacher-v2 rankers."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from twirl.vetting.adjudication import ADJUDICATION_VET_SHEET_VERSION
from twirl.vetting.label_io import candidate_key


ACTIVE_LEARNING_POLICY_VERSION = "s56_s64_existing_teacher_multi_ranker_v1"
MIXED_ACTIVE_LEARNING_POLICY_VERSION = "s56_s57_existing_teacher_balanced_v1"


@dataclass(frozen=True)
class EnrichmentQuotas:
    compact_transit: int = 400
    eclipse_contact: int = 300
    smooth_variable: int = 100
    model_disagreement: int = 100
    stratified_control: int = 100

    @property
    def total(self) -> int:
        return int(sum(asdict(self).values()))


DEFAULT_ENRICHMENT_QUOTAS = EnrichmentQuotas()

IDENTITY_COLUMNS: tuple[str, ...] = (
    "review_id",
    "tic",
    "sector",
    "period_d",
    "t0_bjd",
    "duration_min",
)

MORPHOLOGY_SCORE_COLUMNS: tuple[str, ...] = (
    "p_planet_like",
    "std_p_planet_like",
    "p_eclipse_contact",
    "std_p_eclipse_contact",
    "p_smooth_variable",
    "std_p_smooth_variable",
    "p_other",
    "std_p_other",
    "p_preserve",
    "std_p_preserve",
    "p_compact_transit",
    "std_p_compact_transit",
)

PUBLIC_COLUMNS: tuple[str, ...] = (
    "row_id",
    "review_id",
    "tic",
    "sector",
    "cam",
    "ccd",
    "tmag",
    "source_kind",
    "source_bucket",
    "period_d",
    "t0_bjd",
    "duration_min",
    "depth",
    "depth_snr",
    "sde_max",
    "rep_aperture",
    "rep_peak_rank",
    "n_apertures_agree",
    "apertures_agree",
    "aperture_period_rel_delta",
    "aperture_disagreement_flag",
    "candidate_key",
    "vet_sheet_version",
    "twirl_vet_sheet_name",
    "twirl_vet_sheet_pdf_name",
)


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
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _require_columns(frame: pd.DataFrame, columns: Sequence[str], *, name: str) -> None:
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise KeyError(f"{name} is missing columns: {missing}")


def _require_unique_review_ids(frame: pd.DataFrame, *, name: str) -> None:
    review_id = frame["review_id"].fillna("").astype(str)
    if review_id.eq("").any() or review_id.duplicated().any():
        raise ValueError(f"{name} must contain nonempty unique review_id values")


def combine_existing_ranker_scores(
    compact_scores: pd.DataFrame,
    morphology_scores: pd.DataFrame,
) -> pd.DataFrame:
    """Combine chronology compact scores with Shape+BLS morphology scores."""

    _require_columns(
        compact_scores,
        (*IDENTITY_COLUMNS, "p_compact_transit", "std_p_compact_transit"),
        name="compact score table",
    )
    _require_columns(
        morphology_scores,
        (*IDENTITY_COLUMNS, *MORPHOLOGY_SCORE_COLUMNS),
        name="morphology score table",
    )
    _require_unique_review_ids(compact_scores, name="compact score table")
    _require_unique_review_ids(morphology_scores, name="morphology score table")

    right_columns = list(IDENTITY_COLUMNS) + list(MORPHOLOGY_SCORE_COLUMNS)
    for column in ("model_profile", "model_version"):
        if column in morphology_scores:
            right_columns.append(column)
    right = morphology_scores.loc[:, right_columns].copy()
    right = right.rename(
        columns={
            **{
                column: f"{column}_morphology"
                for column in IDENTITY_COLUMNS
                if column != "review_id"
            },
            "p_compact_transit": "p_compact_transit_morphology_profile",
            "std_p_compact_transit": "std_p_compact_transit_morphology_profile",
            "model_profile": "morphology_model_profile",
            "model_version": "morphology_model_version",
        }
    )
    compact_only_scores = {"p_compact_transit", "std_p_compact_transit"}
    left = compact_scores.drop(
        columns=[
            column
            for column in MORPHOLOGY_SCORE_COLUMNS
            if column not in compact_only_scores and column in compact_scores
        ]
    ).copy()
    left = left.rename(
        columns={
            "model_profile": "compact_model_profile",
            "model_version": "compact_model_version",
        }
    )
    combined = left.merge(right, on="review_id", how="inner", validate="one_to_one")
    if len(combined) != len(left) or len(combined) != len(right):
        raise ValueError("ranker score tables do not contain the same review_id set")

    for column in IDENTITY_COLUMNS:
        if column == "review_id":
            continue
        other = f"{column}_morphology"
        if column in {"tic", "sector"}:
            matches = pd.to_numeric(combined[column], errors="coerce").eq(
                pd.to_numeric(combined[other], errors="coerce")
            )
        else:
            left_value = pd.to_numeric(combined[column], errors="coerce").to_numpy(float)
            right_value = pd.to_numeric(combined[other], errors="coerce").to_numpy(float)
            matches = pd.Series(
                np.isclose(left_value, right_value, rtol=0.0, atol=1.0e-10),
                index=combined.index,
            )
        if not bool(matches.all()):
            raise ValueError(f"ranker score tables disagree on {column}")
        combined = combined.drop(columns=other)

    if "source_kind" in combined:
        injected = combined["source_kind"].fillna("").astype(str).str.contains(
            "inject", case=False
        )
        if injected.any():
            raise ValueError("real enrichment inputs contain injected rows")
    if "is_injected_row" in combined:
        injected = combined["is_injected_row"].fillna(False).astype(bool)
        if injected.any():
            raise ValueError("real enrichment inputs contain injected rows")

    probability_columns = (
        "p_planet_like",
        "p_eclipse_contact",
        "p_smooth_variable",
        "p_other",
    )
    probabilities = combined.loc[:, probability_columns].apply(
        pd.to_numeric, errors="coerce"
    ).fillna(0.0).clip(0.0, 1.0)
    entropy = -(probabilities * np.log(probabilities.clip(lower=1.0e-12))).sum(axis=1)
    combined["morphology_entropy"] = entropy / np.log(len(probability_columns))
    sorted_probability = np.sort(probabilities.to_numpy(float), axis=1)
    combined["morphology_margin"] = sorted_probability[:, -1] - sorted_probability[:, -2]
    combined["compact_profile_disagreement"] = np.abs(
        pd.to_numeric(combined["p_compact_transit"], errors="coerce")
        - pd.to_numeric(
            combined["p_compact_transit_morphology_profile"], errors="coerce"
        )
    )
    return combined


def _excluded_tics(tables: Sequence[pd.DataFrame]) -> set[int]:
    excluded: set[int] = set()
    for frame in tables:
        if frame.empty or "tic" not in frame:
            continue
        values = pd.to_numeric(frame["tic"], errors="coerce").dropna().astype(np.int64)
        excluded.update(int(value) for value in values)
    return excluded


def _top_unique_tics(
    candidates: pd.DataFrame,
    *,
    score_column: str,
    count: int,
    used_tics: set[int],
) -> pd.DataFrame:
    available = candidates.loc[~candidates["tic"].isin(used_tics)].copy()
    available = available.sort_values(
        [score_column, "sde_max", "rep_peak_rank", "tic"],
        ascending=[False, False, True, True],
        kind="stable",
        na_position="last",
    ).drop_duplicates("tic", keep="first")
    if len(available) < count:
        raise ValueError(
            f"bucket {score_column} requested {count} unique TICs but found {len(available)}"
        )
    selected = available.head(int(count)).copy()
    selected["selection_rank"] = np.arange(1, len(selected) + 1, dtype=int)
    return selected


def _stable_hash(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _stratified_controls(
    candidates: pd.DataFrame,
    *,
    count: int,
    used_tics: set[int],
    seed: int,
) -> pd.DataFrame:
    available = candidates.loc[~candidates["tic"].isin(used_tics)].copy()
    available = available.sort_values(
        ["rep_peak_rank", "sde_max", "tic"],
        ascending=[True, False, True],
        kind="stable",
        na_position="last",
    ).drop_duplicates("tic", keep="first")
    if len(available) < count:
        raise ValueError(f"control bucket requested {count} TICs but found {len(available)}")
    bin_columns: list[str] = []
    for source in ("tmag", "period_d", "sde_max"):
        output = f"_control_{source}_bin"
        values = pd.to_numeric(available.get(source), errors="coerce")
        try:
            available[output] = pd.qcut(
                values.rank(method="first"), q=4, labels=False, duplicates="drop"
            ).fillna(-1).astype(int)
        except ValueError:
            available[output] = -1
        bin_columns.append(output)
    available["_control_stratum"] = available[bin_columns].astype(str).agg("|".join, axis=1)
    available["_control_hash"] = [
        _stable_hash(f"{seed}|{tic}") for tic in available["tic"].astype(np.int64)
    ]
    available = available.sort_values(
        ["_control_stratum", "_control_hash"], kind="stable"
    )
    available["_control_round"] = available.groupby("_control_stratum").cumcount()
    selected = available.sort_values(
        ["_control_round", "_control_hash"], kind="stable"
    ).head(int(count)).copy()
    selected["selection_rank"] = np.arange(1, len(selected) + 1, dtype=int)
    return selected.drop(columns=[*bin_columns, "_control_stratum", "_control_hash", "_control_round"])


def build_existing_teacher_enrichment_batch(
    scores: pd.DataFrame,
    *,
    sector: int,
    batch_index: int,
    excluded_tables: Sequence[pd.DataFrame] = (),
    quotas: EnrichmentQuotas = DEFAULT_ENRICHMENT_QUOTAS,
    seed: int = 56,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Select one deterministic, blinded real-TIC enrichment batch."""

    if quotas.total <= 0:
        raise ValueError("enrichment quotas must contain at least one row")
    _require_columns(
        scores,
        (*IDENTITY_COLUMNS, *MORPHOLOGY_SCORE_COLUMNS, "tmag", "sde_max", "rep_peak_rank"),
        name="scores",
    )
    candidates = scores.copy()
    candidates["tic"] = pd.to_numeric(candidates["tic"], errors="coerce").astype("Int64")
    candidates["sector"] = pd.to_numeric(candidates["sector"], errors="coerce").astype("Int64")
    candidates = candidates.dropna(subset=["tic", "sector", "period_d", "t0_bjd", "duration_min"])
    candidates["tic"] = candidates["tic"].astype(np.int64)
    candidates["sector"] = candidates["sector"].astype(np.int16)
    candidates = candidates.loc[candidates["sector"].eq(int(sector))].copy()
    if candidates.empty:
        raise ValueError(f"score table contains no Sector {sector} rows")
    excluded = _excluded_tics(excluded_tables)
    candidates = candidates.loc[~candidates["tic"].isin(excluded)].copy()
    if candidates["tic"].nunique() < quotas.total:
        raise ValueError(
            f"Sector {sector} has {candidates['tic'].nunique()} eligible TICs; "
            f"requested={quotas.total}"
        )

    for column in (
        "p_compact_transit",
        "p_eclipse_contact",
        "p_smooth_variable",
        "std_p_compact_transit",
        "std_p_eclipse_contact",
        "std_p_smooth_variable",
        "std_p_other",
    ):
        candidates[column] = pd.to_numeric(candidates.get(column), errors="coerce").fillna(0.0)
    candidates["_score_compact_transit"] = candidates["p_compact_transit"]
    candidates["_score_eclipse_contact"] = candidates["p_eclipse_contact"]
    candidates["_score_smooth_variable"] = candidates["p_smooth_variable"]
    uncertainty = candidates[
        [
            "std_p_compact_transit",
            "std_p_eclipse_contact",
            "std_p_smooth_variable",
            "std_p_other",
        ]
    ].max(axis=1)
    candidates["_score_model_disagreement"] = (
        0.40 * pd.to_numeric(candidates["morphology_entropy"], errors="coerce").fillna(0.0)
        + 0.35
        * pd.to_numeric(candidates["compact_profile_disagreement"], errors="coerce").fillna(0.0)
        + 0.25 * uncertainty
    )

    used_tics: set[int] = set()
    selections: list[pd.DataFrame] = []
    # Remove likely EBs and variables before selecting the compact-transit bucket.
    bucket_specs = (
        ("eclipse_contact", quotas.eclipse_contact, "_score_eclipse_contact"),
        ("smooth_variable", quotas.smooth_variable, "_score_smooth_variable"),
        ("compact_transit", quotas.compact_transit, "_score_compact_transit"),
        ("model_disagreement", quotas.model_disagreement, "_score_model_disagreement"),
    )
    for bucket, count, score_column in bucket_specs:
        if count <= 0:
            continue
        pool_size = int(candidates.loc[~candidates["tic"].isin(used_tics), "tic"].nunique())
        selected = _top_unique_tics(
            candidates,
            score_column=score_column,
            count=int(count),
            used_tics=used_tics,
        )
        selected["selection_bucket"] = bucket
        selected["selection_score"] = selected[score_column]
        selected["selection_pool_size"] = pool_size
        selected["selection_fraction"] = float(count) / max(pool_size, 1)
        used_tics.update(int(value) for value in selected["tic"])
        selections.append(selected)
    if quotas.stratified_control > 0:
        pool_size = int(candidates.loc[~candidates["tic"].isin(used_tics), "tic"].nunique())
        controls = _stratified_controls(
            candidates,
            count=int(quotas.stratified_control),
            used_tics=used_tics,
            seed=int(seed) + int(sector) + 1000 * int(batch_index),
        )
        controls["selection_bucket"] = "stratified_control"
        controls["selection_score"] = np.nan
        controls["selection_pool_size"] = pool_size
        controls["selection_fraction"] = float(quotas.stratified_control) / max(pool_size, 1)
        used_tics.update(int(value) for value in controls["tic"])
        selections.append(controls)

    hidden = pd.concat(selections, ignore_index=True, sort=False)
    if len(hidden) != quotas.total or hidden["tic"].nunique() != quotas.total:
        raise RuntimeError("enrichment selection did not produce one unique TIC per requested row")
    hidden = hidden.sample(
        frac=1.0,
        random_state=int(seed) + int(sector) + int(batch_index),
    ).reset_index(drop=True)
    hidden["source_candidate_review_id"] = hidden["review_id"].astype(str)
    hidden["review_id"] = [
        f"s{int(sector):04d}-a2v1-enrich-b{int(batch_index):02d}-{index:04d}"
        for index in range(len(hidden))
    ]
    hidden["row_id"] = np.arange(len(hidden), dtype=int)
    hidden["source_kind"] = "real_candidate"
    hidden["source_bucket"] = ""
    hidden["candidate_key"] = hidden.apply(candidate_key, axis=1)
    hidden["vet_sheet_version"] = ADJUDICATION_VET_SHEET_VERSION
    hidden["twirl_vet_sheet_name"] = (
        hidden["review_id"] + "_twirl_twoap_current_adp.png"
    )
    hidden["twirl_vet_sheet_pdf_name"] = ""
    hidden["active_learning_policy_version"] = ACTIVE_LEARNING_POLICY_VERSION
    hidden["batch_index"] = int(batch_index)
    hidden["double_review"] = False
    overlap_indices = hidden.sample(
        n=min(100, len(hidden)),
        random_state=int(seed) + int(sector) + 10_000 + int(batch_index),
    ).index
    hidden.loc[overlap_indices, "double_review"] = True

    queue = hidden.loc[:, [column for column in PUBLIC_COLUMNS if column in hidden]].copy()
    overlap = queue.loc[overlap_indices].reset_index(drop=True).copy()
    overlap["row_id"] = np.arange(len(overlap), dtype=int)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "policy_version": ACTIVE_LEARNING_POLICY_VERSION,
        "vet_sheet_version": ADJUDICATION_VET_SHEET_VERSION,
        "sector": int(sector),
        "batch_index": int(batch_index),
        "n_rows": int(len(queue)),
        "n_unique_tics": int(queue["tic"].nunique()),
        "n_double_review": int(len(overlap)),
        "n_excluded_tics": int(len(excluded)),
        "quotas": asdict(quotas),
        "selection_bucket_counts": {
            str(key): int(value)
            for key, value in hidden["selection_bucket"].value_counts().sort_index().items()
        },
        "training_performed": False,
        "compact_ranker": "existing shape_plus_raw_chronology ensemble",
        "morphology_ranker": "existing shape_plus_periodogram_bls ensemble",
        "scores_hidden_from_browser": True,
    }
    verify_enrichment_batch(queue, overlap, hidden, expected_quotas=quotas)
    return queue, overlap, hidden, summary


def verify_enrichment_batch(
    queue: pd.DataFrame,
    overlap: pd.DataFrame,
    hidden: pd.DataFrame,
    *,
    expected_quotas: EnrichmentQuotas = DEFAULT_ENRICHMENT_QUOTAS,
    sheet_dir: Path | None = None,
) -> dict[str, Any]:
    if len(queue) != expected_quotas.total or queue["tic"].nunique() != len(queue):
        raise ValueError("public queue must contain exactly one unique TIC per requested row")
    if queue["row_id"].tolist() != list(range(len(queue))):
        raise ValueError("public queue row_id is not normalized display order")
    expected_key = queue.apply(candidate_key, axis=1)
    if not queue["candidate_key"].astype(str).equals(expected_key.astype(str)):
        raise ValueError("public queue candidate_key mismatch")
    hidden_tokens = ("p_", "std_p_", "member_", "selection_", "model_", "source_candidate")
    exposed = [column for column in queue if column.startswith(hidden_tokens)]
    if exposed:
        raise ValueError(f"public queue exposes hidden score/provenance columns: {exposed}")
    if len(hidden) != len(queue) or set(hidden["candidate_key"]) != set(queue["candidate_key"]):
        raise ValueError("hidden provenance does not match the public queue")
    counts = hidden["selection_bucket"].value_counts().to_dict()
    expected = {key: value for key, value in asdict(expected_quotas).items() if value > 0}
    if counts != expected:
        raise ValueError(f"selection bucket counts differ: observed={counts}, expected={expected}")
    if len(overlap) != min(100, len(queue)):
        raise ValueError("double-review queue must contain 100 rows or the complete short queue")
    if not set(overlap["candidate_key"]).issubset(set(queue["candidate_key"])):
        raise ValueError("double-review queue contains rows absent from the public queue")
    sheet_status = "not_checked"
    if sheet_dir is not None:
        sheet_dir = Path(sheet_dir)
        missing = [name for name in queue["twirl_vet_sheet_name"] if not (sheet_dir / name).is_file()]
        if missing:
            raise ValueError(f"missing {len(missing)} vet sheets; first={missing[:3]}")
        if list(sheet_dir.glob("*.pdf")):
            raise ValueError("enrichment sheet directory contains PDFs")
        sheet_status = "passed"
    return {
        "passed": True,
        "n_rows": int(len(queue)),
        "n_unique_tics": int(queue["tic"].nunique()),
        "n_overlap": int(len(overlap)),
        "sheet_status": sheet_status,
    }


def write_existing_teacher_enrichment_batch(
    *,
    compact_scores_path: Path,
    morphology_scores_path: Path,
    out_dir: Path,
    sector: int,
    batch_index: int,
    exclude_paths: Sequence[Path] = (),
    quotas: EnrichmentQuotas = DEFAULT_ENRICHMENT_QUOTAS,
) -> dict[str, Any]:
    compact = _read_table(compact_scores_path)
    morphology = _read_table(morphology_scores_path)
    combined = combine_existing_ranker_scores(compact, morphology)
    exclusions = [_read_table(path) for path in exclude_paths if Path(path).is_file()]
    queue, overlap, hidden, summary = build_existing_teacher_enrichment_batch(
        combined,
        sector=int(sector),
        batch_index=int(batch_index),
        excluded_tables=exclusions,
        quotas=quotas,
    )
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    queue_path = _write_table(queue, out_dir / "review_queue_1k.csv")
    overlap_path = _write_table(overlap, out_dir / "double_review_queue_100.csv")
    hidden_path = _write_table(hidden, out_dir / "hidden_selection_provenance.parquet")
    verify_enrichment_batch(
        pd.read_csv(queue_path),
        pd.read_csv(overlap_path),
        pd.read_parquet(hidden_path),
        expected_quotas=quotas,
    )
    summary.update(
        {
            "compact_scores_path": str(compact_scores_path),
            "compact_scores_sha256": _sha256(compact_scores_path),
            "morphology_scores_path": str(morphology_scores_path),
            "morphology_scores_sha256": _sha256(morphology_scores_path),
            "exclude_paths": [str(path) for path in exclude_paths],
            "outputs": {
                "review_queue": str(queue_path),
                "double_review_queue": str(overlap_path),
                "hidden_selection_provenance": str(hidden_path),
                "summary": str(out_dir / "summary.json"),
            },
        }
    )
    (out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    return summary


def _scaled_quotas(quotas: EnrichmentQuotas, factor: int) -> EnrichmentQuotas:
    values = {key: int(value) * int(factor) for key, value in asdict(quotas).items()}
    return EnrichmentQuotas(**values)


def build_mixed_existing_teacher_enrichment_batch(
    score_tables: Mapping[int, pd.DataFrame],
    *,
    batch_index: int,
    excluded_tables: Sequence[pd.DataFrame] = (),
    per_sector_quotas: EnrichmentQuotas,
    double_review_count: int = 100,
    seed: int = 56,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Build a sector-balanced queue with globally unique host TICs."""

    sectors = tuple(sorted(int(sector) for sector in score_tables))
    if len(sectors) < 2:
        raise ValueError("mixed enrichment requires at least two sectors")
    if per_sector_quotas.total <= 0:
        raise ValueError("per-sector quotas must contain at least one row")

    selected_hidden: list[pd.DataFrame] = []
    cumulative_exclusions = list(excluded_tables)
    per_sector_summaries: dict[str, Any] = {}
    for sector in sectors:
        _, _, hidden, sector_summary = build_existing_teacher_enrichment_batch(
            score_tables[sector],
            sector=sector,
            batch_index=batch_index,
            excluded_tables=cumulative_exclusions,
            quotas=per_sector_quotas,
            seed=seed,
        )
        selected_hidden.append(hidden)
        cumulative_exclusions.append(hidden.loc[:, ["tic"]].copy())
        per_sector_summaries[str(sector)] = sector_summary

    hidden = pd.concat(selected_hidden, ignore_index=True, sort=False)
    hidden = hidden.sample(
        frac=1.0,
        random_state=int(seed) + int(batch_index) + sum(sectors),
    ).reset_index(drop=True)
    hidden["row_id"] = np.arange(len(hidden), dtype=int)
    hidden["candidate_key"] = hidden.apply(candidate_key, axis=1)
    hidden["active_learning_policy_version"] = MIXED_ACTIVE_LEARNING_POLICY_VERSION
    hidden["double_review"] = False

    n_overlap = min(max(0, int(double_review_count)), len(hidden))
    overlap_indices = hidden.sample(
        n=n_overlap,
        random_state=int(seed) + int(batch_index) + 10_000 + sum(sectors),
    ).index
    hidden.loc[overlap_indices, "double_review"] = True
    queue = hidden.loc[:, [column for column in PUBLIC_COLUMNS if column in hidden]].copy()
    overlap = queue.loc[overlap_indices].reset_index(drop=True).copy()
    overlap["row_id"] = np.arange(len(overlap), dtype=int)

    global_quotas = _scaled_quotas(per_sector_quotas, len(sectors))
    verify_enrichment_batch(
        queue,
        overlap,
        hidden,
        expected_quotas=global_quotas,
    )
    sector_counts = {
        str(key): int(value)
        for key, value in queue["sector"].value_counts().sort_index().items()
    }
    expected_sector_count = int(per_sector_quotas.total)
    if sector_counts != {str(sector): expected_sector_count for sector in sectors}:
        raise ValueError(
            f"mixed queue sector balance differs: observed={sector_counts}, "
            f"expected_per_sector={expected_sector_count}"
        )
    if queue["tic"].nunique() != len(queue):
        raise ValueError("mixed queue contains a TIC repeated across sectors")

    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "policy_version": MIXED_ACTIVE_LEARNING_POLICY_VERSION,
        "vet_sheet_version": ADJUDICATION_VET_SHEET_VERSION,
        "sectors": list(sectors),
        "batch_index": int(batch_index),
        "n_rows": int(len(queue)),
        "n_unique_tics": int(queue["tic"].nunique()),
        "n_double_review": int(len(overlap)),
        "n_initial_excluded_tics": int(len(_excluded_tics(excluded_tables))),
        "per_sector_quotas": asdict(per_sector_quotas),
        "global_quotas": asdict(global_quotas),
        "sector_counts": sector_counts,
        "selection_bucket_counts": {
            str(key): int(value)
            for key, value in hidden["selection_bucket"].value_counts().sort_index().items()
        },
        "per_sector_summaries": per_sector_summaries,
        "cross_sector_tic_deduplication": True,
        "training_performed": False,
        "compact_ranker": "existing shape_plus_raw_chronology ensemble",
        "morphology_ranker": "existing shape_plus_periodogram_bls ensemble",
        "scores_hidden_from_browser": True,
    }
    return queue, overlap, hidden, summary


def write_mixed_existing_teacher_enrichment_batch(
    *,
    sector_score_paths: Mapping[int, tuple[Path, Path]],
    out_dir: Path,
    batch_index: int,
    exclude_paths: Sequence[Path] = (),
    per_sector_quotas: EnrichmentQuotas,
    double_review_count: int = 100,
) -> dict[str, Any]:
    """Read frozen score products and write one blinded mixed-sector batch."""

    scores: dict[int, pd.DataFrame] = {}
    score_provenance: dict[str, Any] = {}
    for sector, (compact_path, morphology_path) in sorted(sector_score_paths.items()):
        compact_path = Path(compact_path)
        morphology_path = Path(morphology_path)
        scores[int(sector)] = combine_existing_ranker_scores(
            _read_table(compact_path),
            _read_table(morphology_path),
        )
        score_provenance[str(int(sector))] = {
            "compact_scores_path": str(compact_path),
            "compact_scores_sha256": _sha256(compact_path),
            "morphology_scores_path": str(morphology_path),
            "morphology_scores_sha256": _sha256(morphology_path),
        }
    exclusions = [_read_table(path) for path in exclude_paths if Path(path).is_file()]
    queue, overlap, hidden, summary = build_mixed_existing_teacher_enrichment_batch(
        scores,
        batch_index=int(batch_index),
        excluded_tables=exclusions,
        per_sector_quotas=per_sector_quotas,
        double_review_count=int(double_review_count),
    )

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    queue_path = _write_table(queue, out_dir / "review_queue_1k.csv")
    overlap_path = _write_table(overlap, out_dir / "double_review_queue_100.csv")
    hidden_path = _write_table(hidden, out_dir / "hidden_selection_provenance.parquet")
    sector_queue_paths: dict[str, str] = {}
    for sector in sorted(scores):
        sector_queue = queue.loc[queue["sector"].eq(int(sector))].copy()
        path = _write_table(
            sector_queue,
            out_dir / f"sector_{int(sector):04d}_review_queue_{len(sector_queue)}.csv",
        )
        sector_queue_paths[str(int(sector))] = str(path)

    global_quotas = _scaled_quotas(per_sector_quotas, len(scores))
    verify_enrichment_batch(
        pd.read_csv(queue_path),
        pd.read_csv(overlap_path),
        pd.read_parquet(hidden_path),
        expected_quotas=global_quotas,
    )
    summary.update(
        {
            "score_provenance": score_provenance,
            "exclude_paths": [str(path) for path in exclude_paths],
            "outputs": {
                "review_queue": str(queue_path),
                "double_review_queue": str(overlap_path),
                "hidden_selection_provenance": str(hidden_path),
                "sector_queues": sector_queue_paths,
                "summary": str(out_dir / "summary.json"),
            },
        }
    )
    (out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    return summary


__all__ = [
    "ACTIVE_LEARNING_POLICY_VERSION",
    "DEFAULT_ENRICHMENT_QUOTAS",
    "EnrichmentQuotas",
    "build_existing_teacher_enrichment_batch",
    "build_mixed_existing_teacher_enrichment_batch",
    "combine_existing_ranker_scores",
    "verify_enrichment_batch",
    "write_existing_teacher_enrichment_batch",
    "write_mixed_existing_teacher_enrichment_batch",
]
