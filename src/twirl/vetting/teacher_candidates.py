"""Build A2v1 real-candidate rows for harmonic-teacher inference."""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

from twirl.io.compact_export import read_compact_lc_export
from twirl.vetting.adp_only import (
    ADP_ONLY_APERTURES,
    ADP_ONLY_CONTRACT_VERSION,
    assert_adp_only_search_frame,
    classify_period_relation,
)
from twirl.vetting.harmonic_inputs import native_group_path
from twirl.vetting.label_io import candidate_key
from twirl.vetting.teacher_active_learning import A2V1_TEACHER_INPUT_CONTRACT
from twirl.vetting.two_aperture import measure_two_aperture_candidate_metadata


PEAK_FIELDS: tuple[str, ...] = (
    "peak_rank",
    "period_d",
    "t0_bjd",
    "duration_min",
    "depth",
    "depth_snr",
    "sde",
    "log_power",
)


def normalize_a2v1_peak_candidates(
    peaks: pd.DataFrame,
    *,
    small_peaks_per_tic: int = 3,
    sector: int | None = None,
) -> pd.DataFrame:
    """Use ADP-small top-N peaks and attach the ADP-primary rank-1 context."""

    frame = peaks.copy()
    if sector is None:
        if "sector" not in frame:
            sector = 56
        else:
            sectors = (
                pd.to_numeric(frame["sector"], errors="coerce")
                .dropna()
                .astype(int)
                .unique()
            )
            if len(sectors) != 1:
                raise ValueError(
                    "teacher scoring requires exactly one sector; "
                    f"got {sorted(sectors.tolist())}"
                )
            sector = int(sectors[0])
    sector = int(sector)
    assert_adp_only_search_frame(frame)
    if "source_product_tag" not in frame:
        raise KeyError("A2v1 BLS peaks are missing source_product_tag")
    tags = set(frame["source_product_tag"].fillna("").astype(str))
    if tags != {"A2v1"}:
        raise ValueError(f"teacher scoring requires source_product_tag=A2v1; got {sorted(tags)}")
    rank = pd.to_numeric(frame["peak_rank"], errors="coerce")
    small = frame.loc[
        frame["aperture"].astype(str).eq(ADP_ONLY_APERTURES[0])
        & rank.between(1, int(small_peaks_per_tic))
        & frame["status"].astype(str).eq("ok")
    ].copy()
    if small.empty:
        raise ValueError("A2v1 peak table contains no valid ADP-small candidates")
    primary = frame.loc[
        frame["aperture"].astype(str).eq(ADP_ONLY_APERTURES[1])
        & rank.eq(1)
        & frame["status"].astype(str).eq("ok")
    ].copy()
    primary = primary.sort_values("tic", kind="stable").drop_duplicates("tic", keep="first")
    primary_fields = [column for column in PEAK_FIELDS if column in primary]
    primary = primary.loc[:, ["tic", *primary_fields]].rename(
        columns={column: f"adp_{column}" for column in primary_fields}
    )
    out = small.merge(primary, on="tic", how="left", validate="many_to_one")
    out = out.dropna(subset=["adp_period_d", "adp_t0_bjd", "adp_duration_min"])
    out["tic"] = pd.to_numeric(out["tic"], errors="coerce").astype(np.int64)
    out["sector"] = int(sector)
    out["rep_aperture"] = ADP_ONLY_APERTURES[0]
    out["rep_peak_rank"] = pd.to_numeric(out["peak_rank"], errors="coerce").astype(int)
    out["sde_max"] = pd.to_numeric(out["sde"], errors="coerce")
    out["anchor_aperture"] = ADP_ONLY_APERTURES[0]
    out["anchor_period_d"] = pd.to_numeric(out["period_d"], errors="coerce")
    out["anchor_t0_bjd"] = pd.to_numeric(out["t0_bjd"], errors="coerce")
    out["anchor_duration_min"] = pd.to_numeric(out["duration_min"], errors="coerce")
    out["anchor_sde"] = out["sde_max"]
    for column in PEAK_FIELDS:
        if column in out:
            out[f"adp_sml_{column}"] = pd.to_numeric(out[column], errors="coerce")
    relation = classify_period_relation(out["adp_sml_period_d"], out["adp_period_d"])
    out["aperture_period_relation"] = relation
    out["aperture_period_rel_delta"] = (
        np.abs(out["adp_period_d"] - out["adp_sml_period_d"]) / out["adp_sml_period_d"]
    )
    out["aperture_depth_ratio_primary_over_small"] = (
        pd.to_numeric(out.get("adp_depth"), errors="coerce")
        / pd.to_numeric(out.get("adp_sml_depth"), errors="coerce").replace(0, np.nan)
    )
    out["aperture_disagreement_flag"] = relation.eq("unrelated")
    out["n_apertures_agree"] = np.where(relation.isin({"exact", "harmonic"}), 2, 1)
    out["apertures_agree"] = np.where(
        relation.isin({"exact", "harmonic"}),
        ",".join(ADP_ONLY_APERTURES),
        ADP_ONLY_APERTURES[0],
    )
    out["review_id"] = [
        f"s{int(sector):04d}-A2v1-adp-small-{int(tic):016d}-r{int(peak_rank)}"
        for tic, peak_rank in zip(out["tic"], out["rep_peak_rank"])
    ]
    if out["review_id"].duplicated().any():
        raise ValueError("A2v1 candidate review_id values are not unique")
    out["source_kind"] = "real_candidate"
    out["source_product_tag"] = "A2v1"
    out["bls_search_branch"] = "A2v1_current_adp"
    out["adp_only_contract_version"] = ADP_ONLY_CONTRACT_VERSION
    out["input_contract_version"] = A2V1_TEACHER_INPUT_CONTRACT
    out["native_input_include"] = True
    out["native_group_path"] = [native_group_path(row) for row in out.to_dict("records")]
    out["candidate_key"] = out.apply(candidate_key, axis=1)
    return out.sort_values(["tic", "rep_peak_rank"], kind="stable").reset_index(drop=True)


def _peak_from_record(record: dict[str, Any], prefix: str) -> dict[str, Any]:
    return {
        key: record.get(f"{prefix}{key}", np.nan)
        for key in ("period_d", "t0_bjd", "duration_min", "depth", "depth_snr", "sde")
    }


def _measure_tic(payload: tuple[int, list[dict[str, Any]], str]) -> list[dict[str, Any]]:
    tic, records, compact_lc_path = payload
    lc = read_compact_lc_export(
        Path(compact_lc_path), tic=int(tic), columns=ADP_ONLY_APERTURES
    )
    if lc is None:
        return [
            {"review_id": str(record["review_id"]), "metadata_status": "missing_compact_lc"}
            for record in records
        ]
    output: list[dict[str, Any]] = []
    for record in records:
        anchor = _peak_from_record(record, "")
        own = {
            ADP_ONLY_APERTURES[0]: _peak_from_record(record, "adp_sml_"),
            ADP_ONLY_APERTURES[1]: _peak_from_record(record, "adp_"),
        }
        try:
            measured = measure_two_aperture_candidate_metadata(
                lc,
                anchor_peak=anchor,
                own_peaks=own,
                apertures=ADP_ONLY_APERTURES,
            )
            measured["review_id"] = str(record["review_id"])
            measured["metadata_status"] = "ok"
        except Exception as exc:
            measured = {
                "review_id": str(record["review_id"]),
                "metadata_status": "error",
                "metadata_error": str(exc),
            }
        output.append(measured)
    return output


def enrich_candidate_metadata(
    candidates: pd.DataFrame,
    *,
    compact_lc_path: Path,
    workers: int = 1,
    progress_every: int = 500,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Measure training-compatible odd/even, trend, and aperture metadata."""

    work = candidates.copy()
    payloads = [
        (int(tic), group.to_dict("records"), str(compact_lc_path))
        for tic, group in work.groupby("tic", sort=True)
    ]
    rows: list[dict[str, Any]] = []
    if int(workers) <= 1:
        iterator = map(_measure_tic, payloads)
        executor = None
    else:
        executor = ProcessPoolExecutor(max_workers=int(workers))
        iterator = executor.map(_measure_tic, payloads, chunksize=1)
    try:
        for index, batch in enumerate(iterator, start=1):
            rows.extend(batch)
            if progress_every > 0 and index % int(progress_every) == 0:
                print(
                    f"[teacher-candidates] metadata {index:,}/{len(payloads):,} TICs",
                    flush=True,
                )
    finally:
        if executor is not None:
            executor.shutdown(wait=True)
    metrics = pd.DataFrame(rows)
    if metrics["review_id"].duplicated().any():
        raise RuntimeError("candidate metadata contains duplicate review_id values")
    overlap = [
        column for column in metrics.columns if column != "review_id" and column in work.columns
    ]
    work = work.drop(columns=overlap).merge(
        metrics, on="review_id", how="left", validate="one_to_one"
    )
    status = work["metadata_status"].fillna("missing").astype(str)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "n_candidates": int(len(work)),
        "n_unique_tics": int(work["tic"].nunique()),
        "compact_lc_path": str(compact_lc_path),
        "workers": int(workers),
        "metadata_status_counts": {
            str(key): int(value) for key, value in status.value_counts().sort_index().items()
        },
        "passed": bool(status.eq("ok").all()),
    }
    return work, summary


__all__ = [
    "enrich_candidate_metadata",
    "normalize_a2v1_peak_candidates",
]
