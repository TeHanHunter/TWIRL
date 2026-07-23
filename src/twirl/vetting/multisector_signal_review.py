"""Final Planet-like/EB morphology review across multiple TESS sectors.

This module keeps source labels immutable, collapses only exact candidate
duplicates under an explicit precedence rule, and prepares a pre-filled local
review queue.  Cross-sector observations of the same TIC remain separate.
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
from typing import Any, Iterable, Mapping

import numpy as np
import pandas as pd

from twirl.vetting.adjudication_audit import (
    HARMONIC_CNN_TARGET_POLICY,
    MORPHOLOGY_TARGET_MAP_V1,
    PRESERVE_LABELS_V1,
    REJECT_LABELS_V1,
    TRANSIT_HARMONIC_LABELS_V1,
)
from twirl.vetting.label_io import candidate_key, canonicalize_candidate_key


SIGNAL_REREVIEW_POLICY_VERSION = "s56_s62_planet_eb_rereview_v1"
SIGNAL_PRECEDENCE_POLICY_VERSION = "exact_candidate_later_review_wins_v1"
SIGNAL_REREVIEW_LABELS: frozenset[str] = frozenset(
    {"planet_like", "eclipsing_binary_or_pceb"}
)
FINAL_LABELS: frozenset[str] = frozenset(
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

_COMMON_REQUIRED: frozenset[str] = frozenset(
    {"review_id", "tic", "sector", "period_d", "t0_bjd", "duration_min"}
)


def _text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except TypeError:
        pass
    return str(value).strip()


def _truth(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return _text(value).lower() in {"1", "1.0", "true", "t", "yes", "y"}


def _require(frame: pd.DataFrame, columns: Iterable[str], *, name: str) -> None:
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise KeyError(f"{name} is missing columns: {missing}")


def _factor_text(value: Any) -> str:
    text = _text(value).lower()
    if not text or text in {"nan", "none"}:
        return "1"
    if text in {"unresolved", "?"}:
        return "unresolved"
    try:
        numeric = float(text)
    except ValueError:
        return "1"
    factors = (
        (0.25, "0.25"),
        (1.0 / 3.0, "0.3333333333333333"),
        (0.5, "0.5"),
        (1.0, "1"),
        (2.0, "2"),
        (3.0, "3"),
        (4.0, "4"),
    )
    for expected, canonical in factors:
        if np.isclose(numeric, expected, rtol=0.0, atol=1.0e-8):
            return canonical
    return "1"


def _canonical_observation_key(row: Mapping[str, Any]) -> str:
    payload = pd.Series(
        {
            "tic": row.get("tic", ""),
            "sector": row.get("sector", ""),
            "period_d": row.get("period_d", ""),
            "t0_bjd": row.get("t0_bjd", ""),
            "source_bucket": "",
        }
    )
    return candidate_key(payload)


def _source_record(
    row: Mapping[str, Any],
    *,
    source_batch_id: str,
    source_uid: str,
    source_label: str,
    source_labeler: str,
    source_notes: str,
    source_updated_utc: str,
    source_period_factor: Any,
    source_period_status: str,
    source_priority: int,
    original_harmonic_supervision_verified: bool,
    original_harmonic_target: str,
    original_harmonic_include: bool,
    original_factor_review_status: str,
) -> dict[str, Any]:
    if source_label not in SIGNAL_REREVIEW_LABELS:
        raise ValueError(f"non-signal label entered signal re-review: {source_label!r}")
    record = {column: row.get(column, "") for column in _COMMON_REQUIRED}
    for column in (
        "cam",
        "ccd",
        "tmag",
        "depth",
        "depth_snr",
        "sde_max",
        "rep_aperture",
        "rep_peak_rank",
        "n_apertures_agree",
        "apertures_agree",
        "aperture_period_rel_delta",
        "aperture_disagreement_flag",
        "twirl_vet_sheet_name",
        "twirl_vet_sheet_pdf_name",
        "adp_only_review_ok",
        "adp_only_exclusion_reason",
    ):
        record[column] = row.get(column, "")
    record.update(
        {
            "source_batch_id": source_batch_id,
            "source_uid": source_uid,
            "source_review_id": _text(row.get("review_id")),
            "source_candidate_key": _text(
                row.get("pipeline_candidate_key", row.get("candidate_key", ""))
            )
            or _canonical_observation_key(row),
            "source_label": source_label,
            "source_labeler": source_labeler,
            "source_notes": source_notes,
            "source_updated_utc": source_updated_utc,
            "source_period_factor": _factor_text(source_period_factor),
            "source_period_status": source_period_status,
            "source_priority": int(source_priority),
            "original_harmonic_supervision_verified": bool(
                original_harmonic_supervision_verified
            ),
            "original_harmonic_target": original_harmonic_target,
            "original_harmonic_include": bool(original_harmonic_include),
            "original_factor_review_status": original_factor_review_status,
        }
    )
    record["observation_candidate_key"] = _canonical_observation_key(record)
    return record


def normalize_s56_adjudicated_signals(frame: pd.DataFrame) -> pd.DataFrame:
    """Normalize the accepted S56 real table without changing its supervision."""

    _require(
        frame,
        {
            *_COMMON_REQUIRED,
            "human_label",
            "human_labeler",
            "is_injected_row",
            "source_uid",
        },
        name="S56 adjudicated table",
    )
    rows: list[dict[str, Any]] = []
    for source_record in frame.to_dict("records"):
        record = dict(source_record)
        if not _text(record.get("sector")):
            record["sector"] = 56
        if _truth(record.get("is_injected_row")):
            continue
        adjudicated = _text(record.get("human_label_adjudicated"))
        is_final = _truth(record.get("adjudication_final"))
        label = adjudicated if is_final and adjudicated else _text(record.get("human_label"))
        if label not in SIGNAL_REREVIEW_LABELS:
            continue
        labeler = (
            _text(record.get("adjudicated_labeler"))
            if is_final and adjudicated
            else _text(record.get("human_labeler"))
        )
        notes = (
            _text(record.get("adjudicated_notes"))
            if is_final and adjudicated
            else _text(record.get("human_notes"))
        )
        updated = (
            _text(record.get("adjudicated_updated_utc"))
            if is_final and adjudicated
            else _text(record.get("human_updated_utc"))
        )
        harmonic_include = _truth(record.get("harmonic_include_v1"))
        rows.append(
            _source_record(
                record,
                source_batch_id="s56_adjudicated_real343",
                source_uid=_text(record.get("source_uid")),
                source_label=label,
                source_labeler=labeler,
                source_notes=notes,
                source_updated_utc=updated,
                source_period_factor=record.get("effective_period_factor", 1),
                source_period_status=_text(record.get("adjudicated_period_status")),
                source_priority=100,
                original_harmonic_supervision_verified=harmonic_include,
                original_harmonic_target=_text(record.get("harmonic_target_v1")),
                original_harmonic_include=harmonic_include,
                original_factor_review_status=(
                    "verified_existing" if harmonic_include else "not_verified"
                ),
            )
        )
    return pd.DataFrame(rows)


def normalize_browser_signal_rows(
    queue: pd.DataFrame,
    labels: pd.DataFrame,
    *,
    source_batch_id: str,
    source_priority: int = 200,
    require_complete: bool = True,
) -> pd.DataFrame:
    """Normalize one complete repository-browser queue and label file."""

    _require(queue, {*_COMMON_REQUIRED, "row_id", "candidate_key"}, name="review queue")
    _require(
        labels,
        {
            "row_id",
            "candidate_key",
            "tic",
            "sector",
            "label",
            "labeler",
            "notes",
            "period_factor",
            "period_status",
            "updated_utc",
        },
        name="review labels",
    )
    public = queue.copy()
    returned = labels.copy()
    for name, table in (("queue", public), ("labels", returned)):
        table["row_id"] = table["row_id"].map(_text)
        if table["row_id"].eq("").any() or table["row_id"].duplicated().any():
            raise ValueError(f"{name} has blank or duplicate row_id values")
    public_ids = set(public["row_id"])
    returned_ids = set(returned["row_id"])
    if require_complete and public_ids != returned_ids:
        raise ValueError("review labels do not cover the exact source queue")
    if not returned_ids.issubset(public_ids):
        raise ValueError("review labels contain rows absent from the source queue")
    if not require_complete:
        public = public.loc[public["row_id"].isin(returned_ids)].copy()
    returned = returned.rename(
        columns={
            "candidate_key": "returned_candidate_key",
            "tic": "returned_tic",
            "sector": "returned_sector",
            "label": "returned_label",
            "labeler": "returned_labeler",
            "notes": "returned_notes",
            "period_factor": "returned_period_factor",
            "period_status": "returned_period_status",
            "updated_utc": "returned_updated_utc",
        }
    )
    joined = public.merge(returned, on="row_id", how="left", validate="one_to_one")
    expected_key = joined["candidate_key"].map(canonicalize_candidate_key)
    returned_key = joined["returned_candidate_key"].map(canonicalize_candidate_key)
    mismatch = (
        expected_key.ne(returned_key)
        | joined["tic"].map(_text).ne(joined["returned_tic"].map(_text))
        | joined["sector"].map(_text).ne(joined["returned_sector"].map(_text))
    )
    if mismatch.any():
        raise ValueError("review labels do not match their exact source queue")
    rows: list[dict[str, Any]] = []
    for record in joined.to_dict("records"):
        label = _text(record.get("returned_label"))
        if label not in SIGNAL_REREVIEW_LABELS:
            continue
        source_key = canonicalize_candidate_key(record.get("candidate_key"))
        rows.append(
            _source_record(
                record,
                source_batch_id=source_batch_id,
                source_uid=f"{source_batch_id}:{source_key}",
                source_label=label,
                source_labeler=_text(record.get("returned_labeler")),
                source_notes=_text(record.get("returned_notes")),
                source_updated_utc=_text(record.get("returned_updated_utc")),
                source_period_factor=record.get("returned_period_factor", 1),
                source_period_status=_text(record.get("returned_period_status")),
                source_priority=source_priority,
                original_harmonic_supervision_verified=False,
                original_harmonic_target="",
                original_harmonic_include=False,
                original_factor_review_status="not_explicitly_reviewed",
            )
        )
    return pd.DataFrame(rows)


def normalize_accepted_franklin_signals(frame: pd.DataFrame) -> pd.DataFrame:
    """Normalize one strict Franklin return while keeping its factors masked."""

    _require(
        frame,
        {
            *_COMMON_REQUIRED,
            "human_label",
            "human_labeler",
            "human_notes",
            "human_updated_utc",
            "source_batch_id",
            "source_uid",
        },
        name="accepted Franklin table",
    )
    rows: list[dict[str, Any]] = []
    for record in frame.to_dict("records"):
        label = _text(record.get("human_label"))
        if label not in SIGNAL_REREVIEW_LABELS:
            continue
        harmonic_include = _truth(record.get("harmonic_include_v1"))
        if harmonic_include or _truth(record.get("harmonic_supervision_verified")):
            raise ValueError("accepted Franklin morphology unexpectedly contains harmonic truth")
        rows.append(
            _source_record(
                record,
                source_batch_id=_text(record.get("source_batch_id")),
                source_uid=_text(record.get("source_uid")),
                source_label=label,
                source_labeler=_text(record.get("human_labeler")),
                source_notes=_text(record.get("human_notes")),
                source_updated_utc=_text(record.get("human_updated_utc")),
                source_period_factor=record.get("reported_period_factor", 1),
                source_period_status=_text(record.get("reported_period_status")),
                source_priority=100,
                original_harmonic_supervision_verified=False,
                original_harmonic_target="",
                original_harmonic_include=False,
                original_factor_review_status=_text(
                    record.get("factor_review_status")
                )
                or "not_explicitly_reviewed",
            )
        )
    return pd.DataFrame(rows)


def build_signal_rereview_queue(
    sources: Iterable[pd.DataFrame],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Build the public queue, source manifest, and required-sheet manifest."""

    nonempty = [frame.copy() for frame in sources if not frame.empty]
    if not nonempty:
        raise ValueError("no Planet-like or EB sources were supplied")
    provenance = pd.concat(nonempty, ignore_index=True, sort=False)
    if provenance["source_uid"].map(_text).eq("").any():
        raise ValueError("signal source contains a blank source_uid")
    if provenance["source_uid"].duplicated().any():
        raise ValueError("signal source_uid values are not unique")
    if not set(provenance["source_label"]).issubset(SIGNAL_REREVIEW_LABELS):
        raise ValueError("signal source contains an unsupported class")
    provenance["tic"] = pd.to_numeric(provenance["tic"], errors="raise").astype(np.int64)
    provenance["sector"] = pd.to_numeric(
        provenance["sector"], errors="raise"
    ).astype(np.int16)
    provenance["observation_candidate_key"] = provenance.apply(
        _canonical_observation_key, axis=1
    )
    if provenance["observation_candidate_key"].eq("").any():
        raise ValueError("signal source contains a blank observation identity")

    selected_rows: list[dict[str, Any]] = []
    for observation_key, group in provenance.groupby(
        "observation_candidate_key", sort=True, dropna=False
    ):
        ordered = group.sort_values(
            ["source_priority", "source_updated_utc", "source_uid"], kind="stable"
        )
        selected = ordered.iloc[-1].to_dict()
        verified = group.loc[
            group["original_harmonic_supervision_verified"].map(_truth)
            & group["original_harmonic_include"].map(_truth)
        ].copy()
        if not verified.empty:
            targets = {
                _text(value)
                for value in verified["original_harmonic_target"]
                if _text(value)
            }
            if len(targets) != 1:
                raise ValueError(
                    "exact candidate duplicate has conflicting verified "
                    f"harmonic targets: {sorted(targets)}"
                )
            harmonic_source = verified.sort_values(
                ["source_priority", "source_updated_utc", "source_uid"],
                kind="stable",
            ).iloc[-1]
            selected["original_harmonic_supervision_verified"] = True
            selected["original_harmonic_include"] = True
            selected["original_harmonic_target"] = next(iter(targets))
            selected["original_factor_review_status"] = _text(
                harmonic_source["original_factor_review_status"]
            )
            selected["source_period_factor"] = _text(
                harmonic_source["source_period_factor"]
            )
            selected["source_period_status"] = _text(
                harmonic_source["source_period_status"]
            )
        digest = hashlib.sha256(str(observation_key).encode()).hexdigest()[:16]
        selected["review_id"] = (
            f"s{int(selected['sector']):04d}-signal-final-{digest}"
        )
        selected["prior_label"] = selected["source_label"]
        selected["prior_labeler"] = selected["source_labeler"]
        selected["prior_notes"] = selected["source_notes"]
        selected["prior_period_factor"] = selected["source_period_factor"]
        selected["factor_review_status"] = selected[
            "original_factor_review_status"
        ]
        selected["initial_label"] = selected["source_label"]
        selected["initial_notes"] = selected["source_notes"]
        selected["initial_period_factor"] = selected["source_period_factor"]
        selected["initial_period_status"] = selected["source_period_status"]
        selected["prior_label_conflict"] = group["source_label"].nunique() > 1
        selected["n_source_records"] = int(len(group))
        selected["selected_source_uid"] = selected["source_uid"]
        selected["candidate_key"] = selected["observation_candidate_key"]
        selected_rows.append(selected)

    queue = pd.DataFrame(selected_rows).sort_values(
        ["sector", "prior_label", "tic", "period_d"], kind="stable"
    ).reset_index(drop=True)
    queue.insert(0, "row_id", np.arange(len(queue), dtype=int))
    same_sector_count = queue.groupby(["sector", "tic"])["tic"].transform("size")
    queue["same_sector_tic_count"] = same_sector_count.astype(int)
    queue["rereview_policy_version"] = SIGNAL_REREVIEW_POLICY_VERSION
    queue["precedence_policy_version"] = SIGNAL_PRECEDENCE_POLICY_VERSION

    row_by_key = queue.set_index("observation_candidate_key")["row_id"].to_dict()
    provenance["selected_queue_row_id"] = provenance[
        "observation_candidate_key"
    ].map(row_by_key)
    provenance["selected_for_queue"] = False
    selected_uid = set(queue["source_uid"].astype(str))
    provenance.loc[
        provenance["source_uid"].astype(str).isin(selected_uid),
        "selected_for_queue",
    ] = True

    public_columns = (
        "row_id",
        "review_id",
        "tic",
        "sector",
        "source_batch_id",
        "prior_label",
        "prior_labeler",
        "factor_review_status",
        "same_sector_tic_count",
        "cam",
        "ccd",
        "tmag",
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
        "twirl_vet_sheet_name",
        "initial_label",
        "initial_notes",
        "initial_period_factor",
        "initial_period_status",
        "prior_label_conflict",
        "n_source_records",
        "selected_source_uid",
        "observation_candidate_key",
        "original_harmonic_supervision_verified",
        "original_harmonic_target",
        "original_harmonic_include",
        "rereview_policy_version",
        "precedence_policy_version",
        "candidate_key",
    )
    public = queue.loc[:, [column for column in public_columns if column in queue]].copy()
    asset_manifest = public.loc[
        :,
        [
            "row_id",
            "sector",
            "tic",
            "source_batch_id",
            "twirl_vet_sheet_name",
        ],
    ].copy()
    asset_manifest["sheet_required"] = True
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "policy_version": SIGNAL_REREVIEW_POLICY_VERSION,
        "precedence_policy_version": SIGNAL_PRECEDENCE_POLICY_VERSION,
        "n_source_records": int(len(provenance)),
        "n_queue_rows": int(len(public)),
        "n_unique_tics": int(public["tic"].nunique()),
        "n_same_sector_multi_candidate_rows": int(
            public["same_sector_tic_count"].gt(1).sum()
        ),
        "n_prior_label_conflicts": int(public["prior_label_conflict"].sum()),
        "sector_counts": {
            str(key): int(value)
            for key, value in public["sector"].value_counts().sort_index().items()
        },
        "prior_label_counts": {
            str(key): int(value)
            for key, value in public["prior_label"].value_counts().sort_index().items()
        },
        "source_batch_counts": {
            str(key): int(value)
            for key, value in provenance["source_batch_id"].value_counts().sort_index().items()
        },
    }
    return public, provenance.reset_index(drop=True), asset_manifest, summary


def standalone_app_candidate_key(row: Mapping[str, Any]) -> str:
    """Reproduce the exact identity used by ``franklin_vetting_app.py``."""

    return "|".join(
        _text(row.get(column, ""))
        for column in ("review_id", "tic", "sector", "period_d", "t0_bjd")
    )


def finalize_signal_rereview(
    queue: pd.DataFrame,
    labels: pd.DataFrame,
    *,
    adjudicator: str,
    accepted_utc: str,
) -> pd.DataFrame:
    """Strictly join a complete final pass without changing harmonic authority."""

    _require(
        queue,
        {
            "row_id",
            "review_id",
            "tic",
            "sector",
            "period_d",
            "t0_bjd",
            "prior_label",
            "original_harmonic_supervision_verified",
            "original_harmonic_target",
            "original_harmonic_include",
        },
        name="signal re-review queue",
    )
    _require(
        labels,
        {
            "row_id",
            "candidate_key",
            "tic",
            "sector",
            "label",
            "labeler",
            "notes",
            "period_factor",
            "period_status",
            "updated_utc",
        },
        name="signal re-review labels",
    )
    if not _text(adjudicator) or not _text(accepted_utc):
        raise ValueError("adjudicator and accepted_utc must be nonblank")
    public = queue.copy()
    returned = labels.copy()
    for name, table in (("queue", public), ("labels", returned)):
        table["row_id"] = table["row_id"].map(_text)
        if table["row_id"].eq("").any() or table["row_id"].duplicated().any():
            raise ValueError(f"{name} has blank or duplicate row_id values")
    if set(public["row_id"]) != set(returned["row_id"]):
        raise ValueError("final signal review is incomplete or contains extra rows")
    expected = public.apply(standalone_app_candidate_key, axis=1)
    expected_by_id = dict(zip(public["row_id"], expected))
    observed = returned["row_id"].map(expected_by_id)
    if returned["candidate_key"].map(_text).ne(observed.map(_text)).any():
        raise ValueError("final signal labels do not match the exact frozen queue")
    if returned["tic"].map(_text).ne(
        returned["row_id"].map(public.set_index("row_id")["tic"]).map(_text)
    ).any():
        raise ValueError("final signal labels contain a TIC mismatch")
    if returned["sector"].map(_text).ne(
        returned["row_id"].map(public.set_index("row_id")["sector"]).map(_text)
    ).any():
        raise ValueError("final signal labels contain a sector mismatch")
    final_labels = returned["label"].map(_text)
    invalid = sorted(set(final_labels) - FINAL_LABELS)
    if invalid or final_labels.eq("").any():
        raise ValueError(f"final signal review contains invalid labels: {invalid}")
    observed_labelers = set(returned["labeler"].map(_text))
    if observed_labelers != {adjudicator}:
        raise ValueError(
            "final signal review labeler differs from the declared adjudicator: "
            f"{sorted(observed_labelers)}"
        )

    returned = returned.rename(
        columns={
            "candidate_key": "final_app_candidate_key",
            "label": "final_human_label",
            "labeler": "final_labeler",
            "notes": "final_notes",
            "period_factor": "rereview_period_factor",
            "period_status": "rereview_period_status",
            "updated_utc": "final_updated_utc",
            "tic": "returned_tic",
            "sector": "returned_sector",
        }
    )
    joined = public.merge(returned, on="row_id", how="left", validate="one_to_one")
    joined["morphology_adjudicator"] = adjudicator
    joined["morphology_accepted_utc"] = accepted_utc
    joined["morphology_review_status"] = "explicit_row_review"
    joined["human_label"] = joined["final_human_label"]
    joined["human_label_source"] = "human_adjudication"
    joined["human_labeler"] = joined["final_labeler"]
    joined["human_notes"] = joined["final_notes"]
    joined["human_updated_utc"] = joined["final_updated_utc"]
    joined["morphology_target_v1"] = (
        joined["human_label"].map(MORPHOLOGY_TARGET_MAP_V1).fillna("")
    )
    joined["morphology_include_v1"] = joined["morphology_target_v1"].ne("")
    joined["preserve_target_v1"] = ""
    joined.loc[
        joined["human_label"].isin(PRESERVE_LABELS_V1), "preserve_target_v1"
    ] = "preserve"
    joined.loc[
        joined["human_label"].isin(REJECT_LABELS_V1), "preserve_target_v1"
    ] = "reject"
    joined["preserve_include_v1"] = joined["preserve_target_v1"].ne("")
    original_verified = joined[
        "original_harmonic_supervision_verified"
    ].map(_truth)
    keep_harmonic = (
        original_verified
        & joined["human_label"].isin(TRANSIT_HARMONIC_LABELS_V1)
        & joined["original_harmonic_target"].map(_text).ne("")
    )
    joined["harmonic_supervision_verified"] = keep_harmonic
    joined["harmonic_target_v1"] = joined["original_harmonic_target"].where(
        keep_harmonic, ""
    )
    joined["harmonic_include_v1"] = keep_harmonic
    joined["broad_preserve_only"] = joined["human_label"].eq("wide_transit_like")
    joined["model_target_policy_version"] = HARMONIC_CNN_TARGET_POLICY
    joined["rereview_period_is_audit_only"] = True
    return joined.drop(columns=["returned_tic", "returned_sector"])


__all__ = [
    "FINAL_LABELS",
    "SIGNAL_PRECEDENCE_POLICY_VERSION",
    "SIGNAL_REREVIEW_LABELS",
    "SIGNAL_REREVIEW_POLICY_VERSION",
    "build_signal_rereview_queue",
    "finalize_signal_rereview",
    "normalize_accepted_franklin_signals",
    "normalize_browser_signal_rows",
    "normalize_s56_adjudicated_signals",
    "standalone_app_candidate_key",
]
