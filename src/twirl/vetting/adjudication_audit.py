"""Post-review audit and training-table assembly for S56 adjudication."""
from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from twirl.vetting.adjudication import (
    ADJUDICATION_LEAKAGE_COLUMNS,
    ADP_APERTURES,
    SIGNAL_LABELS,
    join_browser_labels,
    read_table,
    reanchor_deprecated_sources,
)
from twirl.vetting.adp_only import canonical_det_flux_columns
from twirl.vetting.label_io import (
    ADJUDICATION_LABEL_COLUMNS,
    latest_label_records,
    normalize_review_queue,
    validate_label_records,
)
from twirl.vetting.recovery50_teacher import add_deterministic_splits, add_label_roles


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except TypeError:
        pass
    return str(value)


def _as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).astype(bool)
    return series.fillna("").astype(str).str.strip().str.lower().isin({"1", "true", "yes", "y"})


def normalize_period_factor(value: Any) -> tuple[str, str, float]:
    text = _text(value).strip().lower() or "1"
    aliases = {"0.250": "0.25", "0.50": "0.5", "1.0": "1", "2.0": "2", "4.0": "4", "?": "unresolved"}
    text = aliases.get(text, text)
    if text == "unresolved":
        return text, "unresolved", float("nan")
    try:
        factor = float(text)
    except ValueError as exc:
        raise ValueError(f"invalid adjudication period factor: {value}") from exc
    if factor not in {0.25, 0.5, 1.0, 2.0, 4.0}:
        raise ValueError(f"invalid adjudication period factor: {value}")
    return f"{factor:g}", "review_period" if factor == 1.0 else "refolded", factor


def load_adjudication_reviews(
    *,
    queue_csv: Path,
    labels_csv: Path,
    manifest_csv: Path,
) -> pd.DataFrame:
    """Load and strictly validate row-level adjudication decisions."""

    queue = normalize_review_queue(read_table(queue_csv))
    labels = latest_label_records(labels_csv)
    validate_label_records(queue, labels)
    manifest = read_table(manifest_csv)
    if manifest["row_id"].duplicated().any():
        raise ValueError("adjudication manifest has duplicate normalized row IDs")
    if len(manifest) != len(queue):
        raise ValueError(f"manifest rows={len(manifest)} queue rows={len(queue)}")
    manifest = manifest.set_index("row_id").reindex(queue["row_id"]).reset_index()
    if manifest["source_uid"].isna().any():
        raise ValueError("adjudication manifest is missing normalized queue rows")
    expected_keys = queue.set_index("row_id")["candidate_key"].astype(str)
    manifest_keys = manifest.set_index("row_id")["candidate_key"].fillna("").astype(str)
    if not manifest_keys.equals(expected_keys):
        raise ValueError("adjudication manifest candidate keys do not match the normalized queue")
    label_columns = {
        "label": "adjudicated_row_label",
        "label_source": "adjudicated_row_label_source",
        "labeler": "adjudicated_row_labeler",
        "notes": "adjudicated_row_notes",
        "updated_utc": "adjudicated_row_updated_utc",
        "period_factor": "adjudicated_row_period_factor",
        "period_status": "adjudicated_row_period_status_saved",
    }
    labels = labels.rename(columns=label_columns)
    keep = ["row_id", *label_columns.values()]
    rows = manifest.merge(labels.loc[:, keep], on="row_id", how="left", validate="one_to_one")
    for column in label_columns.values():
        if column not in rows:
            rows[column] = ""
        rows[column] = rows[column].fillna("").astype(str)

    factors: list[str] = []
    statuses: list[str] = []
    periods: list[float] = []
    for _, row in rows.iterrows():
        if not _text(row["adjudicated_row_label"]):
            factors.append("")
            statuses.append("")
            periods.append(float("nan"))
            continue
        factor_text, status, factor = normalize_period_factor(row["adjudicated_row_period_factor"])
        review_period = float(row["review_period_d"])
        factors.append(factor_text)
        statuses.append(status)
        periods.append(review_period * factor if np.isfinite(factor) else float("nan"))
    rows["adjudicated_row_period_factor"] = factors
    rows["adjudicated_row_period_status"] = statuses
    rows["adjudicated_row_period_d"] = periods
    return rows


def _cohen_kappa(left: pd.Series, right: pd.Series) -> float:
    left = left.fillna("").astype(str)
    right = right.fillna("").astype(str)
    if len(left) == 0:
        return float("nan")
    observed = float(left.eq(right).mean())
    labels = sorted(set(left) | set(right))
    expected = sum(float(left.eq(label).mean()) * float(right.eq(label).mean()) for label in labels)
    if np.isclose(expected, 1.0):
        return 1.0 if np.isclose(observed, 1.0) else float("nan")
    return (observed - expected) / (1.0 - expected)


def repeat_agreement(rows: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    repeated = rows.loc[rows["repeat_group"].fillna("").astype(str).ne("")].copy()
    records: list[dict[str, Any]] = []
    for group, part in repeated.groupby("repeat_group", sort=True):
        if len(part) != 2:
            raise ValueError(f"repeat group {group} has {len(part)} rows instead of 2")
        original = part.loc[part["repeat_occurrence"].eq("original")]
        blind = part.loc[part["repeat_occurrence"].eq("blind_repeat")]
        if len(original) != 1 or len(blind) != 1:
            raise ValueError(f"repeat group {group} lacks one original and one blind repeat")
        left = original.iloc[0]
        right = blind.iloc[0]
        left_label = _text(left["adjudicated_row_label"])
        right_label = _text(right["adjudicated_row_label"])
        left_factor = _text(left["adjudicated_row_period_factor"])
        right_factor = _text(right["adjudicated_row_period_factor"])
        complete = bool(left_label and right_label)
        label_agreement = complete and left_label == right_label
        factor_agreement = complete and left_factor == right_factor
        records.append(
            {
                "repeat_group": group,
                "source_uid": left["source_uid"],
                "reference_label": left.get("repeat_reference_label", ""),
                "original_row_id": int(left["row_id"]),
                "repeat_row_id": int(right["row_id"]),
                "original_label": left_label,
                "repeat_label": right_label,
                "original_period_factor": left_factor,
                "repeat_period_factor": right_factor,
                "complete": complete,
                "label_agreement": label_agreement,
                "factor_agreement": factor_agreement,
                "exact_agreement": label_agreement and factor_agreement,
                "discordant": complete and not (label_agreement and factor_agreement),
            }
        )
    pairs = pd.DataFrame(records)
    complete_pairs = pairs.loc[_as_bool(pairs["complete"])].copy() if len(pairs) else pairs
    per_class: dict[str, dict[str, Any]] = {}
    if len(complete_pairs):
        for label, part in complete_pairs.groupby("reference_label", sort=True):
            per_class[str(label)] = {
                "n": int(len(part)),
                "label_agreement": float(_as_bool(part["label_agreement"]).mean()),
                "exact_agreement": float(_as_bool(part["exact_agreement"]).mean()),
            }
    summary = {
        "n_repeat_pairs": int(len(pairs)),
        "n_complete_pairs": int(len(complete_pairs)),
        "n_discordant_pairs": int(_as_bool(complete_pairs["discordant"]).sum()) if len(complete_pairs) else 0,
        "label_exact_agreement": float(_as_bool(complete_pairs["label_agreement"]).mean()) if len(complete_pairs) else float("nan"),
        "label_and_factor_exact_agreement": float(_as_bool(complete_pairs["exact_agreement"]).mean()) if len(complete_pairs) else float("nan"),
        "cohen_kappa": _cohen_kappa(complete_pairs["original_label"], complete_pairs["repeat_label"]),
        "per_reference_class": per_class,
    }
    return pairs, summary


def write_repeat_resolution_queue(
    *,
    queue_csv: Path,
    manifest_csv: Path,
    repeat_pairs: pd.DataFrame,
    out_dir: Path,
) -> dict[str, Any]:
    queue = normalize_review_queue(read_table(queue_csv))
    manifest = read_table(manifest_csv)
    discordant = repeat_pairs.loc[_as_bool(repeat_pairs["discordant"])].copy()
    out_dir.mkdir(parents=True, exist_ok=True)
    queue_path = out_dir / "review_queue_repeat_resolution.csv"
    manifest_path = out_dir / "repeat_resolution_manifest_private.csv"
    labels_path = out_dir / "human_labels_repeat_resolution.csv"
    if discordant.empty:
        resolution_queue = queue.head(0).copy()
        resolution_manifest = manifest.head(0).copy()
    else:
        source_rows: list[pd.Series] = []
        hidden_rows: list[pd.Series] = []
        queue_by_id = queue.set_index("row_id")
        manifest_by_id = manifest.set_index("row_id")
        for resolution_index, pair in discordant.reset_index(drop=True).iterrows():
            source_row_id = int(pair["original_row_id"])
            public = queue_by_id.loc[source_row_id].copy()
            public["review_id"] = f"resolution:{resolution_index:03d}"
            public["source_bucket"] = "real_adjudication_resolution"
            source_rows.append(public)
            hidden = manifest_by_id.loc[source_row_id].copy()
            hidden["resolution_index"] = resolution_index
            hidden["repeat_group"] = pair["repeat_group"]
            hidden_rows.append(hidden)
        resolution_queue = normalize_review_queue(pd.DataFrame(source_rows).reset_index(drop=True))
        resolution_manifest = pd.DataFrame(hidden_rows).reset_index(drop=True)
        resolution_manifest.insert(0, "row_id", resolution_queue["row_id"].to_numpy())
        resolution_manifest["candidate_key"] = resolution_queue["candidate_key"].to_numpy()
    resolution_queue.to_csv(queue_path, index=False)
    resolution_manifest.to_csv(manifest_path, index=False)
    if not labels_path.exists():
        pd.DataFrame(columns=ADJUDICATION_LABEL_COLUMNS).to_csv(labels_path, index=False)
    return {
        "n_resolution_rows": int(len(resolution_queue)),
        "queue": str(queue_path),
        "manifest": str(manifest_path),
        "labels": str(labels_path),
    }


def _load_resolution_decisions(resolution_dir: Path) -> pd.DataFrame:
    queue_path = resolution_dir / "review_queue_repeat_resolution.csv"
    manifest_path = resolution_dir / "repeat_resolution_manifest_private.csv"
    labels_path = resolution_dir / "human_labels_repeat_resolution.csv"
    columns = [
        "source_uid",
        "resolution_label",
        "resolution_period_factor",
        "resolution_period_status",
        "resolution_period_d",
        "resolution_labeler",
        "resolution_updated_utc",
    ]
    if not queue_path.exists() or not manifest_path.exists() or not labels_path.exists():
        return pd.DataFrame(columns=columns)
    queue = normalize_review_queue(read_table(queue_path))
    labels = latest_label_records(labels_path)
    validate_label_records(queue, labels)
    manifest = read_table(manifest_path)
    joined = manifest.merge(labels, on="row_id", how="left", validate="one_to_one")
    records: list[dict[str, Any]] = []
    for _, row in joined.iterrows():
        label = _text(row.get("label"))
        if not label:
            continue
        factor_text, status, factor = normalize_period_factor(row.get("period_factor"))
        review_period = float(row["review_period_d"])
        records.append(
            {
                "source_uid": row["source_uid"],
                "resolution_label": label,
                "resolution_period_factor": factor_text,
                "resolution_period_status": status,
                "resolution_period_d": review_period * factor if np.isfinite(factor) else np.nan,
                "resolution_labeler": _text(row.get("labeler")),
                "resolution_updated_utc": _text(row.get("updated_utc")),
            }
        )
    return pd.DataFrame(records, columns=columns)


def final_unique_decisions(rows: pd.DataFrame, resolution: pd.DataFrame) -> pd.DataFrame:
    if len(resolution) and resolution["source_uid"].duplicated().any():
        raise ValueError("repeat resolution contains duplicate source UIDs")
    resolution_by_uid = resolution.set_index("source_uid") if len(resolution) else None
    records: list[dict[str, Any]] = []
    for source_uid, part in rows.groupby("source_uid", sort=False):
        base = part.sort_values("row_id", kind="stable").iloc[0].to_dict()
        labeled = part.loc[part["adjudicated_row_label"].fillna("").astype(str).ne("")].copy()
        decision: dict[str, Any] | None = None
        decision_source = ""
        if len(part) == 1 and len(labeled) == 1:
            decision = labeled.iloc[0].to_dict()
            decision_source = "single_review"
        elif len(part) == 2 and len(labeled) == 2:
            same_label = labeled["adjudicated_row_label"].nunique(dropna=False) == 1
            same_factor = labeled["adjudicated_row_period_factor"].nunique(dropna=False) == 1
            if same_label and same_factor:
                decision = labeled.iloc[0].to_dict()
                decision_source = "blind_repeat_agreement"
            elif resolution_by_uid is not None and source_uid in resolution_by_uid.index:
                resolved = resolution_by_uid.loc[source_uid]
                base.update(
                    {
                        "adjudicated_row_label": resolved["resolution_label"],
                        "adjudicated_row_period_factor": resolved["resolution_period_factor"],
                        "adjudicated_row_period_status": resolved["resolution_period_status"],
                        "adjudicated_row_period_d": resolved["resolution_period_d"],
                        "adjudicated_row_labeler": resolved["resolution_labeler"],
                        "adjudicated_row_updated_utc": resolved["resolution_updated_utc"],
                    }
                )
                decision = base
                decision_source = "third_pass_resolution"
        base["human_label_adjudicated"] = _text(decision.get("adjudicated_row_label")) if decision else ""
        base["adjudicated_period_factor"] = _text(decision.get("adjudicated_row_period_factor")) if decision else ""
        base["adjudicated_period_status"] = _text(decision.get("adjudicated_row_period_status")) if decision else ""
        base["adjudicated_period_d"] = decision.get("adjudicated_row_period_d", np.nan) if decision else np.nan
        base["adjudicated_labeler"] = _text(decision.get("adjudicated_row_labeler")) if decision else ""
        base["adjudicated_updated_utc"] = _text(decision.get("adjudicated_row_updated_utc")) if decision else ""
        base["adjudicated_decision_source"] = decision_source
        base["adjudication_final"] = bool(decision)
        records.append(base)
    return pd.DataFrame(records)


def _transition_table(frame: pd.DataFrame) -> pd.DataFrame:
    complete = frame.loc[_as_bool(frame["adjudication_final"])].copy()
    return pd.crosstab(
        complete["raw_human_label"].fillna("").astype(str),
        complete["human_label_adjudicated"].fillna("").astype(str),
        rownames=["raw_human_label"],
        colnames=["human_label_adjudicated"],
        dropna=False,
    )


def write_transition_audit(final: pd.DataFrame, out_dir: Path) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    _transition_table(final).to_csv(out_dir / "transition_matrix_overall.csv")
    for column, prefix in (
        ("source_cohort", "source_cohort"),
        ("origin_queue", "origin_queue"),
        ("labeling_era", "labeling_era"),
    ):
        for value, part in final.groupby(column, dropna=False, sort=True):
            safe = _text(value).replace("/", "_").replace(" ", "_") or "missing"
            _transition_table(part).to_csv(out_dir / f"transition_matrix_{prefix}_{safe}.csv")
    old_label = final["raw_human_label"].fillna("").astype(str)
    new_label = final["human_label_adjudicated"].fillna("").astype(str)
    final_mask = _as_bool(final["adjudication_final"])

    def count(old: str, new: str) -> int:
        return int((final_mask & old_label.eq(old) & new_label.eq(new)).sum())

    dedicated = {
        "planet_like_to_wide_transit_like": count("planet_like", "wide_transit_like"),
        "planet_like_to_eclipsing_binary_or_pceb": count("planet_like", "eclipsing_binary_or_pceb"),
        "wide_transit_like_to_eclipsing_binary_or_pceb": count("wide_transit_like", "eclipsing_binary_or_pceb"),
        "eclipsing_binary_or_pceb_to_stellar_variability": count("eclipsing_binary_or_pceb", "stellar_variability"),
        "uncertain_to_instrumental_or_systematic": count("uncertain", "instrumental_or_systematic"),
        "instrumental_or_systematic_to_uncertain": count("instrumental_or_systematic", "uncertain"),
    }
    harmonic = pd.crosstab(
        final.loc[final_mask, "human_label_adjudicated"].fillna("").astype(str),
        final.loc[final_mask, "adjudicated_period_factor"].fillna("").astype(str),
        rownames=["human_label_adjudicated"],
        colnames=["adjudicated_period_factor"],
        dropna=False,
    )
    harmonic.to_csv(out_dir / "harmonic_factor_by_adjudicated_class.csv")
    (out_dir / "dedicated_transitions.json").write_text(json.dumps(dedicated, indent=2, sort_keys=True) + "\n")
    return dedicated


def _set_current_adp_ephemeris(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    for prefix in ("model_", "display_", "anchor_", "queue_"):
        out[f"{prefix}period_d"] = out["period_d"]
        out[f"{prefix}t0_bjd"] = out["t0_bjd"]
        out[f"{prefix}duration_min"] = out["duration_min"]
    out["model_ephemeris_source"] = "current_adp_sml_bls"
    out["display_ephemeris_source"] = "current_adp_sml_bls"
    out["display_anchor_aperture"] = ADP_APERTURES[0]
    out["anchor_aperture"] = ADP_APERTURES[0]
    out["tensor_apertures"] = ",".join(ADP_APERTURES)
    out["adp_only_contract_version"] = "s56_adp_pair_v1"
    return out


def build_adjudicated_training_table(
    *,
    active_training_table: Path,
    eb_pilot_queue: Path,
    eb_pilot_labels: Path,
    deprecated_eb_queue: Path,
    deprecated_eb_labels: Path,
    current_adp_bls: Path,
    final_decisions: pd.DataFrame,
) -> pd.DataFrame:
    """Overlay final decisions once and retain untouched real/injected labels."""

    active = read_table(active_training_table).copy()
    active["source_uid"] = "teacher2k:" + active["review_id"].fillna("").astype(str)
    active["human_label_raw"] = active["human_label"].fillna("").astype(str)
    active["human_notes_raw"] = active.get("human_notes", pd.Series("", index=active.index)).fillna("").astype(str)
    active["human_updated_utc_raw"] = active.get(
        "human_updated_utc", pd.Series("", index=active.index)
    ).fillna("").astype(str)

    pilot = join_browser_labels(eb_pilot_queue, eb_pilot_labels)
    pilot["source_uid"] = "eb_pilot_adp:" + pilot["review_id"].fillna("").astype(str)
    pilot["human_label"] = pilot["browser_label"].fillna("").astype(str)
    pilot["human_label_source"] = pilot["browser_label_source"].fillna("").astype(str)
    pilot["human_labeler"] = pilot["browser_labeler"].fillna("").astype(str)
    pilot["human_notes"] = pilot["browser_notes"].fillna("").astype(str)
    pilot["human_updated_utc"] = pilot["browser_updated_utc"].fillna("").astype(str)
    pilot["human_label_raw"] = pilot["human_label"]
    pilot["human_notes_raw"] = pilot["human_notes"]
    pilot["human_updated_utc_raw"] = pilot["human_updated_utc"]
    pilot["source_kind"] = "real_candidate"
    pilot["is_injected_row"] = False
    pilot = _set_current_adp_ephemeris(pilot)

    deprecated = join_browser_labels(deprecated_eb_queue, deprecated_eb_labels)
    deprecated["ephemeris_source"] = "deprecated_queue_pending_rebuild"
    deprecated = reanchor_deprecated_sources(deprecated, read_table(current_adp_bls))
    deprecated = _set_current_adp_ephemeris(deprecated)
    deprecated["source_uid"] = "eb_deprecated:" + deprecated["review_id"].fillna("").astype(str)
    deprecated["human_label"] = deprecated["browser_label"].fillna("").astype(str)
    deprecated["human_label_source"] = deprecated["browser_label_source"].fillna("").astype(str)
    deprecated["human_labeler"] = deprecated["browser_labeler"].fillna("").astype(str)
    deprecated["human_notes"] = deprecated["browser_notes"].fillna("").astype(str)
    deprecated["human_updated_utc"] = deprecated["browser_updated_utc"].fillna("").astype(str)
    deprecated["human_label_raw"] = deprecated["human_label"]
    deprecated["human_notes_raw"] = deprecated["human_notes"]
    deprecated["human_updated_utc_raw"] = deprecated["human_updated_utc"]
    deprecated["source_kind"] = "real_candidate"
    deprecated["is_injected_row"] = False

    master = pd.concat([active, pilot, deprecated], ignore_index=True, sort=False)
    master = master.loc[master["human_label_raw"].fillna("").astype(str).ne("")].copy()
    if master["source_uid"].duplicated().any():
        raise ValueError("adjudicated training master contains duplicate source UIDs")
    canonical_columns = canonical_det_flux_columns(master.columns)
    if canonical_columns:
        master = master.drop(columns=canonical_columns)

    decision_columns = [
        "source_uid",
        "human_label_adjudicated",
        "adjudicated_period_factor",
        "adjudicated_period_status",
        "adjudicated_period_d",
        "adjudicated_labeler",
        "adjudicated_updated_utc",
        "adjudication_final",
        "adjudicated_decision_source",
    ]
    decisions = final_decisions.loc[:, decision_columns].drop_duplicates("source_uid", keep="last")
    master = master.merge(decisions, on="source_uid", how="left", validate="one_to_one")
    for column in decision_columns[1:]:
        if column not in master:
            master[column] = ""
    final_mask = _as_bool(master["adjudication_final"])
    master.loc[final_mask, "human_label"] = master.loc[final_mask, "human_label_adjudicated"]
    resolved = final_mask & master["adjudicated_period_status"].fillna("").astype(str).ne("unresolved")
    master.loc[resolved, "model_period_d"] = pd.to_numeric(
        master.loc[resolved, "adjudicated_period_d"], errors="coerce"
    )
    master.loc[resolved, "model_t0_bjd"] = master.loc[resolved, "t0_bjd"]
    master.loc[resolved, "model_duration_min"] = master.loc[resolved, "duration_min"]
    master.loc[resolved, "model_ephemeris_source"] = "human_adjudicated_harmonic"
    factor = master["adjudicated_period_factor"].fillna("").astype(str)
    master["harmonic_suspect"] = final_mask & (factor.ne("1") | master["adjudicated_period_status"].eq("unresolved"))
    master["refold_factor"] = pd.to_numeric(factor, errors="coerce")
    master["refold_period_d"] = pd.to_numeric(master["adjudicated_period_d"], errors="coerce")
    master["ephemeris_status"] = master["adjudicated_period_status"].fillna("").astype(str)

    master = add_label_roles(master)
    label = master["human_label"].fillna("").astype(str)
    unresolved_signal = master["adjudicated_period_status"].fillna("").astype(str).eq("unresolved") & label.isin(
        SIGNAL_LABELS
    )
    master["preserve_training_include"] = ~label.isin({"", "skip"})
    master["compact_morphology_include"] = master["main_teacher_include"].fillna(False).astype(bool) & ~unresolved_signal
    master.loc[unresolved_signal, "main_teacher_target"] = ""
    master.loc[unresolved_signal, "main_teacher_include"] = False
    master = add_deterministic_splits(master, validation_fraction=0.20, test_fraction=0.20, random_state=56)
    real = master["source_kind"].fillna("").astype(str).eq("real_candidate")
    master["eb_cnn_include"] = real & ~label.isin({"", "skip"}) & ~unresolved_signal
    master["eb_miner_target"] = np.where(
        master["eb_cnn_include"],
        np.where(label.eq("eclipsing_binary_or_pceb"), "eb_pceb", "not_eb"),
        "",
    )
    master["adjudication_leakage_columns"] = ",".join(sorted(ADJUDICATION_LEAKAGE_COLUMNS))
    master["source_row_id"] = master.get("row_id", pd.Series(np.nan, index=master.index))
    master["row_id"] = np.arange(len(master), dtype=int)
    return master.reset_index(drop=True)


def run_adjudication_audit(
    *,
    root: Path,
    active_training_table: Path,
    eb_pilot_queue: Path,
    eb_pilot_labels: Path,
    deprecated_eb_queue: Path,
    deprecated_eb_labels: Path,
    current_adp_bls: Path,
) -> dict[str, Any]:
    root = Path(root)
    audit_dir = root / "audit"
    audit_dir.mkdir(parents=True, exist_ok=True)
    rows = load_adjudication_reviews(
        queue_csv=root / "review_queue_real343.csv",
        labels_csv=root / "human_labels_adjudicated.csv",
        manifest_csv=root / "adjudication_manifest_private.csv",
    )
    rows.to_csv(audit_dir / "adjudication_rows_joined_private.csv", index=False)
    pairs, repeat_summary = repeat_agreement(rows)
    pairs.to_csv(audit_dir / "blind_repeat_pairs.csv", index=False)
    resolution_summary = write_repeat_resolution_queue(
        queue_csv=root / "review_queue_real343.csv",
        manifest_csv=root / "adjudication_manifest_private.csv",
        repeat_pairs=pairs,
        out_dir=root / "repeat_resolution",
    )
    resolution = _load_resolution_decisions(root / "repeat_resolution")
    final = final_unique_decisions(rows, resolution)
    final.to_csv(audit_dir / "adjudicated_unique_sources.csv", index=False)
    dedicated = write_transition_audit(final, audit_dir)
    training = build_adjudicated_training_table(
        active_training_table=active_training_table,
        eb_pilot_queue=eb_pilot_queue,
        eb_pilot_labels=eb_pilot_labels,
        deprecated_eb_queue=deprecated_eb_queue,
        deprecated_eb_labels=deprecated_eb_labels,
        current_adp_bls=current_adp_bls,
        final_decisions=final,
    )
    training_dir = root / "adjudicated_training_table"
    training_dir.mkdir(parents=True, exist_ok=True)
    training_path = training_dir / "human_vetting_training_table_adjudicated.csv"
    training.to_csv(training_path, index=False)
    leakage_path = training_dir / "model_input_leakage_columns.json"
    leakage_path.write_text(
        json.dumps({"columns": sorted(ADJUDICATION_LEAKAGE_COLUMNS)}, indent=2, sort_keys=True) + "\n"
    )

    labeled_rows = rows["adjudicated_row_label"].fillna("").astype(str).ne("")
    final_sources = _as_bool(final["adjudication_final"])
    pending_repeat_resolution = int(resolution_summary["n_resolution_rows"] - len(resolution))
    cleanup_complete = bool(
        int(labeled_rows.sum()) == len(rows)
        and int(final_sources.sum()) == len(final)
        and pending_repeat_resolution == 0
    )
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "n_display_rows": int(len(rows)),
        "n_labeled_rows": int(labeled_rows.sum()),
        "n_unique_sources": int(len(final)),
        "n_final_unique_sources": int(final_sources.sum()),
        "repeat_agreement": repeat_summary,
        "repeat_resolution": resolution_summary,
        "n_completed_resolution_rows": int(len(resolution)),
        "n_pending_resolution_rows": pending_repeat_resolution,
        "dedicated_transitions": dedicated,
        "cleanup_complete": cleanup_complete,
        "training_rows": int(len(training)),
        "training_human_label_counts": {
            str(key): int(value)
            for key, value in training["human_label"].fillna("").astype(str).value_counts().sort_index().items()
        },
        "outputs": {
            "joined_rows": str(audit_dir / "adjudication_rows_joined_private.csv"),
            "repeat_pairs": str(audit_dir / "blind_repeat_pairs.csv"),
            "unique_sources": str(audit_dir / "adjudicated_unique_sources.csv"),
            "training_table": str(training_path),
            "leakage_columns": str(leakage_path),
        },
    }
    (audit_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n"
    )
    return summary


__all__ = [
    "build_adjudicated_training_table",
    "final_unique_decisions",
    "load_adjudication_reviews",
    "normalize_period_factor",
    "repeat_agreement",
    "run_adjudication_audit",
    "write_repeat_resolution_queue",
]
