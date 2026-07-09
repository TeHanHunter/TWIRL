#!/usr/bin/env python3
"""Build a focused revisit queue for real S56 rows labeled planet-like."""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))


DEFAULT_TRAINING_TABLE = (
    REPO_ROOT
    / "reports/stage5_validation/s56_recovery50_teacher_queue_2k/"
    / "human_training_table/human_vetting_training_table.csv"
)
DEFAULT_OUTPUT_DIR = (
    REPO_ROOT
    / "reports/stage5_validation/s56_recovery50_teacher_queue_2k/"
    / "real_planet_revisit_wide"
)
DEFAULT_HALF_PERIOD_REVIEW_IDS = (
    "real:1400739415",
    "real:1961800465",
    "real:160520962",
    "real:1400855795",
)
DEFAULT_SHEET_BRANCH_NAME = "current_adp"


def _clean(value: Any) -> Any:
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


def _finite_float(value: Any) -> float | None:
    value = _clean(value)
    if value is None:
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if np.isfinite(out) else None


def _display_ephemeris(row: pd.Series) -> tuple[float | None, float | None, float | None]:
    period = _finite_float(row.get("display_period_d"))
    t0 = _finite_float(row.get("display_t0_bjd"))
    duration = _finite_float(row.get("display_duration_min"))
    if period is None:
        period = _finite_float(row.get("period_d"))
    if t0 is None:
        t0 = _finite_float(row.get("t0_bjd"))
    if duration is None:
        duration = _finite_float(row.get("duration_min"))
    return period, t0, duration


def _set_if_present_or_new(row: dict[str, Any], key: str, value: Any) -> None:
    row[key] = "" if value is None else value


def _base_revisit_row(
    source: pd.Series,
    *,
    variant: str,
    period_d: float | None,
    t0_bjd: float | None,
    duration_min: float | None,
    twirl_vet_sheet_name: str,
    notes: str,
    sheet_branch_name: str = DEFAULT_SHEET_BRANCH_NAME,
) -> dict[str, Any]:
    out = source.to_dict()
    original_review_id = str(_clean(source.get("review_id")) or "")
    source_bucket = f"real_planet_revisit_{variant}"
    review_suffix = "current_fold" if variant == "current" else variant

    out["review_id"] = (
        f"{original_review_id}:{review_suffix}" if original_review_id else f"real_planet_revisit:{review_suffix}"
    )
    safe_review_id = out["review_id"].replace("/", "_").replace(":", "_")
    out["revisit_source_review_id"] = original_review_id
    out["revisit_variant"] = variant
    out["revisit_created_utc"] = datetime.now(timezone.utc).isoformat()
    out["previous_human_label"] = str(_clean(source.get("human_label")) or "")
    out["previous_human_notes"] = str(_clean(source.get("human_notes")) or "")
    out["source_bucket"] = source_bucket
    out["source_kind"] = "real_candidate"
    out["vet_class"] = (
        "revisit: current fold"
        if variant == "current"
        else "revisit: half-period fold"
    )
    out["label"] = ""
    out["label_source"] = ""
    out["labeler"] = ""
    out["updated_utc"] = ""
    out["notes"] = notes
    out["twirl_vet_sheet_name"] = twirl_vet_sheet_name or f"{safe_review_id}_twirl_twoap_{sheet_branch_name}.png"
    out["twirl_vet_sheet_pdf_name"] = ""
    _set_if_present_or_new(out, "period_d", period_d)
    _set_if_present_or_new(out, "t0_bjd", t0_bjd)
    _set_if_present_or_new(out, "duration_min", duration_min)
    out["candidate_key"] = ""
    out.pop("row_id", None)
    return out


def build_revisit_queue(
    training_table: Path,
    *,
    half_period_review_ids: tuple[str, ...] = DEFAULT_HALF_PERIOD_REVIEW_IDS,
    include_half_period_variants: bool = True,
    sheet_branch_name: str = DEFAULT_SHEET_BRANCH_NAME,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    table = pd.read_csv(training_table)
    source_kind = table.get("source_kind", pd.Series("", index=table.index)).fillna("").astype(str)
    human_label = table.get("human_label", pd.Series("", index=table.index)).fillna("").astype(str)
    review_id = table.get("review_id", pd.Series("", index=table.index)).fillna("").astype(str)
    selected = table[source_kind.eq("real_candidate") & human_label.eq("planet_like")].copy()
    selected = selected.sort_values(["tic", "review_id"], kind="stable").reset_index(drop=True)
    half_ids = set(half_period_review_ids)

    pair_rows: list[dict[str, Any]] = []
    other_rows: list[dict[str, Any]] = []
    selected_review_ids = set(selected.get("review_id", pd.Series(dtype=str)).fillna("").astype(str))
    missing_half_ids = sorted(half_ids - selected_review_ids)
    for _, source in selected.iterrows():
        period, t0, duration = _display_ephemeris(source)
        source_review_id = str(_clean(source.get("review_id")) or "")
        current = _base_revisit_row(
            source,
            variant="current",
            period_d=period,
            t0_bjd=t0,
            duration_min=duration,
            twirl_vet_sheet_name="",
            notes=(
                "Revisit real row previously labeled planet_like. "
                "Use wide_transit_like for broad/long-duration transit-like events that should not train "
                "the compact planet class."
            ),
            sheet_branch_name=sheet_branch_name,
        )
        if include_half_period_variants and source_review_id in half_ids:
            pair_rows.append(current)
            half_period = period / 2.0 if period and period > 0 else None
            period_text = f"{period:.8g} d" if period is not None else "unknown"
            pair_rows.append(
                _base_revisit_row(
                    source,
                    variant="half_period",
                    period_d=half_period,
                    t0_bjd=t0,
                    duration_min=duration,
                    twirl_vet_sheet_name="",
                    notes=(
                        "Half-period refold of a real row previously labeled planet_like. "
                        f"Original displayed period was {period_text}. "
                        "Use EB/PCEB if alternating/eclipse-like behavior is clearer; use wide_transit_like "
                        "for broad transit-like events; use Planet only if this remains a compact planet-like signal."
                    ),
                    sheet_branch_name=sheet_branch_name,
                )
            )
        else:
            other_rows.append(current)

    queue = pd.DataFrame(pair_rows + other_rows)

    summary = {
        "training_table": str(training_table),
        "n_real_planet_like_source_rows": int(len(selected)),
        "n_half_period_source_rows": int(sum(selected["review_id"].astype(str).isin(half_ids))),
        "n_revisit_rows": int(len(queue)),
        "n_current_fold_rows": int((queue.get("revisit_variant", pd.Series(dtype=str)) == "current").sum()),
        "n_half_period_fold_rows": int((queue.get("revisit_variant", pd.Series(dtype=str)) == "half_period").sum()),
        "half_period_review_ids": sorted(half_ids),
        "missing_half_period_review_ids": missing_half_ids,
        "include_half_period_variants": bool(include_half_period_variants),
        "label_policy_note": (
            "wide_transit_like is a separate raw human label; the main teacher target should continue "
            "to keep compact planet_like separate from wide/long-duration signals."
        ),
        "sheet_branch_name": sheet_branch_name,
    }
    return queue, summary


def _write_summary_md(summary: dict[str, Any], path: Path) -> None:
    lines = [
        "# S56 Real Planet-Like Revisit Queue",
        "",
        f"- Source rows previously labeled `planet_like` and real: `{summary['n_real_planet_like_source_rows']}`",
        f"- Revisit rows written: `{summary['n_revisit_rows']}`",
        f"- Current-fold rows: `{summary['n_current_fold_rows']}`",
        f"- Half-period refold rows: `{summary['n_half_period_fold_rows']}`",
        "",
        "Use `wide_transit_like` for broad/long-duration transit-like signals that should be preserved for audit but should not train the compact planet-like class.",
        "Use `eclipsing_binary_or_pceb` when the half-period fold makes the signal look EB/PCEB-like.",
        "The original human labels are not overwritten; new labels belong in `human_labels_revisit.csv`.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def _build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--training-table", type=Path, default=DEFAULT_TRAINING_TABLE)
    ap.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    ap.add_argument(
        "--half-period-review-id",
        action="append",
        default=None,
        help="Review id to duplicate with period_d/2. May be repeated.",
    )
    ap.add_argument(
        "--no-half-period-variants",
        action="store_true",
        help="Write only the current-fold revisit rows, without duplicated period/2 rows.",
    )
    ap.add_argument("--sheet-branch-name", default=DEFAULT_SHEET_BRANCH_NAME)
    return ap


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    half_ids = tuple(args.half_period_review_id or DEFAULT_HALF_PERIOD_REVIEW_IDS)
    queue, summary = build_revisit_queue(
        args.training_table,
        half_period_review_ids=half_ids,
        include_half_period_variants=not bool(args.no_half_period_variants),
        sheet_branch_name=str(args.sheet_branch_name),
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    queue_path = args.output_dir / "review_queue_real_planet_revisit.csv"
    labels_path = args.output_dir / "human_labels_revisit.csv"
    queue.to_csv(queue_path, index=False)
    if not labels_path.exists():
        pd.DataFrame(
            columns=["row_id", "candidate_key", "tic", "sector", "label", "label_source", "labeler", "notes", "updated_utc"]
        ).to_csv(labels_path, index=False)
    summary["review_queue"] = str(queue_path)
    summary["labels_out"] = str(labels_path)
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    _write_summary_md(summary, args.output_dir / "summary.md")
    print(f"[revisit-queue] wrote {len(queue)} rows: {queue_path}")
    print(f"[revisit-queue] labels: {labels_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
