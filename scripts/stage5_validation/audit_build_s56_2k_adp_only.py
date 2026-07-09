#!/usr/bin/env python3
"""Audit the labeled S56 2k ephemerides and build the active ADP-only table."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.vetting.adp_only import (  # noqa: E402
    ADP_ONLY_APERTURES,
    build_adp_only_training_frame,
    canonical_det_flux_columns,
    classify_period_relation,
)
from twirl.vetting.recovery50_teacher import json_default  # noqa: E402


DEFAULT_ROOT = REPO_ROOT / "reports/stage5_validation/s56_recovery50_teacher_queue_2k"
DEFAULT_TRAINING_TABLE = (
    DEFAULT_ROOT / "human_training_table_revisit_harmonic/human_vetting_training_table.csv"
)
DEFAULT_METRICS = (
    REPO_ROOT / "reports/stage5_validation/s56_recovery50_teacher_queue/twirl_vet_metrics_real.csv",
    REPO_ROOT
    / "reports/stage5_validation/s56_recovery50_teacher_queue/"
    "twirl_vet_metrics_injected_fullphase_binmatch.csv",
    REPO_ROOT
    / "reports/stage5_validation/s56_recovery50_teacher_queue_next4k/twirl_vet_metrics_real.csv",
    REPO_ROOT
    / "reports/stage5_validation/s56_recovery50_teacher_queue_next4k/twirl_vet_metrics_injected.csv",
)
DEFAULT_OUT_DIR = DEFAULT_ROOT / "human_training_table_adp_only"


def _counts(series: pd.Series) -> dict[str, int]:
    return {
        str(key): int(value)
        for key, value in series.fillna("missing").astype(str).value_counts().sort_index().items()
    }


def _load_audit_metrics(paths: tuple[Path, ...]) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for path in paths:
        if not path.exists():
            continue
        frame = pd.read_csv(path, low_memory=False)
        if "review_id" not in frame:
            continue
        frame["metric_table"] = str(path)
        frame["supplement_kind"] = (
            "adp_primary"
            if "adp_period_d" in frame
            else "canonical_sml" if "sml_period_d" in frame else "missing"
        )
        frames.append(frame)
    if not frames:
        raise FileNotFoundError("no vet-metric tables were available")
    metrics = pd.concat(frames, ignore_index=True, sort=False)
    metrics["review_id"] = metrics["review_id"].fillna("").astype(str)
    return metrics.drop_duplicates("review_id", keep="last")


def _markdown(summary: dict[str, Any]) -> str:
    lines = [
        "# S56 Labeled 2k ADP-Only Audit",
        "",
        f"- rows: `{summary['n_rows']}`",
        f"- human display exactly matched ADP small: `{summary['display_vs_adp_small']['exact']}`",
        f"- human display missing ADP small: `{summary['display_vs_adp_small'].get('missing', 0)}`",
        f"- rows with ADP primary supplement: `{summary['supplement_kind_counts'].get('adp_primary', 0)}`",
        f"- rows with historical canonical supplement: `{summary['supplement_kind_counts'].get('canonical_sml', 0)}`",
        f"- active ADP-only rows: `{summary['adp_only']['n_adp_only_review_ok']}`",
        f"- excluded pending ADP-only review: `{summary['adp_only']['n_excluded']}`",
        "",
        "The original queue period frequently differs from the ADP-small BLS anchor, but the",
        "vet sheets used the ADP-small anchor for 1,999/2,000 rows. The lone missing row was",
        "already labeled skip. Historical canonical columns are removed from the active table.",
        "",
    ]
    return "\n".join(lines)


def build_audit(
    *,
    training_table: Path,
    metrics_tables: tuple[Path, ...],
    out_dir: Path,
) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    training = pd.read_csv(training_table, low_memory=False)
    metrics = _load_audit_metrics(metrics_tables)
    audit_columns = [
        column
        for column in (
            "review_id",
            "anchor_aperture",
            "anchor_period_d",
            "anchor_t0_bjd",
            "anchor_duration_min",
            "adp_sml_period_d",
            "adp_sml_t0_bjd",
            "adp_sml_duration_min",
            "adp_period_d",
            "sml_period_d",
            "supplement_kind",
            "metric_table",
        )
        if column in metrics
    ]
    audit = training.merge(
        metrics.loc[:, audit_columns],
        on="review_id",
        how="left",
        validate="one_to_one",
    ).copy()
    audit["display_vs_adp_small"] = classify_period_relation(
        audit["display_period_d"], audit["adp_sml_period_d"]
    )
    audit["legacy_queue_vs_adp_small"] = classify_period_relation(
        audit["queue_period_d"], audit["adp_sml_period_d"]
    )
    supplement_period = pd.to_numeric(audit.get("adp_period_d"), errors="coerce")
    canonical_supplement = pd.to_numeric(audit.get("sml_period_d"), errors="coerce")
    use_canonical = audit["supplement_kind"].fillna("").eq("canonical_sml")
    supplement_period.loc[use_canonical] = canonical_supplement.loc[use_canonical]
    audit["display_vs_supplement"] = classify_period_relation(
        audit["display_period_d"], supplement_period
    )

    active, active_summary = build_adp_only_training_frame(
        training,
        metrics_tables=metrics_tables,
    )
    active_path = out_dir / "human_vetting_training_table.csv"
    active.to_csv(active_path, index=False)
    audit_path = out_dir / "ephemeris_audit_rows.csv"
    audit.to_csv(audit_path, index=False)
    excluded = active.loc[~active["adp_only_review_ok"].astype(bool)].copy()
    excluded_path = out_dir / "excluded_pending_adp_review.csv"
    excluded.to_csv(excluded_path, index=False)

    existing_metadata: list[str] = []
    model_summary = DEFAULT_ROOT / "cnn_teacher_harmonic_raw_multiview/cnn_shape_plus_bls/summary.json"
    if model_summary.exists():
        payload = json.loads(model_summary.read_text())
        existing_metadata = [str(value) for value in payload.get("metadata_columns", [])]
    summary = {
        "training_table": str(training_table),
        "metrics_tables": [str(path) for path in metrics_tables],
        "out_dir": str(out_dir),
        "n_rows": int(len(audit)),
        "display_vs_adp_small": _counts(audit["display_vs_adp_small"]),
        "legacy_queue_vs_adp_small": _counts(audit["legacy_queue_vs_adp_small"]),
        "display_vs_supplement": _counts(audit["display_vs_supplement"]),
        "supplement_kind_counts": _counts(audit["supplement_kind"]),
        "display_source_counts": _counts(audit["display_ephemeris_source"]),
        "existing_combined_model_canonical_metadata": canonical_det_flux_columns(existing_metadata),
        "adp_only": active_summary,
        "active_apertures": list(ADP_ONLY_APERTURES),
        "injection_source_h5_counts": _counts(
            active.loc[
                active["source_kind"].fillna("").astype(str).eq("injection_recovery"),
                "source_h5",
            ]
        ),
        "outputs": {
            "active_training_table": str(active_path),
            "audit_rows": str(audit_path),
            "excluded_rows": str(excluded_path),
            "summary": str(out_dir / "summary.json"),
        },
    }
    (out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, default=json_default) + "\n"
    )
    (out_dir / "summary.md").write_text(_markdown(summary))
    return summary


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training-table", type=Path, default=DEFAULT_TRAINING_TABLE)
    parser.add_argument("--metrics-table", type=Path, action="append", default=list(DEFAULT_METRICS))
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    summary = build_audit(
        training_table=args.training_table,
        metrics_tables=tuple(args.metrics_table or ()),
        out_dir=args.out_dir,
    )
    print(json.dumps(summary, indent=2, sort_keys=True, default=json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
