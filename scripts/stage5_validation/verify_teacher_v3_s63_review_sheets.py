#!/usr/bin/env python3
"""Fail-closed coverage verification for the frozen S63 reviewer sheets."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import sys
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.vetting.s63_preprocessing import (  # noqa: E402
    file_sha256,
    validate_producer_git_sha,
    validate_sha256,
)


CONTRACT_VERSION = "twirl_teacher_v3_s63_review_sheet_coverage_v1"
APERTURES = ("DET_FLUX_ADP_SML", "DET_FLUX_ADP")
QUALITY_COUNT_COLUMNS = (
    "quality_n_cad_total",
    "quality_n_cad_internal_bad",
    "quality_n_cad_external_bad",
    "quality_n_cad_external_only_bad",
    "quality_n_cad_authority_excluded",
    "quality_n_cad_effective_bad",
)


def _load_json(path: Path, *, label: str) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{label} must contain one JSON object")
    return payload


def _exact_float(left: pd.Series, right: pd.Series, *, label: str) -> None:
    observed = pd.to_numeric(left, errors="raise").to_numpy(dtype=float)
    expected = pd.to_numeric(right, errors="raise").to_numpy(dtype=float)
    if not np.isfinite(observed).all() or not np.isfinite(expected).all():
        raise ValueError(f"{label} must be finite")
    if not np.allclose(observed, expected, rtol=0.0, atol=1.0e-12):
        raise ValueError(f"rendered {label} differs from the frozen queue")


def verify_review_sheets(
    *,
    queue_csv: Path,
    public_summary_path: Path,
    bundle_completion_marker: Path,
    metrics_csv: Path,
    render_summary_path: Path,
    sheet_dir: Path,
    compact_h5: Path,
    cadence_reference: Path,
    cadence_reference_manifest: Path,
    expected_rows: int,
    producer_git_sha: str,
    launch_manifest_sha256: str,
) -> dict[str, Any]:
    queue_csv = Path(queue_csv).resolve(strict=True)
    public_summary_path = Path(public_summary_path).resolve(strict=True)
    bundle_completion_marker = Path(bundle_completion_marker).resolve(strict=True)
    metrics_csv = Path(metrics_csv).resolve(strict=True)
    render_summary_path = Path(render_summary_path).resolve(strict=True)
    sheet_dir = Path(sheet_dir).resolve(strict=True)
    compact_h5 = Path(compact_h5).resolve(strict=True)
    cadence_reference = Path(cadence_reference).resolve(strict=True)
    cadence_reference_manifest = Path(cadence_reference_manifest).resolve(
        strict=True
    )
    producer = validate_producer_git_sha(producer_git_sha)
    launch_sha = validate_sha256(
        launch_manifest_sha256, label="launch-manifest SHA-256"
    )
    if queue_csv.name != "review_queue_1100.csv":
        raise ValueError("S63 public queue filename must be review_queue_1100.csv")
    if public_summary_path.name != "public_summary.json":
        raise ValueError("S63 public summary filename must be public_summary.json")

    queue = pd.read_csv(queue_csv, dtype=str, keep_default_na=False)
    metrics = pd.read_csv(metrics_csv, dtype=str, keep_default_na=False)
    if len(queue) != int(expected_rows) or len(metrics) != int(expected_rows):
        raise ValueError("queue and metrics must both have the exact frozen row count")
    required_queue = {
        "review_id",
        "tic",
        "sector",
        "period_d",
        "t0_bjd",
        "duration_min",
        "twirl_vet_sheet_name",
        "twirl_vet_sheet_pdf_name",
    }
    required_metrics = {
        *required_queue - {"period_d", "t0_bjd", "duration_min"},
        "twirl_vet_status",
        "anchor_aperture",
        "anchor_period_d",
        "anchor_t0_bjd",
        "anchor_duration_min",
        "cadence_reference_sha256",
        "cadence_reference_manifest_sha256",
        "orbitid_policy",
        *QUALITY_COUNT_COLUMNS,
    }
    missing_queue = sorted(required_queue - set(queue))
    missing_metrics = sorted(required_metrics - set(metrics))
    if missing_queue or missing_metrics:
        raise ValueError(
            f"review coverage inputs lack columns: queue={missing_queue}, "
            f"metrics={missing_metrics}"
        )
    if queue["review_id"].eq("").any() or queue["review_id"].duplicated().any():
        raise ValueError("queue review IDs must be nonempty and unique")
    if metrics["review_id"].eq("").any() or metrics["review_id"].duplicated().any():
        raise ValueError("metrics review IDs must be nonempty and unique")
    expected_names = queue["review_id"] + "_twirl_twoap_current_adp.png"
    if not queue["twirl_vet_sheet_name"].equals(expected_names):
        raise ValueError("queue sheet filenames are not the exact current-ADP names")
    if queue["twirl_vet_sheet_name"].duplicated().any():
        raise ValueError("queue declares duplicate sheet filenames")
    if queue["twirl_vet_sheet_pdf_name"].ne("").any():
        raise ValueError("S63 queue must declare no PDF sheets")

    merged = queue.merge(
        metrics,
        on="review_id",
        how="outer",
        validate="one_to_one",
        suffixes=("_queue", "_metric"),
        indicator=True,
        sort=False,
    )
    if set(merged["_merge"]) != {"both"}:
        raise ValueError("each frozen queue row must have exactly one metrics row")
    for field in ("tic", "sector", "twirl_vet_sheet_name"):
        if not merged[f"{field}_queue"].equals(merged[f"{field}_metric"]):
            raise ValueError(f"metrics {field} differs from the queue")
    if merged["twirl_vet_sheet_pdf_name_metric"].ne("").any():
        raise ValueError("render metrics unexpectedly declare PDFs")
    if set(merged["twirl_vet_status"]) != {"ok"}:
        raise ValueError("every S63 reviewer sheet must render with status=ok")
    if set(merged["anchor_aperture"]) != {APERTURES[0]}:
        raise ValueError("review sheets are not anchored on ADP-small")
    period_column = "period_d_queue" if "period_d_queue" in merged else "period_d"
    epoch_column = "t0_bjd_queue" if "t0_bjd_queue" in merged else "t0_bjd"
    duration_column = (
        "duration_min_queue" if "duration_min_queue" in merged else "duration_min"
    )
    _exact_float(merged["anchor_period_d"], merged[period_column], label="period")
    _exact_float(merged["anchor_t0_bjd"], merged[epoch_column], label="epoch")
    _exact_float(
        merged["anchor_duration_min"],
        merged[duration_column],
        label="duration",
    )
    cadence_sha = file_sha256(cadence_reference)
    cadence_manifest_sha = file_sha256(cadence_reference_manifest)
    if set(merged["cadence_reference_sha256"]) != {cadence_sha} or set(
        merged["cadence_reference_manifest_sha256"]
    ) != {cadence_manifest_sha}:
        raise ValueError("sheet metrics do not bind the exact cadence authority")
    if set(merged["orbitid_policy"]) != {"strict"}:
        raise ValueError("sheet rendering did not use strict orbit IDs")
    for column in QUALITY_COUNT_COLUMNS:
        values = pd.to_numeric(merged[column], errors="raise")
        if (values < 0).any() or (values != np.floor(values)).any():
            raise ValueError(f"invalid external-quality audit column {column}")

    declared_pngs = set(queue["twirl_vet_sheet_name"])
    observed_pngs = {
        path.name for path in sheet_dir.iterdir() if path.is_file() and path.suffix == ".png"
    }
    if observed_pngs != declared_pngs:
        raise ValueError("sheet directory does not exactly equal the declared PNG set")
    if any((sheet_dir / name).stat().st_size <= 0 for name in declared_pngs):
        raise ValueError("one or more declared PNG sheets is empty")
    if any(sheet_dir.rglob("*.pdf")):
        raise ValueError("S63 reviewer render unexpectedly contains PDFs")

    public_summary = _load_json(public_summary_path, label="public queue summary")
    if public_summary.get("schema_version") != (
        "twirl_teacher_v3_s63_public_review_queue_v1"
    ) or int(public_summary.get("n_queue_rows", -1)) != int(expected_rows):
        raise ValueError("public queue summary contract/count mismatch")
    if validate_producer_git_sha(public_summary.get("producer_git_sha")) != producer:
        raise ValueError("public queue producer Git SHA differs from render")
    queue_declaration = public_summary.get("outputs", {}).get("review_queue", {})
    if validate_sha256(
        queue_declaration.get("sha256"), label="public queue summary queue SHA-256"
    ) != file_sha256(queue_csv):
        raise ValueError("public summary does not bind the exact queue")
    marker = _load_json(bundle_completion_marker, label="queue bundle marker")
    for field, expected in {
        "schema_version": "twirl_teacher_v3_s63_queue_bundle_complete_v1",
        "sector": 63,
        "passed": True,
        "publication_complete": True,
        "visibility": "reviewer_safe_unpublished",
        "private_completion_marker_required": True,
    }.items():
        if marker.get(field) != expected:
            raise ValueError(f"queue bundle marker has incompatible {field}")
    marker_artifacts = marker.get("artifacts")
    if not isinstance(marker_artifacts, dict):
        raise ValueError("queue bundle marker lacks artifact records")
    if set(marker_artifacts) != {"public_review_queue", "public_summary"}:
        raise ValueError("public bundle marker contains non-public artifact records")
    if re.fullmatch(r"[0-9a-f]{32}", str(marker.get("bundle_id", ""))) is None:
        raise ValueError("public queue bundle marker has invalid bundle_id")
    for name, path in (
        ("public_review_queue", queue_csv),
        ("public_summary", public_summary_path),
    ):
        declaration = marker_artifacts.get(name)
        if not isinstance(declaration, dict):
            raise ValueError(f"queue bundle marker lacks {name}")
        if declaration.get("filename") != path.name or validate_sha256(
            declaration.get("sha256"), label=f"queue bundle {name} SHA-256"
        ) != file_sha256(path):
            raise ValueError(f"queue bundle marker does not bind {name}")

    render_summary = _load_json(render_summary_path, label="render summary")
    required_summary = {
        "branch_name": "current_adp",
        "apertures": list(APERTURES),
        "n_rows": int(expected_rows),
        "use_row_ephemeris": True,
        "write_pdf": False,
        "orbitid_policy": "strict",
        "producer_git_sha": producer,
        "launch_manifest_sha256": launch_sha,
        "queue_csv_sha256": file_sha256(queue_csv),
        "lc_export_h5_sha256": file_sha256(compact_h5),
    }
    for field, expected in required_summary.items():
        if render_summary.get(field) != expected:
            raise ValueError(f"render summary {field} differs from {expected!r}")
    external = render_summary.get("external_quality")
    if not isinstance(external, dict) or external.get(
        "cadence_reference_table_sha256"
    ) != cadence_sha or external.get(
        "cadence_reference_manifest_sha256"
    ) != cadence_manifest_sha:
        raise ValueError("render summary does not bind the cadence authority pair")
    if render_summary.get("status_counts") != {"ok": int(expected_rows)}:
        raise ValueError("render summary does not declare complete successful coverage")

    return {
        "contract_version": CONTRACT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "producer_git_sha": producer,
        "launch_manifest_sha256": launch_sha,
        "queue_sha256": file_sha256(queue_csv),
        "public_summary_sha256": file_sha256(public_summary_path),
        "queue_bundle_completion_marker_sha256": file_sha256(
            bundle_completion_marker
        ),
        "compact_h5_sha256": file_sha256(compact_h5),
        "cadence_reference_sha256": cadence_sha,
        "cadence_reference_manifest_sha256": cadence_manifest_sha,
        "metrics_sha256": file_sha256(metrics_csv),
        "render_summary_sha256": file_sha256(render_summary_path),
        "n_queue_rows": int(expected_rows),
        "n_metrics_rows": int(expected_rows),
        "n_unique_pngs": len(observed_pngs),
        "all_rows_exactly_once": True,
        "all_row_ephemerides_exact": True,
        "external_quality_identical_to_bls_authority": True,
        "pdfs_written": False,
        "passed": True,
    }


def _publish_json(payload: dict[str, Any], out: Path) -> None:
    out = Path(out).resolve()
    pending = out.with_suffix(out.suffix + ".pending")
    if out.exists() or pending.exists():
        raise FileExistsError(f"refusing review-verification overwrite: {out}")
    out.parent.mkdir(parents=True, exist_ok=True)
    pending.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    os.replace(pending, out)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queue-csv", type=Path, required=True)
    parser.add_argument("--public-summary", type=Path, required=True)
    parser.add_argument("--bundle-completion-marker", type=Path, required=True)
    parser.add_argument("--metrics-csv", type=Path, required=True)
    parser.add_argument("--render-summary", type=Path, required=True)
    parser.add_argument("--sheet-dir", type=Path, required=True)
    parser.add_argument("--compact-h5", type=Path, required=True)
    parser.add_argument("--cadence-reference", type=Path, required=True)
    parser.add_argument("--cadence-reference-manifest", type=Path, required=True)
    parser.add_argument("--expected-rows", type=int, default=1100)
    parser.add_argument("--producer-git-sha", required=True)
    parser.add_argument("--launch-manifest-sha256", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = verify_review_sheets(
        queue_csv=args.queue_csv,
        public_summary_path=args.public_summary,
        bundle_completion_marker=args.bundle_completion_marker,
        metrics_csv=args.metrics_csv,
        render_summary_path=args.render_summary,
        sheet_dir=args.sheet_dir,
        compact_h5=args.compact_h5,
        cadence_reference=args.cadence_reference,
        cadence_reference_manifest=args.cadence_reference_manifest,
        expected_rows=args.expected_rows,
        producer_git_sha=args.producer_git_sha,
        launch_manifest_sha256=args.launch_manifest_sha256,
    )
    _publish_json(result, args.out)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
