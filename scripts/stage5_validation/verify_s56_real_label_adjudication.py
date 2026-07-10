#!/usr/bin/env python3
"""Verify the fixed S56 real-label adjudication queue and rendered sheets."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.vetting.adjudication import (  # noqa: E402
    ADP_APERTURES,
    HIDDEN_PUBLIC_SUBSTRINGS,
    join_browser_labels,
    verify_queue_contract,
)
from twirl.vetting.label_io import latest_label_records, normalize_review_queue, validate_label_records  # noqa: E402


DEFAULT_ROOT = REPO_ROOT / "reports/stage5_validation/s56_label_adjudication_real343"
DEFAULT_PILOT_ROOT = REPO_ROOT / "reports/stage5_validation/s56_eb_miner_adp_only/review_queue_eb_priority_pilot100"


def verify(root: Path, *, require_sheets: bool) -> dict:
    queue_path = root / "review_queue_real343.csv"
    manifest_path = root / "adjudication_manifest_private.csv"
    labels_path = root / "human_labels_adjudicated.csv"
    sheets_dir = root / "twirl_vet_sheets"
    failures: list[str] = []
    queue = pd.read_csv(queue_path)
    manifest = pd.read_csv(manifest_path)
    contract = verify_queue_contract(queue, manifest)
    failures.extend(contract["failures"])

    exposed = [
        column
        for column in queue.columns
        if any(token in column.lower() for token in HIDDEN_PUBLIC_SUBSTRINGS)
    ]
    if exposed:
        failures.append(f"public queue exposes hidden columns: {exposed}")
    if set(queue["source_kind"].fillna("").astype(str)) != {"real_candidate"}:
        failures.append("public queue is not real-only")
    expected_apertures = ",".join(ADP_APERTURES)
    if set(queue["tensor_apertures"].fillna("").astype(str)) != {expected_apertures}:
        failures.append("public queue does not carry the exact ADP-only aperture signature")
    deprecated = manifest["source_cohort"].fillna("").astype(str).eq("deprecated_eb_signal")
    deprecated_sources = manifest.loc[deprecated].drop_duplicates("source_uid")
    if len(deprecated_sources) != 12:
        failures.append(f"deprecated source count={len(deprecated_sources)} expected=12")
    if not deprecated_sources["ephemeris_source"].fillna("").astype(str).eq(
        "current_adp_sml_bls_top1"
    ).all():
        failures.append("one or more deprecated rows retain a stale ephemeris")

    normalized = normalize_review_queue(queue)
    labels = latest_label_records(labels_path)
    try:
        validate_label_records(normalized, labels)
    except ValueError as exc:
        failures.append(str(exc))

    pilot_joined = join_browser_labels(
        DEFAULT_PILOT_ROOT / "review_queue_eb_priority_100.csv",
        DEFAULT_PILOT_ROOT / "human_labels_vetted.csv",
    )
    if len(pilot_joined) != 100 or pilot_joined["browser_label"].fillna("").astype(str).eq("").any():
        failures.append("completed EB pilot did not recover all 100 normalized labels")

    pngs = sorted(sheets_dir.glob("*.png")) if sheets_dir.exists() else []
    pdfs = sorted(sheets_dir.glob("*.pdf")) if sheets_dir.exists() else []
    expected_names = set(queue["twirl_vet_sheet_name"].fillna("").astype(str))
    actual_names = {path.name for path in pngs}
    if require_sheets and actual_names != expected_names:
        failures.append(
            f"sheet name mismatch: expected={len(expected_names)} actual={len(actual_names)} "
            f"missing={len(expected_names - actual_names)} extra={len(actual_names - expected_names)}"
        )
    if pdfs:
        failures.append(f"found {len(pdfs)} PDFs; adjudication sheets must be PNG-only")

    result = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "passed": not failures,
        "failures": failures,
        "contract": contract,
        "n_labels": int(labels["label"].fillna("").astype(str).ne("").sum()),
        "n_png": int(len(pngs)),
        "n_pdf": int(len(pdfs)),
        "require_sheets": bool(require_sheets),
    }
    (root / "verification.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--require-sheets", action="store_true")
    args = parser.parse_args(argv)
    result = verify(args.root, require_sheets=args.require_sheets)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
