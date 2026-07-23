#!/usr/bin/env python3
"""Verify queue identity, completion state, and local vet-sheet coverage."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import pandas as pd

from twirl.vetting.multisector_signal_review import standalone_app_candidate_key


def verify(bundle_dir: Path) -> dict[str, object]:
    bundle_dir = Path(bundle_dir)
    queue_path = bundle_dir / "review_queue_planet_eb.csv"
    labels_path = bundle_dir / "human_labels_final.csv"
    manifest_path = bundle_dir / "required_sheet_manifest.csv"
    sheet_root = bundle_dir / "vet_sheets"
    queue = pd.read_csv(queue_path, dtype=str, keep_default_na=False)
    labels = pd.read_csv(labels_path, dtype=str, keep_default_na=False)
    manifest = pd.read_csv(manifest_path, dtype=str, keep_default_na=False)
    required_queue = {
        "row_id",
        "candidate_key",
        "observation_candidate_key",
        "selected_source_uid",
        "twirl_vet_sheet_name",
    }
    missing_queue = sorted(required_queue - set(queue.columns))
    if missing_queue:
        raise KeyError(f"review queue is missing columns: {missing_queue}")
    if queue["row_id"].duplicated().any() or manifest["row_id"].duplicated().any():
        raise ValueError("queue or sheet manifest has duplicate row_id values")
    if (
        queue["candidate_key"].eq("").any()
        or queue["selected_source_uid"].eq("").any()
        or queue["candidate_key"].ne(queue["observation_candidate_key"]).any()
    ):
        raise ValueError(
            "review queue has a blank or inconsistent immutable candidate/source "
            "identity"
        )
    if set(queue["row_id"]) != set(manifest["row_id"]):
        raise ValueError("sheet manifest does not cover the exact queue")
    if labels["row_id"].duplicated().any():
        raise ValueError("final label file has duplicate row_id values")
    if not set(labels["row_id"]).issubset(set(queue["row_id"])):
        raise ValueError("final label file contains unknown row_id values")
    if not labels.empty:
        expected = queue.set_index("row_id").apply(
            standalone_app_candidate_key, axis=1
        )
        mismatch = labels["candidate_key"].ne(labels["row_id"].map(expected))
        if mismatch.any():
            raise ValueError("final label file candidate_key mismatch")
    missing: list[str] = []
    invalid_png: list[str] = []
    for name in manifest["twirl_vet_sheet_name"]:
        path = sheet_root / name
        if not path.is_file():
            missing.append(name)
            continue
        with path.open("rb") as handle:
            signature = handle.read(8)
        if path.stat().st_size < 8 or signature != b"\x89PNG\r\n\x1a\n":
            invalid_png.append(name)
    missing_path = bundle_dir / "missing_sheet_names.txt"
    missing_path.write_text("".join(f"{name}\n" for name in sorted(missing)))
    summary: dict[str, object] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "n_queue_rows": int(len(queue)),
        "n_reviewed_rows": int(len(labels)),
        "n_pending_rows": int(len(queue) - len(labels)),
        "n_required_sheets": int(len(manifest)),
        "n_present_sheets": int(len(manifest) - len(missing)),
        "n_missing_sheets": int(len(missing)),
        "n_invalid_png": int(len(invalid_png)),
        "missing_sheet_names": missing,
        "invalid_png_names": invalid_png,
        "review_complete": bool(len(labels) == len(queue)),
        "assets_complete": bool(not missing and not invalid_png),
    }
    (bundle_dir / "verification.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle-dir", type=Path, required=True)
    parser.add_argument(
        "--allow-incomplete",
        action="store_true",
        help="Return success while sheets or row reviews are still pending.",
    )
    args = parser.parse_args()
    summary = verify(args.bundle_dir)
    print(json.dumps(summary, indent=2, sort_keys=True))
    complete = bool(summary["assets_complete"] and summary["review_complete"])
    if not complete and not args.allow_incomplete:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
