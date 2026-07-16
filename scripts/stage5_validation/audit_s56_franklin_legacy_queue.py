#!/usr/bin/env python3
"""Snapshot and quarantine the legacy canonical-flux Franklin S56 queue."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.vetting.teacher_active_learning import run_legacy_franklin_audit  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--queue",
        type=Path,
        default=ROOT / "reports/stage5_validation/s56_franklin_real5k_handoff/franklin_review_queue_5k_real.csv",
    )
    parser.add_argument(
        "--labels",
        type=Path,
        default=ROOT / "reports/stage5_validation/s56_franklin_real5k_handoff/franklin_labels_vetted.csv",
    )
    parser.add_argument("--adp-peaks", type=Path, required=True)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "reports/stage5_validation/s56_franklin_real5k_handoff/legacy_adp_audit",
    )
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    snapshot = None
    snapshot_sha256 = None
    if args.labels.exists():
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        snapshot_dir = args.out_dir / "label_snapshots"
        snapshot_dir.mkdir(parents=True, exist_ok=True)
        snapshot = snapshot_dir / f"franklin_labels_vetted_{timestamp}.csv"
        shutil.copy2(args.labels, snapshot)
        snapshot_sha256 = hashlib.sha256(snapshot.read_bytes()).hexdigest()
    summary = run_legacy_franklin_audit(
        queue_path=args.queue,
        adp_peaks_path=args.adp_peaks,
        labels_path=args.labels,
        out_dir=args.out_dir,
    )
    summary["label_snapshot"] = str(snapshot) if snapshot is not None else None
    summary["label_snapshot_sha256"] = snapshot_sha256
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
