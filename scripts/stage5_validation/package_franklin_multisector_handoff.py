#!/usr/bin/env python3
"""Validate and ZIP a Franklin multisector vetting handoff."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import zipfile

import pandas as pd


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--handoff-dir", type=Path, required=True)
    parser.add_argument("--expected-rows", type=int, required=True)
    parser.add_argument("--out-zip", type=Path)
    args = parser.parse_args()

    root = args.handoff_dir
    queue_paths = sorted(root.glob("franklin_review_queue_*_real.csv"))
    if len(queue_paths) != 1:
        raise ValueError(f"expected one public queue under {root}, found {queue_paths}")
    queue_path = queue_paths[0]
    required = (
        root / "franklin_labels_vetted.csv",
        root / "franklin_vetting_app.py",
        root / "run_franklin_vetting.sh",
        root / "README_Franklin_vetting.md",
        root / "summary.json",
        root / "reference_examples.csv",
        root / "reference_examples",
        root / "vet_sheets",
    )
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"handoff is missing required paths: {missing}")

    queue = pd.read_csv(queue_path, low_memory=False)
    if len(queue) != int(args.expected_rows):
        raise ValueError(
            f"queue rows differ: observed={len(queue)}, expected={args.expected_rows}"
        )
    if queue["tic"].nunique() != len(queue):
        raise ValueError("public queue repeats a TIC")
    if set(queue["source_kind"].fillna("").astype(str)) != {"real_candidate"}:
        raise ValueError("public queue is not real-only")
    if not pd.to_numeric(queue["rep_peak_rank"], errors="coerce").eq(1).all():
        raise ValueError("public queue includes a non-rank-1 ephemeris")
    sheet_names = queue["twirl_vet_sheet_name"].fillna("").astype(str)
    missing_sheets = [
        name for name in sheet_names if not name or not (root / "vet_sheets" / name).is_file()
    ]
    if missing_sheets:
        raise FileNotFoundError(
            f"handoff is missing {len(missing_sheets)} sheets; first={missing_sheets[:3]}"
        )
    if list((root / "vet_sheets").glob("*.pdf")):
        raise ValueError("handoff contains PDF vet sheets")
    forbidden = (
        "hidden_selection_provenance",
        "teacher_scores",
        "model_scores",
    )
    leaked = [
        str(path)
        for path in root.rglob("*")
        if path.is_file() and any(token in path.name for token in forbidden)
    ]
    if leaked:
        raise ValueError(f"handoff exposes hidden model provenance: {leaked}")

    include = sorted(path for path in root.rglob("*") if path.is_file())
    out_zip = args.out_zip or root.with_suffix(".zip")
    out_zip.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(
        out_zip,
        "w",
        compression=zipfile.ZIP_STORED,
        allowZip64=True,
    ) as archive:
        for index, path in enumerate(include, start=1):
            archive.write(path, root.name / path.relative_to(root))
            if index % 250 == 0:
                print(f"[franklin-package] added {index:,}/{len(include):,} files", flush=True)

    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "handoff_dir": str(root),
        "out_zip": str(out_zip),
        "zip_size_bytes": int(out_zip.stat().st_size),
        "zip_sha256": _sha256(out_zip),
        "n_queue_rows": int(len(queue)),
        "n_unique_tics": int(queue["tic"].nunique()),
        "sector_counts": {
            str(key): int(value)
            for key, value in queue["sector"].value_counts().sort_index().items()
        },
        "n_png_sheets": int(len(sheet_names)),
        "n_reference_pngs": int(len(list((root / "reference_examples").rglob("*.png")))),
        "n_included_files": int(len(include)),
        "compression": "stored (PNG assets are already compressed)",
    }
    manifest_path = out_zip.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
