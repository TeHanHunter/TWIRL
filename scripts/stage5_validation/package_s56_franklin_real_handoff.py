#!/usr/bin/env python3
"""Create a self-contained ZIP for the Franklin S56 real-only vetting handoff."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import zipfile

import pandas as pd


DEFAULT_QUEUE_DIR = Path("reports/stage5_validation/s56_franklin_real5k_handoff")


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _require(path: Path, missing: list[str]) -> None:
    if not path.exists():
        missing.append(str(path))


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queue-dir", type=Path, default=DEFAULT_QUEUE_DIR)
    parser.add_argument("--expected-rows", type=int, default=5000)
    parser.add_argument("--out-zip", type=Path, default=None)
    return parser


def main() -> int:
    args = _build_arg_parser().parse_args()
    queue_dir = args.queue_dir
    out_zip = args.out_zip or queue_dir.with_name(f"{queue_dir.name}.zip")

    queue_csv = queue_dir / "franklin_review_queue_5k_real.csv"
    labels_csv = queue_dir / "franklin_labels_vetted.csv"
    app_py = queue_dir / "franklin_vetting_app.py"
    launcher = queue_dir / "run_franklin_vetting.sh"
    readme = queue_dir / "README_Franklin_vetting.md"
    summary_json = queue_dir / "summary.json"
    reference_csv = queue_dir / "reference_examples.csv"
    sheet_root = queue_dir / "vet_sheets"
    examples_root = queue_dir / "reference_examples"

    missing: list[str] = []
    for path in (queue_csv, labels_csv, app_py, launcher, readme, summary_json, reference_csv, sheet_root, examples_root):
        _require(path, missing)
    if missing:
        raise FileNotFoundError("missing required handoff files: " + ", ".join(missing))

    queue = pd.read_csv(queue_csv)
    if len(queue) != int(args.expected_rows):
        raise ValueError(f"queue rows {len(queue)} != expected {args.expected_rows}")
    source_counts = queue["source_kind"].fillna("").astype(str).value_counts().to_dict()
    if source_counts != {"real_candidate": int(args.expected_rows)}:
        raise ValueError(f"queue is not real-only: {source_counts}")

    missing_sheets: list[str] = []
    for value in queue["twirl_vet_sheet_name"].fillna("").astype(str):
        if not value:
            missing_sheets.append("<blank>")
        elif not (sheet_root / value).exists():
            missing_sheets.append(value)
    if missing_sheets:
        raise FileNotFoundError(f"missing {len(missing_sheets)} PNG sheets; first={missing_sheets[0]}")

    pdf_count = len(list(sheet_root.glob("*.pdf")))
    if pdf_count:
        raise ValueError(f"unexpected PDF sheets remain in package tree: {pdf_count}")

    png_count = len(list(sheet_root.glob("*.png")))
    ref_png_count = len(list(examples_root.rglob("*.png")))

    include_paths: list[Path] = [
        queue_csv,
        labels_csv,
        app_py,
        launcher,
        readme,
        summary_json,
        reference_csv,
    ]
    for optional in (
        queue_dir / "twirl_vet_metrics_current_adp.csv",
        queue_dir / "twirl_two_aperture_vet_summary_current_adp.json",
    ):
        if optional.exists():
            include_paths.append(optional)
    include_paths.extend(sorted(sheet_root.glob("*.png")))
    include_paths.extend(sorted(examples_root.rglob("*.png")))

    out_zip.parent.mkdir(parents=True, exist_ok=True)
    base = queue_dir.name
    if out_zip.exists():
        out_zip.unlink()
    with zipfile.ZipFile(out_zip, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=4) as zf:
        for path in include_paths:
            zf.write(path, Path(base) / path.relative_to(queue_dir))

    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "queue_dir": str(queue_dir),
        "out_zip": str(out_zip),
        "zip_size_bytes": out_zip.stat().st_size,
        "zip_sha256": _sha256(out_zip),
        "n_queue_rows": int(len(queue)),
        "source_kind_counts": source_counts,
        "n_png_sheets": int(png_count),
        "n_pdf_sheets_included": 0,
        "n_reference_pngs": int(ref_png_count),
        "included_files": int(len(include_paths)),
    }
    manifest_path = out_zip.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
