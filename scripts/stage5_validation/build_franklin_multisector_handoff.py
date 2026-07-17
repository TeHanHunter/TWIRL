#!/usr/bin/env python3
"""Build a self-contained real-only multisector vetting handoff for Franklin."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import stat
from typing import Any

import pandas as pd


LABEL_HEADER: tuple[str, ...] = (
    "row_id",
    "candidate_key",
    "tic",
    "sector",
    "label",
    "label_source",
    "labeler",
    "notes",
    "period_factor",
    "period_status",
    "updated_utc",
)
HIDDEN_PREFIXES: tuple[str, ...] = (
    "p_",
    "std_p_",
    "member_",
    "selection_",
    "model_",
    "source_candidate",
    "active_learning_",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_empty_labels(path: Path) -> None:
    with path.open("w", newline="") as handle:
        csv.DictWriter(handle, fieldnames=LABEL_HEADER).writeheader()


def _parse_sector_counts(values: list[list[str]]) -> dict[int, int]:
    out: dict[int, int] = {}
    for sector_text, count_text in values:
        sector = int(sector_text)
        if sector in out:
            raise ValueError(f"duplicate expected count for Sector {sector}")
        out[sector] = int(count_text)
    return out


def _copy_sheet(source_root: Path, destination_root: Path, name: str) -> None:
    source = source_root / name
    if not source.is_file():
        raise FileNotFoundError(source)
    destination = destination_root / name
    if source.resolve() == destination.resolve():
        return
    shutil.copy2(source, destination)


def _write_launcher(out_dir: Path, *, queue_name: str, port: int) -> None:
    path = out_dir / "run_franklin_vetting.sh"
    path.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "cd \"$(dirname \"$0\")\"\n"
        "exec python3 franklin_vetting_app.py \\\n"
        f"  --queue {queue_name} \\\n"
        "  --labels-out franklin_labels_vetted.csv \\\n"
        "  --sheet-root vet_sheets \\\n"
        "  --labeler franklin \\\n"
        "  --host 127.0.0.1 \\\n"
        f"  --port {int(port)}\n"
    )
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def _write_readme(
    out_dir: Path,
    *,
    queue_name: str,
    n_rows: int,
    sector_counts: dict[int, int],
    port: int,
) -> None:
    counts = ", ".join(
        f"S{sector}: {count:,}" for sector, count in sorted(sector_counts.items())
    )
    text = f"""# TWIRL S57-S59 Enriched Vetting Handoff

This package contains {n_rows:,} real candidates selected for human review.
Sector allocation: {counts}.

The queue is enriched by a frozen teacher model, but model scores and
selection buckets are intentionally absent. Every displayed ephemeris is the
rank-1 `DET_FLUX_ADP_SML` BLS candidate. The sheets use only
`DET_FLUX_ADP_SML` and `DET_FLUX_ADP`.

## Start

From this directory:

```bash
python3 franklin_vetting_app.py --check-only \\
  --queue {queue_name} \\
  --labels-out franklin_labels_vetted.csv \\
  --sheet-root vet_sheets
./run_franklin_vetting.sh
```

Open `http://127.0.0.1:{int(port)}/`.

Every click is saved immediately to `franklin_labels_vetted.csv`. Return that
CSV; the PNGs do not need to be returned.

## Six Labels

- **Planet-like**: compact isolated transit-like event without a convincing
  secondary or binary morphology.
- **Eclipse/contact**: secondary eclipse, alternating depths, or convincing
  detached/contact-binary morphology.
- **Smooth variable**: continuous or sinusoidal modulation without discrete
  eclipses.
- **Systematic/artifact**: window-edge leakage, cadence aliases, detrending
  structure, or another recognizable artifact.
- **Flat/no signal**: no obvious useful event. This is not an ambiguity label.
- **Broad isolated dip**: broad primary-like event without a convincing
  secondary or continuous variability.

`Skip` is only for broken or unusable evidence.

Use Eclipse/contact instead of Broad isolated dip when a secondary, odd/even
difference, or contact-binary pattern is credible.

## Period Control

Leave the period at `P` unless the displayed fold is visibly at the wrong
harmonic. Choose `P/4`, `P/2`, `2P`, or `4P` when one is clearly better.
Choose `Unresolved` and leave a short note for another factor such as `3P`.
The label and period choice are saved atomically.

## Keyboard

- `1`: Planet-like
- `2`: Eclipse/contact
- `3`: Smooth variable
- `4`: Systematic/artifact
- `5`: Flat/no signal
- `6`: Broad isolated dip
- `0`: Skip
- left/right arrow: previous/next

Reference examples are under `reference_examples/`.
"""
    (out_dir / "README_Franklin_vetting.md").write_text(text)


def build_handoff(
    *,
    queue_path: Path,
    sheet_root: Path,
    out_dir: Path,
    app_source: Path,
    expected_sector_counts: dict[int, int],
    reference_root: Path | None,
    reference_csv: Path | None,
    port: int,
) -> dict[str, Any]:
    queue = pd.read_csv(queue_path, low_memory=False)
    if queue.empty:
        raise ValueError("Franklin queue is empty")
    required = {
        "row_id",
        "candidate_key",
        "tic",
        "sector",
        "source_kind",
        "rep_peak_rank",
        "twirl_vet_sheet_name",
        "twirl_vet_sheet_pdf_name",
    }
    missing = sorted(required - set(queue.columns))
    if missing:
        raise KeyError(f"Franklin queue is missing columns: {missing}")
    if queue["tic"].nunique() != len(queue):
        raise ValueError("Franklin queue repeats a TIC")
    if set(queue["source_kind"].fillna("").astype(str)) != {"real_candidate"}:
        raise ValueError("Franklin queue is not real-only")
    if not pd.to_numeric(queue["rep_peak_rank"], errors="coerce").eq(1).all():
        raise ValueError("Franklin queue includes a non-rank-1 ephemeris")
    if queue["twirl_vet_sheet_pdf_name"].fillna("").astype(str).ne("").any():
        raise ValueError("Franklin queue requests PDF sheets")
    exposed = [column for column in queue if column.startswith(HIDDEN_PREFIXES)]
    if exposed:
        raise ValueError(f"Franklin queue exposes hidden provenance: {exposed}")
    observed_sector_counts = {
        int(key): int(value)
        for key, value in queue["sector"].value_counts().sort_index().items()
    }
    if observed_sector_counts != expected_sector_counts:
        raise ValueError(
            f"sector counts differ: observed={observed_sector_counts}, "
            f"expected={expected_sector_counts}"
        )

    if out_dir.exists() and any(out_dir.iterdir()):
        raise FileExistsError(
            f"handoff destination must be absent or empty: {out_dir}"
        )
    out_dir.mkdir(parents=True, exist_ok=True)
    destination_sheets = out_dir / "vet_sheets"
    destination_sheets.mkdir(parents=True, exist_ok=True)
    names = queue["twirl_vet_sheet_name"].fillna("").astype(str).tolist()
    if not all(names):
        raise ValueError("Franklin queue contains a blank sheet name")
    if len(names) != len(set(names)):
        raise ValueError("Franklin queue contains duplicate sheet names")
    for index, name in enumerate(names, start=1):
        _copy_sheet(Path(sheet_root), destination_sheets, name)
        if index % 100 == 0:
            print(f"[franklin-handoff] copied {index:,}/{len(names):,} sheets", flush=True)

    count_token = (
        f"{len(queue) // 1000}k" if len(queue) % 1000 == 0 else str(len(queue))
    )
    queue_name = f"franklin_review_queue_{count_token}_real.csv"
    destination_queue = out_dir / queue_name
    if Path(queue_path).resolve() != destination_queue.resolve():
        queue.to_csv(destination_queue, index=False, float_format="%.15g")
    labels_path = out_dir / "franklin_labels_vetted.csv"
    _write_empty_labels(labels_path)
    shutil.copy2(app_source, out_dir / "franklin_vetting_app.py")
    _write_launcher(out_dir, queue_name=queue_name, port=port)
    _write_readme(
        out_dir,
        queue_name=queue_name,
        n_rows=len(queue),
        sector_counts=observed_sector_counts,
        port=port,
    )

    n_reference_pngs = 0
    if reference_root is not None:
        reference_destination = out_dir / "reference_examples"
        shutil.copytree(reference_root, reference_destination)
        n_reference_pngs = len(list(reference_destination.rglob("*.png")))
    if reference_csv is not None:
        shutil.copy2(reference_csv, out_dir / "reference_examples.csv")

    missing_sheets = [
        name for name in names if not (destination_sheets / name).is_file()
    ]
    if missing_sheets:
        raise FileNotFoundError(
            f"handoff is missing {len(missing_sheets)} sheets; first={missing_sheets[:3]}"
        )
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "n_rows": int(len(queue)),
        "n_unique_tics": int(queue["tic"].nunique()),
        "sector_counts": observed_sector_counts,
        "rank_policy": "ADP-small representative BLS peak rank 1 only",
        "source_kind_counts": {"real_candidate": int(len(queue))},
        "n_png_sheets": int(len(names)),
        "n_pdf_sheets": 0,
        "n_reference_pngs": int(n_reference_pngs),
        "queue_sha256": _sha256(destination_queue),
        "app_sha256": _sha256(out_dir / "franklin_vetting_app.py"),
        "scores_in_package": False,
        "selection_provenance_in_package": False,
        "outputs": {
            "queue": str(destination_queue),
            "labels": str(labels_path),
            "sheet_root": str(destination_sheets),
            "readme": str(out_dir / "README_Franklin_vetting.md"),
            "launcher": str(out_dir / "run_franklin_vetting.sh"),
        },
    }
    (out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queue", type=Path, required=True)
    parser.add_argument("--sheet-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument(
        "--app-source",
        type=Path,
        default=Path(__file__).with_name("franklin_vetting_app.py"),
    )
    parser.add_argument(
        "--expected-sector-count",
        nargs=2,
        action="append",
        metavar=("SECTOR", "COUNT"),
        required=True,
    )
    parser.add_argument("--reference-root", type=Path)
    parser.add_argument("--reference-csv", type=Path)
    parser.add_argument("--port", type=int, default=5003)
    args = parser.parse_args()
    summary = build_handoff(
        queue_path=args.queue,
        sheet_root=args.sheet_root,
        out_dir=args.out_dir,
        app_source=args.app_source,
        expected_sector_counts=_parse_sector_counts(args.expected_sector_count),
        reference_root=args.reference_root,
        reference_csv=args.reference_csv,
        port=args.port,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
