#!/usr/bin/env python3
"""Build one pre-filled Planet-like/EB review queue for S56--S62."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shlex

import pandas as pd

from twirl.vetting.multisector_signal_review import (
    build_signal_rereview_queue,
    normalize_accepted_franklin_signals,
    normalize_browser_signal_rows,
    normalize_s56_adjudicated_signals,
)


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_S56 = (
    ROOT
    / "reports/stage5_validation/s56_label_adjudication_real343/"
    "adjudicated_training_table/human_vetting_training_table_adjudicated.csv"
)
DEFAULT_REVISIT_ROOT = (
    ROOT
    / "reports/stage5_validation/s56_s64_existing_teacher_enrichment/"
    "sector_0056/batch_00/compact_planet_revisit_current"
)
DEFAULT_ENRICHMENT_ROOT = (
    ROOT
    / "reports/stage5_validation/s56_s64_existing_teacher_enrichment/"
    "sector_0056/batch_00"
)
DEFAULT_FRANKLIN = (
    ROOT
    / "reports/stage5_validation/franklin_s57_s59_label_return_20260721/"
    "accepted_morphology_labels.csv"
)
DEFAULT_OUT = ROOT / "data_local/label_reviews/s56_s62_signal_rereview"

LABEL_HEADER = (
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


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_readme(out_dir: Path, *, n_rows: int, port: int) -> None:
    text = f"""# S56--S62 final signal review

This queue contains the union of every observation previously labeled
Planet-like or Eclipse/contact in the accepted inputs. It currently has
{n_rows:,} rows. Cross-sector observations remain separate. Exact duplicate
candidates are collapsed by the recorded precedence policy, while the private
source manifest preserves every contributing decision.

The buttons are pre-filled from the prior decision, but a row is not final
until it has been explicitly saved in this pass. Period selections from this
morphology pass remain audit metadata: they do not create new harmonic truth
or erase previously verified S56 harmonic supervision.

Put the exact PNG files named in `required_sheet_manifest.csv` under
`vet_sheets/`, then run:

```bash
./run_local_review.sh
```

Open `http://127.0.0.1:{port}/`. Completion requires `reviewed == count` in
`http://127.0.0.1:{port}/api/summary`.
"""
    (out_dir / "README.md").write_text(text)


def build(
    *,
    s56_adjudicated: Path,
    s56_enrichment_queue: Path,
    s56_enrichment_labels: Path,
    s56_revisit_queue: Path,
    s56_revisit_labels: Path,
    accepted_franklin: list[Path],
    out_dir: Path,
    port: int,
) -> dict[str, object]:
    source_paths = [
        Path(s56_adjudicated),
        Path(s56_enrichment_queue),
        Path(s56_enrichment_labels),
        Path(s56_revisit_queue),
        Path(s56_revisit_labels),
        *(Path(path) for path in accepted_franklin),
    ]
    missing = [str(path) for path in source_paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing signal review inputs: {missing}")
    normalized = [
        normalize_s56_adjudicated_signals(
            pd.read_csv(s56_adjudicated, low_memory=False)
        ),
        normalize_browser_signal_rows(
            pd.read_csv(s56_enrichment_queue, low_memory=False),
            pd.read_csv(s56_enrichment_labels, low_memory=False),
            source_batch_id="s56_enrichment_checkpoint_177",
            source_priority=150,
            require_complete=False,
        ),
        normalize_browser_signal_rows(
            pd.read_csv(s56_revisit_queue, low_memory=False),
            pd.read_csv(s56_revisit_labels, low_memory=False),
            source_batch_id="s56_compact_planet_revisit_407",
        ),
    ]
    normalized.extend(
        normalize_accepted_franklin_signals(
            pd.read_csv(path, low_memory=False)
        )
        for path in accepted_franklin
    )
    public, provenance, assets, summary = build_signal_rereview_queue(normalized)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "vet_sheets").mkdir(exist_ok=True)
    queue_path = out_dir / "review_queue_planet_eb.csv"
    provenance_path = out_dir / "source_label_provenance.csv"
    assets_path = out_dir / "required_sheet_manifest.csv"
    labels_path = out_dir / "human_labels_final.csv"
    if labels_path.exists():
        prior_labels = pd.read_csv(labels_path, dtype=str, keep_default_na=False)
        if not prior_labels.empty:
            raise RuntimeError(
                "refusing to rebuild a signal queue after row-level final review "
                f"started: {labels_path}"
            )
    public.to_csv(queue_path, index=False, float_format="%.15g")
    provenance.to_csv(provenance_path, index=False, float_format="%.15g")
    assets.to_csv(assets_path, index=False, float_format="%.15g")
    sheet_lists: dict[str, dict[str, str]] = {}
    all_names_path = out_dir / "required_sheet_names_all.txt"
    all_names = sorted(
        {
            str(value).strip()
            for value in assets["twirl_vet_sheet_name"]
            if str(value).strip()
        }
    )
    all_names_path.write_text("".join(f"{name}\n" for name in all_names))
    for sector, sector_rows in assets.groupby("sector", sort=True):
        names_path = out_dir / f"required_sheet_names_s{int(sector):04d}.txt"
        names = sorted(
            {
                str(value).strip()
                for value in sector_rows["twirl_vet_sheet_name"]
                if str(value).strip()
            }
        )
        names_path.write_text("".join(f"{name}\n" for name in names))
        sheet_lists[str(int(sector))] = {
            "path": str(names_path),
            "sha256": _sha256(names_path),
        }
    sheet_lists["all"] = {
        "path": str(all_names_path),
        "sha256": _sha256(all_names_path),
    }
    if not labels_path.exists():
        pd.DataFrame(columns=LABEL_HEADER).to_csv(labels_path, index=False)

    app_path = ROOT / "scripts/stage5_validation/franklin_vetting_app.py"
    command = (
        f"#!/usr/bin/env bash\nset -euo pipefail\n"
        f"python {shlex.quote(str(app_path))} \\\n"
        f"  --queue {shlex.quote(str(queue_path))} \\\n"
        f"  --labels-out {shlex.quote(str(labels_path))} \\\n"
        f"  --sheet-root {shlex.quote(str(out_dir / 'vet_sheets'))} \\\n"
        f"  --exact-sheets-only \\\n"
        f"  --labeler tehan --host 127.0.0.1 --port {int(port)}\n"
    )
    launcher = out_dir / "run_local_review.sh"
    launcher.write_text(command)
    launcher.chmod(0o755)
    _write_readme(out_dir, n_rows=len(public), port=int(port))

    summary.update(
        {
            "inputs": {
                str(path): {"sha256": _sha256(path)} for path in source_paths
            },
            "outputs": {
                "review_queue": {
                    "path": str(queue_path),
                    "sha256": _sha256(queue_path),
                },
                "source_label_provenance": {
                    "path": str(provenance_path),
                    "sha256": _sha256(provenance_path),
                },
                "required_sheet_manifest": {
                    "path": str(assets_path),
                    "sha256": _sha256(assets_path),
                },
                "required_sheet_name_lists": sheet_lists,
                "human_labels_final": {"path": str(labels_path)},
            },
            "review_port": int(port),
        }
    )
    summary_path = out_dir / "summary.json"
    summary["outputs"]["summary"] = {"path": str(summary_path)}  # type: ignore[index]
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--s56-adjudicated", type=Path, default=DEFAULT_S56)
    parser.add_argument(
        "--s56-enrichment-queue",
        type=Path,
        default=DEFAULT_ENRICHMENT_ROOT / "review_queue_1k.csv",
    )
    parser.add_argument(
        "--s56-enrichment-labels",
        type=Path,
        default=DEFAULT_ENRICHMENT_ROOT / "human_labels_vetted.csv",
    )
    parser.add_argument(
        "--s56-revisit-queue",
        type=Path,
        default=DEFAULT_REVISIT_ROOT / "review_queue_compact_planet.csv",
    )
    parser.add_argument(
        "--s56-revisit-labels",
        type=Path,
        default=DEFAULT_REVISIT_ROOT / "human_labels_revisit.csv",
    )
    parser.add_argument(
        "--accepted-franklin",
        type=Path,
        action="append",
        default=None,
        help="Accepted morphology table; repeat for S57--S59 and S60--S62.",
    )
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--port", type=int, default=5014)
    args = parser.parse_args()
    accepted = args.accepted_franklin or [DEFAULT_FRANKLIN]
    summary = build(
        s56_adjudicated=args.s56_adjudicated,
        s56_enrichment_queue=args.s56_enrichment_queue,
        s56_enrichment_labels=args.s56_enrichment_labels,
        s56_revisit_queue=args.s56_revisit_queue,
        s56_revisit_labels=args.s56_revisit_labels,
        accepted_franklin=accepted,
        out_dir=args.out_dir,
        port=args.port,
    )
    summary["command_completed_utc"] = datetime.now(timezone.utc).isoformat()
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
