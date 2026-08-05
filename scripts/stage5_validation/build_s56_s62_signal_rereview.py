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
DEFAULT_FRANKLIN_S57_S59 = (
    ROOT
    / "reports/stage5_validation/franklin_s57_s59_label_return_20260721/"
    "accepted_morphology_labels.csv"
)
DEFAULT_FRANKLIN_S60_S62 = (
    ROOT
    / "reports/stage5_validation/franklin_s60_s62_label_return_20260724/"
    "accepted_morphology_labels.csv"
)
DEFAULT_FRANKLIN = (
    DEFAULT_FRANKLIN_S57_S59,
    DEFAULT_FRANKLIN_S60_S62,
)
DEFAULT_OUT = ROOT / "data_local/label_reviews/s56_s62_signal_rereview"
SHEET_SET_VERSION = "s56_s62_a2v1_current_adp_v1"
SHEET_BRANCH_NAME = "current_adp"
SHEET_RENDERER_VERSION = "S56-ADP-HV2"
SHEET_APERTURES = ("DET_FLUX_ADP_SML", "DET_FLUX_ADP")
SHEET_HARMONIC_FACTORS = (0.25, 0.5, 1.0, 2.0, 4.0)
SHEET_N_PERIODS = 20_000
SHEET_N_PEAKS = 10
SHEET_LC_EXPORT_H5_BY_SECTOR = {
    56: (
        "/pdo/users/tehan/twirl_stage5/s56_A2v1_teacher_search_v1/"
        "inputs/s56_A2v1_adp_pair.h5"
    ),
    57: (
        "/pdo/users/tehan/twirl_stage5/s57_A2v1_teacher_search_v1/"
        "inputs/s57_A2v1_adp_pair.h5"
    ),
    58: (
        "/pdo/users/tehan/twirl_stage5/franklin_s56_s59_20260717/"
        "sector_0058/inputs/s58_A2v1_adp_pair.h5"
    ),
    59: (
        "/pdo/users/tehan/twirl_stage5/franklin_s56_s59_20260717/"
        "sector_0059/inputs/s59_A2v1_adp_pair.h5"
    ),
    60: (
        "/pdo/users/tehan/twirl_stage5/s60_A2v1_teacher_search_v1/"
        "inputs/s60_A2v1_adp_pair.h5"
    ),
    61: (
        "/pdo/users/tehan/twirl_stage5/s61_A2v1_teacher_search_v1/"
        "inputs/s61_A2v1_adp_pair.h5"
    ),
    62: (
        "/pdo/users/tehan/twirl_stage5/s62_A2v1_teacher_search_v1/"
        "inputs/s62_A2v1_adp_pair.h5"
    ),
}
SHEET_ROOT_NAME = f"vet_sheets_{SHEET_SET_VERSION}"
SHEET_PROVENANCE_ROOT_NAME = f"render_provenance_{SHEET_SET_VERSION}"

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


def _latest_sheet_name(review_id: object) -> str:
    safe = str(review_id).replace("/", "_").replace(":", "_")
    return f"{safe}_twirl_twoap_{SHEET_BRANCH_NAME}.png"


def _write_readme(out_dir: Path, *, n_rows: int, port: int) -> None:
    text = f"""# S56--S62 final signal review

This queue contains the union of every observation previously labeled
Planet-like or Eclipse/contact in the accepted inputs. It currently has
{n_rows:,} rows. Cross-sector observations remain separate. Exact duplicate
candidates are collapsed by the recorded precedence policy, while the private
source manifest preserves every contributing decision.

The standard TeHan vetter shows the prior decision as a suggestion, but a row
is not final until it has been explicitly saved in this pass. All Planet-like
suggestions appear first, followed by all Eclipse/contact suggestions. Period
selections from this morphology pass remain audit metadata: they do not create
new harmonic truth or erase previously verified S56 harmonic supervision.

This review is locked to the uniform `{SHEET_SET_VERSION}` sheet set:
`DET_FLUX_ADP_SML + DET_FLUX_ADP`, the frozen row ephemeris, and the current
`{SHEET_RENDERER_VERSION}` renderer with `{SHEET_N_PERIODS:,}` trial periods
and `{SHEET_N_PEAKS}` retained peaks. Put the exact PNG files named in
`required_sheet_manifest.csv` under `{SHEET_ROOT_NAME}/` and the renderer
metrics/summaries under `{SHEET_PROVENANCE_ROOT_NAME}/`, then run:

```bash
./run_local_review.sh
```

Open `http://127.0.0.1:{port}/`. The launcher uses the TeHan website with
exact-sheet matching: it cannot fall back to an older TIC-matched image.
Completion requires `reviewed == count` in
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
    expected_sectors = set(range(56, 63))
    observed_sectors = set(
        pd.to_numeric(public["sector"], errors="raise").astype(int)
    )
    if observed_sectors != expected_sectors:
        raise ValueError(
            "S56--S62 final review requires all seven sectors; "
            f"observed={sorted(observed_sectors)}"
        )
    latest_names = public["review_id"].map(_latest_sheet_name)
    public["twirl_vet_sheet_name"] = latest_names
    public["vet_sheet_set_version"] = SHEET_SET_VERSION
    latest_name_by_row = dict(zip(public["row_id"], latest_names))
    assets["twirl_vet_sheet_name"] = assets["row_id"].map(latest_name_by_row)
    assets["sheet_set_version"] = SHEET_SET_VERSION
    assets["branch_name"] = SHEET_BRANCH_NAME
    assets["renderer_version"] = SHEET_RENDERER_VERSION
    assets["apertures"] = ",".join(SHEET_APERTURES)
    assets["use_row_ephemeris"] = True
    assets["write_pdf"] = False
    assets["harmonic_factors"] = ",".join(
        format(value, "g") for value in SHEET_HARMONIC_FACTORS
    )

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sheet_root = out_dir / SHEET_ROOT_NAME
    provenance_root = out_dir / SHEET_PROVENANCE_ROOT_NAME
    render_queue_root = out_dir / "render_queues"
    sheet_root.mkdir(exist_ok=True)
    provenance_root.mkdir(exist_ok=True)
    render_queue_root.mkdir(exist_ok=True)
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
    render_queues: dict[str, dict[str, str]] = {}
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
        render_queue_path = (
            render_queue_root / f"review_queue_s{int(sector):04d}.csv"
        )
        public.loc[
            pd.to_numeric(public["sector"], errors="raise").eq(int(sector))
        ].to_csv(render_queue_path, index=False, float_format="%.15g")
        render_queues[str(int(sector))] = {
            "path": str(render_queue_path),
            "sha256": _sha256(render_queue_path),
        }
    sheet_lists["all"] = {
        "path": str(all_names_path),
        "sha256": _sha256(all_names_path),
    }
    if not labels_path.exists():
        pd.DataFrame(columns=LABEL_HEADER).to_csv(labels_path, index=False)

    app_path = ROOT / "scripts/stage5_validation/run_lightcurve_vetting_app.py"
    python_path = ROOT / ".venv/bin/python"
    command = (
        f"#!/usr/bin/env bash\nset -euo pipefail\n"
        f"cd {shlex.quote(str(ROOT))}\n"
        f"{shlex.quote(str(python_path))} {shlex.quote(str(app_path))} \\\n"
        f"  --candidates {shlex.quote(str(queue_path))} \\\n"
        f"  --labels-out {shlex.quote(str(labels_path))} \\\n"
        f"  --hlsp-root {shlex.quote(str(out_dir / '_sheet_only_no_hlsp'))} \\\n"
        f"  --twirl-vet-root {shlex.quote(str(sheet_root))} \\\n"
        f"  --exact-twirl-vet-sheets \\\n"
        f"  --labeler tehan --host 127.0.0.1 --port {int(port)}\n"
    )
    launcher = out_dir / "run_local_review.sh"
    launcher.write_text(command)
    launcher.chmod(0o755)
    _write_readme(out_dir, n_rows=len(public), port=int(port))
    sheet_set_path = out_dir / "sheet_set.json"
    sheet_set = {
        "sheet_set_version": SHEET_SET_VERSION,
        "sheet_root": SHEET_ROOT_NAME,
        "render_provenance_root": SHEET_PROVENANCE_ROOT_NAME,
        "branch_name": SHEET_BRANCH_NAME,
        "renderer_version": SHEET_RENDERER_VERSION,
        "apertures": list(SHEET_APERTURES),
        "harmonic_factors": list(SHEET_HARMONIC_FACTORS),
        "n_periods": SHEET_N_PERIODS,
        "n_peaks": SHEET_N_PEAKS,
        "lc_export_h5_by_sector": {
            str(sector): path
            for sector, path in SHEET_LC_EXPORT_H5_BY_SECTOR.items()
        },
        "use_row_ephemeris": True,
        "write_pdf": False,
        "require_render_provenance": True,
        "app": "twirl.vetting.lightcurve_label_app.LightCurveVettingApp",
        "exact_sheet_names_only": True,
    }
    sheet_set_path.write_text(
        json.dumps(sheet_set, indent=2, sort_keys=True) + "\n"
    )

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
                "render_queues": render_queues,
                "sheet_set": {
                    "path": str(sheet_set_path),
                    "sha256": _sha256(sheet_set_path),
                },
                "human_labels_final": {"path": str(labels_path)},
            },
            "review_port": int(port),
            "sheet_set": sheet_set,
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
    accepted = args.accepted_franklin or list(DEFAULT_FRANKLIN)
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
