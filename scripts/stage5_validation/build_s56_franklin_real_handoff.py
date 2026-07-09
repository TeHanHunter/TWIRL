#!/usr/bin/env python3
"""Build the Franklin real-only S56 vetting handoff scaffold."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sys
from typing import Any

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIR = Path(__file__).resolve().parent
SRC_ROOT = REPO_ROOT / "src"
for path in (str(SCRIPT_DIR), str(SRC_ROOT)):
    if path not in sys.path:
        sys.path.insert(0, path)

from build_s56_mixed_teacher_queue import (  # noqa: E402
    DEFAULT_REAL_CANDIDATES,
    LABEL_HEADER,
    _blind_for_browser,
    select_real_pool,
)


DEFAULT_OUT_DIR = REPO_ROOT / "reports/stage5_validation/s56_franklin_real5k_handoff"
DEFAULT_LABELED_JOINED = (
    REPO_ROOT
    / "reports/stage5_validation/s56_recovery50_teacher_queue_2k/human_label_audit/joined_human_labels.csv"
)
DEFAULT_REVISIT_QUEUE = (
    REPO_ROOT
    / "reports/stage5_validation/s56_recovery50_teacher_queue_2k/real_planet_revisit_wide/"
    / "review_queue_real_planet_revisit.csv"
)
DEFAULT_REVISIT_LABELS = (
    REPO_ROOT
    / "reports/stage5_validation/s56_recovery50_teacher_queue_2k/real_planet_revisit_wide/"
    / "human_labels_revisit.csv"
)
DEFAULT_EXCLUDE_QUEUE_CSVS = (
    REPO_ROOT / "reports/stage5_validation/s56_recovery50_teacher_queue/review_queue_1k.csv",
    REPO_ROOT / "reports/stage5_validation/s56_recovery50_teacher_queue_next1k/review_queue_1k.csv",
    REPO_ROOT / "reports/stage5_validation/s56_recovery50_teacher_queue_next4k/review_queue_4k.csv",
    REPO_ROOT / "reports/stage5_validation/s56_mixed_teacher_queue_pdo/review_queue_1k.csv",
)

EXAMPLE_SHEET_ROOTS = (
    REPO_ROOT / "reports/stage5_validation/s56_recovery50_teacher_queue/twirl_vet_sheets_fullphase_binmatch",
    REPO_ROOT / "reports/stage5_validation/s56_recovery50_teacher_queue_next4k/twirl_vet_sheets",
    REPO_ROOT
    / "reports/stage5_validation/s56_recovery50_teacher_queue_2k/real_planet_revisit_wide/"
    / "twirl_vet_sheets_secondary_phase",
)

REFERENCE_LABELS = (
    "planet_like",
    "wide_transit_like",
    "eclipsing_binary_or_pceb",
    "stellar_variability",
    "instrumental_or_systematic",
    "uncertain",
    "skip",
)


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=_json_default) + "\n")


def _safe_sheet_name(review_id: Any, branch_name: str = "current_adp") -> str:
    safe = str(review_id).replace("/", "_").replace(":", "_")
    return f"{safe}_twirl_twoap_{branch_name}.png"


def _select_stratified_remaining(
    frame: pd.DataFrame,
    *,
    n: int,
    random_state: int,
) -> pd.DataFrame:
    if len(frame) < n:
        raise ValueError(f"only {len(frame):,} available real rows after exclusions; need {n:,}")
    rng = np.random.default_rng(random_state)
    counts = frame["selection_bucket"].fillna("missing").astype(str).value_counts()
    quota_float = counts / counts.sum() * n
    quotas = np.floor(quota_float).astype(int).to_dict()
    remaining = n - sum(quotas.values())
    if remaining:
        fractions = (quota_float - np.floor(quota_float)).sort_values(ascending=False)
        for bucket in fractions.index[:remaining]:
            quotas[str(bucket)] = quotas.get(str(bucket), 0) + 1

    pieces: list[pd.DataFrame] = []
    used: set[int] = set()
    for offset, (bucket, quota) in enumerate(quotas.items()):
        if quota <= 0:
            continue
        available = frame.loc[frame["selection_bucket"].fillna("missing").astype(str).eq(bucket)]
        take = min(int(quota), len(available))
        if take:
            piece = available.sample(n=take, random_state=random_state + offset)
            pieces.append(piece)
            used.update(piece.index.tolist())
    selected = pd.concat(pieces, ignore_index=False) if pieces else frame.head(0)
    if len(selected) < n:
        fill = frame.loc[[idx for idx in frame.index if idx not in used]].sample(
            n=n - len(selected),
            random_state=random_state + 10_000,
        )
        selected = pd.concat([selected, fill], ignore_index=False)
    elif len(selected) > n:
        selected = selected.sample(n=n, random_state=random_state + 20_000)
    return selected.sample(frac=1.0, random_state=random_state + 30_000).reset_index(drop=True)


def _find_sheet(name: str, review_id: str, tic: Any, roots: tuple[Path, ...]) -> Path | None:
    names = [name] if name else []
    if review_id:
        names.append(_safe_sheet_name(review_id))
    if tic is not None and str(tic):
        names.append(f"*{int(float(tic))}*_twirl_twoap_*.png")
    for root in roots:
        if not root.exists():
            continue
        for candidate in names:
            if "*" in candidate:
                matches = sorted(root.glob(candidate))
                if matches:
                    return matches[0]
            else:
                path = root / candidate
                if path.exists():
                    return path
    return None


def _load_reference_rows(
    *,
    labeled_joined: Path,
    revisit_queue: Path,
    revisit_labels: Path,
) -> pd.DataFrame:
    pieces: list[pd.DataFrame] = []
    if labeled_joined.exists():
        labeled = pd.read_csv(labeled_joined)
        labeled = labeled.loc[labeled.get("source_kind", "").eq("real_candidate")].copy()
        labeled["reference_label"] = labeled["human_label"].fillna("").astype(str)
        labeled["reference_source"] = "main_2k_labels"
        pieces.append(labeled)
    if revisit_queue.exists() and revisit_labels.exists():
        queue = pd.read_csv(revisit_queue)
        labels = pd.read_csv(revisit_labels)
        joined = queue.merge(
            labels[["row_id", "label", "notes"]],
            left_index=True,
            right_on="row_id",
            how="inner",
            suffixes=("", "_revisit"),
        )
        joined["reference_label"] = joined["label_revisit"].fillna("").astype(str)
        joined["reference_source"] = "real_planet_revisit_wide"
        joined["human_label"] = joined["reference_label"]
        joined["notes"] = joined["notes_revisit"].fillna(joined.get("notes", ""))
        pieces.append(joined)
    if not pieces:
        return pd.DataFrame()
    out = pd.concat(pieces, ignore_index=True, sort=False)
    out = out.loc[out["reference_label"].isin(REFERENCE_LABELS)].copy()
    return out


def _copy_reference_examples(
    *,
    out_dir: Path,
    labeled_joined: Path,
    revisit_queue: Path,
    revisit_labels: Path,
    examples_per_label: int,
    random_state: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    ref = _load_reference_rows(
        labeled_joined=labeled_joined,
        revisit_queue=revisit_queue,
        revisit_labels=revisit_labels,
    )
    examples_dir = out_dir / "reference_examples"
    examples_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    missing: list[dict[str, Any]] = []
    for offset, label in enumerate(REFERENCE_LABELS):
        subset = ref.loc[ref["reference_label"].eq(label)].copy()
        if subset.empty:
            continue
        take = min(examples_per_label, len(subset))
        subset = subset.sample(n=take, random_state=random_state + offset)
        label_dir = examples_dir / label
        label_dir.mkdir(parents=True, exist_ok=True)
        for _, row in subset.iterrows():
            source_sheet = _find_sheet(
                str(row.get("twirl_vet_sheet_name", "") or ""),
                str(row.get("review_id", "") or ""),
                row.get("tic", ""),
                EXAMPLE_SHEET_ROOTS,
            )
            rec = {
                "label": label,
                "tic": row.get("tic", ""),
                "review_id": row.get("review_id", ""),
                "source": row.get("reference_source", ""),
                "notes": row.get("notes", ""),
                "source_sheet": str(source_sheet) if source_sheet else "",
                "example_sheet": "",
            }
            if source_sheet is None:
                missing.append(rec)
                rows.append(rec)
                continue
            dest = label_dir / f"{label}_{int(float(row.get('tic')))}.png"
            shutil.copy2(source_sheet, dest)
            rec["example_sheet"] = str(dest.relative_to(out_dir))
            rows.append(rec)
    examples = pd.DataFrame(rows)
    examples.to_csv(out_dir / "reference_examples.csv", index=False)
    summary = {
        "examples_per_label_requested": int(examples_per_label),
        "example_counts": examples["label"].value_counts().sort_index().to_dict() if len(examples) else {},
        "missing_example_sheets": missing,
    }
    return examples, summary


def _write_readme(out_dir: Path, *, n_rows: int, summary: dict[str, Any]) -> None:
    text = f"""# TWIRL S56 Real-Candidate Vetting Handoff

This package is a real-only S56 review queue for Franklin.

## What Is Included

- `franklin_review_queue_5k_real.csv`: {n_rows:,} real S56 candidate rows.
- `vet_sheets/`: pre-rendered TWIRL two-aperture PNG sheets for the browser app.
- `franklin_vetting_app.py`: standalone local browser app.
- `run_franklin_vetting.sh`: one-command launcher.
- `reference_examples/`: examples from TeHan's real-data labels.
- `reference_examples.csv`: table linking examples to labels and TICs.
- `franklin_labels_vetted.csv`: created by the app as labels are saved.

The main queue intentionally contains no injected rows. Injection truth is not
part of Franklin's label task.

## Start Labeling

From the unzipped package directory:

```bash
python3 franklin_vetting_app.py --check-only
./run_franklin_vetting.sh
```

Then open:

```text
http://127.0.0.1:5003/
```

Labels are saved immediately to `franklin_labels_vetted.csv`.

## Label Meanings

- `planet_like`: preserve a compact transit-like signal.
- `wide_transit_like`: preserve a broad/long-duration transit-like signal, but
  do not merge it blindly with compact planet morphology.
- `eclipsing_binary_or_pceb`: EB/PCEB-like shape, odd/even mismatch, secondary
  eclipse, or stellar-companion morphology.
- `stellar_variability`: astrophysical variability or repeating non-transit
  structure.
- `instrumental_or_systematic`: BLS is catching window-edge leakage, cadence
  artifacts, extraction artifacts, or other non-astrophysical structure.
- `uncertain`: flat/no obvious useful event. TeHan used this heavily for flats.
- `skip`: unusable row or app/sheet problem.

If a row looks period-folded at the wrong harmonic, keep the preserve/reject
label you think is scientifically correct and write a short note such as
`half period`, `refold at P/2`, or `possible harmonic`.

## Keyboard Shortcuts

- `1`: planet_like
- `6`: wide_transit_like
- `2`: eclipsing_binary_or_pceb
- `3`: stellar_variability
- `4`: instrumental_or_systematic
- `5`: uncertain
- `0`: skip
- left/right arrows: previous/next row

## Provenance

Created UTC: `{datetime.now(timezone.utc).isoformat()}`

Summary:

```json
{json.dumps(summary, indent=2, sort_keys=True, default=_json_default)}
```
"""
    (out_dir / "README_Franklin_vetting.md").write_text(text)


def _write_launcher(out_dir: Path) -> None:
    path = out_dir / "run_franklin_vetting.sh"
    path.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "cd \"$(dirname \"$0\")\"\n"
        "exec python3 franklin_vetting_app.py \\\n"
        "  --queue franklin_review_queue_5k_real.csv \\\n"
        "  --labels-out franklin_labels_vetted.csv \\\n"
        "  --sheet-root vet_sheets \\\n"
        "  --labeler franklin \\\n"
        "  --host 127.0.0.1 \\\n"
        "  --port 5003\n"
    )
    path.chmod(0o755)


def _review_ids_from_queue(path: Path) -> set[str]:
    if not path.exists():
        return set()
    table = pd.read_csv(path, usecols=lambda col: col in {"review_id", "source_kind"})
    if "review_id" not in table:
        return set()
    if "source_kind" in table:
        table = table.loc[table["source_kind"].fillna("").astype(str).eq("real_candidate")]
    return set(table["review_id"].fillna("").astype(str).replace("", np.nan).dropna().tolist())


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--real-candidates", type=Path, default=DEFAULT_REAL_CANDIDATES)
    parser.add_argument("--labeled-joined", type=Path, default=DEFAULT_LABELED_JOINED)
    parser.add_argument("--revisit-queue", type=Path, default=DEFAULT_REVISIT_QUEUE)
    parser.add_argument("--revisit-labels", type=Path, default=DEFAULT_REVISIT_LABELS)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--n-review", type=int, default=5000)
    parser.add_argument("--n-source-pool", type=int, default=9000)
    parser.add_argument("--random-state", type=int, default=561208)
    parser.add_argument("--cadence-alias-tolerance", type=float, default=0.02)
    parser.add_argument("--examples-per-label", type=int, default=2)
    parser.add_argument(
        "--exclude-queue-csv",
        type=Path,
        action="append",
        default=list(DEFAULT_EXCLUDE_QUEUE_CSVS),
        help=(
            "Prior queue CSV whose real review_id rows should be excluded from "
            "Franklin's fresh target set. May be repeated."
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    real_pool, real_summary = select_real_pool(
        args.real_candidates,
        n_real=args.n_source_pool,
        random_state=args.random_state,
        alias_tolerance=args.cadence_alias_tolerance,
    )
    labeled = pd.read_csv(args.labeled_joined) if args.labeled_joined.exists() else pd.DataFrame()
    exclude_review_ids = set(
        labeled.loc[
            labeled.get("source_kind", pd.Series("", index=labeled.index)).fillna("").astype(str).eq("real_candidate"),
            "review_id",
        ]
        .fillna("")
        .astype(str)
        .tolist()
    ) if len(labeled) and "review_id" in labeled else set()
    exclude_queue_counts: dict[str, int] = {}
    for path in args.exclude_queue_csv or []:
        ids = _review_ids_from_queue(path)
        exclude_queue_counts[str(path)] = len(ids)
        exclude_review_ids.update(ids)
    available = real_pool.loc[~real_pool["review_id"].fillna("").astype(str).isin(exclude_review_ids)].copy()
    review_unblinded = _select_stratified_remaining(
        available,
        n=args.n_review,
        random_state=args.random_state + 1000,
    )
    review = _blind_for_browser(review_unblinded)
    review["source_kind"] = "real_candidate"
    review["truth_source_kind"] = "real_candidate"
    review["twirl_vet_sheet_name"] = review["review_id"].map(_safe_sheet_name)
    review["twirl_vet_sheet_pdf_name"] = review["twirl_vet_sheet_name"].str.replace(".png", ".pdf", regex=False)
    for column in LABEL_HEADER:
        if column not in review:
            review[column] = ""
    for column in ("label", "label_source", "labeler", "notes", "updated_utc"):
        review[column] = ""
    queue_csv = args.out_dir / "franklin_review_queue_5k_real.csv"
    review.to_csv(queue_csv, index=False)

    labels_path = args.out_dir / "franklin_labels_vetted.csv"
    if not labels_path.exists():
        pd.DataFrame(columns=LABEL_HEADER).to_csv(labels_path, index=False)

    app_src = SCRIPT_DIR / "franklin_vetting_app.py"
    shutil.copy2(app_src, args.out_dir / "franklin_vetting_app.py")
    _write_launcher(args.out_dir)

    examples, example_summary = _copy_reference_examples(
        out_dir=args.out_dir,
        labeled_joined=args.labeled_joined,
        revisit_queue=args.revisit_queue,
        revisit_labels=args.revisit_labels,
        examples_per_label=args.examples_per_label,
        random_state=args.random_state + 2000,
    )

    verification_failures: list[str] = []
    if len(review) != args.n_review:
        verification_failures.append(f"queue rows {len(review)} != {args.n_review}")
    if review["source_kind"].fillna("").astype(str).ne("real_candidate").any():
        verification_failures.append("queue contains non-real rows")
    dup = int(review["review_id"].fillna("").astype(str).duplicated().sum())
    if dup:
        verification_failures.append(f"queue has {dup} duplicate review_id rows")
    overlap = int(review["review_id"].fillna("").astype(str).isin(exclude_review_ids).sum())
    if overlap:
        verification_failures.append(f"queue overlaps {overlap} already-labeled real rows")

    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "queue_csv": str(queue_csv),
        "labels_csv": str(labels_path),
        "n_review": int(len(review)),
        "n_source_pool": int(len(real_pool)),
        "n_available_after_excluding_prior_real_rows": int(len(available)),
        "n_excluded_prior_real_review_ids": int(len(exclude_review_ids)),
        "exclude_queue_counts": exclude_queue_counts,
        "real_summary": real_summary,
        "selection_bucket_counts": review["selection_bucket"].fillna("").astype(str).value_counts().sort_index().to_dict(),
        "source_kind_counts": review["source_kind"].fillna("").astype(str).value_counts().sort_index().to_dict(),
        "example_summary": example_summary,
        "n_reference_examples": int(len(examples)),
        "verification_passed": not verification_failures,
        "verification_failures": verification_failures,
    }
    _write_json(args.out_dir / "summary.json", summary)
    _write_readme(args.out_dir, n_rows=len(review), summary=summary)
    if verification_failures:
        for failure in verification_failures:
            print(f"[franklin-handoff] failure: {failure}", file=sys.stderr)
        return 1
    print(json.dumps(summary, indent=2, sort_keys=True, default=_json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
