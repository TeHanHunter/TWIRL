#!/usr/bin/env python3
"""Freeze and normalize one completed Franklin multisector label return."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil

import pandas as pd

from twirl.vetting.franklin_multisector import normalize_franklin_label_return


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _parse_mapping(values: list[str], *, value_kind: str) -> dict[int, str]:
    parsed: dict[int, str] = {}
    for value in values:
        try:
            sector_text, mapped = value.split("=", 1)
            sector = int(sector_text)
        except ValueError as exc:
            raise ValueError(
                f"invalid {value_kind} mapping {value!r}; expected SECTOR=VALUE"
            ) from exc
        if sector in parsed or not mapped:
            raise ValueError(f"invalid or duplicate {value_kind} mapping: {value!r}")
        parsed[sector] = mapped
    return parsed


def _snapshot(source: Path, destination: Path) -> Path:
    source = Path(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if source.resolve() == destination.resolve():
        return destination
    if destination.exists():
        if _sha256(source) == _sha256(destination):
            return destination
        raise FileExistsError(
            f"refusing to overwrite frozen input with different bytes: {destination}"
        )
    temporary = destination.with_name(destination.name + ".tmp")
    shutil.copy2(source, temporary)
    if _sha256(source) != _sha256(temporary):
        temporary.unlink(missing_ok=True)
        raise OSError(f"snapshot hash mismatch while copying {source}")
    temporary.replace(destination)
    return destination


def ingest(
    *,
    queue_path: Path,
    labels_path: Path,
    out_dir: Path,
    source_batch_id: str,
    morphology_adjudicator: str,
    morphology_accepted_utc: str,
    expected_sector_counts: dict[int, int],
    native_h5_by_sector: dict[int, str],
    expected_labeler: str = "franklin",
    expected_label_source: str = "human",
) -> dict[str, object]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    # Validate original sources before writing or replacing any frozen input.
    queue = pd.read_csv(queue_path, dtype=str, keep_default_na=False)
    labels = pd.read_csv(labels_path, dtype=str, keep_default_na=False)
    normalized = normalize_franklin_label_return(
        queue,
        labels,
        source_batch_id=source_batch_id,
        morphology_adjudicator=morphology_adjudicator,
        morphology_accepted_utc=morphology_accepted_utc,
        expected_sector_counts=expected_sector_counts,
        native_h5_by_sector=native_h5_by_sector,
        expected_labeler=expected_labeler,
        expected_label_source=expected_label_source,
    )
    queue_snapshot = _snapshot(
        Path(queue_path), out_dir / "frozen_review_queue.csv"
    )
    labels_snapshot = _snapshot(
        Path(labels_path), out_dir / "franklin_labels_returned.csv"
    )
    normalized_path = out_dir / "accepted_morphology_labels.csv"
    normalized.to_csv(normalized_path, index=False, float_format="%.15g")
    summary: dict[str, object] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_batch_id": source_batch_id,
        "label_unit": "sector_observation",
        "morphology_review_status": "accepted_batch_level",
        "morphology_adjudicator": morphology_adjudicator,
        "morphology_accepted_utc": morphology_accepted_utc,
        "expected_original_labeler": expected_labeler,
        "expected_label_source": expected_label_source,
        "period_supervision_policy": (
            "reported factor/status retained for audit; all real-row harmonic "
            "targets masked until an explicit factor-only review"
        ),
        "n_rows": int(len(normalized)),
        "n_unique_tics": int(pd.to_numeric(normalized["tic"], errors="raise").nunique()),
        "sector_counts": {
            str(key): int(value)
            for key, value in pd.to_numeric(
                normalized["sector"], errors="raise"
            ).value_counts().sort_index().items()
        },
        "human_label_counts": {
            str(key): int(value)
            for key, value in normalized["human_label"].value_counts().sort_index().items()
        },
        "morphology_target_counts": {
            str(key): int(value)
            for key, value in normalized.loc[
                normalized["morphology_include_v1"].astype(bool),
                "morphology_target_v1",
            ].value_counts().sort_index().items()
        },
        "n_harmonic_targets": int(normalized["harmonic_include_v1"].astype(bool).sum()),
        "n_native_h5_paths_resolved": int(
            normalized["native_h5_path"].fillna("").astype(str).ne("").sum()
        ),
        "inputs": {
            "frozen_review_queue": {
                "path": str(queue_snapshot),
                "sha256": _sha256(queue_snapshot),
            },
            "franklin_labels_returned": {
                "path": str(labels_snapshot),
                "sha256": _sha256(labels_snapshot),
            },
        },
        "outputs": {
            "accepted_morphology_labels": {
                "path": str(normalized_path),
                "sha256": _sha256(normalized_path),
            }
        },
    }
    summary_path = out_dir / "summary.json"
    summary["outputs"]["summary"] = {"path": str(summary_path)}  # type: ignore[index]
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queue", type=Path, required=True)
    parser.add_argument("--labels", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--source-batch-id", required=True)
    parser.add_argument("--morphology-adjudicator", required=True)
    parser.add_argument("--expected-labeler", default="franklin")
    parser.add_argument("--expected-label-source", default="human")
    parser.add_argument(
        "--morphology-accepted-utc",
        required=True,
        help="Frozen batch-level acceptance timestamp.",
    )
    parser.add_argument(
        "--expected-sector-count",
        action="append",
        default=[],
        metavar="SECTOR=N",
    )
    parser.add_argument(
        "--native-h5",
        action="append",
        default=[],
        metavar="SECTOR=PATH",
    )
    args = parser.parse_args()
    expected_raw = _parse_mapping(
        args.expected_sector_count, value_kind="expected-sector-count"
    )
    expected = {sector: int(value) for sector, value in expected_raw.items()}
    native = _parse_mapping(args.native_h5, value_kind="native-h5")
    summary = ingest(
        queue_path=args.queue,
        labels_path=args.labels,
        out_dir=args.out_dir,
        source_batch_id=args.source_batch_id,
        morphology_adjudicator=args.morphology_adjudicator,
        morphology_accepted_utc=args.morphology_accepted_utc,
        expected_sector_counts=expected,
        native_h5_by_sector=native,
        expected_labeler=args.expected_labeler,
        expected_label_source=args.expected_label_source,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
