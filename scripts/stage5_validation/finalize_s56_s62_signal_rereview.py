#!/usr/bin/env python3
"""Freeze a complete S56--S62 Planet-like/EB morphology re-review."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import pandas as pd

from twirl.vetting.multisector_signal_review import finalize_signal_rereview


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def finalize(
    *,
    queue_path: Path,
    labels_path: Path,
    provenance_path: Path,
    out_dir: Path,
    adjudicator: str,
    accepted_utc: str,
) -> dict[str, object]:
    queue = pd.read_csv(queue_path, dtype=str, keep_default_na=False)
    labels = pd.read_csv(labels_path, dtype=str, keep_default_na=False)
    provenance = pd.read_csv(
        provenance_path, dtype=str, keep_default_na=False
    )
    required_provenance = {
        "source_uid",
        "source_batch_id",
        "observation_candidate_key",
        "selected_queue_row_id",
        "selected_for_queue",
    }
    missing_provenance = sorted(required_provenance - set(provenance.columns))
    if missing_provenance:
        raise KeyError(
            f"source provenance is missing columns: {missing_provenance}"
        )
    required_queue = {
        "row_id",
        "observation_candidate_key",
        "selected_source_uid",
    }
    missing_queue = sorted(required_queue - set(queue.columns))
    if missing_queue:
        raise KeyError(
            f"review queue is missing provenance columns: {missing_queue}"
        )
    selected = provenance["selected_for_queue"].str.lower().isin(
        {"1", "1.0", "true", "t", "yes", "y"}
    )
    selected_provenance = provenance.loc[
        selected,
        [
            "selected_queue_row_id",
            "observation_candidate_key",
            "source_uid",
        ],
    ].copy()
    if (
        len(selected_provenance) != len(queue)
        or selected_provenance["selected_queue_row_id"].duplicated().any()
    ):
        raise ValueError(
            "source provenance does not bind exactly one selected source to "
            "every final queue row"
        )
    selected_provenance = selected_provenance.rename(
        columns={
            "selected_queue_row_id": "row_id",
            "observation_candidate_key": "provenance_observation_candidate_key",
            "source_uid": "provenance_source_uid",
        }
    )
    binding = queue.loc[
        :,
        ["row_id", "observation_candidate_key", "selected_source_uid"],
    ].merge(
        selected_provenance,
        on="row_id",
        how="left",
        validate="one_to_one",
    )
    binding_mismatch = (
        binding["provenance_observation_candidate_key"].isna()
        | binding["observation_candidate_key"].ne(
            binding["provenance_observation_candidate_key"]
        )
        | binding["selected_source_uid"].ne(binding["provenance_source_uid"])
    )
    if binding_mismatch.any():
        raise ValueError(
            "source provenance does not match the exact queue row, candidate, "
            "and selected source identity"
        )
    final = finalize_signal_rereview(
        queue,
        labels,
        adjudicator=adjudicator,
        accepted_utc=accepted_utc,
    )
    provenance_sha256 = _sha256(provenance_path)
    final["source_label_provenance_sha256"] = provenance_sha256
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    final_path = out_dir / "accepted_signal_rereview.csv"
    if final_path.exists():
        raise FileExistsError(f"refusing to overwrite frozen review: {final_path}")
    final.to_csv(final_path, index=False, float_format="%.15g")
    summary: dict[str, object] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "adjudicator": adjudicator,
        "accepted_utc": accepted_utc,
        "n_rows": int(len(final)),
        "n_unique_tics": int(
            pd.to_numeric(final["tic"], errors="raise").nunique()
        ),
        "sector_counts": {
            str(key): int(value)
            for key, value in pd.to_numeric(
                final["sector"], errors="raise"
            ).value_counts().sort_index().items()
        },
        "final_label_counts": {
            str(key): int(value)
            for key, value in final["final_human_label"]
            .value_counts()
            .sort_index()
            .items()
        },
        "n_preserved_verified_harmonics": int(
            final["harmonic_include_v1"].astype(bool).sum()
        ),
        "inputs": {
            "review_queue": {
                "path": str(queue_path),
                "sha256": _sha256(queue_path),
            },
            "human_labels_final": {
                "path": str(labels_path),
                "sha256": _sha256(labels_path),
            },
            "source_label_provenance": {
                "path": str(provenance_path),
                "sha256": provenance_sha256,
            },
        },
        "outputs": {
            "accepted_signal_rereview": {
                "path": str(final_path),
                "sha256": _sha256(final_path),
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
    parser.add_argument(
        "--source-provenance",
        type=Path,
        default=None,
        help="Defaults to source_label_provenance.csv next to the queue.",
    )
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--adjudicator", default="tehan")
    parser.add_argument(
        "--accepted-utc",
        default=None,
        help="Explicit UTC acceptance timestamp; defaults to the current time.",
    )
    args = parser.parse_args()
    accepted_utc = args.accepted_utc or datetime.now(timezone.utc).isoformat()
    provenance_path = (
        args.source_provenance
        or args.queue.parent / "source_label_provenance.csv"
    )
    summary = finalize(
        queue_path=args.queue,
        labels_path=args.labels,
        provenance_path=provenance_path,
        out_dir=args.out_dir,
        adjudicator=args.adjudicator,
        accepted_utc=accepted_utc,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
