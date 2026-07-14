#!/usr/bin/env python3
"""Build immutable-derived S56 Teacher-v2 and S57 external-test role tables."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil

import pandas as pd

from twirl.vetting.teacher_v2 import (
    TEACHER_V2_ROLE_POLICY,
    build_global_tic_split_registry,
    build_franklin_a2v1_rereview_queue,
    build_franklin_current_a2v1_audit_candidates,
    build_real_human_training_rows,
    build_s56_injection_training_rows,
    build_s56_real_workload_candidates,
    build_s57_external_role_manifest,
    join_franklin_labels_to_queue,
    normalize_franklin_labels,
    transfer_human_labels_to_a2v1_candidates,
)


def _read(path: Path) -> pd.DataFrame:
    return pd.read_parquet(path) if path.suffix.lower() == ".parquet" else pd.read_csv(path, low_memory=False)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _write(frame: pd.DataFrame, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(path, compression="zstd", index=False)
    return path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--s56-schedule", type=Path, required=True)
    parser.add_argument("--s56-candidates", type=Path, required=True)
    parser.add_argument("--s56-injection-native-h5", type=Path, required=True)
    parser.add_argument("--human-table", type=Path, required=True)
    parser.add_argument("--human-native-h5", type=Path, required=True)
    parser.add_argument("--human-prior-splits", type=Path)
    parser.add_argument("--human-a2v1-compatibility", type=Path)
    parser.add_argument("--s56-real-candidates", type=Path, required=True)
    parser.add_argument("--franklin-labels", type=Path, required=True)
    parser.add_argument("--franklin-queue", type=Path)
    parser.add_argument(
        "--prepare-franklin-transfer-audit",
        action="store_true",
        help=(
            "Opt in to the deferred Franklin A2v1 transfer/re-review audit. "
            "Initial Teacher-v2 training only snapshots his labels and excludes "
            "their TICs from evaluation partitions."
        ),
    )
    parser.add_argument("--s57-schedule", type=Path, required=True)
    parser.add_argument("--s57-candidates", type=Path, required=True)
    parser.add_argument("--s57-native-h5", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--allow-incomplete-s57", action="store_true")
    args = parser.parse_args()

    sources = {
        "s56_schedule": args.s56_schedule,
        "s56_candidates": args.s56_candidates,
        "s56_injection_native_h5": args.s56_injection_native_h5,
        "human_table": args.human_table,
        "human_native_h5": args.human_native_h5,
        "s56_real_candidates": args.s56_real_candidates,
        "franklin_labels": args.franklin_labels,
        "s57_schedule": args.s57_schedule,
        "s57_candidates": args.s57_candidates,
        "s57_native_h5": args.s57_native_h5,
    }
    if args.prepare_franklin_transfer_audit:
        if args.franklin_queue is None:
            parser.error(
                "--prepare-franklin-transfer-audit requires --franklin-queue"
            )
        sources["franklin_queue"] = args.franklin_queue
    if args.human_prior_splits:
        sources["human_prior_splits"] = args.human_prior_splits
    if args.human_a2v1_compatibility:
        sources["human_a2v1_compatibility"] = args.human_a2v1_compatibility
    missing = [f"{name}={path}" for name, path in sources.items() if not path.exists()]
    if missing:
        raise FileNotFoundError(f"missing Teacher-v2 inputs: {missing}")

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    snapshot_dir = out_dir / "source_snapshots"
    snapshot_dir.mkdir(parents=True, exist_ok=True)
    franklin_snapshot = snapshot_dir / "franklin_labels_vetted.csv"
    if args.franklin_labels.resolve() != franklin_snapshot.resolve():
        shutil.copy2(args.franklin_labels, franklin_snapshot)

    s56_schedule = _read(args.s56_schedule)
    s56_candidates = _read(args.s56_candidates)
    injection_rows, injection_roles = build_s56_injection_training_rows(
        s56_schedule,
        s56_candidates,
        native_h5=args.s56_injection_native_h5,
    )

    human = _read(args.human_table)
    if args.human_prior_splits:
        prior = _read(args.human_prior_splits)
        split_columns = prior[["review_id", "fixed_split", "cv_fold"]].drop_duplicates(
            "review_id", keep="last"
        )
        human = human.drop(columns=["fixed_split", "cv_fold"], errors="ignore")
        if "source_review_id" in human:
            split_columns = split_columns.rename(columns={"review_id": "source_review_id"})
            human = human.merge(
                split_columns,
                on="source_review_id",
                how="left",
                validate="one_to_one",
            )
        else:
            human = human.merge(
                split_columns, on="review_id", how="left", validate="one_to_one"
            )
        if human[["fixed_split", "cv_fold"]].isna().any().any():
            raise RuntimeError(
                "A2v1 human transfers are missing prior Teacher-v1 grouped splits"
            )
    compatibility = (
        _read(args.human_a2v1_compatibility)
        if args.human_a2v1_compatibility
        else None
    )
    human_rows = build_real_human_training_rows(
        human,
        native_h5=args.human_native_h5,
        compatibility=compatibility,
    )
    franklin = normalize_franklin_labels(_read(args.franklin_labels))
    s57_schedule = _read(args.s57_schedule)
    registry = build_global_tic_split_registry(
        s56_injection_roles=injection_roles,
        s56_human_rows=human_rows,
        franklin_rows=franklin,
        s57_schedule=s57_schedule,
    )
    s57_external, s57_summary = build_s57_external_role_manifest(
        s57_schedule, registry
    )
    real_candidates = _read(args.s56_real_candidates)
    real_workload, real_workload_summary = build_s56_real_workload_candidates(
        real_candidates, registry
    )
    franklin_outputs: dict[str, Path] = {}
    franklin_a2v1 = pd.DataFrame()
    franklin_compatibility = pd.DataFrame()
    franklin_rereview_summary: dict[str, object] = {
        "status": "deferred_not_on_teacher_v2_critical_path"
    }
    franklin_current_a2v1_summary: dict[str, object] = {
        "status": "deferred_not_on_teacher_v2_critical_path"
    }
    if args.prepare_franklin_transfer_audit:
        franklin = join_franklin_labels_to_queue(
            _read(args.franklin_labels),
            pd.read_csv(args.franklin_queue, dtype=str, keep_default_na=False),
        )
        franklin_for_transfer = franklin.copy()
        franklin_for_transfer["source_uid"] = franklin_for_transfer[
            "franklin_identity"
        ]
        franklin_for_transfer["review_id"] = franklin_for_transfer[
            "franklin_queue_review_id"
        ]
        franklin_for_transfer["human_label"] = franklin_for_transfer["label"]
        franklin_for_transfer["is_injected_row"] = False
        franklin_a2v1, franklin_compatibility = (
            transfer_human_labels_to_a2v1_candidates(
                franklin_for_transfer,
                real_candidates,
            )
        )
        franklin_a2v1["teacher_v2_role"] = "franklin_a2v1_audit_only"
        franklin_a2v1["teacher_v2_training_include"] = False
        franklin_public, franklin_private, franklin_rereview_summary = (
            build_franklin_a2v1_rereview_queue(franklin, real_candidates)
        )
        franklin_current_a2v1, franklin_current_a2v1_summary = (
            build_franklin_current_a2v1_audit_candidates(
                franklin,
                real_candidates,
                franklin_compatibility,
            )
        )
        franklin_rereview_csv = out_dir / "franklin_a2v1_rereview_queue.csv"
        franklin_public.to_csv(franklin_rereview_csv, index=False)
        franklin_outputs = {
            "franklin_a2v1_audit_candidates": _write(
                franklin_a2v1,
                out_dir / "franklin_a2v1_audit_candidates.parquet",
            ),
            "franklin_a2v1_compatibility": _write(
                franklin_compatibility,
                out_dir / "franklin_a2v1_compatibility.parquet",
            ),
            "franklin_current_a2v1_audit_candidates": _write(
                franklin_current_a2v1,
                out_dir / "franklin_current_a2v1_audit_candidates.parquet",
            ),
            "franklin_a2v1_rereview_queue": _write(
                franklin_public,
                out_dir / "franklin_a2v1_rereview_queue.parquet",
            ),
            "franklin_a2v1_rereview_queue_csv": franklin_rereview_csv,
            "franklin_a2v1_rereview_private_manifest": _write(
                franklin_private,
                out_dir / "franklin_a2v1_rereview_private_manifest.parquet",
            ),
        }
    if not s57_summary["passed"] and not args.allow_incomplete_s57:
        raise RuntimeError(
            "existing S57 schedule fails the external-test support gate; regenerate it"
        )

    training = pd.concat([human_rows, injection_rows], ignore_index=True, sort=False)
    if training["review_id"].fillna("").astype(str).duplicated().any():
        raise RuntimeError("combined Teacher-v2 training review IDs are not unique")
    training_tics = set(training.loc[training["fixed_split"].eq("development"), "tic"])
    holdout_tics = set(training.loc[training["fixed_split"].eq("test"), "tic"])
    if training_tics & holdout_tics:
        raise RuntimeError("Teacher-v2 development and holdout TICs overlap")
    if set(s57_external["tic"]) & (training_tics | holdout_tics):
        raise RuntimeError("Teacher-v2 S57 external TICs overlap S56 partitions")

    outputs = {
        "training_rows": _write(training, out_dir / "training_rows.parquet"),
        "s56_injection_role_manifest": _write(
            injection_roles, out_dir / "s56_injection_role_manifest.parquet"
        ),
        "s56_injection_development_manifest": _write(
            injection_roles.loc[
                injection_roles["teacher_v2_partition"].astype(str).eq("development")
            ].copy(),
            out_dir / "s56_injection_development_manifest.parquet",
        ),
        "s56_injection_locked_holdout_manifest": _write(
            injection_roles.loc[
                injection_roles["teacher_v2_partition"].astype(str).eq(
                    "locked_holdout"
                )
            ].copy(),
            out_dir / "s56_injection_locked_holdout_manifest.parquet",
        ),
        "s56_human_role_manifest": _write(
            human_rows, out_dir / "s56_human_role_manifest.parquet"
        ),
        "franklin_audit_rows": _write(
            franklin, out_dir / "franklin_audit_rows.parquet"
        ),
        "global_tic_split_registry": _write(
            registry, out_dir / "global_tic_split_registry.parquet"
        ),
        "s57_external_role_manifest": _write(
            s57_external, out_dir / "s57_external_role_manifest.parquet"
        ),
        "s56_real_workload_candidates": _write(
            real_workload, out_dir / "s56_real_workload_candidates.parquet"
        ),
        **franklin_outputs,
    }
    source_sha256 = {name: _sha256(path) for name, path in sources.items() if path.is_file()}
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "role_policy": TEACHER_V2_ROLE_POLICY,
        "n_training_rows": int(len(training)),
        "n_human_rows": int(len(human_rows)),
        "n_injection_rows": int(len(injection_rows)),
        "injection_role_counts": {
            str(key): int(value)
            for key, value in injection_rows["teacher_v2_role"].value_counts().items()
        },
        "fixed_split_counts": {
            str(key): int(value) for key, value in training["fixed_split"].value_counts().items()
        },
        "n_franklin_rows": int(len(franklin)),
        "franklin_transfer_audit_status": (
            "prepared"
            if args.prepare_franklin_transfer_audit
            else "deferred_not_on_teacher_v2_critical_path"
        ),
        "n_franklin_a2v1_compatible": int(len(franklin_a2v1)),
        "franklin_a2v1_transfer_fraction": (
            float(franklin_compatibility["a2v1_transfer_ok"].mean())
            if len(franklin_compatibility)
            else 0.0
        ),
        "franklin_rereview": franklin_rereview_summary,
        "franklin_current_a2v1_score_audit": franklin_current_a2v1_summary,
        "franklin_label_counts": {
            str(key): int(value) for key, value in franklin["label"].value_counts().items()
        },
        "s57_external": s57_summary,
        "s56_real_workload": real_workload_summary,
        "source_sha256": source_sha256,
        "franklin_snapshot_sha256": _sha256(franklin_snapshot),
        "outputs": {name: str(path) for name, path in outputs.items()},
    }
    (out_dir / "role_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
