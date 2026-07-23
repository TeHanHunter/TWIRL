#!/usr/bin/env python3
"""Score A2v1 recovery candidates with one manifest-pinned native-v2 ensemble."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import pandas as pd

from twirl.vetting.harmonic_inference import (
    RECOVERY_SCORE_ARTIFACT_CONTRACT,
    score_harmonic_teacher_ensemble,
)
from twirl.vetting.harmonic_training import (
    verify_deployed_teacher_checkpoint_manifest,
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _input_hashes(paths: tuple[Path, ...]) -> dict[str, str]:
    return {str(path): _sha256(path) for path in paths}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--native-h5", type=Path, required=True)
    parser.add_argument("--checkpoints", type=Path, nargs=5, required=True)
    parser.add_argument("--checkpoint-manifest", type=Path, required=True)
    parser.add_argument("--out-scores", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--filter-manifest", type=Path)
    parser.add_argument("--filter-key", default="injection_id")
    parser.add_argument("--allow-cpu", action="store_true")
    args = parser.parse_args()

    candidates_path = args.candidates.resolve()
    native_h5 = args.native_h5.resolve()
    checkpoint_paths = tuple(path.resolve() for path in args.checkpoints)
    checkpoint_manifest = args.checkpoint_manifest.resolve()
    filter_manifest = args.filter_manifest.resolve() if args.filter_manifest else None
    checkpoint_root = checkpoint_paths[0].parent.parent
    expected_checkpoints = tuple(
        checkpoint_root / f"fold_{fold}" / "teacher.pt" for fold in range(5)
    )
    if checkpoint_paths != expected_checkpoints:
        raise ValueError(
            "checkpoint arguments must be the ordered fold_0..fold_4 ensemble "
            f"under {checkpoint_root}"
        )
    input_paths = (
        candidates_path,
        native_h5,
        checkpoint_manifest,
        *checkpoint_paths,
        *((filter_manifest,) if filter_manifest is not None else ()),
    )
    before = _input_hashes(input_paths)
    checkpoint_manifest_payload = verify_deployed_teacher_checkpoint_manifest(
        manifest_path=checkpoint_manifest,
        checkpoint_root=checkpoint_root,
    )
    if _input_hashes(input_paths) != before:
        raise RuntimeError(
            "teacher recovery inputs changed during checkpoint-manifest verification"
        )
    with h5py.File(native_h5, "r") as h5:
        declared_training_table = Path(
            str(h5.attrs.get("training_table", ""))
        ).resolve()
        declared_training_sha256 = str(
            h5.attrs.get("training_table_sha256", "")
        )
    if declared_training_table != candidates_path:
        raise ValueError(
            "native recovery input is not bound to the exact candidate-table path"
        )
    if declared_training_sha256 != before[str(candidates_path)]:
        raise ValueError(
            "native recovery input training_table_sha256 does not match candidates"
        )

    candidates = (
        pd.read_parquet(candidates_path)
        if candidates_path.suffix.lower() == ".parquet"
        else pd.read_csv(candidates_path, low_memory=False)
    )
    n_source_rows = len(candidates)
    if filter_manifest is not None:
        filter_frame = (
            pd.read_parquet(filter_manifest)
            if filter_manifest.suffix.lower() == ".parquet"
            else pd.read_csv(filter_manifest, low_memory=False)
        )
        if args.filter_key not in candidates or args.filter_key not in filter_frame:
            raise KeyError(
                f"filter key {args.filter_key!r} is absent from candidates or manifest"
            )
        allowed = frozenset(filter_frame[args.filter_key].fillna("").astype(str))
        candidates = candidates.loc[
            candidates[args.filter_key].fillna("").astype(str).isin(allowed)
        ].copy()
    if args.limit is not None:
        candidates = candidates.head(int(args.limit)).copy()
    scores, summary = score_harmonic_teacher_ensemble(
        candidates=candidates,
        native_h5=native_h5,
        checkpoint_paths=checkpoint_paths,
        batch_size=args.batch_size,
        workers=args.workers,
        require_cuda=not args.allow_cpu,
        allow_injections=True,
    )
    out_scores = args.out_scores.resolve()
    if out_scores.suffix.lower() != ".parquet":
        raise ValueError("recovery teacher scores must use a .parquet output")
    out_scores.parent.mkdir(parents=True, exist_ok=True)
    score_temporary = out_scores.with_suffix(out_scores.suffix + ".tmp")
    summary_path = out_scores.with_suffix(".summary.json")
    summary_temporary = summary_path.with_suffix(summary_path.suffix + ".tmp")
    score_temporary.unlink(missing_ok=True)
    summary_temporary.unlink(missing_ok=True)
    publication_started = False
    try:
        scores.to_parquet(score_temporary, compression="zstd", index=False)
        after = _input_hashes(input_paths)
        if after != before:
            raise RuntimeError(
                "teacher recovery inputs changed before outputs were published"
            )
        expected_checkpoint_hashes = {
            str(path): before[str(path)] for path in checkpoint_paths
        }
        if summary.get("checkpoint_sha256") != expected_checkpoint_hashes:
            raise RuntimeError(
                "ensemble checkpoint hashes disagree with pre-scoring hashes"
            )
        summary.update(
            {
                "artifact_contract_version": RECOVERY_SCORE_ARTIFACT_CONTRACT,
                "strict_provenance_passed": True,
                "out_scores": str(out_scores),
                "out_scores_sha256": _sha256(score_temporary),
                "n_source_candidate_rows": int(n_source_rows),
                "filter_manifest": str(filter_manifest or ""),
                "filter_key": str(args.filter_key),
                "limit": int(args.limit) if args.limit is not None else None,
                "checkpoint_manifest": {
                    "path": str(checkpoint_manifest),
                    "sha256": before[str(checkpoint_manifest)],
                    "selected_profile": str(
                        checkpoint_manifest_payload["selected_profile"]
                    ),
                    "checkpoint_namespace": str(
                        checkpoint_manifest_payload["checkpoint_namespace"]
                    ),
                    "input_contract_version": str(
                        checkpoint_manifest_payload["input_contract_version"]
                    ),
                },
                "input_sha256_before": before,
                "input_sha256_after": after,
            }
        )
        summary_temporary.write_text(
            json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n"
        )
        if _input_hashes(input_paths) != before:
            raise RuntimeError(
                "teacher recovery inputs changed during output staging"
            )
        publication_started = True
        score_temporary.replace(out_scores)
        summary_temporary.replace(summary_path)
        if _input_hashes(input_paths) != before:
            out_scores.unlink(missing_ok=True)
            summary_path.unlink(missing_ok=True)
            raise RuntimeError(
                "teacher recovery inputs changed during output publication"
            )
    except Exception:
        score_temporary.unlink(missing_ok=True)
        summary_temporary.unlink(missing_ok=True)
        if publication_started:
            out_scores.unlink(missing_ok=True)
            summary_path.unlink(missing_ok=True)
        raise
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
