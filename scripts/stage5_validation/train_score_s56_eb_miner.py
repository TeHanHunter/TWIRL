#!/usr/bin/env python3
"""Train and apply the S56 EB/PCEB active-learning CNN miner."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.vetting.eb_miner import (  # noqa: E402
    EB_TARGET_POSITIVE,
    aggregate_ensemble_scores,
)
from twirl.vetting.recovery50_cnn import (  # noqa: E402
    CnnTrainConfig,
    TensorConfig,
    build_recovery50_cnn_tensors,
    score_recovery50_cnn_teacher,
    train_recovery50_cnn_teacher,
)
from twirl.vetting.recovery50_teacher import DEFAULT_APERTURES, json_default, read_table, write_table  # noqa: E402


DEFAULT_ROOT = REPO_ROOT / "reports/stage5_validation/s56_eb_miner_adp_only"
DEFAULT_TRAINING_TABLE = DEFAULT_ROOT / "training/eb_miner_training_table.csv"
DEFAULT_CANDIDATE_POOL = DEFAULT_ROOT / "candidate_pool/eb_miner_candidate_pool.csv"
DEFAULT_COMPACT_LC = (
    REPO_ROOT / "data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_adp_lc_export_pdo.h5"
)
DEFAULT_HLSP_ROOT = REPO_ROOT / "data_local/stage1_lightcurves/hlsp_s0056_twirl_fs_v2_compare"


def _average_precision(y_true: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y_true, dtype=bool)
    s = np.asarray(score, dtype=float)
    good = np.isfinite(s)
    y = y[good]
    s = s[good]
    if y.sum() == 0:
        return float("nan")
    order = np.argsort(-s, kind="mergesort")
    y = y[order]
    cum_pos = np.cumsum(y)
    precision = cum_pos / np.arange(1, len(y) + 1)
    return float(np.sum(precision[y]) / max(int(y.sum()), 1))


def _binary_eval(predictions: pd.DataFrame, *, positive_label: str = EB_TARGET_POSITIVE) -> dict[str, Any]:
    out: dict[str, Any] = {}
    target = predictions["main_teacher_target"].fillna("").astype(str)
    pred = predictions["cnn_label"].fillna("").astype(str)
    score_col = f"cnn_p_{positive_label}"
    score = pd.to_numeric(predictions.get(score_col, np.nan), errors="coerce").to_numpy(dtype=float)
    y_true = target.eq(positive_label).to_numpy()
    y_pred = pred.eq(positive_label).to_numpy()
    for split_name, part in [("all", predictions)] + list(predictions.groupby("cnn_training_split", dropna=False)):
        idx = part.index.to_numpy()
        yt = y_true[idx]
        yp = y_pred[idx]
        sc = score[idx]
        tp = int(np.sum(yt & yp))
        fp = int(np.sum(~yt & yp))
        fn = int(np.sum(yt & ~yp))
        tn = int(np.sum(~yt & ~yp))
        precision = tp / (tp + fp) if tp + fp else float("nan")
        recall = tp / (tp + fn) if tp + fn else float("nan")
        specificity = tn / (tn + fp) if tn + fp else float("nan")
        topk: dict[str, Any] = {}
        order = np.argsort(-np.nan_to_num(sc, nan=-np.inf), kind="mergesort")
        for k in (5, 10, 20, 50):
            kk = min(k, len(order))
            if kk <= 0:
                continue
            top = yt[order[:kk]]
            topk[str(k)] = {"n": int(kk), "n_positive": int(top.sum()), "positive_fraction": float(top.mean())}
        out[str(split_name)] = {
            "n": int(len(part)),
            "n_positive": int(yt.sum()),
            "n_negative": int((~yt).sum()),
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "tn": tn,
            "precision": float(precision),
            "recall": float(recall),
            "specificity": float(specificity),
            "average_precision": _average_precision(yt, sc),
            "topk": topk,
        }
    return out


def _training_eval_for_profile(profile_dir: Path) -> dict[str, Any]:
    pred_path = profile_dir / "cnn_predictions.parquet"
    if not pred_path.exists():
        pred_path = profile_dir / "cnn_predictions.csv"
    pred = read_table(pred_path)
    return _binary_eval(pred)


def _chunks(frame: pd.DataFrame, chunk_size: int) -> list[pd.DataFrame]:
    return [frame.iloc[start:start + chunk_size].copy() for start in range(0, len(frame), chunk_size)]


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text()) if path.exists() else {}


def _resolve_score_profile(summary: dict[str, Any], requested: str) -> str:
    profile = (
        str(summary.get("best_profile_by_validation_balanced_accuracy", ""))
        if requested == "auto"
        else str(requested)
    )
    if profile not in {"cnn_shape_only", "cnn_shape_plus_bls"}:
        raise ValueError(f"invalid or unavailable EB score profile: {profile!r}")
    return profile


def _score_prediction_path(chunk_root: Path, member_idx: int) -> Path | None:
    score_dir = chunk_root / "scores" / f"member_{member_idx:02d}"
    for name in ("cnn_scored_candidates.parquet", "cnn_scored_candidates.csv"):
        path = score_dir / name
        if path.exists() and path.stat().st_size > 0:
            return path
    return None


def _chunk_score_paths(chunk_root: Path, n_members: int) -> list[Path]:
    paths: list[Path] = []
    for member_idx in range(n_members):
        path = _score_prediction_path(chunk_root, member_idx)
        if path is not None:
            paths.append(path)
    return paths


def _tensor_summary_if_complete(chunk_root: Path) -> dict[str, Any] | None:
    summary_path = chunk_root / "tensors" / "summary.json"
    if not summary_path.exists():
        return None
    summary = _read_json(summary_path)
    outputs = summary.get("outputs", {})
    if not outputs.get("npz") or not outputs.get("rows"):
        return None
    npz = Path(outputs["npz"])
    rows = Path(outputs["rows"])
    if npz.exists() and rows.exists() and npz.stat().st_size > 0 and rows.stat().st_size > 0:
        return summary
    return None


def _build_score_chunk_tensors_task(
    *,
    chunk_idx: int,
    chunk_csv: str,
    chunk_root: str,
    compact_lc_h5: str,
    hlsp_root: str,
    tensor_config: dict[str, Any],
) -> tuple[int, dict[str, Any]]:
    summary = build_recovery50_cnn_tensors(
        training_table=Path(chunk_csv),
        out_dir=Path(chunk_root) / "tensors",
        compact_lc_h5=Path(compact_lc_h5),
        hlsp_root=Path(hlsp_root),
        injection_h5_override=None,
        config=TensorConfig(**tensor_config),
        progress_every=100,
    )
    return int(chunk_idx), summary


def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--training-table", type=Path, default=DEFAULT_TRAINING_TABLE)
    ap.add_argument("--candidate-pool", type=Path, default=DEFAULT_CANDIDATE_POOL)
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_ROOT / "eb_cnn_miner")
    ap.add_argument("--compact-lc-h5", type=Path, default=DEFAULT_COMPACT_LC)
    ap.add_argument("--hlsp-root", type=Path, default=DEFAULT_HLSP_ROOT)
    ap.add_argument("--apertures", default=",".join(DEFAULT_APERTURES))
    ap.add_argument("--ensemble-size", type=int, default=5)
    ap.add_argument("--seed", type=int, default=56017)
    ap.add_argument("--epochs", type=int, default=80)
    ap.add_argument("--batch-size", type=int, default=32)
    ap.add_argument("--score-batch-size", type=int, default=256)
    ap.add_argument("--candidate-chunk-size", type=int, default=512)
    ap.add_argument("--score-chunk-workers", type=int, default=1)
    ap.add_argument("--learning-rate", type=float, default=5.0e-4)
    ap.add_argument("--weight-decay", type=float, default=1.0e-4)
    ap.add_argument("--dropout", type=float, default=0.25)
    ap.add_argument("--early-stop-patience", type=int, default=16)
    ap.add_argument("--validation-fraction", type=float, default=0.20)
    ap.add_argument("--test-fraction", type=float, default=0.20)
    ap.add_argument(
        "--score-profile",
        choices=("auto", "cnn_shape_only", "cnn_shape_plus_bls"),
        default="auto",
        help="Profile to score with; auto selects each ensemble member by validation balanced accuracy.",
    )
    ap.add_argument("--allow-cpu", action="store_true")
    ap.add_argument("--reuse-existing-models", action="store_true")
    ap.add_argument("--no-resume-scoring", dest="resume_scoring", action="store_false")
    ap.add_argument("--skip-scoring", action="store_true")
    ap.set_defaults(resume_scoring=True)
    return ap


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    apertures = tuple(part.strip() for part in args.apertures.split(",") if part.strip())
    tensor_cfg = TensorConfig(apertures=apertures)

    train_tensor_dir = args.out_dir / "training_tensors"

    member_summaries: list[dict[str, Any]] = []
    model_paths: list[Path] = []
    selected_profiles: list[str] = []
    if args.reuse_existing_models:
        tensor_summary = _read_json(train_tensor_dir / "summary.json") or {
            "status": "reused_existing_models",
            "out_dir": str(train_tensor_dir),
        }
        for member in range(args.ensemble_size):
            seed = int(args.seed + member)
            member_dir = args.out_dir / "ensemble" / f"seed_{seed}"
            profile_summary = _read_json(member_dir / "summary.json")
            selected_profile = _resolve_score_profile(profile_summary, args.score_profile)
            profile_dir = member_dir / selected_profile
            model_path = profile_dir / "cnn_teacher.pt"
            if not model_path.exists():
                raise FileNotFoundError(f"missing reused EB model: {model_path}")
            selected_profiles.append(selected_profile)
            model_paths.append(model_path)
            profile_dir = model_path.parent
            member_summaries.append(
                {
                    "status": "reused_existing_model",
                    "seed": seed,
                    "score_profile": selected_profile,
                    "model_path": str(model_path),
                    "profile_summary": _read_json(profile_dir / "summary.json"),
                    "eb_binary_eval": _training_eval_for_profile(profile_dir),
                }
            )
    else:
        tensor_summary = build_recovery50_cnn_tensors(
            training_table=args.training_table,
            out_dir=train_tensor_dir,
            compact_lc_h5=args.compact_lc_h5,
            hlsp_root=args.hlsp_root,
            injection_h5_override=None,
            config=tensor_cfg,
            progress_every=25,
        )
        tensor_npz = Path(tensor_summary["outputs"]["npz"])
        tensor_rows = Path(tensor_summary["outputs"]["rows"])

        for member in range(args.ensemble_size):
            seed = int(args.seed + member)
            member_dir = args.out_dir / "ensemble" / f"seed_{seed}"
            cfg = CnnTrainConfig(
                epochs=args.epochs,
                batch_size=args.batch_size,
                learning_rate=args.learning_rate,
                weight_decay=args.weight_decay,
                dropout=args.dropout,
                early_stop_patience=args.early_stop_patience,
                min_class_count=1,
                validation_fraction=args.validation_fraction,
                test_fraction=args.test_fraction,
                seed=seed,
                require_cuda=not args.allow_cpu,
            )
            summary = train_recovery50_cnn_teacher(
                tensor_npz=tensor_npz,
                tensor_rows=tensor_rows,
                training_table=args.training_table,
                metrics_tables=(),
                out_dir=member_dir,
                cfg=cfg,
            )
            selected_profile = _resolve_score_profile(summary, args.score_profile)
            selected_profiles.append(selected_profile)
            profile_dir = member_dir / selected_profile
            model_paths.append(profile_dir / "cnn_teacher.pt")
            summary["score_profile"] = selected_profile
            summary["eb_binary_eval"] = _training_eval_for_profile(profile_dir)
            member_summaries.append(summary)

    def score_chunk(
        *,
        chunk_idx: int,
        chunk_root: Path,
        chunk_csv: Path,
        tensor_summary_chunk: dict[str, Any],
    ) -> list[Path]:
        paths: list[Path] = []
        for member_idx, model_path in enumerate(model_paths):
            existing = _score_prediction_path(chunk_root, member_idx) if args.resume_scoring else None
            if existing is not None:
                print(
                    f"[eb-score] chunk={chunk_idx:04d} member={member_idx:02d} reused={existing}",
                    flush=True,
                )
                paths.append(existing)
                continue
            score_summary = score_recovery50_cnn_teacher(
                tensor_npz=Path(tensor_summary_chunk["outputs"]["npz"]),
                tensor_rows=Path(tensor_summary_chunk["outputs"]["rows"]),
                model_path=model_path,
                feature_table=chunk_csv,
                out_dir=chunk_root / "scores" / f"member_{member_idx:02d}",
                batch_size=args.score_batch_size,
                require_cuda=not args.allow_cpu,
            )
            path = Path(score_summary["outputs"]["predictions"])
            paths.append(path)
            print(
                f"[eb-score] chunk={chunk_idx:04d} member={member_idx:02d} scored={path}",
                flush=True,
            )
        return paths

    score_tables: list[Path] = []
    scoring_summary: dict[str, Any] = {"skipped": bool(args.skip_scoring)}
    if not args.skip_scoring:
        candidate_pool = read_table(args.candidate_pool)
        chunk_dir = args.out_dir / "score_chunks"
        chunk_infos: list[tuple[int, Path, Path]] = []
        for chunk_idx, chunk in enumerate(_chunks(candidate_pool, args.candidate_chunk_size)):
            one = chunk.copy().reset_index(drop=True)
            one["row_id"] = np.arange(len(one), dtype=int)
            chunk_root = chunk_dir / f"chunk_{chunk_idx:04d}"
            chunk_root.mkdir(parents=True, exist_ok=True)
            chunk_csv = chunk_root / "candidate_chunk.csv"
            one.to_csv(chunk_csv, index=False)
            chunk_infos.append((chunk_idx, chunk_root, chunk_csv))

        score_chunk_workers = max(1, int(args.score_chunk_workers))
        futures: dict[Any, tuple[int, Path, Path]] = {}
        ready: list[tuple[int, Path, Path, dict[str, Any]]] = []
        print(
            f"[eb-score] chunks={len(chunk_infos)} workers={score_chunk_workers} resume={args.resume_scoring}",
            flush=True,
        )
        if score_chunk_workers > 1:
            with ProcessPoolExecutor(max_workers=score_chunk_workers) as executor:
                for chunk_idx, chunk_root, chunk_csv in chunk_infos:
                    existing_paths = _chunk_score_paths(chunk_root, len(model_paths)) if args.resume_scoring else []
                    if len(existing_paths) == len(model_paths):
                        score_tables.extend(existing_paths)
                        print(f"[eb-score] chunk={chunk_idx:04d} reused_complete_scores", flush=True)
                        continue
                    summary = _tensor_summary_if_complete(chunk_root)
                    if summary is not None:
                        ready.append((chunk_idx, chunk_root, chunk_csv, summary))
                        continue
                    future = executor.submit(
                        _build_score_chunk_tensors_task,
                        chunk_idx=chunk_idx,
                        chunk_csv=str(chunk_csv),
                        chunk_root=str(chunk_root),
                        compact_lc_h5=str(args.compact_lc_h5),
                        hlsp_root=str(args.hlsp_root),
                        tensor_config=tensor_cfg.__dict__,
                    )
                    futures[future] = (chunk_idx, chunk_root, chunk_csv)
                for chunk_idx, chunk_root, chunk_csv, summary in ready:
                    score_tables.extend(
                        score_chunk(
                            chunk_idx=chunk_idx,
                            chunk_root=chunk_root,
                            chunk_csv=chunk_csv,
                            tensor_summary_chunk=summary,
                        )
                    )
                for future in as_completed(futures):
                    chunk_idx, chunk_root, chunk_csv = futures[future]
                    _, summary = future.result()
                    print(f"[eb-score] chunk={chunk_idx:04d} tensors_ready", flush=True)
                    score_tables.extend(
                        score_chunk(
                            chunk_idx=chunk_idx,
                            chunk_root=chunk_root,
                            chunk_csv=chunk_csv,
                            tensor_summary_chunk=summary,
                        )
                    )
        else:
            for chunk_idx, chunk_root, chunk_csv in chunk_infos:
                existing_paths = _chunk_score_paths(chunk_root, len(model_paths)) if args.resume_scoring else []
                if len(existing_paths) == len(model_paths):
                    score_tables.extend(existing_paths)
                    print(f"[eb-score] chunk={chunk_idx:04d} reused_complete_scores", flush=True)
                    continue
                tensor_summary_chunk = _tensor_summary_if_complete(chunk_root)
                if tensor_summary_chunk is None:
                    tensor_summary_chunk = build_recovery50_cnn_tensors(
                        training_table=chunk_csv,
                        out_dir=chunk_root / "tensors",
                        compact_lc_h5=args.compact_lc_h5,
                        hlsp_root=args.hlsp_root,
                        injection_h5_override=None,
                        config=tensor_cfg,
                        progress_every=100,
                    )
                score_tables.extend(
                    score_chunk(
                        chunk_idx=chunk_idx,
                        chunk_root=chunk_root,
                        chunk_csv=chunk_csv,
                        tensor_summary_chunk=tensor_summary_chunk,
                    )
                )
        ensemble = aggregate_ensemble_scores(score_tables)
        ensemble_path = write_table(ensemble, args.out_dir / "scored_candidates_ensemble.csv")
        scoring_summary = {
            "skipped": False,
            "candidate_pool": str(args.candidate_pool),
            "n_input_candidates": int(len(candidate_pool)),
            "n_score_tables": int(len(score_tables)),
            "n_ensemble_rows": int(len(ensemble)),
            "score_profile_policy": args.score_profile,
            "selected_profiles": selected_profiles,
            "score_chunk_workers": int(score_chunk_workers),
            "resume_scoring": bool(args.resume_scoring),
            "reuse_existing_models": bool(args.reuse_existing_models),
            "outputs": {"ensemble_scores": str(ensemble_path)},
        }

    top_summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "training_table": str(args.training_table),
        "candidate_pool": str(args.candidate_pool),
        "out_dir": str(args.out_dir),
        "tensor_summary": tensor_summary,
        "ensemble_size": int(args.ensemble_size),
        "score_profile_policy": args.score_profile,
        "selected_profiles": selected_profiles,
        "member_summaries": member_summaries,
        "scoring": scoring_summary,
    }
    (args.out_dir / "summary.json").write_text(
        json.dumps(top_summary, indent=2, sort_keys=True, default=json_default) + "\n"
    )
    print(json.dumps(top_summary, indent=2, sort_keys=True, default=json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
