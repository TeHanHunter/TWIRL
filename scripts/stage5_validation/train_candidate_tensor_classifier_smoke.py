#!/usr/bin/env python3
"""Train a tiny candidate-centered tensor classifier smoke model.

This is the first H200 *training path* smoke for TWIRL. It is intentionally
small and bounded: it consumes tensors from
``build_s56_candidate_tensor_smoke.py`` and teacher rows from the human-label
training table. If real teacher labels are not available, it can run an explicit
synthetic-label smoke that validates the GPU/data plumbing without producing a
scientific model.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_TENSOR_NPZ = (
    REPO_ROOT
    / "reports/stage5_validation/s56_candidate_tensor_smoke_h200_torch/candidate_tensor_smoke.npz"
)
DEFAULT_TENSOR_ROWS = DEFAULT_TENSOR_NPZ.with_name("candidate_tensor_smoke_rows.csv")
DEFAULT_LABELS = (
    REPO_ROOT
    / "reports/stage5_validation/s56_allhost_ranker_selected_real_leo_queue_pdo/"
    "human_training_table/teacher_labeled_rows.csv"
)
DEFAULT_OUT_DIR = REPO_ROOT / "reports/stage5_validation/s56_candidate_tensor_train_smoke_h200"


@dataclass
class TrainConfig:
    epochs: int = 20
    batch_size: int = 32
    learning_rate: float = 1.0e-3
    weight_decay: float = 1.0e-4
    validation_fraction: float = 0.20
    test_fraction: float = 0.20
    min_teacher_rows: int = 20
    seed: int = 56


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix == ".parquet":
        return pd.read_parquet(path)
    if suffix in {".json", ".jsonl"}:
        return pd.read_json(path, lines=suffix == ".jsonl")
    raise ValueError(f"unsupported table format: {path}")


def _load_tensor_product(tensor_npz: Path, tensor_rows: Path) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    if not tensor_npz.exists():
        raise FileNotFoundError(f"missing tensor NPZ: {tensor_npz}")
    if not tensor_rows.exists():
        raise FileNotFoundError(f"missing tensor row table: {tensor_rows}")
    with np.load(tensor_npz, allow_pickle=True) as npz:
        tensor = np.asarray(npz["tensor"], dtype=np.float32)
        mask = np.asarray(npz["mask"], dtype=np.float32) if "mask" in npz else np.isfinite(tensor).astype(np.float32)
    rows = _read_table(tensor_rows)
    if len(rows) != tensor.shape[0]:
        raise ValueError(f"tensor rows mismatch: {len(rows)} rows for tensor shape {tensor.shape}")
    return tensor, mask, rows


def _prepare_inputs(tensor: np.ndarray, mask: np.ndarray) -> np.ndarray:
    values = np.nan_to_num(tensor, nan=1.0, posinf=1.0, neginf=1.0)
    centered = values - 1.0
    return np.concatenate([centered, mask.astype(np.float32)], axis=1).astype(np.float32)


def _synthetic_labels(rows: pd.DataFrame) -> pd.DataFrame:
    out = rows.copy()
    metric = pd.to_numeric(out.get("tmag", pd.Series(np.arange(len(out)), index=out.index)), errors="coerce")
    if metric.notna().sum() < 2:
        metric = pd.Series(np.arange(len(out)), index=out.index, dtype=float)
    threshold = float(metric.median())
    out["teacher_target"] = np.where(metric <= threshold, "synthetic_bright", "synthetic_faint")
    out["training_split"] = ""
    return out.loc[:, ["review_id", "teacher_target", "training_split"]]


def _load_labels(
    *,
    rows: pd.DataFrame,
    labels: Path,
    label_column: str,
    synthetic_label_smoke: bool,
) -> tuple[pd.DataFrame, str]:
    if synthetic_label_smoke:
        return _synthetic_labels(rows), "synthetic_label_smoke"
    if not labels.exists() or labels.stat().st_size == 0:
        raise SystemExit(
            f"[tensor-train] missing teacher labels: {labels}\n"
            "Run the human-label audit first, or use --synthetic-label-smoke for an environment-only test."
        )
    label_table = _read_table(labels)
    if "review_id" not in label_table:
        raise KeyError(f"label table missing review_id: {labels}")
    if label_column not in label_table:
        raise KeyError(f"label table missing {label_column}: {labels}")
    keep = ["review_id", label_column]
    if "training_split" in label_table:
        keep.append("training_split")
    out = label_table.loc[:, keep].copy()
    if label_column != "teacher_target":
        out = out.rename(columns={label_column: "teacher_target"})
    if "training_split" not in out:
        out["training_split"] = ""
    return out, "human_teacher_labels"


def _build_training_rows(
    *,
    rows: pd.DataFrame,
    labels_df: pd.DataFrame,
    cfg: TrainConfig,
    label_mode: str,
) -> pd.DataFrame:
    work = rows.reset_index(names="tensor_index").merge(labels_df, on="review_id", how="inner")
    work["teacher_target"] = work["teacher_target"].fillna("").astype(str)
    work = work[work["teacher_target"].ne("")].copy()
    if len(work) < cfg.min_teacher_rows and label_mode != "synthetic_label_smoke":
        raise SystemExit(
            f"[tensor-train] only {len(work)} teacher rows joined; "
            f"need at least {cfg.min_teacher_rows} for a real training smoke"
        )
    if work["teacher_target"].nunique() < 2:
        raise SystemExit("[tensor-train] need at least two teacher classes")
    work["training_split"] = work.get("training_split", "").fillna("").astype(str)
    if work["training_split"].eq("").all() or work["training_split"].eq("unlabeled_or_audit").all():
        work["training_split"] = _assign_splits(work, cfg)
    return work.reset_index(drop=True)


def _assign_splits(work: pd.DataFrame, cfg: TrainConfig) -> pd.Series:
    rng = np.random.default_rng(cfg.seed)
    split = pd.Series("train", index=work.index, dtype=object)
    for _, idx_values in work.groupby("teacher_target").groups.items():
        idx = np.asarray(list(idx_values), dtype=int)
        rng.shuffle(idx)
        n = len(idx)
        if n < 3:
            continue
        n_test = max(1, int(round(cfg.test_fraction * n))) if cfg.test_fraction > 0 else 0
        n_val = max(1, int(round(cfg.validation_fraction * n))) if cfg.validation_fraction > 0 and n - n_test >= 3 else 0
        n_test = min(n_test, max(0, n - 1))
        n_val = min(n_val, max(0, n - n_test - 1))
        if n_test:
            split.iloc[idx[:n_test]] = "test"
        if n_val:
            split.iloc[idx[n_test:n_test + n_val]] = "validation"
    return split


def _accuracy(pred: np.ndarray, truth: np.ndarray) -> float:
    if truth.size == 0:
        return float("nan")
    return float(np.mean(pred == truth))


def train_tensor_classifier(
    *,
    tensor_npz: Path,
    tensor_rows: Path,
    labels: Path,
    out_dir: Path,
    label_column: str,
    synthetic_label_smoke: bool,
    require_cuda: bool,
    cfg: TrainConfig,
) -> dict[str, Any]:
    import torch
    from torch import nn
    from torch.utils.data import DataLoader, TensorDataset

    out_dir.mkdir(parents=True, exist_ok=True)
    torch.manual_seed(cfg.seed)
    np.random.seed(cfg.seed)

    tensor, mask, rows = _load_tensor_product(tensor_npz, tensor_rows)
    x_all = _prepare_inputs(tensor, mask)
    labels_df, label_mode = _load_labels(
        rows=rows,
        labels=labels,
        label_column=label_column,
        synthetic_label_smoke=synthetic_label_smoke,
    )
    train_rows = _build_training_rows(rows=rows, labels_df=labels_df, cfg=cfg, label_mode=label_mode)

    classes = tuple(sorted(train_rows["teacher_target"].unique()))
    class_to_id = {label: idx for idx, label in enumerate(classes)}
    y = train_rows["teacher_target"].map(class_to_id).to_numpy(dtype=np.int64)
    x = x_all[train_rows["tensor_index"].to_numpy(dtype=int)]

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if require_cuda and device.type != "cuda":
        raise RuntimeError("CUDA was required but torch.cuda.is_available() is false")

    model = nn.Sequential(
        nn.Conv1d(x.shape[1], 16, kernel_size=5, padding=2),
        nn.ReLU(),
        nn.Conv1d(16, 16, kernel_size=5, padding=2),
        nn.ReLU(),
        nn.AdaptiveAvgPool1d(1),
        nn.Flatten(),
        nn.Linear(16, len(classes)),
    ).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=cfg.learning_rate, weight_decay=cfg.weight_decay)
    loss_fn = nn.CrossEntropyLoss()

    split = train_rows["training_split"].fillna("train").astype(str).to_numpy()
    if not np.any(split == "train"):
        split[:] = "train"
    train_idx = np.where(split == "train")[0]
    train_data = TensorDataset(
        torch.from_numpy(x[train_idx]).float(),
        torch.from_numpy(y[train_idx]).long(),
    )
    generator = torch.Generator().manual_seed(cfg.seed)
    loader = DataLoader(train_data, batch_size=cfg.batch_size, shuffle=True, generator=generator)

    history: list[dict[str, float]] = []
    for epoch in range(1, cfg.epochs + 1):
        model.train()
        total_loss = 0.0
        n_seen = 0
        for xb, yb in loader:
            xb = xb.to(device)
            yb = yb.to(device)
            opt.zero_grad(set_to_none=True)
            loss = loss_fn(model(xb), yb)
            loss.backward()
            opt.step()
            total_loss += float(loss.detach().cpu()) * len(yb)
            n_seen += len(yb)
        metrics = {"epoch": float(epoch), "train_loss": total_loss / max(n_seen, 1)}
        metrics.update(_evaluate_model(model, x, y, split, device))
        history.append(metrics)
        print(
            "[tensor-train] "
            f"epoch={epoch:03d} loss={metrics['train_loss']:.5f} "
            f"train_acc={metrics.get('train_accuracy', float('nan')):.3f} "
            f"val_acc={metrics.get('validation_accuracy', float('nan')):.3f}",
            flush=True,
        )

    model.eval()
    with torch.no_grad():
        logits = model(torch.from_numpy(x).float().to(device)).cpu().numpy()
    pred_id = np.argmax(logits, axis=1)
    pred_label = np.asarray(classes, dtype=object)[pred_id]
    predictions = train_rows.copy()
    predictions["predicted_label"] = pred_label
    predictions["predicted_correct"] = predictions["predicted_label"].eq(predictions["teacher_target"])
    for idx, label in enumerate(classes):
        predictions[f"p_{label}"] = _softmax_np(logits)[:, idx]

    model_path = out_dir / "tensor_classifier_smoke.pt"
    torch.save(
        {
            "model_state_dict": model.state_dict(),
            "classes": classes,
            "input_channels": int(x.shape[1]),
            "config": asdict(cfg),
            "label_mode": label_mode,
        },
        model_path,
    )
    predictions.to_csv(out_dir / "predictions.csv", index=False)
    pd.DataFrame(history).to_csv(out_dir / "training_history.csv", index=False)

    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "tensor_npz": str(tensor_npz),
        "tensor_rows": str(tensor_rows),
        "labels": str(labels),
        "out_dir": str(out_dir),
        "label_mode": label_mode,
        "synthetic_label_smoke": bool(synthetic_label_smoke),
        "n_tensor_rows": int(tensor.shape[0]),
        "n_training_rows": int(len(train_rows)),
        "classes": list(classes),
        "class_counts": {str(k): int(v) for k, v in train_rows["teacher_target"].value_counts().sort_index().items()},
        "split_counts": {str(k): int(v) for k, v in train_rows["training_split"].value_counts().sort_index().items()},
        "device": str(device),
        "torch_version": str(torch.__version__),
        "torch_cuda_built": str(torch.version.cuda),
        "torch_cuda_available": bool(torch.cuda.is_available()),
        "torch_device_count": int(torch.cuda.device_count()) if torch.cuda.is_available() else 0,
        "final_metrics": history[-1] if history else {},
        "outputs": {
            "model": str(model_path),
            "predictions": str(out_dir / "predictions.csv"),
            "history": str(out_dir / "training_history.csv"),
            "summary": str(out_dir / "summary.json"),
        },
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n")
    (out_dir / "summary.md").write_text(_summary_markdown(summary))
    print(json.dumps(summary, indent=2, sort_keys=True, default=_json_default))
    return summary


def _evaluate_model(model: Any, x: np.ndarray, y: np.ndarray, split: np.ndarray, device: Any) -> dict[str, float]:
    import torch

    model.eval()
    with torch.no_grad():
        logits = model(torch.from_numpy(x).float().to(device)).cpu().numpy()
    pred = np.argmax(logits, axis=1)
    out: dict[str, float] = {}
    for name in ("train", "validation", "test"):
        idx = split == name
        out[f"{name}_accuracy"] = _accuracy(pred[idx], y[idx])
        out[f"{name}_n"] = float(np.sum(idx))
    return out


def _softmax_np(logits: np.ndarray) -> np.ndarray:
    shifted = logits - np.max(logits, axis=1, keepdims=True)
    exp = np.exp(shifted)
    return exp / np.sum(exp, axis=1, keepdims=True)


def _summary_markdown(summary: dict[str, Any]) -> str:
    lines = [
        "# S56 Candidate Tensor Training Smoke",
        "",
        f"- label mode: `{summary['label_mode']}`",
        f"- synthetic label smoke: `{summary['synthetic_label_smoke']}`",
        f"- training rows: `{summary['n_training_rows']}`",
        f"- classes: `{', '.join(summary['classes'])}`",
        f"- split counts: `{summary['split_counts']}`",
        f"- device: `{summary['device']}`",
        f"- torch: `{summary['torch_version']}` / CUDA `{summary['torch_cuda_built']}`",
        f"- final metrics: `{summary['final_metrics']}`",
        "",
        "This is a bounded infrastructure smoke. Synthetic-label runs are not scientific models.",
        "",
    ]
    return "\n".join(lines)


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tensor-npz", type=Path, default=DEFAULT_TENSOR_NPZ)
    parser.add_argument("--tensor-rows", type=Path, default=DEFAULT_TENSOR_ROWS)
    parser.add_argument("--labels", type=Path, default=DEFAULT_LABELS)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--label-column", default="teacher_target")
    parser.add_argument("--synthetic-label-smoke", action="store_true")
    parser.add_argument("--require-cuda", action="store_true")
    parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--learning-rate", type=float, default=1.0e-3)
    parser.add_argument("--weight-decay", type=float, default=1.0e-4)
    parser.add_argument("--validation-fraction", type=float, default=0.20)
    parser.add_argument("--test-fraction", type=float, default=0.20)
    parser.add_argument("--min-teacher-rows", type=int, default=20)
    parser.add_argument("--seed", type=int, default=56)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    cfg = TrainConfig(
        epochs=args.epochs,
        batch_size=args.batch_size,
        learning_rate=args.learning_rate,
        weight_decay=args.weight_decay,
        validation_fraction=args.validation_fraction,
        test_fraction=args.test_fraction,
        min_teacher_rows=args.min_teacher_rows,
        seed=args.seed,
    )
    train_tensor_classifier(
        tensor_npz=args.tensor_npz,
        tensor_rows=args.tensor_rows,
        labels=args.labels,
        out_dir=args.out_dir,
        label_column=args.label_column,
        synthetic_label_smoke=args.synthetic_label_smoke,
        require_cuda=args.require_cuda,
        cfg=cfg,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
