#!/usr/bin/env python3
"""Build a small candidate-centered light-curve tensor smoke product.

This is the data-preparation smoke for the later H200 light-curve model. It
does not train the production classifier. It verifies that a compact S56 LC
export plus ranker-selected ephemerides can be converted into fixed-shape,
multi-aperture tensors with explicit masks and provenance.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import importlib.util
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

try:
    from twirl.io.hlsp import BJDREFI  # noqa: E402
except ModuleNotFoundError:
    BJDREFI = 2457000.0


DEFAULT_LC_H5 = (
    REPO_ROOT
    / "data_local/stage3_injections/s56_twirlfs_v2_lc_export/s56_twirlfs_v2_lc_export_pdo.h5"
)
DEFAULT_CANDIDATES = (
    REPO_ROOT
    / "reports/stage5_validation/s56_allhost_ranker_selected_real_leo_queue_pdo/review_queue.csv"
)
DEFAULT_OUT_DIR = REPO_ROOT / "reports/stage5_validation/s56_candidate_tensor_smoke"
DEFAULT_APERTURES = ("DET_FLUX_SML", "DET_FLUX", "DET_FLUX_LAG")


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


def _safe_float(value: Any, default: float = float("nan")) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if np.isfinite(out) else default


def _time_to_bjd(time: np.ndarray) -> np.ndarray:
    finite = time[np.isfinite(time)]
    if finite.size and np.nanmedian(finite) < 1.0e5:
        return time + BJDREFI
    return time


def _bin_channel(
    *,
    time_bjd: np.ndarray,
    flux: np.ndarray,
    quality: np.ndarray,
    period_d: float,
    t0_bjd: float,
    duration_min: float,
    x_edges: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    duration_hr = duration_min / 60.0
    if duration_hr <= 0 or period_d <= 0:
        n = len(x_edges) - 1
        return np.full(n, np.nan, dtype=np.float32), np.zeros(n, dtype=np.bool_), np.zeros(n, dtype=np.int16)

    phase_hr = (((time_bjd - t0_bjd + 0.5 * period_d) % period_d) - 0.5 * period_d) * 24.0
    x = phase_hr / duration_hr
    good = (quality == 0) & np.isfinite(x) & np.isfinite(flux)
    if not np.any(good):
        n = len(x_edges) - 1
        return np.full(n, np.nan, dtype=np.float32), np.zeros(n, dtype=np.bool_), np.zeros(n, dtype=np.int16)

    x_good = x[good]
    flux_good = flux[good].astype(np.float64)
    oot = np.abs(x_good) > 1.5
    baseline = float(np.nanmedian(flux_good[oot])) if np.any(oot) else float(np.nanmedian(flux_good))
    if not np.isfinite(baseline) or abs(baseline) < 1.0e-8:
        baseline = 1.0
    y = flux_good / baseline

    bin_id = np.digitize(x_good, x_edges) - 1
    n_bins = len(x_edges) - 1
    out = np.full(n_bins, np.nan, dtype=np.float32)
    counts = np.zeros(n_bins, dtype=np.int16)
    for idx in range(n_bins):
        values = y[bin_id == idx]
        values = values[np.isfinite(values)]
        if values.size:
            out[idx] = float(np.nanmedian(values))
            counts[idx] = min(int(values.size), np.iinfo(np.int16).max)
    mask = np.isfinite(out)
    return out, mask, counts


def build_tensor_smoke(
    *,
    lc_h5: Path,
    candidate_table: Path,
    out_dir: Path,
    apertures: tuple[str, ...],
    max_candidates: int,
    n_points: int,
    window_durations: float,
    require_labels: bool,
    torch_smoke: bool,
    require_torch: bool,
    progress_every: int,
) -> dict[str, Any]:
    import h5py

    if not lc_h5.exists():
        raise FileNotFoundError(f"missing compact LC HDF5: {lc_h5}")
    if not candidate_table.exists():
        raise FileNotFoundError(f"missing candidate table: {candidate_table}")
    out_dir.mkdir(parents=True, exist_ok=True)

    candidates = _read_table(candidate_table).copy()
    if require_labels:
        labeled = candidates.get("human_label", candidates.get("label", pd.Series("", index=candidates.index)))
        candidates = candidates[labeled.fillna("").astype(str).ne("")].copy()
    if max_candidates > 0:
        candidates = candidates.head(max_candidates).copy()
    print(
        f"[tensor-smoke] candidate rows={len(candidates):,}; "
        f"apertures={','.join(apertures)}",
        flush=True,
    )

    x_grid = np.linspace(-window_durations, window_durations, n_points, dtype=np.float32)
    step = float(x_grid[1] - x_grid[0]) if n_points > 1 else float(window_durations)
    x_edges = np.concatenate(
        ([x_grid[0] - 0.5 * step], 0.5 * (x_grid[1:] + x_grid[:-1]), [x_grid[-1] + 0.5 * step])
    )

    tensors: list[np.ndarray] = []
    masks: list[np.ndarray] = []
    count_arrays: list[np.ndarray] = []
    records: list[dict[str, Any]] = []
    skipped = {"missing_lc": 0, "missing_aperture": 0, "invalid_ephemeris": 0}
    missing_apertures = {aperture: 0 for aperture in apertures}

    print(f"[tensor-smoke] opening compact LC HDF5: {lc_h5}", flush=True)
    with h5py.File(lc_h5, "r") as h5:
        targets = h5["targets"]
        print("[tensor-smoke] compact LC HDF5 opened", flush=True)
        for processed, (row_index, row) in enumerate(candidates.iterrows(), start=1):
            tic = int(row["tic"])
            key = f"{tic:016d}"
            if key not in targets:
                skipped["missing_lc"] += 1
                if progress_every > 0 and (processed % progress_every == 0 or processed == len(candidates)):
                    print(
                        f"[tensor-smoke] processed={processed:,}/{len(candidates):,}; "
                        f"accepted={len(tensors):,}; skipped={skipped}",
                        flush=True,
                    )
                continue
            group = targets[key]
            missing = [ap for ap in apertures if ap not in group]
            if missing:
                skipped["missing_aperture"] += 1
                for aperture in missing:
                    missing_apertures[aperture] += 1
                if progress_every > 0 and (processed % progress_every == 0 or processed == len(candidates)):
                    print(
                        f"[tensor-smoke] processed={processed:,}/{len(candidates):,}; "
                        f"accepted={len(tensors):,}; skipped={skipped}",
                        flush=True,
                    )
                continue
            period_d = _safe_float(row.get("period_d"))
            t0_bjd = _safe_float(row.get("t0_bjd"))
            duration_min = _safe_float(row.get("duration_min"))
            if not all(np.isfinite([period_d, t0_bjd, duration_min])) or period_d <= 0 or duration_min <= 0:
                skipped["invalid_ephemeris"] += 1
                if progress_every > 0 and (processed % progress_every == 0 or processed == len(candidates)):
                    print(
                        f"[tensor-smoke] processed={processed:,}/{len(candidates):,}; "
                        f"accepted={len(tensors):,}; skipped={skipped}",
                        flush=True,
                    )
                continue

            time_bjd = _time_to_bjd(np.asarray(group["time"], dtype=np.float64))
            quality = np.asarray(group["quality"], dtype=np.int32)
            channels = []
            channel_masks = []
            channel_counts = []
            for aperture in apertures:
                flux = np.asarray(group[aperture], dtype=np.float64)
                values, mask, counts = _bin_channel(
                    time_bjd=time_bjd,
                    flux=flux,
                    quality=quality,
                    period_d=period_d,
                    t0_bjd=t0_bjd,
                    duration_min=duration_min,
                    x_edges=x_edges,
                )
                channels.append(values)
                channel_masks.append(mask)
                channel_counts.append(counts)

            tensors.append(np.stack(channels, axis=0))
            masks.append(np.stack(channel_masks, axis=0))
            count_arrays.append(np.stack(channel_counts, axis=0))
            records.append(
                {
                    "row_index": int(row_index),
                    "review_id": str(row.get("review_id", "")),
                    "tic": tic,
                    "period_d": period_d,
                    "t0_bjd": t0_bjd,
                    "duration_min": duration_min,
                    "tmag": _safe_float(row.get("tmag")),
                    "source_kind": str(row.get("source_kind", "")),
                    "leo_class": str(row.get("leo_class", "")),
                }
            )
            if progress_every > 0 and (processed % progress_every == 0 or processed == len(candidates)):
                print(
                    f"[tensor-smoke] processed={processed:,}/{len(candidates):,}; "
                    f"accepted={len(tensors):,}; skipped={skipped}",
                    flush=True,
                )

    if tensors:
        tensor = np.stack(tensors, axis=0).astype(np.float32)
        mask = np.stack(masks, axis=0).astype(np.bool_)
        counts = np.stack(count_arrays, axis=0).astype(np.int16)
    else:
        tensor = np.empty((0, len(apertures), n_points), dtype=np.float32)
        mask = np.empty((0, len(apertures), n_points), dtype=np.bool_)
        counts = np.empty((0, len(apertures), n_points), dtype=np.int16)

    out_npz = out_dir / "candidate_tensor_smoke.npz"
    np.savez_compressed(
        out_npz,
        tensor=tensor,
        mask=mask,
        counts=counts,
        x_duration=x_grid,
        aperture=np.asarray(apertures, dtype=object),
        tic=np.asarray([r["tic"] for r in records], dtype=np.int64),
        row_index=np.asarray([r["row_index"] for r in records], dtype=np.int64),
    )
    record_table = pd.DataFrame(records)
    record_table.to_csv(out_dir / "candidate_tensor_smoke_rows.csv", index=False)

    torch_info: dict[str, Any] = {"requested": bool(torch_smoke), "available": False}
    if torch_smoke or require_torch:
        torch_info = _torch_smoke(tensor, require_torch=require_torch)

    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "lc_h5": str(lc_h5),
        "candidate_table": str(candidate_table),
        "out_npz": str(out_npz),
        "n_input_rows": int(len(candidates)),
        "n_tensor_rows": int(tensor.shape[0]),
        "shape": list(tensor.shape),
        "apertures": list(apertures),
        "n_points": int(n_points),
        "window_durations": float(window_durations),
        "finite_fraction": float(np.isfinite(tensor).sum() / tensor.size) if tensor.size else 0.0,
        "observed_fraction": float(mask.sum() / mask.size) if mask.size else 0.0,
        "skipped": skipped,
        "missing_apertures": missing_apertures,
        "torch": torch_info,
        "outputs": {
            "npz": str(out_npz),
            "rows": str(out_dir / "candidate_tensor_smoke_rows.csv"),
            "summary": str(out_dir / "summary.json"),
        },
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n")
    (out_dir / "summary.md").write_text(_summary_markdown(summary))
    return summary


def _torch_smoke(tensor: np.ndarray, *, require_torch: bool) -> dict[str, Any]:
    info: dict[str, Any] = {
        "requested": True,
        "available": importlib.util.find_spec("torch") is not None,
        "cuda_available": False,
        "device_count": 0,
    }
    if not info["available"]:
        if require_torch:
            raise RuntimeError("torch is required but is not importable")
        return info
    import torch

    info["version"] = str(torch.__version__)
    info["cuda_available"] = bool(torch.cuda.is_available())
    info["device_count"] = int(torch.cuda.device_count()) if torch.cuda.is_available() else 0
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    x = np.nan_to_num(tensor[: min(len(tensor), 32)], nan=1.0).astype(np.float32)
    if x.size:
        xt = torch.from_numpy(x).to(device)
        centered = xt - xt.mean(dim=-1, keepdim=True)
        score = centered.square().mean().sqrt()
        info["device"] = str(device)
        info["rms"] = float(score.detach().cpu())
    return info


def _summary_markdown(summary: dict[str, Any]) -> str:
    lines = [
        "# S56 Candidate Tensor Smoke",
        "",
        f"- input rows: `{summary['n_input_rows']}`",
        f"- tensor rows: `{summary['n_tensor_rows']}`",
        f"- shape: `{tuple(summary['shape'])}`",
        f"- apertures: `{', '.join(summary['apertures'])}`",
        f"- observed fraction: `{summary['observed_fraction']:.3f}`",
        f"- finite fraction: `{summary['finite_fraction']:.3f}`",
        f"- skipped: `{summary['skipped']}`",
        f"- missing apertures: `{summary['missing_apertures']}`",
        f"- torch: `{summary['torch']}`",
        "",
        "This is a data-shape and environment smoke only, not a trained model.",
        "",
    ]
    return "\n".join(lines)


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lc-h5", type=Path, default=DEFAULT_LC_H5)
    parser.add_argument("--candidate-table", type=Path, default=DEFAULT_CANDIDATES)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--apertures", nargs="+", default=list(DEFAULT_APERTURES))
    parser.add_argument("--max-candidates", type=int, default=128)
    parser.add_argument("--n-points", type=int, default=257)
    parser.add_argument("--window-durations", type=float, default=8.0)
    parser.add_argument("--require-labels", action="store_true")
    parser.add_argument("--torch-smoke", action="store_true")
    parser.add_argument("--require-torch", action="store_true")
    parser.add_argument(
        "--progress-every",
        type=int,
        default=25,
        help="Print a progress line every N candidate rows; use 0 to disable.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    summary = build_tensor_smoke(
        lc_h5=args.lc_h5,
        candidate_table=args.candidate_table,
        out_dir=args.out_dir,
        apertures=tuple(args.apertures),
        max_candidates=args.max_candidates,
        n_points=args.n_points,
        window_durations=args.window_durations,
        require_labels=args.require_labels,
        torch_smoke=args.torch_smoke,
        require_torch=args.require_torch,
        progress_every=args.progress_every,
    )
    print("[tensor-smoke] complete")
    print(f"  rows: {summary['n_tensor_rows']:,}/{summary['n_input_rows']:,}")
    print(f"  shape: {tuple(summary['shape'])}")
    print(f"  out: {summary['out_npz']}")
    print(f"  torch: {summary['torch']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
