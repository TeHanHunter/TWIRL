"""Lightweight HTTP vetting app for human light-curve labels."""
from __future__ import annotations

from dataclasses import dataclass
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from io import BytesIO
import json
import os
from pathlib import Path
import tempfile
from typing import Any
from urllib.parse import parse_qs, urlparse

import numpy as np
import pandas as pd

from twirl.io.hlsp import APERTURES, BJDREFI, read_hlsp
from twirl.search.injections import box_transit_mask
from twirl.vetting.label_schema import LABEL_BUTTONS, LABEL_KEY_ALIASES, LABEL_OPTIONS


REVIEW_APERTURES: tuple[str, ...] = (
    *APERTURES,
    "DET_FLUX_ADP",
    "DET_FLUX_ADP_SML",
    "DET_FLUX_ADP_LAG",
    "DET_FLUX_ADP015",
    "DET_FLUX_ADP015_SML",
    "DET_FLUX_ADP015_LAG",
)

DEFAULT_DISPLAY_COLUMNS: tuple[str, ...] = (
    "tic",
    "sector",
    "cam",
    "ccd",
    "tmag",
    "vet_class",
    "class_rank",
    "blind_rank",
    "period_d",
    "t0_bjd",
    "duration_min",
    "depth",
    "depth_snr",
    "sde_max",
    "rep_aperture",
    "n_apertures_agree",
    "apertures_agree",
    "n_in_transit",
    "n_oot_band",
)

PROVENANCE_DISPLAY_COLUMNS: frozenset[str] = frozenset(
    {
        "source_bucket",
        "source_kind",
        "recovery_status",
        "injection_id",
        "signal_family",
        "truth_source_kind",
        "truth_source_bucket",
    }
)


def tic_shard_path(hlsp_root: Path, tic: int) -> Path:
    """Return the standard 4x4-digit TIC shard directory under an HLSP root."""
    tic16 = f"{int(tic):016d}"
    return Path(hlsp_root).joinpath(tic16[0:4], tic16[4:8], tic16[8:12], tic16[12:16])


def find_hlsp_path(hlsp_root: Path, tic: int, sector: int | None = None) -> Path | None:
    """Find the HLSP FITS for ``tic`` using the standard shard path first."""
    tic16 = f"{int(tic):016d}"
    sector_part = f"s{int(sector):04d}-" if sector is not None and np.isfinite(sector) else "s"
    shard = tic_shard_path(hlsp_root, int(tic))
    patterns = [
        f"hlsp_*_tess_ffi_{sector_part}{tic16}_tess_v*_llc.fits",
        f"hlsp_*_tess_ffi_*{tic16}*_llc.fits",
    ]
    for pattern in patterns:
        matches = sorted(shard.glob(pattern))
        if matches:
            return matches[0]
    fallback = sorted(Path(hlsp_root).rglob(f"hlsp_*_tess_ffi_*{tic16}*_llc.fits"))
    return fallback[0] if fallback else None


def find_leo_report(leo_report_roots: tuple[Path, ...], tic: int) -> Path | None:
    """Find a pre-rendered LEO-Vetter PDF report for ``tic``."""
    tic10 = f"{int(tic):010d}"
    patterns = (f"*_tic{tic10}_*.pdf", f"*_tic{int(tic)}_*.pdf")
    for root in leo_report_roots:
        root = Path(root)
        if not root.exists():
            continue
        for pattern in patterns:
            matches = sorted(root.glob(pattern))
            if matches:
                return matches[0]
    return None


def find_leo_report_for_row(leo_report_roots: tuple[Path, ...], row: pd.Series) -> Path | None:
    """Find the row-specific LEO report when the queue records its filename."""
    name = _clean_value(row.get("leo_report_name")) or ""
    if name:
        for root in leo_report_roots:
            path = Path(root) / str(name)
            if path.exists():
                return path
    tic = _safe_float(row.get("tic"))
    return None if tic is None else find_leo_report(leo_report_roots, int(tic))


def find_twirl_vet_sheet_for_row(twirl_vet_roots: tuple[Path, ...], row: pd.Series) -> Path | None:
    """Find a pre-rendered TWIRL two-aperture PNG vet sheet for ``row``."""
    names: list[str] = []
    for col in ("twirl_vet_sheet_name", "twirl_vet_sheet"):
        value = _clean_value(row.get(col))
        if value:
            names.append(str(value))
    review_id = str(_clean_value(row.get("review_id")) or "")
    tic = _safe_float(row.get("tic"))
    if review_id:
        safe = review_id.replace("/", "_").replace(":", "_")
        names.append(f"{safe}_twirl_twoap_twirl_fs_v2_adp015q.png")
        names.append(f"{safe}_twirl_twoap_adpplus_0p20.png")
    if tic is not None:
        names.append(f"*tic{int(tic)}*_twirl_twoap_*.png")
    for root in twirl_vet_roots:
        root = Path(root)
        if not root.exists():
            continue
        for name in names:
            if "*" in name:
                matches = sorted(root.glob(name))
                if matches:
                    return matches[0]
            else:
                path = root / name
                if path.exists():
                    return path
    return None


def leo_class_from_report(path: Path | None) -> str:
    """Return the LEO filename class prefix (PC/FA/FP) when available."""
    if path is None:
        return ""
    prefix = Path(path).name.split("_", 1)[0]
    return prefix if prefix in {"PC", "FA", "FP"} else ""


def load_leo_metric_maps(candidates_path: Path) -> tuple[dict[str, dict[str, Any]], dict[int, dict[str, Any]]]:
    """Load lightweight LEO report metadata, when available next to the review queue."""
    metrics_path = Path(candidates_path).parent / "leo_metrics.csv"
    if not metrics_path.exists():
        return {}, {}
    wanted = {"review_id", "tic", "error", "plot_error", "leo_class", "leo_report_name", "leo_report_path"}
    try:
        metrics = pd.read_csv(metrics_path, usecols=lambda col: col in wanted)
    except Exception:
        return {}, {}

    by_review_id: dict[str, dict[str, Any]] = {}
    by_tic: dict[int, dict[str, Any]] = {}
    for rec in metrics.to_dict("records"):
        review_id = str(_clean_value(rec.get("review_id")) or "")
        if review_id:
            by_review_id[review_id] = rec
        tic = _safe_float(rec.get("tic"))
        if tic is not None:
            by_tic[int(tic)] = rec
    return by_review_id, by_tic


def _clean_value(value: Any) -> Any:
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except TypeError:
        pass
    if isinstance(value, np.generic):
        return value.item()
    return value


def _safe_float(value: Any) -> float | None:
    value = _clean_value(value)
    if value is None:
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if np.isfinite(out) else None


def _candidate_key(row: pd.Series) -> str:
    return "|".join(
        str(_clean_value(row.get(col, "")) or "")
        for col in ("tic", "sector", "period_d", "t0_bjd", "source_bucket")
    )


@dataclass
class CandidateStore:
    candidates_path: Path
    labels_out: Path
    shuffle_order: bool = False
    random_seed: int = 56
    unlabeled_first: bool = True

    def __post_init__(self) -> None:
        self.candidates_path = Path(self.candidates_path)
        self.labels_out = Path(self.labels_out)
        self.frame = pd.read_csv(self.candidates_path)
        self.frame.insert(0, "row_id", np.arange(len(self.frame), dtype=int))
        self.frame["candidate_key"] = self.frame.apply(_candidate_key, axis=1)
        self.labels = self._load_labels()
        self._apply_labels()
        self.review_order = self._build_review_order()

    def _build_review_order(self) -> np.ndarray:
        order = np.arange(len(self.frame), dtype=int)
        if not self.shuffle_order or len(order) <= 1:
            return order
        rng = np.random.default_rng(int(self.random_seed))
        if self.unlabeled_first:
            labels = self.frame["label"].fillna("").astype(str).to_numpy()
            unlabeled = order[labels == ""]
            labeled = order[labels != ""]
            rng.shuffle(unlabeled)
            rng.shuffle(labeled)
            order = np.concatenate([unlabeled, labeled])
        else:
            rng.shuffle(order)
        if len(order) > 1 and int(order[0]) == 0:
            order = np.roll(order, -1)
        return order

    def _load_labels(self) -> pd.DataFrame:
        if self.labels_out.exists():
            labels = pd.read_csv(self.labels_out)
            if "row_id" in labels.columns:
                labels["row_id"] = labels["row_id"].astype(int)
            return labels
        return pd.DataFrame(columns=["row_id", "candidate_key", "label", "label_source", "labeler", "notes", "updated_utc"])

    def _apply_labels(self) -> None:
        for col in ("label", "label_source", "labeler", "notes", "updated_utc"):
            if col not in self.frame.columns:
                self.frame[col] = ""
            else:
                self.frame[col] = self.frame[col].fillna("").astype(object)
        if self.labels.empty:
            return
        labels = self.labels.drop_duplicates("row_id", keep="last").set_index("row_id")
        for row_id, label_row in labels.iterrows():
            if row_id not in self.frame.index:
                continue
            for col in ("label", "label_source", "labeler", "notes", "updated_utc"):
                if col in label_row:
                    self.frame.loc[int(row_id), col] = label_row[col]

    @property
    def count(self) -> int:
        return int(len(self.frame))

    def row(self, index: int) -> pd.Series:
        if self.count == 0:
            raise IndexError("candidate table is empty")
        index = max(0, min(int(index), self.count - 1))
        return self.frame.iloc[int(self.review_order[index])]

    def current_index_for_row_id(self, row_id: int) -> int:
        hits = np.flatnonzero(self.frame["row_id"].to_numpy(dtype=int) == int(row_id))
        if hits.size == 0:
            raise KeyError(f"unknown row_id: {row_id}")
        return int(hits[0])

    def save_label(self, row_id: int, label: str, labeler: str, notes: str) -> dict[str, Any]:
        if label not in LABEL_OPTIONS and label != "":
            raise ValueError(f"unsupported label: {label}")
        index = self.current_index_for_row_id(row_id)
        row = self.frame.iloc[index]
        from datetime import datetime, timezone

        record = {
            "row_id": int(row_id),
            "candidate_key": str(row["candidate_key"]),
            "tic": int(row["tic"]) if _safe_float(row.get("tic")) is not None else "",
            "sector": int(row["sector"]) if _safe_float(row.get("sector")) is not None else "",
            "label": label,
            "label_source": "human" if label else "",
            "labeler": labeler,
            "notes": notes,
            "updated_utc": datetime.now(timezone.utc).isoformat(),
        }
        labels = self.labels[self.labels.get("row_id", pd.Series(dtype=int)).astype(str) != str(row_id)].copy()
        labels = pd.concat([labels, pd.DataFrame([record])], ignore_index=True)
        self.labels = labels
        for key, value in record.items():
            if key in self.frame.columns:
                self.frame.loc[index, key] = value
        self.labels_out.parent.mkdir(parents=True, exist_ok=True)
        self.labels.to_csv(self.labels_out, index=False)
        return record

    def summary(self) -> dict[str, Any]:
        labels = self.frame["label"].fillna("").astype(str)
        labeled = labels.ne("")
        return {
            "count": self.count,
            "labeled": int(labeled.sum()),
            "unlabeled": int((~labeled).sum()),
            "label_counts": labels[labeled].value_counts().sort_index().to_dict(),
            "candidates_path": str(self.candidates_path),
            "labels_out": str(self.labels_out),
        }


class LightCurveVettingApp:
    def __init__(
        self,
        *,
        candidates_path: Path,
        labels_out: Path,
        hlsp_root: Path,
        leo_report_roots: tuple[Path, ...] = (),
        twirl_vet_roots: tuple[Path, ...] = (),
        default_aperture: str = "DET_FLUX",
        labeler: str = "",
        shuffle_order: bool = False,
        random_seed: int = 56,
        unlabeled_first: bool = True,
    ) -> None:
        self.store = CandidateStore(
            candidates_path=Path(candidates_path),
            labels_out=Path(labels_out),
            shuffle_order=shuffle_order,
            random_seed=random_seed,
            unlabeled_first=unlabeled_first,
        )
        self.hlsp_root = Path(hlsp_root)
        self.leo_report_roots = tuple(Path(p) for p in leo_report_roots)
        self.twirl_vet_roots = tuple(Path(p) for p in twirl_vet_roots)
        self.default_aperture = default_aperture
        self.labeler = labeler
        self.leo_metrics_by_review_id, self.leo_metrics_by_tic = load_leo_metric_maps(self.store.candidates_path)

    def _leo_metric_for_row(self, row: pd.Series) -> dict[str, Any]:
        review_id = str(_clean_value(row.get("review_id")) or "")
        if review_id and review_id in self.leo_metrics_by_review_id:
            return self.leo_metrics_by_review_id[review_id]
        tic = _safe_float(row.get("tic"))
        if tic is not None:
            return self.leo_metrics_by_tic.get(int(tic), {})
        return {}

    def candidate_payload(self, index: int) -> dict[str, Any]:
        row = self.store.row(index)
        label = _clean_value(row.get("label")) or ""
        display = {}
        for col in DEFAULT_DISPLAY_COLUMNS:
            if col not in row.index or col in PROVENANCE_DISPLAY_COLUMNS:
                continue
            value = _clean_value(row.get(col))
            if col == "vet_class" and isinstance(value, str) and value.startswith("injected_") and not label:
                value = "review_candidate"
            display[col] = value
        sector = _safe_float(row.get("sector"))
        path = None
        if _safe_float(row.get("tic")) is not None:
            path = find_hlsp_path(
                self.hlsp_root,
                int(float(row["tic"])),
                int(sector) if sector is not None else None,
            )
        leo_report = find_leo_report_for_row(self.leo_report_roots, row)
        twirl_vet_sheet = find_twirl_vet_sheet_for_row(self.twirl_vet_roots, row)
        leo_metric = self._leo_metric_for_row(row)
        leo_plot_error = _clean_value(leo_metric.get("plot_error")) or ""
        leo_error = _clean_value(leo_metric.get("error")) or ""
        leo_report_kind = ""
        if leo_report is not None:
            leo_report_kind = "fallback_plot" if leo_plot_error else "leo_vetter_report"
        return {
            "index": max(0, min(int(index), self.store.count - 1)),
            "count": self.store.count,
            "row_id": int(row["row_id"]),
            "candidate_key": str(row["candidate_key"]),
            "display": display,
            "label": label,
            "label_source": (_clean_value(row.get("label_source")) or "") if label else "",
            "labeler": _clean_value(row.get("labeler")) or self.labeler,
            "notes": _clean_value(row.get("notes")) or "",
            "updated_utc": _clean_value(row.get("updated_utc")) or "",
            "hlsp_path": str(path) if path else None,
            "leo_report_path": str(leo_report) if leo_report else None,
            "leo_report_name": leo_report.name if leo_report else "",
            "twirl_vet_sheet_path": str(twirl_vet_sheet) if twirl_vet_sheet else None,
            "twirl_vet_sheet_name": twirl_vet_sheet.name if twirl_vet_sheet else "",
            "leo_class": leo_class_from_report(leo_report),
            "leo_report_kind": leo_report_kind,
            "leo_plot_error": leo_plot_error,
            "leo_error": leo_error,
            "summary": self.store.summary(),
            "label_options": LABEL_OPTIONS,
            "apertures": REVIEW_APERTURES,
            "default_aperture": self.default_aperture,
        }

    def plot_png(self, index: int, aperture: str | None = None) -> bytes:
        mpl_cache = Path(tempfile.gettempdir()) / "twirl_mplconfig"
        mpl_cache.mkdir(parents=True, exist_ok=True)
        xdg_cache = Path(tempfile.gettempdir()) / "twirl_cache"
        xdg_cache.mkdir(parents=True, exist_ok=True)
        os.environ.setdefault("MPLCONFIGDIR", str(mpl_cache))
        os.environ.setdefault("XDG_CACHE_HOME", str(xdg_cache))
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        aperture = aperture or self.default_aperture
        row = self.store.row(index)
        tic = _safe_float(row.get("tic"))
        sector = _safe_float(row.get("sector"))
        path = None if tic is None else find_hlsp_path(
            self.hlsp_root,
            int(tic),
            int(sector) if sector is not None else None,
        )

        fig, axes = plt.subplots(2, 1, figsize=(10.5, 6.0), constrained_layout=True)
        if path is None:
            for ax in axes:
                ax.axis("off")
            axes[0].text(0.5, 0.5, f"No HLSP FITS found for TIC {tic}", ha="center", va="center")
            return _figure_png(fig)

        lc = read_hlsp(path, columns=(aperture,))
        if lc is None or aperture not in lc.flux:
            for ax in axes:
                ax.axis("off")
            axes[0].text(0.5, 0.5, f"Could not read {aperture} from {path}", ha="center", va="center")
            return _figure_png(fig)

        time = np.asarray(lc.time, dtype=float)
        flux = np.asarray(lc.flux[aperture], dtype=float)
        quality = np.asarray(lc.quality)
        good = (quality == 0) & np.isfinite(time) & np.isfinite(flux)
        flagged = (quality != 0) & np.isfinite(time) & np.isfinite(flux)

        period = _safe_float(row.get("period_d"))
        t0_bjd = _safe_float(row.get("t0_bjd"))
        duration = _safe_float(row.get("duration_min"))
        t0_d = None if t0_bjd is None else t0_bjd - BJDREFI
        in_transit = (
            box_transit_mask(time, period_d=period, t0_d=t0_d, duration_min=duration)
            if period and t0_d is not None and duration
            else np.zeros_like(time, dtype=bool)
        )

        ax = axes[0]
        ax.scatter(time[good], flux[good], s=4, c="#1b1f24", alpha=0.75, linewidths=0, label="QUALITY=0")
        if np.any(flagged):
            ax.scatter(time[flagged], flux[flagged], s=6, c="#c43c2f", alpha=0.75, linewidths=0, label="flagged")
        if np.any(in_transit):
            ax.scatter(time[in_transit & good], flux[in_transit & good], s=18, facecolors="none", edgecolors="#1f77b4", linewidths=0.8, label="ephemeris")
        ax.set_xlabel(f"Time (BJD - {BJDREFI})")
        ax.set_ylabel(aperture)
        ax.set_title(
            f"TIC {int(lc.tic)}  S{int(lc.sector)} cam{int(lc.cam)}/ccd{int(lc.ccd)}  T={lc.tmag:.2f}"
        )
        ax.legend(loc="best", fontsize=8)

        ax = axes[1]
        if period and t0_d is not None:
            phase_d = ((time - t0_d + 0.5 * period) % period) - 0.5 * period
            phase_hr = phase_d * 24.0
            keep = np.abs(phase_hr) <= max(2.0, (duration or 10.0) / 60.0 * 8.0)
            ax.scatter(phase_hr[good & keep], flux[good & keep], s=6, c="#1b1f24", alpha=0.75, linewidths=0)
            if np.any(flagged & keep):
                ax.scatter(phase_hr[flagged & keep], flux[flagged & keep], s=8, c="#c43c2f", alpha=0.75, linewidths=0)
            if duration:
                half_hr = duration / 120.0
                ax.axvspan(-half_hr, half_hr, color="#1f77b4", alpha=0.15)
            ax.set_xlabel("Phase from candidate epoch (hr)")
            ax.set_ylabel(aperture)
        else:
            ax.axis("off")
            ax.text(0.5, 0.5, "No ephemeris columns available", ha="center", va="center")

        for ax in axes:
            ax.grid(alpha=0.2)

        return _figure_png(fig)

    def make_handler(self) -> type[BaseHTTPRequestHandler]:
        app = self

        class Handler(BaseHTTPRequestHandler):
            def do_GET(self) -> None:  # noqa: N802
                parsed = urlparse(self.path)
                query = parse_qs(parsed.query)
                try:
                    if parsed.path == "/":
                        self._send_html(_index_html())
                    elif parsed.path == "/api/candidate":
                        index = int(query.get("index", ["0"])[0])
                        self._send_json(app.candidate_payload(index))
                    elif parsed.path == "/api/summary":
                        self._send_json(app.store.summary())
                    elif parsed.path == "/plot.png":
                        index = int(query.get("index", ["0"])[0])
                        aperture = query.get("aperture", [app.default_aperture])[0]
                        self._send_png(app.plot_png(index, aperture))
                    elif parsed.path == "/leo_report.pdf":
                        index = int(query.get("index", ["0"])[0])
                        row = app.store.row(index)
                        report = find_leo_report_for_row(app.leo_report_roots, row)
                        if report is None:
                            self.send_error(HTTPStatus.NOT_FOUND, "LEO report not found")
                        else:
                            self._send_file(report, "application/pdf", cache_seconds=3600)
                    elif parsed.path == "/twirl_vet_sheet.png":
                        index = int(query.get("index", ["0"])[0])
                        row = app.store.row(index)
                        sheet = find_twirl_vet_sheet_for_row(app.twirl_vet_roots, row)
                        if sheet is None:
                            self.send_error(HTTPStatus.NOT_FOUND, "TWIRL vet sheet not found")
                        else:
                            self._send_file(sheet, "image/png", cache_seconds=3600)
                    elif parsed.path == "/labels.csv":
                        self._send_file(app.store.labels_out, "text/csv")
                    else:
                        self.send_error(HTTPStatus.NOT_FOUND, "not found")
                except Exception as exc:
                    self._send_json({"error": str(exc)}, status=HTTPStatus.INTERNAL_SERVER_ERROR)

            def do_POST(self) -> None:  # noqa: N802
                parsed = urlparse(self.path)
                if parsed.path != "/api/label":
                    self.send_error(HTTPStatus.NOT_FOUND, "not found")
                    return
                try:
                    length = int(self.headers.get("Content-Length", "0"))
                    payload = json.loads(self.rfile.read(length).decode("utf-8"))
                    record = app.store.save_label(
                        row_id=int(payload.get("row_id")),
                        label=str(payload.get("label", "")),
                        labeler=str(payload.get("labeler", "")),
                        notes=str(payload.get("notes", "")),
                    )
                    self._send_json({"ok": True, "record": record, "summary": app.store.summary()})
                except Exception as exc:
                    self._send_json({"ok": False, "error": str(exc)}, status=HTTPStatus.BAD_REQUEST)

            def log_message(self, fmt: str, *args: Any) -> None:
                return

            def _send_json(self, payload: dict[str, Any], status: HTTPStatus = HTTPStatus.OK) -> None:
                body = json.dumps(payload, default=str).encode("utf-8")
                self.send_response(int(status))
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)

            def _send_html(self, html: str) -> None:
                body = html.encode("utf-8")
                self.send_response(HTTPStatus.OK)
                self.send_header("Content-Type", "text/html; charset=utf-8")
                self.send_header("Cache-Control", "no-store")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)

            def _send_png(self, body: bytes) -> None:
                self.send_response(HTTPStatus.OK)
                self.send_header("Content-Type", "image/png")
                self.send_header("Cache-Control", "no-store")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)

            def _send_file(self, path: Path, content_type: str, cache_seconds: int = 0) -> None:
                if not path.exists():
                    self.send_error(HTTPStatus.NOT_FOUND, "file not found")
                    return
                body = path.read_bytes()
                self.send_response(HTTPStatus.OK)
                self.send_header("Content-Type", content_type)
                if cache_seconds > 0:
                    self.send_header("Cache-Control", f"private, max-age={int(cache_seconds)}")
                else:
                    self.send_header("Cache-Control", "no-store")
                self.send_header("Content-Length", str(len(body)))
                filename = "leo_report.pdf" if content_type == "application/pdf" else path.name
                self.send_header("Content-Disposition", f'inline; filename="{filename}"')
                self.end_headers()
                self.wfile.write(body)

        return Handler

    def serve(self, host: str = "127.0.0.1", port: int = 5000) -> None:
        server = ThreadingHTTPServer((host, int(port)), self.make_handler())
        print(f"[vet-app] serving http://{host}:{port}/", flush=True)
        print(f"[vet-app] candidates: {self.store.candidates_path}", flush=True)
        print(f"[vet-app] labels: {self.store.labels_out}", flush=True)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass
        finally:
            server.server_close()


def _figure_png(fig) -> bytes:
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=150)
    import matplotlib.pyplot as plt

    plt.close(fig)
    return buf.getvalue()


def _index_html() -> str:
    labels = "".join(
        (
            f'<button class="label" data-label="{label}" data-shortcut="{key}" '
            f'title="{key}: {text}" aria-keyshortcuts="{key}">'
            f'<span class="key">{key}</span><span>{text}</span></button>'
        )
        for key, label, text in LABEL_BUTTONS
    )
    apertures = "".join(f'<option value="{ap}">{ap}</option>' for ap in REVIEW_APERTURES)
    shortcut_labels = {key: label for key, label, _ in LABEL_BUTTONS}
    shortcut_labels.update(LABEL_KEY_ALIASES)
    shortcut_json = json.dumps(shortcut_labels, sort_keys=True)
    return f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>TWIRL S56 vetter</title>
<style>
body {{ margin: 0; font-family: -apple-system, BlinkMacSystemFont, Segoe UI, sans-serif; color: #17202a; background: #f6f7f9; }}
header {{ display: grid; grid-template-columns: minmax(360px, 1fr) auto; gap: 8px 14px; align-items: center; padding: 9px 14px; background: #ffffff; border-bottom: 1px solid #d9dee7; position: sticky; top: 0; z-index: 2; box-shadow: 0 1px 8px rgba(23, 32, 42, 0.06); }}
button, select, input, textarea {{ font: inherit; }}
button {{ border: 1px solid #bac3d0; background: #fff; border-radius: 6px; padding: 6px 10px; cursor: pointer; }}
button:hover {{ background: #eef3f8; }}
button:disabled {{ opacity: 0.55; cursor: wait; }}
button.active {{ background: #17324d; color: #fff; border-color: #17324d; }}
.toolbar {{ display: flex; align-items: center; gap: 8px; min-width: 0; }}
.toolbar input[type="number"] {{ width: 82px; }}
.quick-labels {{ display: flex; flex-wrap: wrap; gap: 6px; justify-content: flex-end; }}
button.label {{ display: inline-flex; align-items: center; gap: 6px; min-height: 34px; padding: 5px 9px 5px 6px; border-color: #aeb8c5; }}
button.label .key {{ display: inline-grid; place-items: center; min-width: 22px; height: 22px; padding: 0 2px; border-radius: 4px; background: #edf1f5; color: #233040; font-weight: 700; font-size: 12px; }}
button.label.active .key {{ background: rgba(255,255,255,0.18); color: #fff; }}
main {{ display: grid; grid-template-columns: minmax(480px, 1fr) 390px; gap: 12px; padding: 12px; }}
.panel {{ background: #fff; border: 1px solid #d9dee7; border-radius: 8px; padding: 12px; }}
#plot {{ width: 100%; min-height: 420px; border: 1px solid #e1e5eb; border-radius: 6px; background: #fff; display: block; }}
#plotMessage {{ margin-top: 8px; color: #9b2f25; font-size: 13px; min-height: 18px; }}
#leoFrame {{ width: 100%; height: min(78vh, 900px); border: 1px solid #d9dee7; border-radius: 6px; background: #fff; }}
#leoEmpty {{ display: none; min-height: 180px; align-items: center; justify-content: center; color: #607080; border: 1px dashed #bac3d0; border-radius: 6px; }}
#twirlVetImage {{ width: 100%; border: 1px solid #d9dee7; border-radius: 6px; background: #fff; display: none; }}
#preloadBin {{ position: fixed; left: -10000px; top: -10000px; width: 1px; height: 1px; overflow: hidden; opacity: 0; pointer-events: none; }}
#preloadBin iframe {{ width: 1px; height: 1px; border: 0; }}
.lc-panel {{ grid-column: 1 / -1; }}
.lc-panel summary {{ cursor: pointer; font-weight: 600; margin-bottom: 8px; }}
.meta {{ display: grid; grid-template-columns: 150px 1fr; gap: 4px 8px; font-size: 13px; }}
.meta div:nth-child(odd) {{ color: #536170; }}
.notes {{ border-top: 1px solid #e1e5eb; margin-top: 14px; padding-top: 12px; }}
.notes-row {{ display: grid; grid-template-columns: minmax(0, 1fr) auto; gap: 8px; align-items: end; margin-top: 8px; }}
textarea {{ width: 100%; min-height: 90px; resize: vertical; box-sizing: border-box; }}
.muted {{ color: #607080; font-size: 13px; }}
.status {{ color: #536170; font-size: 13px; white-space: nowrap; }}
@media (max-width: 1100px) {{ header {{ grid-template-columns: 1fr; }} .quick-labels {{ justify-content: flex-start; }} main {{ grid-template-columns: 1fr; }} .notes-row {{ grid-template-columns: 1fr; }} }}
</style>
</head>
<body tabindex="-1">
<header>
  <div class="toolbar">
    <button id="prev">Prev</button>
    <input id="idx" type="number" min="0" value="0">
    <button id="next">Next</button>
    <select id="aperture">{apertures}</select>
    <button id="reload">Reload</button>
    <a href="/labels.csv" target="_blank">labels.csv</a>
    <span class="muted" id="summary"></span>
    <span class="status" id="status"></span>
  </div>
  <div class="quick-labels">{labels}</div>
</header>
<main>
  <section class="panel">
    <div style="display: flex; align-items: center; gap: 10px; margin-bottom: 8px;">
      <h2 id="reportTitle" style="margin: 0; font-size: 16px;">TWIRL Vet Sheet</h2>
      <a id="leoLink" target="_blank" style="display: none;">open report</a>
      <span class="muted" id="leoStatus"></span>
    </div>
    <img id="twirlVetImage" alt="TWIRL two-aperture vet sheet">
    <iframe id="leoFrame" title="LEO-Vetter report"></iframe>
    <div id="leoEmpty">No pre-rendered TWIRL or LEO report for this candidate yet.</div>
  </section>
  <aside class="panel">
    <div class="muted" id="counter"></div>
    <h2 id="title" style="margin: 6px 0 10px 0; font-size: 18px;"></h2>
    <div class="meta" id="meta"></div>
    <div class="notes">
      <textarea id="notes" placeholder="notes"></textarea>
      <div class="notes-row">
        <input id="labeler" placeholder="labeler" style="width: 100%; box-sizing: border-box;">
        <button id="save">Save notes</button>
      </div>
    </div>
  </aside>
  <details class="panel lc-panel" id="lcDetails">
    <summary>TWIRL Light Curve</summary>
    <img id="plot" alt="light curve">
    <div id="plotMessage"></div>
  </details>
</main>
<div id="preloadBin" aria-hidden="true"></div>
<script>
const shortcutLabels = {shortcut_json};
const preloadAhead = 4;
const preloadedPdfUrls = new Set();
const preloadedCandidateIndices = new Set();
const preloadFrames = new Map();
let state = {{ index: 0, row_id: null, label: "" }};
let isSaving = false;
function esc(x) {{ return String(x ?? "").replace(/[&<>"']/g, c => ({{'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}}[c])); }}
function currentAperture() {{
  return document.getElementById("aperture").value || state.default_aperture || "DET_FLUX";
}}
function setActiveLabel(label) {{
  state.label = label || "";
  document.querySelectorAll("button.label").forEach(b => b.classList.toggle("active", b.dataset.label === state.label));
}}
function setSavingUi(saving) {{
  isSaving = saving;
  document.querySelectorAll("button.label").forEach(b => b.disabled = saving);
  document.getElementById("save").disabled = saving;
}}
function showError(err) {{
  document.getElementById("status").textContent = err.message || String(err);
}}
function isTextEntry(target) {{
  if (!target) return false;
  const tag = target.tagName;
  return target.isContentEditable || tag === "INPUT" || tag === "TEXTAREA" || tag === "SELECT";
}}
function leoReportUrl(index) {{
  return `/leo_report.pdf?index=${{index}}`;
}}
function twirlVetUrl(index) {{
  return `/twirl_vet_sheet.png?index=${{index}}`;
}}
function preloadPdf(url) {{
  if (preloadedPdfUrls.has(url)) return;
  preloadedPdfUrls.add(url);
  const link = document.createElement("link");
  link.rel = "prefetch";
  link.as = "document";
  link.href = url;
  document.head.appendChild(link);
  fetch(url, {{ cache: "force-cache" }}).catch(() => preloadedPdfUrls.delete(url));
  const bin = document.getElementById("preloadBin");
  if (!bin) return;
  const frame = document.createElement("iframe");
  frame.src = url;
  frame.loading = "eager";
  frame.tabIndex = -1;
  frame.setAttribute("aria-hidden", "true");
  bin.appendChild(frame);
  preloadFrames.set(url, frame);
  while (preloadFrames.size > preloadAhead) {{
    const oldestUrl = preloadFrames.keys().next().value;
    const oldestFrame = preloadFrames.get(oldestUrl);
    if (oldestFrame) oldestFrame.remove();
    preloadFrames.delete(oldestUrl);
  }}
}}
function preloadNextCandidates(index) {{
  for (let offset = 1; offset <= preloadAhead; offset += 1) {{
    const nextIndex = Number(index) + offset;
    if (!Number.isFinite(nextIndex) || nextIndex >= (state.count || 0) || preloadedCandidateIndices.has(nextIndex)) continue;
    preloadedCandidateIndices.add(nextIndex);
    fetch(`/api/candidate?index=${{nextIndex}}`, {{ cache: "no-store" }})
      .then(resp => resp.ok ? resp.json() : null)
      .then(data => {{
        if (data && data.twirl_vet_sheet_path) {{
          const img = new Image();
          img.src = twirlVetUrl(data.index);
        }} else if (data && data.leo_report_path) {{
          preloadPdf(leoReportUrl(data.index));
        }}
      }})
      .catch(() => preloadedCandidateIndices.delete(nextIndex));
  }}
}}
async function loadPlot(index, aperture) {{
  const img = document.getElementById("plot");
  const msg = document.getElementById("plotMessage");
  msg.textContent = "";
  const resp = await fetch(`/plot.png?index=${{index}}&aperture=${{encodeURIComponent(aperture)}}&t=${{Date.now()}}`);
  if (!resp.ok) {{
    let detail = `${{resp.status}} ${{resp.statusText}}`;
    try {{
      const data = await resp.json();
      if (data.error) detail = data.error;
    }} catch (_) {{}}
    img.removeAttribute("src");
    msg.textContent = `Plot failed: ${{detail}}`;
    return;
  }}
  const blob = await resp.blob();
  if (!blob.type.startsWith("image/")) {{
    img.removeAttribute("src");
    msg.textContent = `Plot failed: unexpected content type ${{blob.type || "unknown"}}`;
    return;
  }}
  img.src = URL.createObjectURL(blob);
}}
async function loadCandidate(index) {{
  const resp = await fetch(`/api/candidate?index=${{index}}`);
  const data = await resp.json();
  if (data.error) throw new Error(data.error);
  state = data;
  document.getElementById("idx").value = data.index;
  document.getElementById("counter").textContent = `${{data.index + 1}} / ${{data.count}}`;
  document.getElementById("title").textContent = `TIC ${{data.display.tic ?? ""}}`;
  document.getElementById("labeler").value = data.labeler || "";
  document.getElementById("notes").value = data.notes || "";
  setActiveLabel(data.label || "");
  const rows = Object.entries(data.display).map(([k, v]) => `<div>${{esc(k)}}</div><div>${{esc(v)}}</div>`).join("");
  document.getElementById("meta").innerHTML = rows
    + `<div>leo_class</div><div>${{esc(data.leo_class || "")}}</div>`
    + (data.leo_plot_error ? `<div>leo_plot_error</div><div>${{esc(data.leo_plot_error)}}</div>` : "")
    + (data.leo_error ? `<div>leo_error</div><div>${{esc(data.leo_error)}}</div>` : "")
    + `<div>hlsp_path</div><div>${{esc(data.hlsp_path || "")}}</div>`;
  document.getElementById("summary").textContent = `${{data.summary.labeled}} labeled, ${{data.summary.unlabeled}} open`;
  const leoFrame = document.getElementById("leoFrame");
  const leoEmpty = document.getElementById("leoEmpty");
  const leoLink = document.getElementById("leoLink");
  const leoStatus = document.getElementById("leoStatus");
  const reportTitle = document.getElementById("reportTitle");
  const twirlVetImage = document.getElementById("twirlVetImage");
  if (data.twirl_vet_sheet_path) {{
    const url = twirlVetUrl(data.index);
    twirlVetImage.style.display = "block";
    twirlVetImage.src = url;
    leoFrame.removeAttribute("src");
    leoFrame.style.display = "none";
    leoEmpty.style.display = "none";
    leoLink.href = url;
    leoLink.textContent = "open PNG";
    leoLink.style.display = "inline";
    reportTitle.textContent = "TWIRL Two-Aperture Vet Sheet";
    leoStatus.textContent = data.twirl_vet_sheet_name || "";
  }} else if (data.leo_report_path) {{
    const url = leoReportUrl(data.index);
    twirlVetImage.removeAttribute("src");
    twirlVetImage.style.display = "none";
    leoFrame.style.display = "block";
    leoEmpty.style.display = "none";
    leoFrame.src = url;
    leoLink.href = url;
    leoLink.textContent = "open PDF";
    leoLink.style.display = "inline";
    if (data.leo_report_kind === "fallback_plot") {{
      reportTitle.textContent = "LEO Fallback Plot";
      leoStatus.textContent = data.leo_class ? `class: ${{data.leo_class}} | full plot failed` : "full plot failed";
    }} else {{
      reportTitle.textContent = "LEO-Vetter Report";
      leoStatus.textContent = data.leo_class ? `class: ${{data.leo_class}}` : "";
    }}
  }} else {{
    twirlVetImage.removeAttribute("src");
    twirlVetImage.style.display = "none";
    leoFrame.removeAttribute("src");
    leoFrame.style.display = "none";
    leoEmpty.style.display = "flex";
    leoLink.style.display = "none";
    reportTitle.textContent = "TWIRL Vet Sheet";
    leoStatus.textContent = "";
  }}
  const apertureSelect = document.getElementById("aperture");
  if (!apertureSelect.dataset.initialized) {{
    apertureSelect.value = data.default_aperture;
    apertureSelect.dataset.initialized = "1";
  }}
  const lcDetails = document.getElementById("lcDetails");
  document.getElementById("plot").removeAttribute("src");
  document.getElementById("plotMessage").textContent = "";
  lcDetails.open = !(data.twirl_vet_sheet_path || data.leo_report_path);
  if (lcDetails.open) {{
    await loadPlot(data.index, currentAperture());
  }}
  document.getElementById("status").textContent = data.label ? `saved: ${{data.label}}` : "";
  preloadNextCandidates(data.index);
  document.body.focus({{ preventScroll: true }});
}}
async function saveLabel(options = {{}}) {{
  if (isSaving) return;
  const label = options.label ?? state.label ?? "";
  const advance = Boolean(options.advance);
  const currentIndex = Number(state.index || 0);
  const payload = {{
    row_id: state.row_id,
    label,
    labeler: document.getElementById("labeler").value,
    notes: document.getElementById("notes").value
  }};
  setSavingUi(true);
  document.getElementById("status").textContent = "saving...";
  try {{
    const resp = await fetch("/api/label", {{ method: "POST", headers: {{ "Content-Type": "application/json" }}, body: JSON.stringify(payload) }});
    const data = await resp.json();
    if (!data.ok) throw new Error(data.error || "save failed");
    const nextIndex = advance ? Math.min(state.count - 1, currentIndex + 1) : currentIndex;
    await loadCandidate(nextIndex);
    const action = advance && nextIndex !== currentIndex ? "advanced" : "saved";
    document.getElementById("status").textContent = payload.label ? `${{action}}: ${{payload.label}}` : "saved notes";
  }} finally {{
    setSavingUi(false);
  }}
}}
async function labelAndNext(label) {{
  setActiveLabel(label);
  await saveLabel({{ label, advance: true }});
}}
document.getElementById("prev").onclick = () => loadCandidate(Math.max(0, state.index - 1));
document.getElementById("next").onclick = () => loadCandidate(Math.min(state.count - 1, state.index + 1));
document.getElementById("reload").onclick = () => loadCandidate(Number(document.getElementById("idx").value || 0));
document.getElementById("idx").onchange = e => loadCandidate(Number(e.target.value || 0));
document.getElementById("aperture").onchange = () => {{
  if (document.getElementById("lcDetails").open) loadPlot(state.index, currentAperture());
}};
document.getElementById("save").onclick = () => saveLabel({{ advance: false }}).catch(showError);
document.getElementById("lcDetails").addEventListener("toggle", () => {{
  if (document.getElementById("lcDetails").open && state.count) {{
    loadPlot(state.index, currentAperture()).catch(err => document.getElementById("plotMessage").textContent = err.message);
  }}
}});
document.querySelectorAll("button.label").forEach(b => b.onclick = () => labelAndNext(b.dataset.label).catch(showError));
document.addEventListener("keydown", event => {{
  if (event.defaultPrevented || event.metaKey || event.ctrlKey || event.altKey || isSaving || isTextEntry(event.target)) return;
  const key = event.key.toLowerCase();
  const label = shortcutLabels[key];
  if (label) {{
    event.preventDefault();
    labelAndNext(label).catch(showError);
    return;
  }}
  if (key === "arrowleft") {{
    event.preventDefault();
    loadCandidate(Math.max(0, state.index - 1)).catch(showError);
  }} else if (key === "arrowright") {{
    event.preventDefault();
    loadCandidate(Math.min(state.count - 1, state.index + 1)).catch(showError);
  }}
}});
loadCandidate(0).catch(err => document.getElementById("status").textContent = err.message);
</script>
</body>
</html>"""


__all__ = [
    "CandidateStore",
    "LABEL_OPTIONS",
    "LightCurveVettingApp",
    "find_leo_report_for_row",
    "find_hlsp_path",
    "tic_shard_path",
]
