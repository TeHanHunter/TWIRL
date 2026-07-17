#!/usr/bin/env python3
"""Standalone browser vetting app for a Franklin TWIRL handoff package."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, urlparse


LABEL_OPTIONS: tuple[str, ...] = (
    "planet_like",
    "wide_transit_like",
    "eclipsing_binary_or_pceb",
    "stellar_variability",
    "instrumental_or_systematic",
    "uncertain",
    "skip",
)

LABEL_BUTTONS: tuple[tuple[str, str, str], ...] = (
    ("1", "planet_like", "Planet-like"),
    ("6", "wide_transit_like", "Broad isolated dip"),
    ("2", "eclipsing_binary_or_pceb", "Eclipse/contact"),
    ("3", "stellar_variability", "Smooth variable"),
    ("4", "instrumental_or_systematic", "Systematic/artifact"),
    ("5", "uncertain", "Flat/no signal"),
    ("0", "skip", "Skip"),
)

PERIOD_FACTOR_OPTIONS: tuple[tuple[str, str], ...] = (
    ("0.25", "P/4"),
    ("0.5", "P/2"),
    ("1", "P"),
    ("2", "2P"),
    ("4", "4P"),
    ("unresolved", "Unresolved"),
)

DISPLAY_COLUMNS: tuple[str, ...] = (
    "tic",
    "sector",
    "cam",
    "ccd",
    "tmag",
    "period_d",
    "t0_bjd",
    "duration_min",
    "depth",
    "depth_snr",
    "sde_max",
    "rep_aperture",
    "n_apertures_agree",
    "apertures_agree",
    "centroid_status",
    "centroid_delta_pix",
    "centroid_z",
    "n_in_transit",
)

LABEL_HEADER: tuple[str, ...] = (
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


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return [dict(row) for row in csv.DictReader(fh)]


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: tuple[str, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def _candidate_key(row: dict[str, Any]) -> str:
    return "|".join(str(row.get(col, "") or "") for col in ("review_id", "tic", "sector", "period_d", "t0_bjd"))


def _safe_sheet_names(row: dict[str, Any]) -> list[str]:
    names: list[str] = []
    value = str(row.get("twirl_vet_sheet_name", "") or "")
    if value:
        names.append(value)
    review_id = str(row.get("review_id", "") or "")
    if review_id:
        safe = review_id.replace("/", "_").replace(":", "_")
        names.append(f"{safe}_twirl_twoap_current_adp.png")
        names.append(f"{safe}_twirl_twoap_twirl_fs_v2_adp015q.png")
    tic = str(row.get("tic", "") or "")
    if tic:
        names.append(f"*{tic}*_twirl_twoap_*.png")
    return names


class Store:
    def __init__(self, queue: Path, labels_out: Path, sheet_root: Path, labeler: str) -> None:
        self.queue = Path(queue)
        self.labels_out = Path(labels_out)
        self.sheet_root = Path(sheet_root)
        self.labeler = labeler
        self.rows = _read_csv(self.queue)
        for idx, row in enumerate(self.rows):
            row["row_id"] = str(idx)
            row["candidate_key"] = _candidate_key(row)
        self.labels = self._load_labels()
        self._apply_labels()

    def _load_labels(self) -> list[dict[str, str]]:
        if self.labels_out.exists():
            rows = _read_csv(self.labels_out)
            row_ids = [str(row.get("row_id", "")) for row in rows]
            if len(row_ids) != len(set(row_ids)):
                raise ValueError("label CSV contains duplicate row_id values")
            queue_by_id = {str(row["row_id"]): row for row in self.rows}
            for label in rows:
                row_id = str(label.get("row_id", ""))
                if row_id not in queue_by_id:
                    raise ValueError(f"label CSV row_id {row_id!r} is absent from the queue")
                expected = str(queue_by_id[row_id]["candidate_key"])
                observed = str(label.get("candidate_key", "") or "")
                if observed and observed != expected:
                    raise ValueError(
                        f"candidate_key mismatch for row_id={row_id}: "
                        f"labels={observed!r}, queue={expected!r}"
                    )
            return rows
        return []

    def _apply_labels(self) -> None:
        by_id = {str(row.get("row_id", "")): row for row in self.labels}
        for row in self.rows:
            label = by_id.get(str(row["row_id"]))
            row["label"] = str(label.get("label", "") if label else "")
            row["labeler"] = str(label.get("labeler", self.labeler) if label else self.labeler)
            row["notes"] = str(label.get("notes", "") if label else "")
            row["period_factor"] = str(
                label.get("period_factor", "1") if label else "1"
            ) or "1"
            row["period_status"] = str(
                label.get("period_status", "") if label else ""
            )
            row["updated_utc"] = str(label.get("updated_utc", "") if label else "")

    @property
    def count(self) -> int:
        return len(self.rows)

    def row(self, index: int) -> dict[str, str]:
        if not self.rows:
            raise IndexError("empty queue")
        return self.rows[max(0, min(int(index), len(self.rows) - 1))]

    def sheet_path(self, index: int) -> Path | None:
        row = self.row(index)
        for name in _safe_sheet_names(row):
            if "*" in name:
                matches = sorted(self.sheet_root.glob(name))
                if matches:
                    return matches[0]
            else:
                path = self.sheet_root / name
                if path.exists():
                    return path
        return None

    def payload(self, index: int) -> dict[str, Any]:
        row = self.row(index)
        display = {col: row.get(col, "") for col in DISPLAY_COLUMNS if col in row}
        sheet = self.sheet_path(index)
        return {
            "index": max(0, min(int(index), self.count - 1)),
            "count": self.count,
            "row_id": int(row["row_id"]),
            "candidate_key": row["candidate_key"],
            "display": display,
            "label": row.get("label", ""),
            "labeler": row.get("labeler", self.labeler) or self.labeler,
            "notes": row.get("notes", ""),
            "period_factor": row.get("period_factor", "1") or "1",
            "period_status": row.get("period_status", ""),
            "updated_utc": row.get("updated_utc", ""),
            "sheet_name": sheet.name if sheet else "",
            "sheet_missing": sheet is None,
            "label_options": LABEL_OPTIONS,
            "period_factor_options": PERIOD_FACTOR_OPTIONS,
        }

    def save_label(
        self,
        row_id: int,
        label: str,
        labeler: str,
        notes: str,
        period_factor: str = "1",
    ) -> dict[str, Any]:
        if label not in LABEL_OPTIONS and label != "":
            raise ValueError(f"unsupported label: {label}")
        valid_factors = {value for value, _ in PERIOD_FACTOR_OPTIONS}
        period_factor = str(period_factor or "1")
        if period_factor not in valid_factors:
            raise ValueError(f"unsupported period factor: {period_factor}")
        row = self.rows[int(row_id)]
        record = {
            "row_id": int(row_id),
            "candidate_key": row["candidate_key"],
            "tic": row.get("tic", ""),
            "sector": row.get("sector", ""),
            "label": label,
            "label_source": "human" if label else "",
            "labeler": labeler or self.labeler,
            "notes": notes,
            "period_factor": period_factor,
            "period_status": (
                "unresolved" if period_factor == "unresolved" else "resolved"
            ),
            "updated_utc": datetime.now(timezone.utc).isoformat(),
        }
        kept = [item for item in self.labels if str(item.get("row_id", "")) != str(row_id)]
        kept.append(record)
        self.labels = kept
        _write_csv(self.labels_out, self.labels, LABEL_HEADER)
        self._apply_labels()
        return record

    def summary(self) -> dict[str, Any]:
        counts: dict[str, int] = {}
        factor_counts: dict[str, int] = {}
        for row in self.rows:
            label = row.get("label", "")
            if label:
                counts[label] = counts.get(label, 0) + 1
                factor = str(row.get("period_factor", "1") or "1")
                factor_counts[factor] = factor_counts.get(factor, 0) + 1
        return {
            "queue": str(self.queue),
            "labels_out": str(self.labels_out),
            "sheet_root": str(self.sheet_root),
            "count": self.count,
            "labeled": sum(counts.values()),
            "unlabeled": self.count - sum(counts.values()),
            "label_counts": counts,
            "period_factor_counts": factor_counts,
        }


def _index_html() -> str:
    buttons = "".join(
        f'<button class="label-btn" data-key="{key}" data-label="{label}">{key}: {text}</button>'
        for key, label, text in LABEL_BUTTONS
    )
    period_buttons = "".join(
        f'<button class="period-btn" data-factor="{value}">{text}</button>'
        for value, text in PERIOD_FACTOR_OPTIONS
    )
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>TWIRL Franklin Vetting</title>
<style>
body {{ margin: 0; font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; color: #161616; background: #f5f5f2; }}
.layout {{ display: grid; grid-template-columns: minmax(620px, 1fr) 360px; min-height: 100vh; }}
.sheet {{ padding: 12px; display: flex; align-items: center; justify-content: center; background: #ffffff; border-right: 1px solid #bdbdb8; }}
.sheet img {{ max-width: 100%; max-height: calc(100vh - 24px); object-fit: contain; border: 1px solid #3b3b3b; background: white; }}
.side {{ padding: 16px; overflow: auto; }}
.topline {{ display: flex; align-items: baseline; justify-content: space-between; gap: 12px; margin-bottom: 12px; }}
.counter {{ font-variant-numeric: tabular-nums; font-weight: 700; }}
.nav {{ display: flex; gap: 8px; margin-bottom: 12px; }}
button {{ border: 1px solid #222; background: white; padding: 8px 10px; font: inherit; cursor: pointer; }}
button:hover {{ background: #ecece6; }}
.label-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 8px; margin: 12px 0; }}
.label-btn.active {{ background: #183a59; color: white; }}
.period-title {{ margin-top: 14px; font-size: 13px; font-weight: 700; }}
.period-grid {{ display: grid; grid-template-columns: repeat(3, 1fr); gap: 6px; margin: 6px 0 12px; }}
.period-btn {{ min-width: 0; padding: 6px 4px; font-size: 12px; }}
.period-btn.active {{ background: #183a59; color: white; }}
textarea {{ width: 100%; min-height: 90px; box-sizing: border-box; font: inherit; padding: 8px; border: 1px solid #888; background: white; }}
table {{ width: 100%; border-collapse: collapse; margin-top: 12px; font-size: 13px; }}
td {{ border-bottom: 1px solid #ddd; padding: 4px 2px; vertical-align: top; }}
td:first-child {{ color: #555; width: 42%; }}
.status {{ margin-top: 10px; min-height: 22px; font-size: 13px; color: #333; }}
.missing {{ padding: 40px; border: 1px solid #a33; color: #a33; background: #fff7f7; }}
@media (max-width: 980px) {{ .layout {{ grid-template-columns: 1fr; }} .sheet {{ border-right: 0; border-bottom: 1px solid #bdbdb8; }} .sheet img {{ max-height: 68vh; }} }}
</style>
</head>
<body>
<div class="layout">
  <main class="sheet" id="sheet"></main>
  <aside class="side">
    <div class="topline"><div class="counter" id="counter"></div><div id="saved"></div></div>
    <div class="nav"><button id="prev">Previous</button><button id="next">Next</button><input id="jump" type="number" min="1" style="width:86px" aria-label="row"><button id="go">Go</button></div>
    <div class="label-grid">{buttons}</div>
    <div class="period-title">Best displayed period</div>
    <div class="period-grid">{period_buttons}</div>
    <textarea id="notes" placeholder="Notes, e.g. half period, refold at P/2, possible harmonic"></textarea>
    <div class="status" id="status"></div>
    <table id="meta"></table>
  </aside>
</div>
<script>
let index = Number(new URLSearchParams(location.search).get("i") || 0);
let current = null;
let periodFactor = "1";
async function load(i) {{
  const r = await fetch(`/api/candidate?index=${{i}}`);
  current = await r.json();
  index = current.index;
  history.replaceState(null, "", `?i=${{index}}`);
  document.getElementById("counter").textContent = `${{index + 1}} / ${{current.count}}`;
  document.getElementById("jump").value = index + 1;
  document.getElementById("notes").value = current.notes || "";
  document.getElementById("saved").textContent = current.label ? `saved: ${{current.label}}` : "";
  document.querySelectorAll(".label-btn").forEach(b => b.classList.toggle("active", b.dataset.label === current.label));
  periodFactor = current.period_factor || "1";
  document.querySelectorAll(".period-btn").forEach(b => b.classList.toggle("active", b.dataset.factor === periodFactor));
  const sheet = document.getElementById("sheet");
  if (current.sheet_missing) sheet.innerHTML = '<div class="missing">Missing vet sheet</div>';
  else sheet.innerHTML = `<img src="/vet_sheet.png?index=${{index}}&v=${{Date.now()}}" alt="vet sheet">`;
  const meta = document.getElementById("meta");
  meta.innerHTML = Object.entries(current.display).map(([k,v]) => `<tr><td>${{k}}</td><td>${{v ?? ""}}</td></tr>`).join("");
  document.getElementById("status").textContent = "";
}}
async function save(label) {{
  if (!current) return;
  const payload = {{ row_id: current.row_id, label, labeler: current.labeler, notes: document.getElementById("notes").value, period_factor: periodFactor }};
  const r = await fetch("/api/label", {{ method: "POST", headers: {{ "Content-Type": "application/json" }}, body: JSON.stringify(payload) }});
  const out = await r.json();
  if (!out.ok) {{ document.getElementById("status").textContent = out.error || "save failed"; return; }}
  document.getElementById("status").textContent = `Saved ${{label}}`;
  await load(Math.min(index + 1, current.count - 1));
}}
document.getElementById("prev").onclick = () => load(Math.max(index - 1, 0));
document.getElementById("next").onclick = () => load(Math.min(index + 1, current.count - 1));
document.getElementById("go").onclick = () => load(Math.max(0, Number(document.getElementById("jump").value || 1) - 1));
document.querySelectorAll(".label-btn").forEach(b => b.onclick = () => save(b.dataset.label));
document.querySelectorAll(".period-btn").forEach(b => b.onclick = () => {{
  periodFactor = b.dataset.factor;
  document.querySelectorAll(".period-btn").forEach(x => x.classList.toggle("active", x.dataset.factor === periodFactor));
}});
document.addEventListener("keydown", e => {{
  if (e.target.tagName === "TEXTAREA" || e.target.tagName === "INPUT") return;
  if (e.key === "ArrowLeft") load(Math.max(index - 1, 0));
  else if (e.key === "ArrowRight") load(Math.min(index + 1, current.count - 1));
  else {{
    const btn = Array.from(document.querySelectorAll(".label-btn")).find(b => b.dataset.key === e.key);
    if (btn) save(btn.dataset.label);
  }}
}});
load(index);
</script>
</body>
</html>"""


class Handler(BaseHTTPRequestHandler):
    store: Store

    def do_GET(self) -> None:  # noqa: N802
        parsed = urlparse(self.path)
        query = parse_qs(parsed.query)
        try:
            if parsed.path == "/":
                self._send(_index_html().encode(), "text/html; charset=utf-8")
            elif parsed.path == "/api/candidate":
                index = int(query.get("index", ["0"])[0])
                self._send_json(self.store.payload(index))
            elif parsed.path == "/api/summary":
                self._send_json(self.store.summary())
            elif parsed.path == "/vet_sheet.png":
                index = int(query.get("index", ["0"])[0])
                path = self.store.sheet_path(index)
                if path is None:
                    self.send_error(HTTPStatus.NOT_FOUND, "sheet not found")
                else:
                    self._send(path.read_bytes(), "image/png")
            elif parsed.path == "/labels.csv":
                if not self.store.labels_out.exists():
                    _write_csv(self.store.labels_out, [], LABEL_HEADER)
                self._send(self.store.labels_out.read_bytes(), "text/csv")
            else:
                self.send_error(HTTPStatus.NOT_FOUND, "not found")
        except Exception as exc:
            self._send_json({"error": str(exc)}, status=HTTPStatus.INTERNAL_SERVER_ERROR)

    def do_POST(self) -> None:  # noqa: N802
        if urlparse(self.path).path != "/api/label":
            self.send_error(HTTPStatus.NOT_FOUND, "not found")
            return
        try:
            length = int(self.headers.get("Content-Length", "0"))
            payload = json.loads(self.rfile.read(length).decode())
            record = self.store.save_label(
                row_id=int(payload["row_id"]),
                label=str(payload.get("label", "")),
                labeler=str(payload.get("labeler", "")),
                notes=str(payload.get("notes", "")),
                period_factor=str(payload.get("period_factor", "1")),
            )
            self._send_json({"ok": True, "record": record, "summary": self.store.summary()})
        except Exception as exc:
            self._send_json({"ok": False, "error": str(exc)}, status=HTTPStatus.BAD_REQUEST)

    def log_message(self, fmt: str, *args: Any) -> None:
        return

    def _send_json(self, payload: dict[str, Any], status: HTTPStatus = HTTPStatus.OK) -> None:
        self._send(json.dumps(payload, default=str).encode(), "application/json", status=status)

    def _send(self, body: bytes, content_type: str, status: HTTPStatus = HTTPStatus.OK) -> None:
        self.send_response(int(status))
        self.send_header("Content-Type", content_type)
        self.send_header("Cache-Control", "no-store")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queue", type=Path, default=Path("franklin_review_queue_5k_real.csv"))
    parser.add_argument("--labels-out", type=Path, default=Path("franklin_labels_vetted.csv"))
    parser.add_argument("--sheet-root", type=Path, default=Path("vet_sheets"))
    parser.add_argument("--labeler", default="franklin")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=5003)
    parser.add_argument("--check-only", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    store = Store(args.queue, args.labels_out, args.sheet_root, args.labeler)
    missing = [idx for idx in range(store.count) if store.sheet_path(idx) is None]
    print(json.dumps({**store.summary(), "missing_sheets": len(missing)}, indent=2))
    if args.check_only:
        return 0 if not missing else 1
    Handler.store = store
    server = ThreadingHTTPServer((args.host, int(args.port)), Handler)
    print(f"Open http://{args.host}:{args.port}/")
    server.serve_forever()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
