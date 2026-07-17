#!/usr/bin/env python3
"""Build a balanced TWIRL-FS light-curve QA gallery from an HLSP tree."""
from __future__ import annotations

import argparse
import csv
import random
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from astropy.io import fits  # noqa: E402


def _mad_sigma(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return float("nan")
    med = np.nanmedian(x)
    return float(1.4826 * np.nanmedian(np.abs(x - med)))


def _ylim_percentile(x: np.ndarray, lo_pct: float, hi_pct: float) -> tuple[float, float] | None:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 3:
        return None
    lo, hi = np.nanpercentile(x, [lo_pct, hi_pct])
    if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
        return None
    pad = 0.08 * (hi - lo)
    return float(lo - pad), float(hi + pad)


def _tmag_bin(tmag: float) -> str:
    if not np.isfinite(tmag):
        return "unknown"
    if tmag < 16.5:
        return "bright"
    if tmag < 18.5:
        return "mid"
    return "dim"


def _header_rows(root: Path) -> list[dict[str, object]]:
    rows = []
    for path in sorted(root.rglob("*.fits")):
        try:
            h = fits.getheader(path, 0)
        except Exception as exc:
            print(f"[warn] cannot read header {path}: {type(exc).__name__}: {exc}")
            continue
        rows.append(
            {
                "path": path,
                "tic": int(h.get("TICID", -1)),
                "tmag": float(h.get("TESSMAG", np.nan)),
                "camera": int(h.get("CAMERA", -1)),
                "ccd": int(h.get("CCD", -1)),
                "method": str(h.get("METHOD", "")),
                "gapsplit": float(h.get("GAPSPLIT", np.nan)),
                "nseg": int(h.get("NSEG", -1)),
                "fitcnt": int(h.get("FITCNT", -1)),
                "scalesrc": str(h.get("SCALESRC", "")),
                "cotstat": str(h.get("COTSTAT", "")),
                "hasadp": bool(h.get("HASADP", False)),
                "adpmeth": str(h.get("ADPMETH", "")),
                "adpbksp": float(h.get("ADPBKSP", np.nan)),
                "adpcots": str(h.get("ADPCOTS", "")),
            }
        )
    return rows


def _hlsp_path(root: Path, sector: int, tic: int) -> Path:
    s = f"{tic:016d}"
    return (
        root
        / s[0:4]
        / s[4:8]
        / s[8:12]
        / s[12:16]
        / f"hlsp_twirlfs_tess_ffi_s{sector:04d}-{s}_tess_v01_llc.fits"
    )


def _candidate_rows_from_orbits(root: Path, sector: int, orbit_roots: list[Path]) -> list[dict[str, object]]:
    by_tic: dict[int, dict[str, object]] = {}
    for orbit_root in orbit_roots:
        for h5 in sorted(orbit_root.rglob("LC/*.h5")):
            try:
                tic = int(h5.stem)
            except ValueError:
                continue
            if tic in by_tic:
                continue
            try:
                ccd = int(h5.parent.parent.name.replace("ccd", ""))
                camera = int(h5.parent.parent.parent.name.replace("cam", ""))
            except ValueError:
                continue
            by_tic[tic] = {
                "path": _hlsp_path(root, sector, tic),
                "tic": tic,
                "tmag": float("nan"),
                "camera": camera,
                "ccd": ccd,
                "method": "",
                "gapsplit": float("nan"),
                "nseg": -1,
                "fitcnt": -1,
                "scalesrc": "",
                "cotstat": "",
                "hasadp": False,
                "adpmeth": "",
                "adpbksp": float("nan"),
                "adpcots": "",
            }
    return list(by_tic.values())


def _read_header_for_candidate(row: dict[str, object]) -> dict[str, object] | None:
    path = Path(row["path"])
    if not path.exists():
        return None
    try:
        h = fits.getheader(path, 0)
    except Exception as exc:
        print(f"[warn] cannot read header {path}: {type(exc).__name__}: {exc}")
        return None
    out = dict(row)
    out.update(
        {
            "tmag": float(h.get("TESSMAG", np.nan)),
            "camera": int(h.get("CAMERA", out["camera"])),
            "ccd": int(h.get("CCD", out["ccd"])),
            "method": str(h.get("METHOD", "")),
            "gapsplit": float(h.get("GAPSPLIT", np.nan)),
            "nseg": int(h.get("NSEG", -1)),
            "fitcnt": int(h.get("FITCNT", -1)),
            "scalesrc": str(h.get("SCALESRC", "")),
            "cotstat": str(h.get("COTSTAT", "")),
            "hasadp": bool(h.get("HASADP", False)),
            "adpmeth": str(h.get("ADPMETH", "")),
            "adpbksp": float(h.get("ADPBKSP", np.nan)),
            "adpcots": str(h.get("ADPCOTS", "")),
        }
    )
    return out


def _sample_header_rows_from_candidates(
    candidates: list[dict[str, object]],
    *,
    seed: int,
    include_tic: int | None,
    per_ccd: int,
) -> list[dict[str, object]]:
    rng = random.Random(seed)
    by_ccd: dict[tuple[int, int], list[dict[str, object]]] = defaultdict(list)
    for row in candidates:
        by_ccd[(int(row["camera"]), int(row["ccd"]))].append(row)
    rows_to_read: list[dict[str, object]] = []
    for ccd_key in sorted(by_ccd):
        bucket = by_ccd[ccd_key]
        rng.shuffle(bucket)
        rows_to_read.extend(bucket[:per_ccd])
    if include_tic is not None:
        rows_to_read.extend([r for r in candidates if int(r["tic"]) == include_tic])

    seen = set()
    out = []
    for row in rows_to_read:
        tic = int(row["tic"])
        if tic in seen:
            continue
        seen.add(tic)
        header_row = _read_header_for_candidate(row)
        if header_row is not None:
            out.append(header_row)
    return out


def _choose_balanced(
    rows: list[dict[str, object]],
    *,
    n: int,
    seed: int,
    include_tic: int | None,
) -> list[dict[str, object]]:
    rng = random.Random(seed)
    selected: list[dict[str, object]] = []
    selected_tics: set[int] = set()

    by_ccd: dict[tuple[int, int], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_ccd[(int(row["camera"]), int(row["ccd"]))].append(row)
    for bucket in by_ccd.values():
        rng.shuffle(bucket)

    if include_tic is not None:
        for row in rows:
            if int(row["tic"]) == include_tic:
                selected.append(row)
                selected_tics.add(include_tic)
                break

    per_ccd = max(1, n // max(1, len(by_ccd)))
    # First pass: make each CCD visible and emphasize dim targets.
    for ccd_key in sorted(by_ccd):
        bucket = [r for r in by_ccd[ccd_key] if int(r["tic"]) not in selected_tics]
        dim = [r for r in bucket if _tmag_bin(float(r["tmag"])) == "dim"]
        mid = [r for r in bucket if _tmag_bin(float(r["tmag"])) == "mid"]
        bright = [r for r in bucket if _tmag_bin(float(r["tmag"])) == "bright"]
        ordered = dim[: max(1, per_ccd - 2)] + mid[:1] + bright[:1]
        if len(ordered) < per_ccd:
            ordered += [r for r in bucket if r not in ordered][: per_ccd - len(ordered)]
        for row in ordered[:per_ccd]:
            if len(selected) >= n:
                break
            selected.append(row)
            selected_tics.add(int(row["tic"]))

    # Fill remaining slots with the dimmest random targets not already selected.
    remaining = [r for r in rows if int(r["tic"]) not in selected_tics]
    remaining.sort(key=lambda r: (0 if _tmag_bin(float(r["tmag"])) == "dim" else 1, -float(r["tmag"])))
    rng.shuffle(remaining[: min(len(remaining), 500)])
    for row in remaining:
        if len(selected) >= n:
            break
        selected.append(row)
        selected_tics.add(int(row["tic"]))

    return selected[:n]


def _read_metrics(
    row: dict[str, object],
    *,
    det_column: str,
    alt_det_column: str | None,
) -> dict[str, object]:
    path = Path(row["path"])
    with fits.open(path, memmap=False) as hdul:
        data = hdul[1].data
        time = np.asarray(data["TIME"], dtype=float)
        sap = np.asarray(data["SAP_FLUX"], dtype=float)
        det = np.asarray(data[det_column], dtype=float)
        quality = np.asarray(data["QUALITY"], dtype=int)
        q0 = quality == 0
        neg_q0 = q0 & np.isfinite(sap) & (sap < 0)
        row = dict(row)
        row.update(
            {
                "n_cadences": int(len(time)),
                "n_q0": int(q0.sum()),
                "neg_sap_q0": int(neg_q0.sum()),
                "finite_det_q0": int(np.isfinite(det[q0]).sum()),
                "det_median_q0": float(np.nanmedian(det[q0])),
                "det_mad_q0": _mad_sigma(det[q0]),
            }
        )
        if alt_det_column and alt_det_column in data.columns.names:
            alt = np.asarray(data[alt_det_column], dtype=float)
            row.update(
                {
                    "finite_alt_q0": int(np.isfinite(alt[q0]).sum()),
                    "alt_median_q0": float(np.nanmedian(alt[q0])),
                    "alt_mad_q0": _mad_sigma(alt[q0]),
                }
            )
    return row


def _plot_cell(
    axes: list[plt.Axes],
    row: dict[str, object],
    *,
    det_column: str,
    alt_det_column: str | None,
) -> None:
    path = Path(row["path"])
    with fits.open(path, memmap=False) as hdul:
        data = hdul[1].data
        time = np.asarray(data["TIME"], dtype=float)
        sap = np.asarray(data["SAP_FLUX"], dtype=float)
        det = np.asarray(data[det_column], dtype=float)
        alt = None
        if alt_det_column and alt_det_column in data.columns.names:
            alt = np.asarray(data[alt_det_column], dtype=float)
        quality = np.asarray(data["QUALITY"], dtype=int)
    q0 = quality == 0
    ax_raw = axes[0]
    ax_det = axes[1]

    if (~q0).any():
        ax_raw.scatter(time[~q0], sap[~q0], s=0.8, c="#b0b0b0", alpha=0.3, linewidths=0)
        ax_det.scatter(time[~q0], det[~q0], s=0.8, c="#b0b0b0", alpha=0.3, linewidths=0)
        if alt is not None and len(axes) > 2:
            axes[2].scatter(time[~q0], alt[~q0], s=0.8, c="#b0b0b0", alpha=0.3, linewidths=0)
    ax_raw.scatter(time[q0], sap[q0], s=0.9, c="#262626", alpha=0.6, linewidths=0)
    ax_det.scatter(time[q0], det[q0], s=0.9, c="#2459a6", alpha=0.6, linewidths=0)
    if alt is not None and len(axes) > 2:
        axes[2].scatter(time[q0], alt[q0], s=0.9, c="#2f7d32", alpha=0.6, linewidths=0)
    ax_raw.axhline(0.0, color="#b22222", lw=0.5, alpha=0.8)
    ax_det.axhline(1.0, color="#b22222", lw=0.5, alpha=0.8)
    if alt is not None and len(axes) > 2:
        axes[2].axhline(1.0, color="#b22222", lw=0.5, alpha=0.8)

    raw_ylim = _ylim_percentile(sap[q0], 1, 99)
    det_ylim = _ylim_percentile(det[q0], 1, 99)
    if raw_ylim is not None:
        ax_raw.set_ylim(*raw_ylim)
    if det_ylim is not None:
        ax_det.set_ylim(*det_ylim)
    if alt is not None and len(axes) > 2:
        alt_ylim = _ylim_percentile(alt[q0], 1, 99)
        if alt_ylim is not None:
            axes[2].set_ylim(*alt_ylim)

    ax_raw.set_title(
        f"TIC {int(row['tic'])} T={float(row['tmag']):.2f} "
        f"C{int(row['camera'])}-{int(row['ccd'])} "
        f"neg={int(row['neg_sap_q0'])} {row['cotstat']}",
        fontsize=5.2,
        pad=1.5,
    )
    for ax in axes:
        ax.set_xticks([])
        ax.tick_params(axis="y", labelsize=4, length=1)
        ax.grid(True, lw=0.25, color="#d0d0d0", alpha=0.7)
        for spine in ax.spines.values():
            spine.set_linewidth(0.3)


def _write_gallery(
    rows: list[dict[str, object]],
    out_dir: Path,
    per_page: int,
    *,
    det_column: str,
    alt_det_column: str | None,
    title: str,
    file_prefix: str,
) -> list[Path]:
    pages = []
    ncols = 5
    nrows = int(np.ceil(per_page / ncols))
    n_panel_rows = 3 if alt_det_column else 2
    height_ratios = [1] * n_panel_rows
    for page_idx, start in enumerate(range(0, len(rows), per_page), 1):
        chunk = rows[start : start + per_page]
        fig = plt.figure(figsize=(17, 12), constrained_layout=False)
        outer = fig.add_gridspec(nrows, ncols, hspace=0.42, wspace=0.18)
        for i, row in enumerate(chunk):
            r, c = divmod(i, ncols)
            inner = outer[r, c].subgridspec(
                n_panel_rows, 1, hspace=0.02, height_ratios=height_ratios
            )
            axes = [fig.add_subplot(inner[j, 0]) for j in range(n_panel_rows)]
            _plot_cell(
                axes,
                row,
                det_column=det_column,
                alt_det_column=alt_det_column,
            )
            if c == 0:
                axes[0].set_ylabel("SAP", fontsize=5)
                axes[1].set_ylabel("DET", fontsize=5)
                if alt_det_column:
                    axes[2].set_ylabel("ADP", fontsize=5)
        fig.suptitle(
            f"{title}, page {page_idx}",
            fontsize=13,
        )
        out = out_dir / f"{file_prefix}_page{page_idx:02d}.png"
        fig.savefig(out, dpi=180, bbox_inches="tight")
        plt.close(fig)
        pages.append(out)
    return pages


def _write_outputs(
    rows: list[dict[str, object]],
    pages: list[Path],
    out_dir: Path,
    seed: int,
    *,
    det_column: str,
    alt_det_column: str | None,
    summary_title: str,
    csv_name: str,
) -> None:
    csv_path = out_dir / csv_name
    fieldnames = [
        "tic", "tmag", "camera", "ccd", "method", "gapsplit", "nseg",
        "fitcnt", "scalesrc", "cotstat", "hasadp", "adpmeth", "adpbksp",
        "adpcots", "n_cadences", "n_q0", "neg_sap_q0", "finite_det_q0",
        "det_median_q0", "det_mad_q0", "finite_alt_q0", "alt_median_q0",
        "alt_mad_q0", "path",
    ]
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})

    ccd_counts = Counter((int(r["camera"]), int(r["ccd"])) for r in rows)
    mag_counts = Counter(_tmag_bin(float(r["tmag"])) for r in rows)
    method_counts = Counter(str(r["method"]) for r in rows)
    cot_counts = Counter(str(r["cotstat"]) for r in rows)
    adp_counts = Counter(str(r.get("adpmeth", "")) for r in rows)
    adp_cot_counts = Counter(str(r.get("adpcots", "")) for r in rows)

    md_path = out_dir / "summary.md"
    with md_path.open("w") as f:
        f.write(f"# {summary_title}\n\n")
        f.write(f"Seed: `{seed}`\n\n")
        f.write(f"Sample size: `{len(rows)}` FITS light curves.\n\n")
        if alt_det_column:
            f.write(
                "Each gallery cell shows raw `SAP_FLUX` on top, "
                f"`{det_column}` in the middle, and `{alt_det_column}` on bottom.\n\n"
            )
        else:
            f.write(
                f"Each gallery cell shows raw `SAP_FLUX` on top and `{det_column}` on bottom.\n\n"
            )
        f.write("## Contact Sheets\n\n")
        for page in pages:
            f.write(f"- [{page.name}]({page.name})\n")
        f.write("\n## Sample Balance\n\n")
        f.write(f"- Method counts: `{dict(method_counts)}`\n")
        if alt_det_column:
            f.write(f"- Adaptive method counts: `{dict(adp_counts)}`\n")
            f.write(f"- Adaptive cotrend status counts: `{dict(adp_cot_counts)}`\n")
        f.write(f"- Tmag-bin counts: `{dict(mag_counts)}`\n")
        f.write(f"- Cotrend status counts: `{dict(cot_counts)}`\n")
        f.write("- Camera/CCD counts:\n")
        for key in sorted(ccd_counts):
            f.write(f"  - cam{key[0]}/ccd{key[1]}: `{ccd_counts[key]}`\n")
        f.write("\n## CSV\n\n")
        f.write(f"- [{csv_path.name}]({csv_path.name})\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--sector", type=int, default=56)
    parser.add_argument(
        "--orbit-roots",
        type=Path,
        nargs="*",
        default=None,
        help=(
            "Optional orbit-*/ffi roots. When provided, sample candidates from "
            "the HDF5 LC layout and read only a bounded candidate set of FITS "
            "headers instead of indexing the full HLSP tree."
        ),
    )
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--n", type=int, default=100)
    parser.add_argument("--seed", type=int, default=5602)
    parser.add_argument("--include-tic", type=int, default=267574918)
    parser.add_argument("--per-page", type=int, default=25)
    parser.add_argument("--candidate-per-ccd", type=int, default=50)
    parser.add_argument("--det-column", default="DET_FLUX")
    parser.add_argument("--alt-det-column", default=None)
    parser.add_argument(
        "--title",
        default="S56 TWIRL-FS v2 quick QA: raw SAP vs DET_FLUX",
    )
    parser.add_argument(
        "--summary-title",
        default="S56 TWIRL-FS v2 Quick QA Gallery",
    )
    parser.add_argument(
        "--file-prefix",
        default="s56_twirl_fs_v2_qa_gallery",
    )
    parser.add_argument(
        "--csv-name",
        default="s56_twirl_fs_v2_qa_sample.csv",
    )
    args = parser.parse_args(argv)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    if args.orbit_roots:
        print(f"[gallery] sampling candidates from {len(args.orbit_roots)} orbit roots")
        candidates = _candidate_rows_from_orbits(args.root, args.sector, args.orbit_roots)
        print(f"[gallery] HDF5-derived candidates: {len(candidates)}")
        headers = _sample_header_rows_from_candidates(
            candidates,
            seed=args.seed,
            include_tic=args.include_tic,
            per_ccd=args.candidate_per_ccd,
        )
    else:
        print(f"[gallery] indexing {args.root}")
        headers = _header_rows(args.root)
    print(f"[gallery] readable FITS headers: {len(headers)}")
    sample = _choose_balanced(headers, n=args.n, seed=args.seed, include_tic=args.include_tic)
    print(f"[gallery] selected {len(sample)} targets")
    rows = [
        _read_metrics(
            row,
            det_column=args.det_column,
            alt_det_column=args.alt_det_column,
        )
        for row in sample
    ]
    pages = _write_gallery(
        rows,
        args.out_dir,
        args.per_page,
        det_column=args.det_column,
        alt_det_column=args.alt_det_column,
        title=args.title,
        file_prefix=args.file_prefix,
    )
    _write_outputs(
        rows,
        pages,
        args.out_dir,
        args.seed,
        det_column=args.det_column,
        alt_det_column=args.alt_det_column,
        summary_title=args.summary_title,
        csv_name=args.csv_name,
    )
    for page in pages:
        print(page)
    print(args.out_dir / "summary.md")
    print(args.out_dir / args.csv_name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
