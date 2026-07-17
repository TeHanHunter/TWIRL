#!/usr/bin/env python3
"""Audit how well an S56 injection set covers the exported light-curve hosts."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_EXPORT_MANIFEST = (
    REPO_ROOT
    / "data_local/stage3_injections/s56_twirlfs_v2_lc_export/"
    / "s56_twirlfs_v2_lc_export_pdo.manifest.json"
)
DEFAULT_INJECTION_MANIFEST = (
    REPO_ROOT
    / "data_local/stage3_injections/s56_twirlfs_v2_injection_training/"
    / "pdo_20k_predetrend_batman_periodradius_bright_balanced/injection_manifest.csv"
)
DEFAULT_INJECTION_SUMMARY = DEFAULT_INJECTION_MANIFEST.with_name("summary.json")
DEFAULT_OUT_DIR = (
    REPO_ROOT
    / "reports/stage5_validation/s56_20k_predetrend_physical_bright_bls_map_pdo/"
    / "injection_host_coverage_orcd"
)
DEFAULT_TMAG_BINS = (0.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 30.0)


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text())


def _records_from_export_manifest(manifest: dict[str, Any]) -> pd.DataFrame:
    records = manifest.get("records", [])
    if not isinstance(records, list):
        raise TypeError("export manifest 'records' must be a list")
    df = pd.DataFrame(records)
    if df.empty:
        return pd.DataFrame(columns=["tic", "tessmag"])
    if "tic" not in df:
        raise KeyError("export manifest records are missing tic")
    df = df.copy()
    df["tic"] = pd.to_numeric(df["tic"], errors="coerce").astype("Int64")
    if "tessmag" in df:
        df["tessmag"] = pd.to_numeric(df["tessmag"], errors="coerce")
    else:
        df["tessmag"] = np.nan
    return df.dropna(subset=["tic"])


def _read_injection_manifest(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "tic" not in df:
        raise KeyError(f"injection manifest is missing tic: {path}")
    df = df.copy()
    df["tic"] = pd.to_numeric(df["tic"], errors="coerce").astype("Int64")
    if "tessmag" in df:
        df["tessmag"] = pd.to_numeric(df["tessmag"], errors="coerce")
    else:
        df["tessmag"] = np.nan
    return df.dropna(subset=["tic"])


def _bin_labels(edges: tuple[float, ...]) -> list[str]:
    labels = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        if lo <= 0:
            labels.append(f"<{hi:g}")
        elif hi >= 30:
            labels.append(f">={lo:g}")
        else:
            labels.append(f"[{lo:g},{hi:g})")
    return labels


def _tmag_bin_summary(
    export_records: pd.DataFrame,
    injections: pd.DataFrame,
    *,
    bins: tuple[float, ...],
) -> list[dict[str, Any]]:
    labels = _bin_labels(bins)
    export = export_records[["tic", "tessmag"]].drop_duplicates("tic").copy()
    injected = injections[["tic", "tessmag"]].copy()
    injected_unique = injected.drop_duplicates("tic").copy()
    for frame in (export, injected, injected_unique):
        frame["tmag_bin"] = pd.cut(
            pd.to_numeric(frame["tessmag"], errors="coerce"),
            bins=bins,
            labels=labels,
            right=False,
            include_lowest=True,
        )
    rows: list[dict[str, Any]] = []
    for label in labels:
        export_bin = export.loc[export["tmag_bin"].astype(str) == label]
        injected_rows = injected.loc[injected["tmag_bin"].astype(str) == label]
        injected_tics = injected_unique.loc[injected_unique["tmag_bin"].astype(str) == label]
        n_export = int(len(export_bin))
        n_injected_tics = int(len(injected_tics))
        rows.append(
            {
                "tmag_bin": label,
                "n_exported_targets": n_export,
                "n_injected_unique_tics": n_injected_tics,
                "n_injections": int(len(injected_rows)),
                "unique_host_fraction": float(n_injected_tics / n_export) if n_export else None,
            }
        )
    return rows


def _value_counts(df: pd.DataFrame, col: str) -> dict[str, int]:
    if col not in df:
        return {}
    return {
        str(key): int(value)
        for key, value in df[col].fillna("").astype(str).value_counts().sort_index().items()
        if str(key) != ""
    }


def _grid_summary(injections: pd.DataFrame, injection_summary: dict[str, Any] | None) -> dict[str, Any]:
    out: dict[str, Any] = {
        "grid_cell_counts_present": "grid_cell_id" in injections,
        "n_grid_cells_observed": 0,
    }
    if injection_summary:
        for key in (
            "sampling_mode",
            "grid_period_bins",
            "grid_radius_bins",
            "grid_depth_bins",
            "grid_period_range_d",
            "grid_radius_range_rearth",
            "grid_depth_range",
        ):
            if key in injection_summary:
                out[key] = injection_summary[key]
    if "grid_cell_id" not in injections:
        return out
    counts = injections["grid_cell_id"].fillna("").astype(str)
    counts = counts[counts != ""].value_counts()
    out.update(
        {
            "n_grid_cells_observed": int(len(counts)),
            "grid_cell_min_count": int(counts.min()) if len(counts) else 0,
            "grid_cell_median_count": float(counts.median()) if len(counts) else None,
            "grid_cell_max_count": int(counts.max()) if len(counts) else 0,
        }
    )
    return out


def build_coverage_audit(
    *,
    export_manifest: Path,
    injection_manifest: Path,
    injection_summary: Path | None = None,
    expected_total_targets: int | None = 19072,
    tmag_bins: tuple[float, ...] = DEFAULT_TMAG_BINS,
) -> dict[str, Any]:
    export_payload = _load_json(export_manifest)
    export_records = _records_from_export_manifest(export_payload)
    injections = _read_injection_manifest(injection_manifest)
    injection_summary_payload = _load_json(injection_summary) if injection_summary and injection_summary.exists() else None

    exported_tics = set(export_records["tic"].dropna().astype(int).tolist())
    injected_tics = set(injections["tic"].dropna().astype(int).tolist())
    missing_from_export = sorted(injected_tics - exported_tics)
    exported_missing_from_injections = sorted(exported_tics - injected_tics)
    n_exported = int(export_payload.get("n_exported_targets", len(exported_tics)))
    n_discovered = int(export_payload.get("n_discovered_files", n_exported))
    skipped = export_payload.get("skipped", {})
    if not isinstance(skipped, dict):
        skipped = {}

    n_injections = int(len(injections))
    n_injected_tics = int(len(injected_tics))
    injections_per_host = injections.groupby("tic").size() if n_injections else pd.Series(dtype=int)
    duplicate_injections_per_host = injections_per_host[injections_per_host > 1].sort_values(ascending=False)
    expected_total = int(expected_total_targets) if expected_total_targets else None
    expected_fraction = float(n_injected_tics / expected_total) if expected_total else None

    warnings = []
    if n_injected_tics < n_exported:
        warnings.append("injection hosts are a subset of exported S56 targets")
    if exported_missing_from_injections:
        warnings.append("some exported TICs are absent from the injection hosts")
    if len(duplicate_injections_per_host):
        warnings.append("some TIC hosts have multiple injections")
    if expected_total and n_injected_tics < expected_total:
        warnings.append("injection hosts are a subset of the expected full S56 HLSP target count")
    if missing_from_export:
        warnings.append("some injected TICs are absent from the export manifest")

    return {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "export_manifest": str(export_manifest),
        "injection_manifest": str(injection_manifest),
        "injection_summary": str(injection_summary) if injection_summary else None,
        "expected_total_targets": expected_total,
        "n_discovered_files": n_discovered,
        "n_exported_targets": n_exported,
        "n_export_manifest_records": int(len(export_records)),
        "export_skipped": {str(key): int(value) for key, value in skipped.items()},
        "n_injections": n_injections,
        "n_unique_injected_tics": n_injected_tics,
        "unique_injected_fraction_of_export": float(n_injected_tics / n_exported) if n_exported else None,
        "unique_injected_fraction_of_expected_total": expected_fraction,
        "n_injected_tics_missing_from_export": int(len(missing_from_export)),
        "sample_injected_tics_missing_from_export": missing_from_export[:20],
        "n_exported_tics_missing_from_injections": int(len(exported_missing_from_injections)),
        "sample_exported_tics_missing_from_injections": exported_missing_from_injections[:20],
        "n_duplicate_injected_tics": int(len(duplicate_injections_per_host)),
        "sample_duplicate_injected_tics": {
            str(int(tic)): int(count) for tic, count in duplicate_injections_per_host.head(20).items()
        },
        "injections_per_host_min": int(injections_per_host.min()) if len(injections_per_host) else 0,
        "injections_per_host_median": float(injections_per_host.median()) if len(injections_per_host) else None,
        "injections_per_host_max": int(injections_per_host.max()) if len(injections_per_host) else 0,
        "split_counts": _value_counts(injections, "split"),
        "signal_family_counts": _value_counts(injections, "signal_family"),
        "aperture_counts": _value_counts(injections, "aperture"),
        "apertures_counts": _value_counts(injections, "apertures"),
        "tmag_bins": _tmag_bin_summary(export_records, injections, bins=tmag_bins),
        "grid_summary": _grid_summary(injections, injection_summary_payload),
        "is_literal_full_s56_host_coverage": bool(
            n_injected_tics >= n_exported and not missing_from_export and not exported_missing_from_injections
        ),
        "warnings": warnings,
    }


def markdown_summary(payload: dict[str, Any]) -> str:
    def fmt_frac(value: Any) -> str:
        if value is None:
            return "n/a"
        return f"{100.0 * float(value):.1f}%"

    lines = [
        "# S56 Injection Host Coverage Audit",
        "",
        f"- injections: `{payload['n_injections']:,}`",
        f"- unique injected TIC hosts: `{payload['n_unique_injected_tics']:,}`",
        f"- exported S56 targets: `{payload['n_exported_targets']:,}`",
        f"- expected full S56 targets: `{payload['expected_total_targets']:,}`"
        if payload.get("expected_total_targets")
        else "- expected full S56 targets: `n/a`",
        f"- unique host coverage of export: `{fmt_frac(payload['unique_injected_fraction_of_export'])}`",
        f"- literal full-S56 host coverage: `{payload['is_literal_full_s56_host_coverage']}`",
        f"- exported TICs absent from injections: `{payload.get('n_exported_tics_missing_from_injections', 0):,}`",
        f"- injected TICs absent from export: `{payload.get('n_injected_tics_missing_from_export', 0):,}`",
        f"- duplicate injected TIC hosts: `{payload.get('n_duplicate_injected_tics', 0):,}`",
        "",
    ]
    if payload["warnings"]:
        lines.append("## Warnings")
        lines.extend(f"- {item}" for item in payload["warnings"])
        lines.append("")
    if payload.get("sample_exported_tics_missing_from_injections") or payload.get("sample_duplicate_injected_tics"):
        lines.append("## Host Gaps")
        missing = payload.get("sample_exported_tics_missing_from_injections", [])
        duplicates = payload.get("sample_duplicate_injected_tics", {})
        if missing:
            lines.append(
                "- sample exported TICs absent from injections: "
                + ", ".join(f"`{tic}`" for tic in missing)
            )
        if duplicates:
            lines.append(
                "- sample duplicate injected TIC hosts: "
                + ", ".join(f"`{tic}` ({count})" for tic, count in duplicates.items())
            )
        lines.append("")

    lines.extend(
        [
            "## Tmag Coverage",
            "",
            "| Tmag bin | exported targets | injected TIC hosts | injections | host coverage |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for row in payload["tmag_bins"]:
        lines.append(
            "| {tmag_bin} | {n_exported_targets:,} | {n_injected_unique_tics:,} | "
            "{n_injections:,} | {fraction} |".format(
                **row,
                fraction=fmt_frac(row["unique_host_fraction"]),
            )
        )
    grid = payload.get("grid_summary", {})
    lines.extend(
        [
            "",
            "## Parameter Grid",
            "",
            f"- sampling mode: `{grid.get('sampling_mode', '')}`",
            f"- observed grid cells: `{grid.get('n_grid_cells_observed', 0)}`",
            f"- per-cell count min/median/max: "
            f"`{grid.get('grid_cell_min_count', 0)}` / "
            f"`{grid.get('grid_cell_median_count', 'n/a')}` / "
            f"`{grid.get('grid_cell_max_count', 0)}`",
            "",
            "## Interpretation",
            "",
            "This audit describes host coverage only. Injection truth can train signal "
            "visibility and peak-ranking behavior, but it should not be treated as real "
            "object taxonomy. Real false-positive taxonomy still comes from LEO-assisted "
            "human labels.",
            "",
        ]
    )
    return "\n".join(lines)


def _parse_float_tuple(raw: str) -> tuple[float, ...]:
    return tuple(float(part.strip()) for part in str(raw).split(",") if part.strip())


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--export-manifest", type=Path, default=DEFAULT_EXPORT_MANIFEST)
    parser.add_argument("--injection-manifest", type=Path, default=DEFAULT_INJECTION_MANIFEST)
    parser.add_argument("--injection-summary", type=Path, default=DEFAULT_INJECTION_SUMMARY)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--expected-total-targets", type=int, default=19072)
    parser.add_argument("--tmag-bins", default=",".join(str(value) for value in DEFAULT_TMAG_BINS))
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    payload = build_coverage_audit(
        export_manifest=args.export_manifest,
        injection_manifest=args.injection_manifest,
        injection_summary=args.injection_summary,
        expected_total_targets=args.expected_total_targets,
        tmag_bins=_parse_float_tuple(args.tmag_bins),
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    summary_json = args.out_dir / "summary.json"
    summary_md = args.out_dir / "summary.md"
    summary_json.write_text(json.dumps(payload, indent=2, sort_keys=True, default=_json_default) + "\n")
    summary_md.write_text(markdown_summary(payload))
    print(json.dumps(payload, indent=2, sort_keys=True, default=_json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
