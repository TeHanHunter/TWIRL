#!/usr/bin/env python3
"""Run TWIRL's S56 semi-supervised candidate triage.

The main mode trains a teacher model on human/injection labels, pseudo-labels
high-confidence unlabeled S56 candidates, and trains a student model on the
combined set. Use ``--write-label-template`` first to produce the human-triage
CSV skeleton from an existing ``vetted_per_tic.parquet`` table.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

DEFAULT_CANDIDATES = (
    REPO_ROOT
    / "data_local/stage2/bls_first_pass_v2/sector_0056/vetted_per_tic_centroid.parquet"
)
DEFAULT_CONFIG = REPO_ROOT / "configs/detection/self_training_s56.yaml"
DEFAULT_OUT = REPO_ROOT / "reports/stage5_validation/self_training_s56"
WD1856_TIC = 267574918


def _read_table(path: Path):
    import pandas as pd

    suffix = path.suffix.lower()
    if suffix == ".parquet":
        return pd.read_parquet(path)
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix in {".json", ".jsonl"}:
        return pd.read_json(path, lines=suffix == ".jsonl")
    raise ValueError(f"unsupported table format: {path}")


def _load_config(path: Path | None):
    from twirl.vetting.self_training import config_from_mapping

    if path is None or not path.exists():
        return config_from_mapping(None)
    import yaml

    with open(path) as fh:
        data = yaml.safe_load(fh) or {}
    return config_from_mapping(data)


def _ordered(df, order):
    sort_cols = [col for col, _ in order if col in df.columns]
    ascending = [asc for col, asc in order if col in df.columns]
    if not sort_cols:
        return df
    return df.sort_values(sort_cols, ascending=ascending, na_position="last")


def _template_key(base):
    if "tic" in base.columns:
        return base["tic"].astype("Int64").astype(str)
    return base.index.astype(str)


def _take_bucket(base, used, bucket: str, mask, quota: int, order):
    if quota <= 0:
        return base.iloc[0:0].copy()
    subset = base.loc[mask & ~base["__template_key"].isin(used)].copy()
    subset = _ordered(subset, order).head(quota)
    used.update(subset["__template_key"].tolist())
    subset["source_bucket"] = bucket
    return subset


def _take_random(base, used, bucket: str, quota: int, random_state: int):
    if quota <= 0:
        return base.iloc[0:0].copy()
    subset = base.loc[~base["__template_key"].isin(used)].copy()
    if len(subset) > quota:
        subset = subset.sample(n=quota, random_state=random_state)
    used.update(subset["__template_key"].tolist())
    subset["source_bucket"] = bucket
    return subset


def _balanced_template(candidates, n: int):
    import pandas as pd

    base = candidates.copy()
    base["__template_key"] = _template_key(base)
    if "vet_class" not in base.columns:
        base["vet_class"] = "unclassified"
    if "centroid_pass" in base.columns:
        centroid_pass = base["centroid_pass"].fillna(False).astype(bool)
    else:
        centroid_pass = pd.Series(False, index=base.index)
    centroid_status = (
        base["centroid_status"].fillna("").astype(str)
        if "centroid_status" in base.columns
        else pd.Series("", index=base.index)
    )

    used: set[str] = set()
    buckets = []
    planet = base["vet_class"].astype(str).eq("planet_candidate")
    pceb_suspect = base["vet_class"].astype(str).eq("sub_roche_pceb_suspect")
    pceb_grid = base["vet_class"].astype(str).eq("pceb_grid_ceiling")
    tmag = base["tmag"] if "tmag" in base.columns else pd.Series(float("nan"), index=base.index)

    rank_order = [
        ("class_rank", True),
        ("blind_rank", True),
        ("sde_max", False),
    ]
    sde_order = [
        ("sde_max", False),
        ("class_rank", True),
        ("blind_rank", True),
    ]

    bucket_plan = [
        (
            "wd1856_benchmark",
            1,
            base["tic"].eq(WD1856_TIC) if "tic" in base.columns else False,
            rank_order,
        ),
        (
            "planet_centroid_pass",
            round(n * 0.14),
            planet & centroid_pass,
            sde_order,
        ),
        (
            "planet_centroid_review",
            round(n * 0.10),
            planet & (~centroid_pass | ~centroid_status.eq("on_target")),
            rank_order,
        ),
        ("planet_candidate_top", round(n * 0.22), planet, rank_order),
        (
            "sub_roche_pceb_suspect",
            round(n * 0.18),
            pceb_suspect,
            rank_order,
        ),
        ("pceb_grid_ceiling", round(n * 0.15), pceb_grid, rank_order),
        ("high_sde_all_classes", round(n * 0.10), pd.Series(True, index=base.index), sde_order),
        (
            "faint_planet_candidate",
            round(n * 0.07),
            planet & tmag.ge(18.5),
            rank_order,
        ),
    ]
    for bucket, quota, mask, order in bucket_plan:
        remaining = n - sum(len(part) for part in buckets)
        if remaining <= 0:
            break
        buckets.append(
            _take_bucket(base, used, bucket, mask, min(int(quota), remaining), order)
        )

    remaining = n - sum(len(part) for part in buckets)
    if remaining > 0:
        buckets.append(_take_random(base, used, "random_background", remaining, 56))

    template = pd.concat(buckets, ignore_index=True) if buckets else base.head(0)
    template = template.drop(columns=[c for c in template.columns if c.startswith("__")])
    return template


def _ranked_template(candidates, n: int):
    base = candidates.copy()
    sort_cols = []
    ascending = []
    if "vet_class" in base.columns:
        class_order = {
            "planet_candidate": 0,
            "sub_roche_pceb_suspect": 1,
            "pceb_grid_ceiling": 2,
            "duration_violator": 3,
            "alias_artifact": 4,
        }
        base["__vet_order"] = base["vet_class"].map(class_order).fillna(99)
        sort_cols.append("__vet_order")
        ascending.append(True)
    if "class_rank" in base.columns:
        sort_cols.append("class_rank")
        ascending.append(True)
    if "sde_max" in base.columns:
        sort_cols.append("sde_max")
        ascending.append(False)
    base = base.sort_values(sort_cols, ascending=ascending) if sort_cols else base
    base = base.head(n).copy()
    base["source_bucket"] = "ranked_candidate"
    return base.drop(columns=[c for c in base.columns if c.startswith("__")])


def _write_labeling_guide(out_path: Path, template) -> None:
    counts = template["source_bucket"].value_counts().to_dict()
    lines = [
        "# S56 Human Labeling Guide",
        "",
        "Fill the `label`, `labeler`, and `notes` columns. Keep `label_source=human` for rows you inspect by eye.",
        "",
        "Suggested label vocabulary:",
        "",
        "- `planet_like`: WD 1856-like or otherwise transit-like and worth preserving.",
        "- `eclipsing_binary_or_pceb`: likely stellar/PCEB eclipse or Roche-limit contaminant.",
        "- `stellar_variability`: variable-star or coherent non-transit behavior.",
        "- `instrumental_or_systematic`: scattered light, aperture, cadence, or detrending artifact.",
        "- `centroid_contaminant`: signal likely belongs to a nearby source.",
        "- `uncertain`: keep in the human-review pool; do not use as a strong class.",
        "- `skip`: unusable row or duplicate inspection target.",
        "",
        "The `source_bucket` column records why each row was selected for this template; it is not a training label.",
        "",
        "Bucket counts:",
        "",
    ]
    lines.extend([f"- `{key}`: {val}" for key, val in counts.items()])
    out_path.write_text("\n".join(lines) + "\n")


def _write_template_summary(out_path: Path, template, mode: str) -> None:
    summary = {
        "mode": mode,
        "rows": int(len(template)),
        "source_bucket_counts": template["source_bucket"].value_counts().to_dict(),
    }
    if "vet_class" in template.columns:
        summary["vet_class_counts"] = template["vet_class"].value_counts().to_dict()
    if "tic" in template.columns:
        summary["contains_wd1856"] = bool(template["tic"].eq(WD1856_TIC).any())
    out_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


def _write_label_template(candidates, out_path: Path, n: int, mode: str) -> None:
    if mode == "ranked":
        base = _ranked_template(candidates, n)
    elif mode == "balanced":
        base = _balanced_template(candidates, n)
    else:
        raise ValueError(f"unsupported template mode: {mode}")

    keep = [
        c for c in [
            "source_bucket",
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
            "centroid_status",
            "centroid_pass",
            "centroid_delta_pix",
            "centroid_z",
            "n_in_transit",
            "n_oot_band",
        ]
        if c in base.columns
    ]
    template = base.loc[:, keep].head(n).copy()
    template.insert(len(template.columns), "label", "")
    template.insert(len(template.columns), "label_source", "human")
    template.insert(len(template.columns), "labeler", "")
    template.insert(len(template.columns), "notes", "")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    template.to_csv(out_path, index=False)
    _write_labeling_guide(out_path.parent / "labeling_guide.md", template)
    _write_template_summary(out_path.parent / "label_template_summary.json", template, mode)
    print(f"[self-train] wrote label template: {out_path} ({len(template):,} rows)")
    print(f"[self-train] wrote labeling guide: {out_path.parent / 'labeling_guide.md'}")


def _build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--candidate-table", type=Path, default=DEFAULT_CANDIDATES)
    ap.add_argument("--labels", type=Path, default=None,
                    help="Human/injection label table. CSV or Parquet with TIC and label columns.")
    ap.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--write-label-template", action="store_true",
                    help="Write a blank human-triage label CSV and exit.")
    ap.add_argument("--template-rows", type=int, default=300)
    ap.add_argument("--template-mode", choices=("balanced", "ranked"), default="balanced")
    return ap


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    if not args.candidate_table.exists():
        print(f"[self-train] missing candidate table: {args.candidate_table}", file=sys.stderr)
        return 2

    candidates = _read_table(args.candidate_table)
    cfg = _load_config(args.config)

    if args.write_label_template:
        out = args.out_dir / "human_labels_template.csv"
        _write_label_template(candidates, out, args.template_rows, args.template_mode)
        return 0

    if args.labels is None or not args.labels.exists():
        print("[self-train] --labels is required for training", file=sys.stderr)
        return 2

    labels = _read_table(args.labels)

    from twirl.vetting.self_training import (
        save_model,
        train_teacher_student,
        write_summary,
    )

    result = train_teacher_student(candidates, labels, cfg)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    result.scored.to_parquet(args.out_dir / "scored_candidates.parquet", compression="zstd")
    result.pseudo_labels.to_parquet(args.out_dir / "pseudo_labels.parquet", compression="zstd")
    result.review_queue.to_csv(args.out_dir / "human_review_queue.csv", index=False)
    save_model(result.teacher, args.out_dir / "teacher_model.npz")
    save_model(result.student, args.out_dir / "student_model.npz")
    write_summary(
        {
            **result.summary,
            "candidate_table": str(args.candidate_table),
            "labels": str(args.labels),
            "config": str(args.config),
            "out_dir": str(args.out_dir),
        },
        args.out_dir / "summary.json",
    )

    print("[self-train] complete")
    print(f"  human labels: {result.summary['n_human_labeled']:,}")
    print(f"  pseudo labels: {result.summary['n_pseudo']:,}")
    print(f"  scored candidates: {len(result.scored):,}")
    print(f"  review queue: {len(result.review_queue):,}")
    print(f"  out: {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
