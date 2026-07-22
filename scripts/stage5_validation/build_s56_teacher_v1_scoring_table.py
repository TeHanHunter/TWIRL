#!/usr/bin/env python3
"""Build training-compatible A2v1 candidate rows for teacher-v1 inference."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.vetting.teacher_candidates import (  # noqa: E402
    enrich_candidate_metadata,
    filter_tier1_eligible_candidates,
    normalize_a2v1_peak_candidates,
    validate_quality_bound_bls_evidence,
    validate_tier1_enrichment_gate,
)
from twirl.vetting.harmonic_inputs import (  # noqa: E402
    CANDIDATE_PROVENANCE_CONTRACT_VERSION,
)


def _read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        return pd.read_parquet(path)
    if suffix in {".csv", ".txt"}:
        return pd.read_csv(path, low_memory=False)
    raise ValueError(f"unsupported table format: {path}")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--adp-peaks", type=Path, required=True)
    parser.add_argument("--adp-peaks-summary", type=Path, required=True)
    parser.add_argument("--compact-lc", type=Path, required=True)
    parser.add_argument("--cadence-reference-table", type=Path, required=True)
    parser.add_argument("--cadence-reference-manifest", type=Path, required=True)
    parser.add_argument("--tier1-target-eligibility", type=Path, required=True)
    parser.add_argument("--tier1-gate-json", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--sector", type=int)
    parser.add_argument("--small-peaks-per-tic", type=int, default=3)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--progress-every", type=int, default=500)
    args = parser.parse_args()
    adp_peaks_sha256 = _sha256(args.adp_peaks)
    peaks = _read_table(args.adp_peaks)
    adp_peaks_summary_bytes = args.adp_peaks_summary.read_bytes()
    adp_peaks_summary_sha256 = hashlib.sha256(adp_peaks_summary_bytes).hexdigest()
    adp_peaks_summary = json.loads(adp_peaks_summary_bytes)
    tier1_target_eligibility_sha256 = _sha256(args.tier1_target_eligibility)
    target_eligibility = _read_table(args.tier1_target_eligibility)
    if _sha256(args.tier1_target_eligibility) != tier1_target_eligibility_sha256:
        raise RuntimeError("Tier-1 target eligibility changed while it was being read")
    tier1_gate_bytes = args.tier1_gate_json.read_bytes()
    tier1_gate_json_sha256 = hashlib.sha256(tier1_gate_bytes).hexdigest()
    tier1_gate_summary = json.loads(tier1_gate_bytes)
    compact_lc_sha256 = _sha256(args.compact_lc)
    cadence_reference_sha256 = _sha256(args.cadence_reference_table)
    cadence_reference_manifest_sha256 = _sha256(
        args.cadence_reference_manifest
    )
    tier1_gate = validate_tier1_enrichment_gate(
        tier1_gate_summary,
        target_eligibility,
        target_eligibility_sha256=tier1_target_eligibility_sha256,
        compact_lc_sha256=compact_lc_sha256,
    )
    bls_evidence = validate_quality_bound_bls_evidence(
        peaks,
        adp_peaks_summary,
        tier1_gate_summary,
        target_eligibility,
        peak_table_sha256=adp_peaks_sha256,
        compact_lc_sha256=compact_lc_sha256,
    )
    if cadence_reference_sha256 != bls_evidence["cadence_reference_sha256"]:
        raise ValueError("candidate metadata cadence table does not match BLS/Tier-1")
    if cadence_reference_manifest_sha256 != bls_evidence[
        "cadence_reference_manifest_sha256"
    ]:
        raise ValueError(
            "candidate metadata cadence manifest does not match BLS/Tier-1"
        )
    candidates = normalize_a2v1_peak_candidates(
        peaks,
        small_peaks_per_tic=args.small_peaks_per_tic,
        sector=args.sector,
    )
    candidates, tier1_filter = filter_tier1_eligible_candidates(
        candidates, target_eligibility
    )
    print(
        "[teacher-candidates] Tier-1 filter: "
        f"{tier1_filter['n_candidates_before']:,} -> "
        f"{tier1_filter['n_candidates_after']:,} candidates; "
        f"{tier1_filter['n_candidate_tics_before']:,} -> "
        f"{tier1_filter['n_candidate_tics_after']:,} TICs",
        flush=True,
    )
    if candidates.empty:
        raise RuntimeError("Tier-1 target eligibility excluded every candidate")
    candidates, summary = enrich_candidate_metadata(
        candidates,
        compact_lc_path=args.compact_lc,
        cadence_reference_table=args.cadence_reference_table,
        cadence_reference_manifest=args.cadence_reference_manifest,
        workers=args.workers,
        progress_every=args.progress_every,
    )
    if not summary["passed"]:
        raise RuntimeError(f"candidate metadata failed: {summary['metadata_status_counts']}")
    if _sha256(args.tier1_target_eligibility) != tier1_target_eligibility_sha256:
        raise RuntimeError("Tier-1 target eligibility changed during candidate build")
    if _sha256(args.adp_peaks) != adp_peaks_sha256:
        raise RuntimeError("A2v1 BLS peak table changed during candidate build")
    if _sha256(args.adp_peaks_summary) != adp_peaks_summary_sha256:
        raise RuntimeError("A2v1 BLS summary changed during candidate build")
    if _sha256(args.tier1_gate_json) != tier1_gate_json_sha256:
        raise RuntimeError("Tier-1 gate JSON changed during candidate build")
    if _sha256(args.compact_lc) != compact_lc_sha256:
        raise RuntimeError("compact LC changed during candidate build")
    if _sha256(args.cadence_reference_table) != cadence_reference_sha256:
        raise RuntimeError("cadence-reference table changed during candidate build")
    if _sha256(
        args.cadence_reference_manifest
    ) != cadence_reference_manifest_sha256:
        raise RuntimeError("cadence-reference manifest changed during candidate build")
    summary.update(
        {
            "provenance_contract_version": (
                CANDIDATE_PROVENANCE_CONTRACT_VERSION
            ),
            "tier1_target_eligibility_path": str(
                args.tier1_target_eligibility.resolve()
            ),
            "tier1_target_eligibility_sha256": tier1_target_eligibility_sha256,
            "tier1_gate_json_path": str(args.tier1_gate_json.resolve()),
            "tier1_gate_json_sha256": tier1_gate_json_sha256,
            "compact_lc_sha256": compact_lc_sha256,
            "cadence_reference_table_path": str(
                args.cadence_reference_table.resolve()
            ),
            "cadence_reference_sha256": cadence_reference_sha256,
            "cadence_reference_manifest_path": str(
                args.cadence_reference_manifest.resolve()
            ),
            "cadence_reference_manifest_sha256": (
                cadence_reference_manifest_sha256
            ),
            "adp_peaks_path": str(args.adp_peaks.resolve()),
            "adp_peaks_sha256": adp_peaks_sha256,
            "adp_peaks_summary_path": str(args.adp_peaks_summary.resolve()),
            "adp_peaks_summary_sha256": adp_peaks_summary_sha256,
            "bls_search_contract_version": bls_evidence[
                "bls_search_contract_version"
            ],
            "bls_config_sha256": bls_evidence["bls_config_sha256"],
            "bls_evidence": bls_evidence,
            "tier1_gate": tier1_gate,
            "tier1_filter": tier1_filter,
        }
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    output = args.out_dir / "teacher_v1_scoring_candidates.parquet"
    try:
        candidates.to_parquet(output, compression="zstd", index=False)
    except (ImportError, ModuleNotFoundError, ValueError):
        output = output.with_suffix(".csv")
        candidates.to_csv(output, index=False)
    candidate_table_sha256 = _sha256(output)
    summary["n_candidate_rows"] = int(len(candidates))
    summary["candidate_table_sha256"] = candidate_table_sha256
    summary["outputs"] = {
        "candidate_table": str(output),
        "summary": str(args.out_dir / "summary.json"),
    }
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
