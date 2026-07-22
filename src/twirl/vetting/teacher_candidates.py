"""Build A2v1 real-candidate rows for harmonic-teacher inference."""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd

from twirl.io.compact_export import read_compact_lc_export
from twirl.lightcurves.a2v1_cadence_reference import (
    S56_EXPECTED_DETECTORS,
    S56_EXPECTED_ORBITS,
)
from twirl.lightcurves.external_quality import (
    ExternalQualityReference,
    load_external_quality_reference,
)
from twirl.search.a2v1_bls_contract import (
    A2V1_TEACHER_BLS_SEARCH_CONTRACT,
    approved_a2v1_teacher_bls_config,
    bls_config_sha256,
)
from twirl.vetting.adp_only import (
    ADP_ONLY_APERTURES,
    ADP_ONLY_CONTRACT_VERSION,
    assert_adp_only_search_frame,
    classify_period_relation,
)
from twirl.vetting.harmonic_inputs import native_group_path
from twirl.vetting.label_io import candidate_key
from twirl.vetting.teacher_active_learning import A2V1_TEACHER_INPUT_CONTRACT
from twirl.vetting.two_aperture import measure_two_aperture_candidate_metadata


PEAK_FIELDS: tuple[str, ...] = (
    "peak_rank",
    "period_d",
    "t0_bjd",
    "duration_min",
    "depth",
    "depth_snr",
    "sde",
    "log_power",
)

TIER1_ELIGIBILITY_KEY: tuple[str, str] = ("sector", "tic")
TIER1_ELIGIBILITY_PROVENANCE_FIELDS: tuple[str, ...] = (
    "tier1_contract_version",
    "tier1_config_name",
    "tier1_scope",
    "sector_tic_key",
    "tier1_target_qa_status",
    "tier1_target_qa_reasons",
    "tier1_target_qa_pass",
)

A2V1_BLS_EXTERNAL_QUALITY_CONTRACT = (
    "a2v1_bls_internal_or_authoritative_external_quality_v1"
)
A2V1_CADENCE_REFERENCE_CONTRACT_TEMPLATE = "s{sector}_a2v1_cadence_reference_v1"
A2V1_CADENCE_AUTHORITY = "qlp_cam_quat"
A2V1_QUALITY_AUTHORITY = "spoc_and_qlp_quality_flags"
TIER1_REQUIRED_ENRICHMENT_GATES: tuple[str, ...] = (
    "cadence_reference_prerequisite",
    "injection_source_parity_prerequisite",
    "tier0_prerequisite",
    "population_scatter",
    "cadence_and_finite_data",
    "aperture_outliers",
    "fixed_injection_preservation",
    "independent_extraction",
)
_METADATA_QUALITY_REFERENCE: ExternalQualityReference | None = None


def _validated_sha256(value: Any, *, label: str) -> str:
    if not isinstance(value, str):
        raise ValueError(f"{label} must be a SHA-256 string")
    normalized = value.lower()
    if len(normalized) != 64 or any(
        character not in "0123456789abcdef" for character in normalized
    ):
        raise ValueError(f"{label} must be a valid SHA-256 string")
    return normalized


def validate_tier1_enrichment_gate(
    gate_summary: Mapping[str, Any],
    target_eligibility: pd.DataFrame,
    *,
    target_eligibility_sha256: str,
    compact_lc_sha256: str,
) -> dict[str, Any]:
    """Validate that target eligibility came from a passed enrichment gate."""

    if not isinstance(gate_summary, Mapping):
        raise TypeError("Tier-1 gate summary must be a JSON object")
    if gate_summary.get("status") != "pass":
        raise ValueError("Tier-1 gate summary status is not pass")
    if gate_summary.get("passed") is not True:
        raise ValueError("Tier-1 gate summary passed is not true")
    if gate_summary.get("enrichment_ready") is not True:
        raise ValueError("Tier-1 gate summary enrichment_ready is not true")
    if gate_summary.get("qa_tier") != "tier1_bounded_enrichment_qa":
        raise ValueError("Tier-1 gate summary has an incompatible qa_tier")
    if gate_summary.get("scope") != "active_search_pair":
        raise ValueError("Tier-1 gate summary is not for active_search_pair")
    if gate_summary.get("science_ready") is not False:
        raise ValueError("bounded Tier-1 gate science_ready must be explicitly false")
    if gate_summary.get("promotion_enabled") is not False:
        raise ValueError("bounded Tier-1 gate promotion_enabled must be false")
    if gate_summary.get("apertures") != list(ADP_ONLY_APERTURES):
        raise ValueError("Tier-1 gate does not cover the locked ADP pair")
    gates = gate_summary.get("gates")
    if not isinstance(gates, Mapping):
        raise ValueError("Tier-1 gate summary is missing its nested gates")
    if set(gates) != set(TIER1_REQUIRED_ENRICHMENT_GATES):
        raise ValueError("Tier-1 nested gate inventory is incompatible")
    for gate_name in TIER1_REQUIRED_ENRICHMENT_GATES:
        nested = gates.get(gate_name)
        if not isinstance(nested, Mapping) or nested.get("status") != "pass":
            raise ValueError(f"Tier-1 nested gate {gate_name} did not pass")

    summary_sector = gate_summary.get("sector")
    if type(summary_sector) is not int or summary_sector <= 0:
        raise ValueError("Tier-1 gate summary has an invalid sector")
    summary_identity = {
        "tier1_contract_version": gate_summary.get("contract_version"),
        "tier1_config_name": gate_summary.get("config_name"),
        "tier1_scope": gate_summary.get("scope"),
    }
    required_columns = {
        "sector",
        "tier1_target_qa_status",
        "tier1_target_qa_pass",
        *summary_identity,
    }
    missing = sorted(required_columns - set(target_eligibility.columns))
    if missing:
        raise KeyError(
            "Tier-1 target eligibility is missing gate identity columns: "
            f"{missing}"
        )
    eligibility_sectors = _normalize_positive_integer_key(
        target_eligibility["sector"],
        column="sector",
        table_name="Tier-1 target eligibility",
    )
    if set(eligibility_sectors.tolist()) != {summary_sector}:
        raise ValueError(
            "Tier-1 gate summary sector does not match target eligibility"
        )
    for column, summary_value in summary_identity.items():
        if not isinstance(summary_value, str) or not summary_value.strip():
            raise ValueError(f"Tier-1 gate summary has invalid {column}")
        values = target_eligibility[column].astype("string")
        observed = set(values.dropna().tolist())
        if values.isna().any() or observed != {summary_value}:
            raise ValueError(
                f"Tier-1 gate summary {column} does not match target eligibility"
            )

    target_qa = gate_summary.get("target_qa")
    if not isinstance(target_qa, Mapping):
        raise ValueError("Tier-1 gate summary is missing target_qa provenance")
    if target_qa.get("observation_key") != list(TIER1_ELIGIBILITY_KEY):
        raise ValueError("Tier-1 gate summary has an incompatible observation key")
    if target_qa.get("candidate_teacher_filter") != (
        "tier1_target_qa_pass == True"
    ):
        raise ValueError("Tier-1 gate summary has an incompatible candidate filter")
    eligibility_status = target_eligibility["tier1_target_qa_status"].astype("string")
    if eligibility_status.isna().any() or not eligibility_status.isin(
        ("pass", "review", "fail")
    ).all():
        raise ValueError("Tier-1 target eligibility has invalid QA status values")
    eligibility_pass = _normalize_strict_boolean(
        target_eligibility["tier1_target_qa_pass"],
        column="tier1_target_qa_pass",
    )
    if not eligibility_pass.eq(eligibility_status.eq("pass")).all():
        raise ValueError("Tier-1 target eligibility status/pass values disagree")
    expected_counts = {
        "n_pass": int(eligibility_status.eq("pass").sum()),
        "n_review": int(eligibility_status.eq("review").sum()),
        "n_fail": int(eligibility_status.eq("fail").sum()),
    }
    for field, expected in expected_counts.items():
        observed = target_qa.get(field)
        if type(observed) is not int or observed != expected:
            raise ValueError(f"Tier-1 gate target_qa {field} is inconsistent")

    provenance = gate_summary.get("provenance")
    output_sha256 = (
        provenance.get("output_sha256") if isinstance(provenance, Mapping) else None
    )
    if not isinstance(output_sha256, Mapping):
        raise ValueError("Tier-1 gate summary is missing output SHA-256 provenance")
    expected_eligibility_sha256 = _validated_sha256(
        output_sha256.get("target_eligibility"),
        label="Tier-1 gate target-eligibility provenance",
    )
    observed_eligibility_sha256 = _validated_sha256(
        target_eligibility_sha256,
        label="observed Tier-1 target eligibility",
    )
    if expected_eligibility_sha256 != observed_eligibility_sha256:
        raise ValueError(
            "Tier-1 target eligibility SHA-256 does not match the gate summary"
        )
    expected_compact_sha256 = _validated_sha256(
        provenance.get("compact_lc_sha256"),
        label="Tier-1 gate compact-LC provenance",
    )
    observed_compact_sha256 = _validated_sha256(
        compact_lc_sha256,
        label="observed compact LC",
    )
    if expected_compact_sha256 != observed_compact_sha256:
        raise ValueError("compact LC SHA-256 does not match the Tier-1 gate summary")

    return {
        "status": "pass",
        "passed": True,
        "enrichment_ready": True,
        "sector": summary_sector,
        "scope": str(gate_summary["scope"]),
        "contract_version": str(gate_summary["contract_version"]),
        "config_name": str(gate_summary["config_name"]),
        "observation_key": list(TIER1_ELIGIBILITY_KEY),
        "target_eligibility_sha256": observed_eligibility_sha256,
        "compact_lc_sha256": observed_compact_sha256,
    }


def validate_quality_bound_bls_evidence(
    peaks: pd.DataFrame,
    summary: Mapping[str, Any],
    gate_summary: Mapping[str, Any],
    target_eligibility: pd.DataFrame,
    *,
    peak_table_sha256: str,
    compact_lc_sha256: str,
) -> dict[str, Any]:
    """Bind a complete ADP peak table to the passed Tier-1 evidence set.

    This check prevents a valid Tier-1 target mask from being joined to peaks
    produced with the older internal-only quality mask, a different compact
    export, a partial BLS shard, or incomplete eligible-target coverage.
    """

    if not isinstance(summary, Mapping):
        raise TypeError("A2v1 BLS summary must be a JSON object")
    if not isinstance(gate_summary, Mapping):
        raise TypeError("Tier-1 gate summary must be a JSON object")
    if peaks.empty:
        raise ValueError("A2v1 BLS peak table is empty")
    observed_peak_sha = _validated_sha256(
        peak_table_sha256, label="observed A2v1 BLS peak table"
    )
    observed_compact_sha = _validated_sha256(
        compact_lc_sha256, label="observed compact LC"
    )
    if _validated_sha256(
        summary.get("peak_table_sha256"), label="A2v1 BLS peak-table provenance"
    ) != observed_peak_sha:
        raise ValueError("A2v1 BLS peak-table SHA-256 does not match its summary")
    if _validated_sha256(
        summary.get("compact_lc_sha256"), label="A2v1 BLS compact-LC provenance"
    ) != observed_compact_sha:
        raise ValueError("A2v1 BLS compact-LC SHA-256 does not match the input")

    sector = summary.get("sector")
    if type(sector) is not int or sector <= 0:
        raise ValueError("A2v1 BLS summary has an invalid sector")
    if sector != gate_summary.get("sector"):
        raise ValueError("A2v1 BLS and Tier-1 sectors do not match")
    if summary.get("contract_version") != ADP_ONLY_CONTRACT_VERSION:
        raise ValueError("A2v1 BLS summary has an incompatible ADP contract")
    if summary.get("external_quality_policy_contract") != (
        A2V1_BLS_EXTERNAL_QUALITY_CONTRACT
    ):
        raise ValueError("A2v1 BLS summary lacks the authoritative quality policy")
    if summary.get("apertures") != list(ADP_ONLY_APERTURES):
        raise ValueError("A2v1 BLS summary does not cover the locked ADP pair")
    if summary.get("source_product_tag") != "A2v1":
        raise ValueError("A2v1 BLS summary has the wrong source product tag")
    expected_bls_config = approved_a2v1_teacher_bls_config()
    if summary.get("bls_search_contract_version") != (
        A2V1_TEACHER_BLS_SEARCH_CONTRACT
    ):
        raise ValueError("A2v1 BLS summary has an incompatible search contract")
    if summary.get("config") != expected_bls_config:
        raise ValueError("A2v1 BLS summary does not use the approved search config")
    expected_bls_config_sha256 = bls_config_sha256(expected_bls_config)
    if _validated_sha256(
        summary.get("bls_config_sha256"), label="A2v1 BLS config provenance"
    ) != expected_bls_config_sha256:
        raise ValueError("A2v1 BLS config SHA-256 is inconsistent")
    if int(summary.get("n_shards", -1)) != 1 or int(
        summary.get("shard_index", -1)
    ) != 0:
        raise ValueError("teacher scoring requires one complete, unsharded BLS table")
    if int(summary.get("n_targets", -1)) != int(
        summary.get("n_targets_total", -2)
    ):
        raise ValueError("A2v1 BLS summary describes a partial target population")
    if int(summary.get("n_rows", -1)) != len(peaks):
        raise ValueError("A2v1 BLS summary row count disagrees with the peak table")

    provenance = gate_summary.get("provenance")
    if not isinstance(provenance, Mapping):
        raise ValueError("Tier-1 gate lacks provenance for BLS binding")
    cadence_gate = gate_summary.get("gates", {}).get(
        "cadence_reference_prerequisite"
    ) if isinstance(gate_summary.get("gates"), Mapping) else None
    if not isinstance(cadence_gate, Mapping) or cadence_gate.get("status") != "pass":
        raise ValueError("Tier-1 cadence-reference prerequisite did not pass")
    binding_fields = (
        ("compact_lc_sha256", "compact_lc_sha256"),
        ("cadence_reference_sha256", "cadence_reference_sha256"),
        (
            "cadence_reference_manifest_sha256",
            "cadence_reference_manifest_sha256",
        ),
    )
    for summary_field, gate_field in binding_fields:
        observed = _validated_sha256(
            summary.get(summary_field), label=f"A2v1 BLS {summary_field}"
        )
        expected = _validated_sha256(
            provenance.get(gate_field), label=f"Tier-1 gate {gate_field}"
        )
        if observed != expected:
            raise ValueError(
                f"A2v1 BLS {summary_field} does not match the Tier-1 gate"
            )
    if summary.get("cadence_reference_contract_version") != (
        A2V1_CADENCE_REFERENCE_CONTRACT_TEMPLATE.format(sector=sector)
    ):
        raise ValueError("A2v1 BLS cadence-reference contract mismatch")
    if summary.get("cadence_reference_cadence_authority") != A2V1_CADENCE_AUTHORITY:
        raise ValueError("A2v1 BLS cadence authority mismatch")
    if summary.get("cadence_reference_quality_authority") != A2V1_QUALITY_AUTHORITY:
        raise ValueError("A2v1 BLS quality authority mismatch")
    if cadence_gate.get("cadence_authority") != A2V1_CADENCE_AUTHORITY or (
        cadence_gate.get("quality_authority") != A2V1_QUALITY_AUTHORITY
    ):
        raise ValueError("Tier-1 cadence-reference authorities are incompatible")
    for summary_field, gate_field in (
        ("cadence_reference_sha256", "table_sha256"),
        ("cadence_reference_manifest_sha256", "manifest_sha256"),
    ):
        if _validated_sha256(
            summary.get(summary_field), label=f"A2v1 BLS {summary_field}"
        ) != _validated_sha256(
            cadence_gate.get(gate_field), label=f"Tier-1 cadence gate {gate_field}"
        ):
            raise ValueError("A2v1 BLS cadence evidence disagrees with Tier-1")
    _validated_sha256(
        summary.get("cadence_reference_source_hashes_sha256"),
        label="A2v1 BLS cadence source declarations",
    )

    required_peak_columns = {
        "sector",
        "tic",
        "aperture",
        "peak_rank",
        "status",
        "source_product_tag",
        "adp_only_contract_version",
        "external_quality_policy_contract",
        "bls_search_contract_version",
        "bls_config_sha256",
        "cadence_reference_sha256",
        "cadence_reference_manifest_sha256",
    }
    missing = sorted(required_peak_columns - set(peaks.columns))
    if missing:
        raise KeyError(f"A2v1 BLS peak table lacks provenance columns: {missing}")
    if peaks.columns.duplicated().any():
        raise ValueError("A2v1 BLS peak table contains duplicate column names")
    peak_sector = _normalize_positive_integer_key(
        peaks["sector"], column="sector", table_name="A2v1 BLS peak table"
    )
    peak_tic = _normalize_positive_integer_key(
        peaks["tic"], column="tic", table_name="A2v1 BLS peak table"
    )
    if set(peak_sector.tolist()) != {sector}:
        raise ValueError("A2v1 BLS peak rows have the wrong sector")
    exact_row_values = {
        "source_product_tag": "A2v1",
        "adp_only_contract_version": ADP_ONLY_CONTRACT_VERSION,
        "external_quality_policy_contract": A2V1_BLS_EXTERNAL_QUALITY_CONTRACT,
        "bls_search_contract_version": A2V1_TEACHER_BLS_SEARCH_CONTRACT,
        "bls_config_sha256": expected_bls_config_sha256,
        "cadence_reference_sha256": str(summary["cadence_reference_sha256"]),
        "cadence_reference_manifest_sha256": str(
            summary["cadence_reference_manifest_sha256"]
        ),
    }
    for column, expected in exact_row_values.items():
        values = peaks[column].astype("string")
        if values.isna().any() or set(values.tolist()) != {expected}:
            raise ValueError(f"A2v1 BLS rows have inconsistent {column}")

    eligibility_sector = _normalize_positive_integer_key(
        target_eligibility["sector"],
        column="sector",
        table_name="Tier-1 target eligibility",
    )
    eligibility_tic = _normalize_positive_integer_key(
        target_eligibility["tic"],
        column="tic",
        table_name="Tier-1 target eligibility",
    )
    eligibility_keys = set(zip(eligibility_sector, eligibility_tic, strict=True))
    peak_keys = set(zip(peak_sector, peak_tic, strict=True))
    if peak_keys != eligibility_keys:
        missing_keys = sorted(eligibility_keys - peak_keys)[:10]
        extra_keys = sorted(peak_keys - eligibility_keys)[:10]
        raise ValueError(
            "A2v1 BLS target coverage does not exactly match Tier-1 eligibility; "
            f"missing={missing_keys}, extra={extra_keys}"
        )
    if int(summary.get("n_targets", -1)) != len(peak_keys) or int(
        summary.get("n_unique_tics", -1)
    ) != len({tic for _, tic in peak_keys}):
        raise ValueError("A2v1 BLS target counts disagree with the peak table")

    eligible_pass = _normalize_strict_boolean(
        target_eligibility["tier1_target_qa_pass"],
        column="tier1_target_qa_pass",
    )
    passing_keys = set(
        zip(
            eligibility_sector.loc[eligible_pass],
            eligibility_tic.loc[eligible_pass],
            strict=True,
        )
    )
    valid_rank = pd.to_numeric(peaks["peak_rank"], errors="coerce").eq(1)
    valid = peaks["status"].astype(str).eq("ok") & valid_rank
    valid_rows = peaks.loc[valid, ["sector", "tic", "aperture"]].copy()
    valid_rows["sector"] = peak_sector.loc[valid].to_numpy()
    valid_rows["tic"] = peak_tic.loc[valid].to_numpy()
    coverage = {
        (int(key[0]), int(key[1])): set(group["aperture"].astype(str))
        for key, group in valid_rows.groupby(["sector", "tic"], sort=False)
    }
    incomplete = sorted(
        key
        for key in passing_keys
        if coverage.get(key, set()) != set(ADP_ONLY_APERTURES)
    )
    if incomplete:
        raise ValueError(
            "Tier-1-passing targets lack rank-1 BLS results in both apertures; "
            f"examples={incomplete[:10]}"
        )
    return {
        "status": "pass",
        "sector": sector,
        "peak_table_sha256": observed_peak_sha,
        "compact_lc_sha256": observed_compact_sha,
        "cadence_reference_sha256": str(summary["cadence_reference_sha256"]),
        "cadence_reference_manifest_sha256": str(
            summary["cadence_reference_manifest_sha256"]
        ),
        "quality_policy_contract": A2V1_BLS_EXTERNAL_QUALITY_CONTRACT,
        "bls_search_contract_version": A2V1_TEACHER_BLS_SEARCH_CONTRACT,
        "bls_config_sha256": expected_bls_config_sha256,
        "n_targets": len(peak_keys),
        "n_tier1_passing_targets": len(passing_keys),
        "n_rows": len(peaks),
    }


def _normalize_positive_integer_key(
    values: pd.Series, *, column: str, table_name: str
) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    finite = numeric.notna() & np.isfinite(numeric)
    integral = finite & numeric.eq(np.floor(numeric))
    positive = integral & numeric.gt(0)
    if not positive.all():
        bad = values.loc[~positive].astype(str).head(5).tolist()
        raise ValueError(
            f"{table_name} has invalid {column} values; expected positive integers, "
            f"examples={bad}"
        )
    return numeric.astype(np.int64)


def _normalize_strict_boolean(values: pd.Series, *, column: str) -> pd.Series:
    if pd.api.types.is_bool_dtype(values.dtype):
        if values.isna().any():
            raise ValueError(f"Tier-1 eligibility has missing {column} values")
        return values.astype(bool)
    tokens = values.astype("string").str.strip().str.lower()
    invalid = tokens.isna() | ~tokens.isin(("true", "false"))
    if invalid.any():
        bad = values.loc[invalid].astype(str).head(5).tolist()
        raise ValueError(
            f"Tier-1 eligibility {column} must contain only true/false values; "
            f"examples={bad}"
        )
    return tokens.eq("true")


def filter_tier1_eligible_candidates(
    candidates: pd.DataFrame,
    target_eligibility: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Attach the exact Tier-1 target decision and retain only passing rows.

    The observation key is exactly ``(sector, tic)``.  Candidate rows may repeat
    an observation (for example, one row per BLS peak), while the eligibility
    input must contain one and only one decision for every candidate
    observation.  Extra eligibility rows are allowed because the Tier-1 product
    describes the full audited sector population.
    """

    candidate_required = set(TIER1_ELIGIBILITY_KEY)
    eligibility_required = candidate_required | set(TIER1_ELIGIBILITY_PROVENANCE_FIELDS)
    missing_candidates = sorted(candidate_required - set(candidates.columns))
    if missing_candidates:
        raise KeyError(
            f"candidate table is missing Tier-1 join columns: {missing_candidates}"
        )
    missing_eligibility = sorted(eligibility_required - set(target_eligibility.columns))
    if missing_eligibility:
        raise KeyError(
            "Tier-1 target eligibility is missing required columns: "
            f"{missing_eligibility}"
        )
    if candidates.empty:
        raise ValueError("candidate table is empty before Tier-1 filtering")
    if target_eligibility.empty:
        raise ValueError("Tier-1 target eligibility table is empty")
    if candidates.columns.duplicated().any():
        raise ValueError("candidate table contains duplicate column names")
    if target_eligibility.columns.duplicated().any():
        raise ValueError("Tier-1 target eligibility contains duplicate column names")

    stale = sorted(set(TIER1_ELIGIBILITY_PROVENANCE_FIELDS) & set(candidates.columns))
    if stale:
        raise ValueError(
            "candidate table already contains Tier-1 eligibility provenance; "
            f"refusing to overwrite stale columns: {stale}"
        )

    work = candidates.copy()
    eligibility = target_eligibility.copy()
    for column in TIER1_ELIGIBILITY_KEY:
        work[column] = _normalize_positive_integer_key(
            work[column], column=column, table_name="candidate table"
        )
        eligibility[column] = _normalize_positive_integer_key(
            eligibility[column],
            column=column,
            table_name="Tier-1 target eligibility",
        )

    duplicate_mask = eligibility.duplicated(list(TIER1_ELIGIBILITY_KEY), keep=False)
    if duplicate_mask.any():
        duplicate_keys = (
            eligibility.loc[duplicate_mask, list(TIER1_ELIGIBILITY_KEY)]
            .drop_duplicates()
            .head(5)
            .to_dict("records")
        )
        raise ValueError(
            "Tier-1 target eligibility requires unique (sector, TIC) rows; "
            f"duplicate keys={duplicate_keys}"
        )

    status = eligibility["tier1_target_qa_status"].astype("string")
    invalid_status = status.isna() | ~status.isin(("pass", "review", "fail"))
    if invalid_status.any():
        bad = (
            eligibility.loc[invalid_status, "tier1_target_qa_status"]
            .astype(str)
            .head(5)
            .tolist()
        )
        raise ValueError(
            f"Tier-1 eligibility has invalid target QA status values; examples={bad}"
        )
    eligibility["tier1_target_qa_status"] = status.astype(str)
    eligibility["tier1_target_qa_pass"] = _normalize_strict_boolean(
        eligibility["tier1_target_qa_pass"], column="tier1_target_qa_pass"
    )
    expected_pass = eligibility["tier1_target_qa_status"].eq("pass")
    if not eligibility["tier1_target_qa_pass"].eq(expected_pass).all():
        raise ValueError(
            "Tier-1 eligibility target QA status and pass flag are inconsistent"
        )

    reasons = eligibility["tier1_target_qa_reasons"].fillna("").astype(str)
    missing_failure_reason = ~expected_pass & reasons.str.strip().eq("")
    if missing_failure_reason.any():
        bad_keys = (
            eligibility.loc[missing_failure_reason, list(TIER1_ELIGIBILITY_KEY)]
            .head(5)
            .to_dict("records")
        )
        raise ValueError(
            f"Tier-1 review/fail decisions require reason provenance; keys={bad_keys}"
        )
    eligibility["tier1_target_qa_reasons"] = reasons

    for column in ("tier1_contract_version", "tier1_config_name", "tier1_scope"):
        values = eligibility[column].astype("string")
        if values.isna().any() or values.str.strip().eq("").any():
            raise ValueError(f"Tier-1 eligibility has missing {column} provenance")

    expected_sector_tic_key = pd.Series(
        [
            f"s{int(sector):04d}-tic{int(tic):016d}"
            for sector, tic in zip(
                eligibility["sector"], eligibility["tic"], strict=True
            )
        ],
        index=eligibility.index,
    )
    if not eligibility["sector_tic_key"].astype(str).eq(expected_sector_tic_key).all():
        raise ValueError(
            "Tier-1 eligibility sector_tic_key does not match its (sector, TIC) columns"
        )

    evidence = eligibility.loc[
        :,
        [*TIER1_ELIGIBILITY_KEY, *TIER1_ELIGIBILITY_PROVENANCE_FIELDS],
    ]
    joined = work.merge(
        evidence,
        on=list(TIER1_ELIGIBILITY_KEY),
        how="left",
        validate="many_to_one",
        indicator="_tier1_eligibility_merge",
        sort=False,
    )
    missing_evidence = joined["_tier1_eligibility_merge"].ne("both")
    if missing_evidence.any():
        missing_keys = (
            joined.loc[missing_evidence, list(TIER1_ELIGIBILITY_KEY)]
            .drop_duplicates()
            .head(10)
            .to_dict("records")
        )
        raise ValueError(
            "Tier-1 target eligibility does not completely cover candidate "
            f"(sector, TIC) keys; missing={missing_keys}"
        )
    joined = joined.drop(columns="_tier1_eligibility_merge")

    before_keys = joined.drop_duplicates(list(TIER1_ELIGIBILITY_KEY))
    passing = joined["tier1_target_qa_pass"]
    filtered = joined.loc[passing].copy().reset_index(drop=True)
    after_keys = filtered.drop_duplicates(list(TIER1_ELIGIBILITY_KEY))
    summary = {
        "join_key": list(TIER1_ELIGIBILITY_KEY),
        "filter_expression": "tier1_target_qa_pass == True",
        "provenance_columns": list(TIER1_ELIGIBILITY_PROVENANCE_FIELDS),
        "n_eligibility_rows": int(len(eligibility)),
        "n_candidates_before": int(len(joined)),
        "n_candidates_after": int(len(filtered)),
        "n_candidate_tics_before": int(joined["tic"].nunique()),
        "n_candidate_tics_after": int(filtered["tic"].nunique()),
        "n_candidate_observations_before": int(len(before_keys)),
        "n_candidate_observations_after": int(len(after_keys)),
        "n_candidates_excluded": int((~passing).sum()),
        "candidate_status_counts": {
            str(key): int(value)
            for key, value in joined["tier1_target_qa_status"]
            .value_counts()
            .sort_index()
            .items()
        },
        "candidate_observation_status_counts": {
            str(key): int(value)
            for key, value in before_keys["tier1_target_qa_status"]
            .value_counts()
            .sort_index()
            .items()
        },
    }
    return filtered, summary


def normalize_a2v1_peak_candidates(
    peaks: pd.DataFrame,
    *,
    small_peaks_per_tic: int = 3,
    sector: int | None = None,
) -> pd.DataFrame:
    """Use ADP-small top-N peaks and attach the ADP-primary rank-1 context."""

    frame = peaks.copy()
    if sector is None:
        if "sector" not in frame:
            sector = 56
        else:
            sectors = (
                pd.to_numeric(frame["sector"], errors="coerce")
                .dropna()
                .astype(int)
                .unique()
            )
            if len(sectors) != 1:
                raise ValueError(
                    "teacher scoring requires exactly one sector; "
                    f"got {sorted(sectors.tolist())}"
                )
            sector = int(sectors[0])
    sector = int(sector)
    assert_adp_only_search_frame(frame)
    if "source_product_tag" not in frame:
        raise KeyError("A2v1 BLS peaks are missing source_product_tag")
    tags = set(frame["source_product_tag"].fillna("").astype(str))
    if tags != {"A2v1"}:
        raise ValueError(f"teacher scoring requires source_product_tag=A2v1; got {sorted(tags)}")
    rank = pd.to_numeric(frame["peak_rank"], errors="coerce")
    small = frame.loc[
        frame["aperture"].astype(str).eq(ADP_ONLY_APERTURES[0])
        & rank.between(1, int(small_peaks_per_tic))
        & frame["status"].astype(str).eq("ok")
    ].copy()
    if small.empty:
        raise ValueError("A2v1 peak table contains no valid ADP-small candidates")
    primary = frame.loc[
        frame["aperture"].astype(str).eq(ADP_ONLY_APERTURES[1])
        & rank.eq(1)
        & frame["status"].astype(str).eq("ok")
    ].copy()
    primary = primary.sort_values("tic", kind="stable").drop_duplicates("tic", keep="first")
    primary_fields = [column for column in PEAK_FIELDS if column in primary]
    primary = primary.loc[:, ["tic", *primary_fields]].rename(
        columns={column: f"adp_{column}" for column in primary_fields}
    )
    out = small.merge(primary, on="tic", how="left", validate="many_to_one")
    out = out.dropna(subset=["adp_period_d", "adp_t0_bjd", "adp_duration_min"])
    out["tic"] = pd.to_numeric(out["tic"], errors="coerce").astype(np.int64)
    out["sector"] = int(sector)
    out["rep_aperture"] = ADP_ONLY_APERTURES[0]
    out["rep_peak_rank"] = pd.to_numeric(out["peak_rank"], errors="coerce").astype(int)
    out["sde_max"] = pd.to_numeric(out["sde"], errors="coerce")
    out["anchor_aperture"] = ADP_ONLY_APERTURES[0]
    out["anchor_period_d"] = pd.to_numeric(out["period_d"], errors="coerce")
    out["anchor_t0_bjd"] = pd.to_numeric(out["t0_bjd"], errors="coerce")
    out["anchor_duration_min"] = pd.to_numeric(out["duration_min"], errors="coerce")
    out["anchor_sde"] = out["sde_max"]
    for column in PEAK_FIELDS:
        if column in out:
            out[f"adp_sml_{column}"] = pd.to_numeric(out[column], errors="coerce")
    relation = classify_period_relation(out["adp_sml_period_d"], out["adp_period_d"])
    out["aperture_period_relation"] = relation
    out["aperture_period_rel_delta"] = (
        np.abs(out["adp_period_d"] - out["adp_sml_period_d"]) / out["adp_sml_period_d"]
    )
    out["aperture_depth_ratio_primary_over_small"] = (
        pd.to_numeric(out.get("adp_depth"), errors="coerce")
        / pd.to_numeric(out.get("adp_sml_depth"), errors="coerce").replace(0, np.nan)
    )
    out["aperture_disagreement_flag"] = relation.eq("unrelated")
    out["n_apertures_agree"] = np.where(relation.isin({"exact", "harmonic"}), 2, 1)
    out["apertures_agree"] = np.where(
        relation.isin({"exact", "harmonic"}),
        ",".join(ADP_ONLY_APERTURES),
        ADP_ONLY_APERTURES[0],
    )
    out["review_id"] = [
        f"s{int(sector):04d}-A2v1-adp-small-{int(tic):016d}-r{int(peak_rank)}"
        for tic, peak_rank in zip(out["tic"], out["rep_peak_rank"])
    ]
    if out["review_id"].duplicated().any():
        raise ValueError("A2v1 candidate review_id values are not unique")
    out["source_kind"] = "real_candidate"
    out["source_product_tag"] = "A2v1"
    out["bls_search_branch"] = "A2v1_current_adp"
    out["adp_only_contract_version"] = ADP_ONLY_CONTRACT_VERSION
    out["input_contract_version"] = A2V1_TEACHER_INPUT_CONTRACT
    out["native_input_include"] = True
    out["native_group_path"] = [native_group_path(row) for row in out.to_dict("records")]
    out["candidate_key"] = out.apply(candidate_key, axis=1)
    return out.sort_values(["tic", "rep_peak_rank"], kind="stable").reset_index(drop=True)


def _peak_from_record(record: dict[str, Any], prefix: str) -> dict[str, Any]:
    return {
        key: record.get(f"{prefix}{key}", np.nan)
        for key in ("period_d", "t0_bjd", "duration_min", "depth", "depth_snr", "sde")
    }


def _initialize_metadata_quality_reference(
    reference: ExternalQualityReference,
) -> None:
    global _METADATA_QUALITY_REFERENCE
    _METADATA_QUALITY_REFERENCE = reference


def _measure_tic(payload: tuple[int, list[dict[str, Any]], str]) -> list[dict[str, Any]]:
    tic, records, compact_lc_path = payload
    lc = read_compact_lc_export(
        Path(compact_lc_path), tic=int(tic), columns=ADP_ONLY_APERTURES
    )
    if lc is None:
        return [
            {"review_id": str(record["review_id"]), "metadata_status": "missing_compact_lc"}
            for record in records
        ]
    if _METADATA_QUALITY_REFERENCE is None:
        raise RuntimeError("candidate metadata quality reference was not initialized")
    try:
        overlay = _METADATA_QUALITY_REFERENCE.apply(
            sector=int(lc.sector),
            camera=int(lc.cam),
            ccd=int(lc.ccd),
            cadenceno=lc.cadenceno,
            orbitid=lc.orbitid,
            internal_quality=lc.quality,
            context=f"candidate metadata TIC {tic}",
        )
        lc.quality = overlay.quality
    except Exception as exc:
        return [
            {
                "review_id": str(record["review_id"]),
                "metadata_status": "quality_overlay_error",
                "metadata_error": str(exc),
            }
            for record in records
        ]
    quality_audit = {
        f"metadata_quality_{name}": int(value)
        for name, value in overlay.counts.items()
    }
    output: list[dict[str, Any]] = []
    for record in records:
        anchor = _peak_from_record(record, "")
        own = {
            ADP_ONLY_APERTURES[0]: _peak_from_record(record, "adp_sml_"),
            ADP_ONLY_APERTURES[1]: _peak_from_record(record, "adp_"),
        }
        try:
            measured = measure_two_aperture_candidate_metadata(
                lc,
                anchor_peak=anchor,
                own_peaks=own,
                apertures=ADP_ONLY_APERTURES,
            )
            measured["review_id"] = str(record["review_id"])
            measured["metadata_status"] = "ok"
            measured.update(quality_audit)
        except Exception as exc:
            measured = {
                "review_id": str(record["review_id"]),
                "metadata_status": "error",
                "metadata_error": str(exc),
            }
        output.append(measured)
    return output


def enrich_candidate_metadata(
    candidates: pd.DataFrame,
    *,
    compact_lc_path: Path,
    cadence_reference_table: Path,
    cadence_reference_manifest: Path,
    workers: int = 1,
    progress_every: int = 500,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Measure training-compatible odd/even, trend, and aperture metadata."""

    work = candidates.copy()
    sectors = sorted(
        set(pd.to_numeric(work["sector"], errors="raise").astype(int).tolist())
    )
    if sectors != [56]:
        raise ValueError(
            "the current quality-bound candidate metadata contract requires S56"
        )
    quality_reference = load_external_quality_reference(
        table_path=cadence_reference_table,
        manifest_path=cadence_reference_manifest,
        sector=56,
        expected_orbits=S56_EXPECTED_ORBITS,
        expected_detectors=S56_EXPECTED_DETECTORS,
    )
    payloads = [
        (int(tic), group.to_dict("records"), str(compact_lc_path))
        for tic, group in work.groupby("tic", sort=True)
    ]
    rows: list[dict[str, Any]] = []
    if int(workers) <= 1:
        _initialize_metadata_quality_reference(quality_reference)
        iterator = map(_measure_tic, payloads)
        executor = None
    else:
        executor = ProcessPoolExecutor(
            max_workers=int(workers),
            initializer=_initialize_metadata_quality_reference,
            initargs=(quality_reference,),
        )
        iterator = executor.map(_measure_tic, payloads, chunksize=1)
    try:
        for index, batch in enumerate(iterator, start=1):
            rows.extend(batch)
            if progress_every > 0 and index % int(progress_every) == 0:
                print(
                    f"[teacher-candidates] metadata {index:,}/{len(payloads):,} TICs",
                    flush=True,
                )
    finally:
        if executor is not None:
            executor.shutdown(wait=True)
    quality_reference.assert_unchanged()
    metrics = pd.DataFrame(rows)
    if metrics["review_id"].duplicated().any():
        raise RuntimeError("candidate metadata contains duplicate review_id values")
    overlap = [
        column for column in metrics.columns if column != "review_id" and column in work.columns
    ]
    work = work.drop(columns=overlap).merge(
        metrics, on="review_id", how="left", validate="one_to_one"
    )
    status = work["metadata_status"].fillna("missing").astype(str)
    quality_count_columns = [
        f"metadata_quality_{name}"
        for name in (
            "n_cad_total",
            "n_cad_internal_bad",
            "n_cad_external_bad",
            "n_cad_external_only_bad",
            "n_cad_effective_bad",
        )
    ]
    unique_tic_quality = work.drop_duplicates("tic", keep="first")
    aggregate_quality_counts = {
        name.removeprefix("metadata_quality_"): int(
            pd.to_numeric(unique_tic_quality[name], errors="coerce").fillna(0).sum()
        )
        for name in quality_count_columns
        if name in unique_tic_quality
    }
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "n_candidates": int(len(work)),
        "n_unique_tics": int(work["tic"].nunique()),
        "compact_lc_path": str(compact_lc_path),
        "workers": int(workers),
        "external_quality": {
            **quality_reference.provenance,
            "counts_over_unique_tics": aggregate_quality_counts,
        },
        "metadata_status_counts": {
            str(key): int(value) for key, value in status.value_counts().sort_index().items()
        },
        "passed": bool(status.eq("ok").all()),
    }
    return work, summary


__all__ = [
    "enrich_candidate_metadata",
    "filter_tier1_eligible_candidates",
    "normalize_a2v1_peak_candidates",
    "validate_tier1_enrichment_gate",
]
