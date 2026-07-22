#!/usr/bin/env python3
"""Build a full-S56 real-candidate BLS peak table from the ADP pair only."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Mapping

import h5py
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from twirl.io.compact_export import read_compact_lc_export  # noqa: E402
from twirl.search.a2v1_bls_contract import (  # noqa: E402
    A2V1_TEACHER_BLS_SEARCH_CONTRACT,
    approved_a2v1_teacher_bls_config,
    bls_config_sha256,
)
from twirl.search.bls import BLSConfig, run_bls_on_lc  # noqa: E402
from twirl.vetting.adp_only import (  # noqa: E402
    ADP_ONLY_APERTURES,
    ADP_ONLY_CONTRACT_VERSION,
    validate_adp_only_apertures,
)
from twirl.vetting.recovery50_teacher import json_default, write_table  # noqa: E402


DEFAULT_COMPACT_LC = (
    REPO_ROOT
    / "data_local/stage3_injections/s56_twirlfs_v2_lc_export/"
    "s56_twirlfs_v2_adp_lc_export_pdo.h5"
)
DEFAULT_OUT_DIR = REPO_ROOT / "reports/stage5_validation/s56_adp_real_bls_peaks"

EXTERNAL_QUALITY_POLICY_CONTRACT = (
    "a2v1_bls_internal_or_authoritative_external_quality_v1"
)
CADENCE_REFERENCE_CONTRACT_TEMPLATE = "s{sector}_a2v1_cadence_reference_v1"
CADENCE_REFERENCE_COLUMNS = (
    "sector",
    "orbitid",
    "camera",
    "ccd",
    "cadenceno",
    "spoc_quality",
    "qlp_quality",
    "external_quality",
)
EXPECTED_CADENCE_AUTHORITY = "qlp_cam_quat"
EXPECTED_QUALITY_AUTHORITY = "spoc_and_qlp_quality_flags"
EXPECTED_QUALITY_COMPOSITION = {
    "external_quality": "spoc_quality | (qlp_quality << 30)",
    "qlp_quality_raw_values": [0, 1],
    "qlp_quality_external_bit": 30,
}

# Populated once per worker by ``_initialize_external_quality_worker``.  The
# authoritative cadence map is deliberately not repeated in every task
# payload: a full S56 run has tens of thousands of targets but only sixteen
# detector maps.
_EXTERNAL_QUALITY_BY_DETECTOR: dict[
    tuple[int, int, int], tuple[np.ndarray, np.ndarray, np.ndarray]
] | None = None
_EXTERNAL_QUALITY_PROVENANCE: dict[str, str] | None = None


def _sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def _read_table(path: Path) -> pd.DataFrame:
    suffix = Path(path).suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix == ".parquet":
        return pd.read_parquet(path)
    raise ValueError("cadence reference must end in .csv or .parquet")


def _integer_column(frame: pd.DataFrame, column: str) -> pd.Series:
    values = pd.to_numeric(frame[column], errors="raise")
    array = values.to_numpy(dtype=np.float64)
    if not np.all(np.isfinite(array)) or not np.all(array == np.floor(array)):
        raise ValueError(f"cadence-reference {column} must contain finite integers")
    return values.astype(np.int64)


def _canonical_json_sha256(payload: Mapping[str, Any]) -> str:
    text = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _validate_source_hash_declarations(manifest: Mapping[str, Any]) -> str:
    """Validate and fingerprint the authority hashes declared by the manifest.

    The manifest itself and the cadence table are rehashed before and after a
    BLS run.  This canonical fingerprint additionally binds the complete map
    of source-file hashes that the cadence-reference builder verified when it
    published the immutable evidence pair.
    """

    declared = manifest.get("source_file_sha256")
    if not isinstance(declared, Mapping) or not declared:
        raise ValueError("cadence-reference source_file_sha256 is missing or empty")
    normalized: dict[str, str] = {}
    for raw_path, raw_digest in declared.items():
        path = str(raw_path).strip()
        digest = str(raw_digest).lower()
        if not path:
            raise ValueError("cadence-reference source hash has an empty path")
        if len(digest) != 64 or any(
            char not in "0123456789abcdef" for char in digest
        ):
            raise ValueError(
                f"cadence-reference source hash is invalid for {path!r}"
            )
        if path in normalized:
            raise ValueError(f"duplicate cadence-reference source path: {path}")
        normalized[path] = digest

    sources = manifest.get("sources")
    if not isinstance(sources, list) or not sources:
        raise ValueError("cadence-reference sources inventory is missing or empty")
    observed_paths: set[str] = set()
    for index, source in enumerate(sources):
        if not isinstance(source, Mapping):
            raise ValueError(f"cadence-reference source {index} is not an object")
        path = str(source.get("path", "")).strip()
        digest = str(source.get("sha256", "")).lower()
        if not path or path in observed_paths:
            raise ValueError("cadence-reference sources contain an empty/duplicate path")
        observed_paths.add(path)
        if normalized.get(path) != digest:
            raise ValueError(
                "cadence-reference sources inventory disagrees with "
                f"source_file_sha256 for {path!r}"
            )
    if observed_paths != set(normalized):
        raise ValueError(
            "cadence-reference sources inventory does not exactly cover "
            "source_file_sha256"
        )
    return _canonical_json_sha256(normalized)


def load_external_quality_reference(
    *,
    table_path: Path,
    manifest_path: Path,
    sector: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Load and validate one authoritative detector-cadence quality overlay."""

    table_path = Path(table_path)
    manifest_path = Path(manifest_path)
    table_sha256 = _sha256(table_path)
    manifest_sha256 = _sha256(manifest_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(manifest, dict):
        raise ValueError("cadence-reference manifest must contain a JSON object")

    required_manifest = {
        "contract_version",
        "sector",
        "cadence_authority",
        "quality_authority",
        "quality_composition",
        "table_sha256",
        "n_rows",
        "detectors",
        "orbits",
        "source_file_sha256",
        "sources",
        "n_nonzero_spoc_quality",
        "n_nonzero_qlp_quality",
        "n_nonzero_external_quality",
    }
    missing_manifest = sorted(required_manifest - set(manifest))
    if missing_manifest:
        raise ValueError(
            f"cadence-reference manifest is missing fields: {missing_manifest}"
        )
    expected_contract = CADENCE_REFERENCE_CONTRACT_TEMPLATE.format(sector=int(sector))
    if str(manifest["contract_version"]) != expected_contract:
        raise ValueError("cadence-reference manifest contract mismatch")
    if int(manifest["sector"]) != int(sector):
        raise ValueError("cadence-reference manifest sector mismatch")
    if str(manifest["cadence_authority"]) != EXPECTED_CADENCE_AUTHORITY:
        raise ValueError("cadence-reference cadence authority mismatch")
    if str(manifest["quality_authority"]) != EXPECTED_QUALITY_AUTHORITY:
        raise ValueError("cadence-reference quality authority mismatch")
    if manifest["quality_composition"] != EXPECTED_QUALITY_COMPOSITION:
        raise ValueError("cadence-reference external-quality composition mismatch")
    if str(manifest["table_sha256"]) != table_sha256:
        raise ValueError("cadence-reference table hash mismatch")

    frame = _read_table(table_path)
    if tuple(frame.columns) != CADENCE_REFERENCE_COLUMNS:
        raise ValueError(
            "cadence-reference columns must be exactly "
            f"{CADENCE_REFERENCE_COLUMNS}; observed {tuple(frame.columns)}"
        )
    if len(frame) != int(manifest["n_rows"]):
        raise ValueError("cadence-reference row count disagrees with manifest")
    if frame.empty:
        raise ValueError("cadence-reference table is empty")
    for column in CADENCE_REFERENCE_COLUMNS:
        frame[column] = _integer_column(frame, column)
    if set(frame["sector"].astype(int)) != {int(sector)}:
        raise ValueError("cadence-reference table sector mismatch")
    if (
        not frame["camera"].between(1, 4).all()
        or not frame["ccd"].between(1, 4).all()
    ):
        raise ValueError("cadence-reference camera/ccd values must be in 1..4")
    if (frame["orbitid"] <= 0).any() or (frame["cadenceno"] < 0).any():
        raise ValueError("cadence-reference orbit/cadence values are invalid")
    if (frame[["spoc_quality", "qlp_quality", "external_quality"]] < 0).any().any():
        raise ValueError("cadence-reference quality values must be nonnegative")
    if not set(frame["qlp_quality"].astype(int)).issubset({0, 1}):
        raise ValueError("cadence-reference qlp_quality must contain only 0/1")
    expected_external = np.bitwise_or(
        frame["spoc_quality"].to_numpy(dtype=np.int64),
        np.left_shift(frame["qlp_quality"].to_numpy(dtype=np.int64), 30),
    )
    if not np.array_equal(
        expected_external, frame["external_quality"].to_numpy(dtype=np.int64)
    ):
        raise ValueError("cadence-reference external_quality composition is invalid")

    key = ["sector", "camera", "ccd", "cadenceno"]
    duplicate = frame.duplicated(key, keep=False)
    if duplicate.any():
        examples = frame.loc[duplicate, [*key, "orbitid"]].head(10).to_dict("records")
        raise ValueError(
            f"cadence-reference has duplicate detector-cadence mappings: {examples}"
        )
    observed_detectors = sorted(
        {
            f"cam{int(camera)}_ccd{int(ccd)}"
            for camera, ccd in zip(frame["camera"], frame["ccd"], strict=True)
        }
    )
    if sorted(str(value) for value in manifest["detectors"]) != observed_detectors:
        raise ValueError("cadence-reference detector inventory mismatch")
    observed_orbits = sorted(set(frame["orbitid"].astype(int)))
    if sorted(int(value) for value in manifest["orbits"]) != observed_orbits:
        raise ValueError("cadence-reference orbit inventory mismatch")
    count_fields = {
        "n_nonzero_spoc_quality": "spoc_quality",
        "n_nonzero_qlp_quality": "qlp_quality",
        "n_nonzero_external_quality": "external_quality",
    }
    for manifest_field, column in count_fields.items():
        observed = int(np.count_nonzero(frame[column].to_numpy(dtype=np.int64)))
        if int(manifest[manifest_field]) != observed:
            raise ValueError(
                f"cadence-reference {manifest_field} disagrees with the table"
            )

    source_hashes_sha256 = _validate_source_hash_declarations(manifest)
    frame = frame.sort_values(
        ["sector", "camera", "ccd", "cadenceno"], kind="stable"
    ).reset_index(drop=True)
    return frame, {
        "contract_version": expected_contract,
        "cadence_authority": EXPECTED_CADENCE_AUTHORITY,
        "quality_authority": EXPECTED_QUALITY_AUTHORITY,
        "table_sha256": table_sha256,
        "manifest_sha256": manifest_sha256,
        "source_file_sha256_declaration_sha256": source_hashes_sha256,
    }


def _reference_worker_payload(
    frame: pd.DataFrame,
) -> dict[tuple[int, int, int], tuple[np.ndarray, np.ndarray, np.ndarray]]:
    payload: dict[
        tuple[int, int, int], tuple[np.ndarray, np.ndarray, np.ndarray]
    ] = {}
    for (sector, camera, ccd), detector in frame.groupby(
        ["sector", "camera", "ccd"], sort=True
    ):
        ordered = detector.sort_values("cadenceno", kind="stable")
        payload[(int(sector), int(camera), int(ccd))] = (
            ordered["cadenceno"].to_numpy(dtype=np.int64),
            ordered["orbitid"].to_numpy(dtype=np.int64),
            ordered["external_quality"].to_numpy(dtype=np.int64),
        )
    return payload


def _initialize_external_quality_worker(
    reference_by_detector: dict[
        tuple[int, int, int], tuple[np.ndarray, np.ndarray, np.ndarray]
    ],
    provenance: dict[str, str],
) -> None:
    global _EXTERNAL_QUALITY_BY_DETECTOR, _EXTERNAL_QUALITY_PROVENANCE
    _EXTERNAL_QUALITY_BY_DETECTOR = reference_by_detector
    _EXTERNAL_QUALITY_PROVENANCE = provenance


def _apply_external_quality(lc: Any) -> dict[str, int]:
    if _EXTERNAL_QUALITY_BY_DETECTOR is None:
        raise RuntimeError("external-quality worker state was not initialized")
    lengths = {
        "time": len(lc.time),
        "cadenceno": len(lc.cadenceno),
        "orbitid": len(lc.orbitid),
        "quality": len(lc.quality),
    }
    if len(set(lengths.values())) != 1:
        raise ValueError(f"TIC {lc.tic}: compact cadence arrays disagree: {lengths}")
    cadences = np.asarray(lc.cadenceno, dtype=np.int64)
    orbits = np.asarray(lc.orbitid, dtype=np.int64)
    internal_quality = np.asarray(lc.quality, dtype=np.int64)
    if len(cadences) != len(np.unique(cadences)):
        raise ValueError(f"TIC {lc.tic}: compact cadenceno values are not unique")
    detector_key = (int(lc.sector), int(lc.cam), int(lc.ccd))
    reference = _EXTERNAL_QUALITY_BY_DETECTOR.get(detector_key)
    if reference is None:
        raise ValueError(
            f"TIC {lc.tic}: no external-quality detector mapping for {detector_key}"
        )
    reference_cadence, reference_orbit, reference_quality = reference
    positions = np.searchsorted(reference_cadence, cadences)
    in_bounds = positions < len(reference_cadence)
    matched = np.zeros(len(cadences), dtype=bool)
    matched[in_bounds] = reference_cadence[positions[in_bounds]] == cadences[in_bounds]
    if not matched.all():
        missing = cadences[~matched][:20].astype(int).tolist()
        raise ValueError(
            f"TIC {lc.tic}: external-quality coverage is missing compact cadences "
            f"for {detector_key}: {missing}"
        )
    mapped_orbits = reference_orbit[positions]
    orbit_mismatch = mapped_orbits != orbits
    if orbit_mismatch.any():
        examples = [
            {
                "cadenceno": int(cadences[index]),
                "compact_orbitid": int(orbits[index]),
                "reference_orbitid": int(mapped_orbits[index]),
            }
            for index in np.flatnonzero(orbit_mismatch)[:20]
        ]
        raise ValueError(
            f"TIC {lc.tic}: external-quality orbit mapping mismatch: {examples}"
        )
    external_quality = reference_quality[positions]
    internal_bad = internal_quality != 0
    external_bad = external_quality != 0
    effective_bad = internal_bad | external_bad
    lc.quality = effective_bad.astype(np.int32)
    return {
        "n_cad_internal_bad": int(np.count_nonzero(internal_bad)),
        "n_cad_external_bad": int(np.count_nonzero(external_bad)),
        "n_cad_external_only_bad": int(np.count_nonzero(external_bad & ~internal_bad)),
        "n_cad_effective_bad": int(np.count_nonzero(effective_bad)),
    }


def _target_tics(path: Path) -> list[int]:
    with h5py.File(path, "r") as h5:
        if "targets" not in h5:
            raise KeyError(f"compact export has no /targets group: {path}")
        return sorted(int(key) for key in h5["targets"].keys())


def _result_rows(result: Any) -> list[dict[str, Any]]:
    base = {
        "tic": int(result.tic),
        "sector": int(result.sector),
        "cam": int(result.cam),
        "ccd": int(result.ccd),
        "tmag": float(result.tmag),
        "aperture": str(result.aperture),
        "n_cad_total": int(result.n_cad_total),
        "n_cad_quality": int(result.n_cad_quality or 0),
        "n_cad_kept": int(result.n_cad_kept),
        "n_cad_edge_trimmed": int(result.n_cad_edge_trimmed),
        "n_cad_sigma_clipped": int(result.n_cad_sigma_clipped),
        "dropout_frac": float(result.dropout_frac),
        "quality_dropout_frac": float(result.quality_dropout_frac or 0.0),
        "n_orbits": int(result.n_orbits),
        "baseline_d": float(result.baseline_d),
        "status": str(result.status),
        "bls_search_branch": "current_adp",
        "adp_only_contract_version": ADP_ONLY_CONTRACT_VERSION,
    }
    if not result.peaks:
        return [
            {
                **base,
                "peak_rank": 0,
                "period_d": np.nan,
                "t0_bjd": np.nan,
                "duration_min": np.nan,
                "depth": np.nan,
                "depth_snr": np.nan,
                "sde": np.nan,
                "log_power": np.nan,
            }
        ]
    return [{**base, **asdict(peak)} for peak in result.peaks]


def _process_target(payload: tuple[int, str, dict[str, Any]]) -> list[dict[str, Any]]:
    tic, compact_lc_s, cfg_payload = payload
    compact_lc = Path(compact_lc_s)
    cfg = BLSConfig(
        apertures=tuple(cfg_payload.get("apertures", ADP_ONLY_APERTURES)),
        n_periods=int(cfg_payload["n_periods"]),
        n_peaks=int(cfg_payload["n_peaks"]),
        p_min_d=float(cfg_payload["p_min_d"]),
        p_max_cap_d=float(cfg_payload["p_max_cap_d"]),
        max_period_fraction=float(cfg_payload["max_period_fraction"]),
        durations_min=tuple(
            float(value)
            for value in cfg_payload.get("durations_min", BLSConfig.durations_min)
        ),
        period_mask_frac=float(
            cfg_payload.get("period_mask_frac", BLSConfig.period_mask_frac)
        ),
        period_bin_edges=tuple(
            float(value) for value in cfg_payload.get("period_bin_edges", ())
        ),
        max_peaks_per_period_bin=int(
            cfg_payload.get(
                "max_peaks_per_period_bin", BLSConfig.max_peaks_per_period_bin
            )
        ),
        min_cadences=int(cfg_payload.get("min_cadences", BLSConfig.min_cadences)),
        sigma_clip=float(cfg_payload["sigma_clip"]),
        orbit_edge_trim_d=float(cfg_payload["orbit_edge_trim_d"]),
    )
    config_sha256 = bls_config_sha256(cfg_payload)
    search_contract = (
        A2V1_TEACHER_BLS_SEARCH_CONTRACT
        if cfg_payload == approved_a2v1_teacher_bls_config()
        else "custom"
    )
    config_provenance = {
        "bls_n_periods": int(cfg.n_periods),
        "bls_n_peaks": int(cfg.n_peaks),
        "bls_p_min_d": float(cfg.p_min_d),
        "bls_p_max_cap_d": float(cfg.p_max_cap_d),
        "bls_max_period_fraction": float(cfg.max_period_fraction),
        "bls_sigma_clip": float(cfg.sigma_clip),
        "bls_orbit_edge_trim_d": float(cfg.orbit_edge_trim_d),
    }
    lc = read_compact_lc_export(compact_lc, tic=tic, columns=ADP_ONLY_APERTURES)
    if lc is None:
        raise RuntimeError(
            f"TIC {tic}: compact product could not supply the locked ADP pair"
        )
    quality_counts = _apply_external_quality(lc)
    if _EXTERNAL_QUALITY_PROVENANCE is None:
        raise RuntimeError("external-quality provenance was not initialized")
    rows: list[dict[str, Any]] = []
    for aperture in ADP_ONLY_APERTURES:
        result = run_bls_on_lc(lc, cfg, aperture=aperture)
        current = _result_rows(result)
        for row in current:
            row["source_product_tag"] = str(cfg_payload.get("source_product_tag", ""))
            row.update(config_provenance)
            row["bls_search_contract_version"] = search_contract
            row["bls_config_sha256"] = config_sha256
            row.update(quality_counts)
            row["external_quality_policy_contract"] = (
                EXTERNAL_QUALITY_POLICY_CONTRACT
            )
            row["cadence_reference_sha256"] = _EXTERNAL_QUALITY_PROVENANCE[
                "table_sha256"
            ]
            row["cadence_reference_manifest_sha256"] = (
                _EXTERNAL_QUALITY_PROVENANCE["manifest_sha256"]
            )
        rows.extend(current)
    return rows


def build_peak_table(
    *,
    sector: int = 56,
    compact_lc: Path,
    cadence_reference: Path,
    cadence_reference_manifest: Path,
    out_dir: Path,
    workers: int,
    n_periods: int,
    n_peaks: int,
    max_targets: int | None,
    progress_every: int,
    shard_index: int = 0,
    n_shards: int = 1,
    resume: bool = False,
    source_product_tag: str = "",
) -> dict[str, Any]:
    validate_adp_only_apertures(ADP_ONLY_APERTURES)
    compact_lc = Path(compact_lc)
    cadence_reference = Path(cadence_reference)
    cadence_reference_manifest = Path(cadence_reference_manifest)
    input_sha256 = {
        "compact_lc": _sha256(compact_lc),
        "cadence_reference": _sha256(cadence_reference),
        "cadence_reference_manifest": _sha256(cadence_reference_manifest),
    }
    reference, reference_provenance = load_external_quality_reference(
        table_path=cadence_reference,
        manifest_path=cadence_reference_manifest,
        sector=int(sector),
    )
    if reference_provenance["table_sha256"] != input_sha256["cadence_reference"]:
        raise RuntimeError("cadence-reference hash changed while it was validated")
    if (
        reference_provenance["manifest_sha256"]
        != input_sha256["cadence_reference_manifest"]
    ):
        raise RuntimeError("cadence-reference manifest changed while it was validated")
    reference_payload = _reference_worker_payload(reference)
    out_dir.mkdir(parents=True, exist_ok=True)
    tics = _target_tics(compact_lc)
    if max_targets is not None:
        tics = tics[: max(0, int(max_targets))]
    if n_shards < 1 or shard_index < 0 or shard_index >= n_shards:
        raise ValueError("shard_index must satisfy 0 <= shard_index < n_shards")
    n_targets_total = len(tics)
    tics = tics[int(shard_index) :: int(n_shards)]
    suffix = f"_{int(shard_index):03d}" if int(n_shards) > 1 else ""
    output_path = out_dir / f"real_adp_bls_peaks{suffix}.parquet"
    summary_path = out_dir / f"summary{suffix}.json"
    cfg_payload = approved_a2v1_teacher_bls_config()
    cfg_payload.update(
        {
            "n_periods": int(n_periods),
            "n_peaks": int(n_peaks),
            "source_product_tag": str(source_product_tag),
        }
    )
    cfg_sha256 = bls_config_sha256(cfg_payload)
    search_contract = (
        A2V1_TEACHER_BLS_SEARCH_CONTRACT
        if cfg_payload == approved_a2v1_teacher_bls_config()
        else "custom"
    )
    if resume and output_path.exists() and summary_path.exists():
        summary = json.loads(summary_path.read_text())
        if (
            summary.get("contract_version") == ADP_ONLY_CONTRACT_VERSION
            and summary.get("external_quality_policy_contract")
            == EXTERNAL_QUALITY_POLICY_CONTRACT
            and summary.get("compact_lc") == str(compact_lc)
            and summary.get("compact_lc_sha256") == input_sha256["compact_lc"]
            and summary.get("cadence_reference") == str(cadence_reference)
            and summary.get("cadence_reference_sha256")
            == input_sha256["cadence_reference"]
            and summary.get("cadence_reference_manifest")
            == str(cadence_reference_manifest)
            and summary.get("cadence_reference_manifest_sha256")
            == input_sha256["cadence_reference_manifest"]
            and summary.get("cadence_reference_contract_version")
            == reference_provenance["contract_version"]
            and summary.get("cadence_reference_cadence_authority")
            == reference_provenance["cadence_authority"]
            and summary.get("cadence_reference_quality_authority")
            == reference_provenance["quality_authority"]
            and summary.get("cadence_reference_source_hashes_sha256")
            == reference_provenance["source_file_sha256_declaration_sha256"]
            and int(summary.get("shard_index", -1)) == int(shard_index)
            and int(summary.get("n_shards", -1)) == int(n_shards)
            and int(summary.get("n_targets", -1)) == len(tics)
            and int(summary.get("n_targets_total", -1)) == n_targets_total
            and summary.get("source_product_tag") == str(source_product_tag)
            and summary.get("config") == cfg_payload
            and summary.get("bls_config_sha256") == cfg_sha256
            and summary.get("bls_search_contract_version") == search_contract
            and summary.get("peak_table_sha256") == _sha256(output_path)
        ):
            print(json.dumps(summary, indent=2, sort_keys=True))
            return summary
    payloads = [(tic, str(compact_lc), cfg_payload) for tic in tics]
    rows: list[dict[str, Any]] = []
    workers = max(1, int(workers))
    if workers == 1:
        _initialize_external_quality_worker(reference_payload, reference_provenance)
        iterator = map(_process_target, payloads)
        executor = None
    else:
        executor = ProcessPoolExecutor(
            max_workers=workers,
            initializer=_initialize_external_quality_worker,
            initargs=(reference_payload, reference_provenance),
        )
        iterator = executor.map(_process_target, payloads, chunksize=1)
    try:
        for index, batch in enumerate(iterator, start=1):
            rows.extend(batch)
            if progress_every > 0 and index % int(progress_every) == 0:
                print(f"[adp-real-bls] processed {index:,}/{len(payloads):,}", flush=True)
    finally:
        if executor is not None:
            executor.shutdown(wait=True)

    final_input_sha256 = {
        "compact_lc": _sha256(compact_lc),
        "cadence_reference": _sha256(cadence_reference),
        "cadence_reference_manifest": _sha256(cadence_reference_manifest),
    }
    if final_input_sha256 != input_sha256:
        raise RuntimeError(
            "BLS input source changed during the run; refusing to publish peaks"
        )

    peaks = pd.DataFrame(rows)
    output_path = write_table(peaks, output_path)
    output_sha256 = _sha256(output_path)
    status = peaks.get("status", pd.Series(dtype=str)).fillna("").astype(str)
    valid = status.eq("ok") & pd.to_numeric(
        peaks.get("peak_rank"), errors="coerce"
    ).gt(0)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "sector": int(sector),
        "contract_version": ADP_ONLY_CONTRACT_VERSION,
        "bls_search_contract_version": search_contract,
        "bls_config_sha256": cfg_sha256,
        "external_quality_policy_contract": EXTERNAL_QUALITY_POLICY_CONTRACT,
        "compact_lc": str(compact_lc),
        "compact_lc_sha256": input_sha256["compact_lc"],
        "cadence_reference": str(cadence_reference),
        "cadence_reference_sha256": input_sha256["cadence_reference"],
        "cadence_reference_manifest": str(cadence_reference_manifest),
        "cadence_reference_manifest_sha256": input_sha256[
            "cadence_reference_manifest"
        ],
        "cadence_reference_contract_version": reference_provenance[
            "contract_version"
        ],
        "cadence_reference_cadence_authority": reference_provenance[
            "cadence_authority"
        ],
        "cadence_reference_quality_authority": reference_provenance[
            "quality_authority"
        ],
        "cadence_reference_source_hashes_sha256": reference_provenance[
            "source_file_sha256_declaration_sha256"
        ],
        "out_dir": str(out_dir),
        "apertures": list(ADP_ONLY_APERTURES),
        "n_targets": int(len(tics)),
        "n_targets_total": int(n_targets_total),
        "n_rows": int(len(peaks)),
        "n_unique_tics": int(peaks["tic"].nunique()) if "tic" in peaks else 0,
        "n_valid_peak_rows": int(valid.sum()),
        "n_periods": int(n_periods),
        "n_peaks": int(n_peaks),
        "workers": int(workers),
        "shard_index": int(shard_index),
        "n_shards": int(n_shards),
        "source_product_tag": str(source_product_tag),
        "peak_table_sha256": output_sha256,
        "config": cfg_payload,
        "status_counts": {
            str(k): int(v) for k, v in status.value_counts().sort_index().items()
        },
        "aperture_counts": {
            str(k): int(v)
            for k, v in peaks.get("aperture", pd.Series(dtype=str))
            .fillna("")
            .astype(str)
            .value_counts()
            .sort_index()
            .items()
        },
        "outputs": {
            "peak_table": str(output_path),
            "summary": str(summary_path),
        },
    }
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True, default=json_default) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True, default=json_default))
    return summary


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sector", type=int, default=56)
    parser.add_argument("--compact-lc", type=Path, default=DEFAULT_COMPACT_LC)
    parser.add_argument("--cadence-reference", type=Path, required=True)
    parser.add_argument("--cadence-reference-manifest", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--n-periods", type=int, default=50_000)
    parser.add_argument("--n-peaks", type=int, default=10)
    parser.add_argument("--max-targets", type=int, default=None)
    parser.add_argument("--progress-every", type=int, default=100)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--n-shards", type=int, default=1)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--source-product-tag",
        default="",
        help="Set to A2v1 when the compact export came from the production A2v1 FITS tree.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    build_peak_table(
        sector=args.sector,
        compact_lc=args.compact_lc,
        cadence_reference=args.cadence_reference,
        cadence_reference_manifest=args.cadence_reference_manifest,
        out_dir=args.out_dir,
        workers=args.workers,
        n_periods=args.n_periods,
        n_peaks=args.n_peaks,
        max_targets=args.max_targets,
        progress_every=args.progress_every,
        shard_index=args.shard_index,
        n_shards=args.n_shards,
        resume=args.resume,
        source_product_tag=args.source_product_tag,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
