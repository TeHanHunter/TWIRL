#!/usr/bin/env python3
"""Merge and provenance-verify resumable ADP-only real-BLS shards."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd

from twirl.vetting.adp_only import ADP_ONLY_APERTURES


_INVARIANT_SUMMARY_FIELDS = (
    "sector",
    "contract_version",
    "bls_search_contract_version",
    "bls_config_sha256",
    "external_quality_policy_contract",
    "compact_lc",
    "compact_lc_sha256",
    "cadence_reference",
    "cadence_reference_sha256",
    "cadence_reference_manifest",
    "cadence_reference_manifest_sha256",
    "cadence_reference_contract_version",
    "cadence_reference_cadence_authority",
    "cadence_reference_quality_authority",
    "cadence_reference_source_hashes_sha256",
    "authority_exclusion_policy_contract",
    "authority_exclusion_external_bit",
    "authority_exclusions_sha256",
    "n_authority_exclusions",
    "apertures",
    "n_targets_total",
    "n_periods",
    "n_peaks",
    "source_product_tag",
    "config",
)
_TARGET_SELECTION_INVARIANT_FIELDS = (
    "target_selection_contract_version",
    "target_allowlist",
    "target_allowlist_sha256",
    "target_allowlist_count",
    "target_allowlist_tics_sha256",
)
_ORBITID_INVARIANT_FIELDS = (
    "orbitid_policy",
    "orbitid_reconciliation_contract_version",
)
_ORBITID_POLICIES = ("strict", "reference_by_cadence")
_ORBITID_RECONCILIATION_CONTRACT_VERSION = "a2v1_compact_orbitid_reconciliation_v1"
_TARGET_SELECTION_CONTRACT_VERSION = "a2v1_bls_target_allowlist_v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical_json_sha256(payload: Mapping[str, Any]) -> str:
    text = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _tic_inventory_sha256(tics: set[int] | list[int]) -> str:
    normalized = sorted(int(value) for value in tics)
    if any(value <= 0 for value in normalized) or len(normalized) != len(
        set(normalized)
    ):
        raise ValueError("TIC inventory must contain unique positive integers")
    payload = "".join(f"{value}\n" for value in normalized)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def _read_target_allowlist(path: Path) -> list[int]:
    path = Path(path)
    if path.suffix.lower() == ".csv":
        frame = pd.read_csv(path, dtype=str, keep_default_na=False)
        if "tic" not in frame:
            raise ValueError("target allowlist CSV must contain a 'tic' column")
        raw_values = frame["tic"].astype(str).str.strip().tolist()
    elif path.suffix.lower() in {".txt", ".list"}:
        raw_values = []
        for raw_line in path.read_text(encoding="utf-8").splitlines():
            value = raw_line.strip()
            if not value or value.startswith("#"):
                continue
            if not raw_values and value.lower() == "tic":
                continue
            raw_values.append(value)
    else:
        raise ValueError("target allowlist must end in .csv, .txt, or .list")
    if not raw_values:
        raise ValueError("target allowlist is empty")
    if any(
        not value or not value.isascii() or not value.isdigit() for value in raw_values
    ):
        raise ValueError("target allowlist contains a non-integer TIC value")
    tics = [int(value) for value in raw_values]
    if any(value <= 0 or value > np.iinfo(np.int64).max for value in tics) or len(
        tics
    ) != len(set(tics)):
        raise ValueError(
            "target allowlist must contain unique positive int64 TIC values"
        )
    return sorted(tics)


def _orbitid_summary_from_rows(
    frame: pd.DataFrame,
    *,
    orbitid_policy: str,
) -> dict[str, Any]:
    if orbitid_policy not in _ORBITID_POLICIES:
        raise ValueError(f"invalid orbit-ID policy: {orbitid_policy!r}")
    count_columns = (
        "n_cad_orbitid_reference_matched",
        "n_cad_orbitid_mismatch",
        "n_cad_orbitid_corrected",
    )
    required = {
        "tic",
        "orbitid_policy",
        "orbitid_reconciliation_contract_version",
        "orbitid_correction_signature_sha256",
        *count_columns,
    }
    missing = sorted(required - set(frame))
    if missing:
        raise ValueError(
            f"real-BLS shards lack orbit-ID reconciliation columns: {missing}"
        )
    if set(frame["orbitid_policy"].astype(str)) != {orbitid_policy}:
        raise ValueError("real-BLS row orbit-ID policy disagrees with summary")
    if set(frame["orbitid_reconciliation_contract_version"].astype(str)) != {
        _ORBITID_RECONCILIATION_CONTRACT_VERSION
    }:
        raise ValueError("real-BLS rows have the wrong orbit-ID contract")
    if (
        frame.groupby("tic", sort=False)[
            [*count_columns, "orbitid_correction_signature_sha256"]
        ]
        .nunique(dropna=False)
        .gt(1)
        .any()
        .any()
    ):
        raise ValueError(
            "real-BLS orbit-ID reconciliation fields disagree across target rows"
        )
    targets = frame.drop_duplicates("tic", keep="first").copy()
    for column in count_columns:
        values = pd.to_numeric(targets[column], errors="coerce")
        if (
            values.isna().any()
            or (values < 0).any()
            or (values != np.floor(values)).any()
        ):
            raise ValueError(f"real-BLS shards have invalid {column}")
        targets[column] = values.astype(np.int64)
    if (
        targets["n_cad_orbitid_mismatch"] > targets["n_cad_orbitid_reference_matched"]
    ).any():
        raise ValueError("real-BLS orbit-ID mismatch count exceeds matched count")
    if orbitid_policy == "strict":
        if (
            targets[["n_cad_orbitid_mismatch", "n_cad_orbitid_corrected"]].to_numpy(
                dtype=np.int64
            )
            != 0
        ).any():
            raise ValueError("strict orbit-ID shards contain corrections")
    elif (
        not targets["n_cad_orbitid_corrected"]
        .eq(targets["n_cad_orbitid_mismatch"])
        .all()
    ):
        raise ValueError("reference_by_cadence shards did not correct every mismatch")
    signatures = targets["orbitid_correction_signature_sha256"].astype(str)
    if not signatures.str.fullmatch(r"[0-9a-f]{64}").all():
        raise ValueError("real-BLS shards contain invalid correction signatures")
    records = [
        {
            "tic": int(row.tic),
            "n_cad_orbitid_reference_matched": int(row.n_cad_orbitid_reference_matched),
            "n_cad_orbitid_mismatch": int(row.n_cad_orbitid_mismatch),
            "n_cad_orbitid_corrected": int(row.n_cad_orbitid_corrected),
            "orbitid_correction_signature_sha256": str(
                row.orbitid_correction_signature_sha256
            ),
        }
        for row in targets.sort_values("tic", kind="stable").itertuples(index=False)
    ]
    return {
        "n_cad_orbitid_reference_matched": int(
            targets["n_cad_orbitid_reference_matched"].sum()
        ),
        "n_cad_orbitid_mismatch": int(targets["n_cad_orbitid_mismatch"].sum()),
        "n_cad_orbitid_corrected": int(targets["n_cad_orbitid_corrected"].sum()),
        "n_targets_orbitid_mismatch": int(
            targets["n_cad_orbitid_mismatch"].gt(0).sum()
        ),
        "orbitid_corrections_sha256": _canonical_json_sha256(
            {
                "contract_version": (_ORBITID_RECONCILIATION_CONTRACT_VERSION),
                "policy": orbitid_policy,
                "targets": records,
            }
        ),
    }


def merge_shards(
    *,
    shard_dir: Path,
    out_path: Path,
    n_shards: int,
) -> dict[str, object]:
    if int(n_shards) < 1:
        raise ValueError("n_shards must be positive")
    shard_dir = Path(shard_dir)
    out_path = Path(out_path)
    paths = [
        shard_dir / f"real_adp_bls_peaks_{index:03d}.parquet"
        for index in range(n_shards)
    ]
    summaries = [shard_dir / f"summary_{index:03d}.json" for index in range(n_shards)]
    missing = [str(path) for path in [*paths, *summaries] if not path.exists()]
    if missing:
        raise FileNotFoundError(
            f"missing {len(missing)} real-BLS shard products; first={missing[:5]}"
        )
    initial_hashes = {
        str(path.resolve()): _sha256(path) for path in [*paths, *summaries]
    }
    summary_payloads = [json.loads(path.read_text()) for path in summaries]
    baseline = summary_payloads[0]
    invariant_fields = list(_INVARIANT_SUMMARY_FIELDS)
    optional_groups = (
        _TARGET_SELECTION_INVARIANT_FIELDS,
        _ORBITID_INVARIANT_FIELDS,
    )
    for group in optional_groups:
        present = [field in baseline for field in group]
        if any(present) and not all(present):
            missing_group = [
                field for field, is_present in zip(group, present) if not is_present
            ]
            raise ValueError(
                "real-BLS shard summary has incomplete provenance group: "
                f"{missing_group}"
            )
        if all(present):
            invariant_fields.extend(group)
    for field in invariant_fields:
        if field not in baseline:
            raise ValueError(f"real-BLS shard summary is missing {field}")
    frames: list[pd.DataFrame] = []
    target_sets: list[set[int]] = []
    expected_columns: tuple[str, ...] | None = None
    for index, (path, summary_path, payload) in enumerate(
        zip(paths, summaries, summary_payloads, strict=True)
    ):
        if int(payload.get("shard_index", -1)) != index:
            raise ValueError(f"real-BLS shard {index} has the wrong shard_index")
        if int(payload.get("n_shards", -1)) != int(n_shards):
            raise ValueError(f"real-BLS shard {index} has the wrong n_shards")
        for group in optional_groups:
            present = [field in payload for field in group]
            baseline_present = [field in baseline for field in group]
            if present != baseline_present:
                raise ValueError(
                    f"real-BLS shard {index} has inconsistent provenance fields"
                )
        for field in invariant_fields:
            if payload.get(field) != baseline[field]:
                raise ValueError(
                    f"real-BLS shard {index} disagrees on provenance field {field}"
                )
        if payload.get("peak_table_sha256") != initial_hashes[str(path.resolve())]:
            raise ValueError(f"real-BLS shard {index} peak-table SHA256 mismatch")
        outputs = payload.get("outputs", {})
        if (
            not isinstance(outputs, dict)
            or Path(str(outputs.get("peak_table", ""))).resolve() != path.resolve()
        ):
            raise ValueError(f"real-BLS shard {index} output path mismatch")
        frame = pd.read_parquet(path)
        columns = tuple(str(value) for value in frame.columns)
        if expected_columns is None:
            expected_columns = columns
        elif columns != expected_columns:
            raise ValueError("real-BLS shards have different table schemas")
        required_columns = {"tic", "aperture", "peak_rank", "status"}
        missing_columns = sorted(required_columns - set(frame.columns))
        if missing_columns:
            raise ValueError(
                f"real-BLS shard {index} lacks required columns: {missing_columns}"
            )
        if int(payload.get("n_rows", -1)) != len(frame):
            raise ValueError(f"real-BLS shard {index} row-count mismatch")
        if all(field in baseline for field in _ORBITID_INVARIANT_FIELDS):
            orbitid_summary = _orbitid_summary_from_rows(
                frame,
                orbitid_policy=str(baseline["orbitid_policy"]),
            )
            for field, observed in orbitid_summary.items():
                if payload.get(field) != observed:
                    raise ValueError(
                        f"real-BLS shard {index} orbit-ID summary mismatch for {field}"
                    )
        tic = pd.to_numeric(frame.get("tic"), errors="coerce")
        if tic.isna().any() or (tic <= 0).any() or (tic != np.floor(tic)).any():
            raise ValueError(f"real-BLS shard {index} has invalid TIC values")
        targets = set(tic.astype(np.int64).tolist())
        if len(targets) != int(payload.get("n_targets", -1)):
            raise ValueError(f"real-BLS shard {index} target-count mismatch")
        if any(targets & previous for previous in target_sets):
            raise ValueError("real-BLS shard target sets overlap")
        target_sets.append(targets)
        frames.append(frame)

    frame = pd.concat(frames, ignore_index=True)
    apertures = set(frame["aperture"].dropna().astype(str))
    if not apertures.issubset(set(ADP_ONLY_APERTURES)):
        raise ValueError(
            f"non-ADP apertures found in real-BLS shards: {sorted(apertures)}"
        )
    if sorted(apertures) != sorted(str(value) for value in baseline["apertures"]):
        raise ValueError("merged real-BLS apertures disagree with shard summaries")
    key = ["tic", "aperture", "peak_rank"]
    valid_key = frame[key].notna().all(axis=1)
    if frame.loc[valid_key, key].duplicated().any():
        raise ValueError(
            "real-BLS shards contain duplicate TIC/aperture/peak-rank rows"
        )
    expected_targets = int(baseline["n_targets_total"])
    observed_targets = len(set().union(*target_sets))
    if observed_targets != expected_targets:
        raise ValueError(
            f"real-BLS target coverage is incomplete: {observed_targets} != {expected_targets}"
        )
    observed_tics = set().union(*target_sets)
    if all(field in baseline for field in _TARGET_SELECTION_INVARIANT_FIELDS):
        if (
            str(baseline["target_selection_contract_version"])
            != _TARGET_SELECTION_CONTRACT_VERSION
        ):
            raise ValueError("real-BLS target-selection contract mismatch")
        allowlist_count = int(baseline["target_allowlist_count"])
        if allowlist_count != expected_targets:
            raise ValueError(
                "real-BLS target allowlist count disagrees with n_targets_total"
            )
        observed_inventory_sha256 = _tic_inventory_sha256(observed_tics)
        if str(baseline["target_allowlist_tics_sha256"]) != observed_inventory_sha256:
            raise ValueError(
                "real-BLS merged target inventory disagrees with allowlist"
            )
        allowlist_path = baseline["target_allowlist"]
        allowlist_sha256 = baseline["target_allowlist_sha256"]
        if allowlist_path is None:
            if allowlist_sha256 is not None:
                raise ValueError("real-BLS target allowlist hash exists without a path")
        else:
            allowlist = Path(str(allowlist_path))
            if _sha256(allowlist) != str(allowlist_sha256):
                raise ValueError("real-BLS target allowlist SHA256 mismatch")
            allowlist_tics = _read_target_allowlist(allowlist)
            if len(allowlist_tics) != allowlist_count:
                raise ValueError("real-BLS target allowlist row-count mismatch")
            if _tic_inventory_sha256(allowlist_tics) != str(
                baseline["target_allowlist_tics_sha256"]
            ):
                raise ValueError("real-BLS target allowlist content-signature mismatch")
            if set(allowlist_tics) != observed_tics:
                raise ValueError(
                    "real-BLS merged targets do not exactly cover the allowlist"
                )
    frame = frame.sort_values(key, kind="stable").reset_index(drop=True)
    status = frame["status"].fillna("").astype(str)
    valid = status.eq("ok") & pd.to_numeric(frame.get("peak_rank"), errors="coerce").gt(
        0
    )
    quality_count_columns = (
        "n_cad_internal_bad",
        "n_cad_external_bad",
        "n_cad_external_only_bad",
        "n_cad_effective_bad",
        "n_cad_authority_excluded",
    )
    missing_quality_columns = sorted(set(quality_count_columns) - set(frame))
    if missing_quality_columns:
        raise ValueError(
            f"real-BLS shards lack quality-count columns: {missing_quality_columns}"
        )
    if "n_cad_total" not in frame:
        raise ValueError("real-BLS shards lack n_cad_total")
    quality_numeric = frame.loc[
        :, ["tic", "n_cad_total", *quality_count_columns]
    ].copy()
    for column in ("n_cad_total", *quality_count_columns):
        values = pd.to_numeric(quality_numeric[column], errors="coerce")
        if (
            values.isna().any()
            or (values < 0).any()
            or (values != np.floor(values)).any()
        ):
            raise ValueError(f"real-BLS shards have invalid {column}")
        quality_numeric[column] = values.astype(np.int64)
    if (
        quality_numeric.groupby("tic", sort=False)[
            ["n_cad_total", *quality_count_columns]
        ]
        .nunique()
        .gt(1)
        .any()
        .any()
    ):
        raise ValueError("real-BLS quality counts disagree across target rows")
    if (
        quality_numeric["n_cad_external_only_bad"]
        > quality_numeric["n_cad_external_bad"]
    ).any() or (
        quality_numeric["n_cad_authority_excluded"]
        > quality_numeric["n_cad_external_bad"]
    ).any():
        raise ValueError("real-BLS external quality counts are inconsistent")
    if (
        not quality_numeric["n_cad_effective_bad"]
        .eq(
            quality_numeric["n_cad_internal_bad"]
            + quality_numeric["n_cad_external_only_bad"]
        )
        .all()
    ):
        raise ValueError("real-BLS effective quality counts are inconsistent")
    if (
        quality_numeric.loc[:, quality_count_columns]
        .gt(quality_numeric["n_cad_total"], axis=0)
        .any()
        .any()
    ):
        raise ValueError("real-BLS quality counts exceed n_cad_total")
    target_quality = frame.drop_duplicates("tic", keep="first")
    quality_counts_over_unique_targets = {
        column: int(pd.to_numeric(target_quality[column], errors="raise").sum())
        for column in quality_count_columns
    }
    merged_orbitid_summary: dict[str, Any] = {}
    if all(field in baseline for field in _ORBITID_INVARIANT_FIELDS):
        merged_orbitid_summary = _orbitid_summary_from_rows(
            frame,
            orbitid_policy=str(baseline["orbitid_policy"]),
        )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = out_path.with_suffix(out_path.suffix + ".tmp.parquet")
    temporary.unlink(missing_ok=True)
    frame.to_parquet(temporary, compression="zstd", index=False)
    final_hashes = {str(path.resolve()): _sha256(path) for path in [*paths, *summaries]}
    if final_hashes != initial_hashes:
        temporary.unlink(missing_ok=True)
        raise RuntimeError("real-BLS shard inputs changed during merge")
    temporary.replace(out_path)
    summary: dict[str, Any] = {field: baseline[field] for field in invariant_fields}
    summary.update(
        {
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "out_dir": str(out_path.parent),
            "n_shards": 1,
            "shard_index": 0,
            "n_source_shards": int(n_shards),
            "n_targets": observed_targets,
            "n_targets_total": expected_targets,
            "n_rows": int(len(frame)),
            "n_unique_tics": observed_targets,
            "n_valid_peak_rows": int(valid.sum()),
            "workers": 0,
            "peak_table_sha256": _sha256(out_path),
            "status_counts": {
                str(key): int(value)
                for key, value in status.value_counts().sort_index().items()
            },
            "aperture_counts": {
                str(key): int(value)
                for key, value in frame["aperture"]
                .fillna("")
                .astype(str)
                .value_counts()
                .sort_index()
                .items()
            },
            "quality_counts_over_unique_targets": (quality_counts_over_unique_targets),
            **merged_orbitid_summary,
            "source_shards": [
                {
                    "shard_index": index,
                    "peak_table": str(path),
                    "peak_table_sha256": initial_hashes[str(path.resolve())],
                    "summary": str(summary_path),
                    "summary_sha256": initial_hashes[str(summary_path.resolve())],
                    "n_targets": int(payload["n_targets"]),
                    "n_rows": int(payload["n_rows"]),
                }
                for index, (path, summary_path, payload) in enumerate(
                    zip(paths, summaries, summary_payloads, strict=True)
                )
            ],
            "outputs": {
                "peak_table": str(out_path),
                "summary": str(out_path.with_suffix(".summary.json")),
            },
            "passed": True,
        }
    )
    summary_path = out_path.with_suffix(".summary.json")
    summary_tmp = summary_path.with_suffix(summary_path.suffix + ".tmp")
    summary_tmp.write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    summary_tmp.replace(summary_path)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shard-dir", type=Path, required=True)
    parser.add_argument("--out-path", type=Path, required=True)
    parser.add_argument("--n-shards", type=int, required=True)
    args = parser.parse_args()
    summary = merge_shards(
        shard_dir=args.shard_dir,
        out_path=args.out_path,
        n_shards=args.n_shards,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
