#!/usr/bin/env python3
"""Audit one native-v2 shard through the exact full-pool model transform."""
from __future__ import annotations

import argparse
from collections.abc import Mapping
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from twirl.vetting.harmonic_dataset import (  # noqa: E402
    HarmonicNativeDataset,
    collate_native_samples,
)
from twirl.vetting.harmonic_inputs import (  # noqa: E402
    MODEL_INPUT_CONTRACT_VERSION,
)
from twirl.vetting.ssl_full_pool_eligibility import (  # noqa: E402
    PRODUCTION_ELIGIBLE_IDENTITY_SHA256,
    PRODUCTION_ELIGIBLE_OBSERVATIONS,
    PRODUCTION_FULL_OBSERVATIONS,
    observation_identity_sha256,
)
from twirl.vetting.ssl_full_pool_native import (  # noqa: E402
    FULL_POOL_NATIVE_CONTRACT_VERSION,
    load_full_pool_native_registry_release,
)
from twirl.vetting.ssl_full_pool_numeric import (  # noqa: E402
    FULL_POOL_NUMERIC_ENVELOPE_V1,
    MODEL_INPUT_NUMERIC_AUDIT_SCHEMA,
    MODEL_INPUT_NUMERIC_ENVELOPE_SHA256,
    audit_collated_sample,
)
from twirl.vetting.teacher_native_registry import file_sha256  # noqa: E402
from twirl.vetting.teacher_ssl_fullpool import (  # noqa: E402
    FULLPOOL_SSL_PROFILE,
    fullpool_dataset_rows,
    load_fullpool_ssl_registry,
)


AUDIT_BATCH_SIZE = 8
AUTHORITY_NAMES: tuple[str, ...] = (
    "native_registry",
    "native_registry_summary",
    "native_release_summary",
    "ssl_registry",
    "ssl_registry_summary",
)


def _artifact_metadata(path: Path) -> dict[str, Any]:
    resolved = Path(path).expanduser().resolve(strict=True)
    before = resolved.stat()
    if not resolved.is_file() or before.st_size <= 0:
        raise ValueError(f"authority is missing or empty: {resolved}")
    digest = file_sha256(resolved)
    after = resolved.stat()
    identity_before = (
        before.st_size,
        before.st_mtime_ns,
        before.st_dev,
        before.st_ino,
    )
    identity_after = (
        after.st_size,
        after.st_mtime_ns,
        after.st_dev,
        after.st_ino,
    )
    if identity_before != identity_after:
        raise RuntimeError(f"authority changed while hashing: {resolved}")
    return {
        "path": str(resolved),
        "size_bytes": int(after.st_size),
        "sha256": digest,
    }


def _code_revision() -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    revision = completed.stdout.strip()
    if (
        completed.returncode != 0
        or len(revision) != 40
        or any(character not in "0123456789abcdef" for character in revision)
    ):
        raise RuntimeError("cannot bind numeric audit to a Git revision")
    return revision


def _json_bytes(value: Mapping[str, Any]) -> bytes:
    return (
        json.dumps(
            value,
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


def _publish_immutable(path: Path, payload: bytes) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != payload:
            raise FileExistsError(
                f"refusing to replace different numeric-audit bytes: {path}"
            )
        return
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.link(temporary, path)
        except FileExistsError:
            if path.read_bytes() != payload:
                raise
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _overall_max_abs(maxima: Any) -> float | None:
    if not isinstance(maxima, Mapping):
        return None
    values: list[float] = []
    for item in maxima.values():
        if isinstance(item, Mapping):
            value = item.get("max_abs")
            if isinstance(value, (int, float)) and np.isfinite(value):
                values.append(float(value))
    return max(values) if values else None


def _resolved_path_series(values: pd.Series) -> pd.Series:
    """Resolve each distinct path once rather than once per registry row."""

    text = values.fillna("").astype(str)
    resolved = {
        value: str(Path(value).expanduser().resolve())
        for value in text.unique()
    }
    return text.map(resolved)


def _merge_channel_statistics(
    aggregate: dict[str, dict[str, dict[str, Any]]],
    report: Mapping[str, Any],
) -> None:
    maxima = report.get("maxima")
    if not isinstance(maxima, Mapping):
        return
    for tensor_name in (
        "harmonic_values",
        "local_values",
        "periodogram_values",
    ):
        tensor = maxima.get(tensor_name)
        if not isinstance(tensor, Mapping):
            continue
        tensor_out = aggregate.setdefault(tensor_name, {})
        for channel_name, stats in tensor.items():
            if not isinstance(stats, Mapping):
                continue
            out = tensor_out.setdefault(
                str(channel_name),
                {
                    "n_model_values": 0,
                    "n_masked_values": 0,
                    "min": None,
                    "max": None,
                    "max_abs": None,
                },
            )
            out["n_model_values"] += int(stats.get("n_model_values", 0))
            out["n_masked_values"] += int(stats.get("n_masked_values", 0))
            for name, reducer in (
                ("min", min),
                ("max", max),
                ("max_abs", max),
            ):
                value = stats.get(name)
                if not isinstance(value, (int, float)) or not np.isfinite(value):
                    continue
                previous = out[name]
                out[name] = (
                    float(value)
                    if previous is None
                    else float(reducer(float(previous), float(value)))
                )


def audit_numeric_shard(
    *,
    native_registry_path: Path,
    native_registry_summary_path: Path,
    native_release_summary_path: Path,
    registry_path: Path,
    registry_summary_path: Path,
    sector: int,
    shard_index: int,
    report_out: Path,
    expected_code_revision: str,
) -> dict[str, Any]:
    """Audit one exact `(sector, shard)` and publish an immutable report."""

    sector = int(sector)
    shard_index = int(shard_index)
    if sector not in range(56, 63):
        raise ValueError("sector must be in [56, 62]")
    if shard_index not in range(16):
        raise ValueError("shard_index must be in [0, 15]")
    if (
        len(expected_code_revision) != 40
        or any(
            character not in "0123456789abcdef"
            for character in expected_code_revision
        )
    ):
        raise ValueError(
            "expected_code_revision must be a lowercase 40-hex Git SHA"
        )
    if expected_code_revision != _code_revision():
        raise ValueError("numeric audit repository is not the expected revision")

    authority_paths = {
        "native_registry": native_registry_path,
        "native_registry_summary": native_registry_summary_path,
        "native_release_summary": native_release_summary_path,
        "ssl_registry": registry_path,
        "ssl_registry_summary": registry_summary_path,
    }
    authority_bindings = {
        name: _artifact_metadata(path)
        for name, path in authority_paths.items()
    }
    native_registry, native_release = load_full_pool_native_registry_release(
        registry_path=Path(authority_bindings["native_registry"]["path"]),
        registry_summary_path=Path(
            authority_bindings["native_registry_summary"]["path"]
        ),
        release_summary_path=Path(
            authority_bindings["native_release_summary"]["path"]
        ),
        verify_shard_files=False,
    )
    registry, _registry_summary = load_fullpool_ssl_registry(
        registry_path=Path(authority_bindings["ssl_registry"]["path"]),
        summary_path=Path(authority_bindings["ssl_registry_summary"]["path"]),
    )
    include = registry["ssl_pool_include"].astype(bool)
    if (
        len(registry) != PRODUCTION_FULL_OBSERVATIONS
        or int(include.sum()) != PRODUCTION_ELIGIBLE_OBSERVATIONS
        or observation_identity_sha256(registry.loc[include])
        != PRODUCTION_ELIGIBLE_IDENTITY_SHA256
    ):
        raise ValueError("numeric audit registry differs from production")
    included_keys = set(
        zip(
            registry.loc[include, "sector"].astype(int),
            registry.loc[include, "tic"].astype(int),
        )
    )
    native_keys = set(
        zip(
            native_registry["sector"].astype(int),
            native_registry["tic"].astype(int),
        )
    )
    if included_keys != native_keys:
        raise ValueError("numeric audit SSL/native eligible keys differ")
    counts = native_release.get("counts", {})
    if int(counts.get("native_registry_observations", -1)) != (
        PRODUCTION_ELIGIBLE_OBSERVATIONS
    ):
        raise ValueError("numeric audit native release count differs")

    sector_native = native_registry.loc[
        native_registry["sector"].astype(int).eq(sector)
    ].copy()
    sector_native["_resolved_native_h5_path"] = _resolved_path_series(
        sector_native["native_h5_path"]
    )
    candidates = sorted(
        {
            str(path)
            for path in sector_native["_resolved_native_h5_path"]
            if Path(path).name == f"native_{shard_index}.h5"
        }
    )
    if len(candidates) != 1:
        raise ValueError(
            f"expected one native path for S{sector} shard {shard_index}; "
            f"found={candidates}"
    )
    native_path = Path(candidates[0])
    registry["_resolved_native_h5_path"] = _resolved_path_series(
        registry["native_h5_path"]
    )
    selected = registry.loc[
        include
        & registry["sector"].astype(int).eq(sector)
        & registry["_resolved_native_h5_path"].eq(str(native_path))
    ].sort_values(["sector", "tic"], kind="stable")
    native_selected = sector_native.loc[
        sector_native["_resolved_native_h5_path"].eq(str(native_path))
    ].sort_values(["sector", "tic"], kind="stable")
    if (
        selected.empty
        or set(zip(selected["sector"], selected["tic"]))
        != set(zip(native_selected["sector"], native_selected["tic"]))
    ):
        raise ValueError("numeric audit shard selection differs from native registry")
    declared_sha = set(selected["native_h5_sha256"].astype(str))
    if len(declared_sha) != 1 or file_sha256(native_path) != next(iter(declared_sha)):
        raise ValueError("numeric audit native HDF5 hash differs")
    native_binding = _artifact_metadata(native_path)

    rows = fullpool_dataset_rows(selected)
    dataset = HarmonicNativeDataset(
        rows,
        native_h5=None,
        metadata=np.empty((len(rows), 0), dtype=np.float32),
        cache_size=0,
        profile=FULLPOOL_SSL_PROFILE,
    )
    row_reports: list[dict[str, Any]] = []
    channel_statistics: dict[str, dict[str, dict[str, Any]]] = {}
    pending_samples: list[dict[str, Any]] = []
    pending_rows: list[Any] = []

    def flush() -> None:
        if not pending_samples:
            return
        batch = collate_native_samples(pending_samples)
        for sample_index, row in enumerate(pending_rows):
            report = audit_collated_sample(batch, sample_index=sample_index)
            if (
                report.get("review_id") != str(row.ssl_observation_id)
                or int(report.get("tic", -1)) != int(row.tic)
            ):
                raise RuntimeError(
                    "numeric audit collator changed the observation identity"
                )
            _merge_channel_statistics(channel_statistics, report)
            failures = report.get("failures", [])
            row_reports.append(
                {
                    "ssl_observation_id": str(row.ssl_observation_id),
                    "sector": int(row.sector),
                    "tic": int(row.tic),
                    "passed": bool(report["passed"]),
                    "n_failures": int(report["n_failures"]),
                    "failure_codes": sorted(
                        {
                            str(item.get("code"))
                            for item in failures
                            if isinstance(item, Mapping)
                        }
                    ),
                    "failures": failures,
                    "harmonic_max_abs": _overall_max_abs(
                        report["maxima"].get("harmonic_values")
                    ),
                    "local_max_abs": _overall_max_abs(
                        report["maxima"].get("local_values")
                    ),
                    "periodogram_max_abs": _overall_max_abs(
                        report["maxima"].get("periodogram_values")
                    ),
                }
            )
        pending_samples.clear()
        pending_rows.clear()

    try:
        for index, row in enumerate(selected.itertuples(index=False)):
            try:
                sample = dataset[index]
            except Exception as exc:
                flush()
                row_reports.append(
                    {
                        "ssl_observation_id": str(row.ssl_observation_id),
                        "sector": int(row.sector),
                        "tic": int(row.tic),
                        "passed": False,
                        "n_failures": 1,
                        "failure_codes": ["transform_exception"],
                        "failures": [
                            {
                                "code": "transform_exception",
                                "type": type(exc).__name__,
                                "message": str(exc),
                            }
                        ],
                        "harmonic_max_abs": None,
                        "local_max_abs": None,
                        "periodogram_max_abs": None,
                    }
                )
                continue
            pending_samples.append(sample)
            pending_rows.append(row)
            if len(pending_samples) == AUDIT_BATCH_SIZE:
                flush()
        flush()
    finally:
        dataset.close()

    if _artifact_metadata(native_path) != native_binding:
        raise RuntimeError("numeric audit native HDF5 changed during the scan")
    for name, binding in authority_bindings.items():
        if _artifact_metadata(Path(binding["path"])) != binding:
            raise RuntimeError(
                f"numeric audit authority {name} changed during the scan"
            )

    report_frame = pd.DataFrame.from_records(row_reports)
    if (
        len(report_frame) != len(selected)
        or report_frame["ssl_observation_id"].duplicated().any()
        or observation_identity_sha256(report_frame)
        != observation_identity_sha256(selected)
    ):
        raise RuntimeError("numeric shard audit did not preserve exact coverage")
    failed = ~report_frame["passed"].astype(bool)
    report = {
        "schema_version": MODEL_INPUT_NUMERIC_AUDIT_SCHEMA,
        "code_revision": expected_code_revision,
        "model_input_contract_version": MODEL_INPUT_CONTRACT_VERSION,
        "envelope_contract": FULL_POOL_NUMERIC_ENVELOPE_V1.contract_version,
        "envelope_canonical_sha256": MODEL_INPUT_NUMERIC_ENVELOPE_SHA256,
        "envelope": FULL_POOL_NUMERIC_ENVELOPE_V1.as_dict(),
        "sector": sector,
        "shard_index": shard_index,
        "n_shards": 16,
        "passed": not bool(failed.any()),
        "counts": {
            "scanned_observations": int(len(report_frame)),
            "passed_observations": int((~failed).sum()),
            "failed_observations": int(failed.sum()),
        },
        "observation_identity_sha256": observation_identity_sha256(report_frame),
        "native_h5": native_binding,
        "native_contract_version": FULL_POOL_NATIVE_CONTRACT_VERSION,
        "authority_bindings": authority_bindings,
        "channel_statistics": channel_statistics,
        "rows": row_reports,
        "action": "audit_only_no_clip_no_exclusion",
    }
    payload = _json_bytes(report)
    _publish_immutable(report_out, payload)
    digest = hashlib.sha256(payload).hexdigest()
    _publish_immutable(
        Path(str(report_out) + ".sha256"),
        f"{digest}  {Path(report_out).name}\n".encode("ascii"),
    )
    if report["passed"] is not True:
        raise RuntimeError(
            f"numeric shard gate failed for S{sector} shard {shard_index}: "
            f"{int(failed.sum())} observations"
        )
    return report


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native-registry", type=Path, required=True)
    parser.add_argument("--native-registry-summary", type=Path, required=True)
    parser.add_argument("--native-release-summary", type=Path, required=True)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--registry-summary", type=Path, required=True)
    parser.add_argument("--sector", type=int, choices=range(56, 63), required=True)
    parser.add_argument("--shard-index", type=int, choices=range(16), required=True)
    parser.add_argument("--report-out", type=Path, required=True)
    parser.add_argument("--expected-code-revision", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    report = audit_numeric_shard(
        native_registry_path=args.native_registry,
        native_registry_summary_path=args.native_registry_summary,
        native_release_summary_path=args.native_release_summary,
        registry_path=args.registry,
        registry_summary_path=args.registry_summary,
        sector=args.sector,
        shard_index=args.shard_index,
        report_out=args.report_out,
        expected_code_revision=args.expected_code_revision,
    )
    print(
        json.dumps(
            {
                "event": "teacher_ssl_fullpool_numeric_shard_complete",
                "sector": report["sector"],
                "shard_index": report["shard_index"],
                "counts": report["counts"],
                "report": str(args.report_out.expanduser().resolve()),
            },
            sort_keys=True,
        ),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
