#!/usr/bin/env python3
"""Merge all 112 exact model-input audits into a full-pool release gate."""
from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
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
import pyarrow as pa


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from twirl.vetting.harmonic_inputs import (  # noqa: E402
    MODEL_INPUT_CONTRACT_VERSION,
)
from twirl.vetting.ssl_full_pool_eligibility import (  # noqa: E402
    PRODUCTION_ELIGIBLE_IDENTITY_SHA256,
    PRODUCTION_ELIGIBLE_OBSERVATIONS,
    PRODUCTION_EXCLUDED_IDENTITY_SHA256,
    PRODUCTION_EXCLUDED_OBSERVATIONS,
    PRODUCTION_FULL_IDENTITY_SHA256,
    PRODUCTION_FULL_OBSERVATIONS,
    observation_identity_sha256,
)
from twirl.vetting.ssl_full_pool_native import (  # noqa: E402
    load_full_pool_native_registry_release,
)
from twirl.vetting.ssl_full_pool_numeric import (  # noqa: E402
    FULL_POOL_NUMERIC_ENVELOPE_V1,
    MODEL_INPUT_NUMERIC_AUDIT_SCHEMA,
    MODEL_INPUT_NUMERIC_ENVELOPE_SHA256,
    MODEL_INPUT_NUMERIC_RELEASE_SCHEMA,
    normalize_numeric_audit_frame,
    numeric_audit_parquet_bytes,
    read_numeric_audit_parquet,
    validate_numeric_gate_release,
)
from twirl.vetting.teacher_native_registry import file_sha256  # noqa: E402
from twirl.vetting.teacher_ssl_fullpool import (  # noqa: E402
    load_fullpool_ssl_registry,
)


AUTHORITY_NAMES: tuple[str, ...] = (
    "native_registry",
    "native_registry_summary",
    "native_release_summary",
    "ssl_registry",
    "ssl_registry_summary",
)
AUDIT_COLUMNS: tuple[str, ...] = (
    "ssl_observation_id",
    "sector",
    "tic",
    "ssl_pool_include",
    "numeric_status",
    "model_input_numeric_passed",
    "n_failures",
    "failure_codes",
    "failures_json",
    "harmonic_max_abs",
    "local_max_abs",
    "periodogram_max_abs",
)
AUDIT_DTYPES: dict[str, str] = {
    "ssl_observation_id": "object",
    "sector": "int64",
    "tic": "int64",
    "ssl_pool_include": "bool",
    "numeric_status": "object",
    "model_input_numeric_passed": "boolean",
    "n_failures": "int64",
    "failure_codes": "object",
    "failures_json": "object",
    "harmonic_max_abs": "float64",
    "local_max_abs": "float64",
    "periodogram_max_abs": "float64",
}
SHARD_ROW_KEYS: frozenset[str] = frozenset(
    {
        "ssl_observation_id",
        "sector",
        "tic",
        "passed",
        "n_failures",
        "failure_codes",
        "failures",
        "harmonic_max_abs",
        "local_max_abs",
        "periodogram_max_abs",
    }
)


def _artifact_metadata(path: Path) -> dict[str, Any]:
    resolved = Path(path).expanduser().resolve(strict=True)
    before = resolved.stat()
    if not resolved.is_file() or before.st_size <= 0:
        raise ValueError(f"artifact is missing or empty: {resolved}")
    digest = file_sha256(resolved)
    after = resolved.stat()
    if (
        before.st_size,
        before.st_mtime_ns,
        before.st_dev,
        before.st_ino,
    ) != (
        after.st_size,
        after.st_mtime_ns,
        after.st_dev,
        after.st_ino,
    ):
        raise RuntimeError(f"artifact changed while hashing: {resolved}")
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
        raise RuntimeError("cannot bind numeric release to a Git revision")
    return revision


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid numeric shard report: {path}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"numeric shard report must be an object: {path}")
    return value


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


def _canonical_sha256(value: Mapping[str, Any]) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _publish_immutable(path: Path, payload: bytes) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != payload:
            raise FileExistsError(
                f"refusing to replace different numeric-release bytes: {path}"
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


def _parquet_bytes(frame: pd.DataFrame) -> bytes:
    return numeric_audit_parquet_bytes(frame)


def _validate_numeric_audit_frame(frame: pd.DataFrame) -> pd.DataFrame:
    normalized = normalize_numeric_audit_frame(frame)
    observed_dtypes = {
        name: str(normalized[name].dtype) for name in AUDIT_COLUMNS
    }
    if observed_dtypes != AUDIT_DTYPES:
        raise ValueError("numeric audit table has the wrong exact dtypes")
    if (
        len(normalized) != PRODUCTION_FULL_OBSERVATIONS
        or normalized["ssl_observation_id"].duplicated().any()
        or normalized[["sector", "tic"]].duplicated().any()
        or observation_identity_sha256(normalized)
        != PRODUCTION_FULL_IDENTITY_SHA256
    ):
        raise ValueError("numeric audit table differs from the full pool")
    include = normalized["ssl_pool_include"].astype(bool)
    included = normalized.loc[include]
    excluded = normalized.loc[~include]
    if (
        len(included) != PRODUCTION_ELIGIBLE_OBSERVATIONS
        or observation_identity_sha256(included)
        != PRODUCTION_ELIGIBLE_IDENTITY_SHA256
        or not included["numeric_status"].eq("passed").all()
        or included["model_input_numeric_passed"].isna().any()
        or not included["model_input_numeric_passed"].astype(bool).all()
        or not included["n_failures"].eq(0).all()
        or not included["failure_codes"].eq("[]").all()
        or not included["failures_json"].eq("[]").all()
    ):
        raise ValueError("numeric audit eligible partition is invalid")
    if (
        len(excluded) != PRODUCTION_EXCLUDED_OBSERVATIONS
        or observation_identity_sha256(excluded)
        != PRODUCTION_EXCLUDED_IDENTITY_SHA256
        or not excluded["numeric_status"].eq("not_model_eligible").all()
        or not excluded["model_input_numeric_passed"].isna().all()
        or not excluded["n_failures"].eq(0).all()
        or not excluded["failure_codes"].eq("[]").all()
        or not excluded["failures_json"].eq("[]").all()
    ):
        raise ValueError("numeric audit excluded partition is invalid")
    return normalized


def _validate_passed_shard_rows(
    rows: Any,
    *,
    sector: int,
    expected_count: int,
    expected_identity_sha256: Any,
    context: str,
) -> pd.DataFrame:
    """Return exact passed row evidence without truthy coercions."""

    if not isinstance(rows, list) or len(rows) != int(expected_count):
        raise ValueError(f"{context} row inventory is invalid")
    normalized: list[dict[str, Any]] = []
    for index, row in enumerate(rows):
        row_context = f"{context}.rows[{index}]"
        if not isinstance(row, Mapping) or set(row) != SHARD_ROW_KEYS:
            raise ValueError(f"{row_context} has an invalid record shape")
        observation_id = row.get("ssl_observation_id")
        row_sector = row.get("sector")
        tic = row.get("tic")
        n_failures = row.get("n_failures")
        if (
            not isinstance(observation_id, str)
            or not observation_id
            or type(row_sector) is not int
            or int(row_sector) != int(sector)
            or type(tic) is not int
            or type(n_failures) is not int
            or n_failures != 0
            or row.get("passed") is not True
            or row.get("failure_codes") != []
            or row.get("failures") != []
        ):
            raise ValueError(f"{row_context} is not an exact passed row")
        for name in (
            "harmonic_max_abs",
            "local_max_abs",
            "periodogram_max_abs",
        ):
            value = row.get(name)
            if value is not None and (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not np.isfinite(float(value))
                or float(value) < 0.0
            ):
                raise ValueError(
                    f"{row_context}.{name} is not finite nonnegative evidence"
                )
        normalized.append(dict(row))
    frame = pd.DataFrame.from_records(normalized)
    if (
        frame.empty
        or frame["ssl_observation_id"].duplicated().any()
        or frame[["sector", "tic"]].duplicated().any()
        or observation_identity_sha256(frame) != expected_identity_sha256
    ):
        raise ValueError(f"{context} row identity differs")
    return frame


def _normalized_binding(value: Any, *, context: str) -> dict[str, Any]:
    if not isinstance(value, Mapping) or set(value) != {
        "path",
        "size_bytes",
        "sha256",
    }:
        raise ValueError(f"{context} is not an exact artifact binding")
    return {
        "path": str(Path(str(value["path"])).expanduser().resolve()),
        "size_bytes": int(value["size_bytes"]),
        "sha256": str(value["sha256"]),
    }


def _merge_channel_statistics(
    aggregate: dict[str, dict[str, dict[str, Any]]],
    source: Any,
) -> None:
    if not isinstance(source, Mapping):
        raise ValueError("numeric shard channel statistics are invalid")
    for tensor_name, channels in source.items():
        if not isinstance(channels, Mapping):
            raise ValueError("numeric shard channel statistics are invalid")
        tensor_out = aggregate.setdefault(str(tensor_name), {})
        for channel_name, stats in channels.items():
            if not isinstance(stats, Mapping):
                raise ValueError("numeric shard channel statistics are invalid")
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
                out[name] = (
                    float(value)
                    if out[name] is None
                    else float(reducer(float(out[name]), float(value)))
                )


def merge_numeric_release(
    *,
    native_registry_path: Path,
    native_registry_summary_path: Path,
    native_release_summary_path: Path,
    registry_path: Path,
    registry_summary_path: Path,
    shard_report_paths: Sequence[Path],
    release_out: Path,
    expected_code_revision: str,
    audit_out: Path | None = None,
) -> dict[str, Any]:
    """Validate exact shard union and publish the 175,366-row audit release."""

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
        raise ValueError("numeric merge repository is not the expected revision")
    reports = [Path(path).expanduser().resolve(strict=True) for path in shard_report_paths]
    if len(reports) != 112 or len(set(reports)) != 112:
        raise ValueError("numeric merge requires exactly 112 unique shard reports")
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
    native_registry, _native_release = load_full_pool_native_registry_release(
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
    expected_native_sha256: dict[str, str] = {}
    for path_text, group in native_registry.groupby(
        "native_h5_path",
        sort=False,
    ):
        resolved = str(Path(str(path_text)).expanduser().resolve())
        hashes = set(group["native_h5_sha256"].astype(str))
        if len(hashes) != 1:
            raise ValueError(
                f"native registry has conflicting hashes for {resolved}"
            )
        expected_native_sha256[resolved] = next(iter(hashes))

    expected_inventory = {
        (sector, shard) for sector in range(56, 63) for shard in range(16)
    }
    observed_inventory: set[tuple[int, int]] = set()
    observed_native_paths: set[str] = set()
    row_records: list[dict[str, Any]] = []
    report_bindings: list[dict[str, Any]] = []
    channel_statistics: dict[str, dict[str, dict[str, Any]]] = {}
    for path in sorted(reports):
        metadata = _artifact_metadata(path)
        sidecar = Path(str(path) + ".sha256")
        if not sidecar.is_file():
            raise FileNotFoundError(sidecar)
        expected_sidecar = f"{metadata['sha256']}  {path.name}\n"
        if sidecar.read_text(encoding="ascii") != expected_sidecar:
            raise ValueError(f"numeric shard SHA sidecar differs: {path}")
        report = _read_json(path)
        key = (int(report.get("sector", -1)), int(report.get("shard_index", -1)))
        if key not in expected_inventory or key in observed_inventory:
            raise ValueError(f"numeric shard inventory is invalid: {key}")
        observed_inventory.add(key)
        if (
            report.get("schema_version") != MODEL_INPUT_NUMERIC_AUDIT_SCHEMA
            or report.get("code_revision") != expected_code_revision
            or report.get("model_input_contract_version")
            != MODEL_INPUT_CONTRACT_VERSION
            or report.get("envelope_contract")
            != FULL_POOL_NUMERIC_ENVELOPE_V1.contract_version
            or report.get("envelope_canonical_sha256")
            != MODEL_INPUT_NUMERIC_ENVELOPE_SHA256
            or not isinstance(report.get("envelope"), Mapping)
            or _canonical_sha256(report["envelope"])
            != MODEL_INPUT_NUMERIC_ENVELOPE_SHA256
            or report.get("n_shards") != 16
            or report.get("passed") is not True
            or int(report.get("counts", {}).get("failed_observations", -1)) != 0
        ):
            raise ValueError(f"numeric shard report did not pass: {path}")
        declared_authorities = report.get("authority_bindings")
        if not isinstance(declared_authorities, Mapping) or set(
            declared_authorities
        ) != set(AUTHORITY_NAMES):
            raise ValueError(f"numeric shard authorities are invalid: {path}")
        for name, expected in authority_bindings.items():
            if _normalized_binding(
                declared_authorities[name],
                context=f"{path.name}.{name}",
            ) != expected:
                raise ValueError(f"numeric shard authority {name} differs: {path}")
        native_h5 = _normalized_binding(
            report.get("native_h5"),
            context=f"{path.name}.native_h5",
        )
        if native_h5["path"] in observed_native_paths:
            raise ValueError("numeric shard reports reuse a native HDF5")
        if (
            native_h5["path"] not in expected_native_sha256
            or native_h5["sha256"]
            != expected_native_sha256[native_h5["path"]]
        ):
            raise ValueError(
                f"numeric shard native binding differs from registry: {path}"
            )
        observed_native_paths.add(native_h5["path"])
        scanned = report["counts"].get("scanned_observations")
        passed = report["counts"].get("passed_observations")
        if (
            type(scanned) is not int
            or type(passed) is not int
            or int(passed) != int(scanned)
        ):
            raise ValueError(f"numeric shard pass counts differ: {path}")
        report_rows = _validate_passed_shard_rows(
            report.get("rows"),
            sector=key[0],
            expected_count=int(scanned),
            expected_identity_sha256=report.get(
                "observation_identity_sha256"
            ),
            context=str(path),
        )
        row_records.extend(report_rows.to_dict(orient="records"))
        _merge_channel_statistics(channel_statistics, report.get("channel_statistics"))
        report_bindings.append(
            {
                "sector": key[0],
                "shard_index": key[1],
                **metadata,
                "native_h5": native_h5,
                "observation_identity_sha256": report.get(
                    "observation_identity_sha256"
                ),
                "scanned_observations": len(report_rows),
            }
        )
        if _artifact_metadata(path) != metadata:
            raise RuntimeError(
                f"numeric shard report changed during merge: {path}"
            )

    if observed_inventory != expected_inventory:
        raise ValueError("numeric shard inventory is not exact")
    expected_native_paths = set(expected_native_sha256)
    if observed_native_paths != expected_native_paths:
        raise ValueError("numeric shard HDF5 union differs from native release")

    eligible_audit = pd.DataFrame.from_records(row_records)
    if (
        len(eligible_audit) != PRODUCTION_ELIGIBLE_OBSERVATIONS
        or eligible_audit["ssl_observation_id"].duplicated().any()
        or eligible_audit[["sector", "tic"]].duplicated().any()
        or not eligible_audit["passed"].astype(bool).all()
        or int(eligible_audit["n_failures"].sum()) != 0
        or observation_identity_sha256(eligible_audit)
        != PRODUCTION_ELIGIBLE_IDENTITY_SHA256
    ):
        raise ValueError("merged numeric eligible audit differs from production")
    include = registry["ssl_pool_include"].astype(bool)
    included = registry.loc[
        include, ["ssl_observation_id", "sector", "tic", "ssl_pool_include"]
    ].merge(
        eligible_audit,
        on=["ssl_observation_id", "sector", "tic"],
        how="left",
        validate="one_to_one",
    )
    if included["passed"].isna().any():
        raise ValueError("numeric audit is missing eligible registry rows")
    included["numeric_status"] = "passed"
    included["model_input_numeric_passed"] = pd.array(
        np.ones(len(included), dtype=bool),
        dtype="boolean",
    )
    included["failure_codes"] = included["failure_codes"].map(
        lambda values: json.dumps(values, separators=(",", ":"), sort_keys=True)
    )
    included["failures_json"] = included["failures"].map(
        lambda values: json.dumps(values, separators=(",", ":"), sort_keys=True)
    )
    included = included.drop(columns=["passed", "failures"])

    excluded = registry.loc[
        ~include, ["ssl_observation_id", "sector", "tic", "ssl_pool_include"]
    ].copy()
    excluded["numeric_status"] = "not_model_eligible"
    excluded["model_input_numeric_passed"] = pd.array(
        [pd.NA] * len(excluded),
        dtype="boolean",
    )
    excluded["n_failures"] = 0
    excluded["failure_codes"] = "[]"
    excluded["failures_json"] = "[]"
    for name in (
        "harmonic_max_abs",
        "local_max_abs",
        "periodogram_max_abs",
    ):
        excluded[name] = np.nan
    audit = pd.concat([included, excluded], ignore_index=True)
    audit = audit.loc[:, list(AUDIT_COLUMNS)].sort_values(
        ["sector", "tic"],
        kind="stable",
    ).reset_index(drop=True)
    audit = audit.astype(AUDIT_DTYPES)
    excluded_audit = audit.loc[~audit["ssl_pool_include"].astype(bool)]
    if (
        len(audit) != PRODUCTION_FULL_OBSERVATIONS
        or len(excluded_audit) != PRODUCTION_EXCLUDED_OBSERVATIONS
        or observation_identity_sha256(audit) != PRODUCTION_FULL_IDENTITY_SHA256
        or observation_identity_sha256(excluded_audit)
        != PRODUCTION_EXCLUDED_IDENTITY_SHA256
    ):
        raise ValueError("full numeric audit partition differs from production")
    for name, binding in authority_bindings.items():
        if _artifact_metadata(Path(binding["path"])) != binding:
            raise RuntimeError(
                f"numeric release authority {name} changed during merge"
            )

    release_out = Path(release_out).expanduser().resolve()
    audit_out = (
        Path(audit_out).expanduser().resolve()
        if audit_out is not None
        else release_out.parent / "model_input_numeric_audit.parquet"
    )
    audit = _validate_numeric_audit_frame(audit)
    if audit_out.exists():
        existing_audit = _validate_numeric_audit_frame(
            read_numeric_audit_parquet(audit_out)
        )
        try:
            pd.testing.assert_frame_equal(
                existing_audit,
                audit,
                check_dtype=True,
                check_exact=True,
            )
        except AssertionError as exc:
            raise FileExistsError(
                f"existing numeric audit differs from recomputation: {audit_out}"
            ) from exc
        audit_payload = audit_out.read_bytes()
    else:
        audit_payload = _parquet_bytes(audit)
        _validate_numeric_audit_frame(
            read_numeric_audit_parquet(pa.BufferReader(audit_payload))
        )
        _publish_immutable(audit_out, audit_payload)
    audit_sha256 = hashlib.sha256(audit_payload).hexdigest()
    _publish_immutable(
        Path(str(audit_out) + ".sha256"),
        f"{audit_sha256}  {audit_out.name}\n".encode("ascii"),
    )
    top_offenders = {}
    for column in (
        "harmonic_max_abs",
        "local_max_abs",
        "periodogram_max_abs",
    ):
        top_offenders[column] = (
            eligible_audit.nlargest(10, column)[
                ["ssl_observation_id", "sector", "tic", column]
            ].to_dict(orient="records")
        )
    release = {
        "schema_version": MODEL_INPUT_NUMERIC_RELEASE_SCHEMA,
        "code_revision": expected_code_revision,
        "model_input_contract_version": MODEL_INPUT_CONTRACT_VERSION,
        "envelope_contract": FULL_POOL_NUMERIC_ENVELOPE_V1.contract_version,
        "envelope_canonical_sha256": MODEL_INPUT_NUMERIC_ENVELOPE_SHA256,
        "envelope": FULL_POOL_NUMERIC_ENVELOPE_V1.as_dict(),
        "passed": True,
        "counts": {
            "full_observations": PRODUCTION_FULL_OBSERVATIONS,
            "eligible_observations": PRODUCTION_ELIGIBLE_OBSERVATIONS,
            "excluded_observations": PRODUCTION_EXCLUDED_OBSERVATIONS,
            "scanned_observations": PRODUCTION_ELIGIBLE_OBSERVATIONS,
            "failed_observations": 0,
            "native_shards": 112,
        },
        "identity_hashes": {
            "full": PRODUCTION_FULL_IDENTITY_SHA256,
            "eligible": PRODUCTION_ELIGIBLE_IDENTITY_SHA256,
            "excluded": PRODUCTION_EXCLUDED_IDENTITY_SHA256,
        },
        "authority_bindings": authority_bindings,
        "outputs": {
            "numeric_audit": {
                "path": str(audit_out),
                "size_bytes": len(audit_payload),
                "sha256": audit_sha256,
            }
        },
        "shard_reports": report_bindings,
        "channel_statistics": channel_statistics,
        "top_offenders": top_offenders,
        "quality_bad_photometry_policy_verified": True,
        "float32_conversion_verified": True,
        "real_only": True,
        "labels_consumed": False,
        "injections_consumed": False,
        "action": "audit_only_no_clip_no_exclusion",
    }
    release_payload = _json_bytes(release)
    _publish_immutable(release_out, release_payload)
    release_sha256 = hashlib.sha256(release_payload).hexdigest()
    _publish_immutable(
        Path(str(release_out) + ".sha256"),
        f"{release_sha256}  {release_out.name}\n".encode("ascii"),
    )
    validate_numeric_gate_release(
        release_out,
        expected_code_revision=expected_code_revision,
        expected_authority_bindings=authority_bindings,
    )
    return release


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native-registry", type=Path, required=True)
    parser.add_argument("--native-registry-summary", type=Path, required=True)
    parser.add_argument("--native-release-summary", type=Path, required=True)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--registry-summary", type=Path, required=True)
    parser.add_argument(
        "--shard-reports",
        type=Path,
        nargs="+",
        required=True,
    )
    parser.add_argument("--release-out", type=Path, required=True)
    parser.add_argument("--audit-out", type=Path)
    parser.add_argument("--expected-code-revision", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    release = merge_numeric_release(
        native_registry_path=args.native_registry,
        native_registry_summary_path=args.native_registry_summary,
        native_release_summary_path=args.native_release_summary,
        registry_path=args.registry,
        registry_summary_path=args.registry_summary,
        shard_report_paths=args.shard_reports,
        release_out=args.release_out,
        expected_code_revision=args.expected_code_revision,
        audit_out=args.audit_out,
    )
    print(
        json.dumps(
            {
                "event": "teacher_ssl_fullpool_numeric_release_complete",
                "counts": release["counts"],
                "release": str(args.release_out.expanduser().resolve()),
            },
            sort_keys=True,
        ),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
