"""Fail-closed mission-quality source admission for later-sector FM inputs.

The scientific contract is uniform even though the upstream mission-quality
provider changes at Sector 67.  Sectors before 67 use detector-level SPOC
``DQUALITY`` files; Sector 67 onward uses the three-bit TICA header-derived
quality files used by current QLP.  Both are joined with the orbit/detector
QLP qflag authority only after exact cadence coverage is proven.

This gate verifies source files and cadence identities.  It does not open an
FM HDF5 product, build a six-view shard, admit a temporal panel, or promote a
deferred ORCD sector to accepted A2v1.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from .registry import FM0ContractError, sha256_file

MISSION_QUALITY_RECEIPT_SCHEMA_VERSION = (
    "twirl_fm0_mission_quality_source_receipt_v2"
)
MISSION_QUALITY_READY_STATE = "FM_MISSION_QUALITY_READY"
MISSION_QUALITY_POLICY = "spoc_before_s67_tica_from_s67_v1"
MISSION_QUALITY_EXCLUSION_POLICY = "qlp_quaternion_authority_exclusion_v1"
TICA_QUALITY_ALLOWED_MASK = (1 << 2) | (1 << 5) | (1 << 11)
TICA_MATERIALIZATION_SCHEMA_VERSION = "twirl_fm0_tica_quality_materialization_v1"
_SHA40 = re.compile(r"^[0-9a-f]{40}$")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")


def mission_quality_type(sector: int) -> str:
    """Return the current QLP mission-quality provider for one sector."""

    sector = int(sector)
    if sector < 56:
        raise FM0ContractError("FM mission-quality sectors must be S56+")
    return "spoc" if sector < 67 else "tica"


def _expected_orbits(sector: int, supplied: Sequence[int]) -> tuple[int, int]:
    values = tuple(int(value) for value in supplied)
    expected = (2 * int(sector) + 7, 2 * int(sector) + 8)
    if values != expected:
        raise FM0ContractError(
            f"S{sector} expected orbits are {expected}, observed {values}"
        )
    return expected


def _read_quat(path: Path) -> tuple[int, ...]:
    if not path.is_file() or path.stat().st_size <= 0:
        raise FM0ContractError(f"missing quaternion authority: {path}")
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            if "cadence" not in tuple(reader.fieldnames or ()):
                raise FM0ContractError(
                    f"quaternion authority lacks cadence column: {path}"
                )
            values = tuple(int(float(row["cadence"])) for row in reader)
    except (OSError, UnicodeError, TypeError, ValueError) as exc:
        raise FM0ContractError(f"invalid quaternion authority: {path}") from exc
    if not values or len(values) != len(set(values)):
        raise FM0ContractError(f"quaternion cadences are empty or duplicated: {path}")
    return values


def _read_flag_file(path: Path, *, label: str) -> dict[int, int]:
    if not path.is_file() or path.stat().st_size <= 0:
        raise FM0ContractError(f"missing {label}: {path}")
    rows: list[tuple[int, int]] = []
    try:
        for line in path.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            fields = re.split(r"[\s,]+", line.strip())
            if len(fields) != 2:
                raise ValueError("flag row does not contain two fields")
            cadence, value = int(float(fields[0])), int(float(fields[1]))
            if cadence <= 0 or value < 0:
                raise ValueError("flag row is outside integer bounds")
            rows.append((cadence, value))
    except (OSError, UnicodeError, ValueError) as exc:
        raise FM0ContractError(f"invalid {label}: {path}") from exc
    if not rows:
        raise FM0ContractError(f"empty {label}: {path}")
    result = dict(rows)
    if len(result) != len(rows):
        raise FM0ContractError(f"duplicate cadence in {label}: {path}")
    return result


def _assert_exact_coverage(
    *, observed: set[int], expected: set[int], label: str
) -> None:
    missing = sorted(expected - observed)
    extra = sorted(observed - expected)
    if missing or extra:
        raise FM0ContractError(
            f"{label} cadence coverage differs; missing={len(missing)} "
            f"examples={missing[:5]}, extra={len(extra)} examples={extra[:5]}"
        )


def _mission_authority_exclusions(
    *,
    observed: Mapping[int, int],
    expected: set[int],
    label: str,
) -> list[dict[str, int]]:
    missing = sorted(expected - set(observed))
    if missing:
        raise FM0ContractError(
            f"{label} is missing quaternion cadences; missing={len(missing)} "
            f"examples={missing[:5]}"
        )
    return [
        {"cadence": int(cadence), "mission_quality": int(observed[cadence])}
        for cadence in sorted(set(observed) - expected)
    ]


def _canonical_json_sha256(payload: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _load_tica_materialization(
    *, flag_root: Path, sector: int
) -> tuple[dict[tuple[int, int], Mapping[str, Any]], dict[str, Any], list[dict[str, Any]]]:
    summary_path = flag_root / "summary.json"
    ready_path = flag_root / "READY"
    if not summary_path.is_file() or not ready_path.is_file():
        raise FM0ContractError(
            f"S{sector} TICA materialization lacks summary.json or READY"
        )
    try:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise FM0ContractError(
            f"S{sector} TICA materialization summary is invalid"
        ) from exc
    if not isinstance(summary, Mapping):
        raise FM0ContractError(f"S{sector} TICA materialization summary is not an object")
    producer_git_sha = str(summary.get("producer_git_sha", ""))
    qlp_source_sha256 = str(summary.get("qlp_source_sha256", ""))
    if (
        summary.get("schema_version") != TICA_MATERIALIZATION_SCHEMA_VERSION
        or int(summary.get("sector", -1)) != sector
        or summary.get("mission_quality_type") != "tica"
        or summary.get("source") != "qlp.lctools.bin.hlsp.query_ticaflags"
        or int(summary.get("n_detectors", -1)) != 16
        or summary.get("cadence_coverage_verified") is not False
        or _SHA40.fullmatch(producer_git_sha) is None
        or _SHA256.fullmatch(qlp_source_sha256) is None
        or not str(summary.get("qlp_version", "")).strip()
        or not Path(str(summary.get("qlp_source_path", ""))).is_absolute()
    ):
        raise FM0ContractError(
            f"S{sector} TICA materialization summary violates its contract"
        )
    if ready_path.read_text(encoding="utf-8").strip() != producer_git_sha:
        raise FM0ContractError(
            f"S{sector} TICA READY marker disagrees with producer_git_sha"
        )
    detectors = summary.get("detectors")
    if not isinstance(detectors, list) or len(detectors) != 16:
        raise FM0ContractError(
            f"S{sector} TICA materialization detector inventory is incomplete"
        )
    by_detector: dict[tuple[int, int], Mapping[str, Any]] = {}
    for entry in detectors:
        if not isinstance(entry, Mapping):
            raise FM0ContractError(
                f"S{sector} TICA materialization detector entry is invalid"
            )
        camera = int(entry.get("camera", -1))
        ccd = int(entry.get("ccd", -1))
        key = (camera, ccd)
        expected_name = f"ticaffiflag_s{sector}_cam{camera}_ccd{ccd}.txt"
        if (
            camera not in range(1, 5)
            or ccd not in range(1, 5)
            or key in by_detector
            or entry.get("path") != expected_name
            or _SHA256.fullmatch(str(entry.get("sha256", ""))) is None
            or int(entry.get("n_rows", 0)) <= 0
        ):
            raise FM0ContractError(
                f"S{sector} TICA materialization detector inventory is invalid"
            )
        by_detector[key] = entry
    if set(by_detector) != {
        (camera, ccd) for camera in range(1, 5) for ccd in range(1, 5)
    }:
        raise FM0ContractError(
            f"S{sector} TICA materialization detector inventory is incomplete"
        )
    provenance = {
        "producer_git_sha": producer_git_sha,
        "qlp_version": str(summary["qlp_version"]),
        "qlp_source_path": str(summary["qlp_source_path"]),
        "qlp_source_sha256": qlp_source_sha256,
        "summary_path": str(summary_path),
        "summary_sha256": sha256_file(summary_path),
        "ready_path": str(ready_path),
        "ready_sha256": sha256_file(ready_path),
    }
    bindings = [
        {
            "role": "tica_materialization_summary",
            "path": str(summary_path),
            "sha256": provenance["summary_sha256"],
        },
        {
            "role": "tica_materialization_ready",
            "path": str(ready_path),
            "sha256": provenance["ready_sha256"],
        },
    ]
    return by_detector, provenance, bindings


def verify_mission_quality_sources(
    *,
    sector: int,
    expected_orbits: Sequence[int],
    qlp_root: str | Path,
    mission_flag_root: str | Path | None = None,
) -> dict[str, Any]:
    """Verify all QLP and mission-quality source files for one sector."""

    sector = int(sector)
    orbits = _expected_orbits(sector, expected_orbits)
    provider = mission_quality_type(sector)
    root = Path(qlp_root).expanduser().resolve()
    if not root.is_dir():
        raise FM0ContractError(f"QLP authority root is missing: {root}")
    default_flag_dir = root / ("spocflags" if provider == "spoc" else "ticaflags")
    flag_root = Path(mission_flag_root or default_flag_dir).expanduser().resolve()
    if not flag_root.is_dir():
        raise FM0ContractError(
            f"{provider.upper()} mission-quality root is missing: {flag_root}"
        )

    tica_inventory: dict[tuple[int, int], Mapping[str, Any]] = {}
    tica_materialization: dict[str, Any] | None = None
    materialization_bindings: list[dict[str, Any]] = []
    if provider == "tica":
        (
            tica_inventory,
            tica_materialization,
            materialization_bindings,
        ) = _load_tica_materialization(flag_root=flag_root, sector=sector)

    bindings: list[dict[str, Any]] = list(materialization_bindings)
    n_quat_rows = 0
    n_qflag_rows = 0
    n_mission_rows = 0
    n_mission_rows_retained = 0
    n_qflag_nonzero = 0
    n_mission_nonzero = 0
    n_mission_nonzero_retained = 0
    exclusions_by_detector: dict[str, dict[str, Any]] = {}
    for camera in range(1, 5):
        cadence_by_orbit: dict[int, set[int]] = {}
        for orbit in orbits:
            path = root / f"orbit-{orbit}/ffi/run/cam{camera}_quat.txt"
            cadences = _read_quat(path)
            cadence_by_orbit[orbit] = set(cadences)
            n_quat_rows += len(cadences)
            bindings.append(
                {
                    "role": "qlp_quaternion",
                    "orbit": orbit,
                    "camera": camera,
                    "path": str(path),
                    "sha256": sha256_file(path),
                    "n_rows": len(cadences),
                }
            )

        detector_cadences = set().union(*cadence_by_orbit.values())
        if sum(len(value) for value in cadence_by_orbit.values()) != len(
            detector_cadences
        ):
            raise FM0ContractError(
                f"S{sector} camera {camera} quaternion orbits overlap"
            )
        for ccd in range(1, 5):
            for orbit in orbits:
                path = root / f"orbit-{orbit}/ffi/run/cam{camera}ccd{ccd}_qflag.txt"
                flags = _read_flag_file(path, label="QLP qflag authority")
                _assert_exact_coverage(
                    observed=set(flags),
                    expected=cadence_by_orbit[orbit],
                    label=f"S{sector} orbit {orbit} cam{camera}/ccd{ccd} QLP qflag",
                )
                if not set(flags.values()).issubset({0, 1}):
                    raise FM0ContractError(
                        f"S{sector} orbit {orbit} cam{camera}/ccd{ccd} "
                        "QLP qflag is not binary"
                    )
                n_qflag_rows += len(flags)
                n_qflag_nonzero += sum(value != 0 for value in flags.values())
                bindings.append(
                    {
                        "role": "qlp_detector_qflag",
                        "orbit": orbit,
                        "camera": camera,
                        "ccd": ccd,
                        "path": str(path),
                        "sha256": sha256_file(path),
                        "n_rows": len(flags),
                    }
                )

            prefix = "spocffiflag" if provider == "spoc" else "ticaffiflag"
            path = flag_root / f"{prefix}_s{sector}_cam{camera}_ccd{ccd}.txt"
            flags = _read_flag_file(
                path, label=f"{provider.upper()} mission-quality authority"
            )
            observed_hash = sha256_file(path)
            if provider == "tica":
                declared = tica_inventory[(camera, ccd)]
                if (
                    observed_hash != str(declared["sha256"])
                    or len(flags) != int(declared["n_rows"])
                ):
                    raise FM0ContractError(
                        f"S{sector} cam{camera}/ccd{ccd} TICA file disagrees "
                        "with its materialization summary"
                    )
            if provider == "tica" and any(
                int(value) & ~TICA_QUALITY_ALLOWED_MASK for value in flags.values()
            ):
                raise FM0ContractError(
                    f"S{sector} cam{camera}/ccd{ccd} TICA quality contains "
                    "bits outside the current QLP three-bit authority"
                )
            exclusions = _mission_authority_exclusions(
                observed=flags,
                expected=detector_cadences,
                label=f"S{sector} cam{camera}/ccd{ccd} {provider.upper()} quality",
            )
            n_mission_rows += len(flags)
            n_mission_rows_retained += len(detector_cadences)
            n_mission_nonzero += sum(value != 0 for value in flags.values())
            n_mission_nonzero_retained += sum(
                flags[cadence] != 0 for cadence in detector_cadences
            )
            detector_name = f"cam{camera}_ccd{ccd}"
            exclusions_by_detector[detector_name] = {
                "n_rows": len(exclusions),
                "rows": exclusions,
            }
            bindings.append(
                {
                    "role": f"{provider}_mission_quality",
                    "camera": camera,
                    "ccd": ccd,
                    "path": str(path),
                    "sha256": observed_hash,
                    "n_rows": len(flags),
                }
            )

    authority_exclusions = {
        "policy": MISSION_QUALITY_EXCLUSION_POLICY,
        "n_rows": sum(
            int(entry["n_rows"]) for entry in exclusions_by_detector.values()
        ),
        "by_detector": exclusions_by_detector,
    }

    return {
        "schema_version": MISSION_QUALITY_RECEIPT_SCHEMA_VERSION,
        "sector": sector,
        "expected_orbits": list(orbits),
        "quality_policy": MISSION_QUALITY_POLICY,
        "mission_quality_type": provider,
        "quality_state": MISSION_QUALITY_READY_STATE,
        "passed": True,
        "n_detectors": 16,
        "n_quaternion_files": 8,
        "n_qflag_files": 32,
        "n_mission_quality_files": 16,
        "n_quaternion_rows": n_quat_rows,
        "n_qflag_rows": n_qflag_rows,
        "n_mission_quality_rows": n_mission_rows,
        "n_mission_quality_rows_retained": n_mission_rows_retained,
        "n_mission_quality_rows_excluded_by_quat": authority_exclusions[
            "n_rows"
        ],
        "n_nonzero_qflag_rows": n_qflag_nonzero,
        "n_nonzero_mission_quality_rows": n_mission_nonzero,
        "n_nonzero_mission_quality_rows_retained": n_mission_nonzero_retained,
        "mission_quality_authority_exclusions": authority_exclusions,
        "mission_quality_authority_exclusions_sha256": _canonical_json_sha256(
            authority_exclusions
        ),
        "tica_materialization": tica_materialization,
        "source_bindings": bindings,
        "hdf5_quality_join_verified": False,
        "six_view_shards_verified": False,
        "panel_admission_authorized": False,
        "claim_limit": (
            "mission/QLP source and cadence coverage only; HDF5 and shard "
            "quality joins remain required"
        ),
    }


def write_mission_quality_receipt(
    receipt: Mapping[str, Any], output: str | Path
) -> Path:
    """Publish one immutable passing mission-quality receipt."""

    if receipt.get("passed") is not True:
        raise FM0ContractError("refusing to publish a non-passing quality receipt")
    target = Path(output).expanduser().resolve()
    if target.exists():
        raise FM0ContractError(
            f"refusing to overwrite mission-quality receipt: {target}"
        )
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return target


__all__ = [
    "MISSION_QUALITY_EXCLUSION_POLICY",
    "MISSION_QUALITY_POLICY",
    "MISSION_QUALITY_READY_STATE",
    "MISSION_QUALITY_RECEIPT_SCHEMA_VERSION",
    "TICA_QUALITY_ALLOWED_MASK",
    "mission_quality_type",
    "verify_mission_quality_sources",
    "write_mission_quality_receipt",
]
