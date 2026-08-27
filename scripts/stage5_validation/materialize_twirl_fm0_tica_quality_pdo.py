#!/usr/bin/env python3
"""Materialize current-QLP TICA mission-quality flags in a user-owned release."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import inspect
import json
import multiprocessing
import os
import re
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

_SHA40 = re.compile(r"^[0-9a-f]{40}$")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _serialize_rows(rows: list[tuple[int, int]]) -> bytes:
    if not rows:
        raise ValueError("TICA query returned no rows")
    ordered = sorted((int(cadence), int(flag)) for cadence, flag in rows)
    cadences = [cadence for cadence, _flag in ordered]
    if (
        any(cadence <= 0 for cadence in cadences)
        or any(flag < 0 for _cadence, flag in ordered)
        or len(cadences) != len(set(cadences))
    ):
        raise ValueError("TICA query returned invalid or duplicate cadence rows")
    return "".join(f"{cadence}, {flag}\n" for cadence, flag in ordered).encode()


def _query_detector(sector: int, camera: int, ccd: int) -> dict[str, Any]:
    """Run one isolated QLP database query in a spawned worker."""

    from qlp.lctools.bin.hlsp import query_ticaflags

    frame = query_ticaflags(int(sector), int(camera), int(ccd))
    if not {"cadence", "flag"}.issubset(frame.columns):
        raise ValueError("QLP TICA query lacks cadence/flag columns")
    rows = [
        (int(cadence), int(flag))
        for cadence, flag in zip(frame["cadence"], frame["flag"], strict=True)
    ]
    payload = _serialize_rows(rows)
    return {
        "camera": int(camera),
        "ccd": int(ccd),
        "payload": payload,
        "n_rows": len(rows),
        "n_nonzero": sum(flag != 0 for _cadence, flag in rows),
    }


def materialize(
    *, sector: int, output_dir: Path, producer_git_sha: str, workers: int
) -> dict[str, Any]:
    sector = int(sector)
    if sector < 67:
        raise ValueError("TICA mission-quality materialization requires S67+")
    if _SHA40.fullmatch(producer_git_sha) is None:
        raise ValueError("producer_git_sha must be a full lowercase Git SHA")
    if workers not in range(1, 5):
        raise ValueError("workers must be in 1..4 to bound PDO database load")
    final = output_dir.expanduser().resolve()
    partial = final.with_name(final.name + ".partial")
    if final.exists() or partial.exists():
        raise FileExistsError(f"refusing to overwrite TICA quality release: {final}")
    partial.mkdir(parents=True)

    try:
        import qlp.lctools.bin.hlsp as qlp_hlsp

        qlp_version = importlib.metadata.version("qlp")
        qlp_source = Path(inspect.getfile(qlp_hlsp)).resolve()
        tasks = [(sector, camera, ccd) for camera in range(1, 5) for ccd in range(1, 5)]
        results: list[dict[str, Any]] = []
        context = multiprocessing.get_context("spawn")
        with ProcessPoolExecutor(max_workers=workers, mp_context=context) as pool:
            futures = {pool.submit(_query_detector, *task): task for task in tasks}
            for future in as_completed(futures):
                _sector, camera, ccd = futures[future]
                result = future.result()
                name = f"ticaffiflag_s{sector}_cam{camera}_ccd{ccd}.txt"
                path = partial / name
                path.write_bytes(result.pop("payload"))
                result.update(
                    {
                        "path": name,
                        "sha256": _sha256(path),
                    }
                )
                results.append(result)
                print(
                    f"[tica-quality] S{sector} cam{camera}/ccd{ccd}: "
                    f"{result['n_rows']:,} rows",
                    flush=True,
                )
        results.sort(key=lambda row: (int(row["camera"]), int(row["ccd"])))
        if len(results) != 16:
            raise RuntimeError("TICA materializer did not return all 16 detectors")
        summary = {
            "schema_version": "twirl_fm0_tica_quality_materialization_v1",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "producer_git_sha": producer_git_sha,
            "sector": sector,
            "mission_quality_type": "tica",
            "quality_bits": {
                "coarse_pointing": 2,
                "reaction_wheel_desaturation": 5,
                "predicted_straylight": 11,
            },
            "source": "qlp.lctools.bin.hlsp.query_ticaflags",
            "qlp_version": qlp_version,
            "qlp_source_path": str(qlp_source),
            "qlp_source_sha256": _sha256(qlp_source),
            "workers": workers,
            "n_detectors": len(results),
            "n_rows": sum(int(row["n_rows"]) for row in results),
            "n_nonzero": sum(int(row["n_nonzero"]) for row in results),
            "detectors": results,
            "cadence_coverage_verified": False,
            "claim_limit": (
                "immutable TICA flag materialization only; run the separate "
                "mission-quality cadence gate before HDF5 use"
            ),
        }
        (partial / "summary.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n"
        )
        (partial / "READY").write_text(producer_git_sha + "\n")
        os.replace(partial, final)
        return summary
    except Exception:
        shutil.rmtree(partial, ignore_errors=True)
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sector", required=True, type=int)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--producer-git-sha", required=True)
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    summary = materialize(
        sector=args.sector,
        output_dir=args.output_dir,
        producer_git_sha=args.producer_git_sha,
        workers=args.workers,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
