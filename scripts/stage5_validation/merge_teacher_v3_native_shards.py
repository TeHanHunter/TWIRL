#!/usr/bin/env python3
"""Merge and strictly verify one sector's Teacher-v3 native-input shards."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import h5py
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.vetting.harmonic_export import (  # noqa: E402
    merge_raw_pair_shards,
    read_candidate_table,
)
from twirl.vetting.harmonic_inputs import (  # noqa: E402
    native_group_path,
    verify_raw_pair_contract,
)
from twirl.vetting.s63_preprocessing import (  # noqa: E402
    validate_producer_git_sha,
)

S63_PROSPECTIVE_CONTRACT = "s63_teacher_v3_prospective_v1"


def _authorize_sector(*, sector: int, prospective_contract: str | None) -> str:
    if 56 <= int(sector) <= 62:
        if prospective_contract is not None:
            raise ValueError(
                "legacy Teacher-v3 sectors must not set --prospective-contract"
            )
        return "teacher_v3_s56_s62_legacy"
    if int(sector) == 63 and prospective_contract == S63_PROSPECTIVE_CONTRACT:
        return S63_PROSPECTIVE_CONTRACT
    if int(sector) == 63:
        raise ValueError(
            "S63 native merge requires --prospective-contract "
            f"{S63_PROSPECTIVE_CONTRACT}"
        )
    raise ValueError(
        "Teacher-v3 native merge is bounded to legacy S56-S62 plus "
        "explicitly authorized prospective S63"
    )


def _truth(values):
    if values.dtype == bool:
        return values.fillna(False)
    return (
        values.fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
        .isin({"1", "true", "t", "yes", "y"})
    )


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _path_set_sha256(paths: set[str]) -> str:
    payload = ("\n".join(sorted(paths)) + "\n").encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _validate_shard_contract(
    shard_paths: list[Path],
    *,
    training_table_sha256: str,
    expected_producer_git_sha: str | None = None,
) -> dict[str, object]:
    """Require one complete, nonduplicated shard-index sequence."""

    import h5py

    paths = [Path(path).resolve() for path in shard_paths]
    if len(set(paths)) != len(paths):
        raise ValueError("Teacher-v3 native merge received duplicate shard paths")
    expected_n_shards = len(paths)
    shard_indices: list[int] = []
    for path in paths:
        with h5py.File(path, "r") as source:
            for name in ("shard_index", "n_shards", "training_table_sha256"):
                if name not in source.attrs:
                    raise ValueError(f"native shard {path} lacks {name}")
            shard_index = int(source.attrs["shard_index"])
            declared_n_shards = int(source.attrs["n_shards"])
            if declared_n_shards != expected_n_shards:
                raise ValueError(
                    f"native shard {path} declares n_shards={declared_n_shards}; "
                    f"received {expected_n_shards} shard files"
                )
            if str(source.attrs["training_table_sha256"]) != (
                training_table_sha256
            ):
                raise ValueError(
                    f"native shard {path} does not bind the exact training table"
                )
            if expected_producer_git_sha is not None:
                observed_producer = validate_producer_git_sha(
                    source.attrs.get("producer_git_sha"),
                    label=f"native shard {path} producer_git_sha",
                )
                if observed_producer != expected_producer_git_sha:
                    raise ValueError(
                        f"native shard {path} producer Git SHA differs from merge"
                    )
            shard_indices.append(shard_index)
    expected_indices = list(range(expected_n_shards))
    if sorted(shard_indices) != expected_indices:
        raise ValueError(
            "native shard indices are incomplete or duplicated: "
            f"observed={sorted(shard_indices)}, expected={expected_indices}"
        )
    return {
        "n_shards": int(expected_n_shards),
        "shard_indices": sorted(shard_indices),
        "training_table_sha256": training_table_sha256,
    }


def _native_group_paths(path: Path) -> set[str]:
    import h5py

    paths: set[str] = set()
    with h5py.File(path, "r") as source:
        for root in ("targets", "injections"):
            if root not in source:
                raise ValueError(f"merged native input lacks /{root}")
            paths.update(f"{root}/{key}" for key in source[root])
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sector", type=int, required=True)
    parser.add_argument("--shards", type=Path, nargs="+", required=True)
    parser.add_argument("--out-h5", type=Path, required=True)
    parser.add_argument("--training-table", type=Path, required=True)
    parser.add_argument(
        "--prospective-contract",
        choices=(S63_PROSPECTIVE_CONTRACT,),
        help="Required only for the sealed prospective S63 inference lane.",
    )
    parser.add_argument("--producer-git-sha")
    args = parser.parse_args()
    try:
        authorization = _authorize_sector(
            sector=args.sector,
            prospective_contract=args.prospective_contract,
        )
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    producer_git_sha: str | None = None
    if args.sector == 63:
        if args.producer_git_sha is None:
            raise SystemExit("prospective S63 native merge requires --producer-git-sha")
        producer_git_sha = validate_producer_git_sha(args.producer_git_sha)
    elif args.producer_git_sha is not None:
        raise SystemExit("legacy Teacher-v3 native merge must not set --producer-git-sha")
    output = args.out_h5.resolve()
    summary_path = output.with_suffix(".summary.json")
    pending = output.with_suffix(output.suffix + ".merge_pending")
    if output.exists() or summary_path.exists() or pending.exists():
        raise FileExistsError(
            "Teacher-v3 native merge requires fresh output, summary, and "
            f"pending paths: {output}"
        )
    training_table = args.training_table.resolve()
    training_table_sha256 = _file_sha256(training_table)
    source_rows = read_candidate_table(args.training_table)
    if _file_sha256(training_table) != training_table_sha256:
        raise RuntimeError(
            "Teacher-v3 sector table changed while it was initially read"
        )
    if "native_input_include" not in source_rows:
        raise KeyError("Teacher-v3 sector table lacks native_input_include")
    rows = source_rows.loc[_truth(source_rows["native_input_include"])].copy()
    sectors = set(
        pd.to_numeric(rows["sector"], errors="raise").astype(int)
    )
    if sectors != {int(args.sector)}:
        raise ValueError(
            f"active native rows have sectors {sorted(sectors)}; "
            f"expected only {int(args.sector)}"
        )
    if "native_group_path" not in rows:
        rows["native_group_path"] = [
            native_group_path(row) for row in rows.to_dict("records")
        ]
    expected_paths = set(rows["native_group_path"].astype(str))
    injection_paths = {
        path for path in expected_paths if path.startswith("injections/")
    }
    target_paths = expected_paths - injection_paths
    expected = {
        "targets": int(len(target_paths)),
        "injections": int(len(injection_paths)),
    }
    shard_contract = _validate_shard_contract(
        args.shards,
        training_table_sha256=training_table_sha256,
        expected_producer_git_sha=producer_git_sha,
    )
    merge = merge_raw_pair_shards(
        shard_paths=args.shards,
        out_h5=pending,
    )
    if producer_git_sha is not None:
        with h5py.File(pending, "r+") as h5:
            h5.attrs["producer_git_sha"] = producer_git_sha
            h5.flush()
    verification = verify_raw_pair_contract(
        pending,
        require_errors=True,
        require_periodograms=True,
    )
    observed_paths = _native_group_paths(pending)
    missing_paths = sorted(expected_paths - observed_paths)
    unexpected_paths = sorted(observed_paths - expected_paths)
    identity_match = not missing_paths and not unexpected_paths
    count_match = merge["counts"] == expected
    summary = {
        "authorization_contract": authorization,
        "sector": int(args.sector),
        "merge": merge,
        "verification": verification,
        "training_table": str(training_table),
        "training_table_sha256": training_table_sha256,
        "shard_contract": shard_contract,
        "expected_counts": expected,
        "exact_count_match": count_match,
        "exact_group_identity_match": identity_match,
        "expected_group_paths_sha256": _path_set_sha256(expected_paths),
        "observed_group_paths_sha256": _path_set_sha256(observed_paths),
        "missing_group_paths": missing_paths,
        "unexpected_group_paths": unexpected_paths,
    }
    if producer_git_sha is not None:
        summary["producer_git_sha"] = producer_git_sha
    if (
        not verification["passed"]
        or not count_match
        or not identity_match
    ):
        pending.unlink(missing_ok=True)
        summary_path.write_text(
            json.dumps(summary, indent=2, sort_keys=True, allow_nan=False)
            + "\n"
        )
        print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))
        raise SystemExit(2)
    if _file_sha256(training_table) != training_table_sha256:
        pending.unlink(missing_ok=True)
        raise RuntimeError(
            "Teacher-v3 sector table changed during native shard merge"
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    pending.replace(output)
    merge["out_h5"] = str(output)
    summary["merge"] = merge
    summary["out_h5_sha256"] = _file_sha256(output)
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
