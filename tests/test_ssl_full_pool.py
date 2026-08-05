from __future__ import annotations

import hashlib
import json
from pathlib import Path

import h5py
import pandas as pd
import pytest

from twirl.vetting.ssl_full_pool import (
    EXPECTED_FLUX_COLUMNS,
    EXPECTED_SECTORS,
    FULL_POOL_CONTRACT_VERSION,
    POOL_COLUMNS,
    write_ssl_full_pool,
)
from twirl.vetting.teacher_split_registry import build_tic_split_registry


def _write_split_registry(path: Path) -> tuple[int, int, int]:
    rows = pd.DataFrame(
        {
            "sector": [56] * 30,
            "tic": list(range(10_000, 10_030)),
            "candidate_key": [
                f"s0056-tic{tic:016d}" for tic in range(10_000, 10_030)
            ],
            "split_stratum": [
                ("planet_like", "eclipse_contact", "other")[index % 3]
                for index in range(30)
            ],
        }
    )
    registry, _ = build_tic_split_registry(rows, seed=560063)
    registry.to_csv(path, index=False, lineterminator="\n")
    test_tics = (
        registry.loc[registry["fixed_split"].eq("test"), "tic"]
        .astype(int)
        .sort_values()
        .tolist()
    )
    development_tic = int(
        registry.loc[
            registry["fixed_split"].eq("development"), "tic"
        ].min()
    )
    return test_tics[0], test_tics[1], development_tic


def _write_compact_export(
    root: Path,
    *,
    sector: int,
    tics: list[int],
) -> Path:
    h5_path = root / f"s{sector}_A2v1_adp_pair.h5"
    manifest_path = root / f"s{sector}_A2v1_adp_pair.manifest.json"
    records: list[dict[str, object]] = []
    with h5py.File(h5_path, "w") as h5:
        h5.attrs["sector"] = sector
        h5.attrs["n_targets"] = len(tics)
        h5.attrs["time_unit"] = "BJD - 2457000"
        h5.attrs["flux_columns"] = json.dumps(list(EXPECTED_FLUX_COLUMNS))
        targets = h5.create_group("targets")
        for index, tic in enumerate(tics):
            targets.create_group(f"{tic:016d}")
            records.append(
                {
                    "tic": tic,
                    "sector": sector,
                    "camera": index % 4 + 1,
                    "ccd": (index + 1) % 4 + 1,
                    "tessmag": 15.0 + index / 10,
                    "n_cadences": 1_000 + index,
                    "flux_columns": list(EXPECTED_FLUX_COLUMNS),
                    "source_fits": (
                        f"/readonly/hlsp/s{sector}/tic{tic:016d}.fits"
                    ),
                }
            )
    manifest = {
        "created_utc": "2026-07-27T00:00:00+00:00",
        "sector": sector,
        "hlsp_root": f"/readonly/hlsp/s{sector}",
        "out_h5": f"/pdo/users/tehan/{h5_path.name}",
        "time_unit": "BJD - 2457000",
        "requested_columns": list(EXPECTED_FLUX_COLUMNS),
        "n_discovered_files": len(records),
        "n_exported_targets": len(records),
        "skipped": {
            "read_failed": 0,
            "tic_filter": 0,
            "duplicate_tic": 0,
            "no_flux_columns": 0,
        },
        "records": records,
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    return manifest_path


def _write_s63_reservation(path: Path, tics: list[int]) -> None:
    normalized = sorted(tics)
    path.write_text("".join(f"{tic}\n" for tic in normalized))
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    summary = {
        "schema_version": "twirl_teacher_ssl_reserved_sector_tics_v1",
        "created_utc": "2026-07-27T00:00:00+00:00",
        "identity_only": True,
        "light_curves_opened": False,
        "sector": 63,
        "orbits": [133, 134],
        "selection": (
            "sector == 63; orbit in {133,134}; tic_id > 0; unique TIC"
        ),
        "observations": "/readonly/twirl_observations.fits",
        "observations_sha256": "1" * 64,
        "n_selected_observation_rows": len(normalized),
        "n_reserved_tics": len(normalized),
        "reserved_tics": str(path),
        "reserved_tics_sha256": digest,
    }
    path.with_suffix(".summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )


def _inputs(tmp_path: Path) -> dict[str, object]:
    input_dir = tmp_path / "inputs"
    input_dir.mkdir(parents=True)
    registry_path = input_dir / "tic_split_registry.csv"
    both_tic, test_only_tic, s63_only_tic = _write_split_registry(
        registry_path
    )
    manifests = [
        _write_compact_export(
            input_dir,
            sector=sector,
            tics=[
                both_tic,
                test_only_tic,
                s63_only_tic,
                900_000 + sector,
            ],
        )
        for sector in EXPECTED_SECTORS
    ]
    reservation_path = input_dir / "s63_reserved_tics.txt"
    _write_s63_reservation(
        reservation_path,
        [both_tic, s63_only_tic, 999_999],
    )
    return {
        "manifests": manifests,
        "registry": registry_path,
        "reservation": reservation_path,
        "both_tic": both_tic,
        "test_only_tic": test_only_tic,
        "s63_only_tic": s63_only_tic,
    }


def _run(tmp_path: Path, inputs: dict[str, object]) -> dict[str, object]:
    return write_ssl_full_pool(
        compact_manifest_paths=inputs["manifests"],
        split_registry_path=inputs["registry"],
        s63_reserved_tics_path=inputs["reservation"],
        out_dir=tmp_path / "frozen",
    )


def test_full_pool_is_deterministic_observation_keyed_and_leakage_safe(
    tmp_path: Path,
) -> None:
    inputs = _inputs(tmp_path)

    summary = _run(tmp_path, inputs)
    out_dir = tmp_path / "frozen"
    csv_path = out_dir / "teacher_ssl_full_pool_observations.csv"
    parquet_path = out_dir / "teacher_ssl_full_pool_observations.parquet"
    summary_path = out_dir / "teacher_ssl_full_pool_manifest.summary.json"
    csv_bytes = csv_path.read_bytes()
    parquet_bytes = parquet_path.read_bytes()
    summary_bytes = summary_path.read_bytes()
    allowlist_bytes = {
        sector: (out_dir / "allowlists" / f"s{sector}_tics.csv").read_bytes()
        for sector in EXPECTED_SECTORS
    }

    pool_csv = pd.read_csv(csv_path)
    pool_parquet = pd.read_parquet(parquet_path)
    pd.testing.assert_frame_equal(pool_csv, pool_parquet, check_dtype=False)
    assert tuple(pool_csv.columns) == POOL_COLUMNS
    assert pool_csv[["sector", "tic"]].to_dict(orient="records") == [
        {"sector": sector, "tic": 900_000 + sector}
        for sector in EXPECTED_SECTORS
    ]
    assert pool_csv["pool_contract_version"].eq(
        FULL_POOL_CONTRACT_VERSION
    ).all()
    excluded = {
        int(inputs["both_tic"]),
        int(inputs["test_only_tic"]),
        int(inputs["s63_only_tic"]),
    }
    assert excluded.isdisjoint(set(pool_csv["tic"].astype(int)))
    assert pool_csv["observation_id"].is_unique
    assert pool_csv["compact_h5_sha256"].str.fullmatch(
        r"[0-9a-f]{64}"
    ).all()

    assert summary["counts"]["input"] == {
        "n_observations": 28,
        "n_unique_tics": 10,
    }
    assert summary["counts"]["excluded_input"]["n_observations"] == 21
    assert summary["counts"]["retained"] == {
        "n_observations": 7,
        "n_unique_tics": 7,
        "n_multisector_tics": 0,
    }
    categories = summary["counts"]["excluded_input"]["categories"]
    assert {
        key: value["n_observations"] for key, value in categories.items()
    } == {
        "fixed_test_only": 7,
        "s63_reserved_only": 7,
        "fixed_test_and_s63_reserved": 7,
    }
    assert summary["leakage_audit"] == {
        "fixed_test_observations_retained": 0,
        "s63_reserved_observations_retained": 0,
    }
    assert summary["inputs"]["s63_reserved_tics"]["summary_validated"] is True
    assert json.loads(summary_bytes) == summary

    for sector in EXPECTED_SECTORS:
        assert allowlist_bytes[sector] == f"tic\n{900_000 + sector}\n".encode()
        assert summary["per_sector"][str(sector)]["retained"] == {
            "n_observations": 1,
            "n_unique_tics": 1,
        }

    inputs["manifests"] = list(reversed(inputs["manifests"]))
    repeated = _run(tmp_path, inputs)
    assert repeated == summary
    assert csv_path.read_bytes() == csv_bytes
    assert parquet_path.read_bytes() == parquet_bytes
    assert summary_path.read_bytes() == summary_bytes
    assert {
        sector: (out_dir / "allowlists" / f"s{sector}_tics.csv").read_bytes()
        for sector in EXPECTED_SECTORS
    } == allowlist_bytes


def test_full_pool_fails_closed_on_missing_or_duplicate_sectors(
    tmp_path: Path,
) -> None:
    inputs = _inputs(tmp_path)
    manifests = list(inputs["manifests"])

    with pytest.raises(ValueError, match="exactly 7"):
        write_ssl_full_pool(
            compact_manifest_paths=manifests[:-1],
            split_registry_path=inputs["registry"],
            s63_reserved_tics_path=inputs["reservation"],
            out_dir=tmp_path / "missing",
        )

    duplicate_path = manifests[-1].with_name("duplicate_s56.manifest.json")
    duplicate_path.write_bytes(manifests[0].read_bytes())
    with pytest.raises(ValueError, match="duplicate compact export manifest"):
        write_ssl_full_pool(
            compact_manifest_paths=[*manifests[:-1], duplicate_path],
            split_registry_path=inputs["registry"],
            s63_reserved_tics_path=inputs["reservation"],
            out_dir=tmp_path / "duplicate",
        )


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (
            lambda manifest: manifest.update(
                requested_columns=["DET_FLUX_ADP"]
            ),
            "requested_columns must equal",
        ),
        (
            lambda manifest: manifest["records"].append(
                dict(manifest["records"][0])
            ),
            "counts disagree",
        ),
        (
            lambda manifest: manifest["skipped"].update(read_failed=1),
            "not lossless",
        ),
        (
            lambda manifest: manifest.update(unexpected_field=True),
            "wrong fields",
        ),
    ],
)
def test_full_pool_rejects_malformed_compact_manifests(
    tmp_path: Path,
    mutation,
    message: str,
) -> None:
    inputs = _inputs(tmp_path)
    path = inputs["manifests"][0]
    manifest = json.loads(path.read_text())
    mutation(manifest)
    path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )

    with pytest.raises(ValueError, match=message):
        _run(tmp_path, inputs)


def test_full_pool_rejects_h5_inventory_registry_and_s63_corruption(
    tmp_path: Path,
) -> None:
    h5_inputs = _inputs(tmp_path / "h5")
    h5_path = h5_inputs["manifests"][0].with_name(
        h5_inputs["manifests"][0].name.replace(".manifest.json", ".h5")
    )
    with h5py.File(h5_path, "a") as h5:
        del h5["targets"][next(iter(h5["targets"]))]
    with pytest.raises(ValueError, match="/targets inventory disagrees"):
        _run(tmp_path / "h5", h5_inputs)

    registry_inputs = _inputs(tmp_path / "registry")
    registry = pd.read_csv(registry_inputs["registry"])
    registry["unexpected"] = "schema drift"
    registry.to_csv(registry_inputs["registry"], index=False)
    with pytest.raises(ValueError, match="columns/order disagree"):
        _run(tmp_path / "registry", registry_inputs)

    s63_inputs = _inputs(tmp_path / "s63")
    reservation = s63_inputs["reservation"]
    lines = reservation.read_text().splitlines()
    reservation.write_text("\n".join([lines[0], lines[0], *lines[1:]]) + "\n")
    with pytest.raises(ValueError, match="strictly increasing and unique"):
        _run(tmp_path / "s63", s63_inputs)


def test_full_pool_refuses_to_replace_a_changed_artifact(
    tmp_path: Path,
) -> None:
    inputs = _inputs(tmp_path)
    _run(tmp_path, inputs)
    csv_path = (
        tmp_path / "frozen" / "teacher_ssl_full_pool_observations.csv"
    )
    csv_path.write_text(csv_path.read_text() + "corruption\n")

    with pytest.raises(FileExistsError, match="immutable output"):
        _run(tmp_path, inputs)
