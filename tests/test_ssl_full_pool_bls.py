from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

from twirl.search.a2v1_bls_contract import bls_config_sha256
from twirl.vetting.adp_only import (
    ADP_ONLY_APERTURES,
    ADP_ONLY_CONTRACT_VERSION,
)
from twirl.vetting.ssl_full_pool import (
    EXPECTED_SECTORS,
    FULL_POOL_CONTRACT_VERSION,
    FULL_POOL_SELECTION_POLICY_VERSION,
    FULL_POOL_SUMMARY_SCHEMA_VERSION,
    POOL_COLUMNS,
)
from twirl.vetting import ssl_full_pool_bls as module


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _artifact(path: Path, **extra: object) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": _sha256(path),
        **extra,
    }


def _write_pool(tmp_path: Path) -> tuple[Path, dict[int, tuple[int, ...]]]:
    root = tmp_path / "pool"
    root.mkdir(parents=True)
    rows: list[dict[str, object]] = []
    sector_tics: dict[int, tuple[int, ...]] = {}
    exports: list[dict[str, object]] = []
    allowlists: dict[str, dict[str, object]] = {}
    for sector in EXPECTED_SECTORS:
        tics = (777, 10_000 + sector)
        sector_tics[sector] = tics
        compact = root / f"s{sector}.h5"
        compact.write_bytes(f"compact-sector-{sector}\n".encode("ascii"))
        exports.append(
            {
                "sector": sector,
                "compact_h5": _artifact(compact),
            }
        )
        allowlist = root / "allowlists" / f"s{sector}_tics.csv"
        allowlist.parent.mkdir(exist_ok=True)
        allowlist.write_text(
            "tic\n" + "".join(f"{tic}\n" for tic in tics),
            encoding="ascii",
        )
        allowlists[str(sector)] = _artifact(
            allowlist,
            n_tics=len(tics),
            tic_inventory_sha256=module._tic_inventory_sha256(tics),
        )
        for index, tic in enumerate(tics):
            rows.append(
                {
                    "pool_contract_version": FULL_POOL_CONTRACT_VERSION,
                    "observation_id": f"s{sector:04d}-tic{tic:016d}",
                    "sector": sector,
                    "tic": tic,
                    "camera": index + 1,
                    "ccd": index + 1,
                    "tessmag": 16.0 + index,
                    "n_cadences": 1_000,
                    "flux_columns_json": (
                        '["DET_FLUX_ADP_SML","DET_FLUX_ADP"]'
                    ),
                    "compact_h5_path": str(compact.resolve()),
                    "compact_h5_sha256": _sha256(compact),
                    "compact_h5_size_bytes": compact.stat().st_size,
                    "compact_group_path": f"targets/{tic:016d}",
                    "compact_manifest_path": (
                        f"/frozen/s{sector}.manifest.json"
                    ),
                    "compact_manifest_sha256": "a" * 64,
                    "source_fits": f"/readonly/s{sector}/tic{tic}.fits",
                }
            )
    pool = pd.DataFrame.from_records(rows, columns=POOL_COLUMNS)
    pool = pool.sort_values(["sector", "tic"], kind="stable").reset_index(
        drop=True
    )
    csv_path = root / "teacher_ssl_full_pool_observations.csv"
    parquet_path = root / "teacher_ssl_full_pool_observations.parquet"
    pool.to_csv(csv_path, index=False, lineterminator="\n")
    pool.to_parquet(parquet_path, compression="zstd", index=False)
    identities = pool.loc[:, ["sector", "tic"]]
    summary = {
        "passed": True,
        "schema_version": FULL_POOL_SUMMARY_SCHEMA_VERSION,
        "pool_contract_version": FULL_POOL_CONTRACT_VERSION,
        "selection_policy_version": FULL_POOL_SELECTION_POLICY_VERSION,
        "sectors": list(EXPECTED_SECTORS),
        "inputs": {"compact_exports": exports},
        "outputs": {
            "csv": _artifact(csv_path, n_rows=len(pool)),
            "parquet": _artifact(parquet_path, n_rows=len(pool)),
            "sector_allowlists": allowlists,
        },
        "counts": {
            "retained": {
                "n_observations": len(pool),
                "n_unique_tics": int(pool["tic"].nunique()),
            }
        },
        "per_sector": {
            str(sector): {
                "retained": {
                    "n_observations": len(tics),
                    "n_unique_tics": len(tics),
                }
            }
            for sector, tics in sector_tics.items()
        },
        "identity_hashes": {
            "retained_observations_sha256": module._identity_sha256(
                identities
            )
        },
        "leakage_audit": {
            "fixed_test_observations_retained": 0,
            "s63_reserved_observations_retained": 0,
        },
    }
    summary_path = root / "teacher_ssl_full_pool_manifest.summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return summary_path, sector_tics


def _orbit_summary(
    frame: pd.DataFrame,
    *,
    policy: str,
) -> dict[str, object]:
    targets = (
        frame.drop_duplicates("tic")
        .sort_values("tic", kind="stable")
        .reset_index(drop=True)
    )
    records = [
        {
            "tic": int(row.tic),
            "n_cad_orbitid_reference_matched": int(
                row.n_cad_orbitid_reference_matched
            ),
            "n_cad_orbitid_mismatch": int(row.n_cad_orbitid_mismatch),
            "n_cad_orbitid_corrected": int(row.n_cad_orbitid_corrected),
            "orbitid_correction_signature_sha256": str(
                row.orbitid_correction_signature_sha256
            ),
        }
        for row in targets.itertuples(index=False)
    ]
    return {
        "n_cad_orbitid_reference_matched": int(
            targets["n_cad_orbitid_reference_matched"].sum()
        ),
        "n_cad_orbitid_mismatch": int(
            targets["n_cad_orbitid_mismatch"].sum()
        ),
        "n_cad_orbitid_corrected": int(
            targets["n_cad_orbitid_corrected"].sum()
        ),
        "n_targets_orbitid_mismatch": int(
            targets["n_cad_orbitid_mismatch"].gt(0).sum()
        ),
        "orbitid_corrections_sha256": module._canonical_json_sha256(
            {
                "contract_version": (
                    module.ORBITID_RECONCILIATION_CONTRACT_VERSION
                ),
                "policy": policy,
                "targets": records,
            }
        ),
    }


def _write_sector_bls(
    tmp_path: Path,
    *,
    pool_summary_path: Path,
    sector_tics: dict[int, tuple[int, ...]],
) -> dict[int, Path]:
    pool_summary = json.loads(pool_summary_path.read_text(encoding="utf-8"))
    exports = {
        int(item["sector"]): item
        for item in pool_summary["inputs"]["compact_exports"]
    }
    config = module._expected_bls_config()
    config_hash = bls_config_sha256(config)
    outputs: dict[int, Path] = {}
    for sector in EXPECTED_SECTORS:
        policy = module.ORBITID_POLICY_BY_SECTOR[sector]
        rows: list[dict[str, object]] = []
        for tic_index, tic in enumerate(sector_tics[sector]):
            mismatch = int(sector == 62 and tic_index == 1)
            signature = hashlib.sha256(
                (
                    f"sector={sector};tic={tic};policy={policy};"
                    f"mismatch={mismatch}"
                ).encode("ascii")
            ).hexdigest()
            for aperture in ADP_ONLY_APERTURES:
                rows.append(
                    {
                        "sector": sector,
                        "tic": tic,
                        "cam": tic_index + 1,
                        "ccd": tic_index + 1,
                        "tmag": 16.0 + tic_index,
                        "aperture": aperture,
                        "peak_rank": 1,
                        "status": "ok",
                        "period_d": 0.5 + tic_index,
                        "bls_search_branch": "current_adp",
                        "adp_only_contract_version": (
                            ADP_ONLY_CONTRACT_VERSION
                        ),
                        "source_product_tag": (
                            module.FULL_POOL_BLS_SOURCE_PRODUCT_TAG
                        ),
                        "bls_n_periods": int(config["n_periods"]),
                        "bls_n_peaks": int(config["n_peaks"]),
                        "bls_p_min_d": float(config["p_min_d"]),
                        "bls_p_max_cap_d": float(config["p_max_cap_d"]),
                        "bls_max_period_fraction": float(
                            config["max_period_fraction"]
                        ),
                        "bls_sigma_clip": float(config["sigma_clip"]),
                        "bls_orbit_edge_trim_d": float(
                            config["orbit_edge_trim_d"]
                        ),
                        "bls_search_contract_version": "custom",
                        "bls_config_sha256": config_hash,
                        "orbitid_policy": policy,
                        "orbitid_reconciliation_contract_version": (
                            module.ORBITID_RECONCILIATION_CONTRACT_VERSION
                        ),
                        "n_cad_orbitid_reference_matched": 1_000,
                        "n_cad_orbitid_mismatch": mismatch,
                        "n_cad_orbitid_corrected": mismatch,
                        "orbitid_correction_signature_sha256": signature,
                    }
                )
        frame = pd.DataFrame.from_records(rows)
        root = tmp_path / "bls" / f"s{sector}"
        root.mkdir(parents=True)
        table_path = root / "real_adp_bls_peaks.parquet"
        frame.to_parquet(table_path, compression="zstd", index=False)
        allowlist = pool_summary["outputs"]["sector_allowlists"][str(sector)]
        compact = exports[sector]["compact_h5"]
        summary_path = table_path.with_suffix(".summary.json")
        summary = {
            "passed": True,
            "sector": sector,
            "contract_version": ADP_ONLY_CONTRACT_VERSION,
            "bls_search_contract_version": "custom",
            "bls_config_sha256": config_hash,
            "compact_lc": compact["path"],
            "compact_lc_sha256": compact["sha256"],
            "cadence_reference": f"/authority/s{sector}/cadence.csv",
            "cadence_reference_sha256": "b" * 64,
            "cadence_reference_manifest": (
                f"/authority/s{sector}/cadence.manifest.json"
            ),
            "cadence_reference_manifest_sha256": "c" * 64,
            "cadence_reference_contract_version": (
                f"s{sector}_a2v1_cadence_reference_v1"
            ),
            "cadence_reference_cadence_authority": "qlp_cam_quat",
            "cadence_reference_quality_authority": (
                "spoc_and_qlp_quality_flags"
            ),
            "cadence_reference_source_hashes_sha256": "d" * 64,
            "authority_exclusion_policy_contract": (
                "a2v1_quat_absent_spoc_cadence_exclusions_v1"
            ),
            "authority_exclusion_external_bit": 62,
            "authority_exclusions_sha256": "e" * 64,
            "n_authority_exclusions": 0,
            "target_selection_contract_version": (
                module.TARGET_SELECTION_CONTRACT_VERSION
            ),
            "target_allowlist": allowlist["path"],
            "target_allowlist_sha256": allowlist["sha256"],
            "target_allowlist_count": len(sector_tics[sector]),
            "target_allowlist_tics_sha256": (
                allowlist["tic_inventory_sha256"]
            ),
            "orbitid_policy": policy,
            "orbitid_reconciliation_contract_version": (
                module.ORBITID_RECONCILIATION_CONTRACT_VERSION
            ),
            **_orbit_summary(frame, policy=policy),
            "apertures": list(ADP_ONLY_APERTURES),
            "n_targets": len(sector_tics[sector]),
            "n_targets_total": len(sector_tics[sector]),
            "n_unique_tics": len(sector_tics[sector]),
            "n_rows": len(frame),
            "n_periods": int(config["n_periods"]),
            "n_peaks": int(config["n_peaks"]),
            "n_shards": 1,
            "shard_index": 0,
            "n_source_shards": 2,
            "source_product_tag": module.FULL_POOL_BLS_SOURCE_PRODUCT_TAG,
            "config": config,
            "peak_table_sha256": _sha256(table_path),
            "outputs": {
                "peak_table": str(table_path.resolve()),
                "summary": str(summary_path.resolve()),
            },
        }
        summary_path.write_text(
            json.dumps(
                summary,
                indent=2,
                sort_keys=True,
                allow_nan=False,
            )
            + "\n",
            encoding="utf-8",
        )
        outputs[sector] = table_path
    return outputs


def _inputs(
    tmp_path: Path,
) -> tuple[Path, dict[int, Path], dict[int, tuple[int, ...]]]:
    pool_summary, sector_tics = _write_pool(tmp_path)
    sector_bls = _write_sector_bls(
        tmp_path,
        pool_summary_path=pool_summary,
        sector_tics=sector_tics,
    )
    return pool_summary, sector_bls, sector_tics


def _run(
    tmp_path: Path,
    *,
    pool_summary: Path,
    sector_bls: dict[int, Path],
) -> dict[str, object]:
    return module.write_global_full_pool_bls(
        pool_summary_path=pool_summary,
        sector_bls_paths=sector_bls,
        out_parquet_path=tmp_path / "global" / "full_pool_bls.parquet",
    )


def test_global_bls_is_deterministic_and_exactly_covers_frozen_pool(
    tmp_path: Path,
) -> None:
    pool_summary, sector_bls, _ = _inputs(tmp_path)

    summary = _run(
        tmp_path,
        pool_summary=pool_summary,
        sector_bls=sector_bls,
    )
    out_path = tmp_path / "global" / "full_pool_bls.parquet"
    summary_path = out_path.with_suffix(".summary.json")
    out_bytes = out_path.read_bytes()
    summary_bytes = summary_path.read_bytes()
    frame = pd.read_parquet(out_path)

    assert summary["passed"] is True
    assert summary["schema_version"] == (
        module.GLOBAL_BLS_SUMMARY_SCHEMA_VERSION
    )
    assert summary["counts"]["n_rows"] == len(frame) == 28
    assert summary["counts"]["n_observations"] == 14
    assert summary["counts"]["n_unique_tics"] == 8
    assert summary["counts"]["n_multisector_tics"] == 1
    assert summary["coverage_audit"]["missing_frozen_observations"] == 0
    assert summary["coverage_audit"]["unexpected_bls_observations"] == 0
    assert summary["bls_contract"]["orbitid_policy_by_sector"] == {
        **{str(sector): "strict" for sector in range(56, 62)},
        "62": "reference_by_cadence",
    }
    assert Path(summary["output"]["path"]) == out_path.resolve()
    assert summary["output"]["sha256"] == _sha256(out_path)
    assert json.loads(summary_bytes) == summary
    assert frame[["sector", "tic"]].drop_duplicates().shape[0] == 14

    repeated = module.write_global_full_pool_bls(
        pool_summary_path=pool_summary,
        sector_bls_paths=dict(reversed(list(sector_bls.items()))),
        out_parquet_path=out_path,
    )
    assert repeated == summary
    assert out_path.read_bytes() == out_bytes
    assert summary_path.read_bytes() == summary_bytes


def test_global_bls_rejects_missing_sector(tmp_path: Path) -> None:
    pool_summary, sector_bls, _ = _inputs(tmp_path)
    sector_bls.pop(62)

    with pytest.raises(ValueError, match="exactly S56--S62"):
        _run(
            tmp_path,
            pool_summary=pool_summary,
            sector_bls=sector_bls,
        )


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("compact_lc_sha256", "0" * 64, "compact-HDF5 SHA-256"),
        ("target_allowlist_count", 1, "allowlist count"),
        ("orbitid_policy", "reference_by_cadence", "must equal 'strict'"),
    ],
)
def test_global_bls_rejects_sector_summary_contract_drift(
    tmp_path: Path,
    field: str,
    value: object,
    message: str,
) -> None:
    pool_summary, sector_bls, _ = _inputs(tmp_path)
    summary_path = sector_bls[61].with_suffix(".summary.json")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary[field] = value
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match=message):
        _run(
            tmp_path,
            pool_summary=pool_summary,
            sector_bls=sector_bls,
        )


def test_global_bls_rejects_config_and_exact_coverage_drift(
    tmp_path: Path,
) -> None:
    config_root = tmp_path / "config"
    pool_summary, sector_bls, _ = _inputs(config_root)
    summary_path = sector_bls[58].with_suffix(".summary.json")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["config"]["p_min_d"] = 0.2
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="locked full-pool config"):
        _run(
            config_root,
            pool_summary=pool_summary,
            sector_bls=sector_bls,
        )

    coverage_root = tmp_path / "coverage"
    pool_summary, sector_bls, _ = _inputs(coverage_root)
    table_path = sector_bls[59]
    frame = pd.read_parquet(table_path)
    dropped_tic = int(frame["tic"].max())
    frame = frame.loc[~frame["tic"].eq(dropped_tic)].reset_index(drop=True)
    frame.to_parquet(table_path, compression="zstd", index=False)
    summary_path = table_path.with_suffix(".summary.json")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["peak_table_sha256"] = _sha256(table_path)
    summary["n_rows"] = len(frame)
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="coverage disagrees"):
        _run(
            coverage_root,
            pool_summary=pool_summary,
            sector_bls=sector_bls,
        )


def test_global_bls_rejects_changed_frozen_pool_and_output(
    tmp_path: Path,
) -> None:
    pool_root = tmp_path / "pool-drift"
    pool_summary, sector_bls, _ = _inputs(pool_root)
    summary = json.loads(pool_summary.read_text(encoding="utf-8"))
    summary["outputs"]["parquet"]["sha256"] = "0" * 64
    pool_summary.write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="Parquet SHA-256"):
        _run(
            pool_root,
            pool_summary=pool_summary,
            sector_bls=sector_bls,
        )

    immutable_root = tmp_path / "immutable"
    pool_summary, sector_bls, _ = _inputs(immutable_root)
    _run(
        immutable_root,
        pool_summary=pool_summary,
        sector_bls=sector_bls,
    )
    out_path = immutable_root / "global" / "full_pool_bls.parquet"
    out_path.write_bytes(out_path.read_bytes() + b"corruption")
    with pytest.raises(FileExistsError, match="immutable output"):
        _run(
            immutable_root,
            pool_summary=pool_summary,
            sector_bls=sector_bls,
        )
