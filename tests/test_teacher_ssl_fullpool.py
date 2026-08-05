from __future__ import annotations

from pathlib import Path
import hashlib
import json
import random

import numpy as np
import pandas as pd
import pytest

import twirl.vetting.teacher_ssl_fullpool as fullpool
from twirl.vetting.ssl_full_pool_eligibility import (
    ArtifactBinding,
    EligibilityAuthority,
    observation_identity_sha256,
)


def _bls_row(
    *,
    sector: int,
    tic: int,
    aperture: str = fullpool.FULLPOOL_SSL_ANCHOR_APERTURE,
    peak_rank: int = 1,
    status: str = "ok",
) -> dict[str, object]:
    valid = peak_rank == 1 and status == "ok"
    return {
        "sector": sector,
        "tic": tic,
        "aperture": aperture,
        "peak_rank": peak_rank,
        "status": status,
        "period_d": 1.0 + tic / 1.0e6 if valid else np.nan,
        "t0_bjd": 2_457_000.0 + sector if valid else np.nan,
        "duration_min": 30.0 if valid else np.nan,
    }


def _registry_inputs(
    tmp_path: Path,
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
]:
    small_rows: list[dict[str, object]] = []
    for sector in fullpool.FULLPOOL_SSL_SECTORS:
        tic = 10_000 + sector
        small_rows.append(_bls_row(sector=sector, tic=tic))
        small_rows.append(
            _bls_row(
                sector=sector,
                tic=tic,
                aperture="DET_FLUX_ADP",
            )
        )
    # Frozen development host observed in two sectors.
    small_rows.extend(
        [
            _bls_row(sector=56, tic=20_000),
            _bls_row(sector=57, tic=20_000),
        ]
    )
    # One explicit rank-zero status row remains visible in the audit.
    small_rows.append(
        _bls_row(
            sector=61,
            tic=50_000,
            peak_rank=0,
            status="too_few_cadences",
        )
    )
    bls = pd.DataFrame(small_rows)

    native_records: list[dict[str, object]] = []
    anchor = bls.loc[
        bls["aperture"].eq(fullpool.FULLPOOL_SSL_ANCHOR_APERTURE)
        & bls["peak_rank"].eq(1)
    ]
    for row in anchor.itertuples(index=False):
        native_records.append(
            {
                "sector": int(row.sector),
                "tic": int(row.tic),
                "native_h5_path": str(
                    (tmp_path / f"s{int(row.sector)}.h5").resolve()
                ),
                "native_group_path": (
                    f"/observations/s{int(row.sector):04d}/tic{int(row.tic)}"
                ),
                "native_h5_sha256": f"{int(row.sector) % 10}" * 64,
                "native_contract_version": (
                    fullpool.FULL_POOL_NATIVE_CONTRACT_VERSION
                ),
            }
        )
    native = pd.DataFrame(native_records)
    split = pd.DataFrame(
        [
            {"tic": 20_000, "fixed_split": "development", "cv_fold": 2},
            {"tic": 30_000, "fixed_split": "test", "cv_fold": -1},
        ]
    )
    reserved = pd.DataFrame({"tic": [40_000, 99_999]})
    anchor = bls.loc[
        bls["aperture"].eq(fullpool.FULLPOOL_SSL_ANCHOR_APERTURE)
        & bls["peak_rank"].isin([0, 1])
    ].sort_values(["sector", "tic"], kind="stable")
    pool_records: list[dict[str, object]] = []
    for row in anchor.itertuples(index=False):
        sector = int(row.sector)
        tic = int(row.tic)
        pool_records.append(
            {
                "pool_contract_version": fullpool.FULL_POOL_CONTRACT_VERSION,
                "observation_id": f"s{sector:04d}-tic{tic:016d}",
                "sector": sector,
                "tic": tic,
                "camera": 1,
                "ccd": 1,
                "tessmag": 18.0,
                "n_cadences": 1_000,
                "flux_columns_json": (
                    '["DET_FLUX_ADP_SML", "DET_FLUX_ADP"]'
                ),
                "compact_h5_path": str(
                    (tmp_path / f"s{sector}_compact.h5").resolve()
                ),
                "compact_h5_sha256": "e" * 64,
                "compact_h5_size_bytes": 10_000,
                "compact_group_path": f"/targets/{tic:016d}",
                "compact_manifest_path": str(
                    (tmp_path / f"s{sector}.manifest.json").resolve()
                ),
                "compact_manifest_sha256": "f" * 64,
                "source_fits": f"/readonly/s{sector}/tic{tic}.fits",
            }
        )
    pool = pd.DataFrame(
        pool_records,
        columns=fullpool.POOL_COLUMNS,
    )
    return bls, native, split, reserved, pool


def _build(
    tmp_path: Path,
) -> tuple[pd.DataFrame, dict[str, object]]:
    bls, native, split, reserved, pool = _registry_inputs(tmp_path)
    return fullpool.build_fullpool_ssl_registry(
        bls,
        native,
        split,
        reserved,
        frozen_pool=pool,
    )


def _artifact(tmp_path: Path, name: str) -> ArtifactBinding:
    path = tmp_path / name
    path.write_bytes(f"{name}\n".encode("ascii"))
    return ArtifactBinding(
        path=path.resolve(),
        size_bytes=path.stat().st_size,
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
    )


def _training_authority_inputs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> tuple[
    EligibilityAuthority,
    pd.DataFrame,
    dict[str, object],
    pd.DataFrame,
    dict[str, object],
    dict[str, dict[str, object]],
]:
    registry, registry_summary = _build(tmp_path)
    include = registry["ssl_pool_include"]
    full_keys = frozenset(
        zip(registry["sector"].astype(int), registry["tic"].astype(int))
    )
    eligible_keys = frozenset(
        zip(
            registry.loc[include, "sector"].astype(int),
            registry.loc[include, "tic"].astype(int),
        )
    )
    excluded_keys = frozenset(full_keys - eligible_keys)
    full_identity = observation_identity_sha256(full_keys)
    eligible_identity = observation_identity_sha256(eligible_keys)
    excluded_identity = observation_identity_sha256(excluded_keys)
    for name, value in (
        ("PRODUCTION_FULL_OBSERVATIONS", len(full_keys)),
        ("PRODUCTION_ELIGIBLE_OBSERVATIONS", len(eligible_keys)),
        ("PRODUCTION_EXCLUDED_OBSERVATIONS", len(excluded_keys)),
        ("PRODUCTION_FULL_IDENTITY_SHA256", full_identity),
        ("PRODUCTION_ELIGIBLE_IDENTITY_SHA256", eligible_identity),
        ("PRODUCTION_EXCLUDED_IDENTITY_SHA256", excluded_identity),
    ):
        monkeypatch.setattr(fullpool, name, value)

    source_bindings = {
        name: _artifact(tmp_path, f"{name}.dat")
        for name in (
            "frozen_pool",
            "frozen_pool_summary",
            "global_bls",
            "global_bls_summary",
            "exclusions",
            "eligibility_summary",
        )
    }
    eligibility = EligibilityAuthority(
        contract_version="eligibility-test-v2",
        release_binding="eligibility-test-release",
        anchor_aperture=fullpool.FULLPOOL_SSL_ANCHOR_APERTURE,
        full_keys=full_keys,
        eligible_keys=eligible_keys,
        excluded_keys=excluded_keys,
        n_full=len(full_keys),
        n_eligible=len(eligible_keys),
        n_excluded=len(excluded_keys),
        full_observation_identity_sha256=full_identity,
        eligible_observation_identity_sha256=eligible_identity,
        excluded_observation_identity_sha256=excluded_identity,
        bindings=source_bindings,
        exclusions=pd.DataFrame(),
        summary={},
    )
    _, native_registry, _, _, _ = _registry_inputs(tmp_path)
    native_registry = native_registry.sort_values(
        ["sector", "tic"], kind="stable"
    ).reset_index(drop=True)
    native_release = {
        "schema_version": "native-release-test-v2",
        "release_binding": "native-release-test",
        "native_contract_version": fullpool.FULL_POOL_NATIVE_CONTRACT_VERSION,
        "counts": {
            "full_observations": len(full_keys),
            "eligible_observations": len(eligible_keys),
            "excluded_observations": len(excluded_keys),
            "native_registry_observations": len(eligible_keys),
        },
        "identity_hashes": {
            "full_observations_sha256": full_identity,
            "eligible_observations_sha256": eligible_identity,
            "excluded_observations_sha256": excluded_identity,
            "native_registry_observations_sha256": eligible_identity,
        },
    }
    explicit_bindings = {
        name: _artifact(tmp_path, f"{name}.authority").metadata()
        for name in (
            "native_registry",
            "native_registry_summary",
            "native_release_summary",
            "registry",
            "registry_summary",
        )
    }
    explicit_bindings["eligibility_exclusions"] = source_bindings[
        "exclusions"
    ].metadata()
    explicit_bindings["eligibility_summary"] = source_bindings[
        "eligibility_summary"
    ].metadata()
    split_binding = _artifact(tmp_path, "split-registry.csv").metadata()
    reserved_binding = _artifact(tmp_path, "reserved-hosts.txt").metadata()
    registry_summary = dict(registry_summary)
    registry_summary["source_provenance"] = {
        "frozen_pool": source_bindings["frozen_pool"].metadata(),
        "frozen_pool_summary": source_bindings[
            "frozen_pool_summary"
        ].metadata(),
        "bls_summary": source_bindings["global_bls_summary"].metadata(),
        "eligibility_exclusions": explicit_bindings[
            "eligibility_exclusions"
        ],
        "eligibility_summary": explicit_bindings["eligibility_summary"],
        "native_registry": explicit_bindings["native_registry"],
        "native_registry_summary": explicit_bindings[
            "native_registry_summary"
        ],
        "native_release_summary": explicit_bindings[
            "native_release_summary"
        ],
        "frozen_split_registry": split_binding,
        "reserved_hosts": reserved_binding,
        "frozen_pool_authority_bindings": {
            "split_registry_sha256_equal": True,
            "reserved_hosts_sha256_equal": True,
            "split_registry_sha256": split_binding["sha256"],
            "reserved_hosts_sha256": reserved_binding["sha256"],
        },
        "global_bls_authority_binding": {
            "summary_path": str(
                source_bindings["global_bls_summary"].path
            ),
            "summary_sha256": source_bindings["global_bls_summary"].sha256,
            "artifact_path": str(source_bindings["global_bls"].path),
            "artifact_sha256": source_bindings["global_bls"].sha256,
            "artifact_matches_summary": True,
            "summary_schema_version": (
                fullpool.GLOBAL_BLS_SUMMARY_SCHEMA_VERSION
            ),
            "summary_contract_version": fullpool.GLOBAL_BLS_CONTRACT_VERSION,
        },
        "native_model_eligibility_binding": {
            "contract_version": eligibility.contract_version,
            "release_binding": eligibility.release_binding,
            "full_observations": len(full_keys),
            "eligible_observations": len(eligible_keys),
            "excluded_observations": len(excluded_keys),
            "full_observation_identity_sha256": full_identity,
            "eligible_observation_identity_sha256": eligible_identity,
            "excluded_observation_identity_sha256": excluded_identity,
        },
        "native_release_binding": {
            "schema_version": native_release["schema_version"],
            "release_binding": native_release["release_binding"],
            "native_contract_version": native_release[
                "native_contract_version"
            ],
            "release_summary_sha256": explicit_bindings[
                "native_release_summary"
            ]["sha256"],
        },
    }
    return (
        eligibility,
        native_registry,
        native_release,
        registry,
        registry_summary,
        explicit_bindings,
    )


def test_fullpool_registry_is_label_free_and_excludes_hosts(
    tmp_path: Path,
) -> None:
    registry, audit = _build(tmp_path)

    assert not (set(registry) & fullpool._FORBIDDEN_LABEL_COLUMNS)
    assert audit["labels_consumed"] is False
    assert audit["fixed_test_exclusion_scope"] == "whole_tic"
    assert audit["prospective_exclusion_scope"] == "whole_tic"
    assert not registry["tic"].isin([30_000, 40_000]).any()
    assert not registry["fixed_test_member"].any()
    assert not registry["reserved_prospective_member"].any()
    assert registry.loc[
        registry["tic"].eq(50_000), "ssl_pool_exclusion_reason"
    ].eq("bls_unsearchable+native_missing").all()


def test_unseen_tics_are_present_in_all_folds_and_frozen_fold_is_held(
    tmp_path: Path,
) -> None:
    registry, _ = _build(tmp_path)
    unknown = registry[
        registry["tic"].between(10_056, 10_062)
    ]
    assert unknown["ssl_pool_include"].all()
    assert unknown["ssl_held_out_fold"].eq(-1).all()
    assert unknown["fold_assignment_source"].eq("unlabeled_all_folds").all()

    selected_fold_2, audit_2 = fullpool.select_fullpool_ssl_fold(
        registry,
        held_out_fold=2,
    )
    selected_fold_1, audit_1 = fullpool.select_fullpool_ssl_fold(
        registry,
        held_out_fold=1,
    )
    assert not selected_fold_2["tic"].eq(20_000).any()
    assert selected_fold_1["tic"].eq(20_000).sum() == 2
    for selected in (selected_fold_2, selected_fold_1):
        assert set(unknown["tic"]).issubset(set(selected["tic"]))
        assert not selected["fixed_test_member"].any()
        assert not selected["reserved_prospective_member"].any()
    assert audit_2["n_held_tics"] == 1
    assert audit_1["n_held_tics"] == 0
    assert all(audit_2["tic_disjoint"].values())


def test_fullpool_dataset_adapter_uses_no_supervised_targets(
    tmp_path: Path,
) -> None:
    registry, _ = _build(tmp_path)
    selected, _ = fullpool.select_fullpool_ssl_fold(
        registry,
        held_out_fold=2,
    )
    rows = fullpool.fullpool_dataset_rows(selected)

    assert "human_label" not in rows
    assert rows["review_id"].equals(rows["ssl_observation_id"])
    assert rows["morphology_target_index"].eq(-1).all()
    assert rows["preserve_target_index"].eq(-1).all()
    assert rows["harmonic_target_index"].eq(-1).all()
    assert rows["morphology_weight"].eq(0).all()
    assert not rows["is_injected_row"].any()


def test_fullpool_registry_summary_is_hash_bound_and_immutable(
    tmp_path: Path,
) -> None:
    registry, audit = _build(tmp_path)
    registry_path = tmp_path / "fullpool.csv"
    summary_path = tmp_path / "fullpool.summary.json"
    result = fullpool.write_fullpool_ssl_registry(
        registry,
        audit,
        registry_path=registry_path,
        summary_path=summary_path,
        source_provenance={"test": {"sha256": "a" * 64}},
    )
    loaded, summary = fullpool.load_fullpool_ssl_registry(
        registry_path=registry_path,
        summary_path=summary_path,
    )
    assert len(loaded) == len(registry)
    assert summary["registry_sha256"] == result["registry_sha256"]

    # Byte-identical retries are accepted.
    repeated = fullpool.write_fullpool_ssl_registry(
        registry,
        audit,
        registry_path=registry_path,
        summary_path=summary_path,
        source_provenance={"test": {"sha256": "a" * 64}},
    )
    assert repeated["registry_sha256"] == result["registry_sha256"]

    registry_path.write_bytes(registry_path.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="hash does not match"):
        fullpool.load_fullpool_ssl_registry(
            registry_path=registry_path,
            summary_path=summary_path,
        )


def test_training_authority_chain_is_exact_and_provenance_bound(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    (
        eligibility,
        native_registry,
        native_release,
        registry,
        registry_summary,
        authority_bindings,
    ) = _training_authority_inputs(tmp_path, monkeypatch)

    audit = fullpool._validate_training_authority_chain(
        eligibility=eligibility,
        native_registry=native_registry,
        native_release=native_release,
        registry=registry,
        registry_summary=registry_summary,
        authority_bindings=authority_bindings,
    )
    assert audit["production_lock_passed"] is True
    assert audit["source_provenance_verified"] is True
    assert set(audit["authority_bindings"]) == {
        "eligibility_exclusions",
        "eligibility_summary",
        "native_registry",
        "native_registry_summary",
        "native_release_summary",
        "registry",
        "registry_summary",
    }

    changed_mapping = registry.copy()
    first_included = changed_mapping["ssl_pool_include"].idxmax()
    changed_mapping.loc[
        first_included, "native_group_path"
    ] = "targets/wrong"
    with pytest.raises(ValueError, match="native mapping differs"):
        fullpool._validate_training_authority_chain(
            eligibility=eligibility,
            native_registry=native_registry,
            native_release=native_release,
            registry=changed_mapping,
            registry_summary=registry_summary,
            authority_bindings=authority_bindings,
        )

    changed_provenance = json.loads(json.dumps(registry_summary))
    changed_provenance["source_provenance"]["native_registry"][
        "sha256"
    ] = "0" * 64
    with pytest.raises(ValueError, match="provenance native_registry differs"):
        fullpool._validate_training_authority_chain(
            eligibility=eligibility,
            native_registry=native_registry,
            native_release=native_release,
            registry=registry,
            registry_summary=changed_provenance,
            authority_bindings=authority_bindings,
        )


def test_selected_native_verification_requires_actual_v2_group_identity(
    tmp_path: Path,
) -> None:
    h5py = pytest.importorskip("h5py")
    native_path = tmp_path / "native.h5"
    with h5py.File(native_path, "w") as h5:
        h5.attrs["contract_version"] = (
            fullpool.FULL_POOL_NATIVE_CONTRACT_VERSION
        )
        group = h5.create_group("targets/123")
        group.attrs["sector"] = 56
        group.attrs["tic"] = 123
    selected = pd.DataFrame(
        {
            "sector": [56],
            "tic": [123],
            "native_h5_path": [str(native_path)],
            "native_group_path": ["targets/123"],
            "native_h5_sha256": [
                hashlib.sha256(native_path.read_bytes()).hexdigest()
            ],
            "native_contract_version": [
                fullpool.FULL_POOL_NATIVE_CONTRACT_VERSION
            ],
        }
    )
    records = fullpool._verify_native_files(selected)
    assert records[0]["hash_verified_now"] is True
    assert records[0]["root_contract_verified_now"] is True
    assert records[0]["group_identities_verified_now"] is True

    with h5py.File(native_path, "r+") as h5:
        h5.attrs["contract_version"] = "native-v1"
    selected.loc[0, "native_h5_sha256"] = hashlib.sha256(
        native_path.read_bytes()
    ).hexdigest()
    with pytest.raises(ValueError, match="not native-v2"):
        fullpool._verify_native_files(selected)

    with h5py.File(native_path, "r+") as h5:
        h5.attrs["contract_version"] = (
            fullpool.FULL_POOL_NATIVE_CONTRACT_VERSION
        )
        h5["targets/123"].attrs["tic"] = 999
    selected.loc[0, "native_h5_sha256"] = hashlib.sha256(
        native_path.read_bytes()
    ).hexdigest()
    with pytest.raises(ValueError, match="group identity differs"):
        fullpool._verify_native_files(selected)


def test_training_cli_requires_full_authority_chain_and_has_no_hash_bypass() -> None:
    script = (
        Path(__file__).resolve().parents[1]
        / "scripts"
        / "stage5_validation"
        / "train_teacher_ssl_fullpool_fold.py"
    ).read_text(encoding="utf-8")
    for option in (
        "--eligibility-exclusions",
        "--eligibility-summary",
        "--native-registry",
        "--native-registry-summary",
        "--native-release-summary",
        "--registry",
        "--registry-summary",
    ):
        assert option in script
    assert "--skip-native-hash-verification" not in script


def test_fold_run_directory_resumes_only_identical_contract(
    tmp_path: Path,
) -> None:
    fold_dir = tmp_path / "fold_0"
    contract = {"schema_version": "test", "fold": 0, "value": 1}
    digest, existed = fullpool._prepare_fold_directory(
        fold_dir=fold_dir,
        run_contract=contract,
        resume=False,
    )
    assert len(digest) == 64
    assert existed is False

    repeated_digest, existed = fullpool._prepare_fold_directory(
        fold_dir=fold_dir,
        run_contract=contract,
        resume=True,
    )
    assert repeated_digest == digest
    assert existed is True

    with pytest.raises(RuntimeError, match="contract differs"):
        fullpool._prepare_fold_directory(
            fold_dir=fold_dir,
            run_contract={**contract, "value": 2},
            resume=True,
        )
    with pytest.raises(FileExistsError, match="not empty"):
        fullpool._prepare_fold_directory(
            fold_dir=fold_dir,
            run_contract=contract,
            resume=False,
        )


def test_fullpool_builder_rejects_label_columns_and_missing_native(
    tmp_path: Path,
) -> None:
    bls, native, split, reserved, pool = _registry_inputs(tmp_path)
    with pytest.raises(ValueError, match="label-free"):
        fullpool.build_fullpool_ssl_registry(
            bls.assign(human_label="uncertain"),
            native,
            split,
            reserved,
            frozen_pool=pool,
        )
    unknown_key = (
        native["tic"].between(10_056, 10_062)
        & native["sector"].eq(56)
    )
    with pytest.raises(ValueError, match="lack full-pool native mappings"):
        fullpool.build_fullpool_ssl_registry(
            bls,
            native.loc[~unknown_key],
            split,
            reserved,
            frozen_pool=pool,
        )

    with pytest.raises(ValueError, match="coverage differs"):
        fullpool.build_fullpool_ssl_registry(
            bls.loc[~(bls["sector"].eq(56) & bls["tic"].eq(10_056))],
            native,
            split,
            reserved,
            frozen_pool=pool,
        )


def test_text_tic_inventory_is_strict_sorted_and_unique(
    tmp_path: Path,
) -> None:
    path = tmp_path / "reserved.txt"
    path.write_text("10\n20\n30\n", encoding="ascii")
    assert fullpool.read_tic_inventory(path)["tic"].tolist() == [10, 20, 30]

    path.write_text("20\n10\n", encoding="ascii")
    with pytest.raises(ValueError, match="sorted and unique"):
        fullpool.read_tic_inventory(path)


def test_frozen_pool_consumer_validates_summary_and_allowlists(
    tmp_path: Path,
) -> None:
    *_, pool = _registry_inputs(tmp_path)
    bundle = tmp_path / "frozen"
    allowlists = bundle / "allowlists"
    allowlists.mkdir(parents=True)
    pool_path = bundle / "teacher_ssl_full_pool_observations.csv"
    pool.to_csv(pool_path, index=False, lineterminator="\n")
    allowlist_metadata: dict[str, dict[str, object]] = {}
    for sector in fullpool.FULLPOOL_SSL_SECTORS:
        path = allowlists / f"s{sector}_tics.csv"
        tics = (
            pool.loc[pool["sector"].eq(sector), "tic"]
            .astype(int)
            .sort_values()
            .tolist()
        )
        payload = ("tic\n" + "".join(f"{tic}\n" for tic in tics)).encode()
        path.write_bytes(payload)
        allowlist_metadata[str(sector)] = {
            "path": str(path),
            "size_bytes": len(payload),
            "sha256": hashlib.sha256(payload).hexdigest(),
            "n_tics": len(tics),
            "tic_inventory_sha256": "unused_by_consumer",
        }
    summary = {
        "schema_version": fullpool.FULL_POOL_SUMMARY_SCHEMA_VERSION,
        "pool_contract_version": fullpool.FULL_POOL_CONTRACT_VERSION,
        "sectors": list(fullpool.FULLPOOL_SSL_SECTORS),
        "counts": {
            "retained": {
                "n_observations": len(pool),
                "n_unique_tics": pool["tic"].nunique(),
                "n_multisector_tics": int(
                    pool.groupby("tic")["sector"].nunique().gt(1).sum()
                ),
            }
        },
        "identity_hashes": {
            "retained_observations_sha256": fullpool._pool_identity_sha256(
                pool
            )
        },
        "leakage_audit": {
            "fixed_test_observations_retained": 0,
            "s63_reserved_observations_retained": 0,
        },
        "outputs": {
            "csv": {
                "path": str(pool_path),
                "size_bytes": pool_path.stat().st_size,
                "sha256": hashlib.sha256(pool_path.read_bytes()).hexdigest(),
                "n_rows": len(pool),
            },
            "sector_allowlists": allowlist_metadata,
        },
    }
    summary_path = bundle / "teacher_ssl_full_pool_manifest.summary.json"
    summary_path.write_text(
        json.dumps(summary, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    loaded, _ = fullpool.load_frozen_ssl_full_pool(
        pool_path=pool_path,
        summary_path=summary_path,
    )
    assert len(loaded) == len(pool)

    (allowlists / "s56_tics.csv").write_text("tic\n999\n", encoding="ascii")
    with pytest.raises(ValueError, match="allowlist hash mismatch"):
        fullpool.load_frozen_ssl_full_pool(
            pool_path=pool_path,
            summary_path=summary_path,
        )


def test_global_bls_consumer_requires_summary_authorized_parquet(
    tmp_path: Path,
) -> None:
    pytest.importorskip("pyarrow")
    bls, *_, pool = _registry_inputs(tmp_path)
    pool_summary_path = tmp_path / "frozen-pool.summary.json"
    pool_summary_path.write_text("{}\n", encoding="utf-8")
    bls_path = tmp_path / "global.parquet"
    bls.to_parquet(bls_path, index=False)
    identity = fullpool._pool_identity_sha256(pool)
    n_observations = len(pool)
    summary = {
        "passed": True,
        "schema_version": fullpool.GLOBAL_BLS_SUMMARY_SCHEMA_VERSION,
        "contract_version": fullpool.GLOBAL_BLS_CONTRACT_VERSION,
        "sectors": list(fullpool.FULLPOOL_SSL_SECTORS),
        "observation_identity_columns": ["sector", "tic"],
        "frozen_pool": {
            "summary": {
                "path": str(pool_summary_path),
                "size_bytes": pool_summary_path.stat().st_size,
                "sha256": hashlib.sha256(
                    pool_summary_path.read_bytes()
                ).hexdigest(),
            },
            "contract_version": fullpool.FULL_POOL_CONTRACT_VERSION,
            "selection_policy_version": "test",
            "observation_identity_sha256": identity,
            "n_observations": n_observations,
            "n_unique_tics": pool["tic"].nunique(),
        },
        "bls_contract": {},
        "counts": {
            "n_rows": len(bls),
            "n_observations": n_observations,
            "n_unique_tics": pool["tic"].nunique(),
            "n_multisector_tics": 1,
            "status_rows": {},
            "aperture_rows": {},
        },
        "coverage_audit": {
            "missing_frozen_observations": 0,
            "unexpected_bls_observations": 0,
            "observation_identity_sha256": identity,
        },
        "sector_products": [],
        "output": {
            "path": str(bls_path),
            "size_bytes": bls_path.stat().st_size,
            "sha256": hashlib.sha256(bls_path.read_bytes()).hexdigest(),
            "n_rows": len(bls),
            "n_observations": n_observations,
        },
    }
    summary_path = tmp_path / "global.summary.json"
    summary_path.write_text(
        json.dumps(summary, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    loaded, _, authorized_path = fullpool.load_global_full_pool_bls(
        summary_path=summary_path,
        frozen_pool=pool,
        frozen_pool_summary_path=pool_summary_path,
    )
    assert len(loaded) == len(bls)
    assert authorized_path == bls_path.resolve()

    summary["counts"]["n_rows"] += 1
    summary_path.write_text(
        json.dumps(summary, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="row count mismatch"):
        fullpool.load_global_full_pool_bls(
            summary_path=summary_path,
            frozen_pool=pool,
            frozen_pool_summary_path=pool_summary_path,
        )


def test_rng_state_round_trip_supports_exact_epoch_resume() -> None:
    torch = pytest.importorskip("torch")
    random.seed(4)
    np.random.seed(5)
    torch.manual_seed(6)
    state = fullpool._capture_rng_state()
    expected = (
        random.random(),
        np.random.random(),
        torch.rand(3),
    )
    random.random()
    np.random.random()
    torch.rand(3)
    fullpool._restore_rng_state(state)
    observed = (
        random.random(),
        np.random.random(),
        torch.rand(3),
    )
    assert observed[0] == expected[0]
    assert observed[1] == expected[1]
    assert torch.equal(observed[2], expected[2])


def test_training_numerics_gates_reject_nonfinite_loss_and_state() -> None:
    torch = pytest.importorskip("torch")
    finite = torch.tensor(1.0, requires_grad=True)
    fullpool._assert_finite_loss_and_components(
        finite,
        {"invariance": 0.5},
        fold=0,
        epoch=1,
        batch_index=0,
    )
    with pytest.raises(FloatingPointError, match="non-finite.*loss"):
        fullpool._assert_finite_loss_and_components(
            torch.tensor(float("nan"), requires_grad=True),
            {},
            fold=0,
            epoch=1,
            batch_index=0,
        )
    with pytest.raises(FloatingPointError, match="components"):
        fullpool._assert_finite_loss_and_components(
            finite,
            {"variance": float("inf")},
            fold=0,
            epoch=1,
            batch_index=0,
        )

    module = torch.nn.Linear(2, 1)
    fullpool._assert_finite_module_state(
        {"encoder": module},
        fold=0,
        epoch=1,
    )
    with torch.no_grad():
        module.weight[0, 0] = float("inf")
    with pytest.raises(FloatingPointError, match="model state"):
        fullpool._assert_finite_module_state(
            {"encoder": module},
            fold=0,
            epoch=1,
        )


def test_resume_pointer_survives_crash_between_checkpoint_and_state(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fold_dir = tmp_path / "fold_0"
    fold_dir.mkdir()
    contract_sha256 = "a" * 64
    monkeypatch.setattr(
        fullpool,
        "_atomic_torch_save",
        lambda payload, path: fullpool._atomic_write(
            path,
            json.dumps(payload, sort_keys=True).encode("utf-8"),
        ),
    )

    first = {
        "schema_version": fullpool.FULLPOOL_SSL_RESUME_SCHEMA,
        "run_contract_sha256": contract_sha256,
        "epoch": 1,
        "global_step": 10,
        "marker": 1,
    }
    first_path, _ = fullpool._publish_resume_checkpoint(
        fold_dir=fold_dir,
        contract_sha256=contract_sha256,
        checkpoint=first,
        epoch=1,
        global_step=10,
    )
    state = json.loads(
        (fold_dir / "resume_state.json").read_text(encoding="utf-8")
    )
    assert state["checkpoint"] == first_path.name

    # Simulate a process dying after generation 2 is durable but before the
    # authoritative resume_state.json pointer advances.
    second = {
        "schema_version": fullpool.FULLPOOL_SSL_RESUME_SCHEMA,
        "run_contract_sha256": contract_sha256,
        "epoch": 2,
        "global_step": 20,
        "marker": 2,
    }
    orphan_path = fold_dir / "resume_epoch_0002_step_000000000020.pt"
    fullpool._atomic_torch_save(second, orphan_path)
    state = json.loads(
        (fold_dir / "resume_state.json").read_text(encoding="utf-8")
    )
    assert state["checkpoint"] == first_path.name
    assert json.loads(first_path.read_text(encoding="utf-8"))["marker"] == 1

    second_path, _ = fullpool._publish_resume_checkpoint(
        fold_dir=fold_dir,
        contract_sha256=contract_sha256,
        checkpoint=second,
        epoch=2,
        global_step=20,
    )
    assert second_path == orphan_path
    state = json.loads(
        (fold_dir / "resume_state.json").read_text(encoding="utf-8")
    )
    assert state["checkpoint"] == second_path.name
    assert json.loads(second_path.read_text(encoding="utf-8"))["marker"] == 2

    with pytest.raises(RuntimeError, match="authoritative"):
        fullpool._publish_resume_checkpoint(
            fold_dir=fold_dir,
            contract_sha256=contract_sha256,
            checkpoint=second,
            epoch=2,
            global_step=20,
        )
