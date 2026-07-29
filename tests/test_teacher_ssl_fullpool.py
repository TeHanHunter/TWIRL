from __future__ import annotations

import importlib.util
from pathlib import Path
import hashlib
import json
import os
import random
import time

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
                # Keep the synthetic SHA string-valued through CSV round trips.
                # A digits-only 64-character fixture can be inferred as numeric.
                "native_h5_sha256": (
                    "a" + f"{int(row.sector) % 10}" * 63
                ),
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
        **fullpool.MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS,
        "expected_git_sha": "1" * 40,
        "schema_version": (
            fullpool.MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS["schema_version"]
        ),
        "release_binding": fullpool.FULL_POOL_NATIVE_RELEASE_BINDING,
        "native_contract_version": fullpool.FULL_POOL_NATIVE_CONTRACT_VERSION,
        "builder_contract_version": (
            fullpool.FULL_POOL_NATIVE_BUILDER_CONTRACT_VERSION
        ),
        "detrend_contract_version": (
            fullpool.FULL_POOL_NATIVE_DETREND_CONTRACT_VERSION
        ),
        "detrend_config_sha256": (
            fullpool.FULL_POOL_NATIVE_DETREND_CONFIG_SHA256
        ),
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
            "numeric_gate_release",
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
    numeric_gate_release = {
        "schema_version": fullpool.MODEL_INPUT_NUMERIC_RELEASE_SCHEMA,
        "release_binding": fullpool.MODEL_INPUT_NUMERIC_RELEASE_BINDING,
        "native_freshness": dict(
            {
                **fullpool.MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS,
                "expected_git_sha": "1" * 40,
            }
        ),
        "passed": True,
        "model_input_contract_version": (
            fullpool.MODEL_INPUT_CONTRACT_VERSION
        ),
        "envelope_contract": (
            fullpool.MODEL_INPUT_NUMERIC_ENVELOPE_CONTRACT
        ),
        "envelope_canonical_sha256": "9" * 64,
        "counts": {
            "full_observations": len(full_keys),
            "eligible_observations": len(eligible_keys),
            "excluded_observations": len(excluded_keys),
            "scanned_observations": len(eligible_keys),
            "failed_observations": 0,
            "native_shards": native_registry["native_h5_path"].nunique(),
        },
        "identity_hashes": {
            "full": full_identity,
            "eligible": eligible_identity,
            "excluded": excluded_identity,
        },
        "code_revision": "1" * 40,
    }
    return (
        eligibility,
        native_registry,
        native_release,
        registry,
        registry_summary,
        numeric_gate_release,
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


def test_bounded_selection_pins_required_observation_without_changing_default(
    tmp_path: Path,
) -> None:
    registry, _ = _build(tmp_path)
    production_default, default_audit = fullpool.select_fullpool_ssl_fold(
        registry,
        held_out_fold=1,
    )
    production_explicit_empty, empty_audit = (
        fullpool.select_fullpool_ssl_fold(
            registry,
            held_out_fold=1,
            required_observation_ids=[],
        )
    )
    pd.testing.assert_frame_equal(
        production_default,
        production_explicit_empty,
    )
    assert default_audit == empty_audit

    baseline, _ = fullpool.select_fullpool_ssl_fold(
        registry,
        held_out_fold=1,
        max_rows=2,
    )
    forced_id = next(
        value
        for value in production_default["ssl_observation_id"].astype(str)
        if value not in set(baseline["ssl_observation_id"].astype(str))
    )
    selected, audit = fullpool.select_fullpool_ssl_fold(
        registry,
        held_out_fold=1,
        max_rows=2,
        required_observation_ids=[forced_id],
    )

    assert len(selected) == 2
    assert forced_id in set(selected["ssl_observation_id"].astype(str))
    assert audit["required_observation_ids"] == [forced_id]
    assert audit["n_required_observations"] == 1
    assert audit["required_observations_selected"] is True


def test_required_observation_selection_fails_closed(
    tmp_path: Path,
) -> None:
    registry, _ = _build(tmp_path)
    eligible_id = str(
        registry.loc[
            registry["ssl_pool_include"], "ssl_observation_id"
        ].iloc[0]
    )
    held_id = str(
        registry.loc[
            registry["ssl_held_out_fold"].eq(2),
            "ssl_observation_id",
        ].iloc[0]
    )

    with pytest.raises(ValueError, match="only valid with bounded max_rows"):
        fullpool.select_fullpool_ssl_fold(
            registry,
            held_out_fold=1,
            required_observation_ids=[eligible_id],
        )
    with pytest.raises(ValueError, match="contains duplicates"):
        fullpool.select_fullpool_ssl_fold(
            registry,
            held_out_fold=1,
            max_rows=2,
            required_observation_ids=[eligible_id, eligible_id],
        )
    with pytest.raises(ValueError, match="absent from the full-pool registry"):
        fullpool.select_fullpool_ssl_fold(
            registry,
            held_out_fold=1,
            max_rows=2,
            required_observation_ids=["s9999-tic9999999999999999999"],
        )
    with pytest.raises(ValueError, match="held-out fold"):
        fullpool.select_fullpool_ssl_fold(
            registry,
            held_out_fold=2,
            max_rows=2,
            required_observation_ids=[held_id],
        )


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
        numeric_gate_release,
        authority_bindings,
    ) = _training_authority_inputs(tmp_path, monkeypatch)

    audit = fullpool._validate_training_authority_chain(
        eligibility=eligibility,
        native_registry=native_registry,
        native_release=native_release,
        registry=registry,
        registry_summary=registry_summary,
        numeric_gate_release=numeric_gate_release,
        authority_bindings=authority_bindings,
    )
    assert audit["production_lock_passed"] is True
    assert audit["source_provenance_verified"] is True
    assert audit["model_input_contract_version"] == (
        fullpool.MODEL_INPUT_CONTRACT_VERSION
    )
    assert audit["model_input_numeric_envelope_contract"] == (
        fullpool.MODEL_INPUT_NUMERIC_ENVELOPE_CONTRACT
    )
    assert audit["numeric_gate_release"]["binding"] == (
        authority_bindings["numeric_gate_release"]
    )
    assert audit["numeric_gate_release"]["release_binding"] == (
        fullpool.MODEL_INPUT_NUMERIC_RELEASE_BINDING
    )
    assert audit["native_freshness"]["native_contract_version"] == (
        fullpool.FULL_POOL_NATIVE_CONTRACT_VERSION
    )
    assert set(audit["authority_bindings"]) == {
        "eligibility_exclusions",
        "eligibility_summary",
        "native_registry",
        "native_registry_summary",
        "native_release_summary",
        "registry",
        "registry_summary",
        "numeric_gate_release",
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
            numeric_gate_release=numeric_gate_release,
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
            numeric_gate_release=numeric_gate_release,
            authority_bindings=authority_bindings,
        )

    wrong_numeric_contract = dict(numeric_gate_release)
    wrong_numeric_contract["model_input_contract_version"] = "old"
    with pytest.raises(ValueError, match="wrong input contract"):
        fullpool._validate_training_authority_chain(
            eligibility=eligibility,
            native_registry=native_registry,
            native_release=native_release,
            registry=registry,
            registry_summary=registry_summary,
            numeric_gate_release=wrong_numeric_contract,
            authority_bindings=authority_bindings,
        )


def test_selected_native_verification_requires_actual_v3_group_identity(
    tmp_path: Path,
) -> None:
    h5py = pytest.importorskip("h5py")
    native_path = tmp_path / "native.h5"
    with h5py.File(native_path, "w") as h5:
        h5.attrs["contract_version"] = (
            fullpool.FULL_POOL_NATIVE_CONTRACT_VERSION
        )
        h5.attrs["expected_git_sha"] = "1" * 40
        h5.attrs["release_binding"] = (
            fullpool.FULL_POOL_NATIVE_RELEASE_BINDING
        )
        h5.attrs["builder_contract_version"] = (
            fullpool.FULL_POOL_NATIVE_BUILDER_CONTRACT_VERSION
        )
        h5.attrs["builder_code_sha256"] = (
            fullpool.MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS[
                "builder_code_sha256"
            ]
        )
        h5.attrs["detrend_contract_version"] = (
            fullpool.FULL_POOL_NATIVE_DETREND_CONTRACT_VERSION
        )
        h5.attrs["detrend_config_sha256"] = (
            fullpool.FULL_POOL_NATIVE_DETREND_CONFIG_SHA256
        )
        h5.attrs["detrend_quality_source"] = "final_effective_quality"
        h5.attrs["detrend_time_contract_version"] = (
            fullpool.MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS[
                "detrend_time_contract_version"
            ]
        )
        h5.attrs["detrend_time_dataset"] = (
            fullpool.MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS[
                "detrend_time_dataset"
            ]
        )
        h5.attrs["detrend_time_system"] = "BTJD"
        h5.attrs["published_time_system"] = "BJD"
        h5.attrs["btjd_to_bjd_offset_d"] = (
            fullpool.MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS[
                "btjd_to_bjd_offset_d"
            ]
        )
        h5.attrs["warning_capture_policy"] = (
            fullpool.MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS[
                "warning_capture_policy"
            ]
        )
        h5.attrs["rank_warning_publication_policy"] = (
            fullpool.MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS[
                "rank_warning_publication_policy"
            ]
        )
        h5.attrs["rank_warning_count"] = 0
        h5.attrs["rank_warning_ledger_json"] = "[]"
        h5.attrs["rank_warning_ledger_sha256"] = (
            fullpool.MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS[
                "rank_warning_ledger_sha256"
            ]
        )
        h5.attrs["raw_photometry_only"] = 1
        h5.attrs["compact_adp_photometry_reused"] = 0
        h5.attrs["compact_adp_flux_reused"] = 0
        h5.attrs["periodogram_n"] = (
            fullpool.MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS["periodogram_n"]
        )
        group = h5.create_group("targets/123")
        group.attrs["sector"] = 56
        group.attrs["tic"] = 123
        for name in fullpool.PERIODOGRAM_DATASETS:
            group.create_dataset(
                name,
                data=np.ones(
                    fullpool.MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS[
                        "periodogram_n"
                    ],
                    dtype=np.float32,
                ),
            )
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
    records = fullpool._verify_native_files(
        selected,
        expected_git_sha="1" * 40,
    )
    assert records[0]["hash_verified_now"] is True
    assert records[0]["root_contract_verified_now"] is True
    assert records[0]["group_identities_verified_now"] is True

    with h5py.File(native_path, "r+") as h5:
        h5.attrs["contract_version"] = "native-v1"
    selected.loc[0, "native_h5_sha256"] = hashlib.sha256(
        native_path.read_bytes()
    ).hexdigest()
    with pytest.raises(ValueError, match="not fresh native-v3"):
        fullpool._verify_native_files(
            selected,
            expected_git_sha="1" * 40,
        )

    with h5py.File(native_path, "r+") as h5:
        h5.attrs["contract_version"] = (
            fullpool.FULL_POOL_NATIVE_CONTRACT_VERSION
        )
        h5["targets/123"].attrs["tic"] = 999
    selected.loc[0, "native_h5_sha256"] = hashlib.sha256(
        native_path.read_bytes()
    ).hexdigest()
    with pytest.raises(ValueError, match="group identity differs"):
        fullpool._verify_native_files(
            selected,
            expected_git_sha="1" * 40,
        )


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
        "--numeric-gate-release",
        "--required-observation-id",
    ):
        assert option in script
    assert 'action="append"' in script
    assert "--skip-native-hash-verification" not in script
    assert "--skip-numeric-gate" not in script


def test_effective_quality_adp_model_run_uses_fresh_v3r1_namespace() -> None:
    assert fullpool.FULLPOOL_SSL_REGISTRY_SCHEMA.endswith("_v3")
    assert fullpool.FULLPOOL_SSL_REGISTRY_SUMMARY_SCHEMA.endswith("_v3")
    assert fullpool.FULLPOOL_SSL_SELECTION_SCHEMA.endswith("_v5")
    assert fullpool.FULLPOOL_SSL_RUN_CONTRACT_SCHEMA.endswith("_v5")
    assert fullpool.FULLPOOL_SSL_RESUME_SCHEMA.endswith("_v5")
    assert fullpool.FULLPOOL_SSL_CHECKPOINT_SCHEMA.endswith("_v5")
    assert fullpool.FULLPOOL_SSL_SUMMARY_SCHEMA.endswith("_v5")
    assert fullpool.FULLPOOL_SSL_TRAINING_AUTHORITY_SCHEMA.endswith("_v5")
    assert "effective_quality_adp_btjd_v2" in fullpool.FULLPOOL_SSL_RUN_ID
    assert "effective_quality_adp_btjd_v2" in (
        fullpool.FULLPOOL_SSL_ENCODER_NAME
    )
    assert (
        "effective_quality_adp_btjd_v2"
        in fullpool.FULLPOOL_SSL_CHECKPOINT_NAMESPACE
    )


def test_locked_smoke_and_fold_roots_reject_stale_namespaces(
    tmp_path: Path,
) -> None:
    smoke_root = (
        tmp_path
        / "model_runs"
        / fullpool.FULLPOOL_SSL_MODEL_NAMESPACE
        / "smoke"
        / "one_epoch"
    )
    observed = fullpool._validate_locked_training_shape(
        out_root=smoke_root,
        fold=2,
        epochs=1,
        batch_size=64,
        workers=8,
        seed=fullpool.FULLPOOL_SSL_DEFAULT_TRAINING_SEED,
        learning_rate=fullpool.FULLPOOL_SSL_LEARNING_RATE,
        weight_decay=fullpool.FULLPOOL_SSL_WEIGHT_DECAY,
        checkpoint_every=1,
        require_cuda=True,
        max_rows=4096,
        required_observation_ids=[
            fullpool.FULLPOOL_SSL_SMOKE_REQUIRED_OBSERVATION_ID
        ],
    )
    assert observed == smoke_root.resolve()

    stale_root = (
        tmp_path
        / "model_runs"
        / "effective_quality_mask_v1"
        / "smoke"
        / "one_epoch"
    )
    with pytest.raises(ValueError, match="fresh"):
        fullpool._validate_locked_training_shape(
            out_root=stale_root,
            fold=2,
            epochs=1,
            batch_size=64,
            workers=8,
            seed=fullpool.FULLPOOL_SSL_DEFAULT_TRAINING_SEED,
            learning_rate=fullpool.FULLPOOL_SSL_LEARNING_RATE,
            weight_decay=fullpool.FULLPOOL_SSL_WEIGHT_DECAY,
            checkpoint_every=1,
            require_cuda=True,
            max_rows=4096,
            required_observation_ids=[
                fullpool.FULLPOOL_SSL_SMOKE_REQUIRED_OBSERVATION_ID
            ],
        )
    with pytest.raises(ValueError, match="exactly once"):
        fullpool._validate_locked_training_shape(
            out_root=smoke_root,
            fold=2,
            epochs=1,
            batch_size=64,
            workers=8,
            seed=fullpool.FULLPOOL_SSL_DEFAULT_TRAINING_SEED,
            learning_rate=fullpool.FULLPOOL_SSL_LEARNING_RATE,
            weight_decay=fullpool.FULLPOOL_SSL_WEIGHT_DECAY,
            checkpoint_every=1,
            require_cuda=True,
            max_rows=4096,
            required_observation_ids=[],
        )

    production_root = (
        tmp_path
        / "model_runs"
        / fullpool.FULLPOOL_SSL_MODEL_NAMESPACE
        / "training"
        / "five_fold"
    )
    for field, replacement in (
        ("learning_rate", fullpool.FULLPOOL_SSL_LEARNING_RATE * 2.0),
        ("weight_decay", fullpool.FULLPOOL_SSL_WEIGHT_DECAY * 2.0),
    ):
        kwargs = {
            "out_root": production_root,
            "fold": 0,
            "epochs": 20,
            "batch_size": 64,
            "workers": 8,
            "seed": fullpool.FULLPOOL_SSL_DEFAULT_TRAINING_SEED,
            "learning_rate": fullpool.FULLPOOL_SSL_LEARNING_RATE,
            "weight_decay": fullpool.FULLPOOL_SSL_WEIGHT_DECAY,
            "checkpoint_every": 1,
            "require_cuda": True,
            "max_rows": None,
            "required_observation_ids": [],
        }
        kwargs[field] = replacement
        with pytest.raises(ValueError, match=field.replace("_", " ")):
            fullpool._validate_locked_training_shape(**kwargs)


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


def test_epoch_coverage_locks_singleton_skip_rule() -> None:
    assert fullpool._expected_epoch_coverage(64, 64) == {
        "selected_observations": 64,
        "observations_per_epoch": 64,
        "optimizer_steps_per_epoch": 1,
        "singleton_observations_skipped_per_epoch": 0,
    }
    assert fullpool._expected_epoch_coverage(65, 64) == {
        "selected_observations": 65,
        "observations_per_epoch": 64,
        "optimizer_steps_per_epoch": 1,
        "singleton_observations_skipped_per_epoch": 1,
    }
    assert fullpool._expected_epoch_coverage(66, 64) == {
        "selected_observations": 66,
        "observations_per_epoch": 66,
        "optimizer_steps_per_epoch": 2,
        "singleton_observations_skipped_per_epoch": 0,
    }


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


def test_resume_checkpoint_requires_current_exact_native_freshness(
    tmp_path: Path,
) -> None:
    torch = pytest.importorskip("torch")
    fold_dir = tmp_path / "fold_0"
    fold_dir.mkdir()
    contract_sha256 = "a" * 64
    current_freshness = {
        **fullpool.MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS,
        "expected_git_sha": "1" * 40,
    }
    stale_freshness = {
        **current_freshness,
        "expected_git_sha": "2" * 40,
    }
    checkpoint = {
        "schema_version": fullpool.FULLPOOL_SSL_RESUME_SCHEMA,
        "run_id": fullpool.FULLPOOL_SSL_RUN_ID,
        "checkpoint_namespace": fullpool.FULLPOOL_SSL_CHECKPOINT_NAMESPACE,
        "numeric_release_binding": fullpool.MODEL_INPUT_NUMERIC_RELEASE_BINDING,
        "native_freshness": stale_freshness,
        "run_contract_sha256": contract_sha256,
        "epoch": 1,
        "global_step": 10,
    }
    fullpool._publish_resume_checkpoint(
        fold_dir=fold_dir,
        contract_sha256=contract_sha256,
        checkpoint=checkpoint,
        epoch=1,
        global_step=10,
    )
    with pytest.raises(ValueError, match="native freshness"):
        fullpool._load_resume_checkpoint(
            fold_dir=fold_dir,
            contract_sha256=contract_sha256,
            device=torch.device("cpu"),
            expected_native_freshness=current_freshness,
        )


def test_optimizer_tracks_exact_real_ssl_embedding_gradient_path() -> None:
    torch = pytest.importorskip("torch")
    from twirl.vetting.harmonic_cnn import (
        HarmonicModelConfig,
        build_harmonic_cnn,
    )
    from twirl.vetting.harmonic_ssl import vicreg_loss

    config = HarmonicModelConfig(metadata_dim=0)
    model = build_harmonic_cnn(
        config,
        profile=fullpool.FULLPOOL_SSL_PROFILE,
    )
    projector = fullpool._projection_head(2 * config.embedding_dim)
    named = fullpool._optimizer_named_parameters(model, projector)
    optimizer = torch.optim.AdamW(
        [
            {
                "params": [parameter for _, parameter in named],
                "decoupled_weight_decay": True,
            }
        ],
        lr=fullpool.FULLPOOL_SSL_LEARNING_RATE,
        betas=fullpool.FULLPOOL_SSL_ADAMW_BETAS,
        eps=fullpool.FULLPOOL_SSL_ADAMW_EPS,
        weight_decay=fullpool.FULLPOOL_SSL_WEIGHT_DECAY,
        amsgrad=False,
        foreach=None,
        maximize=False,
        capturable=False,
        differentiable=False,
        fused=None,
    )
    batch_size = 4

    def batch() -> dict[str, object]:
        return {
            "harmonic_values": torch.randn(batch_size, 7, 7, 37),
            "harmonic_mask": torch.ones(
                batch_size, 7, 7, 37, dtype=torch.bool
            ),
            "local_values": torch.randn(batch_size, 14, 7, 19),
            "local_mask": torch.ones(
                batch_size, 14, 7, 19, dtype=torch.bool
            ),
            "periodogram_values": torch.randn(batch_size, 4, 31),
            "periodogram_mask": torch.ones(
                batch_size, 4, 31, dtype=torch.bool
            ),
            "metadata": torch.empty(batch_size, 0),
        }

    first = projector(model(batch())["embedding"])
    second = projector(model(batch())["embedding"])
    loss = vicreg_loss(first, second)
    if isinstance(loss, tuple):
        loss = loss[0]
    loss.backward()
    coverage = fullpool._assert_ssl_gradient_coverage(
        model,
        projector,
        optimizer,
        fold=0,
        epoch=1,
        batch_index=0,
    )
    optimizer.step()
    contract = fullpool._optimizer_checkpoint_contract(
        model,
        projector,
        optimizer,
        gradient_coverage_verified=True,
        expected_step=1,
    )

    assert len(list(model.parameters())) + len(
        list(projector.parameters())
    ) == 146
    assert coverage["active_parameter_count"] == 90
    assert coverage["excluded_parameter_count"] == 56
    assert contract["parameter_count"] == 90
    assert contract["state_parameter_count"] == 90
    assert contract["parameter_scope"] == (
        fullpool.FULLPOOL_SSL_OPTIMIZER_PARAMETER_SCOPE
    )
    assert optimizer.state_dict()["param_groups"][0][
        "decoupled_weight_decay"
    ] is True
    active_names = {
        item["name"] for item in contract["parameter_manifest"]
    }
    assert not any(
        name.startswith(
            (
                "encoder.small_encoder.",
                "encoder.supp_encoder.",
                "encoder.supplement_gate.",
                "encoder.morphology_head.",
                "encoder.preserve_head.",
                "encoder.harmonic_head.",
            )
        )
        for name in active_names
    )


def test_optimizer_contract_rejects_inactive_task_head_parameter() -> None:
    torch = pytest.importorskip("torch")
    from twirl.vetting.harmonic_cnn import (
        HarmonicModelConfig,
        build_harmonic_cnn,
    )

    config = HarmonicModelConfig(metadata_dim=0)
    model = build_harmonic_cnn(
        config,
        profile=fullpool.FULLPOOL_SSL_PROFILE,
    )
    projector = fullpool._projection_head(2 * config.embedding_dim)
    optimizer = torch.optim.AdamW(
        [
            {
                "params": (
                    list(model.parameters())
                    + list(projector.parameters())
                ),
                "decoupled_weight_decay": True,
            }
        ],
        lr=fullpool.FULLPOOL_SSL_LEARNING_RATE,
        betas=fullpool.FULLPOOL_SSL_ADAMW_BETAS,
        eps=fullpool.FULLPOOL_SSL_ADAMW_EPS,
        weight_decay=fullpool.FULLPOOL_SSL_WEIGHT_DECAY,
    )
    with pytest.raises(RuntimeError, match="embedding path"):
        fullpool._optimizer_checkpoint_contract(
            model,
            projector,
            optimizer,
            gradient_coverage_verified=True,
            expected_step=1,
        )


def test_training_input_revalidation_detects_late_native_mutation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    h5py = pytest.importorskip("h5py")
    native_path = tmp_path / "native.h5"
    with h5py.File(native_path, "w") as h5:
        h5.attrs["expected_git_sha"] = "1" * 40
        h5.attrs["periodogram_n"] = fullpool.FULL_POOL_NATIVE_PERIODOGRAM_N
    stat = native_path.stat()
    record = {
        "native_h5_path": str(native_path),
        "native_h5_sha256": hashlib.sha256(
            native_path.read_bytes()
        ).hexdigest(),
        "native_h5_size_bytes": stat.st_size,
        "native_h5_mtime_ns": stat.st_mtime_ns,
        "native_h5_ctime_ns": stat.st_ctime_ns,
        "native_h5_device": stat.st_dev,
        "native_h5_inode": stat.st_ino,
    }
    monkeypatch.setattr(
        fullpool,
        "full_pool_native_root_failures",
        lambda attrs: [],
    )
    fullpool._revalidate_training_inputs(
        authority_bindings={},
        native_files=[record],
        expected_git_sha="1" * 40,
        full_native_hash=False,
    )

    # Some ORCD node-local filesystems coalesce ctime updates within one
    # timestamp tick.  Cross a full tick so this test exercises the metadata
    # revalidation path rather than filesystem timestamp granularity.
    time.sleep(1.1)
    with native_path.open("r+b") as handle:
        handle.seek(-1, os.SEEK_END)
        original = handle.read(1)
        handle.seek(-1, os.SEEK_END)
        handle.write(bytes([original[0] ^ 1]))
    os.utime(
        native_path,
        ns=(stat.st_atime_ns, stat.st_mtime_ns),
    )
    with pytest.raises(RuntimeError, match="changed after preflight"):
        fullpool._revalidate_training_inputs(
            authority_bindings={},
            native_files=[record],
            expected_git_sha="1" * 40,
            full_native_hash=False,
        )


def _load_v3_completion_validator() -> object:
    script = (
        Path(__file__).resolve().parents[1]
        / "scripts"
        / "stage5_validation"
        / "validate_teacher_ssl_fullpool_v3_training.py"
    )
    spec = importlib.util.spec_from_file_location(
        "validate_teacher_ssl_fullpool_v3_training_native_test",
        script,
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_completion_validator_real_native_record_path(tmp_path: Path) -> None:
    h5py = pytest.importorskip("h5py")
    validator = _load_v3_completion_validator()
    revision = "1" * 40
    freshness = {
        **fullpool.MODEL_INPUT_NUMERIC_NATIVE_FRESHNESS,
        "expected_git_sha": revision,
    }
    native_path = tmp_path / "native.h5"
    with h5py.File(native_path, "w") as h5:
        h5.attrs["contract_version"] = freshness[
            "native_contract_version"
        ]
        h5.attrs["release_binding"] = freshness["release_binding"]
        h5.attrs["expected_git_sha"] = revision
        h5.attrs["builder_contract_version"] = freshness[
            "builder_contract_version"
        ]
        h5.attrs["builder_code_sha256"] = freshness[
            "builder_code_sha256"
        ]
        h5.attrs["detrend_contract_version"] = freshness[
            "detrend_contract_version"
        ]
        h5.attrs["detrend_config_sha256"] = freshness[
            "detrend_config_sha256"
        ]
        h5.attrs["detrend_quality_source"] = "final_effective_quality"
        h5.attrs["detrend_time_contract_version"] = freshness[
            "detrend_time_contract_version"
        ]
        h5.attrs["detrend_time_dataset"] = freshness[
            "detrend_time_dataset"
        ]
        h5.attrs["detrend_time_system"] = "BTJD"
        h5.attrs["published_time_system"] = "BJD"
        h5.attrs["btjd_to_bjd_offset_d"] = freshness[
            "btjd_to_bjd_offset_d"
        ]
        h5.attrs["warning_capture_policy"] = freshness[
            "warning_capture_policy"
        ]
        h5.attrs["rank_warning_publication_policy"] = freshness[
            "rank_warning_publication_policy"
        ]
        h5.attrs["rank_warning_count"] = 0
        h5.attrs["rank_warning_ledger_json"] = freshness[
            "rank_warning_ledger_json"
        ]
        h5.attrs["rank_warning_ledger_sha256"] = freshness[
            "rank_warning_ledger_sha256"
        ]
        h5.attrs["raw_photometry_only"] = 1
        h5.attrs["compact_adp_photometry_reused"] = 0
        h5.attrs["compact_adp_flux_reused"] = 0
        h5.attrs["periodogram_n"] = freshness["periodogram_n"]
    stat = native_path.stat()
    record = {
        "native_h5_path": str(native_path.resolve()),
        "native_h5_sha256": hashlib.sha256(
            native_path.read_bytes()
        ).hexdigest(),
        "native_h5_size_bytes": stat.st_size,
        "native_h5_mtime_ns": stat.st_mtime_ns,
        "native_h5_ctime_ns": stat.st_ctime_ns,
        "native_h5_device": stat.st_dev,
        "native_h5_inode": stat.st_ino,
        "native_contract_version": freshness["native_contract_version"],
        "native_release_binding": freshness["release_binding"],
        "native_expected_git_sha": revision,
        "native_builder_contract_version": freshness[
            "builder_contract_version"
        ],
        "native_builder_code_sha256": freshness[
            "builder_code_sha256"
        ],
        "native_detrend_contract_version": freshness[
            "detrend_contract_version"
        ],
        "native_detrend_config_sha256": freshness[
            "detrend_config_sha256"
        ],
        "native_detrend_quality_source": "final_effective_quality",
        "native_detrend_time_contract_version": freshness[
            "detrend_time_contract_version"
        ],
        "native_detrend_time_dataset": freshness["detrend_time_dataset"],
        "native_detrend_time_system": "BTJD",
        "native_published_time_system": "BJD",
        "native_btjd_to_bjd_offset_d": freshness[
            "btjd_to_bjd_offset_d"
        ],
        "native_warning_capture_policy": freshness[
            "warning_capture_policy"
        ],
        "native_rank_warning_publication_policy": freshness[
            "rank_warning_publication_policy"
        ],
        "native_rank_warning_count": 0,
        "native_rank_warning_ledger_json": freshness[
            "rank_warning_ledger_json"
        ],
        "native_rank_warning_ledger_sha256": freshness[
            "rank_warning_ledger_sha256"
        ],
        "raw_photometry_only": True,
        "compact_adp_photometry_reused": False,
        "compact_adp_flux_reused": False,
        "periodogram_n": freshness["periodogram_n"],
        "hash_verified_now": True,
        "root_contract_verified_now": True,
        "group_identities_verified_now": True,
    }
    observed = validator._validate_native_record(
        record,
        expected_code_revision=revision,
        index=0,
    )
    assert observed["sha256"] == record["native_h5_sha256"]


def test_completion_metadata_detects_ctime_only_mutation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    validator = _load_v3_completion_validator()
    artifact = tmp_path / "authority.json"
    artifact.write_text("{}\n", encoding="ascii")
    real_sha256 = validator._sha256

    def mutate_after_hash(path: Path) -> str:
        digest = real_sha256(path)
        # See the ctime-granularity note in the late-native-mutation test.
        time.sleep(1.1)
        path.chmod(0o600)
        return digest

    monkeypatch.setattr(validator, "_sha256", mutate_after_hash)
    with pytest.raises(RuntimeError, match="changed while hashing"):
        validator._metadata(artifact)
