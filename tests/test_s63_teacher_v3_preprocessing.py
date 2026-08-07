from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import subprocess

import h5py
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]


def _load_script(relative: str, name: str):
    path = ROOT / relative
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, sort_keys=True) + "\n", encoding="utf-8")


def _accepted_stage1() -> dict:
    return {
        "sector": 63,
        "orbits": [133, 134],
        "ok": True,
        "ok_h5": True,
        "ok_fits": True,
        "n_expected_h5": 4,
        "n_expected_unique_tics": 3,
        "schema": {
            "expected_method": "A2v1",
            "expected_prodtag": "A2v1",
            "schema_only": True,
            "check_h5_open": True,
            "allow_edge_warn_missing": True,
            "required_a2v1_columns": ["DET_FLUX_ADP_SML", "DET_FLUX_ADP"],
        },
        "expected_contract": {
            "ok": True,
            "requested_orbits": [133, 134],
            "missing_requested_orbits": [],
        },
        "h5": {
            "n_missing_h5_non_edge": 0,
            "n_missing_h5_edge_warn": 1,
            "n_zero_byte_h5": 0,
            "n_unreadable_h5": 0,
        },
        "fits": {
            "hlsp_root": "/tmp/twirl-test-s63-hlsp",
            "n_fits": 3,
            "n_found_unique_tics": 3,
            "n_extra_fits_tics": 0,
            "n_missing_fits_non_edge_tics": 0,
            "n_missing_fits_edge_warn_tics": 1,
            "n_bad_checked_fits": 0,
        },
    }


def test_s63_authorization_is_explicit_and_does_not_widen_legacy() -> None:
    native = _load_script(
        "scripts/stage5_validation/build_teacher_v3_native_inputs.py",
        "s63_native_authorization_test",
    )
    merge = _load_script(
        "scripts/stage5_validation/merge_teacher_v3_native_shards.py",
        "s63_merge_authorization_test",
    )
    raw = _load_script(
        "scripts/stage5_validation/export_teacher_v3_raw_sources.py",
        "s63_raw_authorization_test",
    )
    contract = "s63_teacher_v3_prospective_v1"

    assert native._authorize_sector(sector=62, prospective_contract=None) == (
        "teacher_v3_s56_s62_legacy"
    )
    with pytest.raises(ValueError, match="must not"):
        native._authorize_sector(sector=62, prospective_contract=contract)
    with pytest.raises(ValueError, match="requires"):
        merge._authorize_sector(sector=63, prospective_contract=None)
    assert merge._authorize_sector(
        sector=63, prospective_contract=contract
    ) == contract
    with pytest.raises(ValueError, match="133/134"):
        raw._authorize_sector(sector=63, orbits=[132, 133], contract=contract)
    assert raw._authorize_sector(
        sector=63, orbits=[134, 133], contract=contract
    ) == contract
    with pytest.raises(ValueError, match="bounded"):
        native._authorize_sector(sector=64, prospective_contract=contract)


def test_s63_stage1_preflight_is_strict(tmp_path: Path) -> None:
    module = _load_script(
        "scripts/stage5_validation/build_s63_teacher_v3_launch_manifest.py",
        "s63_stage1_preflight_test",
    )
    path = tmp_path / "accepted_stage1.json"
    _write_json(path, _accepted_stage1())
    audit = module.validate_accepted_stage1(path)
    assert audit["passed"] is True
    assert audit["orbits"] == [133, 134]
    assert audit["expected_contract_evidence"] == "embedded_expected_contract"

    legacy = _accepted_stage1()
    legacy.pop("expected_contract")
    legacy["expected_h5_by_orbit"] = {"133": 16, "134": 16}
    legacy["expected_h5_by_ccd"] = {
        f"o{orbit}_cam{camera}_ccd{ccd}": 1
        for orbit in (133, 134)
        for camera in range(1, 5)
        for ccd in range(1, 5)
    }
    legacy["n_expected_h5"] = 32
    _write_json(path, legacy)
    legacy_audit = module.validate_accepted_stage1(path)
    assert legacy_audit["passed"] is True
    assert legacy_audit["expected_contract_evidence"] == (
        "legacy_exact_orbit_detector_counts"
    )

    missing_legacy_evidence = _accepted_stage1()
    missing_legacy_evidence.pop("expected_contract")
    _write_json(path, missing_legacy_evidence)
    with pytest.raises(ValueError, match="expected_h5_by_orbit"):
        module.validate_accepted_stage1(path)

    corrupt_legacy_counts = dict(legacy)
    corrupt_legacy_counts["expected_h5_by_orbit"] = {"133": 16, "134": 17}
    _write_json(path, corrupt_legacy_counts)
    with pytest.raises(ValueError, match="does not sum"):
        module.validate_accepted_stage1(path)

    failed = _accepted_stage1()
    failed["h5"]["n_unreadable_h5"] = 1
    _write_json(path, failed)
    with pytest.raises(ValueError, match="n_unreadable_h5"):
        module.validate_accepted_stage1(path)


def test_s63_preprocessing_shell_assets_are_fail_closed() -> None:
    paths = (
        "scripts/stage1_lightcurves/run_a2v1_cadence_reference_pdo.sh",
        "scripts/stage5_validation/run_teacher_v3_raw_export_pdo.sh",
        "scripts/stage5_validation/run_s63_teacher_v3_preprocessing_pdo.sh",
        "scripts/orcd/slurm_teacher_v3_native_sector_cpu.sbatch",
        "scripts/orcd/slurm_teacher_v3_s63_bls_cpu.sbatch",
        "scripts/orcd/slurm_teacher_v3_s63_bls_merge_cpu.sbatch",
        "scripts/orcd/slurm_teacher_v3_s63_native_cpu.sbatch",
        "scripts/orcd/slurm_teacher_v3_s63_native_merge_cpu.sbatch",
        "scripts/orcd/slurm_teacher_v3_s63_score_h200.sbatch",
        "scripts/orcd/slurm_teacher_v3_s63_launch_cpu.sbatch",
        "scripts/orcd/slurm_teacher_v3_s63_render_review_cpu.sbatch",
        "scripts/orcd/slurm_teacher_v3_s63_queue_cpu.sbatch",
    )
    for relative in paths:
        completed = subprocess.run(
            ["bash", "-n", str(ROOT / relative)],
            check=False,
            capture_output=True,
            text=True,
        )
        assert completed.returncode == 0, f"{relative}: {completed.stderr}"

    controller = (ROOT / paths[2]).read_text(encoding="utf-8")
    assert "50000" in controller
    assert "--n-peaks 10" in controller
    assert "--orbitid-policy strict" in controller
    assert "DET_FLUX_ADP_SML DET_FLUX_ADP" in controller
    assert "TWIRL_EXPECTED_GIT_SHA" in controller
    assert "mkdir \"${ACTIVE_CLAIM}\"" in controller
    assert "PYTHONNOUSERSITE=1" in controller
    assert 'PDO_USER_ROOT="/pdo/users/tehan"' in controller
    assert 'canonical_path="$(realpath -m -- "${!output_var}")"' in controller
    assert '"${canonical_path}" != "${canonical_user_root}"/*' in controller
    assert '"${canonical_path}" == "${REPO}"/*' in controller
    for output_var in (
        "VALIDATION",
        "CADENCE_TABLE",
        "COMPACT_H5",
        "MODEL_READY_ALLOWLIST",
        "BLS_MERGED",
        "CANDIDATES",
        "RAW_SOURCE_H5",
        "NATIVE_H5",
        "LAUNCH_MANIFEST",
    ):
        assert output_var in controller
    assert (
        "/sw/python-versions/python-3.11.9/lib:"
        "/pdo/app/anaconda/anaconda2-4.4.0/lib"
    ) in controller
    assert "No stage is submitted by this controller" in controller

    cadence_wrapper = (ROOT / paths[0]).read_text(encoding="utf-8")
    raw_wrapper = (ROOT / paths[1]).read_text(encoding="utf-8")
    for wrapper in (cadence_wrapper, raw_wrapper):
        assert 'PDO_USER_ROOT="$(realpath -e -- /pdo/users/tehan)"' in wrapper
        assert 'REPO="$(realpath -e -- "${REPO}")"' in wrapper
        assert '"${PDO_USER_ROOT}"/*' in wrapper
        assert '"${REPO}"/*' in wrapper
    assert 'OUTPUT_DIR="$(realpath -m -- "${OUTPUT_DIR}")"' in cadence_wrapper
    assert 'OUT_H5="$(realpath -m -- "${OUT_H5}")"' in raw_wrapper

    bls = (ROOT / paths[4]).read_text(encoding="utf-8")
    assert "#SBATCH --array=0-15%8" in bls
    assert "#SBATCH -c 8" in bls
    assert "--n-periods 50000" in bls
    assert "--n-peaks 10" in bls
    assert "TWIRL_S63_CADENCE_TABLE_SHA256" in bls
    assert "TWIRL_S63_CADENCE_MANIFEST_SHA256" in bls
    assert "scripts/assert_clean_checkout.sh" in bls

    native = (ROOT / paths[6]).read_text(encoding="utf-8")
    assert "#SBATCH --array=0-31%16" in native
    assert "--n-periods 4096" in native
    assert "--orbitid-policy strict" in native
    assert "scripts/assert_clean_checkout.sh" in native

    score = (ROOT / paths[8]).read_text(encoding="utf-8")
    assert score.count("#SBATCH --gres=gpu:h200:1") == 1
    assert "TWIRL_S63_LAUNCH_MANIFEST_SHA256" in score
    assert "TWIRL_TEACHER_V3_CHECKPOINT_MANIFEST_SHA256" in score
    assert "scripts/assert_clean_checkout.sh" in score
    assert '--repo "${REPO}"' in score
    assert '--expected-git-sha "${EXPECTED_GIT_SHA}"' in score
    assert "--allow-cpu" not in score
    assert "umask 077" in score
    assert 'TWIRL_ROOT="$(realpath -e' in score
    assert 'OUT_DIR="$(realpath -m' in score
    assert '"${OUT_DIR}" != "${TWIRL_ROOT}"/*' in score

    for relative in paths[4:8]:
        wrapper = (ROOT / relative).read_text(encoding="utf-8")
        assert "#SBATCH --exclude=node4900" in wrapper
        assert "realpath -m" in wrapper
        assert '"${TWIRL_ROOT}"/*' in wrapper
        assert ".claims" in wrapper
        assert ".pending" in wrapper
        assert "producer-git-sha" in wrapper or "attest_s63" in wrapper

    launch = (ROOT / paths[9]).read_text(encoding="utf-8")
    render = (ROOT / paths[10]).read_text(encoding="utf-8")
    assert "build_s63_teacher_v3_launch_manifest.py" in launch
    assert "verify_teacher_v3_s63_review_sheets.py" in render
    assert "--branch-name current_adp" in render
    assert "--apertures DET_FLUX_ADP_SML DET_FLUX_ADP" in render
    assert "--use-row-ephemeris" in render
    assert "--no-pdf" in render
    assert "--cadence-reference" in render

    legacy_native = (ROOT / paths[3]).read_text(encoding="utf-8")
    assert "SECTOR < 56 || SECTOR > 62" in legacy_native
    assert "#SBATCH -c 1" in legacy_native


def test_full_synthetic_s63_launch_manifest(tmp_path: Path) -> None:
    module = _load_script(
        "scripts/stage5_validation/build_s63_teacher_v3_launch_manifest.py",
        "s63_full_launch_manifest_test",
    )
    accepted = tmp_path / "accepted_stage1.json"
    _write_json(accepted, _accepted_stage1())

    reserved = tmp_path / "reserved.txt"
    teacher_v3_corpus = tmp_path / "teacher_v3_corpus.csv"
    reserved.write_text("1\n2\n3\n", encoding="ascii")
    teacher_v3_corpus.write_text("tic\n2\n", encoding="ascii")
    plan = tmp_path / "prospective_plan.json"
    _write_json(
        plan,
        {
            "schema_version": "twirl_teacher_v3_s63_prospective_plan_v1",
            "status": "plan_frozen_launch_manifest_pending",
            "authorization_contract": module.AUTHORIZATION_CONTRACT,
            "sector": 63,
            "search_contract": {
                "apertures": list(module.ADP_ONLY_APERTURES),
                "anchor_aperture": module.ADP_ONLY_APERTURES[0],
                "context_aperture": module.ADP_ONLY_APERTURES[1],
                "n_periods": 50_000,
                "n_retained_peaks_per_aperture": 10,
                "candidate_peak_rank": 1,
                "orbitid_policy": "strict",
                "absolute_probability_threshold": None,
            },
            "blinded_queue": {
                "selection_seed": 630056,
                "primary": {
                    "n_tics": 1000,
                    "quotas": {
                        "planet_preserve": 300,
                        "eclipse_contact": 200,
                        "smooth_variable": 150,
                        "broad_dip": 100,
                        "disagreement_harmonic": 150,
                        "stratified_control": 100,
                    },
                },
                "repeated_host_side_cohort": {
                    "n_tics": 100,
                    "quotas": {
                        "planet_preserve": 30,
                        "eclipse_contact": 20,
                        "smooth_variable": 15,
                        "broad_dip": 10,
                        "disagreement_harmonic": 15,
                        "stratified_control": 10,
                    },
                },
            },
            "frozen_teacher_v3": {
                "automatic_promotion": False,
                "student_training_authorized": False,
            },
            "frozen_training_identity": {
                "morphology_corpus_sha256": module.file_sha256(
                    teacher_v3_corpus
                ),
                "n_corpus_tics": 1,
            },
            "s63_identity_reservation": {
                "reserved_tics_sha256": module.file_sha256(reserved),
                "n_requested_tics": 3,
            },
        },
    )
    selection_policy = (
        ROOT
        / "reports/stage5_validation/teacher_v3_s63_prospective_v1/"
        "preregistered/selection_policy_v1.json"
    )
    assert selection_policy.is_file()

    cadence = tmp_path / "cadence.csv"
    cadence.write_text("sector,orbitid,camera,ccd,cadenceno\n63,133,1,1,1\n")
    cadence_manifest = tmp_path / "cadence.manifest.json"
    _write_json(
        cadence_manifest,
        {
            "contract_version": "s63_a2v1_cadence_reference_v1",
            "sector": 63,
            "cadence_authority": "qlp_cam_quat",
            "quality_authority": "spoc_and_qlp_quality_flags",
            "orbits": [133, 134],
            "n_spoc_authority_files_verified": 16,
            "n_qlp_qflag_files_verified": 32,
            "n_rows": 1,
            "table_sha256": module.file_sha256(cadence),
        },
    )

    compact = tmp_path / "compact.h5"
    with h5py.File(compact, "w") as h5:
        h5.attrs["sector"] = 63
        h5.attrs["flux_columns"] = json.dumps(list(module.ADP_ONLY_APERTURES))
        h5.attrs["n_targets"] = 2
        targets = h5.create_group("targets")
        for tic in (1, 2):
            group = targets.create_group(f"{tic:016d}")
            group.attrs["tic"] = tic
            group.attrs["sector"] = 63
            group.attrs["camera"] = 1
            group.attrs["ccd"] = 1
            group.create_dataset("time", data=[1.0])
            group.create_dataset("cadenceno", data=[1])
            group.create_dataset("quality", data=[0])
            group.create_dataset("orbitid", data=[133])
            for aperture in module.ADP_ONLY_APERTURES:
                group.create_dataset(aperture, data=[1.0])
    compact_manifest = tmp_path / "compact.manifest.json"
    _write_json(
        compact_manifest,
        {
            "sector": 63,
            "requested_columns": list(module.ADP_ONLY_APERTURES),
            "n_exported_targets": 2,
        },
    )

    model_ready = tmp_path / "model_ready.txt"
    primary = tmp_path / "primary.txt"
    repeated = tmp_path / "repeated.txt"
    excluded = tmp_path / "excluded.txt"
    model_ready.write_text("1\n2\n", encoding="ascii")
    primary.write_text("1\n", encoding="ascii")
    repeated.write_text("2\n", encoding="ascii")
    excluded.write_text("3\n", encoding="ascii")
    model_ready_summary = tmp_path / "model_ready.summary.json"
    _write_json(
        model_ready_summary,
        {
            "schema_version": "twirl_teacher_v3_s63_model_ready_allowlist_v1",
            "sector": 63,
            "orbits": [133, 134],
            "passed": True,
            "light_curve_tensors_read": False,
            "source_hashes": {
                "prospective_plan_sha256": module.file_sha256(plan),
                "accepted_stage1_validation_sha256": module.file_sha256(
                    accepted
                ),
                "compact_h5_sha256": module.file_sha256(compact),
                "compact_manifest_sha256": module.file_sha256(compact_manifest),
                "reserved_tics_sha256": module.file_sha256(reserved),
            },
            "counts": {"n_reserved_tics": 3, "n_model_ready_tics": 2},
            "partition_checks": {
                "model_ready_subset_of_reservation": True,
                "compact_target_count_exact": True,
            },
            "outputs": {
                "model_ready_allowlist": {
                    "sha256": module.file_sha256(model_ready)
                }
            },
        },
    )
    cohort_summary = tmp_path / "cohort_summary.json"
    _write_json(
        cohort_summary,
        {
            "schema_version": "twirl_teacher_v3_s63_cohorts_v1",
            "partition_checks": {
                "model_ready_subset_of_reservation": True,
                "primary_repeated_disjoint": True,
                "primary_union_repeated_equals_model_ready": True,
                "reservation_complement_exact": True,
            },
            "source_hashes": {
                "model_ready_allowlist_sha256": module.file_sha256(model_ready),
                "teacher_v3_corpus_sha256": module.file_sha256(
                    teacher_v3_corpus
                ),
                "reserved_tics_sha256": module.file_sha256(reserved),
            },
            "outputs": {
                "primary_cohort": {"sha256": module.file_sha256(primary)},
                "repeated_host_cohort": {"sha256": module.file_sha256(repeated)},
                "reserved_not_model_ready": {"sha256": module.file_sha256(excluded)},
            },
        },
    )

    peaks = tmp_path / "bls_peaks.parquet"
    peaks.write_bytes(b"synthetic locked BLS table")
    bls_summary = tmp_path / "bls_peaks.summary.json"
    config = module.approved_a2v1_teacher_bls_config()
    _write_json(
        bls_summary,
        {
            "sector": 63,
            "contract_version": module.ADP_ONLY_CONTRACT_VERSION,
            "bls_search_contract_version": module.A2V1_TEACHER_BLS_SEARCH_CONTRACT,
            "bls_config_sha256": module.bls_config_sha256(config),
            "config": config,
            "apertures": list(module.ADP_ONLY_APERTURES),
            "n_periods": 50_000,
            "n_peaks": 10,
            "source_product_tag": "A2v1",
            "orbitid_policy": "strict",
            "n_shards": 1,
            "shard_index": 0,
            "n_source_shards": 2,
            "passed": True,
            "peak_table_sha256": module.file_sha256(peaks),
            "compact_lc_sha256": module.file_sha256(compact),
            "cadence_reference_sha256": module.file_sha256(cadence),
            "cadence_reference_manifest_sha256": module.file_sha256(cadence_manifest),
            "target_allowlist_sha256": module.file_sha256(model_ready),
            "n_targets": 2,
            "n_targets_total": 2,
            "target_allowlist_count": 2,
        },
    )

    candidates = tmp_path / "candidates.csv"
    pd.DataFrame(
        {
            "tic": [1, 2],
            "sector": [63, 63],
            "rep_aperture": [module.ADP_ONLY_APERTURES[0]] * 2,
            "rep_peak_rank": [1, 1],
            "period_d": [1.0, 2.0],
            "t0_bjd": [2459000.0, 2459001.0],
            "duration_min": [10.0, 12.0],
            "native_input_include": [True, True],
            **{
                column: [float(index + 1), float(index + 2)]
                for index, column in enumerate(
                    module.S63_TARGET_METADATA_COLUMNS
                )
            },
        }
    ).to_csv(candidates, index=False)
    candidate_summary = tmp_path / "candidates.summary.json"
    _write_json(
        candidate_summary,
        {
            "schema_version": module.CANDIDATE_CONTRACT_VERSION,
            "sector": 63,
            "science_ready": False,
            "promotion_enabled": False,
            "teacher_v3_metadata_columns": list(
                module.S63_TARGET_METADATA_COLUMNS
            ),
            "n_rows": 2,
            "n_unique_tics": 2,
            "n_model_ready_tics": 2,
            "n_bls_ineligible_tics": 0,
            "bls_ineligible_tics": [],
            "candidate_tics_sha256": module.tic_inventory_sha256({1, 2}),
            "bls_ineligible_tics_sha256": module.tic_inventory_sha256(set()),
            "teacher_scores_opened": False,
            "candidate_table_sha256": module.file_sha256(candidates),
            "artifact_hashes": {
                "prospective_plan_sha256": module.file_sha256(plan),
                "model_ready_allowlist_sha256": module.file_sha256(model_ready),
                "bls_merged_table_sha256": module.file_sha256(peaks),
            },
        },
    )

    raw = tmp_path / "raw.h5"
    with h5py.File(raw, "w") as h5:
        h5.attrs["contract_version"] = module.RAW_SOURCE_CONTRACT_VERSION
        h5.attrs["orbits"] = json.dumps([133, 134])
        targets = h5.create_group("targets")
        targets.create_group(f"{1:016d}")
        targets.create_group(f"{2:016d}")
    raw_summary = tmp_path / "raw.summary.json"
    _write_json(raw_summary, {"n_requested": 2, "n_written": 2})

    native = tmp_path / "native.h5"
    with h5py.File(native, "w") as h5:
        h5.attrs["contract_version"] = module.RAW_PAIR_CONTRACT_VERSION
        h5.attrs["training_table_sha256"] = module.file_sha256(candidates)
        h5.attrs["raw_source_h5_sha256"] = module.file_sha256(raw)
        h5.attrs["compact_adp_h5_sha256"] = module.file_sha256(compact)
        h5.attrs["cadence_reference_table_sha256"] = module.file_sha256(cadence)
        h5.attrs["cadence_reference_manifest_sha256"] = module.file_sha256(
            cadence_manifest
        )
        h5.attrs["orbitid_reconciliation_policy"] = "strict"
        h5.attrs["periodogram_n"] = 4096
        targets = h5.create_group("targets")
        targets.create_group(f"{1:016d}")
        targets.create_group(f"{2:016d}")
        h5.create_group("injections")
    native_summary = tmp_path / "native.summary.json"
    _write_json(
        native_summary,
        {
            "authorization_contract": module.AUTHORIZATION_CONTRACT,
            "sector": 63,
            "verification": {"passed": True},
            "out_h5_sha256": module.file_sha256(native),
            "exact_count_match": True,
            "exact_group_identity_match": True,
        },
    )

    repo = tmp_path / "repo"
    repo.mkdir()
    (repo / "tracked.txt").write_text("x\n")
    subprocess.run(["git", "init", "-q"], cwd=repo, check=True)
    subprocess.run(["git", "add", "tracked.txt"], cwd=repo, check=True)
    subprocess.run(
        [
            "git",
            "-c",
            "user.name=TWIRL Test",
            "-c",
            "user.email=twirl@example.invalid",
            "commit",
            "-qm",
            "test",
        ],
        cwd=repo,
        check=True,
    )
    git_sha = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    upstream_git_sha = "a" * 40

    # Freeze the immutable production receipt, then attest only its exact copy.
    hlsp_root = tmp_path / "accepted_hlsp"
    hlsp_root.mkdir()
    source_fits: dict[int, Path] = {}
    for tic in (1, 2):
        source_fits[tic] = hlsp_root / f"hlsp_twirl_s63_{tic}.fits"
        source_fits[tic].write_bytes(f"fits-{tic}".encode("ascii"))
    source_validation = tmp_path / "source_stage1_validation.json"
    source_payload = _accepted_stage1()
    source_payload["n_expected_unique_tics"] = 2
    source_payload["fits"].update(
        {
            "hlsp_root": str(hlsp_root),
            "n_fits": 2,
            "n_found_unique_tics": 2,
        }
    )
    _write_json(source_validation, source_payload)
    source_sha = module.file_sha256(source_validation)
    accepted_payload = dict(source_payload)
    accepted_payload.update(
        {
            "producer_git_sha": upstream_git_sha,
            "attestation_contract_version": (
                "twirl_teacher_v3_s63_stage1_receipt_attestation_v1"
            ),
            "source_validation_path": str(source_validation),
            "source_validation_sha256": source_sha,
            "pre_attestation_sha256": source_sha,
        }
    )
    _write_json(accepted, accepted_payload)
    accepted_sha = module.file_sha256(accepted)

    with h5py.File(compact, "r+") as h5:
        h5.attrs["hlsp_root"] = str(hlsp_root)
        h5.attrs["producer_git_sha"] = upstream_git_sha
        h5.attrs["attestation_contract_version"] = (
            "twirl_teacher_v3_s63_compact_attestation_v1"
        )
        h5.attrs["accepted_stage1_validation_sha256"] = accepted_sha
        h5.attrs["source_host_out_h5"] = str(compact)
        for tic in (1, 2):
            h5[f"targets/{tic:016d}"].attrs["source_fits"] = str(
                source_fits[tic]
            )
    compact_sha = module.file_sha256(compact)
    _write_json(
        compact_manifest,
        {
            "sector": 63,
            "hlsp_root": str(hlsp_root),
            "out_h5": str(compact),
            "out_h5_sha256": compact_sha,
            "requested_columns": list(module.ADP_ONLY_APERTURES),
            "n_discovered_files": 2,
            "n_exported_targets": 2,
            "skipped": {
                "read_failed": 0,
                "tic_filter": 0,
                "duplicate_tic": 0,
                "no_flux_columns": 0,
            },
            "records": [
                {
                    "tic": tic,
                    "sector": 63,
                    "camera": 1,
                    "ccd": 1,
                    "n_cadences": 1,
                    "flux_columns": list(module.ADP_ONLY_APERTURES),
                    "source_fits": str(source_fits[tic]),
                }
                for tic in (1, 2)
            ],
            "producer_git_sha": upstream_git_sha,
            "attestation_contract_version": (
                "twirl_teacher_v3_s63_compact_attestation_v1"
            ),
            "accepted_stage1_validation_sha256": accepted_sha,
        },
    )
    compact_manifest_sha = module.file_sha256(compact_manifest)

    cadence_payload = json.loads(cadence_manifest.read_text())
    cadence_payload["producer_git_sha"] = upstream_git_sha
    _write_json(cadence_manifest, cadence_payload)
    cadence_manifest_sha = module.file_sha256(cadence_manifest)

    model_payload = json.loads(model_ready_summary.read_text())
    model_payload["producer_git_sha"] = upstream_git_sha
    model_payload["source_hashes"].update(
        {
            "accepted_stage1_validation_sha256": accepted_sha,
            "compact_h5_sha256": compact_sha,
            "compact_manifest_sha256": compact_manifest_sha,
        }
    )
    _write_json(model_ready_summary, model_payload)
    cohort_payload = json.loads(cohort_summary.read_text())
    cohort_payload["producer_git_sha"] = upstream_git_sha
    _write_json(cohort_summary, cohort_payload)
    bls_payload = json.loads(bls_summary.read_text())
    bls_payload.update(
        {
            "producer_git_sha": upstream_git_sha,
            "compact_lc_sha256": compact_sha,
            "cadence_reference_manifest_sha256": cadence_manifest_sha,
            "target_metadata_contract_version": (
                module.S63_TARGET_METADATA_CONTRACT
            ),
            "target_metadata_columns": list(module.S63_TARGET_METADATA_COLUMNS),
            "target_metadata_finite_counts": {
                column: 2 for column in module.S63_TARGET_METADATA_COLUMNS
            },
        }
    )
    _write_json(bls_summary, bls_payload)
    candidate_payload = json.loads(candidate_summary.read_text())
    candidate_payload["producer_git_sha"] = git_sha
    _write_json(candidate_summary, candidate_payload)
    with h5py.File(raw, "r+") as h5:
        h5.attrs["producer_git_sha"] = git_sha
        h5.attrs["training_table_sha256"] = module.file_sha256(candidates)
        h5.attrs["compact_adp_h5_sha256"] = compact_sha
    raw_sha = module.file_sha256(raw)
    _write_json(
        raw_summary,
        {
            "n_requested": 2,
            "n_written": 2,
            "producer_git_sha": git_sha,
            "out_h5_sha256": raw_sha,
            "training_table_sha256": module.file_sha256(candidates),
            "compact_adp_h5_sha256": compact_sha,
        },
    )
    with h5py.File(native, "r+") as h5:
        h5.attrs["producer_git_sha"] = git_sha
        h5.attrs["raw_source_h5_sha256"] = raw_sha
        h5.attrs["compact_adp_h5_sha256"] = compact_sha
        h5.attrs["cadence_reference_manifest_sha256"] = cadence_manifest_sha
    native_sha = module.file_sha256(native)
    _write_json(
        native_summary,
        {
            "authorization_contract": module.AUTHORIZATION_CONTRACT,
            "sector": 63,
            "verification": {"passed": True},
            "out_h5_sha256": native_sha,
            "exact_count_match": True,
            "exact_group_identity_match": True,
            "producer_git_sha": git_sha,
        },
    )

    args = argparse.Namespace(
        preregistered_contract=plan,
        selection_policy=selection_policy,
        accepted_validation=accepted,
        cadence_table=cadence,
        cadence_manifest=cadence_manifest,
        compact_h5=compact,
        compact_manifest=compact_manifest,
        reserved_tics=reserved,
        teacher_v3_corpus=teacher_v3_corpus,
        model_ready_allowlist=model_ready,
        model_ready_summary=model_ready_summary,
        primary_cohort=primary,
        repeated_host_cohort=repeated,
        reserved_not_model_ready=excluded,
        cohort_summary=cohort_summary,
        bls_peaks=peaks,
        bls_summary=bls_summary,
        candidates=candidates,
        candidate_summary=candidate_summary,
        raw_source_h5=raw,
        raw_source_summary=raw_summary,
        native_h5=native,
        native_summary=native_summary,
        repo=repo,
        expected_git_sha=git_sha,
        upstream_producer_git_sha=upstream_git_sha,
    )
    manifest = module.build_launch_manifest(args)
    assert manifest["passed"] is True
    assert manifest["locked_search"] == {
        "apertures": list(module.ADP_ONLY_APERTURES),
        "n_periods": 50_000,
        "n_peaks": 10,
        "orbitid_policy": "strict",
        "native_periodogram_n": 4096,
    }
    assert manifest["checks"]["inventories"]["cohort_counts"] == {
        "primary_corpus_disjoint": 1,
        "repeated_host_side": 1,
        "reserved_not_model_ready": 1,
    }
    output = tmp_path / "launch.json"
    module._publish_json_atomic(manifest, output)
    assert json.loads(output.read_text())["passed"] is True
    with pytest.raises(FileExistsError, match="overwrite"):
        module._publish_json_atomic(manifest, output)
    assert module.file_sha256(source_validation) == source_sha

    compact_manifest_original = compact_manifest.read_text(encoding="utf-8")
    wrong_root_payload = json.loads(compact_manifest_original)
    wrong_root = tmp_path / "different_accepted_root"
    wrong_root.mkdir()
    wrong_root_payload["hlsp_root"] = str(wrong_root)
    _write_json(compact_manifest, wrong_root_payload)
    with pytest.raises(ValueError, match="hlsp_root differs"):
        module._validate_compact(
            compact,
            compact_manifest,
            accepted_validation_path=accepted,
            expected_producer_git_sha=upstream_git_sha,
        )
    compact_manifest.write_text(compact_manifest_original, encoding="utf-8")

    bls_summary_original = bls_summary.read_text(encoding="utf-8")
    wrong_producer_payload = json.loads(bls_summary_original)
    wrong_producer_payload["producer_git_sha"] = "f" * 40
    _write_json(bls_summary, wrong_producer_payload)
    with pytest.raises(ValueError, match="differs from launch"):
        module._validate_bls(
            peaks_path=peaks,
            summary_path=bls_summary,
            compact_sha256=compact_sha,
            cadence_sha256=module.file_sha256(cadence),
            cadence_manifest_sha256=cadence_manifest_sha,
            allowlist_sha256=module.file_sha256(model_ready),
            allowlist_tics=[1, 2],
            expected_producer_git_sha=upstream_git_sha,
        )
    bls_summary.write_text(bls_summary_original, encoding="utf-8")

    subset_candidates = tmp_path / "subset_candidates.csv"
    pd.read_csv(candidates).iloc[:1].to_csv(subset_candidates, index=False)
    subset_summary = tmp_path / "subset_candidates.summary.json"
    _write_json(
        subset_summary,
        {
            "schema_version": module.CANDIDATE_CONTRACT_VERSION,
            "producer_git_sha": git_sha,
            "sector": 63,
            "science_ready": False,
            "promotion_enabled": False,
            "teacher_v3_metadata_columns": list(
                module.S63_TARGET_METADATA_COLUMNS
            ),
            "n_rows": 1,
            "n_unique_tics": 1,
            "n_model_ready_tics": 2,
            "n_bls_ineligible_tics": 1,
            "bls_ineligible_tics": [2],
            "candidate_tics_sha256": module.tic_inventory_sha256({1}),
            "bls_ineligible_tics_sha256": module.tic_inventory_sha256({2}),
            "teacher_scores_opened": False,
            "candidate_table_sha256": module.file_sha256(subset_candidates),
            "artifact_hashes": {
                "prospective_plan_sha256": module.file_sha256(plan),
                "model_ready_allowlist_sha256": module.file_sha256(model_ready),
                "bls_merged_table_sha256": module.file_sha256(peaks),
            },
        },
    )
    subset_result = module._validate_candidates(
        candidates_path=subset_candidates,
        summary_path=subset_summary,
        model_ready_tics={1, 2},
        expected_hashes={
            "prospective_plan": module.file_sha256(plan),
            "bls_peaks": module.file_sha256(peaks),
            "allowlist": module.file_sha256(model_ready),
        },
        expected_producer_git_sha=git_sha,
    )
    assert subset_result["n_candidate_tics"] == 1
    assert subset_result["n_bls_ineligible_tics"] == 1

    (repo / "untracked.txt").write_text("dirty\n", encoding="utf-8")
    with pytest.raises(ValueError, match="fully clean"):
        module._git_identity(repo, git_sha)


def test_s63_review_sheet_verifier_requires_exact_once_coverage(
    tmp_path: Path,
) -> None:
    module = _load_script(
        "scripts/stage5_validation/verify_teacher_v3_s63_review_sheets.py",
        "s63_review_sheet_verifier_test",
    )
    producer = "a" * 40
    launch_sha = "b" * 64
    public_dir = tmp_path / "public"
    sheet_dir = tmp_path / "render" / "sheets"
    public_dir.mkdir()
    sheet_dir.mkdir(parents=True)
    queue = public_dir / "review_queue_1100.csv"
    rows = pd.DataFrame(
        {
            "review_id": ["s63_0001", "s63_0002"],
            "tic": [101, 102],
            "sector": [63, 63],
            "period_d": [1.25, 2.5],
            "t0_bjd": [2_457_001.0, 2_457_002.0],
            "duration_min": [8.0, 12.0],
            "twirl_vet_sheet_name": [
                "s63_0001_twirl_twoap_current_adp.png",
                "s63_0002_twirl_twoap_current_adp.png",
            ],
            "twirl_vet_sheet_pdf_name": ["", ""],
        }
    )
    rows.to_csv(queue, index=False)
    for name in rows["twirl_vet_sheet_name"]:
        (sheet_dir / name).write_bytes(b"png")
    compact = tmp_path / "compact.h5"
    cadence = tmp_path / "cadence.csv"
    cadence_manifest = tmp_path / "cadence.manifest.json"
    compact.write_bytes(b"compact")
    cadence.write_bytes(b"cadence")
    cadence_manifest.write_bytes(b"manifest")
    public_summary = public_dir / "public_summary.json"
    _write_json(
        public_summary,
        {
            "schema_version": "twirl_teacher_v3_s63_public_review_queue_v1",
            "producer_git_sha": producer,
            "n_queue_rows": 2,
            "outputs": {
                "review_queue": {
                    "path": str(queue),
                    "sha256": module.file_sha256(queue),
                }
            },
        },
    )
    bundle_marker = public_dir / "bundle_complete.json"
    _write_json(
        bundle_marker,
        {
            "schema_version": "twirl_teacher_v3_s63_queue_bundle_complete_v1",
            "sector": 63,
            "passed": True,
            "publication_complete": True,
            "visibility": "reviewer_safe_unpublished",
            "private_completion_marker_required": True,
            "bundle_id": "c" * 32,
            "artifacts": {
                "public_review_queue": {
                    "filename": queue.name,
                    "sha256": module.file_sha256(queue),
                },
                "public_summary": {
                    "filename": public_summary.name,
                    "sha256": module.file_sha256(public_summary),
                },
            },
        },
    )
    metrics = tmp_path / "render" / "twirl_vet_metrics.csv"
    metrics_frame = rows.copy()
    metrics_frame["twirl_vet_status"] = "ok"
    metrics_frame["anchor_aperture"] = "DET_FLUX_ADP_SML"
    metrics_frame["anchor_period_d"] = rows["period_d"]
    metrics_frame["anchor_t0_bjd"] = rows["t0_bjd"]
    metrics_frame["anchor_duration_min"] = rows["duration_min"]
    metrics_frame["cadence_reference_sha256"] = module.file_sha256(cadence)
    metrics_frame["cadence_reference_manifest_sha256"] = module.file_sha256(
        cadence_manifest
    )
    metrics_frame["orbitid_policy"] = "strict"
    for column in module.QUALITY_COUNT_COLUMNS:
        metrics_frame[column] = 0
    metrics_frame["quality_n_cad_total"] = 100
    metrics_frame.to_csv(metrics, index=False)
    render_summary = tmp_path / "render" / "render_summary.json"
    _write_json(
        render_summary,
        {
            "branch_name": "current_adp",
            "apertures": list(module.APERTURES),
            "n_rows": 2,
            "use_row_ephemeris": True,
            "write_pdf": False,
            "orbitid_policy": "strict",
            "producer_git_sha": producer,
            "launch_manifest_sha256": launch_sha,
            "queue_csv_sha256": module.file_sha256(queue),
            "lc_export_h5_sha256": module.file_sha256(compact),
            "external_quality": {
                "cadence_reference_table_sha256": module.file_sha256(cadence),
                "cadence_reference_manifest_sha256": module.file_sha256(
                    cadence_manifest
                ),
            },
            "status_counts": {"ok": 2},
        },
    )
    result = module.verify_review_sheets(
        queue_csv=queue,
        public_summary_path=public_summary,
        bundle_completion_marker=bundle_marker,
        metrics_csv=metrics,
        render_summary_path=render_summary,
        sheet_dir=sheet_dir,
        compact_h5=compact,
        cadence_reference=cadence,
        cadence_reference_manifest=cadence_manifest,
        expected_rows=2,
        producer_git_sha=producer,
        launch_manifest_sha256=launch_sha,
    )
    assert result["all_rows_exactly_once"] is True
    assert result["n_unique_pngs"] == 2

    duplicate_metrics = pd.concat(
        [metrics_frame.iloc[:1], metrics_frame.iloc[:1]], ignore_index=True
    )
    duplicate_metrics.to_csv(metrics, index=False)
    with pytest.raises(ValueError, match="nonempty and unique"):
        module.verify_review_sheets(
            queue_csv=queue,
            public_summary_path=public_summary,
            bundle_completion_marker=bundle_marker,
            metrics_csv=metrics,
            render_summary_path=render_summary,
            sheet_dir=sheet_dir,
            compact_h5=compact,
            cadence_reference=cadence,
            cadence_reference_manifest=cadence_manifest,
            expected_rows=2,
            producer_git_sha=producer,
            launch_manifest_sha256=launch_sha,
        )
