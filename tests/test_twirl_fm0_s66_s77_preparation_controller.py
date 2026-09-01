from __future__ import annotations

import hashlib
import json
import os
import subprocess
from pathlib import Path

import pytest

from twirl.models.fm0 import preparation_authority
from twirl.models.fm0.preparation_authority import (
    AUTHORITY_MAP_CAMPAIGN_ID,
    AUTHORITY_MAP_SCHEMA_VERSION,
    PreparationAuthorityMap,
    admit_from_preparation_authorities,
    load_preparation_authority_map,
)
from twirl.models.fm0.registry import FM0ContractError

ROOT = Path(__file__).resolve().parents[1]
MAP = ROOT / "configs/orcd/twirl_fm0_s66_s77_full_visit_preparation_v1.json"
CONTROLLER = ROOT / "scripts/orcd/run_twirl_fm0_s66_s77_preparation_orcd.sh"
SECTORS = list(range(66, 78))


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _controller_text() -> str:
    return CONTROLLER.read_text(encoding="utf-8")


def test_authority_map_is_exactly_s66_s77_and_freezes_s66() -> None:
    payload = json.loads(MAP.read_text(encoding="utf-8"))
    assert payload["schema_version"] == AUTHORITY_MAP_SCHEMA_VERSION
    assert payload["campaign_id"] == AUTHORITY_MAP_CAMPAIGN_ID
    assert payload["scope"] == {
        "ordered_sectors": SECTORS,
        "excluded_sectors": [65],
        "blocked_sectors": [78],
        "untouched_sectors": [79, 80],
        "full_visit_six_view_shards": True,
        "label_blind_sector_wide_hdf5_quality_qa_allowed": True,
        "sealed_aperture_photometry_or_derived_shard_access_forbidden": True,
        "temporal_panel_freeze_authorized": False,
        "model_training_authorized": False,
        "h200_use_authorized": False,
    }
    authorities = payload["authorities"]
    assert sorted(map(int, authorities["source_receipts"])) == SECTORS
    assert sorted(map(int, authorities["hdf5_quality_receipts"])) == SECTORS
    assert authorities["hdf5_quality_receipts"]["66"] == {
        "path": "/orcd/data/mki_aryeh/001/twirl/reports/stage5_validation/"
        "twirl_fm0_2_hdf5_quality_admission_v1/030e3455/s0066/receipt.json",
        "sha256": "c1c4bb613e1cd61fdc9b50b09b35282488e267155feb387220a2981b71c5e985",
    }
    assert authorities["admission_policy"]["sha256"] == (
        "10536c89d2656604c97f6289cb03fc9f5aca459aa5340303324314e4905dfd76"
    )
    for name in ("admission_policy", "exclusion_ledger"):
        binding = authorities[name]
        assert _sha(ROOT / binding["repository_path"]) == binding["sha256"]
    assert payload["outputs"]["six_view_root"] == (
        "/orcd/data/mki_aryeh/001/twirl/exports/fm0/later_sector_releases/"
        "twirl_fm0_s66_s77_full_visit_v1"
    )
    authority = load_preparation_authority_map(MAP, expected_sha256=_sha(MAP))
    assert authority.hdf5_receipt_sha256[66] == (
        "c1c4bb613e1cd61fdc9b50b09b35282488e267155feb387220a2981b71c5e985"
    )


def test_authority_map_rejects_missing_s66_runtime_hash(tmp_path: Path) -> None:
    payload = json.loads(MAP.read_text(encoding="utf-8"))
    payload["authorities"]["hdf5_quality_receipts"]["66"] = {
        "path": payload["authorities"]["hdf5_quality_receipts"]["66"]["path"],
        "sha256": None,
        "required_runtime_env": "TWIRL_FM0_S66_HDF5_QUALITY_RECEIPT_SHA256",
    }
    changed = tmp_path / "authority.json"
    changed.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(FM0ContractError, match="required at runtime"):
        load_preparation_authority_map(changed, expected_sha256=_sha(changed))


def test_controller_is_noninteractive_and_dry_run_by_default() -> None:
    result = subprocess.run(
        ["bash", str(CONTROLLER), "probe"],
        cwd=ROOT,
        env={
            "PATH": "/usr/bin:/bin",
            "TWIRL_EXPECTED_GIT_SHA": "a" * 40,
        },
        check=True,
        capture_output=True,
        text=True,
    )
    assert result.stdout.startswith("+ ssh ")
    script = _controller_text()
    assert "DRY_RUN=1" in script
    for option in (
        "BatchMode=yes",
        "PasswordAuthentication=no",
        "KbdInteractiveAuthentication=no",
        "PubkeyAuthentication=no",
        "HostbasedAuthentication=no",
        "GSSAPIAuthentication=no",
        "NumberOfPasswordPrompts=0",
        "ControlMaster=no",
    ):
        assert option in script
    assert "/Users/tehan/.ssh/cm/tehan@orcd-login.mit.edu:22" in result.stdout


def test_phase_a_dry_run_has_two_chronological_afterok_chains() -> None:
    result = subprocess.run(
        ["bash", str(CONTROLLER), "submit-phase-a"],
        cwd=ROOT,
        env={
            "PATH": "/usr/bin:/bin",
            "TWIRL_EXPECTED_GIT_SHA": "a" * 40,
        },
        check=True,
        capture_output=True,
        text=True,
    )
    output = result.stdout
    assert "slurm_twirl_fm0_mission_quality_reference_cpu.sbatch" in output
    assert "slurm_twirl_fm0_later_source_inventory_cpu.sbatch" in output
    assert output.count("previous") >= 4
    assert "afterok:" in output
    assert "phase_a_quality.tsv" in output
    assert "phase_a_source.tsv" in output
    assert ".submitting" in output
    assert "seq 66 77" in output


def test_deploy_accepts_a_predeployed_standalone_exact_clone_first() -> None:
    script = _controller_text()
    deploy = script[script.index("  deploy)") : script.index("  submit-phase-a)")]
    existing = deploy.index("if [[ -e '${REMOTE_REPO}' ]]")
    fallback_source = deploy.index("test -d '${REMOTE_SOURCE}/.git'")
    detached_gate = deploy.index("${REMOTE_GATE}")
    assert existing < fallback_source < detached_gate
    assert "git -C '${REMOTE_SOURCE}' worktree add --detach" in deploy


def test_controller_requires_phase_a_freeze_before_six_view_and_cpu_admission() -> None:
    script = _controller_text()
    six_view = script[script.index("  submit-six-view)") : script.index("  submit-admission)")]
    assert "submit-phase-a-freeze)" in script
    assert "validate-phase-a" in script
    assert "sacct -nX" in script
    assert "'COMPLETED 0:0'" in script
    assert "slurm_twirl_fm0_later_six_view_cpu.sbatch" in script
    assert "slurm_twirl_fm0_later_admission_v2_cpu.sbatch" in script
    assert "phase_a_record_sha256" in script
    assert "six_view_tail_job" in script
    assert "previous_0=''" in six_view and "previous_3=''" in six_view
    assert "(sector - 66) % 4" in six_view
    assert "six_view_fast.tsv" in script
    assert "phase_a_producer_git_sha" in script
    assert "phase_a_record_sha256" in script
    assert "TWIRL_FM0_PHASE_A_PRODUCER_GIT_SHA" in script
    assert "TWIRL_FM0_PHASE_A_RECORD_SHA256" in script
    assert "admission_v2_retry1.tsv" in script
    assert "'FAILED 1:0'" in script
    assert "submit-training" not in script
    assert "slurm_twirl_fm0_1_real_train_h200.sbatch" not in script


def test_all_new_preparation_wrappers_are_cpu_only_and_fail_closed() -> None:
    wrappers = (
        "slurm_twirl_fm0_mission_quality_reference_cpu.sbatch",
        "slurm_twirl_fm0_later_source_inventory_cpu.sbatch",
        "slurm_twirl_fm0_later_six_view_cpu.sbatch",
        "slurm_twirl_fm0_phase_a_authority_freeze_cpu.sbatch",
        "slurm_twirl_fm0_later_admission_v2_cpu.sbatch",
    )
    for name in wrappers:
        script = (ROOT / "scripts/orcd" / name).read_text(encoding="utf-8")
        assert "#SBATCH -p pg_mki_aryeh" in script
        assert "#SBATCH --exclude=node4900" in script
        assert "#SBATCH --gres" not in script
        assert "set -euo pipefail" in script
        assert "TWIRL_EXPECTED_GIT_SHA" in script
        assert "status --porcelain=v1 --untracked-files=all" in script


def test_controller_shell_syntax() -> None:
    scripts = [
        CONTROLLER,
        ROOT / "scripts/orcd/slurm_twirl_fm0_phase_a_authority_freeze_cpu.sbatch",
        ROOT / "scripts/orcd/slurm_twirl_fm0_later_admission_v2_cpu.sbatch",
    ]
    subprocess.run(
        ["bash", "-n", *map(os.fspath, scripts)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )


def test_admission_keeps_phase_a_and_current_producer_revisions_distinct(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    phase_a_producer = "a" * 40
    current_producer = "b" * 40
    phase_a_digest = "c" * 64
    six_root = tmp_path / "six-view"
    rows = []
    for sector in SECTORS:
        receipt = six_root / f"s{sector:04d}" / "receipt.json"
        receipt.parent.mkdir(parents=True)
        receipt.write_text(json.dumps({"sector": sector}), encoding="utf-8")
        rows.append(
            {
                "sector": sector,
                "mission_quality_reference": {
                    "path": f"/mission/{sector}",
                    "sha256": "1" * 64,
                },
                "hdf5_quality": {
                    "path": f"/hdf5/{sector}",
                    "sha256": "2" * 64,
                },
                "source_inventory": {
                    "path": f"/source/{sector}",
                    "sha256": "3" * 64,
                },
            }
        )
    repo = tmp_path / "repo"
    authority_path = repo / "configs" / "orcd" / "authority.json"
    authority_path.parent.mkdir(parents=True)
    authority_path.write_text("{}", encoding="utf-8")
    output = tmp_path / "admission.json"
    authority = PreparationAuthorityMap(
        path=authority_path,
        sha256="4" * 64,
        payload={
            "authorities": {
                "admission_policy": {
                    "repository_path": "configs/models/policy.yaml",
                    "sha256": "5" * 64,
                },
                "exclusion_ledger": {
                    "repository_path": "configs/models/ledger.yaml",
                    "sha256": "6" * 64,
                },
            },
            "outputs": {
                "six_view_root": str(six_root),
                "admission_receipt": str(output),
            },
        },
        hdf5_receipt_sha256={},
    )
    observed: dict[str, str] = {}

    def fake_load(
        supplied: PreparationAuthorityMap,
        *,
        expected_producer_git_sha: str,
        expected_sha256: str,
    ) -> tuple[dict[str, object], Path, str]:
        assert supplied is authority
        observed["phase_a_producer"] = expected_producer_git_sha
        observed["phase_a_digest"] = expected_sha256
        return {"sectors": rows}, tmp_path / "phase-a.json", phase_a_digest

    def fake_admit(**kwargs: object) -> tuple[dict[str, object], str]:
        observed["current_producer"] = str(kwargs["producer_git_sha"])
        return {"preparation_pool_admitted": True}, "7" * 64

    monkeypatch.setattr(preparation_authority, "load_phase_a_authority_record", fake_load)
    monkeypatch.setattr(preparation_authority, "sha256_file", lambda _path: "8" * 64)
    monkeypatch.setattr(preparation_authority, "admit_preparation_pool", fake_admit)

    result, digest = admit_from_preparation_authorities(
        authority,
        phase_a_producer_git_sha=phase_a_producer,
        phase_a_record_sha256=phase_a_digest,
        producer_git_sha=current_producer,
        output_path=output,
    )
    assert result == {"preparation_pool_admitted": True}
    assert digest == "7" * 64
    assert observed == {
        "phase_a_producer": phase_a_producer,
        "phase_a_digest": phase_a_digest,
        "current_producer": current_producer,
    }
