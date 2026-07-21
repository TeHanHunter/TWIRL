from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).parents[1]
REPRODUCTION_RUNNER = (
    ROOT / "scripts" / "stage1_lightcurves" / "run_a2v1_reproduction_pdo.sh"
)
QUEUE_RUNNER = (
    ROOT / "scripts" / "stage1_lightcurves" / "run_a2v1_sector_queue_pdo.sh"
)
S64_S69_MANIFEST = ROOT / "configs" / "a2v1_production_s64_s69.txt"
CHECKOUT_GUARD = ROOT / "scripts" / "assert_clean_checkout.sh"
TEACHER_RUNNER = (
    ROOT / "scripts" / "stage5_validation" / "run_s56_a2v1_teacher_search_pdo.sh"
)


def _write_executable(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    path.chmod(0o755)


def _write_clean_repo(tmp_path: Path) -> Path:
    repo = tmp_path / "clean-repo"
    assert_script = repo / "scripts" / "assert_clean_checkout.sh"
    assert_script.parent.mkdir(parents=True)
    shutil.copy2(ROOT / "scripts" / "assert_clean_checkout.sh", assert_script)
    assert_script.chmod(0o755)
    subprocess.run(["git", "init", "-q"], cwd=repo, check=True)
    subprocess.run(["git", "add", "."], cwd=repo, check=True)
    subprocess.run(
        [
            "git",
            "-c",
            "user.name=TWIRL Test",
            "-c",
            "user.email=twirl-test@example.invalid",
            "commit",
            "-qm",
            "Initialize clean test checkout",
        ],
        cwd=repo,
        check=True,
    )
    return repo


def test_a2v1_reproduction_stops_when_overlay_setup_fails(tmp_path: Path) -> None:
    fork = tmp_path / "tglc-fork"
    light_curves = fork / "tglc" / "scripts" / "light_curves.py"
    epsfs = fork / "tglc" / "scripts" / "epsfs.py"
    light_curves.parent.mkdir(parents=True)
    light_curves.write_text("# _source_tic_overlay_path\n")
    epsfs.write_text("# source.mask.mask\n")

    failing_python = tmp_path / "failing-python"
    marker = tmp_path / "tglc-was-run"
    tglc_python = tmp_path / "tglc-python"
    clean_repo = _write_clean_repo(tmp_path)
    _write_executable(failing_python, "#!/bin/sh\nexit 37\n")
    _write_executable(tglc_python, '#!/bin/sh\ntouch "$TWIRL_TEST_TGLC_MARKER"\n')

    env = os.environ.copy()
    env.update(
        {
            "TWIRL_REPO": str(clean_repo),
            "TWIRL_SOURCE_ROOT": str(tmp_path / "source"),
            "TWIRL_A2V1_ROOT": str(tmp_path / "a2v1"),
            "TWIRL_TGLC_FORK": str(fork),
            "TWIRL_QLP_PY": str(failing_python),
            "TWIRL_TGLC_PY": str(tglc_python),
            "TWIRL_SCRIPT_DIR": str(ROOT / "scripts" / "stage1_lightcurves"),
            "TWIRL_TEST_TGLC_MARKER": str(marker),
        }
    )
    result = subprocess.run(
        ["bash", str(REPRODUCTION_RUNNER), "56", "119:o1"],
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode != 0
    assert "orbit-119 source_tic overlay failed" in (result.stdout + result.stderr)
    assert not marker.exists()


def test_a2v1_queue_and_teacher_share_validation_report_contract() -> None:
    queue = QUEUE_RUNNER.read_text()
    teacher = TEACHER_RUNNER.read_text()

    assert "set -euo pipefail" in REPRODUCTION_RUNNER.read_text()
    assert (
        "OBSERVATIONS=${TWIRL_OBSERVATIONS:-$REPO/data_local/catalogs/"
        "twirl_master_catalog/twirl_wd_tess_observations_v0.fits}"
    ) in queue
    assert queue.count('--observations "$OBSERVATIONS"') == 2
    assert "s${sector}_A2v1_validation_full_schema.json" in queue
    assert "${SECTOR_SHORT}_A2v1_validation_full_schema.json" in teacher
    assert "s$(printf '%04d' \"$sector\")_A2v1_full_validation.json" not in queue


def test_checkout_guard_uses_legacy_git_compatible_directory_change() -> None:
    guard = CHECKOUT_GUARD.read_text()

    assert 'cd "${REPO}"' in guard
    assert "git -C" not in guard
    assert "--porcelain=v1" not in guard


def test_a2v1_queue_accepts_only_complete_or_absent_legacy_epsfs() -> None:
    queue = QUEUE_RUNNER.read_text()

    assert 'if [ "$kind" = "epsf" ] && [ "${#directories[@]}" -eq 0 ]; then' in queue
    assert 'if [ "$epsf_count" -eq 3136 ]; then' in queue
    assert 'elif [ "$epsf_count" -eq 0 ]; then' in queue
    assert "epsf_mode=refit-all" in queue
    assert "partial old-ePSF preparation ($epsf_count/3136)" in queue


def test_s64_s69_manifest_uses_consecutive_orbit_pairs() -> None:
    rows = [
        line
        for line in S64_S69_MANIFEST.read_text().splitlines()
        if line and not line.startswith("#")
    ]

    assert rows == [
        "64 135:o1 136:o2",
        "65 137:o1 138:o2",
        "66 139:o1 140:o2",
        "67 141:o1 142:o2",
        "68 143:o1 144:o2",
        "69 145:o1 146:o2",
    ]
