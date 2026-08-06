from __future__ import annotations

import importlib.util
from pathlib import Path
import subprocess

import pandas as pd
import pytest

from twirl.vetting.teacher_v3_inference import verify_execution_checkout


ROOT = Path(__file__).resolve().parents[1]
QUEUE_CLI = (
    ROOT / "scripts" / "stage5_validation" / "build_teacher_v3_s63_review_queue.py"
)
QUEUE_WRAPPER = ROOT / "scripts" / "orcd" / "slurm_teacher_v3_s63_queue_cpu.sbatch"
COHORT_CLI = (
    ROOT / "scripts" / "stage5_validation" / "build_teacher_v3_s63_cohorts.py"
)
CANDIDATE_CLI = (
    ROOT / "scripts" / "stage5_validation" / "build_teacher_v3_s63_candidates.py"
)


def _run_git(repo: Path, *arguments: str) -> str:
    completed = subprocess.run(
        ["git", *arguments],
        cwd=repo,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def _clean_repo(path: Path) -> tuple[Path, str, dict[str, object]]:
    path.mkdir()
    _run_git(path, "init")
    _run_git(path, "config", "user.name", "TWIRL test")
    _run_git(path, "config", "user.email", "twirl-test@example.invalid")
    tracked = path / "tracked.txt"
    tracked.write_text("frozen\n", encoding="utf-8")
    _run_git(path, "add", "tracked.txt")
    _run_git(path, "commit", "-m", "frozen checkout")
    sha = _run_git(path, "rev-parse", "HEAD")
    launch_git: dict[str, object] = {
        "repo": str(path.resolve()),
        "sha": sha,
        "checkout_clean": True,
        "untracked_files_checked": True,
    }
    return path, sha, launch_git


def _queue_cli_module():
    spec = importlib.util.spec_from_file_location("s63_queue_cli_test", QUEUE_CLI)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_queue_cli_requires_git_binding_before_score_reads() -> None:
    source = QUEUE_CLI.read_text(encoding="utf-8")
    assert 'parser.add_argument("--repo", type=Path, required=True)' in source
    assert 'parser.add_argument("--expected-git-sha", required=True)' in source
    assert 'parser.add_argument("--producer-git-sha", required=True)' in source
    first_audit = source.index("git_audit = verify_execution_checkout(")
    producer_binding = source.index(
        "producer_git_sha != git_audit[\"observed_git_sha\"]"
    )
    launch_producer_binding = source.index(
        "launch_producer_git_sha != git_audit[\"launch_git_sha\"]"
    )
    scoring_read = source.index("scoring_summary = _bound_json(")
    score_read = source.index("scores = _bound_table(")
    final_audit = source.index("final_git_audit = verify_execution_checkout(")
    publication = source.index("payload = write_s63_prospective_review_queue(")
    assert (
        first_audit
        < producer_binding
        < launch_producer_binding
        < scoring_read
        < score_read
        < final_audit
        < publication
    )


def test_all_queue_precursor_clis_require_producer_git_sha() -> None:
    for path in (COHORT_CLI, CANDIDATE_CLI, QUEUE_CLI):
        source = path.read_text(encoding="utf-8")
        assert 'parser.add_argument("--producer-git-sha", required=True)' in source
        assert "validate_git_sha(" in source


@pytest.mark.parametrize(
    "column,replacement",
    (
        ("sde_max", 999.0),
        ("depth", 0.9),
        ("candidate_key", "tampered-key"),
        ("aperture_period_relation", "unrelated"),
        ("adp_sml_own_even_depth", 42.0),
        ("candidate_custom_context", "tampered"),
    ),
)
def test_candidate_score_equality_covers_every_candidate_derived_field(
    column: str,
    replacement: object,
) -> None:
    module = _queue_cli_module()
    candidates = pd.DataFrame(
        {
            "tic": [11, 12],
            "sector": [63, 63],
            "review_id": ["candidate-11", "candidate-12"],
            "candidate_key": ["key-11", "key-12"],
            "period_d": [1.0, 2.0],
            "t0_bjd": [2460000.0, 2460001.0],
            "duration_min": [10.0, 12.0],
            "depth": [0.1, 0.2],
            "depth_snr": [8.0, 9.0],
            "sde_max": [20.0, 21.0],
            "aperture_period_relation": ["exact", "harmonic"],
            "adp_sml_own_even_depth": [0.1, 0.2],
            "candidate_custom_context": ["alpha", "beta"],
        }
    )
    scores = candidates.copy()
    module._validate_candidate_score_equality(scores, candidates)
    scores.loc[0, column] = replacement
    with pytest.raises(ValueError, match=column):
        module._validate_candidate_score_equality(scores, candidates)


def test_queue_cpu_wrapper_is_hash_bound_and_output_isolated() -> None:
    subprocess.run(["bash", "-n", str(QUEUE_WRAPPER)], check=True)
    source = QUEUE_WRAPPER.read_text(encoding="utf-8")
    assert "#SBATCH -p pg_mki_aryeh" in source
    assert "#SBATCH --exclude=node4900" in source
    assert "--gres" not in source
    assert "umask 077" in source
    assert '"${REPO}/scripts/assert_clean_checkout.sh" "${REPO}"' in source
    assert 'TWIRL_ROOT="$(realpath -e -- "${TWIRL_ROOT}")"' in source
    assert 'PUBLIC_OUT_DIR="$(realpath -m -- "${PUBLIC_OUT_DIR}")"' in source
    assert 'HIDDEN_OUT_DIR="$(realpath -m -- "${HIDDEN_OUT_DIR}")"' in source
    assert '"${output}" != "${TWIRL_ROOT}"/*' in source
    assert "must be distinct and non-nested" in source
    assert "Refusing to overwrite S63 queue output" in source
    assert '--hidden-out-dir "${HIDDEN_OUT_DIR}" >/dev/null' in source
    assert 'public=${PUBLIC_OUT_DIR}"' in source
    assert "hidden=${HIDDEN_OUT_DIR}" not in source
    for environment_binding in (
        "TWIRL_EXPECTED_GIT_SHA",
        "TWIRL_S63_PROSPECTIVE_PLAN_SHA256",
        "TWIRL_S63_SELECTION_POLICY_SHA256",
        "TWIRL_S63_LAUNCH_MANIFEST_SHA256",
        "TWIRL_S63_SCORING_SUMMARY_SHA256",
        "TWIRL_S63_SCORE_TABLE_SHA256",
        "TWIRL_S63_CANDIDATES_SHA256",
        "TWIRL_S63_PRIMARY_COHORT_SHA256",
        "TWIRL_S63_REPEATED_HOST_COHORT_SHA256",
        "TWIRL_S63_COHORT_SUMMARY_SHA256",
    ):
        assert environment_binding in source
    for cli_flag in (
        "--repo",
        "--expected-git-sha",
        "--producer-git-sha",
        "--plan",
        "--expected-plan-sha256",
        "--selection-policy",
        "--expected-selection-policy-sha256",
        "--launch-manifest",
        "--expected-launch-manifest-sha256",
        "--scoring-summary",
        "--expected-scoring-summary-sha256",
        "--scores",
        "--expected-scores-sha256",
        "--candidates",
        "--expected-candidates-sha256",
        "--primary-cohort",
        "--expected-primary-cohort-sha256",
        "--repeated-host-cohort",
        "--expected-repeated-host-cohort-sha256",
        "--cohort-summary",
        "--expected-cohort-summary-sha256",
        "--public-out-dir",
        "--hidden-out-dir",
    ):
        assert cli_flag in source


@pytest.mark.parametrize("dirty_kind", ("tracked", "untracked"))
def test_queue_checkout_audit_rejects_dirty_and_untracked_files(
    tmp_path: Path,
    dirty_kind: str,
) -> None:
    repo, sha, launch_git = _clean_repo(tmp_path / "repo")
    if dirty_kind == "tracked":
        (repo / "tracked.txt").write_text("dirty\n", encoding="utf-8")
    else:
        (repo / "untracked.txt").write_text("untracked\n", encoding="utf-8")
    with pytest.raises(ValueError, match="fully clean"):
        verify_execution_checkout(
            repo=repo,
            expected_git_sha=sha,
            launch_git=launch_git,
        )


def test_queue_checkout_audit_rejects_operator_sha_drift(tmp_path: Path) -> None:
    repo, sha, launch_git = _clean_repo(tmp_path / "repo")
    wrong_sha = ("0" if sha[0] != "0" else "1") + sha[1:]
    with pytest.raises(ValueError, match="Git SHA binding differs"):
        verify_execution_checkout(
            repo=repo,
            expected_git_sha=wrong_sha,
            launch_git=launch_git,
        )
