from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "check_docs.py"
SPEC = importlib.util.spec_from_file_location("check_docs", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
check_docs = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(check_docs)


def test_doc_checker_accepts_local_external_and_historical_links(tmp_path: Path) -> None:
    (tmp_path / "README.md").write_text(
        "# TWIRL\n\n[plan](doc/twirl_plan.md)\n"
        "[site](https://example.org/reference)\n"
    )
    doc = tmp_path / "doc"
    doc.mkdir()
    (doc / "twirl_plan.md").write_text("# Plan\n\n## Immediate priorities\n")
    (doc / "twirl_progress_log.md").write_text(
        "# History\n\nOld `scripts/stage1_lcs/worker.py` deployment.\n"
    )

    assert check_docs.check_repository(tmp_path) == []


def test_doc_checker_reports_broken_links_stale_paths_and_priority_drift(
    tmp_path: Path,
) -> None:
    (tmp_path / "README.md").write_text(
        "# TWIRL\n\n## Current priorities\n\n"
        "[missing](doc/missing.md)\n"
        "Use `scripts/stage1_lcs/worker.py`.\n"
    )
    doc = tmp_path / "doc"
    doc.mkdir()
    (doc / "twirl_plan.md").write_text("# Plan\n")

    errors = check_docs.check_repository(tmp_path)

    assert any("broken local link" in error for error in errors)
    assert any("stale active path" in error for error in errors)
    assert any("priority heading belongs only" in error for error in errors)


def test_doc_checker_ignores_archives_and_fenced_examples(tmp_path: Path) -> None:
    (tmp_path / "README.md").write_text(
        "# TWIRL\n\n```markdown\n[example](missing.md)\n"
        "## Current priorities\n```\n"
    )
    doc = tmp_path / "doc"
    archive = doc / "archive"
    archive.mkdir(parents=True)
    (doc / "twirl_plan.md").write_text("# Plan\n")
    (archive / "old.md").write_text(
        "## Current priorities\n[missing](missing.md)\n"
        "`scripts/stage1_lcs/worker.py`\n"
    )

    assert check_docs.check_repository(tmp_path) == []
