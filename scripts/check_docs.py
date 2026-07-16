#!/usr/bin/env python3
"""Check active TWIRL documentation links and authority boundaries."""
from __future__ import annotations

import argparse
import re
from pathlib import Path
from urllib.parse import unquote, urlsplit


ACTIVE_ROOTS = ("README.md", "AGENTS.md", "doc", "catalogs", "scripts")
EXCLUDED_PARTS = {"archive", "__pycache__"}
HISTORICAL_DOCS = {Path("doc/twirl_progress_log.md")}
AUTHORITATIVE_PLAN = Path("doc/twirl_plan.md")
STALE_STAGE1_PATH = "scripts/stage1_lcs"

INLINE_LINK_RE = re.compile(
    r"!?\[[^\]]*\]\(\s*(?P<target><[^>]+>|[^\s)]+)"
    r"(?:\s+(?:\"[^\"]*\"|'[^']*'))?\s*\)"
)
REFERENCE_LINK_RE = re.compile(
    r"^\s{0,3}\[[^\]]+\]:\s*(?P<target><[^>]+>|\S+)",
    re.MULTILINE,
)
PRIORITY_HEADING_RE = re.compile(
    r"^#{1,6}\s+(?:(?:current|immediate|top|near-term)\s+)"
    r"(?:(?:implementation|project)\s+)?priorit(?:y|ies)"
    r"(?:\s+\([^)]*\))?\s*$",
    re.IGNORECASE | re.MULTILINE,
)


def active_markdown_files(repo_root: Path) -> list[Path]:
    """Return live Markdown sources, excluding archived and generated reports."""
    paths: set[Path] = set()
    for raw_root in ACTIVE_ROOTS:
        root = repo_root / raw_root
        if root.is_file() and root.suffix.lower() == ".md":
            paths.add(root)
            continue
        if not root.is_dir():
            continue
        for path in root.rglob("*.md"):
            relative = path.relative_to(repo_root)
            if any(part in EXCLUDED_PARTS for part in relative.parts):
                continue
            paths.add(path)
    return sorted(paths)


def _without_fenced_code(text: str) -> str:
    """Blank fenced-code content while retaining source line numbers."""
    kept: list[str] = []
    fence: str | None = None
    for line in text.splitlines(keepends=True):
        stripped = line.lstrip()
        marker = stripped[:3]
        if fence is None and marker in {"```", "~~~"}:
            fence = marker
            kept.append("\n" if line.endswith("\n") else "")
        elif fence is not None:
            if stripped.startswith(fence):
                fence = None
            kept.append("\n" if line.endswith("\n") else "")
        else:
            kept.append(line)
    return "".join(kept)


def _line_number(text: str, offset: int) -> int:
    return text.count("\n", 0, offset) + 1


def _local_link_path(source: Path, raw_target: str) -> Path | None:
    target = raw_target.strip()
    if target.startswith("<") and target.endswith(">"):
        target = target[1:-1]
    target = unquote(target)
    if not target or target.startswith("#"):
        return None
    parsed = urlsplit(target)
    if parsed.scheme or parsed.netloc:
        return None
    target_path = parsed.path
    if not target_path or target_path.startswith(("/", "~")):
        return None
    return (source.parent / target_path).resolve()


def _link_errors(repo_root: Path, source: Path, text: str) -> list[str]:
    errors: list[str] = []
    matches = [*INLINE_LINK_RE.finditer(text), *REFERENCE_LINK_RE.finditer(text)]
    for match in sorted(matches, key=lambda item: item.start()):
        target = match.group("target")
        local_path = _local_link_path(source, target)
        if local_path is None or local_path.exists():
            continue
        line = _line_number(text, match.start())
        relative = source.relative_to(repo_root)
        errors.append(f"{relative}:{line}: broken local link: {target}")
    return errors


def check_repository(repo_root: Path) -> list[str]:
    """Return human-readable documentation-policy violations."""
    repo_root = repo_root.resolve()
    errors: list[str] = []
    files = active_markdown_files(repo_root)
    if not files:
        return ["no active Markdown files found"]

    for source in files:
        relative = source.relative_to(repo_root)
        raw_text = source.read_text(encoding="utf-8")
        prose = _without_fenced_code(raw_text)
        errors.extend(_link_errors(repo_root, source, prose))

        if relative not in HISTORICAL_DOCS:
            for line_number, line in enumerate(raw_text.splitlines(), start=1):
                if STALE_STAGE1_PATH in line:
                    errors.append(
                        f"{relative}:{line_number}: stale active path: "
                        f"{STALE_STAGE1_PATH}"
                    )

        if relative != AUTHORITATIVE_PLAN:
            for match in PRIORITY_HEADING_RE.finditer(prose):
                line = _line_number(prose, match.start())
                errors.append(
                    f"{relative}:{line}: priority heading belongs only in "
                    f"{AUTHORITATIVE_PLAN}"
                )

    return sorted(set(errors))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="repository root (default: inferred from this script)",
    )
    args = parser.parse_args(argv)

    errors = check_repository(args.root)
    if errors:
        print(f"documentation checks failed ({len(errors)} issue(s)):")
        for error in errors:
            print(f"- {error}")
        return 1
    print("documentation checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
