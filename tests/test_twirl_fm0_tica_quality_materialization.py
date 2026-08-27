from __future__ import annotations

from pathlib import Path

import pytest

from scripts.stage5_validation.materialize_twirl_fm0_tica_quality_pdo import (
    _serialize_rows,
    materialize,
)


def test_serialize_tica_rows_is_sorted_and_exact() -> None:
    assert _serialize_rows([(12, 32), (10, 0), (11, 4)]) == (
        b"10, 0\n11, 4\n12, 32\n"
    )


@pytest.mark.parametrize(
    "rows",
    ([], [(10, 0), (10, 4)], [(0, 0)], [(10, -1)]),
)
def test_serialize_tica_rows_rejects_invalid_authority(
    rows: list[tuple[int, int]],
) -> None:
    with pytest.raises(ValueError):
        _serialize_rows(rows)


def test_materializer_preflight_rejects_scope_before_importing_qlp(
    tmp_path: Path,
) -> None:
    with pytest.raises(ValueError, match="S67"):
        materialize(
            sector=66,
            output_dir=tmp_path / "out",
            producer_git_sha="a" * 40,
            workers=2,
        )
    with pytest.raises(ValueError, match="workers"):
        materialize(
            sector=67,
            output_dir=tmp_path / "out",
            producer_git_sha="a" * 40,
            workers=5,
        )
