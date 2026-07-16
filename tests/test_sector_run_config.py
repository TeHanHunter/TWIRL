from dataclasses import asdict
from pathlib import Path

import pytest

from twirl.search.bls import BLSConfig
from twirl.search.sector_run import (
    _build_arg_parser,
    _config_from_args,
    resolve_hlsp_root,
)


def test_sector_run_cli_defaults_match_bls_config() -> None:
    args = _build_arg_parser().parse_args(
        ["--sector", "56", "--hlsp-root", "/accepted/s56"]
    )

    assert asdict(_config_from_args(args)) == asdict(BLSConfig())
    assert resolve_hlsp_root(56, args.hlsp_root) == Path("/accepted/s56")


def test_sector_run_yaml_matches_bls_config() -> None:
    config_path = Path(__file__).parents[1] / "configs" / "detection" / "bls_default.yaml"
    args = _build_arg_parser().parse_args(
        [
            "--sector", "56",
            "--hlsp-root", "/accepted/s56",
            "--config", str(config_path),
        ]
    )

    assert asdict(_config_from_args(args)) == asdict(BLSConfig())


def test_sector_run_requires_explicit_hlsp_root() -> None:
    with pytest.raises(SystemExit):
        _build_arg_parser().parse_args(["--sector", "56"])

    with pytest.raises(ValueError, match="explicit --hlsp-root"):
        resolve_hlsp_root(56, None)
