from __future__ import annotations

import importlib.util
import signal
import sys
from pathlib import Path

SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts/stage5_validation/build_twirl_fm0_mission_quality_reference.py"
)


def _load_script():
    spec = importlib.util.spec_from_file_location(
        "build_twirl_fm0_mission_quality_reference", SCRIPT
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_cli_forwards_exact_reference_contract(tmp_path, monkeypatch, capsys):
    module = _load_script()
    transfer = tmp_path / "transfer"
    output = tmp_path / "quality" / "s0066"
    revision = "a" * 40
    transfer_manifest_sha256 = "b" * 64
    observed = {}

    def fake_build(**kwargs):
        observed.update(kwargs)
        return {
            "schema_version": "twirl_fm0_mission_quality_reference_v2",
            "sector": 66,
            "mission_quality_provider": "spoc",
        }

    monkeypatch.setattr(module, "build_mission_quality_reference", fake_build)
    result = module.main(
        [
            "--quality-transfer-root",
            str(transfer),
            "--quality-transfer-manifest-sha256",
            transfer_manifest_sha256,
            "--sector",
            "66",
            "--output-dir",
            str(output),
            "--producer-git-sha",
            revision,
        ]
    )

    assert result == 0
    assert observed == {
        "quality_transfer_root": transfer,
        "expected_quality_transfer_manifest_sha256": transfer_manifest_sha256,
        "sector": 66,
        "output_dir": output,
        "producer_git_sha": revision,
    }
    stdout = capsys.readouterr().out
    assert '"mission_quality_provider": "spoc"' in stdout
    assert '"sector": 66' in stdout


def test_cli_signal_handler_raises_cleanup_friendly_exception():
    module = _load_script()
    try:
        module._interrupt(signal.SIGTERM, None)
    except InterruptedError as exc:
        assert "SIGTERM" in str(exc)
    else:
        raise AssertionError("SIGTERM handler did not interrupt the builder")
