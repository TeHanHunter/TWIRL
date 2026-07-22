from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import h5py
import numpy as np
import pandas as pd
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = (
    REPO_ROOT
    / "scripts"
    / "stage5_validation"
    / "build_s56_adp_real_bls_peaks.py"
)


def _load_module():
    spec = importlib.util.spec_from_file_location(
        "build_s56_adp_real_bls_peaks_test", SCRIPT_PATH
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_compact(path: Path, *, n_cadences: int = 220) -> tuple[int, np.ndarray]:
    tic = 123456789
    cadences = np.arange(1000, 1000 + n_cadences, dtype=np.int32)
    orbitid = np.where(np.arange(n_cadences) < n_cadences // 2, 119, 120)
    internal_quality = np.zeros(n_cadences, dtype=np.int32)
    internal_quality[3] = 7
    time = np.linspace(0.0, 3.0, n_cadences)
    flux = 1.0 + 1.0e-3 * np.sin(2.0 * np.pi * time)
    with h5py.File(path, "w") as h5:
        group = h5.create_group("targets").create_group(f"{tic:016d}")
        group.attrs.update(
            {
                "tic": tic,
                "tessmag": 17.5,
                "sector": 56,
                "camera": 1,
                "ccd": 1,
            }
        )
        group.create_dataset("time", data=time)
        group.create_dataset("cadenceno", data=cadences)
        group.create_dataset("orbitid", data=orbitid.astype(np.int16))
        group.create_dataset("quality", data=internal_quality)
        group.create_dataset("DET_FLUX_ADP_SML", data=flux)
        group.create_dataset("DET_FLUX_ADP", data=flux)
    return tic, cadences


def _reference_frame(cadences: np.ndarray) -> pd.DataFrame:
    n = len(cadences)
    orbitid = np.where(np.arange(n) < n // 2, 119, 120)
    spoc_quality = np.zeros(n, dtype=np.int64)
    qlp_quality = np.zeros(n, dtype=np.int64)
    spoc_quality[10] = 4
    qlp_quality[11] = 1
    external_quality = spoc_quality | (qlp_quality << 30)
    return pd.DataFrame(
        {
            "sector": np.full(n, 56, dtype=np.int64),
            "orbitid": orbitid.astype(np.int64),
            "camera": np.ones(n, dtype=np.int64),
            "ccd": np.ones(n, dtype=np.int64),
            "cadenceno": cadences.astype(np.int64),
            "spoc_quality": spoc_quality,
            "qlp_quality": qlp_quality,
            "external_quality": external_quality,
        }
    )


def _write_reference_pair(
    tmp_path: Path, frame: pd.DataFrame
) -> tuple[Path, Path]:
    table_path = tmp_path / "cadence_reference.csv"
    manifest_path = tmp_path / "cadence_reference.json"
    source_path = tmp_path / "authority_source.txt"
    source_path.write_text("fixed authority input\n", encoding="utf-8")
    frame.to_csv(table_path, index=False, lineterminator="\n")
    source_hash = _sha256(source_path)
    external_quality = frame.get("external_quality", frame.get("quality"))
    assert external_quality is not None
    manifest = {
        "contract_version": "s56_a2v1_cadence_reference_v1",
        "sector": 56,
        "cadence_authority": "qlp_cam_quat",
        "quality_authority": "spoc_and_qlp_quality_flags",
        "quality_composition": {
            "external_quality": "spoc_quality | (qlp_quality << 30)",
            "qlp_quality_raw_values": [0, 1],
            "qlp_quality_external_bit": 30,
        },
        "table_sha256": _sha256(table_path),
        "n_rows": len(frame),
        "detectors": ["cam1_ccd1"],
        "orbits": sorted(set(frame["orbitid"].astype(int))),
        "n_nonzero_spoc_quality": int(np.count_nonzero(frame["spoc_quality"])),
        "n_nonzero_qlp_quality": int(np.count_nonzero(frame["qlp_quality"])),
        "n_nonzero_external_quality": int(
            np.count_nonzero(external_quality)
        ),
        "source_file_sha256": {str(source_path): source_hash},
        "sources": [
            {
                "role": "test_authority",
                "path": str(source_path),
                "sha256": source_hash,
            }
        ],
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return table_path, manifest_path


def _fake_bls_runner(captured_quality: list[np.ndarray]):
    def run(lc, cfg, *, aperture):
        del cfg
        captured_quality.append(np.asarray(lc.quality).copy())
        n_total = len(lc.time)
        n_good = int(np.count_nonzero(lc.quality == 0))
        return SimpleNamespace(
            tic=lc.tic,
            sector=lc.sector,
            cam=lc.cam,
            ccd=lc.ccd,
            tmag=lc.tmag,
            aperture=aperture,
            n_cad_total=n_total,
            n_cad_quality=n_good,
            n_cad_kept=n_good,
            n_cad_edge_trimmed=0,
            n_cad_sigma_clipped=0,
            dropout_frac=(n_total - n_good) / n_total,
            quality_dropout_frac=(n_total - n_good) / n_total,
            n_orbits=2,
            baseline_d=3.0,
            status="ok",
            peaks=[],
        )

    return run


def _cfg_payload() -> dict[str, object]:
    return {
        "n_periods": 50,
        "n_peaks": 2,
        "p_min_d": 0.12,
        "p_max_cap_d": 15.0,
        "max_period_fraction": 0.45,
        "sigma_clip": 5.0,
        "orbit_edge_trim_d": 0.0,
        "source_product_tag": "A2v1-test",
    }


def test_external_bad_cadences_are_combined_before_bls(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module = _load_module()
    compact_path = tmp_path / "compact.h5"
    tic, cadences = _write_compact(compact_path)
    table_path, manifest_path = _write_reference_pair(
        tmp_path, _reference_frame(cadences)
    )
    reference, provenance = module.load_external_quality_reference(
        table_path=table_path, manifest_path=manifest_path, sector=56
    )
    module._initialize_external_quality_worker(
        module._reference_worker_payload(reference), provenance
    )
    captured: list[np.ndarray] = []
    monkeypatch.setattr(module, "run_bls_on_lc", _fake_bls_runner(captured))

    rows = module._process_target((tic, str(compact_path), _cfg_payload()))

    assert len(captured) == 2
    for effective_quality in captured:
        assert np.flatnonzero(effective_quality).tolist() == [3, 10, 11]
    assert {row["n_cad_internal_bad"] for row in rows} == {1}
    assert {row["n_cad_external_bad"] for row in rows} == {2}
    assert {row["n_cad_external_only_bad"] for row in rows} == {2}
    assert {row["n_cad_effective_bad"] for row in rows} == {3}
    assert {row["n_cad_quality"] for row in rows} == {len(cadences) - 3}
    assert {row["external_quality_policy_contract"] for row in rows} == {
        module.EXTERNAL_QUALITY_POLICY_CONTRACT
    }


def test_missing_external_cadence_coverage_fails_closed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module = _load_module()
    compact_path = tmp_path / "compact.h5"
    tic, cadences = _write_compact(compact_path)
    incomplete = _reference_frame(cadences).iloc[:-1].copy()
    table_path, manifest_path = _write_reference_pair(tmp_path, incomplete)
    reference, provenance = module.load_external_quality_reference(
        table_path=table_path, manifest_path=manifest_path, sector=56
    )
    module._initialize_external_quality_worker(
        module._reference_worker_payload(reference), provenance
    )
    monkeypatch.setattr(module, "run_bls_on_lc", _fake_bls_runner([]))

    with pytest.raises(ValueError, match="coverage is missing compact cadences"):
        module._process_target((tic, str(compact_path), _cfg_payload()))


def test_reference_rejects_legacy_generic_quality_column(tmp_path: Path) -> None:
    module = _load_module()
    _, cadences = _write_compact(tmp_path / "compact.h5")
    frame = _reference_frame(cadences).rename(
        columns={"external_quality": "quality"}
    )
    table_path, manifest_path = _write_reference_pair(tmp_path, frame)

    with pytest.raises(ValueError, match="columns must be exactly"):
        module.load_external_quality_reference(
            table_path=table_path, manifest_path=manifest_path, sector=56
        )


def test_summary_and_resume_bind_all_quality_inputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module = _load_module()
    compact_path = tmp_path / "compact.h5"
    _, cadences = _write_compact(compact_path)
    table_path, manifest_path = _write_reference_pair(
        tmp_path, _reference_frame(cadences)
    )
    captured: list[np.ndarray] = []
    monkeypatch.setattr(module, "run_bls_on_lc", _fake_bls_runner(captured))
    out_dir = tmp_path / "output"

    summary = module.build_peak_table(
        compact_lc=compact_path,
        cadence_reference=table_path,
        cadence_reference_manifest=manifest_path,
        out_dir=out_dir,
        workers=1,
        n_periods=50,
        n_peaks=2,
        max_targets=None,
        progress_every=0,
        resume=False,
        source_product_tag="A2v1-test",
    )

    assert summary["compact_lc_sha256"] == _sha256(compact_path)
    assert summary["cadence_reference_sha256"] == _sha256(table_path)
    assert summary["cadence_reference_manifest_sha256"] == _sha256(manifest_path)
    assert (
        summary["external_quality_policy_contract"]
        == module.EXTERNAL_QUALITY_POLICY_CONTRACT
    )
    assert summary["cadence_reference_quality_authority"] == (
        "spoc_and_qlp_quality_flags"
    )
    assert summary["peak_table_sha256"] == _sha256(
        Path(summary["outputs"]["peak_table"])
    )

    def fail_if_called(*args, **kwargs):
        del args, kwargs
        raise AssertionError("valid resume should not rerun BLS")

    monkeypatch.setattr(module, "run_bls_on_lc", fail_if_called)
    resumed = module.build_peak_table(
        compact_lc=compact_path,
        cadence_reference=table_path,
        cadence_reference_manifest=manifest_path,
        out_dir=out_dir,
        workers=1,
        n_periods=50,
        n_peaks=2,
        max_targets=None,
        progress_every=0,
        resume=True,
        source_product_tag="A2v1-test",
    )
    assert resumed == summary
