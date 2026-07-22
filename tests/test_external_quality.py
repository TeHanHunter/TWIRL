from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from twirl.lightcurves.a2v1_cadence_reference import (
    CADENCE_REFERENCE_COLUMNS,
    file_sha256,
)
from twirl.lightcurves.external_quality import (
    EFFECTIVE_QUALITY_POLICY,
    EXTERNAL_QUALITY_POLICY_CONTRACT,
    load_external_quality_reference,
)


def _write_reference(root: Path) -> tuple[Path, Path]:
    table_path = root / "reference.csv"
    manifest_path = root / "reference.json"
    frame = pd.DataFrame(
        {
            "sector": [56, 56, 56, 56],
            "orbitid": [119, 119, 120, 120],
            "camera": [1, 1, 1, 1],
            "ccd": [1, 1, 1, 1],
            "cadenceno": [100, 101, 102, 103],
            "spoc_quality": [0, 0, 0, 8],
            "qlp_quality": [0, 1, 0, 0],
            "external_quality": [0, 1 << 30, 0, 8],
        },
        columns=CADENCE_REFERENCE_COLUMNS,
    )
    frame.to_csv(table_path, index=False)

    sources: list[dict[str, object]] = []
    digest_index = 1

    def add_source(role: str, **values: int) -> None:
        nonlocal digest_index
        path = (root / f"source_{digest_index:02d}.txt").resolve()
        digest = f"{digest_index:064x}"
        sources.append(
            {"role": role, "path": str(path), "sha256": digest, **values}
        )
        digest_index += 1

    for orbit in (119, 120):
        add_source("qlp_cam_quat", orbitid=orbit, camera=1)
        add_source("qlp_detector_qflag", orbitid=orbit, camera=1, ccd=1)
    add_source("spoc_flag_file", camera=1, ccd=1)
    add_source("spoc_quality_table")
    add_source("spoc_quality_provenance")
    source_hashes = {
        str(source["path"]): str(source["sha256"]) for source in sources
    }
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
        "table_sha256": file_sha256(table_path),
        "table_columns": list(CADENCE_REFERENCE_COLUMNS),
        "n_rows": len(frame),
        "detectors": ["cam1_ccd1"],
        "orbits": [119, 120],
        "n_rows_by_detector": {"cam1_ccd1": len(frame)},
        "n_nonzero_spoc_quality": 1,
        "n_nonzero_qlp_quality": 1,
        "n_nonzero_external_quality": 2,
        "n_spoc_authority_files_verified": 1,
        "n_qlp_qflag_files_verified": 2,
        "source_file_sha256": source_hashes,
        "sources": sources,
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return table_path, manifest_path


def test_reference_applies_internal_or_external_quality_exactly(
    tmp_path: Path,
) -> None:
    table, manifest = _write_reference(tmp_path)
    reference = load_external_quality_reference(
        table_path=table,
        manifest_path=manifest,
        sector=56,
        expected_orbits=(119, 120),
        expected_detectors=((1, 1),),
    )

    result = reference.apply(
        sector=56,
        camera=1,
        ccd=1,
        cadenceno=np.asarray([100, 101, 103]),
        orbitid=np.asarray([119, 119, 120]),
        internal_quality=np.asarray([0, 4, 0]),
        context="test target",
    )

    assert result.quality.tolist() == [0, 1, 1]
    assert result.external_quality.tolist() == [0, 1 << 30, 8]
    assert result.counts == {
        "n_cad_total": 3,
        "n_cad_internal_bad": 1,
        "n_cad_external_bad": 2,
        "n_cad_external_only_bad": 1,
        "n_cad_effective_bad": 2,
    }
    assert reference.provenance["policy_contract"] == (
        EXTERNAL_QUALITY_POLICY_CONTRACT
    )
    assert reference.provenance["effective_quality_policy"] == (
        EFFECTIVE_QUALITY_POLICY
    )


def test_reference_fails_closed_on_missing_cadence_or_orbit_mismatch(
    tmp_path: Path,
) -> None:
    table, manifest = _write_reference(tmp_path)
    reference = load_external_quality_reference(
        table_path=table, manifest_path=manifest, sector=56
    )

    with pytest.raises(ValueError, match="coverage is missing cadences"):
        reference.apply(
            sector=56,
            camera=1,
            ccd=1,
            cadenceno=[999],
            orbitid=[119],
            internal_quality=[0],
        )
    with pytest.raises(ValueError, match="orbit mapping mismatch"):
        reference.apply(
            sector=56,
            camera=1,
            ccd=1,
            cadenceno=[100],
            orbitid=[120],
            internal_quality=[0],
        )


def test_reference_rejects_legacy_quality_schema(tmp_path: Path) -> None:
    table, manifest_path = _write_reference(tmp_path)
    frame = pd.read_csv(table).rename(columns={"external_quality": "quality"})
    frame.to_csv(table, index=False)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["table_sha256"] = file_sha256(table)
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(ValueError, match="columns must be exactly"):
        load_external_quality_reference(
            table_path=table, manifest_path=manifest_path, sector=56
        )


def test_reference_detects_evidence_change_after_loading(tmp_path: Path) -> None:
    table, manifest = _write_reference(tmp_path)
    reference = load_external_quality_reference(
        table_path=table, manifest_path=manifest, sector=56
    )
    table.write_text(table.read_text(encoding="utf-8") + "\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match="table changed"):
        reference.assert_unchanged()
