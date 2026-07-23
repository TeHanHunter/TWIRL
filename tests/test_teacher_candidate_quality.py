from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from twirl.lightcurves.a2v1_cadence_reference import (
    AUTHORITY_EXCLUSION_EXTERNAL_BIT,
    AUTHORITY_EXCLUSION_POLICY,
    AUTHORITY_EXCLUSION_POLICY_CONTRACT,
    CADENCE_REFERENCE_COLUMNS,
    authority_exclusions_sha256,
    file_sha256,
)
from twirl.lightcurves.external_quality import load_external_quality_reference
import twirl.vetting.teacher_candidates as teacher_candidates


def _write_reference(root: Path) -> tuple[Path, Path]:
    table = root / "cadence.csv"
    manifest = root / "cadence.json"
    frame = pd.DataFrame(
        {
            "sector": [56, 56, 56, 56],
            "orbitid": [119, 119, 120, 120],
            "camera": [1, 1, 1, 1],
            "ccd": [1, 1, 1, 1],
            "cadenceno": [100, 101, 102, 103],
            "spoc_quality": [0, 8, 0, 0],
            "qlp_quality": [0, 0, 0, 0],
            "external_quality": [0, 8, 0, 0],
        },
        columns=CADENCE_REFERENCE_COLUMNS,
    )
    frame.to_csv(table, index=False)
    sources = []
    index = 1
    for orbit in (119, 120):
        for role, values in (
            ("qlp_cam_quat", {"orbitid": orbit, "camera": 1}),
            (
                "qlp_detector_qflag",
                {"orbitid": orbit, "camera": 1, "ccd": 1},
            ),
        ):
            path = (root / f"source_{index}.txt").resolve()
            sources.append(
                {
                    "role": role,
                    "path": str(path),
                    "sha256": f"{index:064x}",
                    **values,
                }
            )
            index += 1
    for role, values in (
        ("spoc_flag_file", {"camera": 1, "ccd": 1}),
        ("spoc_quality_table", {}),
        ("spoc_quality_provenance", {}),
    ):
        path = (root / f"source_{index}.txt").resolve()
        sources.append(
            {
                "role": role,
                "path": str(path),
                "sha256": f"{index:064x}",
                **values,
            }
        )
        index += 1
    authority_exclusions = {
        "contract_version": AUTHORITY_EXCLUSION_POLICY_CONTRACT,
        "policy": AUTHORITY_EXCLUSION_POLICY,
        "external_bit": AUTHORITY_EXCLUSION_EXTERNAL_BIT,
        "n_rows": 0,
        "by_detector": {"cam1_ccd1": {"n_rows": 0, "rows": []}},
    }
    payload = {
        "contract_version": "s56_a2v1_cadence_reference_v1",
        "builder_version": "a2v1_cadence_reference_builder_v3",
        "sector": 56,
        "cadence_authority": "qlp_cam_quat",
        "quality_authority": "spoc_and_qlp_quality_flags",
        "quality_composition": {
            "external_quality": "spoc_quality | (qlp_quality << 30)",
            "qlp_quality_raw_values": [0, 1],
            "qlp_quality_external_bit": 30,
        },
        "table_sha256": file_sha256(table),
        "table_columns": list(CADENCE_REFERENCE_COLUMNS),
        "n_rows": 4,
        "detectors": ["cam1_ccd1"],
        "orbits": [119, 120],
        "n_rows_by_detector": {"cam1_ccd1": 4},
        "n_nonzero_spoc_quality": 1,
        "n_nonzero_qlp_quality": 0,
        "n_nonzero_external_quality": 1,
        "n_spoc_authority_files_verified": 1,
        "n_qlp_qflag_files_verified": 2,
        "n_spoc_rows_excluded_by_quat": 0,
        "authority_exclusions": authority_exclusions,
        "authority_exclusions_sha256": authority_exclusions_sha256(
            authority_exclusions
        ),
        "source_file_sha256": {
            source["path"]: source["sha256"] for source in sources
        },
        "sources": sources,
    }
    manifest.write_text(json.dumps(payload), encoding="utf-8")
    return table, manifest


def test_candidate_metadata_uses_effective_quality_overlay(
    tmp_path: Path, monkeypatch
) -> None:
    compact = tmp_path / "compact.h5"
    tic = 123
    with h5py.File(compact, "w") as handle:
        group = handle.create_group("targets").create_group(f"{tic:016d}")
        group.attrs.update(
            {"tic": tic, "sector": 56, "camera": 1, "ccd": 1, "tessmag": 18.0}
        )
        group.create_dataset("time", data=np.arange(4, dtype=float))
        group.create_dataset("cadenceno", data=[100, 101, 102, 103])
        group.create_dataset("orbitid", data=[119, 119, 120, 120])
        group.create_dataset("quality", data=[0, 0, 4, 0])
        group.create_dataset("DET_FLUX_ADP_SML", data=np.ones(4))
        group.create_dataset("DET_FLUX_ADP", data=np.ones(4))
    table, manifest = _write_reference(tmp_path)
    reference = load_external_quality_reference(
        table_path=table,
        manifest_path=manifest,
        sector=56,
        expected_orbits=(119, 120),
        expected_detectors=((1, 1),),
    )
    teacher_candidates._initialize_metadata_quality_reference(reference)
    captured: list[np.ndarray] = []

    def fake_measure(lc, **kwargs):
        del kwargs
        captured.append(np.asarray(lc.quality).copy())
        return {"fake_metric": 1.0}

    monkeypatch.setattr(
        teacher_candidates, "measure_two_aperture_candidate_metadata", fake_measure
    )
    rows = teacher_candidates._measure_tic(
        (tic, [{"review_id": "row-1"}], str(compact))
    )

    assert captured[0].tolist() == [0, 1, 1, 0]
    assert rows[0]["metadata_status"] == "ok"
    assert rows[0]["metadata_quality_n_cad_internal_bad"] == 1
    assert rows[0]["metadata_quality_n_cad_authority_excluded"] == 0
    assert rows[0]["metadata_quality_n_cad_external_only_bad"] == 1
    assert rows[0]["metadata_quality_n_cad_effective_bad"] == 2
