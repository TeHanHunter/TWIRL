"""Release-bound orbit-ID reconciliation contracts for native model inputs."""
from __future__ import annotations

from typing import Any, Mapping


S62_TEACHER_V3_ORBITID_RECONCILIATION_RELEASE = (
    "teacher_v3_s56_s62_a2v1_current_adp_s62_boundary_v1"
)
S62_TEACHER_V3_ORBITID_RECONCILIATION_BINDING: Mapping[str, Any] = {
    "training_table_sha256": (
        "47497f068adec26c8a209f47ab17d9601c06411fc83bd32a7fb6fdbbdb9b5422"
    ),
    "raw_source_h5_sha256": (
        "6ce62385794f4cf5fb7b446b76504159c0c03f6a492b900d239fe379420a8563"
    ),
    "compact_adp_h5_sha256": (
        "2c76f33a2f3b104df7ac34e01f6f275984d426a126d7b51faccfe4af739191f2"
    ),
    "cadence_reference_table_sha256": (
        "ee13371211ee3065cde5173da783cacd407e5142212b4cda267e8c4941408cb5"
    ),
    "cadence_reference_manifest_sha256": (
        "72a1419d5d4ed865ff6afc0af6bc407c5aaf42d3fb3435980cb45554877c0ce0"
    ),
    "cadence_reference_source_declaration_sha256": (
        "b1b39374d0a131c678e7000f428cb326207979f1ed3a68996586d13ce6367ac1"
    ),
    "n_real_targets": 997,
    "n_groups_corrected": 996,
    "n_groups_unmodified": 1,
    "n_cad_corrected": 65_478,
    "min_cadenceno_corrected": 766_048,
    "max_cadenceno_corrected": 766_136,
    "input_orbitid": 132,
    "reference_orbitid": 131,
    "n_cad_corrected_by_camera": {
        "cam1": 4_095,
        "cam2": 10_920,
        "cam3": 23_585,
        "cam4": 26_878,
    },
    "n_groups_corrected_by_camera": {
        "cam1": 117,
        "cam2": 312,
        "cam3": 265,
        "cam4": 302,
    },
    "n_cad_corrected_by_detector": {
        "cam1_ccd1": 1_015,
        "cam1_ccd2": 1_260,
        "cam1_ccd3": 1_120,
        "cam1_ccd4": 700,
        "cam2_ccd1": 2_345,
        "cam2_ccd2": 4_025,
        "cam2_ccd3": 2_485,
        "cam2_ccd4": 2_065,
        "cam3_ccd1": 8_455,
        "cam3_ccd2": 5_162,
        "cam3_ccd3": 2_047,
        "cam3_ccd4": 7_921,
        "cam4_ccd1": 3_471,
        "cam4_ccd2": 2_848,
        "cam4_ccd3": 6_942,
        "cam4_ccd4": 13_617,
    },
}


__all__ = [
    "S62_TEACHER_V3_ORBITID_RECONCILIATION_BINDING",
    "S62_TEACHER_V3_ORBITID_RECONCILIATION_RELEASE",
]
