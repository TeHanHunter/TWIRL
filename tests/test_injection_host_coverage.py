import json

import pandas as pd

from scripts.stage5_validation.audit_s56_injection_host_coverage import build_coverage_audit


def test_build_coverage_audit_flags_partial_host_coverage(tmp_path):
    export_manifest = tmp_path / "export.manifest.json"
    injection_manifest = tmp_path / "injection_manifest.csv"
    injection_summary = tmp_path / "summary.json"

    export_manifest.write_text(
        json.dumps(
            {
                "n_discovered_files": 4,
                "n_exported_targets": 4,
                "skipped": {},
                "records": [
                    {"tic": 101, "tessmag": 14.2},
                    {"tic": 102, "tessmag": 15.2},
                    {"tic": 103, "tessmag": 16.2},
                    {"tic": 104, "tessmag": 18.2},
                ],
            }
        )
    )
    pd.DataFrame(
        [
            {"injection_id": "inj_0", "tic": 101, "tessmag": 14.2, "split": "train", "grid_cell_id": "p0_r0"},
            {"injection_id": "inj_1", "tic": 101, "tessmag": 14.2, "split": "train", "grid_cell_id": "p0_r0"},
            {"injection_id": "inj_2", "tic": 104, "tessmag": 18.2, "split": "test", "grid_cell_id": "p1_r0"},
        ]
    ).to_csv(injection_manifest, index=False)
    injection_summary.write_text(
        json.dumps({"sampling_mode": "period_radius_grid", "grid_period_bins": 2, "grid_radius_bins": 1})
    )

    result = build_coverage_audit(
        export_manifest=export_manifest,
        injection_manifest=injection_manifest,
        injection_summary=injection_summary,
        expected_total_targets=4,
        tmag_bins=(0, 15, 17, 30),
    )

    assert result["n_injections"] == 3
    assert result["n_unique_injected_tics"] == 2
    assert result["unique_injected_fraction_of_export"] == 0.5
    assert result["n_exported_tics_missing_from_injections"] == 2
    assert result["sample_exported_tics_missing_from_injections"] == [102, 103]
    assert result["n_duplicate_injected_tics"] == 1
    assert result["sample_duplicate_injected_tics"] == {"101": 2}
    assert result["is_literal_full_s56_host_coverage"] is False
    assert "injection hosts are a subset of exported S56 targets" in result["warnings"]
    assert "some exported TICs are absent from the injection hosts" in result["warnings"]
    assert "some TIC hosts have multiple injections" in result["warnings"]
    assert result["grid_summary"]["n_grid_cells_observed"] == 2
    assert result["tmag_bins"][0]["n_exported_targets"] == 1
    assert result["tmag_bins"][0]["n_injected_unique_tics"] == 1


def test_build_coverage_audit_detects_injected_tic_missing_from_export(tmp_path):
    export_manifest = tmp_path / "export.manifest.json"
    injection_manifest = tmp_path / "injection_manifest.csv"
    export_manifest.write_text(
        json.dumps(
            {
                "n_exported_targets": 1,
                "records": [{"tic": 101, "tessmag": 14.2}],
            }
        )
    )
    pd.DataFrame([{"injection_id": "inj_0", "tic": 999, "tessmag": 14.2}]).to_csv(
        injection_manifest, index=False
    )

    result = build_coverage_audit(
        export_manifest=export_manifest,
        injection_manifest=injection_manifest,
        injection_summary=None,
        expected_total_targets=1,
    )

    assert result["n_injected_tics_missing_from_export"] == 1
    assert result["sample_injected_tics_missing_from_export"] == [999]
    assert "some injected TICs are absent from the export manifest" in result["warnings"]
