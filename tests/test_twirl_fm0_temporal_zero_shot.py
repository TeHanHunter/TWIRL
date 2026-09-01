from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from twirl.models.fm0.temporal_panel import (
    TEMPORAL_PANEL_FIELDS,
    TEMPORAL_PANEL_READY_STATE,
    TEMPORAL_PANEL_RECEIPT_SCHEMA_VERSION,
    TEMPORAL_PANEL_SCHEMA_VERSION,
    TEMPORAL_PANEL_SECTOR_BINDING_FIELDS,
    TEMPORAL_PANEL_SECTOR_BINDING_SCHEMA_VERSION,
)
from twirl.models.fm0.temporal_zero_shot import (
    _sample_set_sha256,
    load_temporal_panel,
    paired_checkpoint_delta,
    select_temporal_rows,
)


def _selection_row(
    component: str,
    observation: str,
    cohort: str,
    sector: int,
) -> dict[str, str]:
    return {
        "leakage_component_id": component,
        "observation_key": observation,
        "temporal_cohort": cohort,
        "source_partition": "poc_development",
        "sector": str(sector),
    }


def test_bounded_selection_is_order_invariant_and_retains_all_visits() -> None:
    rows = [
        _selection_row("repeat-a", "ra-66", "repeated", 66),
        _selection_row("repeat-a", "ra-70", "repeated", 70),
        _selection_row("repeat-b", "rb-67", "repeated", 67),
        _selection_row("repeat-c", "rc-68", "repeated", 68),
        _selection_row("new-a", "na-66", "new", 66),
        _selection_row("new-b", "nb-67", "new", 67),
        _selection_row("new-b", "nb-71", "new", 71),
        _selection_row("new-c", "nc-68", "new", 68),
    ]
    forward = select_temporal_rows(
        rows, max_repeated_components=2, max_new_components=2
    )
    reverse = select_temporal_rows(
        list(reversed(rows)), max_repeated_components=2, max_new_components=2
    )

    assert forward == reverse
    for cohort in ("repeated", "new"):
        selected_components = {row["leakage_component_id"] for row in forward[cohort]}
        assert len(selected_components) == 2
        for component in selected_components:
            expected = {
                row["observation_key"]
                for row in rows
                if row["leakage_component_id"] == component
            }
            observed = {
                row["observation_key"]
                for row in forward[cohort]
                if row["leakage_component_id"] == component
            }
            assert observed == expected


def test_bounded_selection_rejects_component_crossing_cohorts() -> None:
    rows = [
        _selection_row("same", "a", "repeated", 66),
        _selection_row("same", "b", "new", 67),
        _selection_row("repeat-b", "c", "repeated", 68),
        _selection_row("new-b", "d", "new", 69),
    ]
    with pytest.raises(ValueError, match="crosses cohorts"):
        select_temporal_rows(rows, max_repeated_components=2, max_new_components=2)


def _raw_checkpoint(
    representative: np.ndarray,
    visits: np.ndarray,
    *,
    huber: float,
    ratio: float,
) -> dict[str, object]:
    embeddings = {
        "h_window": representative,
        "z_window": representative * np.asarray([1.0, 2.0]),
    }
    visit_embeddings = {
        "h_window": visits,
        "z_window": visits * np.asarray([1.0, 2.0]),
    }
    return {
        "left": embeddings,
        "right": {name: value.copy() for name, value in embeddings.items()},
        "visits": visit_embeddings,
        "reconstruction": {
            "masked_valid_target_count": 24,
            "masked_huber_mean": huber,
            "zero_prediction_masked_huber_mean": 1.0,
            "model_to_zero_baseline_ratio": ratio,
        },
    }


def test_checkpoint_delta_is_source_and_query_paired() -> None:
    step0_representatives = np.asarray(
        [[1.0, 0.0], [1.0, 0.0], [1.0, 0.0]], dtype=np.float64
    )
    step2000_representatives = np.asarray(
        [[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0]], dtype=np.float64
    )
    step0_visits = np.repeat(step0_representatives, 2, axis=0)
    step2000_visits = np.repeat(step2000_representatives, 2, axis=0)
    initial = _raw_checkpoint(step0_representatives, step0_visits, huber=0.8, ratio=0.8)
    trained = _raw_checkpoint(
        step2000_representatives, step2000_visits, huber=0.3, ratio=0.3
    )

    result = paired_checkpoint_delta(
        initial,
        trained,
        representative_components=("a", "b", "c"),
        visit_components=("a", "a", "b", "b", "c", "c"),
        visit_sectors=(66, 67, 66, 67, 66, 67),
    )

    assert result["delta_definition"] == "step2000_minus_step0"
    assert (
        result["representations"]["h_window"]["effective_rank"]["step2000_minus_step0"]
        > 0.0
    )
    retrieval = result["representations"]["h_window"][
        "query_sector_excluded_cross_visit_retrieval"
    ]
    assert retrieval["status"] == "available"
    assert retrieval["n_paired_queries"] == 6
    assert result["masked_reconstruction"][
        "step2000_minus_step0_masked_huber_mean"
    ] == pytest.approx(-0.5)
    assert result["masked_reconstruction"][
        "step2000_minus_step0_model_to_zero_baseline_ratio"
    ] == pytest.approx(-0.5)


def test_sample_hash_binds_rows_windows_and_masks() -> None:
    rows = [{"observation_key": "obs-a"}]
    sample = {
        "flux": np.zeros((2, 4), dtype=np.float32),
        "temporal_mask": np.asarray([False, True, False, True]),
    }
    first = _sample_set_sha256(rows, [sample])
    changed = {key: value.copy() for key, value in sample.items()}
    changed["temporal_mask"][0] = True
    second = _sample_set_sha256(rows, [changed])
    assert first != second


def _write_csv(path: Path, fields: tuple[str, ...], rows: list[dict[str, str]]) -> str:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_panel_loader_verifies_ready_hash_and_scientific_boundaries(
    tmp_path: Path,
) -> None:
    sector_roots: dict[int, Path] = {}
    bindings: list[dict[str, str]] = []
    for sector in range(66, 78):
        root = tmp_path / f"sector-{sector}"
        root.mkdir()
        sector_roots[sector] = root.resolve()
        bindings.append(
            {
                "binding_schema_version": (
                    TEMPORAL_PANEL_SECTOR_BINDING_SCHEMA_VERSION
                ),
                "sector": str(sector),
                "sector_release_root": str(root.resolve()),
                "sector_receipt_path": str(root / "receipt.json"),
                "sector_receipt_sha256": f"{sector - 65:064x}",
                "sector_manifest_path": str(root / "manifest.csv"),
                "sector_manifest_sha256": f"{sector - 64:064x}",
                "n_manifest_rows": "1",
                "n_development_rows": "1" if sector in {66, 67} else "0",
                "n_train_rows_excluded": "0",
                "n_sealed_rows_emitted": "0",
                "scientific_training_eligible": "false",
                "development_evaluation_eligible": "true",
            }
        )

    panel_rows: list[dict[str, str]] = []
    for sector, cohort in ((66, "repeated"), (67, "new")):
        row = {field: "" for field in TEMPORAL_PANEL_FIELDS}
        row.update(
            {
                "panel_schema_version": TEMPORAL_PANEL_SCHEMA_VERSION,
                "observation_key": f"obs-{sector}",
                "product_instance_id": f"product-{sector}",
                "physical_source_id": f"gaia_dr3:{sector}",
                "gaia_dr3_source_id": str(sector),
                "tic_id": str(sector + 1_000),
                "sector": str(sector),
                "leakage_component_id": f"component-{sector}",
                "source_partition": "poc_development",
                "temporal_cohort": cohort,
                "sector_release_root": str(sector_roots[sector]),
                "sector_receipt_sha256": f"{sector - 65:064x}",
                "sector_manifest_sha256": f"{sector - 64:064x}",
                "relative_path": f"shards/obs-{sector}.npz",
                "sha256": f"{sector:064x}",
                "input_source_sha256": "1" * 64,
                "source_row_sha256": "2" * 64,
                "n_cadences": "4096",
                "n_segments": "2",
                "view_present_json": "[1,1,1,1,1,1]",
                "input_adapter": "provider-neutral",
                "mission_quality_provider": "tica",
                "mission_quality_reference_manifest_sha256": "3" * 64,
                "hdf5_quality_receipt_sha256": "4" * 64,
                "full_visit_shard": "true",
                "model_context_length_bound": "false",
                "scientific_training_eligible": "false",
                "development_evaluation_eligible": "true",
            }
        )
        panel_rows.append(row)

    panel_dir = tmp_path / "panel"
    panel_dir.mkdir()
    panel_hash = _write_csv(panel_dir / "panel.csv", TEMPORAL_PANEL_FIELDS, panel_rows)
    binding_hash = _write_csv(
        panel_dir / "sector_bindings.csv",
        TEMPORAL_PANEL_SECTOR_BINDING_FIELDS,
        bindings,
    )
    receipt = {
        "schema_version": TEMPORAL_PANEL_RECEIPT_SCHEMA_VERSION,
        "ready_state": TEMPORAL_PANEL_READY_STATE,
        "passed": True,
        "panel_sectors": list(range(66, 78)),
        "baseline_manifest_sectors": list(range(56, 65)),
        "excluded_sectors": [65],
        "label_blind": True,
        "identity_only": True,
        "shard_payloads_opened": False,
        "temporal_panel_frozen": True,
        "development_evaluation_eligible": True,
        "scientific_training_eligible": False,
        "model_training_authorized": False,
        "sealed_test_access_authorized": False,
        "event_retention_authorized": False,
        "formal_model_gate_authorized": False,
        "production_model_claim": False,
        "prospective_test_claim": False,
        "panel_path": "panel.csv",
        "panel_sha256": panel_hash,
        "sector_bindings_path": "sector_bindings.csv",
        "sector_bindings_sha256": binding_hash,
        "n_panel_rows": 2,
        "n_panel_components": 2,
        "n_repeated_rows": 1,
        "n_repeated_components": 1,
        "n_new_rows": 1,
        "n_new_components": 1,
        "n_sealed_rows_emitted": 0,
        "baseline_manifest_path": str(tmp_path / "manifest.csv"),
        "baseline_manifest_sha256": "5" * 64,
        "config_path": str(tmp_path / "config.yaml"),
        "config_sha256": "6" * 64,
        "admission_receipt_path": str(tmp_path / "admission.json"),
        "admission_receipt_sha256": "7" * 64,
    }
    receipt_path = panel_dir / "receipt.json"
    receipt_path.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    receipt_hash = hashlib.sha256(receipt_path.read_bytes()).hexdigest()
    (panel_dir / "READY").write_text(receipt_hash + "\n", encoding="utf-8")
    for path in panel_dir.iterdir():
        path.chmod(0o444)
    panel_dir.chmod(0o555)

    loaded, summary = load_temporal_panel(panel_dir, receipt_sha256=receipt_hash)
    assert len(loaded) == 2
    assert summary["receipt_sha256"] == receipt_hash
    assert summary["sealed_rows"] == 0

    (panel_dir / "READY").chmod(0o644)
    (panel_dir / "READY").write_text("0" * 64 + "\n", encoding="utf-8")
    (panel_dir / "READY").chmod(0o444)
    with pytest.raises(ValueError, match="READY"):
        load_temporal_panel(panel_dir, receipt_sha256=receipt_hash)
    panel_dir.chmod(0o755)
