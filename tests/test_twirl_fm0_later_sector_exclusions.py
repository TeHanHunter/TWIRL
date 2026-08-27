from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "configs/models/twirl_fm0_later_sector_exclusions_v1.yaml"


def test_s65_exclusion_is_pre_model_and_interval_is_self_consistent() -> None:
    payload = yaml.safe_load(LEDGER.read_text(encoding="utf-8"))

    assert payload["schema_version"] == "twirl_fm0_later_sector_exclusion_ledger_v1"
    assert payload["status"] == "frozen_pre_model_evaluation"
    assert payload["scope"]["selection_is_label_blind"] is True
    assert payload["scope"]["model_training_authorized"] is False
    assert payload["scope"]["temporal_panel_frozen"] is False

    assert len(payload["excluded_sectors"]) == 1
    exclusion = payload["excluded_sectors"][0]
    assert exclusion["sector"] == 65
    assert exclusion["reason_code"] == "incomplete_mission_quality_authority"
    assert exclusion["exclusion_is_model_outcome_based"] is False
    assert exclusion["missing_cadence_intervals_inclusive"] == [[802061, 803522]]
    start, stop = exclusion["missing_cadence_intervals_inclusive"][0]
    assert stop - start + 1 == exclusion["missing_cadence_count"] == 1462
