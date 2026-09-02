from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

import twirl.models.fm0.composite_release as composite
from twirl.models.fm0.input_release import (
    INPUT_RELEASE_SCHEMA_VERSION,
    MANIFEST_COLUMNS,
    ObservationRelease,
    deterministic_npz_bytes,
)
from twirl.models.fm0.later_sector_admission_v2 import POOL_RECEIPT_SCHEMA_VERSION
from twirl.models.fm0.later_sector_release import (
    LATER_SIX_VIEW_MANIFEST_FIELDS,
    LATER_SIX_VIEW_MANIFEST_SCHEMA_VERSION,
)
from twirl.models.fm0.registry import FM0ContractError

PRODUCER = "1" * 40


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _identity(prefix: str, value: str) -> str:
    return f"{prefix}_{hashlib.sha256(value.encode()).hexdigest()}"


def _release() -> ObservationRelease:
    n = 160
    flux = np.zeros((n, 6), dtype=np.float32)
    flux[:, 2] = np.linspace(-0.01, 0.01, n, dtype=np.float32)
    flux[:, 3] = flux[:, 2] + 0.002
    return ObservationRelease(
        flux=flux,
        flux_valid=np.ones((n, 6), dtype=bool),
        flux_error=np.full((n, 2), 0.002, dtype=np.float32),
        error_valid=np.ones((n, 2), dtype=bool),
        local_time_cadences=np.arange(n, dtype=np.float32),
        delta_time_cadences=np.r_[np.float32(0), np.ones(n - 1, np.float32)],
        time_valid=np.ones(n, dtype=bool),
        segment_boundary=np.r_[True, np.zeros(n - 1, dtype=bool)],
        segment_id=np.zeros(n, dtype=np.int32),
        view_present=np.ones(6, dtype=bool),
    )


def _write_csv(path: Path, fields: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_shard(root: Path, observation: str) -> tuple[str, str]:
    relative = f"shards/{observation}.npz"
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(deterministic_npz_bytes(_release()))
    path.chmod(0o444)
    return relative, _sha(path)


def _legacy_row(
    root: Path,
    *,
    name: str,
    component: str,
    partition: str,
    materialize_shard: bool,
) -> dict[str, str]:
    observation = _identity("observation", name)
    if materialize_shard:
        relative, shard_sha = _write_shard(root, observation)
    else:
        relative = f"shards/{observation}.npz"
        shard_sha = "a" * 64
    return {
        "input_release_schema_version": INPUT_RELEASE_SCHEMA_VERSION,
        "observation_key": observation,
        "product_instance_id": _identity("product", name),
        "source_sha256": "b" * 64,
        "leakage_component_id": component,
        "source_partition": partition,
        "product_state": "A2V1_ACCEPTED",
        "relative_path": relative,
        "sha256": shard_sha,
        "input_source_sha256": "b" * 64,
        "n_cadences": "160",
        "n_segments": "1",
        "view_present_json": "[1,1,1,1,1,1]",
        "host_visit_offset_cadences": "0.0",
        "host_visit_gap_cadences": "0.0",
        "host_visit_overlaps_previous": "False",
        "input_adapter": "a2v1_hdf5_quality_aware_v1",
        "scientific_training_eligible": "True",
    }


def _later_row(
    root: Path,
    *,
    sector: int,
    component: str,
    partition: str = "poc_train",
) -> dict[str, str]:
    observation = _identity("observation", f"s{sector}")
    relative, shard_sha = _write_shard(root, observation)
    return {
        "manifest_schema_version": LATER_SIX_VIEW_MANIFEST_SCHEMA_VERSION,
        "observation_key": observation,
        "product_instance_id": _identity("product", f"s{sector}"),
        "physical_source_id": f"gaia_dr3:{sector}",
        "gaia_dr3_source_id": str(sector),
        "tic_id": str(10_000 + sector),
        "sector": str(sector),
        "leakage_component_id": component,
        "source_partition": partition,
        "product_state": "ORCD_COMPLETE_DEFERRED",
        "relative_path": relative,
        "sha256": shard_sha,
        "input_source_sha256": "c" * 64,
        "source_row_sha256": hashlib.sha256(f"source-{sector}".encode()).hexdigest(),
        "n_cadences": "160",
        "n_segments": "1",
        "view_present_json": "[1,1,1,1,1,1]",
        "input_adapter": "a2v1_hdf5_mission_quality_v1",
        "mission_quality_provider": "spoc" if sector == 66 else "tica",
        "mission_quality_reference_manifest_sha256": "d" * 64,
        "hdf5_quality_receipt_sha256": "e" * 64,
        "full_visit_shard": "true",
        "model_context_length_bound": "false",
        "scientific_training_eligible": "false",
        "panel_admission_authorized": "false",
    }


def _fixture(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> tuple[
    composite.SourceManifestBinding,
    tuple[composite.SourceManifestBinding, ...],
    Path,
    str,
]:
    holdout_component = _identity("leakage", "holdout")
    legacy = tmp_path / "legacy"
    legacy_rows = [
        _legacy_row(
            legacy,
            name="legacy-overlap",
            component=holdout_component,
            partition="poc_train",
            materialize_shard=False,
        ),
        _legacy_row(
            legacy,
            name="legacy-train",
            component=_identity("leakage", "legacy-train"),
            partition="poc_train",
            materialize_shard=True,
        ),
        _legacy_row(
            legacy,
            name="legacy-development",
            component=_identity("leakage", "legacy-development"),
            partition="poc_development",
            materialize_shard=False,
        ),
        _legacy_row(
            legacy,
            name="legacy-sealed",
            component=_identity("leakage", "legacy-sealed"),
            partition="poc_sealed_test",
            materialize_shard=False,
        ),
    ]
    _write_csv(legacy / "manifest.csv", MANIFEST_COLUMNS, legacy_rows)
    legacy_receipt = tmp_path / "legacy.receipt.json"
    legacy_receipt.write_text(
        json.dumps(
            {
                "schema_version": "twirl_fm0_1_input_release_receipt_v1",
                "passed": True,
                "scientific_training_eligible": True,
                "release": {"manifest_sha256": _sha(legacy / "manifest.csv")},
            }
        )
        + "\n",
        encoding="utf-8",
    )
    legacy_binding = composite.SourceManifestBinding(
        source_kind="legacy",
        sector_min=56,
        sector_max=64,
        release_root=legacy,
        manifest_sha256=_sha(legacy / "manifest.csv"),
        receipt_path=legacy_receipt,
        receipt_sha256=_sha(legacy_receipt),
    )

    later_bindings: list[composite.SourceManifestBinding] = []
    later_rows_by_root: dict[Path, list[dict[str, str]]] = {}
    for sector in range(66, 78):
        root = tmp_path / f"s{sector}"
        component = (
            holdout_component
            if sector in {66, 77}
            else _identity("leakage", f"s{sector}")
        )
        rows = [_later_row(root, sector=sector, component=component)]
        _write_csv(root / "manifest.csv", LATER_SIX_VIEW_MANIFEST_FIELDS, rows)
        receipt = root / "receipt.json"
        receipt.write_text("{}\n", encoding="utf-8")
        (root / "READY").write_text(_sha(receipt) + "\n", encoding="utf-8")
        later_rows_by_root[root.resolve()] = rows
        later_bindings.append(
            composite.SourceManifestBinding(
                source_kind="later",
                sector_min=sector,
                sector_max=sector,
                release_root=root,
                manifest_sha256=_sha(root / "manifest.csv"),
                receipt_path=receipt,
                receipt_sha256=_sha(receipt),
            )
        )

    def fake_validate(
        root: str | Path,
        *,
        expected_receipt_sha256: str | None,
        require_read_only: bool,
        verify_shard_payloads: bool,
    ):
        root = Path(root).resolve()
        assert verify_shard_payloads is False
        assert expected_receipt_sha256 == _sha(root / "receipt.json")
        rows = later_rows_by_root[root]
        sector = int(rows[0]["sector"])
        return (
            root / "receipt.json",
            {"manifest_sha256": _sha(root / "manifest.csv")},
            {"sector": sector, "manifest_rows": tuple(rows)},
        )

    monkeypatch.setattr(composite, "validate_later_sector_release", fake_validate)
    admission = tmp_path / "admission.json"
    admission.write_text(
        json.dumps(
            {
                "schema_version": POOL_RECEIPT_SCHEMA_VERSION,
                "preparation_pool_sectors": list(range(66, 78)),
                "excluded_sectors": [65],
                "preparation_pool_admitted": True,
                "n_nonsealed_preparation_rows": 12,
                "n_six_view_shards": 12,
                "scientific_training_eligible": False,
                "model_training_authorized": False,
                "sealed_test_access_authorized": False,
                "sector_receipts": [
                    {
                        "sector": item.sector_min,
                        "evidence": {
                            "six_view_release": {
                                "sha256": item.receipt_sha256,
                                "manifest_sha256": item.manifest_sha256,
                            }
                        },
                    }
                    for item in later_bindings
                ],
            },
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return (
        legacy_binding,
        tuple(later_bindings),
        admission,
        _sha(admission),
    )


def _freeze(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> composite.CompositeReleaseResult:
    legacy, later, admission, admission_sha = _fixture(tmp_path, monkeypatch)
    return composite.freeze_composite_release(
        tmp_path / "composite",
        legacy_binding=legacy,
        later_bindings=later,
        admission_receipt_path=admission,
        admission_receipt_sha256=admission_sha,
        producer_git_sha=PRODUCER,
        require_source_read_only=False,
    )


def test_freezer_quarantines_all_earlier_s77_components_without_touching_shards(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    result = _freeze(tmp_path, monkeypatch)
    selection = result.receipt["selection"]

    assert selection["role_counts"] == {
        composite.TRAIN_ROLE: 11,
        composite.HOLDOUT_ROLE: 1,
        composite.EXCLUDED_OVERLAP_ROLE: 2,
    }
    assert selection["pre_quarantine_training_rows"] == 13
    assert selection["n_holdout_components"] == 1
    assert selection["n_excluded_overlap_components"] == 1
    assert selection["role_cadences"] == {
        composite.TRAIN_ROLE: 1760,
        composite.HOLDOUT_ROLE: 160,
        composite.EXCLUDED_OVERLAP_ROLE: 320,
    }
    assert result.receipt["sources"]["excluded_sectors"] == [65]
    assert result.receipt["limits"]["source_shards_opened"] is False
    assert result.receipt["limits"]["shards_rewritten"] is False
    assert result.receipt["limits"]["development_rows_selected"] == 0
    assert result.receipt["limits"]["sealed_rows_selected"] == 0

    with (result.root / "role_index.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert {row["role"] for row in rows} == set(composite.ROLE_ORDER)
    excluded = [row for row in rows if row["role"] == composite.EXCLUDED_OVERLAP_ROLE]
    assert {row["source_id"] for row in excluded} == {
        "legacy_s56_s64",
        "later_s0066",
    }
    assert "development" not in (result.root / "role_index.csv").read_text()
    assert "sealed" not in (result.root / "role_index.csv").read_text()


def test_composite_dataset_resolves_original_train_and_holdout_shards(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    result = _freeze(tmp_path, monkeypatch)
    common = {
        "composite_root": str(result.root),
        "receipt_sha256": result.receipt_sha256,
        "source_bindings_sha256": result.source_bindings_sha256,
        "role_index_sha256": result.role_index_sha256,
        "variant": "TWIRL-FM0.3.1",
        "windows_per_epoch": 4,
        "require_read_only": False,
    }
    train = composite.FM0CompositeReleaseDataset(
        composite.FM0CompositeDatasetConfig(**common)
    )
    holdout = composite.FM0CompositeReleaseDataset(
        composite.FM0CompositeDatasetConfig(**common, role=composite.HOLDOUT_ROLE)
    )

    assert train.contract["n_observations"] == 11
    assert holdout.contract["n_observations"] == 1
    for key in (
        "receipt_sha256",
        "source_bindings_sha256",
        "role_index_sha256",
    ):
        assert train.contract[key] == common[key]
    assert train.sample(0)["flux"].shape == (2, 128)
    assert holdout.sample(0)["flux"].shape == (2, 128)


def test_composite_validation_fails_on_role_index_drift(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    result = _freeze(tmp_path, monkeypatch)
    role_index = result.root / "role_index.csv"
    role_index.chmod(0o644)
    role_index.write_bytes(role_index.read_bytes() + b"\n")
    role_index.chmod(0o444)

    with pytest.raises(FM0ContractError, match="role index hash drifted"):
        composite.validate_composite_release(
            result.root,
            expected_receipt_sha256=result.receipt_sha256,
            expected_source_bindings_sha256=result.source_bindings_sha256,
            expected_role_index_sha256=result.role_index_sha256,
            require_source_read_only=False,
        )


def test_composite_dataset_rejects_non_native_fm03_geometry() -> None:
    with pytest.raises(ValueError, match="L128"):
        composite.FM0CompositeDatasetConfig(
            composite_root="/unused",
            receipt_sha256="a" * 64,
            source_bindings_sha256="b" * 64,
            role_index_sha256="c" * 64,
            variant="TWIRL-FM0.3.2",
            window_length=256,
        )


def _uninitialized_composite_dataset() -> composite.FM0CompositeReleaseDataset:
    dataset = object.__new__(composite.FM0CompositeReleaseDataset)
    dataset.config = composite.FM0CompositeDatasetConfig(
        composite_root="/unused",
        receipt_sha256="a" * 64,
        source_bindings_sha256="b" * 64,
        role_index_sha256="c" * 64,
        variant="TWIRL-FM0.3.1",
        windows_per_epoch=1,
        require_read_only=False,
    )
    dataset.epoch = 0
    return dataset


def test_composite_sampler_rejects_segments_shorter_than_window() -> None:
    dataset = _uninitialized_composite_dataset()
    release = _release()
    short_segment = np.arange(dataset.config.window_length - 1)

    assert (
        dataset._window_start_if_eligible(
            release,
            short_segment,
            sample_index=0,
        )
        is None
    )
    # The inherited FM0.1 sampler deliberately retains its padded-window
    # behavior; only the FM0.3 composite subclass rejects the short segment.
    assert (
        composite.FM0ReleaseDataset._window_start_if_eligible(
            dataset,
            release,
            short_segment,
            sample_index=0,
        )
        == 0
    )
    assert (
        dataset._window_start_if_eligible(
            release,
            np.arange(dataset.config.window_length),
            sample_index=0,
        )
        == 0
    )


def test_composite_selection_fails_if_parent_returns_a_short_slice(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    dataset = _uninitialized_composite_dataset()
    release = _release()
    short_indices = np.arange(dataset.config.window_length - 1)
    row = {"observation_key": _identity("observation", "short")}
    monkeypatch.setattr(
        composite.FM0ReleaseDataset,
        "_selection",
        lambda self, index: (row, release, short_indices, 0),
    )

    with pytest.raises(ValueError, match="shorter than window_length"):
        dataset._selection(0)
