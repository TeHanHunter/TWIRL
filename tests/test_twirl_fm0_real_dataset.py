from __future__ import annotations

import csv
from dataclasses import asdict
import hashlib
from pathlib import Path

import numpy as np
import pytest

from twirl.models.fm0.dataset import (
    FM0ReleaseDataset,
    FM0ReleaseDatasetConfig,
    deterministic_window_start,
)
from twirl.models.fm0.input_release import (
    INPUT_RELEASE_SCHEMA_VERSION,
    MANIFEST_COLUMNS,
    ObservationRelease,
    deterministic_npz_bytes,
)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _release(n: int = 3000) -> ObservationRelease:
    flux = np.zeros((n, 6), dtype=np.float32)
    flux[:, 2] = np.arange(n, dtype=np.float32) / 10_000
    flux[:, 3] = flux[:, 2] + 0.01
    valid = np.ones((n, 6), dtype=bool)
    errors = np.full((n, 2), 0.002, dtype=np.float32)
    return ObservationRelease(
        flux=flux,
        flux_valid=valid,
        flux_error=errors,
        error_valid=np.ones((n, 2), dtype=bool),
        local_time_cadences=np.arange(n, dtype=np.float32),
        delta_time_cadences=np.r_[np.float32(0), np.ones(n - 1, np.float32)],
        time_valid=np.ones(n, dtype=bool),
        segment_boundary=np.r_[True, np.zeros(n - 1, dtype=bool)],
        segment_id=np.zeros(n, dtype=np.int32),
        view_present=np.ones(6, dtype=bool),
    )


def _write_release(
    root: Path,
    *,
    partition: str = "poc_train",
    release: ObservationRelease | None = None,
) -> str:
    release = _release() if release is None else release
    shard = root / "shards" / "observation_test.npz"
    shard.parent.mkdir(parents=True)
    shard.write_bytes(deterministic_npz_bytes(release))
    row = {
        "input_release_schema_version": INPUT_RELEASE_SCHEMA_VERSION,
        "observation_key": "observation_test",
        "product_instance_id": "product_test",
        "source_sha256": "a" * 64,
        "leakage_component_id": "leakage_test",
        "source_partition": partition,
        "product_state": "A2V1_ACCEPTED",
        "relative_path": "shards/observation_test.npz",
        "sha256": _sha(shard),
        "input_source_sha256": "a" * 64,
        "n_cadences": str(release.flux.shape[0]),
        "n_segments": str(np.unique(release.segment_id).size),
        "view_present_json": "[1,1,1,1,1,1]",
        "host_visit_offset_cadences": "0.0",
        "host_visit_gap_cadences": "0.0",
        "host_visit_overlaps_previous": "False",
        "input_adapter": "a2v1_hdf5_quality_aware_v1",
        "scientific_training_eligible": "True",
    }
    manifest = root / "manifest.csv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)
    return _sha(manifest)


def test_release_dataset_is_deterministic_source_first_and_model_safe(tmp_path: Path) -> None:
    manifest_sha = _write_release(tmp_path)
    config = FM0ReleaseDatasetConfig(
        release_root=str(tmp_path),
        manifest_sha256=manifest_sha,
        variant="TWIRL-FM0.1.1",
        windows_per_epoch=8,
    )
    dataset = FM0ReleaseDataset(config)
    first = dataset.sample(3, mask_view=0)
    repeated = dataset.sample(3, mask_view=0)
    second_mask = dataset.sample(3, mask_view=1)

    assert dataset.contract["n_sources"] == 1
    assert dataset.contract["source_partition"] == "poc_train"
    assert set(first) == {
        "flux", "flux_valid", "flux_error", "error_valid",
        "local_time_cadences", "delta_time_cadences", "time_valid",
        "segment_boundary", "view_present", "temporal_mask",
        "reconstruction_mask",
    }
    assert first["flux"].shape == (2, 2048)
    assert np.array_equal(first["flux"], repeated["flux"])
    assert np.array_equal(first["reconstruction_mask"], repeated["reconstruction_mask"])
    for name in (
        "flux",
        "flux_valid",
        "flux_error",
        "error_valid",
        "local_time_cadences",
        "delta_time_cadences",
        "time_valid",
        "segment_boundary",
        "view_present",
    ):
        assert np.array_equal(first[name], second_mask[name]), name
    assert not np.array_equal(first["temporal_mask"], second_mask["temporal_mask"])
    assert not np.array_equal(first["reconstruction_mask"], second_mask["reconstruction_mask"])
    for sample in (first, second_mask):
        assert np.array_equal(
            sample["reconstruction_mask"],
            sample["flux_valid"] & sample["temporal_mask"][None, :],
        )
    assert first["local_time_cadences"][np.flatnonzero(first["time_valid"])[0]] == 0


def test_release_dataset_rejects_manifest_or_partition_drift(tmp_path: Path) -> None:
    manifest_sha = _write_release(tmp_path)
    with pytest.raises(ValueError, match="manifest differs"):
        FM0ReleaseDataset(
            FM0ReleaseDatasetConfig(
                release_root=str(tmp_path),
                manifest_sha256="0" * 64,
                variant="TWIRL-FM0.1.1",
            )
        )
    with pytest.raises(ValueError, match="no rows"):
        FM0ReleaseDataset(
            FM0ReleaseDatasetConfig(
                release_root=str(tmp_path),
                manifest_sha256=manifest_sha,
                variant="TWIRL-FM0.1.1",
                source_partition="poc_development",
            )
        )


def test_release_dataset_excludes_missing_variant_views_without_dropping_healthy_visits(
    tmp_path: Path,
) -> None:
    manifest_sha = _write_release(tmp_path)
    manifest = tmp_path / "manifest.csv"
    with manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    missing = _release()
    missing.flux[:, (0, 2, 4)] = 0.0
    missing.flux_valid[:, (0, 2, 4)] = False
    missing.view_present[[0, 2, 4]] = False
    missing_shard = tmp_path / "shards" / "observation_missing_1x1.npz"
    missing_shard.write_bytes(deterministic_npz_bytes(missing))
    missing_row = dict(rows[0])
    missing_row.update(
        {
            "observation_key": "observation_missing_1x1",
            "relative_path": "shards/observation_missing_1x1.npz",
            "sha256": _sha(missing_shard),
            "view_present_json": "[0,1,0,1,0,1]",
        }
    )
    rows.append(missing_row)
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    manifest_sha = _sha(manifest)

    dataset = FM0ReleaseDataset(
        FM0ReleaseDatasetConfig(
            release_root=str(tmp_path),
            manifest_sha256=manifest_sha,
            variant="TWIRL-FM0.1.1",
            windows_per_epoch=2,
        )
    )

    assert dataset.contract["n_sources"] == 1
    assert dataset.contract["n_observations"] == 1
    assert dataset.contract["n_excluded_missing_required_views"] == 1
    assert dataset.sample(0)["flux"].shape == (2, 2048)


def test_release_dataset_replaces_an_ineligible_window_with_a_valid_one(
    tmp_path: Path,
) -> None:
    release = _release()
    # This observation has healthy ADP data, but its early quality-masked
    # stretch is longer than a training window.  The ordinary deterministic
    # start below lands wholly inside it; the dataset must retain the host and
    # select an eligible window later in the same segment.
    release.flux[:2500, (2, 3)] = 0.0
    release.flux_valid[:2500, (2, 3)] = False
    manifest_sha = _write_release(tmp_path, release=release)
    config = FM0ReleaseDatasetConfig(
        release_root=str(tmp_path),
        manifest_sha256=manifest_sha,
        variant="TWIRL-FM0.1.1",
        windows_per_epoch=1024,
    )
    bad_index = next(
        index
        for index in range(1024)
        if deterministic_window_start(
            3000, sample_index=index, epoch=0, seed=config.seed
        ) <= 452
    )
    dataset = FM0ReleaseDataset(config)

    sample = dataset.sample(bad_index)

    assert np.all(np.any(sample["flux_valid"], axis=1))


def test_one_step_real_training_and_strict_validation(tmp_path: Path) -> None:
    torch = pytest.importorskip("torch")
    from twirl.models.fm0.model import (
        FM0ModelConfig,
        build_fm0_model,
        count_trainable_parameters,
    )
    from twirl.models.fm0.training import FM0OptimizationConfig, run_real_training
    from twirl.models.fm0.validation import (
        REAL_RUN_CONTRACT_SCHEMA_VERSION,
        REAL_RUN_SUMMARY_SCHEMA_VERSION,
        sha256_file,
        validate_real_run_release,
        write_json_with_sha256,
        write_sha256_sidecar,
    )

    release_root = tmp_path / "release"
    manifest_sha = _write_release(release_root)
    dataset = FM0ReleaseDataset(
        FM0ReleaseDatasetConfig(
            release_root=str(release_root),
            manifest_sha256=manifest_sha,
            variant="TWIRL-FM0.1.1",
            windows_per_epoch=2,
            window_length=64,
        )
    )
    receipt = tmp_path / "input.receipt.json"
    write_json_with_sha256(
        receipt,
        {
            "schema_version": "twirl_fm0_1_input_release_receipt_v1",
            "passed": True,
            "scientific_training_eligible": True,
            "release": {"manifest_sha256": manifest_sha},
        },
    )
    optimization = FM0OptimizationConfig(
        warmup_steps=0,
        max_optimizer_steps=1,
        effective_batch_windows=2,
    )
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    contract = {
        "schema_version": REAL_RUN_CONTRACT_SCHEMA_VERSION,
        "variant": "TWIRL-FM0.1.1",
        "architecture": "tcn",
        "seed": 560067,
        "synthetic_smoke": False,
        "real_data_consumed": True,
        "target_step": 1,
        "optimization": asdict(optimization),
        "input_release": {
            "receipt_path": str(receipt.resolve()),
            "receipt_sha256": sha256_file(receipt),
            "manifest_path": str((release_root / "manifest.csv").resolve()),
            "manifest_sha256": manifest_sha,
        },
    }
    contract_sha = write_json_with_sha256(run_dir / "run_contract.json", contract)
    model = build_fm0_model(
        "TWIRL-FM0.1.1",
        config_override=FM0ModelConfig(
            architecture="tcn",
            n_flux_views=2,
            window_length=64,
            d_model=16,
            embedding_dim=16,
            stem_kernel=5,
            dropout=0.0,
            tcn_blocks=2,
            tcn_dilation_cycle=(1, 2),
            conformer_heads=4,
            conformer_conv_kernel=7,
            minimum_parameters=1,
            maximum_parameters=10_000_000,
        ),
        enforce_parameter_budget=False,
    )
    result = run_real_training(
        model=model,
        dataset=dataset,
        output_dir=run_dir,
        run_contract=contract,
        optimization=optimization,
        target_step=1,
        micro_batch_windows=1,
        device="cpu",
        precision="fp32",
    )
    checkpoint_sha = write_sha256_sidecar(run_dir / "checkpoint.pt")
    write_json_with_sha256(
        run_dir / "summary.json",
        {
            "schema_version": REAL_RUN_SUMMARY_SCHEMA_VERSION,
            "passed": True,
            "synthetic_only": False,
            "real_data_consumed": True,
            "variant": "TWIRL-FM0.1.1",
            "architecture": "tcn",
            "parameter_count": count_trainable_parameters(model),
            "global_step": 1,
            "final_metrics": result["final_metrics"],
            "run_contract_sha256": contract_sha,
            "checkpoint_sha256": checkpoint_sha,
        },
    )
    validated = validate_real_run_release(run_dir)
    assert validated["passed"] is True
    assert validated["real_data_consumed"] is True
    checkpoint = torch.load(run_dir / "checkpoint.pt", weights_only=False)
    assert checkpoint["dataset_contract"]["kind"] == "fm0_input_release"
