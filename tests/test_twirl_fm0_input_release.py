from __future__ import annotations

import csv
import hashlib

import numpy as np
import pytest

from twirl.models.fm0.input_release import (
    FLUX_VIEW_NAMES,
    FM0ContractError,
    build_observation_release,
    deterministic_npz_bytes,
    deterministic_training_window,
    evaluation_windows,
    extract_window,
    load_input_release,
    validate_model_tensors,
    write_input_release,
)
from twirl.models.fm0.registry import (
    build_alias_registry,
    build_observation_registry,
    write_registry_release,
)


_RAW_MANIFEST_FIELDS = (
    "observation_key",
    "product_instance_id",
    "source_sha256",
    "raw_npz_path",
    "raw_npz_sha256",
)


def _write_raw_manifest(path, rows) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_RAW_MANIFEST_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def _raw_fixture(n: int = 420) -> dict[str, np.ndarray]:
    cadence = np.arange(700000, 700000 + n, dtype=np.int64)
    time = np.arange(n, dtype=np.float64) * 200.0 / 86400.0 + 2900.0
    split = n // 2
    time[split:] += 1.0
    orbit = np.where(np.arange(n) < split, 119, 120).astype(np.int64)
    phase = np.arange(n) / 43.0
    rng = np.random.default_rng(560067)
    raw_small = 1000.0 + 8.0 * np.sin(phase) + rng.normal(0, 1.0, n)
    raw_large = 2400.0 + 14.0 * np.sin(phase + 0.2) + rng.normal(0, 1.5, n)
    internal = np.zeros(n, dtype=np.uint64)
    spoc = np.zeros(n, dtype=np.uint64)
    qlp = np.zeros(n, dtype=np.uint64)
    excluded = np.zeros(n, dtype=bool)
    internal[4] = 1
    spoc[5] = 2
    qlp[6] = 1
    excluded[7] = True
    time[8] = np.nan
    error_small = np.full(n, 1.0)
    error_large = np.full(n, 1.5)
    error_small[9] = 0.0
    return {
        "time": time,
        "cadence": cadence,
        "orbit": orbit,
        "internal_quality": internal,
        "spoc_quality": spoc,
        "qlp_quality": qlp,
        "authority_excluded": excluded,
        "raw_flux_1x1": raw_small,
        "raw_flux_error_1x1": error_small,
        "raw_flux_3x3": raw_large,
        "raw_flux_error_3x3": error_large,
    }


def test_six_view_quality_error_time_and_window_contract(tmp_path) -> None:
    release = build_observation_release(_raw_fixture())
    assert release.flux.shape == (420, 6)
    assert len(FLUX_VIEW_NAMES) == 6
    assert release.flux_error.shape == (420, 2)
    assert not release.flux_valid[4:9].any(axis=1).any()
    assert not release.error_valid[9, 0]
    assert release.flux_valid[9, 0]
    assert release.segment_boundary[0]
    assert release.segment_boundary[210]
    assert release.n_segments == 2
    assert release.audit["external_quality_formula"] == "spoc_quality | (qlp_quality << 30)"

    specs = evaluation_windows(release)
    assert [(spec.segment_id, spec.start_offset) for spec in specs] == [(0, 0), (1, 0)]
    training = deterministic_training_window(
        release, observation_key="observation_test", epoch=2, draw_index=3
    )
    assert training == deterministic_training_window(
        release, observation_key="observation_test", epoch=2, draw_index=3
    )
    tensors = extract_window(
        release, segment_id=training.segment_id, start_offset=training.start_offset
    )
    assert tensors["flux"].shape == (2048, 6)
    assert tensors["local_time_cadences"][np.flatnonzero(tensors["time_valid"])[0]] == 0
    assert not tensors["time_valid"][training.n_observed :].any()
    with pytest.raises(FM0ContractError, match="allowlist"):
        validate_model_tensors(dict(tensors, sector=np.array(56)))

    payload = deterministic_npz_bytes(release)
    assert payload == deterministic_npz_bytes(release)
    shard = tmp_path / "shard.npz"
    shard.write_bytes(payload)
    loaded = load_input_release(shard)
    assert np.array_equal(loaded.flux, release.flux)
    assert np.array_equal(loaded.local_time_cadences, release.local_time_cadences)


def test_raw_array_allowlist_rejects_bls() -> None:
    raw = _raw_fixture(120)
    raw["bls_period"] = np.ones(120)
    with pytest.raises(FM0ContractError, match="keys mismatch"):
        build_observation_release(raw)


def test_absolute_raw_time_must_be_monotonic_on_finite_cadences() -> None:
    raw = _raw_fixture(120)
    raw["time"][20] = raw["time"][19] - 0.01
    with pytest.raises(FM0ContractError, match="strictly increasing"):
        build_observation_release(raw)


def test_partial_input_release_writer_is_hash_bound_not_certification(tmp_path) -> None:
    aliases = build_alias_registry(
        [{"gaia_dr3_source_id": "123456789012345678", "tic_id": "42"}]
    )
    observations = build_observation_registry(
        [
            {
                "gaia_dr3_source_id": "123456789012345678",
                "tic_id": "42",
                "sector": 56,
                "a2v1_product_version": "A2v1",
                "source_sha256": "a" * 64,
                "product_state": "A2V1_ACCEPTED",
            }
        ],
        aliases,
    )
    registry_dir = tmp_path / "registry"
    write_registry_release(registry_dir, aliases, observations)

    raw_path = tmp_path / "raw.npz"
    np.savez(raw_path, **_raw_fixture(160))
    raw_hash = hashlib.sha256(raw_path.read_bytes()).hexdigest()
    raw_manifest = tmp_path / "raw_manifest.csv"
    with raw_manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "observation_key",
                "product_instance_id",
                "source_sha256",
                "raw_npz_path",
                "raw_npz_sha256",
            ),
        )
        writer.writeheader()
        writer.writerow(
            {
                "observation_key": observations[0]["observation_key"],
                "product_instance_id": observations[0]["product_instance_id"],
                "source_sha256": observations[0]["source_sha256"],
                "raw_npz_path": raw_path,
                "raw_npz_sha256": raw_hash,
            }
        )
    output = tmp_path / "release"
    summary = write_input_release(
        registry_dir=registry_dir,
        raw_manifest_path=raw_manifest,
        out_dir=output,
    )
    assert summary["n_observations"] == 1
    assert not summary["certifies_full_campaign"]
    assert not summary["scientific_training_eligible"]
    assert summary["input_adapter"] == "strict_npz_fixture_v1"
    assert (output / "manifest.csv").is_file()
    rows = list(csv.DictReader((output / "manifest.csv").open(newline="")))
    assert rows[0]["product_instance_id"] == observations[0]["product_instance_id"]
    assert rows[0]["source_sha256"] == observations[0]["source_sha256"]
    assert rows[0]["scientific_training_eligible"] == "False"
    assert float(rows[0]["host_visit_offset_cadences"]) == 0.0
    assert float(rows[0]["host_visit_gap_cadences"]) == 0.0
    assert write_input_release(
        registry_dir=registry_dir,
        raw_manifest_path=raw_manifest,
        out_dir=output,
    ) == summary


@pytest.mark.parametrize(
    ("field", "wrong_value", "message"),
    (
        ("product_instance_id", "product_wrong", "product_instance_id"),
        ("source_sha256", "b" * 64, "source_sha256"),
    ),
)
def test_raw_manifest_must_bind_registry_product_and_source(
    tmp_path, field: str, wrong_value: str, message: str
) -> None:
    aliases = build_alias_registry(
        [{"gaia_dr3_source_id": "123456789012345678", "tic_id": "42"}]
    )
    observations = build_observation_registry(
        [
            {
                "gaia_dr3_source_id": "123456789012345678",
                "tic_id": "42",
                "sector": 56,
                "a2v1_product_version": "A2v1",
                "source_sha256": "a" * 64,
                "product_state": "A2V1_ACCEPTED",
            }
        ],
        aliases,
    )
    registry_dir = tmp_path / "registry"
    write_registry_release(registry_dir, aliases, observations)
    raw_path = tmp_path / "raw.npz"
    np.savez(raw_path, **_raw_fixture(120))
    row = {
        "observation_key": observations[0]["observation_key"],
        "product_instance_id": observations[0]["product_instance_id"],
        "source_sha256": observations[0]["source_sha256"],
        "raw_npz_path": raw_path,
        "raw_npz_sha256": hashlib.sha256(raw_path.read_bytes()).hexdigest(),
    }
    row[field] = wrong_value
    raw_manifest = tmp_path / "raw_manifest.csv"
    _write_raw_manifest(raw_manifest, [row])
    with pytest.raises(FM0ContractError, match=message):
        write_input_release(
            registry_dir=registry_dir,
            raw_manifest_path=raw_manifest,
            out_dir=tmp_path / "release",
        )


def test_identical_raw_hashes_cannot_cross_leakage_components(tmp_path) -> None:
    aliases = build_alias_registry(
        [
            {"gaia_dr3_source_id": "100", "tic_id": "10"},
            {"gaia_dr3_source_id": "200", "tic_id": "20"},
        ]
    )
    observations = build_observation_registry(
        [
            {
                "gaia_dr3_source_id": gaia,
                "tic_id": tic,
                "sector": 56,
                "a2v1_product_version": "A2v1",
                "source_sha256": source_hash * 64,
                "product_state": "A2V1_ACCEPTED",
            }
            for gaia, tic, source_hash in (("100", "10", "a"), ("200", "20", "b"))
        ],
        aliases,
    )
    registry_dir = tmp_path / "registry"
    write_registry_release(registry_dir, aliases, observations)
    raw_path = tmp_path / "raw.npz"
    np.savez(raw_path, **_raw_fixture(120))
    raw_hash = hashlib.sha256(raw_path.read_bytes()).hexdigest()
    rows = [
        {
            "observation_key": observation["observation_key"],
            "product_instance_id": observation["product_instance_id"],
            "source_sha256": observation["source_sha256"],
            "raw_npz_path": raw_path,
            "raw_npz_sha256": raw_hash,
        }
        for observation in observations
    ]
    raw_manifest = tmp_path / "raw_manifest.csv"
    _write_raw_manifest(raw_manifest, rows)
    with pytest.raises(FM0ContractError, match="different leakage components"):
        write_input_release(
            registry_dir=registry_dir,
            raw_manifest_path=raw_manifest,
            out_dir=tmp_path / "release",
        )


def test_host_visit_offset_and_gap_are_derived_from_absolute_raw_time(tmp_path) -> None:
    aliases = build_alias_registry(
        [{"gaia_dr3_source_id": "123456789012345678", "tic_id": "42"}]
    )
    observations = build_observation_registry(
        [
            {
                "gaia_dr3_source_id": "123456789012345678",
                "tic_id": "42",
                "sector": sector,
                "a2v1_product_version": "A2v1",
                "source_sha256": source_hash * 64,
                "product_state": "A2V1_ACCEPTED",
            }
            for sector, source_hash in ((56, "a"), (57, "b"))
        ],
        aliases,
    )
    observations_by_sector = {int(row["sector"]): row for row in observations}
    registry_dir = tmp_path / "registry"
    write_registry_release(registry_dir, aliases, observations)

    raw_56 = _raw_fixture(160)
    raw_57 = _raw_fixture(160)
    raw_57["time"] = raw_57["time"] + 30.0
    paths = {56: tmp_path / "s56.npz", 57: tmp_path / "s57.npz"}
    np.savez(paths[56], **raw_56)
    np.savez(paths[57], **raw_57)
    rows = []
    for sector in (57, 56):  # manifest order is deliberately not chronological
        observation = observations_by_sector[sector]
        rows.append(
            {
                "observation_key": observation["observation_key"],
                "product_instance_id": observation["product_instance_id"],
                "source_sha256": observation["source_sha256"],
                "raw_npz_path": paths[sector],
                "raw_npz_sha256": hashlib.sha256(paths[sector].read_bytes()).hexdigest(),
            }
        )
    raw_manifest = tmp_path / "raw_manifest.csv"
    _write_raw_manifest(raw_manifest, rows)
    output = tmp_path / "release"
    write_input_release(
        registry_dir=registry_dir,
        raw_manifest_path=raw_manifest,
        out_dir=output,
    )
    manifest = {
        row["observation_key"]: row
        for row in csv.DictReader((output / "manifest.csv").open(newline=""))
    }
    first = manifest[observations_by_sector[56]["observation_key"]]
    second = manifest[observations_by_sector[57]["observation_key"]]
    finite_56 = raw_56["time"][np.isfinite(raw_56["time"])]
    finite_57 = raw_57["time"][np.isfinite(raw_57["time"])]
    units_per_day = 86400.0 / 200.0
    assert float(first["host_visit_offset_cadences"]) == 0.0
    assert float(first["host_visit_gap_cadences"]) == 0.0
    assert float(second["host_visit_offset_cadences"]) == pytest.approx(
        (finite_57[0] - finite_56[0]) * units_per_day
    )
    assert float(second["host_visit_gap_cadences"]) == pytest.approx(
        (finite_57[0] - finite_56[-1]) * units_per_day
    )
