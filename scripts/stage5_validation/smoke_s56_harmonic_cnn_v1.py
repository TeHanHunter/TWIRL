#!/usr/bin/env python3
"""Run a deterministic CUDA smoke of the complete harmonic-teacher trainer."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from twirl.lightcurves.a2v1_cadence_reference import (  # noqa: E402
    AUTHORITY_EXCLUSION_EXTERNAL_BIT,
    AUTHORITY_EXCLUSION_POLICY_CONTRACT,
)
from twirl.lightcurves.external_quality import (  # noqa: E402
    EFFECTIVE_QUALITY_POLICY,
    EXTERNAL_QUALITY_POLICY_CONTRACT,
)
from twirl.vetting.adjudication_audit import HARMONIC_CNN_TARGET_POLICY  # noqa: E402
from twirl.vetting.harmonic_cnn import HarmonicTrainConfig  # noqa: E402
from twirl.vetting.harmonic_inputs import (  # noqa: E402
    CHANNEL_CONTRACT,
    NATIVE_DATASETS,
    PERIODOGRAM_DATASETS,
    RAW_PAIR_CONTRACT_VERSION,
)
from twirl.vetting.harmonic_training import run_harmonic_teacher_training  # noqa: E402


def _light_curve(kind: str, *, seed: int, n_cadences: int = 256) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(seed)
    time = 2459825.0 + np.arange(n_cadences) * 200.0 / 86400.0
    period = 0.25
    t0 = 2459825.08
    phase = ((time - t0 + 0.5 * period) % period) / period - 0.5
    small = np.ones(n_cadences)
    primary = np.ones(n_cadences)
    if kind == "planet_like":
        event = np.abs(phase) < 0.025
        small[event] -= 0.12
        primary[event] -= 0.07
    elif kind == "eclipse_contact":
        primary_event = np.abs(phase) < 0.045
        secondary_event = np.abs(np.abs(phase) - 0.5) < 0.035
        small[primary_event] -= 0.20
        primary[primary_event] -= 0.15
        small[secondary_event] -= 0.08
        primary[secondary_event] -= 0.06
    elif kind == "smooth_variable":
        small += 0.08 * np.sin(2.0 * np.pi * phase)
        primary += 0.06 * np.sin(2.0 * np.pi * phase + 0.1)
    else:
        small += 0.02 * np.sin(2.0 * np.pi * np.arange(n_cadences) / 41.0)
        primary -= 0.015 * np.sin(2.0 * np.pi * np.arange(n_cadences) / 41.0)
    small += rng.normal(0.0, 0.008, n_cadences)
    primary += rng.normal(0.0, 0.010, n_cadences)
    raw_small = 100.0 * small + rng.normal(0.0, 1.0, n_cadences)
    raw_primary = 170.0 * primary + rng.normal(0.0, 1.4, n_cadences)
    periods = np.geomspace(0.12, 1.5, 64)
    peak = np.exp(-0.5 * np.square((np.log(periods) - np.log(period)) / 0.08))
    return {
        "time": time.astype(np.float64),
        "cadenceno": np.arange(1000, 1000 + n_cadences, dtype=np.int64),
        "orbitid": np.where(np.arange(n_cadences) < n_cadences // 2, 119, 120).astype(np.int32),
        "quality": np.zeros(n_cadences, dtype=np.int32),
        "raw_flux_small": raw_small.astype(np.float64),
        "raw_flux_err_small": np.full(n_cadences, 1.2, dtype=np.float64),
        "raw_flux_primary": raw_primary.astype(np.float64),
        "raw_flux_err_primary": np.full(n_cadences, 1.8, dtype=np.float64),
        "det_flux_adp_sml": small.astype(np.float64),
        "det_flux_adp": primary.astype(np.float64),
        "bls_log_period_grid": np.log10(periods).astype(np.float32),
        "bls_power_small": (peak + 0.02).astype(np.float32),
        "bls_sde_small": (8.0 * peak).astype(np.float32),
        "bls_power_primary": (0.8 * peak + 0.02).astype(np.float32),
        "bls_sde_primary": (6.0 * peak).astype(np.float32),
    }


def _write_group(group: h5py.Group, payload: dict[str, np.ndarray]) -> None:
    for name in NATIVE_DATASETS + PERIODOGRAM_DATASETS:
        group.create_dataset(name, data=payload[name], compression="lzf", shuffle=True)


def _write_paired_original(group: h5py.Group, payload: dict[str, np.ndarray]) -> None:
    original = _light_curve("other", seed=10_000 + int(group.attrs["tic"]))
    for name in NATIVE_DATASETS[4:]:
        group.create_dataset(
            f"paired_original_{name}", data=original[name], compression="lzf", shuffle=True
        )


def build_smoke_fixture(out_dir: Path, *, seed: int = 56) -> tuple[Path, Path]:
    fixture = out_dir / "fixture"
    fixture.mkdir(parents=True, exist_ok=False)
    table_path = fixture / "training.csv"
    native_path = fixture / "native.h5"
    rows: list[dict[str, object]] = []
    specifications = (
        ("planet_like", "planet_like", "preserve", "p", 12),
        ("eclipsing_binary_or_pceb", "eclipse_contact", "preserve", "2p", 12),
        ("stellar_variability", "smooth_variable", "preserve", "", 12),
        ("instrumental_or_systematic", "other", "reject", "", 8),
        ("uncertain", "other", "reject", "", 8),
    )
    index = 0
    for human_label, morphology, preserve, harmonic, count in specifications:
        for _ in range(count):
            tic = 100_000 + index
            rows.append(
                {
                    "review_id": f"real:{index:04d}",
                    "tic": tic,
                    "source_kind": "real_candidate",
                    "is_injected_row": False,
                    "human_label": human_label,
                    "period_d": 0.25,
                    "t0_bjd": 2459825.08,
                    "duration_min": 12.0,
                    "depth": 0.1,
                    "depth_snr": 12.0,
                    "sde_max": 9.0,
                    "tmag": 18.0 + 0.02 * index,
                    "morphology_target_v1": morphology,
                    "morphology_include_v1": True,
                    "preserve_target_v1": preserve,
                    "preserve_include_v1": True,
                    "harmonic_target_v1": harmonic,
                    "harmonic_include_v1": bool(harmonic),
                    "broad_preserve_only": False,
                    "model_target_policy_version": HARMONIC_CNN_TARGET_POLICY,
                }
            )
            index += 1
    for injection_index in range(10):
        tic = 200_000 + injection_index
        rows.append(
            {
                "review_id": f"injection:{injection_index:04d}",
                "tic": tic,
                "source_kind": "injection_recovery",
                "is_injected_row": True,
                "injection_id": f"smoke_{injection_index:04d}",
                "human_label": "planet_like",
                "period_d": 0.25,
                "t0_bjd": 2459825.08,
                "duration_min": 12.0,
                "depth": 0.1,
                "depth_snr": 12.0,
                "sde_max": 10.0,
                "adp_sml_sde": 10.0,
                "tmag": 18.5,
                "truth_period_d": 0.25,
                "truth_radius_rearth": 4.0,
                "truth_model_depth": 0.12,
                "boundary_period_bin_key": "smoke_p",
                "boundary_radius_bin_key": "smoke_r",
                "boundary_tmag_bin_key": "smoke_m",
                "morphology_target_v1": "planet_like",
                "morphology_include_v1": True,
                "preserve_target_v1": "preserve",
                "preserve_include_v1": True,
                "harmonic_target_v1": "p",
                "harmonic_include_v1": True,
                "broad_preserve_only": False,
                "model_target_policy_version": HARMONIC_CNN_TARGET_POLICY,
            }
        )
    table = pd.DataFrame(rows).sample(frac=1.0, random_state=seed).reset_index(drop=True)
    table.to_csv(table_path, index=False)

    with h5py.File(native_path, "w") as h5:
        h5.attrs["contract_version"] = RAW_PAIR_CONTRACT_VERSION
        h5.attrs["training_table_sha256"] = hashlib.sha256(
            table_path.read_bytes()
        ).hexdigest()
        h5.attrs["time_system"] = "BJD"
        h5.attrs["external_quality_policy_contract"] = (
            EXTERNAL_QUALITY_POLICY_CONTRACT
        )
        h5.attrs["effective_quality_policy"] = EFFECTIVE_QUALITY_POLICY
        h5.attrs["cadence_reference_contract_version"] = (
            "synthetic_smoke_cadence_reference_v1"
        )
        h5.attrs["cadence_reference_cadence_authority"] = "qlp_cam_quat"
        h5.attrs["cadence_reference_quality_authority"] = (
            "spoc_and_qlp_quality_flags"
        )
        h5.attrs["cadence_reference_table"] = "synthetic-smoke://reference"
        h5.attrs["cadence_reference_manifest"] = "synthetic-smoke://manifest"
        h5.attrs["cadence_reference_table_sha256"] = "1" * 64
        h5.attrs["cadence_reference_manifest_sha256"] = "2" * 64
        h5.attrs["cadence_reference_source_declaration_sha256"] = "3" * 64
        h5.attrs["authority_exclusion_policy_contract"] = (
            AUTHORITY_EXCLUSION_POLICY_CONTRACT
        )
        h5.attrs["authority_exclusion_external_bit"] = (
            AUTHORITY_EXCLUSION_EXTERNAL_BIT
        )
        h5.attrs["authority_exclusions_sha256"] = "4" * 64
        h5.attrs["n_authority_exclusions"] = 0
        for name, channels in CHANNEL_CONTRACT.items():
            h5.attrs[name] = json.dumps(channels)
        h5.create_group("targets")
        h5.create_group("injections")
        for row_index, row in table.iterrows():
            morphology = str(row["morphology_target_v1"])
            payload = _light_curve(morphology, seed=seed + row_index)
            if bool(row["is_injected_row"]):
                group = h5.create_group(f"injections/{row['injection_id']}")
            else:
                group = h5.create_group(f"targets/{int(row['tic']):016d}")
            group.attrs["tic"] = int(row["tic"])
            _write_group(group, payload)
            group.attrs["quality_policy_contract"] = (
                EXTERNAL_QUALITY_POLICY_CONTRACT
            )
            group.attrs["n_cad_total"] = len(payload["quality"])
            group.attrs["n_cad_internal_bad"] = 0
            group.attrs["n_cad_external_bad"] = 0
            group.attrs["n_cad_external_only_bad"] = 0
            group.attrs["n_cad_authority_excluded"] = 0
            group.attrs["n_cad_effective_bad"] = 0
            if bool(row["is_injected_row"]):
                _write_paired_original(group, payload)
        h5.attrs["quality_overlay_n_cad_total"] = int(256 * len(table))
        for name in (
            "n_cad_internal_bad",
            "n_cad_external_bad",
            "n_cad_external_only_bad",
            "n_cad_authority_excluded",
            "n_cad_effective_bad",
        ):
            h5.attrs[f"quality_overlay_{name}"] = 0
    return table_path, native_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--profiles", default="metadata_only,full_combined")
    parser.add_argument("--epochs", type=int, default=1)
    parser.add_argument("--pretrain-epochs", type=int, default=1)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=False)
    table_path, native_path = build_smoke_fixture(args.out_dir)
    profiles = tuple(value for value in args.profiles.split(",") if value)
    summary = run_harmonic_teacher_training(
        training_table=table_path,
        native_h5=native_path,
        out_dir=args.out_dir / "training",
        profiles=profiles,
        train_config=HarmonicTrainConfig(
            epochs=args.epochs,
            batch_size=args.batch_size,
            patience=max(1, args.epochs),
        ),
        workers=args.workers,
        pretrain_epochs=args.pretrain_epochs,
        require_cuda=True,
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=True))


if __name__ == "__main__":
    main()
