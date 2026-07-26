"""Seven-sector Teacher v3 training on observation-addressed native inputs.

Teacher v3 is a new data release, not a new neural-network architecture.  It
retains :data:`s56_harmonic_cnn_v1`, fixes the primary and metadata-baseline
profiles before opening the test partition, and fits one temperature per
profile from the concatenated five-fold development OOF logits.
"""
from __future__ import annotations

from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from twirl.vetting.harmonic_cnn import (
    HARMONIC_CLASSES,
    MODEL_VERSION,
    MORPHOLOGY_CLASSES,
    PRESERVE_CLASSES,
    HarmonicModelConfig,
    HarmonicTrainConfig,
    build_harmonic_cnn,
)
from twirl.vetting.harmonic_dataset import (
    HarmonicNativeDataset,
    MetadataNormalization,
    build_metadata_matrix,
    prepare_harmonic_training_rows,
)
from twirl.vetting.harmonic_inputs import (
    CHANNEL_CONTRACT,
    RAW_PAIR_CONTRACT_VERSION,
    verify_raw_pair_contract,
)
from twirl.vetting.harmonic_training import (
    _calibration_by_source,
    _evaluate_loader,
    _file_sha256,
    _loader,
    _softmax,
    _subset_metrics,
    _train_one_fold,
    classification_metrics,
    expected_calibration_error,
    fit_temperature,
    injection_truth_human_audit,
    validate_teacher_input_provenance,
)
from twirl.vetting.teacher_native_registry import (
    file_sha256,
    read_table,
    validate_native_input_registry,
    validate_native_input_registry_path,
)
from twirl.vetting.teacher_v3 import (
    TEACHER_V3_CORPUS_VERSION,
    TEACHER_V3_RUN_NAME,
)


TEACHER_V3_RUN_ID = "teacher_v3_s56_s62_a2v1_current_adp"
TEACHER_V3_CHECKPOINT_NAMESPACE = (
    "twirl_teacher_v3_s56_s62_native_v1"
)
TEACHER_V3_METADATA_CHECKPOINT_NAMESPACE = (
    "twirl_teacher_v3_s56_s62_metadata_baseline_v1"
)
TEACHER_V3_NATIVE_MANIFEST_SCHEMA = (
    "twirl_teacher_v3_native_input_manifest_v1"
)
TEACHER_V3_CALIBRATION_SCHEMA = (
    "twirl_teacher_v3_pooled_oof_calibration_v1"
)
TEACHER_V3_CHECKPOINT_MANIFEST_SCHEMA = (
    "twirl_teacher_v3_selected_checkpoint_manifest_v1"
)
TEACHER_V3_PRIMARY_PROFILE = "shape_plus_periodogram_bls"
TEACHER_V3_BASELINE_PROFILE = "metadata_only"
TEACHER_V3_PROFILES: tuple[str, str] = (
    TEACHER_V3_BASELINE_PROFILE,
    TEACHER_V3_PRIMARY_PROFILE,
)
TEACHER_V3_SECTORS: tuple[int, ...] = tuple(range(56, 63))


def _json_payload(payload: Mapping[str, Any]) -> bytes:
    return (
        json.dumps(
            payload,
            indent=2,
            sort_keys=True,
            allow_nan=True,
        )
        + "\n"
    ).encode("utf-8")


def _write_json(path: Path, payload: Mapping[str, Any]) -> str:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    serialized = _json_payload(payload)
    path.write_bytes(serialized)
    return hashlib.sha256(serialized).hexdigest()


def _truth(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values.fillna(False)
    return (
        values.fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
        .isin({"1", "1.0", "true", "t", "yes", "y"})
    )


def validate_teacher_v3_training_table(source: pd.DataFrame) -> None:
    """Validate the frozen release identity before touching native inputs."""

    required = {
        "review_id",
        "teacher_v3_observation_id",
        "teacher_run_name",
        "teacher_v3_corpus_version",
        "teacher_architecture_version",
        "teacher_v3_training_include",
        "sector",
        "tic",
        "fixed_split",
        "cv_fold",
        "native_h5_path",
        "native_group_path",
    }
    missing = sorted(required - set(source.columns))
    if missing:
        raise KeyError(f"Teacher v3 training table is missing columns: {missing}")
    if source.empty:
        raise ValueError("Teacher v3 training table is empty")
    for name in ("review_id", "teacher_v3_observation_id"):
        values = source[name].fillna("").astype(str).str.strip()
        if values.eq("").any() or values.duplicated().any():
            raise ValueError(f"{name} must be nonblank and globally unique")
    identities = {
        "teacher_run_name": TEACHER_V3_RUN_NAME,
        "teacher_v3_corpus_version": TEACHER_V3_CORPUS_VERSION,
        "teacher_architecture_version": MODEL_VERSION,
    }
    for name, expected in identities.items():
        observed = set(source[name].fillna("").astype(str).str.strip())
        if observed != {expected}:
            raise ValueError(
                f"{name} must be exactly {expected!r}; observed={sorted(observed)}"
            )
    sectors = {
        int(value)
        for value in pd.to_numeric(source["sector"], errors="raise")
    }
    if sectors != set(TEACHER_V3_SECTORS):
        raise ValueError(
            f"Teacher v3 sectors must be {list(TEACHER_V3_SECTORS)}; "
            f"observed={sorted(sectors)}"
        )
    active = _truth(source["teacher_v3_training_include"])
    if not active.any():
        raise ValueError("Teacher v3 has no active training rows")
    for name in ("native_h5_path", "native_group_path"):
        blank = source.loc[active, name].fillna("").astype(str).str.strip().eq("")
        if blank.any():
            raise ValueError(
                f"{name} is blank for {int(blank.sum())} active Teacher v3 rows"
            )
    split = source["fixed_split"].fillna("").astype(str).str.strip().str.lower()
    if set(split) != {"development", "test"}:
        raise ValueError(
            "fixed_split must contain exactly development and test"
        )
    folds = pd.to_numeric(source["cv_fold"], errors="raise").astype(int)
    if set(folds[split.eq("development")]) != set(range(5)):
        raise ValueError("development rows must cover all five folds 0..4")
    if not folds[split.eq("test")].eq(-1).all():
        raise ValueError("fixed test rows must have cv_fold=-1")


def _native_mapping_rows(rows: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "sector",
        "tic",
        "native_h5_path",
        "native_group_path",
        "is_injected_row",
    ]
    missing = sorted(set(columns) - set(rows.columns))
    if missing:
        raise KeyError(f"prepared Teacher v3 rows are missing columns: {missing}")
    mappings = rows.loc[:, columns].copy()
    mappings["sector"] = pd.to_numeric(
        mappings["sector"], errors="raise"
    ).astype(np.int64)
    mappings["tic"] = pd.to_numeric(
        mappings["tic"], errors="raise"
    ).astype(np.int64)
    mappings["native_h5_path"] = (
        mappings["native_h5_path"]
        .fillna("")
        .astype(str)
        .str.strip()
        .map(lambda value: str(Path(value).expanduser().resolve()))
    )
    mappings["native_group_path"] = (
        mappings["native_group_path"].fillna("").astype(str).str.strip()
    )
    mappings["is_injected_row"] = _truth(mappings["is_injected_row"])
    if mappings["native_h5_path"].eq("").any():
        raise ValueError("prepared Teacher v3 rows contain blank native paths")
    if mappings["native_group_path"].eq("").any():
        raise ValueError("prepared Teacher v3 rows contain blank native groups")
    identity = mappings.loc[
        :,
        ["native_h5_path", "native_group_path", "sector", "tic"],
    ].drop_duplicates()
    collisions = identity.duplicated(
        ["native_h5_path", "native_group_path"],
        keep=False,
    )
    if collisions.any():
        examples = identity.loc[collisions].head(5).to_dict("records")
        raise ValueError(
            "one native HDF5 group maps to more than one (sector,TIC); "
            f"first={examples}"
        )
    return mappings.drop_duplicates().reset_index(drop=True)


def _validate_real_registry(
    *,
    mappings: pd.DataFrame,
    registry_path: Path,
    registry_summary_path: Path,
) -> dict[str, Any]:
    audit = validate_native_input_registry_path(
        registry_path=registry_path,
        summary_path=registry_summary_path,
        expected_contract_version=RAW_PAIR_CONTRACT_VERSION,
    )
    registry = validate_native_input_registry(
        read_table(registry_path),
        path_base=Path(registry_path).parent,
        expected_contract_version=RAW_PAIR_CONTRACT_VERSION,
    )
    expected = (
        mappings.loc[
            ~mappings["is_injected_row"],
            ["sector", "tic", "native_h5_path", "native_group_path"],
        ]
        .drop_duplicates()
        .sort_values(["sector", "tic"], kind="stable")
        .reset_index(drop=True)
    )
    observed = (
        registry.loc[
            :,
            ["sector", "tic", "native_h5_path", "native_group_path"],
        ]
        .sort_values(["sector", "tic"], kind="stable")
        .reset_index(drop=True)
    )
    if not expected.equals(observed):
        comparison = expected.merge(
            observed,
            on=["sector", "tic", "native_h5_path", "native_group_path"],
            how="outer",
            indicator=True,
        )
        examples = comparison.loc[comparison["_merge"].ne("both")].head(5)
        raise ValueError(
            "real Teacher v3 rows do not exactly match the frozen native "
            f"registry; first={examples.to_dict('records')}"
        )
    return audit


def build_teacher_v3_native_manifest(
    *,
    rows: pd.DataFrame,
    registry_path: Path,
    registry_summary_path: Path,
) -> dict[str, Any]:
    """Fully validate and hash all row-addressed Teacher v3 HDF5 inputs."""

    import h5py

    mappings = _native_mapping_rows(rows)
    registry_path = Path(registry_path).resolve()
    registry_summary_path = Path(registry_summary_path).resolve()
    registry_audit = _validate_real_registry(
        mappings=mappings,
        registry_path=registry_path,
        registry_summary_path=registry_summary_path,
    )
    records: list[dict[str, Any]] = []
    for path_text, expected_rows in mappings.groupby(
        "native_h5_path",
        sort=True,
    ):
        path = Path(path_text)
        if not path.is_file():
            raise FileNotFoundError(f"missing Teacher v3 native HDF5: {path}")
        stat_before = path.stat()
        verification = verify_raw_pair_contract(
            path,
            require_errors=True,
            require_periodograms=True,
        )
        if not verification["passed"]:
            raise RuntimeError(
                f"Teacher v3 native contract failed for {path}: "
                f"{verification['failures'][:10]}"
            )
        expected_groups = set(expected_rows["native_group_path"])
        with h5py.File(path, "r") as h5:
            contract = str(h5.attrs.get("contract_version", ""))
            observed_groups = {
                f"{root}/{key}"
                for root in ("targets", "injections")
                if root in h5
                for key in h5[root]
            }
            if observed_groups != expected_groups:
                missing = sorted(expected_groups - observed_groups)[:10]
                extra = sorted(observed_groups - expected_groups)[:10]
                raise ValueError(
                    f"native group set mismatch in {path}: "
                    f"missing={missing}, extra={extra}"
                )
            for row in (
                expected_rows.loc[
                    :,
                    ["native_group_path", "sector", "tic"],
                ]
                .drop_duplicates()
                .to_dict("records")
            ):
                group = h5[str(row["native_group_path"])]
                try:
                    observed_sector = int(group.attrs["sector"])
                    observed_tic = int(group.attrs["tic"])
                except (KeyError, TypeError, ValueError) as exc:
                    raise ValueError(
                        f"{path}:{row['native_group_path']} lacks scalar "
                        "sector/TIC identity attributes"
                    ) from exc
                if (
                    observed_sector,
                    observed_tic,
                ) != (int(row["sector"]), int(row["tic"])):
                    raise ValueError(
                        f"{path}:{row['native_group_path']} declares "
                        f"(sector,TIC)=({observed_sector},{observed_tic}), "
                        f"expected=({row['sector']},{row['tic']})"
                    )
            source_table_sha256 = str(
                h5.attrs.get("training_table_sha256", "")
            )
        digest = file_sha256(path)
        stat_after = path.stat()
        if (
            stat_before.st_size != stat_after.st_size
            or stat_before.st_mtime_ns != stat_after.st_mtime_ns
        ):
            raise RuntimeError(
                f"native HDF5 changed while it was validated: {path}"
            )
        records.append(
            {
                "path": str(path),
                "sha256": digest,
                "size_bytes": int(stat_after.st_size),
                "mtime_ns": int(stat_after.st_mtime_ns),
                "contract_version": contract,
                "source_training_table_sha256": source_table_sha256,
                "n_training_rows": int(len(expected_rows)),
                "n_unique_groups": int(len(expected_groups)),
                "group_counts": verification["counts"],
                "external_quality_counts": verification[
                    "external_quality_counts"
                ],
            }
        )
    observed_sectors = set(
        pd.to_numeric(mappings["sector"], errors="raise").astype(int)
    )
    if observed_sectors != set(TEACHER_V3_SECTORS):
        raise ValueError(
            "native manifest does not cover exactly sectors 56--62"
        )
    if len(records) != len(TEACHER_V3_SECTORS):
        raise ValueError(
            "Teacher v3 requires exactly one native HDF5 per sector; "
            f"found {len(records)}"
        )
    manifest: dict[str, Any] = {
        "schema_version": TEACHER_V3_NATIVE_MANIFEST_SCHEMA,
        "run_id": TEACHER_V3_RUN_ID,
        "model_version": MODEL_VERSION,
        "native_contract_version": RAW_PAIR_CONTRACT_VERSION,
        "sectors": list(TEACHER_V3_SECTORS),
        "n_training_rows": int(len(rows)),
        "n_unique_native_groups": int(
            mappings[
                ["native_h5_path", "native_group_path"]
            ].drop_duplicates().shape[0]
        ),
        "registry": {
            "path": str(registry_path),
            "sha256": _file_sha256(registry_path),
            "summary_path": str(registry_summary_path),
            "summary_sha256": _file_sha256(registry_summary_path),
            "audit": registry_audit,
        },
        "native_files": records,
    }
    return manifest


def _manifest_digest(manifest: Mapping[str, Any]) -> str:
    return hashlib.sha256(_json_payload(manifest)).hexdigest()


def build_teacher_v3_input_provenance(
    *,
    training_table: Path,
    native_manifest: Mapping[str, Any],
) -> dict[str, str]:
    """Bind checkpoints to the table and canonical seven-file manifest."""

    manifest_sha256 = _manifest_digest(native_manifest)
    return {
        "checkpoint_namespace": TEACHER_V3_CHECKPOINT_NAMESPACE,
        "input_contract_version": RAW_PAIR_CONTRACT_VERSION,
        # Kept for compatibility with the v1 checkpoint validator.  For v3
        # this is the digest of the canonical multi-file manifest, not one H5.
        "native_h5_sha256": manifest_sha256,
        "native_manifest_sha256": manifest_sha256,
        "training_table_sha256": _file_sha256(training_table),
        "native_registry_sha256": str(
            native_manifest["registry"]["sha256"]
        ),
        "native_registry_summary_sha256": str(
            native_manifest["registry"]["summary_sha256"]
        ),
    }


def _oof_frame_hash(frame: pd.DataFrame) -> str:
    columns = [
        "review_id",
        "morphology_target",
        *[f"logit_{label}" for label in MORPHOLOGY_CLASSES],
    ]
    ordered = frame.loc[:, columns].sort_values(
        "review_id",
        kind="stable",
    )
    payload = ordered.to_csv(
        index=False,
        lineterminator="\n",
        float_format="%.17g",
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def calibrate_teacher_v3_profile_oof(
    *,
    rows: pd.DataFrame,
    out_dir: Path,
    profile: str,
    input_provenance: Mapping[str, str],
) -> dict[str, Any]:
    """Fit and apply one temperature to concatenated five-fold OOF logits."""

    import torch

    out_dir = Path(out_dir)
    parts: list[pd.DataFrame] = []
    for fold in range(5):
        path = out_dir / profile / f"fold_{fold}" / "validation_predictions.csv"
        part = pd.read_csv(path)
        part["cv_fold"] = fold
        parts.append(part)
    predictions = pd.concat(parts, ignore_index=True)
    if predictions["review_id"].duplicated().any():
        raise RuntimeError(
            f"duplicate development OOF predictions for profile {profile}"
        )
    expected_ids = set(
        rows.loc[rows["fixed_split"].eq("development"), "review_id"].astype(str)
    )
    observed_ids = set(predictions["review_id"].astype(str))
    if observed_ids != expected_ids:
        raise RuntimeError(
            f"{profile} OOF predictions do not exactly cover development rows"
        )
    logit_columns = [f"logit_{label}" for label in MORPHOLOGY_CLASSES]
    missing = sorted(set(logit_columns) - set(predictions.columns))
    if missing:
        raise RuntimeError(
            f"{profile} OOF predictions lack raw logits: {missing}"
        )
    logits = predictions.loc[:, logit_columns].to_numpy(dtype=float)
    truth = predictions["morphology_target"].to_numpy(dtype=int)
    temperature = fit_temperature(logits, truth)
    probability = _softmax(logits, temperature=temperature)
    predictions["morphology_prediction"] = probability.argmax(axis=1)
    for index, label in enumerate(MORPHOLOGY_CLASSES):
        predictions[f"p_{label}"] = probability[:, index]
    oof_hash = _oof_frame_hash(predictions)
    active = truth >= 0
    if not np.any(active):
        raise RuntimeError(f"{profile} OOF predictions have no morphology targets")
    uncalibrated = _softmax(logits)
    calibration = {
        "schema_version": TEACHER_V3_CALIBRATION_SCHEMA,
        "run_id": TEACHER_V3_RUN_ID,
        "profile": profile,
        "scope": "concatenated_five_fold_development_oof_logits",
        "n_oof_rows": int(len(predictions)),
        "oof_logits_sha256": oof_hash,
        "temperature": float(temperature),
        "uncalibrated_cross_entropy": float(
            -np.mean(
                np.log(
                    np.clip(
                        uncalibrated[
                            np.flatnonzero(active),
                            truth[active],
                        ],
                        1.0e-12,
                        1.0,
                    )
                )
            )
        ),
        "calibrated_cross_entropy": float(
            -np.mean(
                np.log(
                    np.clip(
                        probability[
                            np.flatnonzero(active),
                            truth[active],
                        ],
                        1.0e-12,
                        1.0,
                    )
                )
            )
        ),
        **dict(input_provenance),
    }
    calibration_path = out_dir / profile / "pooled_oof_calibration.json"
    calibration_sha256 = _write_json(calibration_path, calibration)

    row_lookup = rows.set_index("review_id", drop=False)
    for fold in range(5):
        fold_mask = predictions["cv_fold"].eq(fold)
        fold_predictions = predictions.loc[fold_mask].drop(
            columns=["cv_fold"]
        )
        fold_path = (
            out_dir
            / profile
            / f"fold_{fold}"
            / "validation_predictions.csv"
        )
        fold_predictions.to_csv(fold_path, index=False)
        checkpoint_path = (
            out_dir / profile / f"fold_{fold}" / "teacher.pt"
        )
        checkpoint = torch.load(
            checkpoint_path,
            map_location="cpu",
            weights_only=False,
        )
        validate_teacher_input_provenance(
            checkpoint,
            expected=input_provenance,
            artifact=f"Teacher v3 {profile} fold-{fold} checkpoint",
        )
        if checkpoint.get("temperature_calibration_scope") != (
            "pending_pooled_oof_development"
        ):
            raise RuntimeError(
                f"{checkpoint_path} was not held for pooled OOF calibration"
            )
        checkpoint["temperature"] = float(temperature)
        checkpoint["temperature_calibration_scope"] = (
            "pooled_oof_development"
        )
        checkpoint["pooled_oof_calibration_sha256"] = calibration_sha256
        checkpoint["run_id"] = TEACHER_V3_RUN_ID
        checkpoint["release_name"] = TEACHER_V3_RUN_NAME
        temporary = checkpoint_path.with_suffix(".pt.tmp")
        temporary.unlink(missing_ok=True)
        torch.save(checkpoint, temporary)
        temporary.replace(checkpoint_path)

        fold_truth = fold_predictions[
            "morphology_target"
        ].to_numpy(dtype=int)
        fold_probability = fold_predictions.loc[
            :,
            [f"p_{label}" for label in MORPHOLOGY_CLASSES],
        ].to_numpy(dtype=float)
        fold_rows = row_lookup.loc[
            fold_predictions["review_id"].astype(str)
        ].reset_index(drop=True)
        metrics_path = (
            out_dir / profile / f"fold_{fold}" / "metrics.json"
        )
        metrics = json.loads(metrics_path.read_text())
        metrics["morphology"] = classification_metrics(
            fold_truth,
            fold_probability,
            classes=MORPHOLOGY_CLASSES,
        )
        metrics["morphology_by_source"] = _subset_metrics(
            fold_rows,
            fold_truth,
            fold_probability,
        )
        metrics["calibration"] = expected_calibration_error(
            fold_truth,
            fold_probability,
        )
        metrics["temperature"] = float(temperature)
        metrics["temperature_calibration_scope"] = (
            "pooled_oof_development"
        )
        metrics["pooled_oof_calibration_sha256"] = calibration_sha256
        _write_json(metrics_path, metrics)
    development_rows = row_lookup.loc[
        predictions["review_id"].astype(str)
    ].reset_index(drop=True)
    source_metrics = _subset_metrics(
        development_rows,
        truth,
        probability,
    )
    return {
        "profile": profile,
        "temperature": float(temperature),
        "calibration_path": str(calibration_path),
        "calibration_sha256": calibration_sha256,
        "oof_logits_sha256": oof_hash,
        "metrics": {
            "morphology_by_source": source_metrics,
            "calibration": _calibration_by_source(
                development_rows,
                truth,
                probability,
            ),
        },
    }


def _evaluate_fixed_profile(
    *,
    rows: pd.DataFrame,
    out_dir: Path,
    profile: str,
    train_config: HarmonicTrainConfig,
    workers: int,
    input_provenance: Mapping[str, str],
    require_cuda: bool,
) -> dict[str, Any]:
    import torch

    test_mask = rows["fixed_split"].eq("test").to_numpy()
    test_rows = rows.loc[test_mask].reset_index(drop=True)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if require_cuda and device.type != "cuda":
        raise RuntimeError("CUDA was required but is unavailable")
    morphology_members: list[np.ndarray] = []
    preserve_members: list[np.ndarray] = []
    harmonic_members: list[np.ndarray] = []
    truth: np.ndarray | None = None
    preserve_truth: np.ndarray | None = None
    harmonic_truth: np.ndarray | None = None
    review_ids: list[str] | None = None
    temperatures: set[float] = set()
    for fold in range(5):
        checkpoint_path = (
            Path(out_dir) / profile / f"fold_{fold}" / "teacher.pt"
        )
        checkpoint = torch.load(
            checkpoint_path,
            map_location="cpu",
            weights_only=False,
        )
        validate_teacher_input_provenance(
            checkpoint,
            expected=input_provenance,
            artifact=f"Teacher v3 {profile} fold-{fold} checkpoint",
        )
        if checkpoint.get("temperature_calibration_scope") != (
            "pooled_oof_development"
        ):
            raise RuntimeError(
                f"{checkpoint_path} lacks pooled OOF calibration"
            )
        temperatures.add(float(checkpoint["temperature"]))
        norm = checkpoint["metadata_normalization"]
        normalization = MetadataNormalization(
            columns=tuple(norm["columns"]),
            center=tuple(norm["center"]),
            scale=tuple(norm["scale"]),
        )
        metadata, _ = build_metadata_matrix(
            rows,
            fit_mask=np.zeros(len(rows), dtype=bool),
            normalization=normalization,
        )
        dataset = HarmonicNativeDataset(
            rows,
            native_h5=None,
            metadata=metadata,
            profile=profile,
        )
        loader = _loader(
            dataset,
            np.flatnonzero(test_mask),
            batch_size=train_config.batch_size,
            shuffle=False,
            workers=workers,
            seed=train_config.seed,
        )
        model = build_harmonic_cnn(
            HarmonicModelConfig(**checkpoint["model_config"]),
            profile=profile,
        ).to(device)
        model.load_state_dict(checkpoint["model_state_dict"])
        scored = _evaluate_loader(model, loader, device=device)
        morphology_members.append(
            _softmax(
                scored["morphology_logits"],
                temperature=float(checkpoint["temperature"]),
            )
        )
        preserve_members.append(_softmax(scored["preserve_logits"]))
        harmonic_members.append(_softmax(scored["harmonic_logits"]))
        if truth is None:
            truth = scored["morphology_target"]
            preserve_truth = scored["preserve_target"]
            harmonic_truth = scored["harmonic_target"]
            review_ids = scored["review_id"]
        elif review_ids != scored["review_id"]:
            raise RuntimeError(
                f"{profile} test ordering changed between ensemble members"
            )
    if len(temperatures) != 1:
        raise RuntimeError(
            f"{profile} fold checkpoints disagree on pooled temperature"
        )
    assert truth is not None
    assert preserve_truth is not None
    assert harmonic_truth is not None
    assert review_ids is not None
    morphology_probability = np.mean(morphology_members, axis=0)
    preserve_probability = np.mean(preserve_members, axis=0)
    harmonic_probability = np.mean(harmonic_members, axis=0)
    metrics: dict[str, Any] = {
        "morphology_by_source": _subset_metrics(
            test_rows,
            truth,
            morphology_probability,
        ),
        "calibration": _calibration_by_source(
            test_rows,
            truth,
            morphology_probability,
        ),
        "preserve": classification_metrics(
            preserve_truth,
            preserve_probability,
            classes=PRESERVE_CLASSES,
        ),
        "harmonic": classification_metrics(
            harmonic_truth,
            harmonic_probability,
            classes=HARMONIC_CLASSES,
        ),
        "temperature": float(next(iter(temperatures))),
        "temperature_calibration_scope": "pooled_oof_development",
    }
    uncertain = (
        test_rows["human_label"].fillna("").astype(str).eq("uncertain")
    ).to_numpy()
    metrics["uncertain_masked_sensitivity"] = _subset_metrics(
        test_rows.loc[~uncertain].reset_index(drop=True),
        truth[~uncertain],
        morphology_probability[~uncertain],
    )
    predictions = test_rows.loc[
        :,
        [
            "review_id",
            "sector",
            "tic",
            "human_label",
            "is_injected_row",
        ],
    ].copy()
    predictions["morphology_target_index"] = truth
    predictions["morphology_prediction_index"] = (
        morphology_probability.argmax(axis=1)
    )
    for index, label in enumerate(MORPHOLOGY_CLASSES):
        predictions[f"p_{label}"] = morphology_probability[:, index]
    predictions["p_preserve"] = preserve_probability[:, 1]
    predictions["predicted_period_factor"] = np.asarray(
        [0.25, 1.0 / 3.0, 0.5, 1.0, 2.0, 3.0, 4.0]
    )[harmonic_probability.argmax(axis=1)]
    predictions.to_csv(
        Path(out_dir) / profile / "fixed_test_predictions.csv",
        index=False,
    )
    _write_json(
        Path(out_dir) / profile / "fixed_test_metrics.json",
        metrics,
    )
    return metrics


def _profile_development_record(
    result: Mapping[str, Any],
) -> dict[str, Any]:
    metrics = result["metrics"]["morphology_by_source"]
    all_metrics = metrics["all"]
    real_metrics = metrics["real"]
    per_class = real_metrics["per_class"]
    return {
        "profile": str(result["profile"]),
        "fixed_role": (
            "primary"
            if result["profile"] == TEACHER_V3_PRIMARY_PROFILE
            else "metadata_baseline"
        ),
        "temperature": float(result["temperature"]),
        "validation_macro_f1": float(all_metrics["macro_f1"]),
        "validation_balanced_accuracy": float(
            all_metrics["balanced_accuracy"]
        ),
        "validation_real_planet_recall": float(
            per_class.get("planet_like", {}).get("recall", float("nan"))
        ),
        "validation_real_eb_recall": float(
            per_class.get("eclipse_contact", {}).get(
                "recall",
                float("nan"),
            )
        ),
        "validation_real_variable_recall": float(
            per_class.get("smooth_variable", {}).get(
                "recall",
                float("nan"),
            )
        ),
        "validation_real_other_recall": float(
            per_class.get("other", {}).get("recall", float("nan"))
        ),
        "validation_ece": float(
            result["metrics"]["calibration"]["all"]["ece"]
        ),
    }


def _write_selected_checkpoint_manifest(
    *,
    out_dir: Path,
    input_provenance: Mapping[str, str],
    selected_profile: str = TEACHER_V3_PRIMARY_PROFILE,
    filename: str = "selected_checkpoint_manifest.json",
    selection_policy: str = "fixed_before_test",
) -> Path:
    import torch

    records: list[dict[str, Any]] = []
    for fold in range(5):
        path = (
            Path(out_dir)
            / selected_profile
            / f"fold_{fold}"
            / "teacher.pt"
        )
        checkpoint = torch.load(
            path,
            map_location="cpu",
            weights_only=False,
        )
        validate_teacher_input_provenance(
            checkpoint,
            expected=input_provenance,
            artifact=f"Teacher v3 selected fold-{fold} checkpoint",
        )
        if checkpoint.get("temperature_calibration_scope") != (
            "pooled_oof_development"
        ):
            raise RuntimeError(
                f"selected Teacher v3 checkpoint lacks pooled calibration: {path}"
            )
        records.append(
            {
                "fold": fold,
                "path": str(path.relative_to(out_dir)),
                "sha256": _file_sha256(path),
                "pooled_oof_calibration_sha256": checkpoint[
                    "pooled_oof_calibration_sha256"
                ],
            }
        )
    payload = {
        "schema_version": TEACHER_V3_CHECKPOINT_MANIFEST_SCHEMA,
        "run_id": TEACHER_V3_RUN_ID,
        "release_name": TEACHER_V3_RUN_NAME,
        "model_version": MODEL_VERSION,
        "selected_profile": selected_profile,
        "selection_policy": selection_policy,
        **dict(input_provenance),
        "checkpoints": records,
    }
    path = Path(out_dir) / str(filename)
    _write_json(path, payload)
    return path


def build_teacher_v3_metadata_provenance(
    *,
    training_table: Path,
) -> dict[str, str]:
    """Bind the native-independent metadata baseline to the frozen table."""

    not_applicable = hashlib.sha256(
        b"teacher_v3_metadata_only_no_native_input_v1"
    ).hexdigest()
    return {
        "checkpoint_namespace": TEACHER_V3_METADATA_CHECKPOINT_NAMESPACE,
        "input_contract_version": "teacher_v3_metadata_only_v1",
        "native_h5_sha256": not_applicable,
        "native_manifest_sha256": not_applicable,
        "native_input_mode": "metadata_only_no_hdf5",
        "training_table_sha256": _file_sha256(training_table),
    }


def run_teacher_v3_metadata_baseline(
    *,
    training_table: Path,
    out_dir: Path,
    train_config: HarmonicTrainConfig = HarmonicTrainConfig(seed=560062),
    workers: int = 8,
    require_cuda: bool = True,
) -> dict[str, Any]:
    """Train a genuine seven-sector metadata baseline without native HDF5s."""

    training_table = Path(training_table).resolve()
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    source = pd.read_csv(training_table, low_memory=False)
    validate_teacher_v3_training_table(source)
    rows = prepare_harmonic_training_rows(
        source,
        seed=train_config.seed,
    )
    expected_active = int(
        _truth(source["teacher_v3_training_include"]).sum()
    )
    if len(rows) != expected_active:
        raise RuntimeError(
            "Teacher v3 active-row contract disagrees with model targets: "
            f"{len(rows)} != {expected_active}"
        )
    provenance = build_teacher_v3_metadata_provenance(
        training_table=training_table,
    )
    rows.to_csv(out_dir / "training_rows_with_fixed_splits.csv", index=False)
    injection_audit = injection_truth_human_audit(
        source,
        out_dir=out_dir / "injection_truth_human_audit",
    )
    for fold in range(5):
        _train_one_fold(
            rows=rows,
            native_h5=None,
            out_dir=out_dir,
            profile=TEACHER_V3_BASELINE_PROFILE,
            fold=fold,
            train_config=train_config,
            workers=workers,
            pretrain_epochs=0,
            require_cuda=require_cuda,
            input_provenance=provenance,
            defer_temperature_calibration=True,
            pretraining_cache_namespace="not_applicable_metadata_only",
            skip_encoder_pretraining=True,
        )
    calibration = calibrate_teacher_v3_profile_oof(
        rows=rows,
        out_dir=out_dir,
        profile=TEACHER_V3_BASELINE_PROFILE,
        input_provenance=provenance,
    )
    development = _profile_development_record(calibration)
    pd.DataFrame([development]).to_csv(
        out_dir / "development_fixed_profile_comparison.csv",
        index=False,
    )
    freeze = {
        "schema_version": "twirl_teacher_v3_model_freeze_v1",
        "run_id": TEACHER_V3_RUN_ID,
        "model_version": MODEL_VERSION,
        "profile": TEACHER_V3_BASELINE_PROFILE,
        "fixed_role": "metadata_baseline",
        "native_inputs_used": False,
        "encoder_pretraining_used": False,
        "temperature_calibration": (
            "one_temperature_from_concatenated_five_fold_"
            "development_oof_logits"
        ),
        "test_rows_used_for_selection_or_calibration": False,
        **provenance,
    }
    freeze_path = out_dir / "pretest_model_freeze.json"
    freeze_sha256 = _write_json(freeze_path, freeze)
    if _file_sha256(training_table) != provenance[
        "training_table_sha256"
    ]:
        raise RuntimeError(
            "Teacher v3 training table changed during metadata training"
        )
    manifest = _write_selected_checkpoint_manifest(
        out_dir=out_dir,
        input_provenance=provenance,
        selected_profile=TEACHER_V3_BASELINE_PROFILE,
        filename="metadata_baseline_checkpoint_manifest.json",
        selection_policy=(
            "fixed_metadata_baseline; fixed_test_sealed_until_primary_freeze"
        ),
    )
    summary: dict[str, Any] = {
        "schema_version": (
            "twirl_teacher_v3_metadata_baseline_training_summary_v1"
        ),
        "run_id": TEACHER_V3_RUN_ID,
        "release_name": TEACHER_V3_RUN_NAME,
        "model_version": MODEL_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "profiles": [TEACHER_V3_BASELINE_PROFILE],
        "fixed_role": "metadata_baseline",
        "native_inputs_used": False,
        "encoder_pretraining_used": False,
        "train_config": asdict(train_config),
        "n_training_rows": int(len(rows)),
        "n_development_rows": int(
            rows["fixed_split"].eq("development").sum()
        ),
        "n_fixed_test_rows": int(rows["fixed_split"].eq("test").sum()),
        **provenance,
        "calibration": calibration,
        "development_fixed_profile_comparison": [development],
        "fixed_test_status": "sealed_pending_primary_profile_freeze",
        "test_metrics": {},
        "pretest_model_freeze": str(freeze_path),
        "pretest_model_freeze_sha256": freeze_sha256,
        "metadata_baseline_checkpoint_manifest": str(manifest),
        "metadata_baseline_checkpoint_manifest_sha256": _file_sha256(
            manifest
        ),
        "injection_truth_human_audit": injection_audit,
        "student_training_blocked": True,
        "automatic_production_promotion": False,
    }
    _write_json(out_dir / "summary.json", summary)
    return summary


def _assert_native_inputs_unchanged(
    manifest: Mapping[str, Any],
) -> None:
    for record in manifest["native_files"]:
        path = Path(record["path"])
        if (
            path.stat().st_size != int(record["size_bytes"])
            or path.stat().st_mtime_ns != int(record["mtime_ns"])
            or file_sha256(path) != str(record["sha256"])
        ):
            raise RuntimeError(
                f"Teacher v3 native input changed during training: {path}"
            )
    registry = manifest["registry"]
    for path_name, digest_name in (
        ("path", "sha256"),
        ("summary_path", "summary_sha256"),
    ):
        path = Path(registry[path_name])
        if _file_sha256(path) != str(registry[digest_name]):
            raise RuntimeError(
                f"Teacher v3 native registry input changed: {path}"
            )


def run_teacher_v3_training(
    *,
    training_table: Path,
    native_registry: Path,
    native_registry_summary: Path,
    out_dir: Path,
    train_config: HarmonicTrainConfig = HarmonicTrainConfig(seed=560062),
    workers: int = 8,
    pretrain_epochs: int = 20,
    require_cuda: bool = True,
) -> dict[str, Any]:
    """Train the fixed Teacher v3 primary and baseline, then open the test."""

    training_table = Path(training_table).resolve()
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    source = pd.read_csv(training_table, low_memory=False)
    validate_teacher_v3_training_table(source)
    rows = prepare_harmonic_training_rows(
        source,
        seed=train_config.seed,
    )
    expected_active = int(
        _truth(source["teacher_v3_training_include"]).sum()
    )
    if len(rows) != expected_active:
        raise RuntimeError(
            "Teacher v3 active-row contract disagrees with model targets: "
            f"{len(rows)} != {expected_active}"
        )
    native_manifest = build_teacher_v3_native_manifest(
        rows=rows,
        registry_path=native_registry,
        registry_summary_path=native_registry_summary,
    )
    native_manifest_path = out_dir / "native_input_manifest.json"
    native_manifest_sha256 = _write_json(
        native_manifest_path,
        native_manifest,
    )
    input_provenance = build_teacher_v3_input_provenance(
        training_table=training_table,
        native_manifest=native_manifest,
    )
    if input_provenance["native_manifest_sha256"] != (
        native_manifest_sha256
    ):
        raise RuntimeError("Teacher v3 native manifest serialization changed")
    rows.to_csv(out_dir / "training_rows_with_fixed_splits.csv", index=False)
    injection_audit = injection_truth_human_audit(
        source,
        out_dir=out_dir / "injection_truth_human_audit",
    )

    fold_results: list[dict[str, Any]] = []
    calibration_results: list[dict[str, Any]] = []
    for profile in TEACHER_V3_PROFILES:
        metadata_only = profile == TEACHER_V3_BASELINE_PROFILE
        for fold in range(5):
            fold_results.append(
                _train_one_fold(
                    rows=rows,
                    native_h5=None,
                    out_dir=out_dir,
                    profile=profile,
                    fold=fold,
                    train_config=train_config,
                    workers=workers,
                    pretrain_epochs=(
                        0 if metadata_only else pretrain_epochs
                    ),
                    require_cuda=require_cuda,
                    input_provenance=input_provenance,
                    defer_temperature_calibration=True,
                    pretraining_cache_namespace=(
                        "teacher_v3_s56_s62_native_v1"
                    ),
                    skip_encoder_pretraining=metadata_only,
                )
            )
        calibration_results.append(
            calibrate_teacher_v3_profile_oof(
                rows=rows,
                out_dir=out_dir,
                profile=profile,
                input_provenance=input_provenance,
            )
        )

    development = pd.DataFrame(
        [
            _profile_development_record(result)
            for result in calibration_results
        ]
    )
    development.to_csv(
        out_dir / "development_fixed_profile_comparison.csv",
        index=False,
    )
    freeze = {
        "schema_version": "twirl_teacher_v3_model_freeze_v1",
        "run_id": TEACHER_V3_RUN_ID,
        "model_version": MODEL_VERSION,
        "primary_profile": TEACHER_V3_PRIMARY_PROFILE,
        "baseline_profile": TEACHER_V3_BASELINE_PROFILE,
        "profile_selection_policy": "fixed_before_test",
        "temperature_calibration": (
            "one_temperature_per_profile_from_concatenated_"
            "five_fold_development_oof_logits"
        ),
        "test_rows_used_for_selection_or_calibration": False,
        **input_provenance,
    }
    freeze_path = out_dir / "pretest_model_freeze.json"
    freeze_sha256 = _write_json(freeze_path, freeze)

    test_metrics = {
        profile: _evaluate_fixed_profile(
            rows=rows,
            out_dir=out_dir,
            profile=profile,
            train_config=train_config,
            workers=workers,
            input_provenance=input_provenance,
            require_cuda=require_cuda,
        )
        for profile in TEACHER_V3_PROFILES
    }
    _assert_native_inputs_unchanged(native_manifest)
    if _file_sha256(training_table) != input_provenance[
        "training_table_sha256"
    ]:
        raise RuntimeError(
            "Teacher v3 training table changed during training"
        )
    selected_manifest = _write_selected_checkpoint_manifest(
        out_dir=out_dir,
        input_provenance=input_provenance,
    )
    summary: dict[str, Any] = {
        "schema_version": "twirl_teacher_v3_training_summary_v1",
        "run_id": TEACHER_V3_RUN_ID,
        "release_name": TEACHER_V3_RUN_NAME,
        "model_version": MODEL_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "profiles": list(TEACHER_V3_PROFILES),
        "primary_profile": TEACHER_V3_PRIMARY_PROFILE,
        "baseline_profile": TEACHER_V3_BASELINE_PROFILE,
        "profile_selection_policy": "fixed_before_test",
        "train_config": asdict(train_config),
        "pretraining_epochs": int(pretrain_epochs),
        "n_training_rows": int(len(rows)),
        "n_development_rows": int(
            rows["fixed_split"].eq("development").sum()
        ),
        "n_fixed_test_rows": int(rows["fixed_split"].eq("test").sum()),
        "input_channel_contract": {
            name: list(channels)
            for name, channels in CHANNEL_CONTRACT.items()
        },
        **input_provenance,
        "native_input_manifest": str(native_manifest_path),
        "pretest_model_freeze": str(freeze_path),
        "pretest_model_freeze_sha256": freeze_sha256,
        "calibration": calibration_results,
        "development_fixed_profile_comparison": (
            development.to_dict("records")
        ),
        "test_metrics": test_metrics,
        "selected_checkpoint_manifest": str(selected_manifest),
        "selected_checkpoint_manifest_sha256": _file_sha256(
            selected_manifest
        ),
        "injection_truth_human_audit": injection_audit,
        "student_training_blocked": True,
        "automatic_production_promotion": False,
    }
    _write_json(out_dir / "summary.json", summary)
    return summary


__all__ = [
    "TEACHER_V3_BASELINE_PROFILE",
    "TEACHER_V3_CHECKPOINT_NAMESPACE",
    "TEACHER_V3_METADATA_CHECKPOINT_NAMESPACE",
    "TEACHER_V3_PRIMARY_PROFILE",
    "TEACHER_V3_PROFILES",
    "TEACHER_V3_RUN_ID",
    "build_teacher_v3_input_provenance",
    "build_teacher_v3_metadata_provenance",
    "build_teacher_v3_native_manifest",
    "calibrate_teacher_v3_profile_oof",
    "run_teacher_v3_training",
    "run_teacher_v3_metadata_baseline",
    "validate_teacher_v3_training_table",
]
