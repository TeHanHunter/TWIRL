"""Development-only evaluation of the completed Teacher v4-SSL encoders.

The frozen fixed-test population and prospective Sector 63 remain unopened.
Every supervised estimate is five-fold out-of-fold on the existing real
development labels.  Fold ``k`` always uses the full-pool encoder whose VICReg
pretraining excluded the complete TIC groups assigned to fold ``k``.
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from twirl.vetting.harmonic_cnn import (
    MODEL_VERSION,
    MORPHOLOGY_CLASSES,
    HarmonicModelConfig,
    build_harmonic_cnn,
)
from twirl.vetting.harmonic_dataset import HarmonicNativeDataset
from twirl.vetting.harmonic_training import (
    _file_sha256,
    _loader,
    _softmax,
    classification_metrics,
    expected_calibration_error,
    fit_temperature,
)
from twirl.vetting.teacher_ssl_fullpool import (
    FULLPOOL_SSL_CHECKPOINT_NAMESPACE,
    FULLPOOL_SSL_CHECKPOINT_SCHEMA,
    FULLPOOL_SSL_ENCODER_NAME,
    FULLPOOL_SSL_MODEL_FACING_NAME,
    FULLPOOL_SSL_PROFILE,
    FULLPOOL_SSL_RUN_ID,
    load_fullpool_ssl_registry,
)
from twirl.vetting.teacher_ssl_training import (
    TEACHER_V3_EXPECTED_SUMMARY_SHA256,
    _active_release_rows,
    _assert_development_only,
    _evaluate_loader_with_embeddings,
    _pool_development_oof,
    _run_finetune_fold,
)
from twirl.vetting.teacher_v3_training import (
    _read_training_table_with_stable_hash,
    _truth,
    validate_teacher_v3_training_table,
)


FULLPOOL_EVALUATION_SCHEMA = "twirl_teacher_ssl_fullpool_evaluation_v1"
FULLPOOL_EVALUATION_CHECKPOINT_SCHEMA = (
    "twirl_teacher_ssl_fullpool_finetuned_checkpoint_v1"
)
FULLPOOL_EVALUATION_OOF_SCHEMA = (
    "twirl_teacher_ssl_fullpool_development_oof_v1"
)
FULLPOOL_EVALUATION_RUN_ID = (
    "teacher_ssl_fullpool_v4_development_transfer_evaluation_v1"
)
FULLPOOL_EVALUATION_NAMESPACE = (
    "twirl_teacher_ssl_fullpool_v4_development_transfer_v1"
)
FULLPOOL_COMPLETION_SCHEMA = (
    "twirl_teacher_ssl_fullpool_v4_five_fold_completion_release_v1"
)
FULLPOOL_COMPLETION_BINDING = (
    "teacher_ssl_fullpool_v4_detector_consistent_raw_v1_"
    "effective_quality_adp_btjd_five_fold_complete_v1"
)
EXPECTED_REAL_DEVELOPMENT_ROWS = 6_168
EXPECTED_REAL_DEVELOPMENT_TICS = 6_054
EXPECTED_REAL_DEVELOPMENT_UNCERTAIN_ROWS = 3_780
EXPECTED_MATCHED_DEVELOPMENT_ROWS = 5_560
EXPECTED_MATCHED_DEVELOPMENT_TICS = 5_446
EXPECTED_MATCHED_DEVELOPMENT_UNCERTAIN_ROWS = 3_480
EXPECTED_S63_RESERVED_LABEL_ROWS = 608
EXPECTED_S63_RESERVED_LABEL_TICS = 608
EXPECTED_S63_RESERVED_LABEL_SECTORS = {61: 206, 62: 402}
EXPECTED_FULLPOOL_RAW_ROWS = 212_049
EXPECTED_FULLPOOL_RETAINED_ROWS = 175_366
EXPECTED_FULLPOOL_EXCLUDED_ROWS = 36_683
EXPECTED_S63_RESERVED_TICS = 53_512
EXPECTED_DUPLICATE_TARGET_SECTOR_LABEL_ROWS = 4
EXPECTED_DUPLICATE_TARGET_SECTOR_KEYS = 2
LABEL_POLICIES: tuple[str, ...] = (
    "uncertain_as_other",
    "uncertain_masked",
)


def _json_safe(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return [_json_safe(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, float):
        return value if np.isfinite(value) else None
    if isinstance(value, Path):
        return str(value)
    return value


def _write_json(path: Path, payload: Mapping[str, Any]) -> str:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            _json_safe(payload),
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
        + "\n",
        encoding="utf-8",
    )
    return _file_sha256(path)


def _code_revision() -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=Path(__file__).resolve().parents[3],
        check=False,
        capture_output=True,
        text=True,
    )
    revision = completed.stdout.strip()
    if completed.returncode != 0 or len(revision) != 40:
        raise RuntimeError("cannot bind full-pool evaluation to a Git revision")
    return revision


def _string_set_sha256(values: Sequence[str]) -> str:
    payload = "\n".join(sorted(str(value) for value in values)) + "\n"
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def require_fresh_evaluation_output(path: Path) -> Path:
    """Create one fresh evaluation root and reject ambiguous reuse."""

    output = Path(path).expanduser().resolve()
    if output.exists():
        raise FileExistsError(
            f"full-pool evaluation output must be fresh: {output}"
        )
    output.mkdir(parents=True, exist_ok=False)
    return output


def _read_json_object(path: Path, *, artifact: str) -> dict[str, Any]:
    path = Path(path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise RuntimeError(f"{artifact} must be a JSON object")
    return payload


def _validate_resume_run_contract(
    *,
    output_dir: Path,
    expected_contract: Mapping[str, Any],
    current_input_provenance: Mapping[str, str],
) -> dict[str, Any]:
    """Validate an interrupted evaluation root before any artifact reuse."""

    output = Path(output_dir).expanduser().resolve()
    if not output.is_dir():
        raise FileNotFoundError(output)
    for completed in (
        output / "summary.json",
        output / "development_performance.csv",
        output / "development_per_class_metrics.csv",
    ):
        if completed.exists():
            raise RuntimeError(
                "refusing to resume an evaluation with final artifacts: "
                f"{completed}"
            )
    contract_path = output / "run_contract.json"
    contract = _read_json_object(
        contract_path,
        artifact="interrupted evaluation run contract",
    )
    fields = (
        "schema_version",
        "run_id",
        "model_facing_name",
        "scope",
        "label_policies",
        "fold_encoder_rule",
        "fine_tune_epochs",
        "batch_size",
        "workers",
        "seed",
        "linear_probe_steps",
        "bootstrap_draws",
        "fixed_test_container_opened",
        "fixed_test_tensors_constructed",
        "fixed_test_labels_used",
        "sector_63_rows_loaded",
        "injections_loaded",
        "automatic_production_promotion",
    )
    changed = [
        name
        for name in fields
        if contract.get(name) != expected_contract.get(name)
    ]
    if changed:
        raise RuntimeError(
            "interrupted evaluation contract changed: "
            f"{changed}"
        )
    provenance_fields = (
        "training_table_sha256",
        "fullpool_registry_sha256",
        "fullpool_registry_summary_sha256",
        "fullpool_completion_release_sha256",
        "fullpool_training_code_revision",
    )
    provenance_changed = [
        name
        for name in provenance_fields
        if contract.get(name) != current_input_provenance.get(name)
    ]
    if provenance_changed:
        raise RuntimeError(
            "interrupted evaluation input provenance changed: "
            f"{provenance_changed}"
        )
    original_revision = str(contract.get("evaluation_code_revision", ""))
    if len(original_revision) != 40:
        raise RuntimeError("interrupted evaluation lacks a Git revision")
    return {
        "mode": "contract_checked_resume",
        "interrupted_run_contract": str(contract_path),
        "interrupted_run_contract_sha256": _file_sha256(contract_path),
        "interrupted_evaluation_code_revision": original_revision,
        "resume_evaluation_code_revision": str(
            current_input_provenance["evaluation_code_revision"]
        ),
        "fixed_test_status": "sealed_not_loaded_not_evaluated",
        "sector_63_status": "reserved_not_loaded",
    }


def _read_bound_json(
    binding: Mapping[str, Any],
    *,
    artifact: str,
) -> tuple[dict[str, Any], Path, str]:
    if not isinstance(binding, Mapping) or not {"path", "sha256"}.issubset(
        binding
    ):
        raise RuntimeError(f"{artifact} binding is incomplete")
    path = Path(str(binding["path"])).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    if "size_bytes" in binding and path.stat().st_size != int(
        binding["size_bytes"]
    ):
        raise RuntimeError(f"{artifact} size changed")
    observed_sha256 = _file_sha256(path)
    if observed_sha256 != str(binding["sha256"]):
        raise RuntimeError(f"{artifact} hash changed")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise RuntimeError(f"{artifact} must be a JSON object")
    return payload, path, observed_sha256


def _validate_preregistered_label_exclusions(
    *,
    excluded: pd.DataFrame,
    registry_summary: Mapping[str, Any],
) -> dict[str, Any]:
    """Prove unmatched labels are raw inputs held out for prospective S63."""

    if (
        len(excluded) != EXPECTED_S63_RESERVED_LABEL_ROWS
        or excluded["tic"].nunique() != EXPECTED_S63_RESERVED_LABEL_TICS
    ):
        raise RuntimeError(
            "unmatched development-label population changed: "
            f"rows={len(excluded)}, tics={excluded['tic'].nunique()}"
        )
    sector_counts = {
        int(sector): int(count)
        for sector, count in excluded["sector"].value_counts().sort_index().items()
    }
    if sector_counts != EXPECTED_S63_RESERVED_LABEL_SECTORS:
        raise RuntimeError(
            "unmatched development-label sectors changed: "
            f"{sector_counts} != {EXPECTED_S63_RESERVED_LABEL_SECTORS}"
        )
    if excluded["tic"].duplicated().any():
        raise RuntimeError("an unmatched development TIC has multiple label rows")

    provenance = registry_summary.get("source_provenance")
    if not isinstance(provenance, Mapping):
        raise RuntimeError("SSL registry lacks source provenance")
    fullpool_summary, fullpool_summary_path, fullpool_summary_sha256 = (
        _read_bound_json(
            provenance.get("frozen_pool_summary", {}),
            artifact="frozen full-pool summary",
        )
    )
    counts = fullpool_summary.get("counts", {})
    observed_pool_counts = {
        "input": int(counts.get("input", {}).get("n_observations", -1)),
        "retained": int(counts.get("retained", {}).get("n_observations", -1)),
        "excluded": int(
            counts.get("excluded_input", {}).get("n_observations", -1)
        ),
    }
    expected_pool_counts = {
        "input": EXPECTED_FULLPOOL_RAW_ROWS,
        "retained": EXPECTED_FULLPOOL_RETAINED_ROWS,
        "excluded": EXPECTED_FULLPOOL_EXCLUDED_ROWS,
    }
    if observed_pool_counts != expected_pool_counts:
        raise RuntimeError(
            "frozen full-pool inventory counts changed: "
            f"{observed_pool_counts} != {expected_pool_counts}"
        )
    if fullpool_summary.get("exclusion_scope") != (
        "whole host: exclude every S56--S62 observation of any TIC in the "
        "frozen fixed-test registry or prospective S63 reservation"
    ):
        raise RuntimeError("frozen full-pool whole-host exclusion scope changed")

    summary_inputs = fullpool_summary.get("inputs")
    if not isinstance(summary_inputs, Mapping):
        raise RuntimeError("frozen full-pool summary lacks input bindings")
    reserved_binding = summary_inputs.get("s63_reserved_tics")
    reserved_hosts = provenance.get("reserved_hosts")
    if not isinstance(reserved_binding, Mapping) or not isinstance(
        reserved_hosts, Mapping
    ):
        raise RuntimeError("prospective S63 host authority is missing")
    for field in ("path", "sha256", "size_bytes"):
        if str(reserved_binding.get(field)) != str(reserved_hosts.get(field)):
            raise RuntimeError(
                f"prospective S63 host bindings disagree for {field}"
            )
    reserved_path = Path(str(reserved_binding["path"])).expanduser().resolve()
    if not reserved_path.is_file():
        raise FileNotFoundError(reserved_path)
    if reserved_path.stat().st_size != int(reserved_binding["size_bytes"]):
        raise RuntimeError("prospective S63 host inventory size changed")
    reserved_sha256 = _file_sha256(reserved_path)
    if reserved_sha256 != str(reserved_binding["sha256"]):
        raise RuntimeError("prospective S63 host inventory hash changed")
    reserved_tics = {
        int(value.strip())
        for value in reserved_path.read_text(encoding="utf-8").splitlines()
        if value.strip()
    }
    if len(reserved_tics) != EXPECTED_S63_RESERVED_TICS:
        raise RuntimeError(
            "prospective S63 TIC count changed: "
            f"{len(reserved_tics)} != {EXPECTED_S63_RESERVED_TICS}"
        )
    excluded_tics = set(
        pd.to_numeric(excluded["tic"], errors="raise").astype(np.int64)
    )
    if not excluded_tics.issubset(reserved_tics):
        unexpected = sorted(excluded_tics - reserved_tics)
        raise RuntimeError(
            "unmatched labels are not prospective S63 whole-host exclusions: "
            f"{unexpected[:10]}"
        )

    compact_exports = summary_inputs.get("compact_exports")
    if not isinstance(compact_exports, list):
        raise RuntimeError("frozen full-pool summary lacks compact exports")
    export_by_sector = {
        int(record["sector"]): record
        for record in compact_exports
        if isinstance(record, Mapping) and "sector" in record
    }
    raw_manifest_audit: dict[str, Any] = {}
    for sector, expected_count in sector_counts.items():
        export = export_by_sector.get(sector)
        if not isinstance(export, Mapping):
            raise RuntimeError(f"no compact-export authority for Sector {sector}")
        manifest, manifest_path, manifest_sha256 = _read_bound_json(
            export.get("manifest", {}),
            artifact=f"Sector {sector} compact manifest",
        )
        records = manifest.get("records")
        if not isinstance(records, list) or len(records) != int(
            export.get("n_exported_targets", -1)
        ):
            raise RuntimeError(f"Sector {sector} compact manifest count changed")
        raw_tics = {
            int(record["tic"])
            for record in records
            if isinstance(record, Mapping) and "tic" in record
        }
        expected_tics = set(
            pd.to_numeric(
                excluded.loc[excluded["sector"].eq(sector), "tic"],
                errors="raise",
            ).astype(np.int64)
        )
        if len(expected_tics) != expected_count or not expected_tics.issubset(
            raw_tics
        ):
            absent = sorted(expected_tics - raw_tics)
            raise RuntimeError(
                f"Sector {sector} excluded labels are absent from the raw "
                f"compact manifest: {absent[:10]}"
            )
        raw_manifest_audit[str(sector)] = {
            "path": str(manifest_path),
            "sha256": manifest_sha256,
            "n_raw_observations": int(len(records)),
            "n_excluded_label_observations_present": int(len(expected_tics)),
        }

    return {
        "reason": "prospective_s63_whole_tic_holdout",
        "n_rows": int(len(excluded)),
        "n_tics": int(len(excluded_tics)),
        "sector_counts": {str(key): value for key, value in sector_counts.items()},
        "human_label_counts": {
            str(name): int(value)
            for name, value in excluded["human_label"]
            .astype(str)
            .value_counts()
            .sort_index()
            .items()
        },
        "review_ids_sha256": _string_set_sha256(
            excluded["review_id"].astype(str).tolist()
        ),
        "tics_sha256": _string_set_sha256(
            excluded["tic"].astype(str).tolist()
        ),
        "reserved_inventory": str(reserved_path),
        "reserved_inventory_sha256": reserved_sha256,
        "fullpool_summary": str(fullpool_summary_path),
        "fullpool_summary_sha256": fullpool_summary_sha256,
        "raw_compact_inventory_present": True,
        "raw_compact_manifests": raw_manifest_audit,
    }


def load_fullpool_development_labels(
    *,
    training_table_path: Path,
    registry_path: Path,
    registry_summary_path: Path,
    seed: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Bind real development labels to the exact full-pool native inputs."""

    training_table_path = Path(training_table_path).expanduser().resolve()
    source, training_sha256, stable_read = _read_training_table_with_stable_hash(
        training_table_path,
        artifact="Teacher v4-SSL evaluation label table",
    )
    validate_teacher_v3_training_table(source)
    active = _active_release_rows(source, seed=int(seed))
    real = ~_truth(active["is_injected_row"])
    development = active["fixed_split"].astype(str).eq("development")
    rows = active.loc[real & development].copy().reset_index(drop=True)
    _assert_development_only(
        rows,
        artifact="full-pool evaluation labels",
        real_only=True,
    )
    if len(rows) != EXPECTED_REAL_DEVELOPMENT_ROWS:
        raise RuntimeError(
            "real development-label count changed: "
            f"{len(rows)} != {EXPECTED_REAL_DEVELOPMENT_ROWS}"
        )
    if rows["tic"].nunique() != EXPECTED_REAL_DEVELOPMENT_TICS:
        raise RuntimeError(
            "real development TIC count changed: "
            f"{rows['tic'].nunique()} != {EXPECTED_REAL_DEVELOPMENT_TICS}"
        )
    uncertain = rows["human_label"].astype(str).eq("uncertain")
    if int(uncertain.sum()) != EXPECTED_REAL_DEVELOPMENT_UNCERTAIN_ROWS:
        raise RuntimeError(
            "real development uncertain-label count changed: "
            f"{int(uncertain.sum())} != "
            f"{EXPECTED_REAL_DEVELOPMENT_UNCERTAIN_ROWS}"
        )
    duplicate_input = rows.duplicated(["sector", "tic"], keep=False)
    duplicate_rows = rows.loc[duplicate_input].copy()
    duplicate_keys = duplicate_rows.loc[:, ["sector", "tic"]].drop_duplicates()
    if (
        len(duplicate_rows) != EXPECTED_DUPLICATE_TARGET_SECTOR_LABEL_ROWS
        or len(duplicate_keys) != EXPECTED_DUPLICATE_TARGET_SECTOR_KEYS
    ):
        raise RuntimeError(
            "candidate-label/full-pool-input multiplicity changed: "
            f"rows={len(duplicate_rows)}, keys={len(duplicate_keys)}"
        )
    candidate_ephemeris = ["period_d", "t0_bjd", "duration_min"]
    duplicate_ephemerides = duplicate_rows.loc[
        :, ["sector", "tic", *candidate_ephemeris]
    ].drop_duplicates()
    if len(duplicate_ephemerides) != len(duplicate_rows):
        raise RuntimeError(
            "duplicate target-sector label records do not have distinct "
            "candidate ephemerides"
        )

    registry, registry_summary = load_fullpool_ssl_registry(
        registry_path=registry_path,
        summary_path=registry_summary_path,
    )
    columns = [
        "sector",
        "tic",
        "ssl_observation_id",
        "period_d",
        "t0_bjd",
        "duration_min",
        "bls_status",
        "native_h5_path",
        "native_group_path",
        "native_h5_sha256",
        "native_contract_version",
        "ssl_pool_include",
        "fixed_test_member",
        "reserved_prospective_member",
        "ssl_held_out_fold",
        "fold_assignment_source",
    ]
    registry_labels = registry.loc[:, columns].rename(
        columns={
            name: f"fullpool_{name}"
            for name in columns
            if name not in {"sector", "tic"}
        }
    )
    merged = rows.merge(
        registry_labels,
        on=["sector", "tic"],
        how="left",
        validate="many_to_one",
        indicator=True,
    )
    excluded = merged.loc[merged["_merge"].ne("both")].copy()
    exclusion_audit: dict[str, Any] | None = None
    if len(excluded):
        exclusion_audit = _validate_preregistered_label_exclusions(
            excluded=excluded,
            registry_summary=registry_summary,
        )
    merged = merged.loc[merged["_merge"].eq("both")].drop(columns="_merge")
    if (
        len(merged) != EXPECTED_MATCHED_DEVELOPMENT_ROWS
        or merged["tic"].nunique() != EXPECTED_MATCHED_DEVELOPMENT_TICS
    ):
        raise RuntimeError(
            "retained full-pool development support changed: "
            f"rows={len(merged)}, tics={merged['tic'].nunique()}"
        )
    matched_uncertain = merged["human_label"].astype(str).eq("uncertain")
    if int(matched_uncertain.sum()) != EXPECTED_MATCHED_DEVELOPMENT_UNCERTAIN_ROWS:
        raise RuntimeError(
            "retained full-pool uncertain-label count changed: "
            f"{int(matched_uncertain.sum())} != "
            f"{EXPECTED_MATCHED_DEVELOPMENT_UNCERTAIN_ROWS}"
        )
    for column in (
        "ssl_pool_include",
        "fixed_test_member",
        "reserved_prospective_member",
    ):
        merged[f"fullpool_{column}"] = merged[
            f"fullpool_{column}"
        ].astype(bool)
    if not merged["fullpool_ssl_pool_include"].all():
        raise RuntimeError("a labeled development row is excluded from SSL")
    if merged["fullpool_fixed_test_member"].any():
        raise RuntimeError("fixed-test membership reached development evaluation")
    if merged["fullpool_reserved_prospective_member"].any():
        raise RuntimeError("prospective membership reached development evaluation")
    if not merged["fullpool_fold_assignment_source"].eq(
        "frozen_development_split"
    ).all():
        raise RuntimeError("development labels lack frozen fold assignments")
    if not np.array_equal(
        pd.to_numeric(merged["cv_fold"], errors="raise").to_numpy(dtype=int),
        pd.to_numeric(
            merged["fullpool_ssl_held_out_fold"], errors="raise"
        ).to_numpy(dtype=int),
    ):
        raise RuntimeError(
            "Teacher-v3 development folds disagree with SSL held-out folds"
        )

    # Preserve the candidate-specific frozen ephemeris.  Two S56 target-sector
    # light curves each have two reviewed candidate periods, including one
    # genuine label disagreement.  They share one raw native light curve but
    # must remain distinct folded model inputs for matched Teacher-v3 support.
    for name in ("native_h5_path", "native_group_path"):
        merged[f"teacher_v3_{name}"] = merged[name]
        merged[name] = merged[f"fullpool_{name}"]
    merged["ssl_observation_id"] = merged["fullpool_ssl_observation_id"]
    merged["native_h5_sha256"] = merged["fullpool_native_h5_sha256"]
    merged["native_contract_version"] = merged[
        "fullpool_native_contract_version"
    ]
    merged["input_variant"] = "observed"
    merged = merged.sort_values(["cv_fold", "sector", "tic"], kind="stable")
    merged = merged.reset_index(drop=True)
    _assert_development_only(
        merged,
        artifact="bound full-pool evaluation labels",
        real_only=True,
    )
    audit = {
        "scope": "real_development_labels_only",
        "n_rows": int(len(merged)),
        "n_tics": int(merged["tic"].nunique()),
        "n_uncertain_rows": int(matched_uncertain.sum()),
        "n_source_development_rows": int(len(rows)),
        "n_source_development_tics": int(rows["tic"].nunique()),
        "n_source_uncertain_rows": int(uncertain.sum()),
        "support_rule": (
            "intersection_with_retained_fullpool_ssl_registry_after_"
            "preregistered_whole_tic_exclusions"
        ),
        "excluded_label_authority": exclusion_audit,
        "evaluation_unit": "candidate_review_record",
        "n_unique_target_sector_native_inputs": int(
            merged.loc[:, ["sector", "tic"]].drop_duplicates().shape[0]
        ),
        "n_duplicate_target_sector_label_rows": int(len(duplicate_rows)),
        "n_duplicate_target_sector_keys": int(len(duplicate_keys)),
        "duplicate_candidate_ephemerides_distinct": True,
        "candidate_ephemeris_source": (
            "frozen_teacher_v3_review_record; fullpool_registry_ephemeris_"
            "retained_as_audit_columns_only"
        ),
        "review_ids_sha256": _string_set_sha256(
            merged["review_id"].astype(str).tolist()
        ),
        "ssl_observation_ids_sha256": _string_set_sha256(
            merged["ssl_observation_id"].astype(str).tolist()
        ),
        "training_table": str(training_table_path),
        "training_table_sha256": training_sha256,
        "training_table_initial_read": stable_read,
        "registry": str(Path(registry_path).expanduser().resolve()),
        "registry_sha256": _file_sha256(registry_path),
        "registry_summary": str(
            Path(registry_summary_path).expanduser().resolve()
        ),
        "registry_summary_sha256": _file_sha256(registry_summary_path),
        "registry_summary_schema": registry_summary.get(
            "summary_schema_version"
        ),
        "fixed_test_tensors_constructed": False,
        "fixed_test_labels_used": False,
        "sector_63_rows_present": False,
        "injected_rows_present": False,
    }
    return merged, audit


def load_completed_fullpool_encoders(
    completion_release_path: Path,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Load and checksum-verify the five immutable production encoders."""

    import torch

    release_path = Path(completion_release_path).expanduser().resolve()
    release = json.loads(release_path.read_text(encoding="utf-8"))
    if not isinstance(release, dict):
        raise ValueError("completion release must be a JSON object")
    expected = {
        "passed": True,
        "schema_version": FULLPOOL_COMPLETION_SCHEMA,
        "release_binding": FULLPOOL_COMPLETION_BINDING,
        "run_id": FULLPOOL_SSL_RUN_ID,
        "encoder_name": FULLPOOL_SSL_ENCODER_NAME,
        "checkpoint_namespace": FULLPOOL_SSL_CHECKPOINT_NAMESPACE,
    }
    failures = [
        f"{name}={release.get(name)!r}, expected {value!r}"
        for name, value in expected.items()
        if release.get(name) != value
    ]
    if release.get("counts") != {"completed_epochs": 100, "folds": 5}:
        failures.append("completion counts are not five folds x 20 epochs")
    folds = release.get("folds")
    if not isinstance(folds, list) or len(folds) != 5:
        failures.append("completion release lacks exactly five fold records")
        folds = []
    if failures:
        raise ValueError("invalid full-pool completion release: " + "; ".join(failures))

    loaded: list[dict[str, Any]] = []
    for expected_fold, record in enumerate(folds):
        if not isinstance(record, Mapping) or int(record.get("fold", -1)) != expected_fold:
            raise ValueError("completion folds are not ordered 0..4")
        binding = record.get("checkpoint")
        if not isinstance(binding, Mapping) or set(binding) != {
            "path",
            "sha256",
            "size_bytes",
        }:
            raise ValueError(f"fold {expected_fold} checkpoint binding is invalid")
        checkpoint_path = Path(str(binding["path"])).resolve()
        if not checkpoint_path.is_file():
            raise FileNotFoundError(checkpoint_path)
        stat = checkpoint_path.stat()
        if int(stat.st_size) != int(binding["size_bytes"]):
            raise ValueError(f"fold {expected_fold} checkpoint size changed")
        observed_sha256 = _file_sha256(checkpoint_path)
        if observed_sha256 != str(binding["sha256"]):
            raise ValueError(f"fold {expected_fold} checkpoint hash changed")
        checkpoint = torch.load(
            checkpoint_path,
            map_location="cpu",
            weights_only=False,
        )
        checkpoint_expected = {
            "schema_version": FULLPOOL_SSL_CHECKPOINT_SCHEMA,
            "run_id": FULLPOOL_SSL_RUN_ID,
            "encoder_name": FULLPOOL_SSL_ENCODER_NAME,
            "checkpoint_namespace": FULLPOOL_SSL_CHECKPOINT_NAMESPACE,
            "profile": FULLPOOL_SSL_PROFILE,
            "fold": expected_fold,
            "epochs": 20,
            "model_version": MODEL_VERSION,
        }
        mismatch = [
            name
            for name, value in checkpoint_expected.items()
            if checkpoint.get(name) != value
        ]
        if mismatch:
            raise ValueError(
                f"fold {expected_fold} checkpoint fields changed: {mismatch}"
            )
        config = checkpoint.get("model_config")
        if not isinstance(config, Mapping) or int(config.get("metadata_dim", -1)) != 0:
            raise ValueError(f"fold {expected_fold} is not a metadata-free encoder")
        state = checkpoint.get("encoder_state_dict")
        if not isinstance(state, Mapping) or not state:
            raise ValueError(f"fold {expected_fold} lacks encoder weights")
        loaded.append(
            {
                "fold": expected_fold,
                "checkpoint": str(checkpoint_path),
                "checkpoint_sha256": observed_sha256,
                "cache_identity": {
                    "source": "frozen_fullpool_completion_release",
                    "completion_release": str(release_path),
                    "completion_release_sha256": _file_sha256(release_path),
                    "fold": expected_fold,
                    "fullpool_training_code_revision": release[
                        "expected_code_revision"
                    ],
                },
                "encoder_state_dict": state,
            }
        )
    audit = {
        "completion_release": str(release_path),
        "completion_release_sha256": _file_sha256(release_path),
        "completion_schema": release["schema_version"],
        "completion_binding": release["release_binding"],
        "fullpool_training_code_revision": release["expected_code_revision"],
        "checkpoint_sha256": {
            str(item["fold"]): item["checkpoint_sha256"] for item in loaded
        },
    }
    return loaded, audit


def _extract_encoder_embeddings(
    *,
    rows: pd.DataFrame,
    encoder_state: Mapping[str, Any],
    batch_size: int,
    workers: int,
    seed: int,
    require_cuda: bool,
) -> np.ndarray:
    import torch

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if require_cuda and device.type != "cuda":
        raise RuntimeError("CUDA was required for full-pool embedding export")
    dataset = HarmonicNativeDataset(
        rows,
        native_h5=None,
        metadata=np.empty((len(rows), 0), dtype=np.float32),
        profile=FULLPOOL_SSL_PROFILE,
    )
    try:
        loader = _loader(
            dataset,
            np.arange(len(rows), dtype=int),
            batch_size=int(batch_size),
            shuffle=False,
            workers=int(workers),
            seed=int(seed),
        )
        model = build_harmonic_cnn(
            HarmonicModelConfig(metadata_dim=0),
            profile=FULLPOOL_SSL_PROFILE,
        ).to(device)
        model.load_state_dict(dict(encoder_state), strict=True)
        evaluation = _evaluate_loader_with_embeddings(
            model,
            loader,
            device=device,
        )
    finally:
        dataset.close()
    expected_ids = rows["review_id"].astype(str).tolist()
    if evaluation["review_id"] != expected_ids:
        raise RuntimeError("embedding export changed development-row order")
    embedding = np.asarray(evaluation["embedding"], dtype=np.float32)
    if embedding.shape != (len(rows), 128) or not np.isfinite(embedding).all():
        raise RuntimeError(
            f"invalid encoder embedding matrix: {embedding.shape}"
        )
    return embedding


def _embedding_diagnostics(values: np.ndarray) -> dict[str, Any]:
    values = np.asarray(values, dtype=np.float64)
    centered = values - values.mean(axis=0, keepdims=True)
    standard_deviation = centered.std(axis=0, ddof=1)
    singular = np.linalg.svd(centered, full_matrices=False, compute_uv=False)
    variance = singular**2
    probability = variance / max(float(variance.sum()), 1.0e-12)
    active = probability > 0
    effective_rank = float(
        np.exp(-np.sum(probability[active] * np.log(probability[active])))
    )
    return {
        "n_rows": int(len(values)),
        "embedding_dim": int(values.shape[1]),
        "mean_dimension_std": float(standard_deviation.mean()),
        "minimum_dimension_std": float(standard_deviation.min()),
        "maximum_dimension_std": float(standard_deviation.max()),
        "effective_rank": effective_rank,
        "mean_l2_norm": float(np.linalg.norm(values, axis=1).mean()),
        "finite": True,
    }


def export_all_fold_embeddings(
    *,
    rows: pd.DataFrame,
    encoders: Sequence[Mapping[str, Any]],
    output_dir: Path,
    batch_size: int,
    workers: int,
    seed: int,
    require_cuda: bool,
) -> tuple[list[np.ndarray], dict[str, Any]]:
    """Export coherent reference-fold and fold-safe held-row embeddings."""

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    all_embeddings: list[np.ndarray] = []
    diagnostics: list[dict[str, Any]] = []
    for record in encoders:
        fold = int(record["fold"])
        embedding = _extract_encoder_embeddings(
            rows=rows,
            encoder_state=record["encoder_state_dict"],
            batch_size=int(batch_size),
            workers=int(workers),
            seed=int(seed) + 1000 * fold,
            require_cuda=bool(require_cuda),
        )
        all_embeddings.append(embedding)
        diagnostics.append(
            {
                "fold": fold,
                **_embedding_diagnostics(embedding),
            }
        )
        print(
            f"[fullpool-eval] exported encoder fold {fold} embeddings",
            flush=True,
        )

    reference_path = output_dir / "reference_fold_0_embeddings.npz"
    with reference_path.open("wb") as stream:
        np.savez_compressed(
            stream,
            review_id=rows["review_id"].astype(str).to_numpy(dtype=str),
            ssl_observation_id=rows["ssl_observation_id"].astype(str).to_numpy(
                dtype=str
            ),
            tic=pd.to_numeric(rows["tic"], errors="raise").to_numpy(dtype=np.int64),
            sector=pd.to_numeric(rows["sector"], errors="raise").to_numpy(dtype=np.int16),
            cv_fold=pd.to_numeric(rows["cv_fold"], errors="raise").to_numpy(dtype=np.int8),
            human_label=rows["human_label"].astype(str).to_numpy(dtype=str),
            morphology_target=rows["morphology_target_index"].to_numpy(dtype=np.int8),
            period_d=pd.to_numeric(rows["period_d"], errors="raise").to_numpy(dtype=np.float64),
            duration_min=pd.to_numeric(
                rows["duration_min"], errors="raise"
            ).to_numpy(dtype=np.float64),
            embedding=all_embeddings[0],
        )
    held_embedding = np.empty((len(rows), 128), dtype=np.float32)
    for fold, embedding in enumerate(all_embeddings):
        mask = rows["cv_fold"].to_numpy(dtype=int) == fold
        held_embedding[mask] = embedding[mask]
    held_path = output_dir / "fold_safe_held_embeddings.npz"
    with held_path.open("wb") as stream:
        np.savez_compressed(
            stream,
            review_id=rows["review_id"].astype(str).to_numpy(dtype=str),
            tic=pd.to_numeric(rows["tic"], errors="raise").to_numpy(dtype=np.int64),
            cv_fold=pd.to_numeric(rows["cv_fold"], errors="raise").to_numpy(dtype=np.int8),
            embedding=held_embedding,
        )
    audit = {
        "reference_encoder_fold": 0,
        "reference_scope": (
            "exploratory_coherent_embedding_space; labels_not_used_in_"
            "pretraining; not_a_performance_estimate"
        ),
        "fold_safe_scope": (
            "each_row_from_encoder_that_excluded_its_complete_tic_fold"
        ),
        "reference_embeddings": str(reference_path),
        "reference_embeddings_sha256": _file_sha256(reference_path),
        "fold_safe_embeddings": str(held_path),
        "fold_safe_embeddings_sha256": _file_sha256(held_path),
        "diagnostics": diagnostics,
    }
    _write_json(output_dir / "embedding_export.summary.json", audit)
    return all_embeddings, audit


def validate_and_recompute_existing_embeddings(
    *,
    rows: pd.DataFrame,
    encoders: Sequence[Mapping[str, Any]],
    output_dir: Path,
    batch_size: int,
    workers: int,
    seed: int,
    require_cuda: bool,
) -> tuple[list[np.ndarray], dict[str, Any]]:
    """Recompute fold embeddings and prove equality with persisted exports."""

    output_dir = Path(output_dir).expanduser().resolve()
    summary_path = output_dir / "embedding_export.summary.json"
    summary = _read_json_object(
        summary_path,
        artifact="interrupted embedding export summary",
    )
    reference_path = output_dir / "reference_fold_0_embeddings.npz"
    held_path = output_dir / "fold_safe_held_embeddings.npz"
    expected_bindings = (
        (
            "reference_embeddings",
            "reference_embeddings_sha256",
            reference_path,
        ),
        (
            "fold_safe_embeddings",
            "fold_safe_embeddings_sha256",
            held_path,
        ),
    )
    for path_field, hash_field, expected_path in expected_bindings:
        if Path(str(summary.get(path_field, ""))).resolve() != expected_path:
            raise RuntimeError(
                f"interrupted embedding path changed for {path_field}"
            )
        if _file_sha256(expected_path) != str(summary.get(hash_field, "")):
            raise RuntimeError(
                f"interrupted embedding hash changed for {path_field}"
            )

    all_embeddings: list[np.ndarray] = []
    diagnostics: list[dict[str, Any]] = []
    for record in encoders:
        fold = int(record["fold"])
        embedding = _extract_encoder_embeddings(
            rows=rows,
            encoder_state=record["encoder_state_dict"],
            batch_size=int(batch_size),
            workers=int(workers),
            seed=int(seed) + 1000 * fold,
            require_cuda=bool(require_cuda),
        )
        all_embeddings.append(embedding)
        diagnostics.append({"fold": fold, **_embedding_diagnostics(embedding)})
        print(
            f"[fullpool-eval-resume] recomputed encoder fold {fold} embeddings",
            flush=True,
        )

    # These files were created by the trusted, hash-bound interrupted run.
    # Its first revision stored pandas string arrays as object dtype, so pickle
    # is enabled only after validating the complete artifact hash above.
    with np.load(reference_path, allow_pickle=True) as payload:
        if payload["review_id"].astype(str).tolist() != rows[
            "review_id"
        ].astype(str).tolist():
            raise RuntimeError("reference embedding row identities changed")
        persisted_reference = payload["embedding"].astype(np.float32)
    if not np.array_equal(persisted_reference, all_embeddings[0]):
        raise RuntimeError("recomputed reference embeddings changed")

    recomputed_held = np.empty((len(rows), 128), dtype=np.float32)
    for fold, embedding in enumerate(all_embeddings):
        mask = rows["cv_fold"].to_numpy(dtype=int) == fold
        recomputed_held[mask] = embedding[mask]
    with np.load(held_path, allow_pickle=True) as payload:
        if payload["review_id"].astype(str).tolist() != rows[
            "review_id"
        ].astype(str).tolist():
            raise RuntimeError("held embedding row identities changed")
        persisted_held = payload["embedding"].astype(np.float32)
    if not np.array_equal(persisted_held, recomputed_held):
        raise RuntimeError("recomputed fold-safe embeddings changed")

    audit = dict(summary)
    audit.update(
        {
            "resume_reused": True,
            "resume_validation": (
                "all_five_encoders_recomputed; reference_and_fold_safe_"
                "exports_exactly_equal"
            ),
            "resume_diagnostics": diagnostics,
            "summary": str(summary_path),
            "summary_sha256": _file_sha256(summary_path),
        }
    )
    return all_embeddings, audit


def _policy_rows(rows: pd.DataFrame, policy: str) -> pd.DataFrame:
    if policy not in LABEL_POLICIES:
        raise ValueError(f"unknown label policy {policy!r}")
    if policy == "uncertain_as_other":
        selected = rows.copy()
    else:
        selected = rows.loc[
            ~rows["human_label"].astype(str).eq("uncertain")
        ].copy()
    selected = selected.reset_index(drop=True)
    _assert_development_only(
        selected,
        artifact=f"{policy} supervised evaluation",
        real_only=True,
    )
    if policy == "uncertain_masked" and selected["human_label"].astype(str).eq(
        "uncertain"
    ).any():
        raise RuntimeError("uncertain rows reached the masked evaluation")
    return selected


def _average_precision_binary(truth: np.ndarray, score: np.ndarray) -> float:
    truth = np.asarray(truth, dtype=bool)
    score = np.asarray(score, dtype=float)
    positive = int(truth.sum())
    if positive == 0:
        return float("nan")
    order = np.argsort(-score, kind="stable")
    ranked = truth[order]
    precision = np.cumsum(ranked) / np.arange(1, len(ranked) + 1)
    return float(precision[ranked].sum() / positive)


def probability_metrics(
    truth: np.ndarray,
    probability: np.ndarray,
) -> dict[str, Any]:
    """Return classification, calibration, and one-vs-rest AP metrics."""

    truth = np.asarray(truth, dtype=int)
    probability = np.asarray(probability, dtype=float)
    active = truth >= 0
    selected_truth = truth[active]
    selected_probability = probability[active]
    classification = classification_metrics(
        selected_truth,
        selected_probability,
        classes=MORPHOLOGY_CLASSES,
    )
    per_class_ap: dict[str, float] = {}
    for index, label in enumerate(MORPHOLOGY_CLASSES):
        per_class_ap[label] = _average_precision_binary(
            selected_truth == index,
            selected_probability[:, index],
        )
    supported = [value for value in per_class_ap.values() if np.isfinite(value)]
    return {
        **classification,
        "expected_calibration_error": float(
            expected_calibration_error(
                selected_truth,
                selected_probability,
            )["ece"]
        ),
        "macro_average_precision": (
            float(np.mean(supported)) if supported else float("nan")
        ),
        "per_class_average_precision": per_class_ap,
    }


def _fit_linear_probe(
    *,
    train_embedding: np.ndarray,
    train_target: np.ndarray,
    validation_embedding: np.ndarray,
    seed: int,
    steps: int,
) -> np.ndarray:
    import torch
    from torch.nn import functional

    torch.manual_seed(int(seed))
    train_embedding = np.asarray(train_embedding, dtype=np.float32)
    validation_embedding = np.asarray(validation_embedding, dtype=np.float32)
    train_target = np.asarray(train_target, dtype=np.int64)
    center = train_embedding.mean(axis=0, keepdims=True)
    scale = train_embedding.std(axis=0, keepdims=True)
    scale[~np.isfinite(scale) | (scale < 1.0e-6)] = 1.0
    train_values = torch.from_numpy((train_embedding - center) / scale)
    validation_values = torch.from_numpy(
        (validation_embedding - center) / scale
    )
    targets = torch.from_numpy(train_target)
    counts = np.bincount(train_target, minlength=len(MORPHOLOGY_CLASSES))
    if (counts == 0).any():
        raise RuntimeError("a linear-probe training fold lacks a morphology class")
    class_weight = np.sqrt(counts.max() / counts)
    class_weight /= class_weight.mean()
    weights = torch.from_numpy(class_weight.astype(np.float32))
    probe = torch.nn.Linear(train_values.shape[1], len(MORPHOLOGY_CLASSES))
    optimizer = torch.optim.AdamW(
        probe.parameters(),
        lr=3.0e-2,
        weight_decay=1.0e-3,
    )
    probe.train()
    for _ in range(int(steps)):
        optimizer.zero_grad(set_to_none=True)
        loss = functional.cross_entropy(
            probe(train_values),
            targets,
            weight=weights,
        )
        loss.backward()
        optimizer.step()
    probe.eval()
    with torch.no_grad():
        return probe(validation_values).numpy().astype(np.float64)


def run_frozen_linear_probe(
    *,
    all_rows: pd.DataFrame,
    embeddings_by_fold: Sequence[np.ndarray],
    policy: str,
    output_dir: Path,
    seed: int,
    steps: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Run a fold-safe frozen-encoder morphology probe."""

    rows = _policy_rows(all_rows, policy)
    source_positions = all_rows.reset_index().set_index("review_id")["index"]
    positions = source_positions.loc[rows["review_id"].astype(str)].to_numpy(
        dtype=int
    )
    parts: list[pd.DataFrame] = []
    for fold in range(5):
        target = rows["morphology_target_index"].to_numpy(dtype=int)
        train = (rows["cv_fold"].to_numpy(dtype=int) != fold) & (target >= 0)
        validation = (rows["cv_fold"].to_numpy(dtype=int) == fold) & (target >= 0)
        if not train.any() or not validation.any():
            raise RuntimeError(f"linear probe fold {fold} has an empty partition")
        embedding = np.asarray(embeddings_by_fold[fold], dtype=np.float32)[
            positions
        ]
        logits = _fit_linear_probe(
            train_embedding=embedding[train],
            train_target=target[train],
            validation_embedding=embedding[validation],
            seed=int(seed) + 10_000 + fold,
            steps=int(steps),
        )
        selected = rows.loc[
            validation,
            ["review_id", "tic", "sector", "cv_fold", "human_label"],
        ].copy()
        selected["morphology_target"] = target[validation]
        for index, label in enumerate(MORPHOLOGY_CLASSES):
            selected[f"logit_{label}"] = logits[:, index]
        parts.append(selected)
    oof = pd.concat(parts, ignore_index=True)
    if oof["review_id"].duplicated().any():
        raise RuntimeError("linear-probe OOF predictions contain duplicates")
    expected = set(
        rows.loc[
            rows["morphology_target_index"].ge(0), "review_id"
        ].astype(str)
    )
    if set(oof["review_id"].astype(str)) != expected:
        raise RuntimeError("linear-probe OOF support is incomplete")
    logit_columns = [f"logit_{label}" for label in MORPHOLOGY_CLASSES]
    logits = oof.loc[:, logit_columns].to_numpy(dtype=float)
    truth = oof["morphology_target"].to_numpy(dtype=int)
    temperature = fit_temperature(logits, truth)
    probability = _softmax(logits, temperature=temperature)
    for index, label in enumerate(MORPHOLOGY_CLASSES):
        oof[f"p_{label}"] = probability[:, index]
    oof["morphology_prediction"] = probability.argmax(axis=1)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    predictions_path = output_dir / "development_oof_predictions.csv"
    oof.to_csv(
        predictions_path,
        index=False,
        lineterminator="\n",
        float_format="%.17g",
    )
    metrics = probability_metrics(truth, probability)
    summary = {
        "scope": "five_fold_real_development_oof_frozen_encoder_probe",
        "label_policy": policy,
        "encoder_selection": (
            "fold_k_encoder_excluded_complete_held_fold_k_tics"
        ),
        "probe": "fixed_step_class_weighted_linear_softmax",
        "steps": int(steps),
        "temperature": float(temperature),
        "n_rows": int(len(oof)),
        "n_tics": int(oof["tic"].nunique()),
        "review_ids_sha256": _string_set_sha256(
            oof["review_id"].astype(str).tolist()
        ),
        "metrics": metrics,
        "predictions": str(predictions_path),
        "predictions_sha256": _file_sha256(predictions_path),
    }
    _write_json(output_dir / "summary.json", summary)
    return oof, summary


def load_existing_linear_probe(
    *,
    all_rows: pd.DataFrame,
    policy: str,
    output_dir: Path,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Validate and load a completed probe from an interrupted evaluation."""

    rows = _policy_rows(all_rows, policy)
    output_dir = Path(output_dir).expanduser().resolve()
    summary_path = output_dir / "summary.json"
    summary = _read_json_object(
        summary_path,
        artifact=f"interrupted {policy} frozen-probe summary",
    )
    predictions_path = output_dir / "development_oof_predictions.csv"
    expected_summary = {
        "scope": "five_fold_real_development_oof_frozen_encoder_probe",
        "label_policy": policy,
        "encoder_selection": (
            "fold_k_encoder_excluded_complete_held_fold_k_tics"
        ),
        "probe": "fixed_step_class_weighted_linear_softmax",
    }
    changed = [
        name
        for name, value in expected_summary.items()
        if summary.get(name) != value
    ]
    if changed:
        raise RuntimeError(
            f"interrupted {policy} frozen probe changed: {changed}"
        )
    if Path(str(summary.get("predictions", ""))).resolve() != predictions_path:
        raise RuntimeError("interrupted frozen-probe prediction path changed")
    if _file_sha256(predictions_path) != str(
        summary.get("predictions_sha256", "")
    ):
        raise RuntimeError("interrupted frozen-probe prediction hash changed")
    oof = pd.read_csv(predictions_path)
    expected_rows = rows.loc[rows["morphology_target_index"].ge(0)].copy()
    expected_ids = set(expected_rows["review_id"].astype(str))
    if (
        oof["review_id"].astype(str).duplicated().any()
        or set(oof["review_id"].astype(str)) != expected_ids
    ):
        raise RuntimeError("interrupted frozen-probe support changed")
    aligned = expected_rows.set_index("review_id").loc[
        oof["review_id"].astype(str)
    ]
    for column, expected_values in (
        ("tic", pd.to_numeric(aligned["tic"], errors="raise").to_numpy()),
        ("cv_fold", aligned["cv_fold"].to_numpy(dtype=int)),
        (
            "morphology_target",
            aligned["morphology_target_index"].to_numpy(dtype=int),
        ),
    ):
        observed = pd.to_numeric(oof[column], errors="raise").to_numpy()
        if not np.array_equal(observed, expected_values):
            raise RuntimeError(
                f"interrupted frozen-probe {column} identity changed"
            )
    probability_columns = [f"p_{label}" for label in MORPHOLOGY_CLASSES]
    metrics = probability_metrics(
        oof["morphology_target"].to_numpy(dtype=int),
        oof.loc[:, probability_columns].to_numpy(dtype=float),
    )
    for name in (
        "accuracy",
        "balanced_accuracy",
        "macro_f1",
        "macro_average_precision",
        "expected_calibration_error",
    ):
        if not np.isclose(
            float(metrics[name]),
            float(summary["metrics"][name]),
            rtol=0.0,
            atol=1.0e-12,
        ):
            raise RuntimeError(
                f"interrupted frozen-probe metric changed: {name}"
            )
    loaded = dict(summary)
    loaded.update(
        {
            "resume_reused": True,
            "resume_summary": str(summary_path),
            "resume_summary_sha256": _file_sha256(summary_path),
        }
    )
    return oof, loaded


def load_existing_finetuned_policy(
    *,
    rows: pd.DataFrame,
    policy: str,
    fine_tune_root: Path,
    encoders: Sequence[Mapping[str, Any]],
    expected_input_provenance: Mapping[str, str],
    expected_evaluation_code_revision: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Validate pooled OOF artifacts and calibrated checkpoints for reuse."""

    import torch

    profile_dir = Path(fine_tune_root).expanduser().resolve() / FULLPOOL_SSL_PROFILE
    calibration_path = profile_dir / "pooled_oof_calibration.json"
    calibration = _read_json_object(
        calibration_path,
        artifact=f"interrupted {policy} pooled calibration",
    )
    expected_calibration = {
        "schema_version": FULLPOOL_EVALUATION_OOF_SCHEMA,
        "scope": "five_fold_real_development_oof_fine_tuned",
        "label_policy": policy,
        "n_rows": int(len(rows)),
        "n_real_rows": int(len(rows)),
        "n_tics": int(rows["tic"].nunique()),
        "evaluation_code_revision": str(expected_evaluation_code_revision),
    }
    changed = [
        name
        for name, value in expected_calibration.items()
        if calibration.get(name) != value
    ]
    for name, value in expected_input_provenance.items():
        if name == "evaluation_code_revision":
            continue
        if calibration.get(name) != value:
            changed.append(name)
    if changed:
        raise RuntimeError(
            f"interrupted {policy} pooled calibration changed: "
            f"{sorted(set(changed))}"
        )
    calibration_sha256 = _file_sha256(calibration_path)
    oof_path = profile_dir / "development_oof_predictions.csv"
    if Path(str(calibration.get("oof_predictions", ""))).resolve() != oof_path:
        raise RuntimeError(f"interrupted {policy} OOF path changed")
    if _file_sha256(oof_path) != str(
        calibration.get("oof_predictions_sha256", "")
    ):
        raise RuntimeError(f"interrupted {policy} OOF hash changed")
    oof = pd.read_csv(oof_path)
    expected_ids = set(rows["review_id"].astype(str))
    if (
        oof["review_id"].astype(str).duplicated().any()
        or set(oof["review_id"].astype(str)) != expected_ids
    ):
        raise RuntimeError(f"interrupted {policy} OOF support changed")
    aligned = rows.set_index("review_id").loc[oof["review_id"].astype(str)]
    for column, expected_values in (
        ("tic", pd.to_numeric(aligned["tic"], errors="raise").to_numpy()),
        ("cv_fold", aligned["cv_fold"].to_numpy(dtype=int)),
        (
            "morphology_target",
            aligned["morphology_target_index"].to_numpy(dtype=int),
        ),
    ):
        observed = pd.to_numeric(oof[column], errors="raise").to_numpy()
        if not np.array_equal(observed, expected_values):
            raise RuntimeError(f"interrupted {policy} OOF {column} changed")

    embedding_path = profile_dir / "development_oof_embeddings.npz"
    if not embedding_path.is_file():
        raise FileNotFoundError(embedding_path)
    with np.load(embedding_path) as payload:
        embedding_ids = payload["review_id"].astype(str)
        embedding_tics = payload["tic"].astype(np.int64)
        embedding = payload["embedding"].astype(np.float32)
    if (
        set(embedding_ids) != expected_ids
        or embedding.shape != (len(rows), 128)
        or not np.isfinite(embedding).all()
    ):
        raise RuntimeError(f"interrupted {policy} OOF embeddings changed")
    tic_by_review = dict(
        zip(
            rows["review_id"].astype(str),
            pd.to_numeric(rows["tic"], errors="raise").to_numpy(
                dtype=np.int64
            ),
        )
    )
    if not np.array_equal(
        embedding_tics,
        np.asarray([tic_by_review[value] for value in embedding_ids]),
    ):
        raise RuntimeError(f"interrupted {policy} embedding TICs changed")

    checkpoint_hashes: dict[str, str] = {}
    for encoder in encoders:
        fold = int(encoder["fold"])
        checkpoint_path = profile_dir / f"fold_{fold}" / "teacher.pt"
        checkpoint = torch.load(
            checkpoint_path,
            map_location="cpu",
            weights_only=False,
        )
        expected_checkpoint = {
            "schema_version": FULLPOOL_EVALUATION_CHECKPOINT_SCHEMA,
            "run_id": FULLPOOL_EVALUATION_RUN_ID,
            "model_facing_name": FULLPOOL_SSL_MODEL_FACING_NAME,
            "checkpoint_namespace": FULLPOOL_EVALUATION_NAMESPACE,
            "label_policy": policy,
            "profile": FULLPOOL_SSL_PROFILE,
            "fold": fold,
            "ssl_encoder_checkpoint_sha256": encoder["checkpoint_sha256"],
            "temperature_calibration_scope": "pooled_development_oof",
            "pooled_oof_calibration_sha256": calibration_sha256,
            "evaluation_code_revision": str(
                expected_evaluation_code_revision
            ),
        }
        mismatch = [
            name
            for name, value in expected_checkpoint.items()
            if checkpoint.get(name) != value
        ]
        if mismatch:
            raise RuntimeError(
                f"interrupted {policy} fold {fold} checkpoint changed: "
                f"{mismatch}"
            )
        checkpoint_hashes[str(fold)] = _file_sha256(checkpoint_path)

    loaded = dict(calibration)
    loaded.update(
        {
            "resume_reused": True,
            "pooled_oof_calibration": str(calibration_path),
            "pooled_oof_calibration_sha256": calibration_sha256,
            "oof_embeddings": str(embedding_path),
            "oof_embeddings_sha256": _file_sha256(embedding_path),
            "checkpoint_sha256": checkpoint_hashes,
        }
    )
    return oof, loaded


def load_teacher_v3_baseline_oof(
    *,
    baseline_root: Path,
    policy: str,
    expected_review_ids: Sequence[str],
    expected_truth: Mapping[str, int],
    expected_tic: Mapping[str, int],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Load Teacher-v3 OOF predictions on exactly the requested real support."""

    root = Path(baseline_root).expanduser().resolve()
    summary_path = root / "summary.json"
    if _file_sha256(summary_path) != TEACHER_V3_EXPECTED_SUMMARY_SHA256:
        raise RuntimeError("Teacher-v3 baseline summary hash changed")
    if policy == "uncertain_as_other":
        profile_root = root / FULLPOOL_SSL_PROFILE
    elif policy == "uncertain_masked":
        profile_root = root / "sensitivities" / "uncertain_masked" / FULLPOOL_SSL_PROFILE
    else:
        raise ValueError(f"unknown baseline policy {policy!r}")
    parts: list[pd.DataFrame] = []
    artifact_hashes: dict[str, str] = {}
    for fold in range(5):
        path = profile_root / f"fold_{fold}" / "validation_predictions.csv"
        part = pd.read_csv(path)
        part["cv_fold"] = fold
        parts.append(part)
        artifact_hashes[str(fold)] = _file_sha256(path)
    baseline = pd.concat(parts, ignore_index=True)
    if baseline["review_id"].astype(str).duplicated().any():
        raise RuntimeError("Teacher-v3 baseline OOF contains duplicates")
    expected = set(str(value) for value in expected_review_ids)
    observed = set(baseline["review_id"].astype(str))
    missing = sorted(expected - observed)
    if missing:
        raise RuntimeError(
            "Teacher-v3 baseline lacks matched development rows: "
            f"{missing[:10]}"
        )
    baseline = baseline.loc[
        baseline["review_id"].astype(str).isin(expected)
    ].copy()
    baseline = baseline.sort_values("review_id", kind="stable").reset_index(drop=True)
    if set(baseline["review_id"].astype(str)) != expected:
        raise RuntimeError("Teacher-v3 matched support changed during selection")
    tic_by_review = {
        str(review_id): int(tic) for review_id, tic in expected_tic.items()
    }
    if set(tic_by_review) != expected:
        raise RuntimeError(
            "Teacher-v3 TIC authority does not exactly match requested support"
        )
    expected_tics = np.asarray(
        [tic_by_review[str(value)] for value in baseline["review_id"]],
        dtype=np.int64,
    )
    if "tic" in baseline:
        observed_tics = pd.to_numeric(
            baseline["tic"], errors="raise"
        ).to_numpy(dtype=np.int64)
        if not np.array_equal(observed_tics, expected_tics):
            raise RuntimeError("Teacher-v3 baseline TIC identity changed")
    else:
        # Legacy Teacher-v3 fold predictions predate the explicit TIC column.
        # Bind TIC only from the already-authorized matched evaluation rows so
        # whole-host bootstrap resampling remains exact and auditable.
        baseline["tic"] = expected_tics
    expected_target = np.asarray(
        [expected_truth[str(value)] for value in baseline["review_id"]],
        dtype=int,
    )
    if not np.array_equal(
        baseline["morphology_target"].to_numpy(dtype=int),
        expected_target,
    ):
        raise RuntimeError("Teacher-v3 baseline morphology truth changed")
    probability_columns = [f"p_{label}" for label in MORPHOLOGY_CLASSES]
    if not set(probability_columns).issubset(baseline.columns):
        raise RuntimeError("Teacher-v3 baseline lacks calibrated probabilities")
    probability = baseline.loc[:, probability_columns].to_numpy(dtype=float)
    summary = {
        "policy": policy,
        "baseline_root": str(root),
        "baseline_summary": str(summary_path),
        "baseline_summary_sha256": TEACHER_V3_EXPECTED_SUMMARY_SHA256,
        "fold_prediction_sha256": artifact_hashes,
        "n_rows": int(len(baseline)),
        "n_tics": int(baseline["tic"].nunique()),
        "review_ids_sha256": _string_set_sha256(
            baseline["review_id"].astype(str).tolist()
        ),
        "metrics": probability_metrics(expected_target, probability),
    }
    return baseline, summary


def _tic_cluster_bootstrap(
    *,
    predictions: pd.DataFrame,
    draws: int,
    seed: int,
) -> dict[str, Any]:
    probability_columns = [f"p_{label}" for label in MORPHOLOGY_CLASSES]
    active = predictions["morphology_target"].to_numpy(dtype=int) >= 0
    frame = predictions.loc[active].copy().reset_index(drop=True)
    truth = frame["morphology_target"].to_numpy(dtype=int)
    probability = frame.loc[:, probability_columns].to_numpy(dtype=float)
    tics = pd.to_numeric(frame["tic"], errors="raise").to_numpy(dtype=np.int64)
    unique_tics = np.unique(tics)
    clusters = [np.flatnonzero(tics == tic) for tic in unique_tics]
    metric_names = (
        "balanced_accuracy",
        "macro_f1",
        "macro_average_precision",
        "expected_calibration_error",
    )
    values = {
        name: np.full(int(draws), np.nan, dtype=float) for name in metric_names
    }
    rng = np.random.default_rng(int(seed))
    for draw in range(int(draws)):
        sampled = rng.integers(0, len(clusters), size=len(clusters))
        indices = np.concatenate([clusters[int(index)] for index in sampled])
        metrics = probability_metrics(truth[indices], probability[indices])
        for name in metric_names:
            values[name][draw] = float(metrics[name])
    intervals = {
        name: {
            "low": float(np.nanquantile(metric_values, 0.025)),
            "high": float(np.nanquantile(metric_values, 0.975)),
        }
        for name, metric_values in values.items()
    }
    return {
        "unit": "tic",
        "draws": int(draws),
        "seed": int(seed),
        "confidence_level": 0.95,
        "n_tics": int(len(unique_tics)),
        "intervals": intervals,
    }


def _comparison_record(
    *,
    model: str,
    policy: str,
    predictions: pd.DataFrame,
    bootstrap_draws: int,
    bootstrap_seed: int,
) -> tuple[dict[str, Any], list[dict[str, Any]], dict[str, Any]]:
    probability_columns = [f"p_{label}" for label in MORPHOLOGY_CLASSES]
    truth = predictions["morphology_target"].to_numpy(dtype=int)
    probability = predictions.loc[:, probability_columns].to_numpy(dtype=float)
    metrics = probability_metrics(truth, probability)
    bootstrap = _tic_cluster_bootstrap(
        predictions=predictions,
        draws=int(bootstrap_draws),
        seed=int(bootstrap_seed),
    )
    record: dict[str, Any] = {
        "model": model,
        "label_policy": policy,
        "n_rows": int(metrics["n"]),
        "n_tics": int(predictions.loc[truth >= 0, "tic"].nunique()),
    }
    for name in (
        "accuracy",
        "balanced_accuracy",
        "macro_f1",
        "macro_average_precision",
        "expected_calibration_error",
    ):
        record[name] = float(metrics[name])
        if name in bootstrap["intervals"]:
            record[f"{name}_low"] = bootstrap["intervals"][name]["low"]
            record[f"{name}_high"] = bootstrap["intervals"][name]["high"]
    class_records: list[dict[str, Any]] = []
    for label in MORPHOLOGY_CLASSES:
        classification = metrics["per_class"][label]
        class_records.append(
            {
                "model": model,
                "label_policy": policy,
                "class": label,
                "support": int(classification["n"]),
                "precision": float(classification["precision"]),
                "recall": float(classification["recall"]),
                "f1": float(classification["f1"]),
                "average_precision": float(
                    metrics["per_class_average_precision"][label]
                ),
            }
        )
    return record, class_records, bootstrap


def run_fullpool_ssl_evaluation(
    *,
    training_table_path: Path,
    registry_path: Path,
    registry_summary_path: Path,
    completion_release_path: Path,
    baseline_teacher_v3_root: Path,
    output_dir: Path,
    fine_tune_epochs: int = 100,
    batch_size: int = 32,
    workers: int = 8,
    seed: int = 560064,
    linear_probe_steps: int = 500,
    bootstrap_draws: int = 2000,
    require_cuda: bool = True,
    resume: bool = False,
) -> dict[str, Any]:
    """Run representation export, probes, and full five-fold transfer tests."""

    for name, value in (
        ("fine_tune_epochs", fine_tune_epochs),
        ("batch_size", batch_size),
        ("linear_probe_steps", linear_probe_steps),
        ("bootstrap_draws", bootstrap_draws),
    ):
        if int(value) < 1:
            raise ValueError(f"{name} must be positive")
    if int(workers) < 0:
        raise ValueError("workers must be nonnegative")
    output = (
        Path(output_dir).expanduser().resolve()
        if resume
        else require_fresh_evaluation_output(output_dir)
    )
    code_revision = _code_revision()
    rows, label_audit = load_fullpool_development_labels(
        training_table_path=training_table_path,
        registry_path=registry_path,
        registry_summary_path=registry_summary_path,
        seed=int(seed),
    )
    encoders, encoder_audit = load_completed_fullpool_encoders(
        completion_release_path
    )
    input_provenance = {
        "evaluation_code_revision": code_revision,
        "training_table_sha256": label_audit["training_table_sha256"],
        "fullpool_registry_sha256": label_audit["registry_sha256"],
        "fullpool_registry_summary_sha256": label_audit[
            "registry_summary_sha256"
        ],
        "fullpool_completion_release_sha256": encoder_audit[
            "completion_release_sha256"
        ],
        "fullpool_training_code_revision": encoder_audit[
            "fullpool_training_code_revision"
        ],
    }
    run_contract = {
        "schema_version": FULLPOOL_EVALUATION_SCHEMA,
        "run_id": FULLPOOL_EVALUATION_RUN_ID,
        "model_facing_name": FULLPOOL_SSL_MODEL_FACING_NAME,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "real_development_labels_only",
        "label_policies": list(LABEL_POLICIES),
        "fold_encoder_rule": (
            "fold_k_transfer_and_probe_use_encoder_k_whose_ssl_pretraining_"
            "excluded_all_fold_k_tics"
        ),
        "fine_tune_epochs": int(fine_tune_epochs),
        "batch_size": int(batch_size),
        "workers": int(workers),
        "seed": int(seed),
        "linear_probe_steps": int(linear_probe_steps),
        "bootstrap_draws": int(bootstrap_draws),
        "fixed_test_container_opened": False,
        "fixed_test_tensors_constructed": False,
        "fixed_test_labels_used": False,
        "sector_63_rows_loaded": False,
        "injections_loaded": False,
        "automatic_production_promotion": False,
        **input_provenance,
    }
    contract_path = output / "run_contract.json"
    resume_audit: dict[str, Any] | None = None
    resume_contract_path: Path | None = None
    resume_contract_sha256: str | None = None
    if resume:
        resume_audit = _validate_resume_run_contract(
            output_dir=output,
            expected_contract=run_contract,
            current_input_provenance=input_provenance,
        )
        run_contract_sha256 = _file_sha256(contract_path)
        resume_contract_path = output / "resume_contract.json"
        if resume_contract_path.exists():
            raise RuntimeError(
                f"resume contract already exists: {resume_contract_path}"
            )
        resume_contract = {
            "schema_version": FULLPOOL_EVALUATION_SCHEMA,
            "run_id": FULLPOOL_EVALUATION_RUN_ID,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "action": "contract_checked_resume_after_evaluator_failure",
            "failure_repaired": (
                "legacy_teacher_v3_predictions_missing_tic_column; tic_"
                "bound_by_review_id_from_authorized_evaluation_support"
            ),
            **resume_audit,
            **input_provenance,
        }
        resume_contract_sha256 = _write_json(
            resume_contract_path,
            resume_contract,
        )
        embeddings, embedding_audit = (
            validate_and_recompute_existing_embeddings(
                rows=rows,
                encoders=encoders,
                output_dir=output / "representation",
                batch_size=int(batch_size),
                workers=int(workers),
                seed=int(seed),
                require_cuda=bool(require_cuda),
            )
        )
    else:
        run_contract_sha256 = _write_json(contract_path, run_contract)
        embeddings, embedding_audit = export_all_fold_embeddings(
            rows=rows,
            encoders=encoders,
            output_dir=output / "representation",
            batch_size=int(batch_size),
            workers=int(workers),
            seed=int(seed),
            require_cuda=bool(require_cuda),
        )

    comparison_records: list[dict[str, Any]] = []
    class_records: list[dict[str, Any]] = []
    policy_summaries: dict[str, Any] = {}
    for policy_index, policy in enumerate(LABEL_POLICIES):
        policy_dir = output / "policies" / policy
        policy_rows = _policy_rows(rows, policy)
        fine_tune_root = policy_dir / "fine_tuned"
        reuse_interrupted_policy = bool(resume and policy_index == 0)
        if reuse_interrupted_policy:
            if resume_audit is None:
                raise AssertionError("resume audit was not initialized")
            probe_oof, probe_summary = load_existing_linear_probe(
                all_rows=rows,
                policy=policy,
                output_dir=policy_dir / "frozen_linear_probe",
            )
            fine_oof, fine_summary = load_existing_finetuned_policy(
                rows=policy_rows,
                policy=policy,
                fine_tune_root=fine_tune_root,
                encoders=encoders,
                expected_input_provenance=input_provenance,
                expected_evaluation_code_revision=str(
                    resume_audit["interrupted_evaluation_code_revision"]
                ),
            )
            print(
                f"[fullpool-eval-resume] reused completed {policy} probe "
                "and five-fold fine-tune",
                flush=True,
            )
        else:
            policy_dir.mkdir(parents=True, exist_ok=False)
            probe_oof, probe_summary = run_frozen_linear_probe(
                all_rows=rows,
                embeddings_by_fold=embeddings,
                policy=policy,
                output_dir=policy_dir / "frozen_linear_probe",
                seed=int(seed),
                steps=int(linear_probe_steps),
            )
            fine_tune_root.mkdir(parents=True, exist_ok=False)
            fold_results: list[dict[str, Any]] = []
            artifact_context = {
                "checkpoint_schema": FULLPOOL_EVALUATION_CHECKPOINT_SCHEMA,
                "run_id": FULLPOOL_EVALUATION_RUN_ID,
                "model_facing_name": FULLPOOL_SSL_MODEL_FACING_NAME,
                "checkpoint_namespace": FULLPOOL_EVALUATION_NAMESPACE,
                "label_policy": policy,
                "log_name": f"fullpool-eval-{policy}",
            }
            for fold, encoder in enumerate(encoders):
                fold_results.append(
                    _run_finetune_fold(
                        rows=policy_rows,
                        fold=fold,
                        ssl_result=encoder,
                        out_dir=fine_tune_root,
                        input_provenance=input_provenance,
                        fine_tune_epochs=int(fine_tune_epochs),
                        batch_size=int(batch_size),
                        workers=int(workers),
                        seed=int(seed),
                        require_cuda=bool(require_cuda),
                        artifact_context=artifact_context,
                    )
                )
            fine_oof, fine_summary = _pool_development_oof(
                rows=policy_rows,
                fold_results=fold_results,
                out_dir=fine_tune_root,
                input_provenance=input_provenance,
                artifact_context={
                    "oof_schema": FULLPOOL_EVALUATION_OOF_SCHEMA,
                    "scope": "five_fold_real_development_oof_fine_tuned",
                    "label_policy": policy,
                },
            )
        active_fine = fine_oof.loc[
            fine_oof["morphology_target"].to_numpy(dtype=int) >= 0
        ].copy()
        expected_truth = dict(
            zip(
                active_fine["review_id"].astype(str),
                active_fine["morphology_target"].to_numpy(dtype=int),
            )
        )
        baseline, baseline_summary = load_teacher_v3_baseline_oof(
            baseline_root=baseline_teacher_v3_root,
            policy=policy,
            expected_review_ids=list(expected_truth),
            expected_truth=expected_truth,
            expected_tic=dict(
                zip(
                    active_fine["review_id"].astype(str),
                    pd.to_numeric(
                        active_fine["tic"], errors="raise"
                    ).to_numpy(dtype=np.int64),
                )
            ),
        )
        baseline_path = policy_dir / "teacher_v3_matched_oof_predictions.csv"
        if baseline_path.exists():
            raise RuntimeError(
                f"refusing to overwrite baseline predictions: {baseline_path}"
            )
        baseline.to_csv(
            baseline_path,
            index=False,
            lineterminator="\n",
            float_format="%.17g",
        )
        baseline_summary["matched_predictions"] = str(baseline_path)
        baseline_summary["matched_predictions_sha256"] = _file_sha256(
            baseline_path
        )

        model_predictions = (
            ("teacher_v3", baseline),
            ("fullpool_ssl_frozen_linear_probe", probe_oof),
            ("fullpool_ssl_fine_tuned", active_fine),
        )
        bootstraps: dict[str, Any] = {}
        for model_index, (model_name, predictions) in enumerate(
            model_predictions
        ):
            record, per_class, bootstrap = _comparison_record(
                model=model_name,
                policy=policy,
                predictions=predictions,
                bootstrap_draws=int(bootstrap_draws),
                bootstrap_seed=(
                    int(seed)
                    + 100_000 * policy_index
                    + 10_000 * model_index
                ),
            )
            comparison_records.append(record)
            class_records.extend(per_class)
            bootstraps[model_name] = bootstrap
        policy_summaries[policy] = {
            "n_supervised_rows": int(len(policy_rows)),
            "n_supervised_tics": int(policy_rows["tic"].nunique()),
            "n_morphology_rows": int(
                policy_rows["morphology_target_index"].ge(0).sum()
            ),
            "human_label_counts": {
                str(name): int(value)
                for name, value in policy_rows["human_label"]
                .astype(str)
                .value_counts()
                .sort_index()
                .items()
            },
            "frozen_linear_probe": probe_summary,
            "fine_tuned_oof": fine_summary,
            "teacher_v3_matched": baseline_summary,
            "bootstrap": bootstraps,
        }
        _write_json(policy_dir / "summary.json", policy_summaries[policy])

    comparison = pd.DataFrame(comparison_records)
    comparison_path = output / "development_performance.csv"
    comparison.to_csv(
        comparison_path,
        index=False,
        lineterminator="\n",
        float_format="%.17g",
    )
    per_class = pd.DataFrame(class_records)
    per_class_path = output / "development_per_class_metrics.csv"
    per_class.to_csv(
        per_class_path,
        index=False,
        lineterminator="\n",
        float_format="%.17g",
    )

    final_checkpoint_hashes = {
        str(record["fold"]): _file_sha256(record["checkpoint"])
        for record in encoders
    }
    if final_checkpoint_hashes != encoder_audit["checkpoint_sha256"]:
        raise RuntimeError("a frozen full-pool encoder changed during evaluation")
    summary = {
        "schema_version": FULLPOOL_EVALUATION_SCHEMA,
        "run_id": FULLPOOL_EVALUATION_RUN_ID,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "evaluation_code_revision": code_revision,
        "run_contract": str(contract_path),
        "run_contract_sha256": run_contract_sha256,
        "resume": resume_audit,
        "resume_contract": (
            str(resume_contract_path) if resume_contract_path else None
        ),
        "resume_contract_sha256": resume_contract_sha256,
        "label_authority": label_audit,
        "encoder_authority": encoder_audit,
        "embedding_export": embedding_audit,
        "policies": policy_summaries,
        "development_performance": str(comparison_path),
        "development_performance_sha256": _file_sha256(comparison_path),
        "development_per_class_metrics": str(per_class_path),
        "development_per_class_metrics_sha256": _file_sha256(per_class_path),
        "frozen_encoder_hashes_unchanged": True,
        "fixed_test_status": "sealed_not_loaded_not_evaluated",
        "sector_63_status": "reserved_not_loaded",
        "estimate_status": "matched_development_oof_not_untouched_test",
        "automatic_production_promotion": False,
        "student_training_blocked": True,
        **input_provenance,
    }
    summary_path = output / "summary.json"
    _write_json(summary_path, summary)
    return _json_safe(summary)


def run_fullpool_ssl_evaluation_preflight(
    *,
    training_table_path: Path,
    registry_path: Path,
    registry_summary_path: Path,
    completion_release_path: Path,
    batch_size: int = 32,
    workers: int = 0,
    seed: int = 560064,
    rows_per_fold: int = 8,
    require_cuda: bool = True,
) -> dict[str, Any]:
    """Validate authorities and one balanced bounded encoder forward pass."""

    if int(rows_per_fold) < 1:
        raise ValueError("rows_per_fold must be positive")
    rows, label_audit = load_fullpool_development_labels(
        training_table_path=training_table_path,
        registry_path=registry_path,
        registry_summary_path=registry_summary_path,
        seed=int(seed),
    )
    encoders, encoder_audit = load_completed_fullpool_encoders(
        completion_release_path
    )
    sample = pd.concat(
        [
            rows.loc[rows["cv_fold"].eq(fold)].head(int(rows_per_fold))
            for fold in range(5)
        ],
        ignore_index=True,
    )
    embedding = _extract_encoder_embeddings(
        rows=sample,
        encoder_state=encoders[0]["encoder_state_dict"],
        batch_size=int(batch_size),
        workers=int(workers),
        seed=int(seed),
        require_cuda=bool(require_cuda),
    )
    final_hashes = {
        str(record["fold"]): _file_sha256(record["checkpoint"])
        for record in encoders
    }
    if final_hashes != encoder_audit["checkpoint_sha256"]:
        raise RuntimeError("a frozen encoder changed during evaluation preflight")
    return {
        "passed": True,
        "scope": "bounded_authority_and_forward_pass_preflight",
        "n_authorized_development_rows": label_audit["n_rows"],
        "n_authorized_development_tics": label_audit["n_tics"],
        "n_smoke_rows": int(len(sample)),
        "smoke_embedding": _embedding_diagnostics(embedding),
        "fullpool_completion_release_sha256": encoder_audit[
            "completion_release_sha256"
        ],
        "frozen_encoder_hashes_unchanged": True,
        "fixed_test_status": "sealed_not_loaded_not_evaluated",
        "sector_63_status": "reserved_not_loaded",
    }


__all__ = [
    "EXPECTED_MATCHED_DEVELOPMENT_ROWS",
    "EXPECTED_MATCHED_DEVELOPMENT_TICS",
    "EXPECTED_REAL_DEVELOPMENT_ROWS",
    "EXPECTED_REAL_DEVELOPMENT_TICS",
    "FULLPOOL_EVALUATION_RUN_ID",
    "LABEL_POLICIES",
    "load_completed_fullpool_encoders",
    "load_fullpool_development_labels",
    "probability_metrics",
    "run_fullpool_ssl_evaluation",
    "run_fullpool_ssl_evaluation_preflight",
]
