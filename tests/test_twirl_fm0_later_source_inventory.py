from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest

from twirl.models.fm0 import later_source_inventory as later
from twirl.models.fm0.corpus import (
    CORPUS_SELECTION_FIELDS,
    CORPUS_SELECTION_SCHEMA_VERSION,
)
from twirl.models.fm0.orcd_source_admission import verify_retained_sector_source
from twirl.models.fm0.registry import FM0ContractError, build_alias_registry

SECTOR = 66
ORBITS = (139, 140)
TARGET_AUTHORITY_SHA256 = "d" * 64
PRODUCER_GIT_SHA = "e" * 40


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True) + "\n", encoding="utf-8")


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def _write_corpus_rows(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=CORPUS_SELECTION_FIELDS, lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def _write_ecsv(path: Path, pairs: list[tuple[int, int]]) -> None:
    lines = [
        "# %ECSV 1.0",
        "# ---",
        "# datatype:",
        "# - {name: TIC, datatype: int64}",
        "# - {name: gaia3, datatype: int64}",
        "# schema: astropy-2.0",
        "TIC gaia3",
        *(f"{tic} {gaia}" for gaia, tic in pairs),
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _target_identity(camera: int, ccd: int) -> tuple[int, int]:
    suffix = (camera * 10) + ccd
    return 200_000 + suffix, 100_000 + suffix


def _make_sector(
    tmp_path: Path, *, extra_overlay_pair: tuple[int, int] | None = None
) -> tuple[Path, Path]:
    sector_root = tmp_path / f"s{SECTOR:04d}"
    for orbit in ORBITS:
        for camera in range(1, 5):
            for ccd in range(1, 5):
                gaia, tic = _target_identity(camera, ccd)
                cell = f"s{SECTOR:04d}_o{orbit}_cam{camera}_ccd{ccd}"
                input_root = (
                    sector_root
                    / "inputs"
                    / f"s{SECTOR:04d}-o{orbit}-cam{camera}-ccd{ccd}"
                )
                overlay = input_root / "source_tic" / "source_0_0.ecsv"
                pairs = [(gaia, tic)]
                if (
                    extra_overlay_pair is not None
                    and orbit == ORBITS[0]
                    and camera == 1
                    and ccd == 1
                ):
                    pairs.append(extra_overlay_pair)
                _write_ecsv(overlay, pairs)
                requested = input_root / "metadata" / "requested_tics.txt"
                expected = input_root / "metadata" / "expected_h5_tics.txt"
                requested.parent.mkdir(parents=True, exist_ok=True)
                requested.write_text(f"{tic}\n", encoding="utf-8")
                expected.write_text(f"{tic}\n", encoding="utf-8")
                records = [
                    {
                        "path": "source/source_0_0.pkl",
                        "bytes": len(b"cleaned-source"),
                        "sha256": hashlib.sha256(b"cleaned-source").hexdigest(),
                    },
                    {
                        "path": "source_tic/source_0_0.ecsv",
                        "bytes": overlay.stat().st_size,
                        "sha256": _sha(overlay),
                    },
                    {
                        "path": "metadata/requested_tics.txt",
                        "bytes": requested.stat().st_size,
                        "sha256": _sha(requested),
                    },
                    {
                        "path": "metadata/expected_h5_tics.txt",
                        "bytes": expected.stat().st_size,
                        "sha256": _sha(expected),
                    },
                ]
                cell_manifest = {
                    "schema_version": 1,
                    "cell": {
                        "sector": SECTOR,
                        "orbit": orbit,
                        "orbit_tag": "o1" if orbit == ORBITS[0] else "o2",
                        "camera": camera,
                        "ccd": ccd,
                    },
                    "target_authority_sha256": TARGET_AUTHORITY_SHA256,
                    "tglc_commit": "a" * 40,
                    "counts": {"source": 1, "source_tic": 1, "metadata": 2, "total": 4},
                    "total_bytes": sum(int(record["bytes"]) for record in records),
                    "files": records,
                }
                input_manifest = input_root / "cell_manifest.json"
                _write_json(input_manifest, cell_manifest)
                manifest_sha = _sha(input_manifest)
                (input_root / "READY").write_text(manifest_sha + "\n", encoding="utf-8")

                cell_root = sector_root / "outputs" / cell
                retained_root = cell_root / "retained" / "slurm-1-v2"
                (retained_root / "epsf").mkdir(parents=True)
                (retained_root / "LC").mkdir()
                (retained_root / "metadata").mkdir()
                epsf = retained_root / "epsf" / "epsf_0_0.npy"
                hdf5 = retained_root / "LC" / f"{tic}.h5"
                epsf.write_bytes(f"epsf-{cell}".encode("ascii"))
                hdf5.write_bytes(f"hdf5-{cell}".encode("ascii"))
                (retained_root / "cell_manifest.json").write_bytes(
                    input_manifest.read_bytes()
                )
                (retained_root / "READY").write_text(
                    manifest_sha + "\n", encoding="utf-8"
                )
                for source in (requested, expected):
                    (retained_root / "metadata" / source.name).write_bytes(
                        source.read_bytes()
                    )
                output_manifest = retained_root / "output_manifest.sha256"
                output_manifest.write_text(
                    f"{_sha(epsf)}  epsf/{epsf.name}\n{_sha(hdf5)}  LC/{hdf5.name}\n",
                    encoding="utf-8",
                )
                completion_path = cell_root / "complete" / "slurm-1.json"
                completion = {
                    "schema": "twirl-a2v1-orcd-cell-complete-v1",
                    "state": "complete",
                    "cell": cell,
                    "attempt_id": "slurm-1",
                    "input_manifest_sha256": manifest_sha,
                    "output_manifest_sha256": _sha(output_manifest),
                    "environment_manifest_sha256": "c" * 64,
                }
                _write_json(completion_path, completion)
                retention = {
                    "schema": "twirl-a2v1-orcd-retention-v1",
                    "ok": True,
                    "cell": cell,
                    "attempt_id": "slurm-1",
                    "pdo_return_deferred": True,
                    "pdo_sector_accepted": False,
                    "retained_root": str(retained_root),
                    "retained_root_validation": {
                        "ok": True,
                        "checks": {
                            "completion_marker_gate_reused": True,
                            "exact_epsf_count_nonzero": True,
                            "hdf5_count_nonzero_tic_authority": True,
                            "output_manifest_paths": True,
                        },
                        "completion_json_sha256": _sha(completion_path),
                        "input_manifest_sha256": manifest_sha,
                        "output_manifest_sha256": _sha(output_manifest),
                        "environment_manifest_sha256": "c" * 64,
                        "outputs": {
                            "hdf5": {
                                "validated_tics": 1,
                                "missing_expected": [],
                                "extra_unrequested": [],
                            }
                        },
                    },
                }
                _write_json(
                    cell_root / "retained" / "slurm-1-v2.retention.json",
                    retention,
                )
    receipt = verify_retained_sector_source(
        sector_root=sector_root,
        sector=SECTOR,
        expected_orbits=ORBITS,
    )
    receipt_path = tmp_path / "source.json"
    _write_json(receipt_path, receipt)
    return receipt_path, sector_root


def _write_observation_inventory(
    path: Path, *, edge_camera: int = 4, edge_ccd: int = 4
) -> str:
    with path.open("w", encoding="utf-8", newline="") as handle:
        fields = (
            "source_id",
            "tic_id",
            "sector",
            "orbit",
            "camera",
            "ccd",
            "edge_warn",
        )
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for orbit in ORBITS:
            for camera in range(1, 5):
                for ccd in range(1, 5):
                    gaia, tic = _target_identity(camera, ccd)
                    writer.writerow(
                        {
                            "source_id": gaia,
                            "tic_id": tic,
                            "sector": SECTOR,
                            "orbit": orbit,
                            "camera": camera,
                            "ccd": ccd,
                            "edge_warn": camera == edge_camera and ccd == edge_ccd,
                        }
                    )
    return _sha(path)


def _component(gaia: int, tic: int) -> tuple[str, str]:
    registry = build_alias_registry(
        [{"gaia_dr3_source_id": str(gaia), "tic_id": str(tic)}]
    )
    alias = registry.aliases[0]
    return str(alias["leakage_component_id"]), str(alias["source_partition"])


def _write_corpus_inventory(
    path: Path,
    *,
    selected_detectors: set[tuple[int, int]] | None = None,
) -> str:
    selected = selected_detectors or {
        (camera, ccd)
        for camera in range(1, 5)
        for ccd in range(1, 5)
        if (camera, ccd) != (4, 4)
    }
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=CORPUS_SELECTION_FIELDS, lineterminator="\n"
        )
        writer.writeheader()
        for camera, ccd in sorted(selected):
            gaia, tic = _target_identity(camera, ccd)
            component, partition = _component(gaia, tic)
            writer.writerow(
                {
                    "schema_version": CORPUS_SELECTION_SCHEMA_VERSION,
                    "gaia_dr3_source_id": gaia,
                    "tic_id": tic,
                    "sector": SECTOR,
                    "camera": camera,
                    "ccd": ccd,
                    "leakage_component_id": component,
                    "source_partition": partition,
                    "orbits_json": "[139,140]",
                    "hdf5_paths_json": (
                        f'["orbit-139/ffi/cam{camera}/ccd{ccd}/LC/{tic}.h5",'
                        f'"orbit-140/ffi/cam{camera}/ccd{ccd}/LC/{tic}.h5"]'
                    ),
                }
            )
    return _sha(path)


def _write_alias_authority(path: Path) -> str:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("gaia_dr3_source_id", "tic_id"),
            lineterminator="\n",
        )
        writer.writeheader()
        for camera in range(1, 5):
            for ccd in range(1, 5):
                gaia, tic = _target_identity(camera, ccd)
                writer.writerow({"gaia_dr3_source_id": gaia, "tic_id": tic})
    return _sha(path)


def _write_corpus_summary(
    path: Path,
    *,
    selection_path: Path,
    n_s66_rows: int,
    alias_path: Path,
) -> str:
    payload = {
        "schema_version": CORPUS_SELECTION_SCHEMA_VERSION,
        "sectors": list(range(66, 78)),
        "n_observations": n_s66_rows,
        "observations_by_sector": {
            str(sector): n_s66_rows if sector == SECTOR else 0
            for sector in range(66, 78)
        },
        "input_authorities": {
            "observation_fits": {
                "path": "/authority/twirl_wd_tess_observations_v0.fits",
                "sha256": TARGET_AUTHORITY_SHA256,
            },
            "gaia_tic_alias_table": {
                "path": str(alias_path),
                "sha256": _sha(alias_path),
            },
        },
        "outputs_sha256": {"selection.csv": _sha(selection_path)},
    }
    _write_json(path, payload)
    return _sha(path)


def _build_authorities(
    tmp_path: Path,
    *,
    selection_path: Path,
    n_selection_rows: int,
) -> tuple[Path, str, Path, str]:
    aliases = tmp_path / "gaia_tic_aliases.csv"
    aliases_sha = _write_alias_authority(aliases)
    summary = tmp_path / "summary.json"
    summary_sha = _write_corpus_summary(
        summary,
        selection_path=selection_path,
        n_s66_rows=n_selection_rows,
        alias_path=aliases,
    )
    return summary, summary_sha, aliases, aliases_sha


def _inventory_kwargs(
    tmp_path: Path,
    *,
    receipt: Path,
    selection: Path,
    n_selection_rows: int,
) -> dict[str, object]:
    summary, summary_sha, aliases, aliases_sha = _build_authorities(
        tmp_path,
        selection_path=selection,
        n_selection_rows=n_selection_rows,
    )
    return {
        "source_receipt_path": receipt,
        "source_receipt_sha256": _sha(receipt),
        "sector": SECTOR,
        "expected_orbits": ORBITS,
        "expected_target_authority_sha256": TARGET_AUTHORITY_SHA256,
        "target_inventory_path": selection,
        "target_inventory_sha256": _sha(selection),
        "corpus_summary_path": summary,
        "corpus_summary_sha256": summary_sha,
        "gaia_tic_alias_authority_path": aliases,
        "gaia_tic_alias_authority_sha256": aliases_sha,
        "output_dir": tmp_path / "inventory",
        "producer_git_sha": PRODUCER_GIT_SHA,
    }


def _build(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    *,
    extra_overlay_pair: tuple[int, int] | None = None,
) -> tuple[dict[str, object], Path, Path]:
    monkeypatch.setattr(later, "EXPECTED_SOURCE_GRID_FILES", 1)
    receipt, _ = _make_sector(tmp_path, extra_overlay_pair=extra_overlay_pair)
    targets = tmp_path / "selection.csv"
    target_sha = _write_corpus_inventory(targets)
    summary_path, summary_sha, aliases_path, aliases_sha = _build_authorities(
        tmp_path, selection_path=targets, n_selection_rows=15
    )
    output = tmp_path / "inventory"
    summary = later.build_later_source_inventory(
        source_receipt_path=receipt,
        source_receipt_sha256=_sha(receipt),
        sector=SECTOR,
        expected_orbits=ORBITS,
        expected_target_authority_sha256=TARGET_AUTHORITY_SHA256,
        target_inventory_path=targets,
        target_inventory_sha256=target_sha,
        corpus_summary_path=summary_path,
        corpus_summary_sha256=summary_sha,
        gaia_tic_alias_authority_path=aliases_path,
        gaia_tic_alias_authority_sha256=aliases_sha,
        output_dir=output,
        producer_git_sha=PRODUCER_GIT_SHA,
    )
    return summary, output, receipt


def test_checksum_bound_corpus_selection_groups_two_orbits_and_carries_partition(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    summary, output, _ = _build(tmp_path, monkeypatch)
    assert summary["n_hdf5_products_declared"] == 32
    assert summary["n_hdf5_products_selected_unique"] == 30
    assert summary["n_hdf5_products_excluded_not_in_target_inventory"] == 2
    assert summary["n_source_rows"] == 15
    with (output / "sources.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 15
    assert all(json.loads(row["orbits_json"]) == [139, 140] for row in rows)
    assert all(len(json.loads(row["hdf5_sha256_json"])) == 2 for row in rows)
    assert all(row["leakage_component_id"].startswith("leakage_") for row in rows)
    assert all(
        row["source_partition"] in summary["source_rows_by_partition"] for row in rows
    )
    assert all(row["model_outcome_blind"] == "true" for row in rows)
    assert all(row["panel_admission_authorized"] == "false" for row in rows)
    assert (output / "READY").read_text().strip() == _sha(output / "summary.json")
    assert all(
        not path.is_symlink() and not path.stat().st_mode & 0o222
        for path in (output, *output.rglob("*"))
    )


def test_publish_interruption_cleans_only_owned_partial_and_allows_retry(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(later, "EXPECTED_SOURCE_GRID_FILES", 1)
    receipt, _ = _make_sector(tmp_path)
    selection = tmp_path / "selection.csv"
    _write_corpus_inventory(selection, selected_detectors={(1, 1)})
    kwargs = _inventory_kwargs(
        tmp_path, receipt=receipt, selection=selection, n_selection_rows=1
    )
    output = Path(kwargs["output_dir"])
    sibling = tmp_path / "inventory-unrelated"
    sibling.write_text("keep\n", encoding="utf-8")
    original = later.publish_immutable

    def interrupt_summary(path: str | Path, payload: bytes) -> None:
        if Path(path).name == "summary.json":
            raise KeyboardInterrupt("simulated partial-tree interruption")
        original(path, payload)

    monkeypatch.setattr(later, "publish_immutable", interrupt_summary)
    with pytest.raises(KeyboardInterrupt, match="partial-tree interruption"):
        later.build_later_source_inventory(**kwargs)
    assert not output.exists()
    assert not output.with_name(output.name + ".partial").exists()
    assert sibling.read_text(encoding="utf-8") == "keep\n"

    monkeypatch.setattr(later, "publish_immutable", original)
    summary = later.build_later_source_inventory(**kwargs)
    assert summary["n_source_rows"] == 1
    assert output.is_dir()


def test_builder_refuses_broken_output_symlink_without_touching_target(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(later, "EXPECTED_SOURCE_GRID_FILES", 1)
    receipt, _ = _make_sector(tmp_path)
    selection = tmp_path / "selection.csv"
    _write_corpus_inventory(selection, selected_detectors={(1, 1)})
    kwargs = _inventory_kwargs(
        tmp_path, receipt=receipt, selection=selection, n_selection_rows=1
    )
    output = Path(kwargs["output_dir"])
    target = tmp_path / "outside-target"
    output.symlink_to(target, target_is_directory=True)

    with pytest.raises(FM0ContractError, match="refusing to overwrite"):
        later.build_later_source_inventory(**kwargs)
    assert output.is_symlink()
    assert not target.exists()


def test_manifest_bound_overlay_rejects_unrequested_all_tic_row(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(later, "EXPECTED_SOURCE_GRID_FILES", 1)
    receipt, _ = _make_sector(tmp_path, extra_overlay_pair=(999_999, 888_888))
    targets = tmp_path / "selection.csv"
    _write_corpus_inventory(targets)
    kwargs = _inventory_kwargs(
        tmp_path, receipt=receipt, selection=targets, n_selection_rows=15
    )
    with pytest.raises(FM0ContractError, match="outside the manifest-bound requested"):
        later.build_later_source_inventory(**kwargs)


def test_exact_corpus_selection_is_an_allowlist(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(later, "EXPECTED_SOURCE_GRID_FILES", 1)
    receipt, _ = _make_sector(tmp_path)
    gaia, tic = _target_identity(1, 1)
    component, partition = _component(gaia, tic)
    selection = tmp_path / "selection.csv"
    with selection.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=CORPUS_SELECTION_FIELDS, lineterminator="\n"
        )
        writer.writeheader()
        writer.writerow(
            {
                "schema_version": CORPUS_SELECTION_SCHEMA_VERSION,
                "gaia_dr3_source_id": gaia,
                "tic_id": tic,
                "sector": SECTOR,
                "camera": 1,
                "ccd": 1,
                "leakage_component_id": component,
                "source_partition": partition,
                "orbits_json": "[139,140]",
                "hdf5_paths_json": (
                    f'["orbit-139/ffi/cam1/ccd1/LC/{tic}.h5",'
                    f'"orbit-140/ffi/cam1/ccd1/LC/{tic}.h5"]'
                ),
            }
        )
    kwargs = _inventory_kwargs(
        tmp_path, receipt=receipt, selection=selection, n_selection_rows=1
    )
    output = Path(kwargs["output_dir"])
    summary = later.build_later_source_inventory(**kwargs)
    assert (
        summary["selection_identity_source"] == later.IDENTITY_SOURCE_CORPUS_SELECTION
    )
    assert summary["n_source_rows"] == 1
    assert summary["n_hdf5_products_selected_unique"] == 2
    assert summary["n_hdf5_products_excluded_not_in_target_inventory"] == 30
    with (output / "sources.csv").open(newline="") as handle:
        row = next(csv.DictReader(handle))
    assert row["leakage_component_id"] == component
    assert row["source_partition"] == partition


def test_output_manifest_must_remain_completion_bound(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(later, "EXPECTED_SOURCE_GRID_FILES", 1)
    receipt, sector_root = _make_sector(tmp_path)
    output_manifest = next(
        sector_root.glob("outputs/*/retained/*/output_manifest.sha256")
    )
    output_manifest.write_text(
        output_manifest.read_text(encoding="utf-8") + "0" * 64 + "  LC/999.h5\n",
        encoding="utf-8",
    )
    targets = tmp_path / "selection.csv"
    _write_corpus_inventory(targets)
    kwargs = _inventory_kwargs(
        tmp_path, receipt=receipt, selection=targets, n_selection_rows=15
    )
    with pytest.raises(
        FM0ContractError, match="output manifest is not completion-bound"
    ):
        later.build_later_source_inventory(**kwargs)


def test_observation_inventory_is_metadata_only_for_partitioned_sources(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(later, "EXPECTED_SOURCE_GRID_FILES", 1)
    receipt, _ = _make_sector(tmp_path)
    targets = tmp_path / "selection.csv"
    _write_observation_inventory(targets)
    kwargs = _inventory_kwargs(
        tmp_path, receipt=receipt, selection=targets, n_selection_rows=15
    )
    with pytest.raises(
        FM0ContractError, match="observation inventory is metadata-only"
    ):
        later.build_later_source_inventory(**kwargs)


def test_sealed_partition_hdf5_content_is_never_opened(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(later, "EXPECTED_SOURCE_GRID_FILES", 1)
    receipt, _ = _make_sector(tmp_path)
    selection = tmp_path / "selection.csv"
    _write_corpus_inventory(selection, selected_detectors={(2, 4)})
    gaia, tic = _target_identity(2, 4)
    assert _component(gaia, tic)[1] == "poc_sealed_test"

    original_open = Path.open

    def guarded_open(path: Path, *args: object, **kwargs: object):
        if path.suffix == ".h5":
            raise AssertionError(f"sealed HDF5 content was opened: {path}")
        return original_open(path, *args, **kwargs)

    monkeypatch.setattr(Path, "open", guarded_open)
    kwargs = _inventory_kwargs(
        tmp_path, receipt=receipt, selection=selection, n_selection_rows=1
    )
    summary = later.build_later_source_inventory(**kwargs)
    assert summary["source_rows_by_partition"]["poc_sealed_test"] == 1
    assert summary["hdf5_content_opened"] is False
    assert summary["sealed_hdf5_content_opened"] is False


def test_corpus_selection_hdf5_paths_must_match_exact_identity(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(later, "EXPECTED_SOURCE_GRID_FILES", 1)
    receipt, _ = _make_sector(tmp_path)
    selection = tmp_path / "selection.csv"
    _write_corpus_inventory(selection, selected_detectors={(1, 1)})
    rows = _read_csv(selection)
    rows[0]["hdf5_paths_json"] = (
        '["orbit-139/ffi/cam1/ccd1/LC/999999.h5",'
        '"orbit-140/ffi/cam1/ccd1/LC/999999.h5"]'
    )
    _write_corpus_rows(selection, rows)
    kwargs = _inventory_kwargs(
        tmp_path, receipt=receipt, selection=selection, n_selection_rows=1
    )

    with pytest.raises(FM0ContractError, match="hdf5_paths_json differs"):
        later.build_later_source_inventory(**kwargs)


def test_corpus_selection_rejects_one_gaia_mapped_to_two_selected_tics(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(later, "EXPECTED_SOURCE_GRID_FILES", 1)
    receipt, _ = _make_sector(tmp_path)
    selection = tmp_path / "selection.csv"
    _write_corpus_inventory(selection, selected_detectors={(1, 1), (1, 2)})
    rows = _read_csv(selection)
    rows[1]["gaia_dr3_source_id"] = rows[0]["gaia_dr3_source_id"]
    _write_corpus_rows(selection, rows)
    kwargs = _inventory_kwargs(
        tmp_path, receipt=receipt, selection=selection, n_selection_rows=2
    )

    with pytest.raises(FM0ContractError, match="ambiguous selected Gaia--TIC"):
        later.build_later_source_inventory(**kwargs)


def test_selected_tic_rejects_ambiguous_manifest_bound_overlay(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(later, "EXPECTED_SOURCE_GRID_FILES", 1)
    _gaia, tic = _target_identity(1, 1)
    receipt, _ = _make_sector(tmp_path, extra_overlay_pair=(999_999, tic))
    selection = tmp_path / "selection.csv"
    _write_corpus_inventory(selection, selected_detectors={(1, 1)})
    kwargs = _inventory_kwargs(
        tmp_path, receipt=receipt, selection=selection, n_selection_rows=1
    )

    with pytest.raises(FM0ContractError, match="ambiguous Gaia identity"):
        later.build_later_source_inventory(**kwargs)


def test_corpus_summary_must_bind_observation_authority(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(later, "EXPECTED_SOURCE_GRID_FILES", 1)
    receipt, _ = _make_sector(tmp_path)
    selection = tmp_path / "selection.csv"
    _write_corpus_inventory(selection, selected_detectors={(1, 1)})
    kwargs = _inventory_kwargs(
        tmp_path, receipt=receipt, selection=selection, n_selection_rows=1
    )
    summary_path = Path(kwargs["corpus_summary_path"])
    corpus_summary = json.loads(summary_path.read_text(encoding="utf-8"))
    corpus_summary["input_authorities"]["observation_fits"]["sha256"] = "f" * 64
    _write_json(summary_path, corpus_summary)
    kwargs["corpus_summary_sha256"] = _sha(summary_path)

    with pytest.raises(FM0ContractError, match="observation authority differs"):
        later.build_later_source_inventory(**kwargs)


def test_supplied_alias_authority_is_semantically_verified(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(later, "EXPECTED_SOURCE_GRID_FILES", 1)
    receipt, _ = _make_sector(tmp_path)
    selection = tmp_path / "selection.csv"
    _write_corpus_inventory(selection, selected_detectors={(1, 1)})
    kwargs = _inventory_kwargs(
        tmp_path, receipt=receipt, selection=selection, n_selection_rows=1
    )
    alias_path = Path(kwargs["gaia_tic_alias_authority_path"])
    alias_rows = _read_csv(alias_path)
    _gaia, selected_tic = _target_identity(1, 1)
    alias_rows.append({"gaia_dr3_source_id": "999999", "tic_id": str(selected_tic)})
    with alias_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("gaia_dr3_source_id", "tic_id"),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(alias_rows)
    alias_sha = _sha(alias_path)
    summary_path = Path(kwargs["corpus_summary_path"])
    corpus_summary = json.loads(summary_path.read_text(encoding="utf-8"))
    corpus_summary["input_authorities"]["gaia_tic_alias_table"]["sha256"] = alias_sha
    _write_json(summary_path, corpus_summary)
    kwargs["corpus_summary_sha256"] = _sha(summary_path)
    kwargs["gaia_tic_alias_authority_sha256"] = alias_sha

    with pytest.raises(FM0ContractError, match="differs from the supplied"):
        later.build_later_source_inventory(**kwargs)


def test_source_receipt_requires_expected_sha256(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(later, "EXPECTED_SOURCE_GRID_FILES", 1)
    receipt, _ = _make_sector(tmp_path)
    selection = tmp_path / "selection.csv"
    _write_corpus_inventory(selection, selected_detectors={(1, 1)})
    kwargs = _inventory_kwargs(
        tmp_path, receipt=receipt, selection=selection, n_selection_rows=1
    )
    kwargs["source_receipt_sha256"] = "0" * 64

    with pytest.raises(FM0ContractError, match="source-receipt SHA-256 drifted"):
        later.build_later_source_inventory(**kwargs)


@pytest.mark.parametrize("sector", [65, 78])
def test_later_source_boundary_excludes_sectors_outside_s66_s77(
    tmp_path: Path, sector: int
) -> None:
    with pytest.raises(FM0ContractError, match="restricted to S66--S77"):
        later.build_later_source_inventory(
            source_receipt_path=tmp_path / "unused.json",
            source_receipt_sha256="0" * 64,
            sector=sector,
            expected_orbits=(2 * sector + 7, 2 * sector + 8),
            expected_target_authority_sha256=TARGET_AUTHORITY_SHA256,
            target_inventory_path=tmp_path / "selection.csv",
            target_inventory_sha256="1" * 64,
            corpus_summary_path=tmp_path / "summary.json",
            corpus_summary_sha256="2" * 64,
            output_dir=tmp_path / "inventory",
            producer_git_sha=PRODUCER_GIT_SHA,
        )


def test_strict_loader_filters_sealed_rows_without_opening_hdf5(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    summary, output, _ = _build(tmp_path, monkeypatch)
    assert summary["source_rows_by_partition"]["poc_sealed_test"] > 0
    original_open = Path.open

    def guarded_open(path: Path, *args: object, **kwargs: object):
        if path.suffix == ".h5":
            raise AssertionError(f"sealed HDF5 content was opened: {path}")
        return original_open(path, *args, **kwargs)

    monkeypatch.setattr(Path, "open", guarded_open)
    rows = later.load_later_source_rows(
        output, expected_summary_sha256=_sha(output / "summary.json")
    )

    assert rows
    assert all(row["source_partition"] != "poc_sealed_test" for row in rows)
    assert (
        len(rows)
        == summary["n_source_rows"]
        - summary["source_rows_by_partition"]["poc_sealed_test"]
    )
    with pytest.raises(FM0ContractError, match="sealed rows cannot be loaded"):
        later.load_later_source_rows(
            output,
            expected_summary_sha256=_sha(output / "summary.json"),
            allowed_source_partitions=("poc_sealed_test",),
        )
