from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
from astropy.table import Table


PDO_USER_ROOT = Path("/pdo/users/tehan")


def enforce_pdo_user_output(path: Path, *, allowed_root: Path = PDO_USER_ROOT) -> Path:
    """Return an absolute output path, rejecting writes outside the user PDO tree."""
    resolved = Path(path).expanduser().resolve(strict=False)
    allowed = Path(allowed_root).expanduser().resolve(strict=False)
    if resolved != allowed and allowed not in resolved.parents:
        raise ValueError(f"Refusing to write outside {allowed}: {resolved}")
    return resolved


def gaia_source_ids_from_designations(designations: Iterable[object]) -> np.ndarray:
    """Parse Gaia DR3 source IDs from strings like ``Gaia DR3 123``."""
    ids: list[int] = []
    for value in designations:
        text = str(value)
        token = text.rsplit(" ", 1)[-1]
        try:
            ids.append(int(token))
        except ValueError:
            continue
    return np.asarray(ids, dtype=np.int64)


def _find_column(table: Table, candidates: tuple[str, ...]) -> str:
    lower_to_name = {name.lower(): name for name in table.colnames}
    for candidate in candidates:
        found = lower_to_name.get(candidate.lower())
        if found is not None:
            return found
    raise KeyError(f"None of {candidates} found in table columns {table.colnames}")


def _clean_int_column(values: Iterable[object]) -> tuple[np.ndarray, np.ndarray]:
    out: list[int] = []
    good: list[bool] = []
    for value in values:
        if np.ma.is_masked(value) or value is None:
            out.append(-1)
            good.append(False)
            continue
        try:
            if isinstance(value, str) and value.startswith("Gaia DR3 "):
                value = value.rsplit(" ", 1)[-1]
            if isinstance(value, float) and not np.isfinite(value):
                raise ValueError
            out.append(int(value))
            good.append(True)
        except (TypeError, ValueError):
            out.append(-1)
            good.append(False)
    return np.asarray(out, dtype=np.int64), np.asarray(good, dtype=bool)


def normalize_tic_catalog(tic_catalog: Table, *, tmag_max: float | None = None) -> Table:
    """Normalize TGLC/TIC catalog variants to TIC, Gaia DR3 source ID, and Tmag."""
    tic_col = _find_column(tic_catalog, ("TIC", "ID", "tic", "id"))
    tmag_col = _find_column(tic_catalog, ("Tmag", "tmag", "TESSMAG", "tessmag"))
    try:
        gaia_col = _find_column(tic_catalog, ("gaia3", "dr3_source_id", "gaia_designation"))
    except KeyError:
        gaia_col = _find_column(tic_catalog, ("GAIA3", "GAIA"))

    tic, tic_good = _clean_int_column(tic_catalog[tic_col])
    gaia3, gaia_good = _clean_int_column(tic_catalog[gaia_col])
    tmag = np.asarray(tic_catalog[tmag_col], dtype=float)
    good = tic_good & gaia_good & np.isfinite(tmag)
    if tmag_max is not None:
        good &= tmag <= float(tmag_max)

    table = Table()
    table["TIC"] = tic[good]
    table["gaia3"] = gaia3[good]
    table["tmag"] = tmag[good]

    if len(table) == 0:
        return table

    order = np.lexsort((np.asarray(table["gaia3"]), np.asarray(table["TIC"])))
    table = table[order]
    keep = np.ones(len(table), dtype=bool)
    keep[1:] = np.asarray(table["TIC"][1:]) != np.asarray(table["TIC"][:-1])
    return table[keep]


def source_target_table(
    source_gaia: Table,
    tic_catalog: Table,
    *,
    tmag_max: float | None = None,
    max_targets: int | None = None,
) -> Table:
    """Build the two-column ``source.tic`` table expected by TGLC light-curve code."""
    designation_col = _find_column(source_gaia, ("designation", "DESIGNATION"))
    source_ids = gaia_source_ids_from_designations(source_gaia[designation_col])
    if source_ids.size == 0:
        return Table({"TIC": np.array([], dtype=np.int64), "gaia3": np.array([], dtype=np.int64)})

    normalized_cols = {"TIC", "gaia3", "tmag"}.issubset(set(tic_catalog.colnames))
    if normalized_cols:
        targets = tic_catalog
        if tmag_max is not None:
            targets = targets[np.asarray(targets["tmag"], dtype=float) <= float(tmag_max)]
    else:
        targets = normalize_tic_catalog(tic_catalog, tmag_max=tmag_max)
    if len(targets) == 0:
        return Table({"TIC": np.array([], dtype=np.int64), "gaia3": np.array([], dtype=np.int64)})

    in_source = np.isin(np.asarray(targets["gaia3"], dtype=np.int64), source_ids)
    matched = targets[in_source]
    if max_targets is not None:
        matched = matched[: int(max_targets)]

    out = Table()
    out["TIC"] = np.asarray(matched["TIC"], dtype=np.int64)
    out["gaia3"] = np.asarray(matched["gaia3"], dtype=np.int64)
    return out
