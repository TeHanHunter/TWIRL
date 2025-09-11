#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
TIC -> Gaia DR3 crossmatch (fast, verbose, with failure log) and CSV export.

Outputs per TIC:
- dr3_source_id (as integer, no scientific notation)
- ra, dec (deg)
- parallax (mas)
- distance_pc (simple 1/parallax; empty for <=0 or NaN)
- phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag

Speed-ups:
- Exact-ID MAST query for TIC -> DR2 (no radius search)
- Batched DR2 -> DR3 crossmatch using gaiadr3.dr2_neighbourhood (best=1)
- Optional progress bars via tqdm if installed

Usage:
  python tic_to_gaia_dr3_full.py \
    --in tic_ids_unique.txt \
    --out vanderburg_TOIs.csv \
    --log failed_tic_to_gaia.log \
    --batch 2000 --retries 3
"""

import argparse, sys, time, logging, pathlib
from typing import List, Dict, Tuple, Optional

import pandas as pd
from astroquery.mast import Catalogs
from astroquery.utils.tap.core import TapPlus
from astropy.table import Table

# Optional progress bars
try:
    from tqdm import tqdm
except Exception:  # fallback if tqdm not installed
    def tqdm(iterable=None, **kwargs):
        return iterable if iterable is not None else range(0)


DEFAULT_GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap"
DEFAULT_BATCH = 2000
DEFAULT_RETRIES = 3
RETRY_SLEEP_S = 2.0


def read_tic_ids(path: str) -> List[int]:
    p = pathlib.Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    ids = []
    with p.open("r") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            try:
                ids.append(int(s))
            except ValueError:
                logging.warning(f"Skipping non-integer line: {s}")
    return ids


def mast_tic_row_for_id(tic_id: int, retries: int = DEFAULT_RETRIES) -> Optional[pd.DataFrame]:
    """Use MAST TIC by exact ID; faster than radius queries."""
    for attempt in range(1, retries + 1):
        try:
            tbl = Catalogs.query_criteria(catalog="TIC", ID=tic_id)
            if isinstance(tbl, Table):
                df = tbl.to_pandas()
            else:
                df = tbl
            if df.shape[0] == 0:
                return None
            df = df[df["ID"].astype(int) == tic_id].reset_index(drop=True)
            if df.shape[0] == 0:
                return None
            return df.iloc[[0]]
        except Exception as e:
            logging.warning(f"[MAST] TIC {tic_id} attempt {attempt}/{retries} failed: {e}")
            time.sleep(RETRY_SLEEP_S)
    return None


def collect_dr2_ids(tic_ids: List[int], retries: int = DEFAULT_RETRIES
                    ) -> Tuple[Dict[int, int], List[int], List[int]]:
    """
    For each TIC, get the GAIA (DR2) id from the TIC row.
    Returns:
      - tic_to_dr2: mapping TIC -> DR2 source_id
      - tic_fail: TICs that didn't resolve a DR2 id
      - tic_missing: TICs not found in TIC at all
    """
    tic_to_dr2: Dict[int, int] = {}
    tic_fail: List[int] = []
    tic_missing: List[int] = []

    for tic in tqdm(tic_ids, desc="Querying TIC→DR2", unit="tic"):
        row = mast_tic_row_for_id(tic, retries=retries)
        if row is None:
            tic_missing.append(tic)
            continue
        gaia_val = row.iloc[0].get("GAIA", None)  # DR2 id in TIC v8
        if pd.isna(gaia_val):
            tic_fail.append(tic)
            continue
        try:
            dr2 = int(gaia_val)
            tic_to_dr2[tic] = dr2
        except Exception:
            tic_fail.append(tic)

    logging.info(f"[MAST] Found DR2 for {len(tic_to_dr2)} / {len(tic_ids)} TICs "
                 f"(missing:{len(tic_missing)}, no-GAIA:{len(tic_fail)})")
    return tic_to_dr2, tic_fail, tic_missing


def chunked(lst: List[int], n: int) -> List[List[int]]:
    return [lst[i:i+n] for i in range(0, len(lst), n)]


def fetch_dr2_to_dr3_with_coords(dr2_ids: List[int], tap_url: str,
                                 retries: int = DEFAULT_RETRIES,
                                 batch_size: int = DEFAULT_BATCH
                                 ) -> pd.DataFrame:
    """
    Query gaiadr3.dr2_neighbourhood (best=1) JOIN gaiadr3.gaia_source to get
    dr2_source_id -> (dr3_source_id, ra, dec, parallax, phot mags).
    """
    cols = ["dr2_source_id", "dr3_source_id", "ra", "dec", "parallax",
            "phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag"]
    if not dr2_ids:
        return pd.DataFrame(columns=cols)

    tap = TapPlus(url=tap_url)
    out_frames = []

    for batch in tqdm(chunked(dr2_ids, batch_size), desc="Resolving DR2→DR3", unit="batch"):
        in_list = ",".join(str(x) for x in batch)
        adql = f"""
        SELECT n.dr2_source_id, n.dr3_source_id,
               s.ra, s.dec, s.parallax,
               s.phot_g_mean_mag, s.phot_bp_mean_mag, s.phot_rp_mean_mag
        FROM gaiadr3.dr2_neighbourhood AS n
        JOIN gaiadr3.gaia_source AS s
          ON s.source_id = n.dr3_source_id
        WHERE n.dr2_source_id IN ({in_list})
          AND n.best = 1
        """
        df = None
        for attempt in range(1, retries + 1):
            try:
                job = tap.launch_job(adql, dump_to_file=False)
                res = job.get_results()
                df = res.to_pandas()
                break
            except Exception as e:
                logging.warning(f"[GAIA] Batch {batch[0]}..{batch[-1]} attempt {attempt}/{retries} failed: {e}")
                time.sleep(RETRY_SLEEP_S)
        if df is None:
            fb_adql = f"""
            SELECT n.dr2_source_id, n.dr3_source_id, n.angular_distance,
                   s.ra, s.dec, s.parallax,
                   s.phot_g_mean_mag, s.phot_bp_mean_mag, s.phot_rp_mean_mag
            FROM gaiadr3.dr2_neighbourhood AS n
            JOIN gaiadr3.gaia_source AS s
              ON s.source_id = n.dr3_source_id
            WHERE n.dr2_source_id IN ({in_list})
            """
            try:
                job = tap.launch_job(fb_adql, dump_to_file=False)
                res = job.get_results().to_pandas()
                res = res.sort_values(["dr2_source_id", "angular_distance"]).groupby("dr2_source_id", as_index=False).first()
                df = res[["dr2_source_id", "dr3_source_id", "ra", "dec", "parallax",
                          "phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag"]].copy()
            except Exception as e:
                logging.error(f"[GAIA] Fallback failed for batch {batch[0]}..{batch[-1]}: {e}")
                df = pd.DataFrame(columns=cols)
        out_frames.append(df)

    out = pd.concat(out_frames, ignore_index=True) if out_frames else pd.DataFrame(columns=cols)
    for col in ("dr2_source_id", "dr3_source_id"):
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce").astype("Int64")
    return out


def main():
    ap = argparse.ArgumentParser(description="Crossmatch TIC IDs to Gaia DR3 source IDs and export CSV with coords/parallax/distance/photometry.")
    ap.add_argument("--in", dest="infile", required=True, help="Input text file of TIC IDs (one per line)")
    ap.add_argument("--out", dest="outfile", default="vanderburg_TOIs.csv", help="Output CSV path")
    ap.add_argument("--log", dest="logfile", default="failed_tic_to_gaia.log", help="Failure log path")
    ap.add_argument("--batch", dest="batch", type=int, default=DEFAULT_BATCH, help="Gaia TAP IN(...) batch size")
    ap.add_argument("--retries", dest="retries", type=int, default=DEFAULT_RETRIES, help="Network/query retries")
    ap.add_argument("--tap", dest="tap_url", default=DEFAULT_GAIA_TAP, help="Gaia DR3 TAP URL")
    args = ap.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)]
    )

    tic_ids = read_tic_ids(args.infile)
    logging.info(f"Loaded {len(tic_ids)} TIC IDs from {args.infile}")

    tic_to_dr2, tic_no_gaia, tic_missing = collect_dr2_ids(tic_ids, retries=args.retries)
    dr2_ids = sorted(set(tic_to_dr2.values()))
    logging.info(f"Unique DR2 IDs to resolve -> DR3: {len(dr2_ids)}")

    dr2_to_dr3 = fetch_dr2_to_dr3_with_coords(dr2_ids, tap_url=args.tap_url, retries=args.retries, batch_size=args.batch)

    dr2_map: Dict[int, Tuple[Optional[int], Optional[float], Optional[float], Optional[float], Optional[float], Optional[float], Optional[float]]] = {}
    for _, r in dr2_to_dr3.iterrows():
        dr2 = int(r["dr2_source_id"]) if pd.notna(r["dr2_source_id"]) else None
        dr3 = int(r["dr3_source_id"]) if pd.notna(r["dr3_source_id"]) else None
        ra = float(r["ra"]) if pd.notna(r["ra"]) else None
        dec = float(r["dec"]) if pd.notna(r["dec"]) else None
        par = float(r["parallax"]) if pd.notna(r["parallax"]) else None
        gmag  = float(r["phot_g_mean_mag"]) if pd.notna(r["phot_g_mean_mag"]) else None
        bpmag = float(r["phot_bp_mean_mag"]) if pd.notna(r["phot_bp_mean_mag"]) else None
        rpmag = float(r["phot_rp_mean_mag"]) if pd.notna(r["phot_rp_mean_mag"]) else None
        if dr2 is not None:
            dr2_map[dr2] = (dr3, ra, dec, par, gmag, bpmag, rpmag)

    rows = []
    unresolved_tics = []
    for tic, dr2 in tqdm(tic_to_dr2.items(), desc="Assembling output", unit="tic"):
        dr3, ra, dec, par, gmag, bpmag, rpmag = dr2_map.get(dr2, (None, None, None, None, None, None, None))
        if par is not None and par > 0:
            distance_pc = 1000.0 / par
        else:
            distance_pc = None
        if dr3 is None:
            unresolved_tics.append(tic)
        rows.append({
            "tic_id": tic,
            "dr3_source_id": dr3,
            "ra": ra,
            "dec": dec,
            "parallax_mas": par,
            "distance_pc": distance_pc,
            "phot_g_mean_mag": gmag,
            "phot_bp_mean_mag": bpmag,
            "phot_rp_mean_mag": rpmag
        })

    out_df = pd.DataFrame(rows).sort_values("tic_id").reset_index(drop=True)
    out_df["dr3_source_id"] = out_df["dr3_source_id"].astype("Int64")  # keep integer formatting
    out_df.to_csv(args.outfile, index=False)
    logging.info(f"Wrote {len(out_df)} rows to {args.outfile}")

    fail_lines = []
    if tic_missing:
        fail_lines.append("# TIC not found in MAST TIC")
        fail_lines.extend(str(t) for t in sorted(tic_missing))
    if tic_no_gaia:
        fail_lines.append("\n# TIC found but TIC row has no GAIA (DR2) id")
        fail_lines.extend(str(t) for t in sorted(tic_no_gaia))
    if unresolved_tics:
        fail_lines.append("\n# TIC resolved to DR2 but no DR3 match found")
        fail_lines.extend(str(t) for t in sorted(unresolved_tics))

    if fail_lines:
        with open(args.logfile, "w") as f:
            f.write("\n".join(fail_lines) + "\n")
        logging.info(f"Wrote failure log with {len(tic_missing) + len(tic_no_gaia) + len(unresolved_tics)} TICs -> {args.logfile}")
    else:
        logging.info("All TICs resolved through to DR3; no failures logged.")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logging.error("Interrupted by user.")
        sys.exit(130)