# pip install astroquery astropy numpy
from astroquery.mast import Catalogs
from astropy.table import vstack, Table
import numpy as np

# --- input: put your TICs (one per line) in this file, or set tic_ids = [...] directly
infile = "tic_ids_unique.txt"
with open(infile) as f:
    tic_ids = [int(x.strip()) for x in f if x.strip().isdigit()]

chunksize = 1000
tables = []
for i in range(0, len(tic_ids), chunksize):
    chunk = tic_ids[i:i+chunksize]
    t = Catalogs.query_criteria(catalog="TIC", ID=chunk)
    if len(t):
        tables.append(t["ID","ra","dec","Tmag","plx","e_plx"])

if not tables:
    raise RuntimeError("No TIC rows returned.")

tbl = vstack(tables, metadata_conflicts="silent")

# keep one row per TIC (prefer first if duplicates exist)
_, idx = np.unique(tbl["ID"], return_index=True)
tbl = tbl[np.sort(idx)]

# distance from parallax (mas): d_pc = 1000/parallax  (NaN if parallax<=0 or missing)
p = np.array(tbl["plx"], dtype=float)
pe = np.array(tbl["e_plx"], dtype=float)
d = np.where(p > 0, 1000.0/p, np.nan)
de = np.where(p > 0, 1000.0*pe/(p**2), np.nan)

tbl["distance_pc"] = d
tbl["e_distance_pc"] = de

# save
tbl.write("tic_ra_dec_tmag_parallax.fits", overwrite=True)
tbl.to_pandas().to_csv("tic_ra_dec_tmag_parallax.csv", index=False)

print(f"Rows: {len(tbl)}  -> tic_ra_dec_tmag_parallax.(fits|csv)")