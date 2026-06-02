# S56 TWIRL-FS v2 Compare Handoff

Product root on PDO:

```text
/pdo/users/tehan/tglc-gpu-production/hlsp_s0056_twirl_fs_v2_compare
```

Each FITS file has two sets of detrending settings. 

## Which Flux Column To Use

Use `TIME`, `QUALITY`, and `ORBITID` exactly the same way for both versions.
For normal transit-search tests, start with `QUALITY == 0`.

| Column | Meaning | Recommendation |
|---|---|---|
| `DET_FLUX` | Canonical TWIRL-FS v2 detrending: robust subtractive spline with `BKSPACE=0.8 d` and `GAPSPLIT=0.5 d`. | Default / conservative search input. |
| `DET_FLUX_ERR` | Per-cadence uncertainty estimate for `DET_FLUX`. | Use with `DET_FLUX`. |
| `DET_FLUX_SML` | Canonical v2 detrended `1x1` aperture. | Aperture check. |
| `DET_FLUX_LAG` | Canonical v2 detrended `5x5` aperture. | Aperture check / contamination check. |
| `DET_FLUX_ADP` | Experimental adaptive/tighter detrending: same subtractive model, spline spacing `ADPBKSP=0.3 d`, quantile knot placement, and adaptive gap split `ADPGAP=0.2 d`. | Test when broad residual structure remains after v2. |
| `DET_FLUX_ADP_ERR` | Per-cadence uncertainty estimate for `DET_FLUX_ADP`. | Use with `DET_FLUX_ADP`. |
| `DET_FLUX_ADP_SML` | Experimental tighter detrended `1x1` aperture. | Aperture check for the adaptive version. |
| `DET_FLUX_ADP_LAG` | Experimental tighter detrended `5x5` aperture. | Aperture check for the adaptive version. |

The adaptive columns are intentionally opt-in. They remove more broad structure,
but they are more likely to attenuate long or shallow astrophysical features.
For short WD-like events, compare detection statistics on both columns rather
than replacing the canonical column.

## Header Keys

- `METHOD`: canonical `DET_FLUX` method, expected `twirl-fs-v2`.
- `BKSPACE`, `GAPSPLIT`, `NSEG`, `FITCNT`, `SCALESRC`, `COTSTAT`: diagnostics for canonical `DET_FLUX`.
- `HASADP`: whether adaptive compare columns are present.
- `ADPMETH`, `ADPBKSP`, `ADPGAP`, `ADPKNOT`, `ADPNSEG`, `ADPFIT`, `ADPSCAL`, `ADPCOTS`: diagnostics for `DET_FLUX_ADP`.

## Minimal Python Example

```python
from astropy.io import fits
import numpy as np

path = "hlsp_twirlfs_tess_ffi_s0056-0000000267574918_tess_v01_llc.fits"
with fits.open(path) as hdul:
    lc = hdul[1].data
    good = lc["QUALITY"] == 0
    time = lc["TIME"][good]
    flux_v2 = lc["DET_FLUX"][good]
    flux_adp = lc["DET_FLUX_ADP"][good]

# Run the same search on flux_v2 and flux_adp, then compare rankings.
```
