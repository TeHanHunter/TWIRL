# S56 TWIRL-FS v1 BLS Smoke Test

Run date: `2026-05-18`

Input FITS tree:
`/pdo/users/tehan/tglc-gpu-production/hlsp_s0056_twirl_fs_v1`

Output:
`/pdo/users/tehan/tglc-gpu-production/bls_s0056_twirl_fs_v1_smoke`

Local copy:
`reports/stage2_search/bls_s56_twirl_fs_v1_smoke/`

## Setup

- Host: `pdogpu6`
- Runner: `scripts/stage2_search/run_bls_smoke_csv.py`
- Output format: CSV/JSON only, intentionally bypassing the sector runner's
  parquet/pyarrow dependency for a quick smoke test.
- Targets: WD 1856 plus the six random faint TWIRL-FS preview targets.
- Apertures: `DET_FLUX_SML`, `DET_FLUX`, `DET_FLUX_LAG`
- Period grid: `50,000` periods
- Peaks retained: `3` per target/aperture
- Wall time: `193.7 s`

## Key Result

The TWIRL-FS FITS tree is readable by the TWIRL BLS code path. WD 1856 is
recovered at the expected period in the small and medium apertures:

| TIC | aperture | period d | duration min | depth | depth SNR | SDE |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 267574918 | `DET_FLUX_SML` | 1.407883 | 9.0 | 0.0899 | 0.64 | 62.16 |
| 267574918 | `DET_FLUX` | 1.407883 | 9.0 | 0.2872 | 2.05 | 57.29 |
| 267574918 | `DET_FLUX_LAG` | 8.404481 | 30.0 | 0.4834 | 2.51 | 57.94 |

The WD 1856 period is within the quick-grid resolution of the published
`1.40793903 d` period.

## Highest-SDE Smoke Peaks

| TIC | aperture | period d | duration min | depth | depth SNR | SDE |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 2052947887 | `DET_FLUX_LAG` | 11.521548 | 30.0 | 14.3412 | 43.00 | 115.41 |
| 2052947887 | `DET_FLUX` | 11.521548 | 30.0 | 8.5657 | 25.69 | 96.93 |
| 1201143622 | `DET_FLUX_LAG` | 7.274186 | 30.0 | 3.5287 | 21.14 | 83.10 |
| 1201143622 | `DET_FLUX` | 7.274186 | 30.0 | 2.6934 | 16.13 | 67.56 |
| 267574918 | `DET_FLUX_SML` | 1.407883 | 9.0 | 0.0899 | 0.64 | 62.16 |
| 267574918 | `DET_FLUX_LAG` | 8.404481 | 30.0 | 0.4834 | 2.51 | 57.94 |
| 267574918 | `DET_FLUX` | 1.407883 | 9.0 | 0.2872 | 2.05 | 57.29 |

The extreme depths in a few faint random targets are not interpreted here;
they are the expected next-stage vetting problem. This smoke test only checks
that TWIRL-FS products run through BLS and that the benchmark signal is still
recoverable.

## Note

The standard sector runner imports `pyarrow` to write parquet outputs. On the
current `pdogpu6` `twirl-gpu-venv`, `pyarrow` fails to import because of a
`libstdc++` symbol mismatch. Michelle's BLS code may not use this environment,
but our PDO sector runner needs that runtime issue fixed before a full local
parquet-producing sweep on `pdogpu6`.
