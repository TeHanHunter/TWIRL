"""Heuristic vetting features for TWIRL Stage 2 BLS candidates.

The vetter applies physically motivated cuts to the per-target BLS results to
demote false-positive classes that BLS SDE alone cannot separate from real
WD transits:

  * short-period stellar variability and PCEBs whose minima are wider than any
    plausible planet-sized companion could produce at that period
  * background-EB contamination that grows monotonically with aperture size
    (Day 2 feature; requires per-aperture depths, not in consolidated.parquet)
  * sinusoidal flux modulation that lacks a flat out-of-transit baseline
    (Day 2/3 feature; requires the light curve, not just BLS results)

This module provides the feature computation. Cuts and re-ranking are applied
by the driver script (``scripts/stage5_validation/vet_s56.py``).

The duration envelope is one-sided (upper bound only): grazing transits can be
arbitrarily short, so we cannot reject on short durations. We can only reject
durations longer than the chord-crossing maximum for the largest plausible
super-Jupiter (R_comp <= 2 R_jup; above this we are in the brown-dwarf or
stellar regime, a different FP class).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

_GM_SUN_SI = 1.32712440e20  # m^3 / s^2
_R_SUN_M = 6.957e8
_R_JUP_M = 6.9911e7
_DAY_S = 86400.0

DEFAULT_M_WD_MSUN = 0.4  # low-mass tail of the WD population (~1% of catalog
                          # sits below this); gives the most generous duration
                          # envelope and so the most conservative rejection.
                          # The cut barely changes between M=0.4 and M=1.0
                          # because the chord sum is dominated by R_comp_max.
DEFAULT_R_WD_RSUN = 0.018  # paired with M=0.4 M_sun via the WD mass-radius
                            # relation; negligible effect on the envelope
                            # (R_WD << R_comp_max).
DEFAULT_R_COMP_MAX_RJUP = 2.0
DEFAULT_CADENCE_BUFFER_MIN = 3.34  # 200 s ingress+egress smear at TESS FFI cadence

# BLS p_min aliasing wall: with bls_default.yaml's p_min=0.0833 d (2 h), the
# per-TIC max-SDE peak piles up in [0.083, 0.10) d for ~33% of all S56 WDs.
# Above this band is also where the fluid-body Roche limit for a Jupiter
# around a 0.6 M_sun WD sits (~0.45 d for fluid, ~0.10 d for an iron core),
# so cutting here jointly excludes BLS aliasing and sub-Roche-limit orbits
# for any companion denser than an iron-rich rocky body. The right
# long-term fix is to set BLS p_min higher; the cut here is a vetter
# safety net.
DEFAULT_P_ALIAS_MIN_D = 0.10


def duration_envelope_min(
    period_d,
    m_wd_msun: float = DEFAULT_M_WD_MSUN,
    r_wd_rsun: float = DEFAULT_R_WD_RSUN,
    r_comp_max_rjup: float = DEFAULT_R_COMP_MAX_RJUP,
    cadence_buffer_min: float = DEFAULT_CADENCE_BUFFER_MIN,
):
    """Upper-bound observed transit duration (minutes) for a WD candidate at
    orbital period ``period_d`` (days).

    Edge-on central transit (b=0) of the largest plausible super-Jupiter
    companion, plus a constant cadence-smear additive buffer. Above this
    duration, the candidate's chord cannot be produced by a planet-class
    companion at that period and is therefore either a PCEB / BD or stellar
    variability.

        T_dur,max = 2 (R_WD + R_comp_max) / v_orb(P, M_WD) + cadence_buffer
    """
    p_s = np.asarray(period_d, dtype=float) * _DAY_S
    gm = m_wd_msun * _GM_SUN_SI
    v_orb = (2.0 * np.pi * gm / p_s) ** (1.0 / 3.0)  # m/s
    r_sum = r_wd_rsun * _R_SUN_M + r_comp_max_rjup * _R_JUP_M
    t_dur_min = (2.0 * r_sum / v_orb) / 60.0 + cadence_buffer_min
    return t_dur_min


def add_vetting_features(
    df: pd.DataFrame,
    grid_max_duration_min: float = 20.0,
    p_alias_min_d: float = DEFAULT_P_ALIAS_MIN_D,
    **envelope_kwargs,
) -> pd.DataFrame:
    """Append heuristic-vetter features to a BLS consolidated table.

    Adds columns:
      * ``dur_envelope_min`` — upper duration cap at the row's fitted period
      * ``dur_envelope_pass`` — ``duration_min <= dur_envelope_min``;
        auto-True when the envelope itself exceeds the BLS duration grid
        ceiling (the cut is unconstrained from this run)
      * ``dur_grid_ceiling_hit`` — duration_min is at the BLS grid maximum,
        so the true best-fit duration may be larger and the candidate is a
        PCEB / long-duration-EB suspect
      * ``p_alias_pass`` — period above the BLS p_min aliasing band

    Does not drop any rows. The driver decides which features become cuts.
    """
    out = df.copy()
    envelope = duration_envelope_min(out["period_d"].values, **envelope_kwargs)
    out["dur_envelope_min"] = envelope

    out["dur_grid_ceiling_hit"] = out["duration_min"].values >= grid_max_duration_min

    envelope_below_grid = envelope < grid_max_duration_min
    out["dur_envelope_pass"] = np.where(
        envelope_below_grid,
        out["duration_min"].values <= envelope,
        True,
    )

    out["p_alias_pass"] = out["period_d"].values >= p_alias_min_d

    # Period-pileup feature: count how many other TICs share a similar fitted
    # period. Astrophysical periods are essentially uniformly distributed in
    # log-P across unrelated WDs, so a sharp pileup (>50 TICs within a
    # ~0.1% period window) is almost always an instrumental alias wall.
    p_arr = out["period_d"].values
    log_p = np.log10(p_arr)
    bin_width_log10 = np.log10(1.001)  # ~0.1% relative band
    bin_idx = np.floor(log_p / bin_width_log10).astype(np.int64)
    counts = pd.Series(bin_idx).value_counts()
    out["period_cluster_count"] = bin_idx
    out["period_cluster_count"] = out["period_cluster_count"].map(counts).astype(int)
    return out
