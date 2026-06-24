"""Light-curve injection primitives for TWIRL recovery tests."""
from __future__ import annotations

import numpy as np

TESS_200S_EXPOSURE_D = 200.0 / 86400.0


def box_transit_mask(
    time_d,
    *,
    period_d: float,
    t0_d: float,
    duration_min: float,
) -> np.ndarray:
    """Return cadences inside a periodic box transit.

    ``time_d`` and ``t0_d`` must use the same day zero point. For TWIRL HLSP
    FITS, ``TIME`` is ``BJD - 2457000``; convert full-BJD ephemerides before
    calling this helper.
    """
    if period_d <= 0:
        raise ValueError("period_d must be positive")
    if duration_min <= 0:
        raise ValueError("duration_min must be positive")
    time = np.asarray(time_d, dtype=np.float64)
    half_duration_d = 0.5 * duration_min / (24.0 * 60.0)
    phase_d = ((time - float(t0_d) + 0.5 * period_d) % period_d) - 0.5 * period_d
    return np.isfinite(time) & (np.abs(phase_d) <= half_duration_d)


def inject_box_transit(
    time_d,
    flux,
    *,
    period_d: float,
    t0_d: float,
    duration_min: float,
    depth: float,
    baseline: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Inject a box-shaped transit into a normalized or near-normalized flux array.

    The injected signal subtracts ``depth * baseline`` from in-transit cadences.
    When ``baseline`` is not supplied, the finite-flux median is used. This is
    appropriate for the first TWIRL-FS training pilot where detrended fluxes are
    close to unity, but the full occurrence-rate completeness run should still
    validate this approximation against the final adopted light-curve product.
    """
    if depth < 0:
        raise ValueError("depth must be non-negative")
    flux_arr = np.asarray(flux, dtype=np.float64)
    if baseline is None:
        finite = np.isfinite(flux_arr)
        if not np.any(finite):
            raise ValueError("cannot infer baseline from all-NaN flux")
        baseline = float(np.nanmedian(flux_arr[finite]))
    if not np.isfinite(baseline):
        raise ValueError("baseline must be finite")

    mask = box_transit_mask(
        time_d,
        period_d=period_d,
        t0_d=t0_d,
        duration_min=duration_min,
    )
    injected = flux_arr.copy()
    injected[mask] -= float(depth) * float(baseline)
    return injected, mask


def _require_batman():
    try:
        import batman  # type: ignore
    except ImportError as exc:  # pragma: no cover - exercised on PDO if dependency is absent
        raise ImportError(
            "batman-package is required for TWIRL transit injections; install/import "
            "`batman` in the active environment before generating injection products"
        ) from exc
    return batman


def batman_transit_model(
    time_d,
    *,
    period_d: float,
    t0_d: float,
    radius_rstar: float,
    a_over_rstar: float,
    impact_b: float,
    limb_darkening: tuple[float, float] = (0.05, 0.05),
    supersample_factor: int = 7,
    exp_time_d: float = TESS_200S_EXPOSURE_D,
) -> np.ndarray:
    """Return a finite-exposure `batman` transit model at the input cadences."""
    batman = _require_batman()
    time = np.asarray(time_d, dtype=np.float64)
    if period_d <= 0:
        raise ValueError("period_d must be positive")
    if radius_rstar <= 0:
        raise ValueError("radius_rstar must be positive")
    if a_over_rstar <= 0:
        raise ValueError("a_over_rstar must be positive")
    if impact_b < 0:
        raise ValueError("impact_b must be non-negative")
    if impact_b >= a_over_rstar:
        raise ValueError("impact_b must be smaller than a_over_rstar")
    if supersample_factor < 1:
        raise ValueError("supersample_factor must be positive")

    params = batman.TransitParams()
    params.t0 = float(t0_d)
    params.per = float(period_d)
    params.rp = float(radius_rstar)
    params.a = float(a_over_rstar)
    params.inc = float(np.degrees(np.arccos(np.clip(float(impact_b) / float(a_over_rstar), 0.0, 1.0))))
    params.ecc = 0.0
    params.w = 90.0
    params.u = [float(limb_darkening[0]), float(limb_darkening[1])]
    params.limb_dark = "quadratic"

    kwargs = {}
    if supersample_factor > 1:
        kwargs["supersample_factor"] = int(supersample_factor)
        kwargs["exp_time"] = float(exp_time_d)
    model = batman.TransitModel(params, time, **kwargs)
    return np.asarray(model.light_curve(params), dtype=np.float64)


def inject_batman_transit(
    time_d,
    flux,
    *,
    period_d: float,
    t0_d: float,
    duration_min: float,
    radius_rstar: float,
    a_over_rstar: float,
    impact_b: float,
    baseline: float | None = None,
    limb_darkening: tuple[float, float] = (0.05, 0.05),
    supersample_factor: int = 7,
    exp_time_d: float = TESS_200S_EXPOSURE_D,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Inject a limb-darkened batman transit into a flux array.

    The model is applied multiplicatively around the local continuum baseline:
    ``injected = flux + baseline * (model - 1)``. This preserves the observed
    residual structure while injecting a physical transit shape.
    """
    flux_arr = np.asarray(flux, dtype=np.float64)
    if baseline is None:
        finite = np.isfinite(flux_arr)
        if not np.any(finite):
            raise ValueError("cannot infer baseline from all-NaN flux")
        baseline = float(np.nanmedian(flux_arr[finite]))
    if not np.isfinite(baseline):
        raise ValueError("baseline must be finite")

    time = np.asarray(time_d, dtype=np.float64)
    model = np.ones_like(time, dtype=np.float64)
    model_window = box_transit_mask(
        time,
        period_d=period_d,
        t0_d=t0_d,
        duration_min=max(
            3.0 * float(duration_min),
            float(duration_min) + 4.0 * float(exp_time_d) * 1440.0,
        ),
    )
    if np.any(model_window):
        model[model_window] = batman_transit_model(
            time[model_window],
            period_d=period_d,
            t0_d=t0_d,
            radius_rstar=radius_rstar,
            a_over_rstar=a_over_rstar,
            impact_b=impact_b,
            limb_darkening=limb_darkening,
            supersample_factor=supersample_factor,
            exp_time_d=exp_time_d,
        )
    mask = box_transit_mask(
        time,
        period_d=period_d,
        t0_d=t0_d,
        duration_min=duration_min,
    )
    injected = flux_arr + float(baseline) * (model - 1.0)
    return injected, mask, model


def choose_observed_epoch(
    time_d,
    *,
    period_d: float,
    duration_min: float,
    rng: np.random.Generator,
    quality=None,
    finite_mask=None,
    min_in_transit: int = 1,
    max_tries: int = 1000,
) -> tuple[float, np.ndarray]:
    """Choose ``t0`` so the injected transit window overlaps observed cadences.

    This avoids creating nominal injections that fall entirely in gaps or only
    on flagged/non-finite cadences. The returned mask is the full in-transit
    cadence mask; callers can combine it with their own good-cadence mask for
    recovery statistics.
    """
    if min_in_transit <= 0:
        raise ValueError("min_in_transit must be positive")
    time = np.asarray(time_d, dtype=np.float64)
    good = np.isfinite(time)
    if quality is not None:
        good &= np.asarray(quality) == 0
    if finite_mask is not None:
        good &= np.asarray(finite_mask, dtype=bool)
    if not np.any(good):
        raise ValueError("no observed good cadences available for injection")

    good_time = time[good]
    lo = float(np.nanmin(good_time))
    hi = float(np.nanmax(good_time))
    for _ in range(int(max_tries)):
        t0 = float(rng.uniform(lo, hi))
        mask = box_transit_mask(
            time,
            period_d=period_d,
            t0_d=t0,
            duration_min=duration_min,
        )
        if int(np.count_nonzero(mask & good)) >= min_in_transit:
            return t0, mask

    # Guaranteed to hit at least the chosen good cadence. It may still fail the
    # requested minimum for very sparse data or very short events.
    t0 = float(rng.choice(good_time))
    mask = box_transit_mask(
        time,
        period_d=period_d,
        t0_d=t0,
        duration_min=duration_min,
    )
    n_good = int(np.count_nonzero(mask & good))
    if n_good < min_in_transit:
        raise ValueError(
            f"could not place injection with {min_in_transit} good in-transit "
            f"cadences after {max_tries} tries; best fallback has {n_good}"
        )
    return t0, mask
