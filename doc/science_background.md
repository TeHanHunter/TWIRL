# Science Background for TWIRL

**Status:** maintained scientific framing, reviewed `2026-07-13`.

This document records the scientific context that is safe to reuse in TWIRL
plans, talks, and pipeline documentation. It is intentionally cautious about
sample size, completeness, and expected yield. The older
[HTML reading guide](wd_background_review.html) is retained only as an April
2026 snapshot and is not authoritative.

## Scientific Motivation and Signal Families

White-dwarf planetary systems can reveal both intact survivors of stellar
evolution and smaller bodies being disrupted near the star. TWIRL therefore
needs two complementary discovery branches:

- **Periodic, short-duration occultations.** WD 1856 b is the benchmark for a
  large intact planet on a close orbit.
- **Variable or weakly periodic dips.** WD 1145+017 was observed with multiple
  `4.5-4.9 h` signals, variable depths, and asymmetric profiles interpreted as
  disintegrating planetesimals and dusty effluent
  ([Vanderburg et al. 2015](https://doi.org/10.1038/nature15527)). This is a
  distinct search problem from a fixed-depth periodic box.

The compact host makes large occultation depths possible, while minute-scale
events can be strongly diluted by long integrations. TESS `200 s` full-frame
images are consequently well matched to TWIRL's initial large, deep,
short-duration regime. That qualitative advantage is not itself a detection
efficiency: sensitivity still depends on target brightness, crowding, cadence
retention, aperture behavior, detrending, search thresholds, and vetting.

## Gaia Catalogue Counts and the Survey Denominator

The Gentile Fusillo et al. Gaia EDR3 main catalogue contains two numbers that
must not be conflated:

- `1,280,266` Gaia EDR3 sources survived the broad colour, absolute-magnitude,
  and Gaia-quality filtering used to isolate the white-dwarf locus. These rows
  were assigned a white-dwarf probability, `Pwd`. They are **not** 1,280,266
  confirmed or high-confidence white dwarfs.
- Applying the authors' general-purpose recommendation `Pwd > 0.75` selects
  `359,073` high-confidence white-dwarf candidates.

Both values are reported directly in table 1 and the conclusion of
[Gentile Fusillo et al. 2021](https://doi.org/10.1093/mnras/stab2672). The
paper also emphasizes that users can trade completeness against contamination
with `Pwd` and quality flags; the threshold does not turn every selected row
into a spectroscopically confirmed white dwarf.

For TWIRL, `Pwd > 0.75` is the high-confidence **reference selection**, not an
automatic occurrence-rate denominator. The statistical denominator must also
encode the final parent-sample cuts, actual `Sector >= 56` coverage, successful
product generation, usable cadence coverage, and the domain over which
end-to-end recovery has been measured. The broader `1,280,266`-row catalogue
may support exploratory target searches, but exploratory additions do not
silently enter the statistical sample.

## WD 1856 b: Discovery, Confirmation, and the 2026 Update

[Vanderburg et al. 2020](https://doi.org/10.1038/s41586-020-2713-y) reported a
Jupiter-sized **planet candidate** transiting WD 1856+534 every `1.4 d`. At that
time the available mass limit did not exclude every low-mass brown-dwarf
interpretation, so the original paper's candidate terminology was appropriate.

That terminology is now outdated. JWST/MIRI thermal-emission measurements
constrained the companion to no more than six Jupiter masses and confirmed its
planetary nature
([Limbach et al. 2025](https://doi.org/10.3847/2041-8213/adc9ad)). Current TWIRL
documents should therefore call it the **confirmed planet WD 1856 b** while
retaining “planet candidate” only when describing the conclusion of the 2020
discovery paper.

The July 2026 [arXiv posting](https://arxiv.org/abs/2607.01316) of the accepted
Nature manuscript by MacDonald et al. reports a JWST/NIRSpec atmospheric
detection, including retrieval evidence for aerosols and hydrocarbons and an
interpretation involving migration-related reheating
([version of record](https://doi.org/10.1038/s41586-026-10514-7)). This is an
important characterization result, not a new confirmation event. Its retrieved
abundances, mass range, effective temperature, and thermal-history
interpretation are model-dependent results from that analysis; they should be
cited to the paper rather than generalized into TWIRL population assumptions.

WD 1856 b remains TWIRL's mandatory engineering benchmark because it tests
whether the real extraction, detrending, aperture, search, and vetting stack can
recover a known short-duration event. Recovering it demonstrates pipeline
function in one favorable system; it does not by itself establish survey-wide
completeness.

## Closest Prior Population Study

[Robert et al. 2024](https://doi.org/10.1093/mnras/stae1859) searched TESS
light curves for `313` metal-polluted white dwarfs, alongside high-cadence
ULTRACAM and ULTRASPEC observations. Their TESS search blindly recovered two
previously known irregularly transiting systems and found no new detections,
giving a nominal detectable-transit fraction of
`0.8%` (`-0.4%`, `+0.6%`). They also performed injection-recovery tests for all
light curves and reported occurrence upper limits across the parameter space
they tested (`1 h` to `27 d`, from dwarf- to Kronian-planet radii).

This study is the most direct population-level comparison for TWIRL, but its
rate must not be multiplied by the Gaia catalogue size to forecast TWIRL
detections. Robert et al. selected **metal-polluted** white dwarfs and used a
different mixture of cadence, time coverage, target properties, and detection
sensitivity. Conversely, the `359,073` high-confidence Gaia candidates are not
all polluted, are not all observed in the production sectors, and will not all
have equivalent recovery probability.

## Claim Boundaries for TWIRL

Safe current framing:

- TWIRL is building a systematic search of `200 s` TESS full-frame data for
  transiting or occulting objects around Gaia-selected white-dwarf candidates.
- Its first-year emphasis is on large, deep, short-duration events, with both a
  transparent periodic search and a separate dip-search branch.
- The project can compare its eventual sensitivity and results with prior K2,
  TESS, and high-speed surveys after its own end-to-end recovery surface is
  measured.

Claims to defer until the relevant products are locked:

- Do not describe all `359,073` high-confidence candidates as uniformly
  observed or searchable at `200 s` cadence.
- Do not claim an expected number of detections by scaling a literature transit
  fraction, an assumed occurrence rate, or the seed catalogue size.
- Do not call TWIRL uniquely capable, the first survey with statistical power,
  or a thousand-fold sensitivity improvement without a like-for-like selection
  and recovery comparison.
- Do not publish occurrence constraints, final null results, or discovery
  claims before the parent sample, coverage, end-to-end completeness, candidate
  merging, and validation are frozen.

Any future yield statement should state its assumed population model, geometric
transit probability, sector coverage, photometric quality, injection domain,
and full search-to-vetting recovery. Until then, expected yield is deliberately
left unspecified.

## Primary References

- [Gentile Fusillo et al. (2021), *A catalogue of white dwarfs in Gaia
  EDR3*](https://doi.org/10.1093/mnras/stab2672)
- [Vanderburg et al. (2015), *A disintegrating minor planet transiting a white
  dwarf*](https://doi.org/10.1038/nature15527)
- [Vanderburg et al. (2020), *A giant planet candidate transiting a white
  dwarf*](https://doi.org/10.1038/s41586-020-2713-y)
- [Limbach et al. (2025), *Thermal Emission and Confirmation of the Frigid White
  Dwarf Exoplanet WD 1856+534 b*](https://doi.org/10.3847/2041-8213/adc9ad)
- [MacDonald et al. (2026), *Aerosols and hydrocarbons in the atmosphere of a
  white dwarf planet*](https://doi.org/10.1038/s41586-026-10514-7)
- [Robert et al. (2024), *The frequency of transiting planetary systems around
  polluted white dwarfs*](https://doi.org/10.1093/mnras/stae1859)
