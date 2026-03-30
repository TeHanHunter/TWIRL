# plotting_style.md

This document records the authoritative TWIRL plotting setup.

The actual style definition now lives in [`style.py`](/Users/tehan/PycharmProjects/TWIRL/src/twirl/plotting/style.py). Figure scripts should import and use that module instead of defining their own Seaborn `rc` blocks.

## Default Look

The shared TWIRL style is intentionally:

- light
- serif
- low-contrast but crisp
- publication-oriented rather than slide-oriented

The chosen base is:

- `seaborn` theme: `whitegrid`
- `seaborn` context: `paper`
- font family: `DejaVu Serif`
- math font set: `dejavuserif`
- white backgrounds
- light grey grid lines
- dark grey axes and text
- ordered palette helper from `viridis`

This keeps the earlier serif look while making the figure system explicit and reusable.

## Authoritative API

Use:

```python
from twirl.plotting.style import apply_twirl_style, get_ordered_palette

template = apply_twirl_style("full_page")
palette = get_ordered_palette(4)
```

The named templates currently defined in code are:

- `column`
- `full_page`

The current `Pwd` spatial-position figure should use `full_page`.

## Template Intent

- `column`: single-column manuscript figures and compact diagnostics
- `full_page`: two-panel overview plots and figures intended to span the page width

Both templates keep title, label, tick, legend, and annotation sizes internally consistent so different figures still look like they belong to the same paper.

## Figure Rules

- Do not hardcode Seaborn theme settings inside individual plotting scripts unless there is a figure-specific reason.
- Keep shared legends outside the data area when possible.
- The legend border should be black.
- Prefer tight panel spacing for multi-panel figures.
- Save publication figures as both `PDF` and `PNG`.
- Rasterize only the dense scatter layers, not the full figure.

## Maintenance Rule

If the manuscript visual standard changes, update [`style.py`](/Users/tehan/PycharmProjects/TWIRL/src/twirl/plotting/style.py) first and then update this file to match it.
