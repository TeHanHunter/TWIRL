# plotting_style.md

Publication-facing TWIRL figures should use the shared style module at
`src/twirl/plotting/style.py`. Update the code first, then this note if the
style contract changes.

## Required Pattern

- Import and call `apply_twirl_style(template)` instead of defining script-local
  Seaborn or Matplotlib theme blocks.
- Use `get_ordered_palette(n)` for ordered color sets unless a figure needs a
  scientifically meaningful alternative colormap.
- Save publication-facing figures as both PNG and PDF when practical.
- Rasterize only dense scatter layers, not complete figures.

## Default Figure Style

- Use the shared templates such as `column` and `full_page` for figure size.
- Default style is publication-oriented: serif text, white background, light
  grid, dark axes, and restrained annotations.
- Do not add figure titles by default. Add titles only when the figure would
  otherwise be ambiguous.
- For vertically stacked comparison figures, prefer per-panel axis labels,
  concise panel labels when needed, tight inter-panel spacing, and compact
  nearby legends or colorbars.
- Integer-valued color quantities should use a stepped colorbar, not a
  continuous ramp.
- For Aitoff sky maps, use black dotted longitude guides, a black dashed
  Galactic-plane line, and longitude labels with a thin white bezel.
- Legend borders should be black; place legends outside the data area when
  possible.

## Current Plotting Cautions

- For publication figures, do not use script-local style blocks that drift from
  `src/twirl/plotting/style.py`.
- Keep generated reports under `reports/`.
- Avoid figure-level decoration that hides the science comparison, especially
  for QA figures and aperture-to-aperture comparisons.
