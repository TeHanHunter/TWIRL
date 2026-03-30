from __future__ import annotations

from copy import deepcopy

import seaborn as sns


PLOT_TEMPLATES = {
    "column": {
        "figsize": (3.4, 2.65),
        "title_size": 8,
        "label_size": 7,
        "tick_size": 6,
        "legend_size": 6,
        "legend_title_size": 6.5,
        "annotation_size": 6,
        "dense_marker_size": 1.25,
        "grid_linewidth": 0.55,
        "panel_wspace": 0.05,
    },
    "full_page": {
        "figsize": (7.1, 4.1),
        "title_size": 9,
        "label_size": 8,
        "tick_size": 7,
        "legend_size": 7,
        "legend_title_size": 7.5,
        "annotation_size": 7,
        "dense_marker_size": 1.7,
        "grid_linewidth": 0.6,
        "panel_wspace": 0.08,
    },
}

BASE_RC = {
    "font.family": "serif",
    "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif",
    "axes.facecolor": "white",
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "axes.edgecolor": "0.25",
    "axes.labelcolor": "0.15",
    "axes.linewidth": 0.8,
    "axes.titleweight": "regular",
    "axes.labelweight": "regular",
    "text.color": "0.15",
    "xtick.color": "0.2",
    "ytick.color": "0.2",
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "grid.color": "0.88",
    "grid.alpha": 1.0,
    "legend.frameon": True,
    "legend.facecolor": "white",
    "legend.edgecolor": "black",
    "legend.fancybox": False,
}


def get_ordered_palette(n_colors: int, palette: str = "viridis"):
    return sns.color_palette(palette, n_colors)


def apply_twirl_style(template_name: str = "column") -> dict[str, float]:
    if template_name not in PLOT_TEMPLATES:
        raise ValueError(f"Unknown TWIRL plot template: {template_name}")

    template = deepcopy(PLOT_TEMPLATES[template_name])
    rc = {
        **BASE_RC,
        "font.size": template["tick_size"],
        "axes.titlesize": template["title_size"],
        "axes.labelsize": template["label_size"],
        "xtick.labelsize": template["tick_size"],
        "ytick.labelsize": template["tick_size"],
        "legend.fontsize": template["legend_size"],
        "legend.title_fontsize": template["legend_title_size"],
        "grid.linewidth": template["grid_linewidth"],
    }
    sns.set_theme(style="whitegrid", context="paper", rc=rc)
    return template
