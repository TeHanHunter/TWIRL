from __future__ import annotations

import numpy as np
import pytest

from twirl.plotting.style import safe_log_contour_label_position


class _Contour:
    def __init__(self, vertices: np.ndarray) -> None:
        self.allsegs = [[vertices]]


def test_safe_log_contour_label_position_avoids_panel_edges() -> None:
    contour = _Contour(
        np.array(
            [
                [0.121, 0.19],
                [0.2, 0.4],
                [1.0, 1.0],
                [5.0, 5.0],
                [12.9, 17.9],
            ]
        )
    )
    position = safe_log_contour_label_position(
        contour,
        preferred_xy=(13.0, 18.0),
        xlim=(0.12, 13.0),
        ylim=(0.18, 18.0),
    )
    assert position == (5.0, 5.0)


def test_safe_log_contour_label_position_omits_boundary_only_contour() -> None:
    contour = _Contour(np.array([[0.121, 0.19], [12.9, 17.9]]))
    assert (
        safe_log_contour_label_position(
            contour,
            preferred_xy=(1.0, 1.0),
            xlim=(0.12, 13.0),
            ylim=(0.18, 18.0),
        )
        is None
    )
    with pytest.raises(ValueError, match="margins"):
        safe_log_contour_label_position(
            contour,
            preferred_xy=(1.0, 1.0),
            xlim=(0.12, 13.0),
            ylim=(0.18, 18.0),
            x_margin_fraction=0.5,
        )
