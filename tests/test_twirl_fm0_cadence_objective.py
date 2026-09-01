from __future__ import annotations

import math

import pytest

torch = pytest.importorskip("torch")

from twirl.models.fm0.cadence_objective import (
    position_centered_cadence_vicreg_loss,
)


def test_position_centering_does_not_treat_time_offsets_as_batch_variance() -> None:
    # Each fixed cadence position is constant across the batch, even though
    # values differ strongly between positions.  Pooling batch and time would
    # incorrectly treat this as a high-variance representation.
    cadence_values = torch.tensor([-10.0, 10.0]).view(1, 2, 1)
    first = cadence_values.expand(4, -1, -1).clone().requires_grad_()
    second = first.detach().clone().requires_grad_()
    visible = torch.ones((4, 2), dtype=torch.bool)

    loss, diagnostics = position_centered_cadence_vicreg_loss(
        first,
        second,
        visible,
        visible,
        invariance_weight=0.0,
        covariance_weight=0.0,
    )

    torch.testing.assert_close(diagnostics["variance"], torch.tensor(0.99))
    torch.testing.assert_close(loss, torch.tensor(24.75))
    assert diagnostics["invariance"] == 0
    assert diagnostics["statistical_cadence_positions"] == 2
    assert diagnostics["statistical_degrees_of_freedom"] == 6


def test_common_visible_mask_and_fp32_gradient_are_stable() -> None:
    first = torch.tensor(
        [
            [[0.0, 1.0], [float("nan"), float("nan")]],
            [[2.0, 3.0], [4.0, 5.0]],
            [[4.0, 5.0], [6.0, 7.0]],
        ],
        dtype=torch.bfloat16,
        requires_grad=True,
    )
    second = torch.tensor(
        [
            [[0.5, 1.5], [2.0, 3.0]],
            [[2.5, 3.5], [4.5, 5.5]],
            [[4.5, 5.5], [6.5, 7.5]],
        ],
        dtype=torch.bfloat16,
        requires_grad=True,
    )
    first_visible = torch.tensor(
        [[True, False], [True, True], [True, True]], dtype=torch.bool
    )
    second_visible = torch.tensor(
        [[True, True], [True, True], [True, True]], dtype=torch.bool
    )

    with torch.autocast(device_type="cpu", dtype=torch.bfloat16):
        loss, diagnostics = position_centered_cadence_vicreg_loss(
            first, second, first_visible, second_visible
        )
    loss.backward()

    assert loss.dtype == torch.float32
    assert torch.isfinite(loss)
    assert torch.isfinite(first.grad[first_visible]).all()
    assert torch.isfinite(second.grad[second_visible]).all()
    assert diagnostics["position_batch_counts"].tolist() == [3, 2]
    assert diagnostics["common_visible_tokens"] == 5
    assert diagnostics["statistical_degrees_of_freedom"] == 3


def test_pooled_statistics_use_global_position_centered_df() -> None:
    first = torch.tensor(
        [
            [[0.00, 0.00], [0.00, 0.50]],
            [[0.25, 0.50], [0.50, 0.00]],
            [[0.50, 1.00], [9.00, 9.00]],
        ]
    )
    second = first + 0.1
    first_visible = torch.tensor(
        [[True, True], [True, True], [True, False]], dtype=torch.bool
    )
    second_visible = torch.ones((3, 2), dtype=torch.bool)

    loss, diagnostics = position_centered_cadence_vicreg_loss(
        first, second, first_visible, second_visible
    )

    expected_invariance = 0.01
    expected_variance = 0.5 * (
        1.0 - math.sqrt(0.25 / 3.0 + 1.0e-4)
        + 1.0 - math.sqrt(0.625 / 3.0 + 1.0e-4)
    )
    expected_covariance = 2.0 * (0.125 / 3.0) ** 2
    expected_loss = (
        25.0 * expected_invariance
        + 25.0 * expected_variance
        + expected_covariance
    )
    torch.testing.assert_close(
        diagnostics["invariance"], torch.tensor(expected_invariance)
    )
    torch.testing.assert_close(
        diagnostics["variance"], torch.tensor(expected_variance)
    )
    torch.testing.assert_close(
        diagnostics["covariance"], torch.tensor(expected_covariance)
    )
    torch.testing.assert_close(loss, torch.tensor(expected_loss))
    assert diagnostics["statistical_degrees_of_freedom"] == 3


@pytest.mark.parametrize(
    ("first_visible", "second_visible", "message"),
    [
        (
            [[True, False], [False, False]],
            [[True, False], [False, False]],
            "insufficient degrees of freedom",
        ),
        (
            [[False, False], [False, False]],
            [[False, False], [False, False]],
            "no common-visible tokens",
        ),
    ],
)
def test_cadence_vicreg_fails_closed_without_position_df(
    first_visible: list[list[bool]],
    second_visible: list[list[bool]],
    message: str,
) -> None:
    first = torch.zeros((2, 2, 3))
    second = torch.zeros_like(first)

    with pytest.raises(ValueError, match=message):
        position_centered_cadence_vicreg_loss(
            first,
            second,
            torch.tensor(first_visible, dtype=torch.bool),
            torch.tensor(second_visible, dtype=torch.bool),
        )


def test_singleton_position_is_ignored_by_pooled_statistics() -> None:
    first = torch.tensor(
        [
            [[100.0, -100.0], [0.0, 1.0]],
            [[0.0, 0.0], [2.0, 3.0]],
            [[0.0, 0.0], [4.0, 5.0]],
            [[0.0, 0.0], [6.0, 7.0]],
        ]
    )
    second = first.clone()
    visible = torch.tensor(
        [[True, True], [False, True], [False, True], [False, True]],
        dtype=torch.bool,
    )

    loss, diagnostics = position_centered_cadence_vicreg_loss(
        first, second, visible, visible
    )

    assert torch.isfinite(loss)
    assert diagnostics["position_batch_counts"].tolist() == [1, 4]
    assert diagnostics["statistical_cadence_positions"] == 1
    assert diagnostics["statistical_degrees_of_freedom"] == 3


def test_cadence_vicreg_rejects_nonfinite_common_visible_tokens() -> None:
    first = torch.zeros((3, 2, 3))
    first[0, 0, 0] = float("nan")
    second = torch.zeros_like(first)
    visible = torch.ones((3, 2), dtype=torch.bool)

    with pytest.raises(ValueError, match="nonfinite common-visible token"):
        position_centered_cadence_vicreg_loss(
            first, second, visible, visible
        )
