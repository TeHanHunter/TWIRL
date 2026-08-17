"""Small, parameter-matched TWIRL-FM0.1 TCN and Conformer encoders."""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Literal

from .dataset import WINDOW_LENGTH, variant_view_indices

try:  # Torch is an optional platform dependency for the base TWIRL package.
    import torch
    from torch import nn
    from torch.nn import functional as F
except ImportError as exc:  # pragma: no cover - local CPU environment omits Torch
    torch = None  # type: ignore[assignment]
    nn = None  # type: ignore[assignment]
    F = None  # type: ignore[assignment]
    _TORCH_IMPORT_ERROR: ImportError | None = exc
    _ModuleBase = object
else:
    _TORCH_IMPORT_ERROR = None
    _ModuleBase = nn.Module


Architecture = Literal["tcn", "conformer"]


def require_torch() -> None:
    if _TORCH_IMPORT_ERROR is not None:
        raise RuntimeError("PyTorch is required for TWIRL-FM0 models") from (
            _TORCH_IMPORT_ERROR
        )


@dataclass(frozen=True)
class FM0ModelConfig:
    """Exact implementation realization beneath the frozen FM0.1 contract."""

    architecture: Architecture
    n_flux_views: int
    window_length: int = WINDOW_LENGTH
    d_model: int = 256
    embedding_dim: int = 256
    stem_kernel: int = 15
    dropout: float = 0.1
    tcn_blocks: int = 22
    tcn_kernel: int = 3
    tcn_dilation_cycle: tuple[int, ...] = (1, 2, 4, 8, 16, 32, 64)
    conformer_blocks: int = 6
    conformer_heads: int = 8
    conformer_ff_multiplier: int = 4
    conformer_conv_kernel: int = 31
    patch_stride: int = 4
    minimum_parameters: int = 8_000_000
    maximum_parameters: int = 12_000_000

    def __post_init__(self) -> None:
        if self.architecture not in ("tcn", "conformer"):
            raise ValueError(f"unsupported architecture: {self.architecture!r}")
        if self.n_flux_views not in (2, 4, 6):
            raise ValueError("FM0.1 exposes exactly 2, 4, or 6 flux views")
        if self.window_length <= 0 or self.d_model <= 0:
            raise ValueError("window_length and d_model must be positive")
        if self.embedding_dim != 256 and self.d_model == 256:
            # Tiny unit-test configurations may change both; production remains 256.
            raise ValueError("the frozen production embedding dimension is 256")
        if self.d_model % self.conformer_heads:
            raise ValueError("d_model must be divisible by conformer_heads")
        if self.stem_kernel % 2 == 0 or self.conformer_conv_kernel % 2 == 0:
            raise ValueError("stem and Conformer convolution kernels must be odd")
        if self.tcn_kernel % 2 == 0:
            raise ValueError("TCN kernel must be odd")
        if not self.tcn_dilation_cycle:
            raise ValueError("TCN dilation cycle cannot be empty")

    @property
    def input_channels(self) -> int:
        # flux, visible-valid, artificial-mask, and view-present per view;
        # two error values and masks; two time coordinates; time and boundary.
        return 4 * self.n_flux_views + 8


class _ChannelLayerNorm(_ModuleBase):
    def __init__(self, channels: int) -> None:
        require_torch()
        super().__init__()
        self.norm = nn.LayerNorm(channels)

    def forward(self, values: Any) -> Any:
        return self.norm(values.transpose(1, 2)).transpose(1, 2)


class _CadenceStem(_ModuleBase):
    def __init__(self, config: FM0ModelConfig) -> None:
        require_torch()
        super().__init__()
        channels = config.d_model
        padding = config.stem_kernel // 2
        self.project = nn.Conv1d(config.input_channels, channels, kernel_size=1)
        self.depthwise = nn.Conv1d(
            channels,
            channels,
            kernel_size=config.stem_kernel,
            padding=padding,
            groups=channels,
        )
        self.pointwise = nn.Conv1d(channels, channels, kernel_size=1)
        self.norm = _ChannelLayerNorm(channels)
        self.dropout = nn.Dropout(config.dropout)

    def forward(self, values: Any, token_valid: Any) -> Any:
        mask = token_valid[:, None, :].to(dtype=values.dtype)
        hidden = self.project(values) * mask
        local = self.depthwise(hidden)
        local = self.pointwise(F.gelu(local))
        hidden = self.norm(hidden + self.dropout(local))
        return hidden * mask


class _TCNBlock(_ModuleBase):
    def __init__(
        self,
        channels: int,
        *,
        kernel_size: int,
        dilation: int,
        dropout: float,
    ) -> None:
        require_torch()
        super().__init__()
        padding = dilation * (kernel_size - 1) // 2
        self.conv1 = nn.Conv1d(
            channels,
            channels,
            kernel_size,
            padding=padding,
            dilation=dilation,
        )
        self.conv2 = nn.Conv1d(
            channels,
            channels,
            kernel_size,
            padding=padding,
            dilation=dilation,
        )
        self.norm1 = _ChannelLayerNorm(channels)
        self.norm2 = _ChannelLayerNorm(channels)
        self.dropout = nn.Dropout(dropout)

    def forward(self, values: Any, token_valid: Any) -> Any:
        mask = token_valid[:, None, :].to(dtype=values.dtype)
        hidden = self.conv1(values * mask)
        hidden = self.norm1(F.gelu(hidden)) * mask
        hidden = self.conv2(self.dropout(hidden))
        hidden = self.norm2(hidden)
        return F.gelu(values + self.dropout(hidden)) * mask


class _FeedForward(_ModuleBase):
    def __init__(self, d_model: int, multiplier: int, dropout: float) -> None:
        require_torch()
        super().__init__()
        hidden = d_model * multiplier
        self.layers = nn.Sequential(
            nn.Linear(d_model, hidden),
            nn.SiLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden, d_model),
            nn.Dropout(dropout),
        )

    def forward(self, values: Any) -> Any:
        return self.layers(values)


class _ConformerConv(_ModuleBase):
    def __init__(self, d_model: int, kernel_size: int, dropout: float) -> None:
        require_torch()
        super().__init__()
        self.pointwise_in = nn.Conv1d(d_model, 2 * d_model, kernel_size=1)
        self.depthwise = nn.Conv1d(
            d_model,
            d_model,
            kernel_size=kernel_size,
            padding=kernel_size // 2,
            groups=d_model,
        )
        self.norm = _ChannelLayerNorm(d_model)
        self.pointwise_out = nn.Conv1d(d_model, d_model, kernel_size=1)
        self.dropout = nn.Dropout(dropout)

    def forward(self, values: Any, token_valid: Any) -> Any:
        mask = token_valid[:, None, :].to(dtype=values.dtype)
        hidden = self.pointwise_in(values.transpose(1, 2) * mask)
        # Pointwise bias can be nonzero at padded tokens.  Re-mask before the
        # depthwise convolution so those artificial values cannot leak into a
        # neighboring valid cadence.
        hidden = F.glu(hidden, dim=1) * mask
        hidden = self.depthwise(hidden)
        hidden = self.norm(hidden)
        hidden = self.pointwise_out(F.silu(hidden))
        return self.dropout(hidden).transpose(1, 2) * token_valid[:, :, None]


class _ConformerBlock(_ModuleBase):
    def __init__(self, config: FM0ModelConfig) -> None:
        require_torch()
        super().__init__()
        d_model = config.d_model
        self.norm_ff1 = nn.LayerNorm(d_model)
        self.ff1 = _FeedForward(
            d_model, config.conformer_ff_multiplier, config.dropout
        )
        self.norm_attention = nn.LayerNorm(d_model)
        self.attention = nn.MultiheadAttention(
            d_model,
            config.conformer_heads,
            dropout=config.dropout,
            batch_first=True,
        )
        self.norm_conv = nn.LayerNorm(d_model)
        self.conv = _ConformerConv(
            d_model, config.conformer_conv_kernel, config.dropout
        )
        self.norm_ff2 = nn.LayerNorm(d_model)
        self.ff2 = _FeedForward(
            d_model, config.conformer_ff_multiplier, config.dropout
        )
        self.final_norm = nn.LayerNorm(d_model)

    def forward(self, values: Any, token_valid: Any) -> Any:
        mask = token_valid[:, :, None].to(dtype=values.dtype)
        hidden = values + 0.5 * self.ff1(self.norm_ff1(values))
        attention_input = self.norm_attention(hidden)
        key_padding = ~token_valid
        all_masked = torch.all(key_padding, dim=1)
        if bool(torch.any(all_masked)):
            key_padding = key_padding.clone()
            key_padding[all_masked, 0] = False
            attention_input = attention_input.clone()
            attention_input[all_masked, 0] = 0.0
        attended, _ = self.attention(
            attention_input,
            attention_input,
            attention_input,
            key_padding_mask=key_padding,
            need_weights=False,
        )
        hidden = (hidden + attended) * mask
        hidden = hidden + self.conv(self.norm_conv(hidden), token_valid)
        hidden = hidden + 0.5 * self.ff2(self.norm_ff2(hidden))
        return self.final_norm(hidden) * mask


class _ContinuousTimeEncoding(_ModuleBase):
    def __init__(self, config: FM0ModelConfig) -> None:
        require_torch()
        super().__init__()
        self.window_length = float(config.window_length)
        self.layers = nn.Sequential(
            nn.Linear(2, config.d_model),
            nn.SiLU(),
            nn.Linear(config.d_model, config.d_model),
        )

    def forward(self, local_time: Any, delta_time: Any) -> Any:
        local = local_time / max(1.0, self.window_length - 1.0)
        delta = torch.sign(delta_time) * torch.log1p(torch.abs(delta_time))
        delta = delta / math.log(87.0)
        return self.layers(torch.stack((local, delta), dim=-1))


class TWIRLFM0(_ModuleBase):
    """Window encoder and per-view reconstruction head for one FM0.1 variant."""

    def __init__(self, config: FM0ModelConfig) -> None:
        require_torch()
        super().__init__()
        self.config = config
        self.stem = _CadenceStem(config)
        if config.architecture == "tcn":
            self.tcn = nn.ModuleList(
                _TCNBlock(
                    config.d_model,
                    kernel_size=config.tcn_kernel,
                    dilation=config.tcn_dilation_cycle[
                        index % len(config.tcn_dilation_cycle)
                    ],
                    dropout=config.dropout,
                )
                for index in range(config.tcn_blocks)
            )
            self.time_encoding = None
            self.conformer = None
        else:
            self.tcn = None
            self.time_encoding = _ContinuousTimeEncoding(config)
            self.conformer = nn.ModuleList(
                _ConformerBlock(config) for _ in range(config.conformer_blocks)
            )
        self.embedding_projection = nn.Linear(
            config.d_model, config.embedding_dim
        )
        self.reconstruction_head = nn.Conv1d(
            config.d_model, config.n_flux_views, kernel_size=1
        )

    def _pack_input(self, batch: dict[str, Any]) -> tuple[Any, Any, Any, Any]:
        required = (
            "flux",
            "flux_valid",
            "flux_error",
            "error_valid",
            "local_time_cadences",
            "delta_time_cadences",
            "time_valid",
            "segment_boundary",
            "view_present",
            "reconstruction_mask",
        )
        missing = [name for name in required if name not in batch]
        if missing:
            raise ValueError(f"FM0 batch lacks required tensors: {missing}")
        flux = batch["flux"].float()
        flux_valid = batch["flux_valid"].bool()
        reconstruction_mask = batch["reconstruction_mask"].bool()
        view_present = batch["view_present"].bool()
        if flux.ndim != 3 or flux.shape[1] != self.config.n_flux_views:
            raise ValueError("flux must have shape [batch, configured_view, cadence]")
        if flux.shape[-1] != self.config.window_length:
            raise ValueError("flux cadence length differs from model configuration")
        if flux_valid.shape != flux.shape or reconstruction_mask.shape != flux.shape:
            raise ValueError("flux masks must match flux shape")
        if view_present.shape != flux.shape[:2]:
            raise ValueError("view_present must have shape [batch, view]")

        time_valid = batch["time_valid"].bool()
        if time_valid.shape != (flux.shape[0], flux.shape[2]):
            raise ValueError("time_valid must have shape [batch, cadence]")
        visible = flux_valid & ~reconstruction_mask & view_present[:, :, None]
        safe_flux = torch.where(visible, flux, torch.zeros_like(flux))
        errors = batch["flux_error"].float()
        error_valid = batch["error_valid"].bool()
        if errors.shape != (flux.shape[0], 2, flux.shape[2]):
            raise ValueError("flux_error must have shape [batch, 2, cadence]")
        if error_valid.shape != errors.shape:
            raise ValueError("error_valid must match flux_error")
        safe_errors = torch.where(error_valid, errors, torch.zeros_like(errors))
        local_time = batch["local_time_cadences"].float()
        delta_time = batch["delta_time_cadences"].float()
        boundary = batch["segment_boundary"].bool()
        for name, value in (
            ("local_time_cadences", local_time),
            ("delta_time_cadences", delta_time),
            ("segment_boundary", boundary),
        ):
            if value.shape != time_valid.shape:
                raise ValueError(f"{name} must match time_valid")

        local_scaled = local_time / max(1, self.config.window_length - 1)
        delta_scaled = torch.sign(delta_time) * torch.log1p(torch.abs(delta_time))
        delta_scaled = delta_scaled / math.log(87.0)
        channels = (
            safe_flux,
            visible.to(dtype=flux.dtype),
            reconstruction_mask.to(dtype=flux.dtype),
            view_present[:, :, None]
            .expand(-1, -1, flux.shape[2])
            .to(dtype=flux.dtype),
            safe_errors,
            error_valid.to(dtype=flux.dtype),
            local_scaled[:, None, :],
            delta_scaled[:, None, :],
            time_valid[:, None, :].to(dtype=flux.dtype),
            boundary[:, None, :].to(dtype=flux.dtype),
        )
        packed = torch.cat(channels, dim=1)
        return packed, time_valid, local_time, delta_time

    @staticmethod
    def _masked_mean(values: Any, valid: Any) -> Any:
        weights = valid[:, :, None].to(dtype=values.dtype)
        numerator = torch.sum(values * weights, dim=1)
        denominator = torch.sum(weights, dim=1).clamp_min(1.0)
        return numerator / denominator

    def _patch_mean(self, values: Any, valid: Any) -> tuple[Any, Any]:
        stride = self.config.patch_stride
        length = values.shape[-1]
        padding = (-length) % stride
        if padding:
            values = F.pad(values, (0, padding))
            valid = F.pad(valid, (0, padding), value=False)
        batch, channels, padded = values.shape
        values = values.reshape(batch, channels, padded // stride, stride)
        valid_blocks = valid.reshape(batch, padded // stride, stride)
        weights = valid_blocks[:, None, :, :].to(dtype=values.dtype)
        pooled = torch.sum(values * weights, dim=-1)
        pooled /= torch.sum(weights, dim=-1).clamp_min(1.0)
        return pooled, torch.any(valid_blocks, dim=-1)

    def _patch_time(self, values: Any, valid: Any) -> Any:
        stride = self.config.patch_stride
        padding = (-values.shape[-1]) % stride
        if padding:
            values = F.pad(values, (0, padding))
            valid = F.pad(valid, (0, padding), value=False)
        reshaped = values.reshape(values.shape[0], -1, stride)
        valid_blocks = valid.reshape(valid.shape[0], -1, stride)
        weights = valid_blocks.to(dtype=values.dtype)
        return torch.sum(reshaped * weights, dim=-1) / torch.sum(
            weights, dim=-1
        ).clamp_min(1.0)

    def forward(self, batch: dict[str, Any]) -> dict[str, Any]:
        packed, token_valid, local_time, delta_time = self._pack_input(batch)
        hidden = self.stem(packed, token_valid)
        if self.config.architecture == "tcn":
            for block in self.tcn:
                hidden = block(hidden, token_valid)
            token_hidden = hidden.transpose(1, 2)
            output_valid = token_valid
            full_hidden = hidden
        else:
            patched, output_valid = self._patch_mean(hidden, token_valid)
            patched_local = self._patch_time(local_time, token_valid)
            patched_delta = self._patch_time(delta_time, token_valid)
            token_hidden = patched.transpose(1, 2)
            token_hidden = token_hidden + self.time_encoding(
                patched_local, patched_delta
            )
            token_hidden = token_hidden * output_valid[:, :, None]
            for block in self.conformer:
                token_hidden = block(token_hidden, output_valid)
            upsampled = token_hidden.transpose(1, 2).repeat_interleave(
                self.config.patch_stride, dim=-1
            )[:, :, : self.config.window_length]
            # The pre-patch cadence stem is a local reconstruction skip.  The
            # Conformer context therefore cannot force all four cadences in a
            # patch to share one decoded state and erase sub-patch events.
            full_hidden = (upsampled + hidden) * token_valid[:, None, :]

        embedding = self.embedding_projection(
            self._masked_mean(token_hidden, output_valid)
        )
        reconstruction = self.reconstruction_head(full_hidden)
        reconstruction = reconstruction * token_valid[:, None, :]
        return {
            "reconstruction": reconstruction,
            "z_window": embedding,
            "token_valid": output_valid,
        }


def architecture_for_variant(
    variant: str,
    *,
    development_winner: Architecture | None = None,
) -> Architecture:
    if variant == "TWIRL-FM0.1.1":
        return "tcn"
    if variant == "TWIRL-FM0.1.2":
        return "conformer"
    variant_view_indices(variant)
    if development_winner not in ("tcn", "conformer"):
        raise ValueError(
            f"{variant} requires the frozen development winner architecture"
        )
    return development_winner


def build_fm0_model(
    variant: str,
    *,
    development_winner: Architecture | None = None,
    config_override: FM0ModelConfig | None = None,
    enforce_parameter_budget: bool = True,
) -> TWIRLFM0:
    """Build one frozen-ladder model and enforce its production parameter range."""

    indices = variant_view_indices(variant)
    architecture = architecture_for_variant(
        variant, development_winner=development_winner
    )
    if config_override is None:
        config = FM0ModelConfig(
            architecture=architecture,
            n_flux_views=len(indices),
        )
    else:
        config = config_override
        if config.architecture != architecture or config.n_flux_views != len(indices):
            raise ValueError("model override conflicts with frozen variant")
    model = TWIRLFM0(config)
    if enforce_parameter_budget:
        count = count_trainable_parameters(model)
        if not config.minimum_parameters <= count <= config.maximum_parameters:
            raise ValueError(
                f"{architecture} parameter count {count:,} is outside the "
                f"frozen {config.minimum_parameters:,}--{config.maximum_parameters:,} range"
            )
    return model


def count_trainable_parameters(model: Any) -> int:
    require_torch()
    return int(sum(parameter.numel() for parameter in model.parameters() if parameter.requires_grad))


def parameter_match_fraction(first: Any, second: Any) -> float:
    """Return relative parameter-count difference against the larger model."""

    first_count = count_trainable_parameters(first)
    second_count = count_trainable_parameters(second)
    return abs(first_count - second_count) / max(first_count, second_count, 1)
