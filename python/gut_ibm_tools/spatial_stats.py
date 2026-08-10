"""Composition-aware 3-D pair-correlation and Ripley spatial statistics."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from .path_utils import prepare_output_file


@dataclass(frozen=True)
class SpatialStatsConfig:
    """Radius grid and randomized-null settings."""

    radii: np.ndarray | None = None
    r_min: float = 1e-5
    r_max: float = 5e-4
    n_radii: int = 24
    n_null: int = 99
    z_bins: int = 8

    def radius_grid(self) -> np.ndarray:
        return self.radii.copy() if self.radii is not None else np.linspace(self.r_min, self.r_max, self.n_radii)


@dataclass
class SpatialStatsResult:
    """Observed same-mark curves and pointwise envelopes for both null models.

    ``observed_g`` and ``observed_k`` count only pairs with equal labels,
    normalized by whole-population intensity and size.  With one mark they
    are the ordinary unmarked estimators.
    """

    radii: np.ndarray
    observed_g: np.ndarray
    observed_k: np.ndarray
    observed_l_minus_r: np.ndarray
    z_null_mean_g: np.ndarray
    z_null_low_g: np.ndarray
    z_null_high_g: np.ndarray
    z_null_mean_k: np.ndarray
    z_null_low_k: np.ndarray
    z_null_high_k: np.ndarray
    label_null_mean_g: np.ndarray
    label_null_low_g: np.ndarray
    label_null_high_g: np.ndarray
    label_null_mean_k: np.ndarray
    label_null_low_k: np.ndarray
    label_null_high_k: np.ndarray
    label_null_valid: bool

    def to_table(self) -> dict[str, np.ndarray]:
        """Return a tidy columnar representation suitable for CSV."""
        result: dict[str, np.ndarray] = {"r": self.radii, "observed_g": self.observed_g, "observed_k": self.observed_k, "observed_l_minus_r": self.observed_l_minus_r}
        result.update(
            {
                "z_null_mean_g": self.z_null_mean_g,
                "z_null_low_g": self.z_null_low_g,
                "z_null_high_g": self.z_null_high_g,
                "z_null_mean_k": self.z_null_mean_k,
                "z_null_low_k": self.z_null_low_k,
                "z_null_high_k": self.z_null_high_k,
                "label_null_mean_g": self.label_null_mean_g,
                "label_null_low_g": self.label_null_low_g,
                "label_null_high_g": self.label_null_high_g,
                "label_null_mean_k": self.label_null_mean_k,
                "label_null_low_k": self.label_null_low_k,
                "label_null_high_k": self.label_null_high_k,
                "label_null_valid": np.full(len(self.radii), self.label_null_valid),
            }
        )
        return result

    def write_csv(self, path: str | Path) -> None:
        """Write the tidy curve table."""
        target = prepare_output_file(path)
        table = self.to_table()
        with target.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(table)
            writer.writerows(zip(*(table[column] for column in table)))


def compute_spatial_stats(
    positions: np.ndarray,
    labels: np.ndarray | None = None,
    *,
    config: SpatialStatsConfig | None = None,
    rng: np.random.Generator | None = None,
) -> SpatialStatsResult:
    """Compute 3-D g(r), K(r), and L(r)-r with two composition-matched nulls.

    Intensity is estimated as N divided by the occupied axis-aligned bounding
    box.  No edge correction is applied, so boundary bias is expected; observed
    and both null ensembles use this identical estimator, making their
    comparison valid even when the domain is not rectangular.
    """
    points = _points(positions)
    cfg = config or SpatialStatsConfig()
    generator = rng or np.random.default_rng()
    radii = np.asarray(cfg.radius_grid(), dtype=float)
    if np.any(radii <= 0) or np.any(np.diff(radii) <= 0):
        raise ValueError("radii must be strictly increasing and positive")
    composition = _labels(labels, len(points))
    observed = _curves(points, composition, radii)
    z_curves = np.array(
        [
            _curves(z_stratified_null(points, cfg.z_bins, generator), composition, radii)
            for _ in range(cfg.n_null)
        ]
    )
    label_valid = len(np.unique(composition)) > 1
    label_curves = (
        np.array(
            [
                _curves(
                    points,
                    label_permutation_null(points, composition, generator)[1],
                    radii,
                )
                for _ in range(cfg.n_null)
            ]
        )
        if label_valid
        else None
    )
    z_mean, z_low, z_high = _summaries(z_curves)
    label_summary = _label_null_summaries(label_curves, radii)
    return SpatialStatsResult(
        radii, observed[0], observed[1], observed[2],
        z_mean[0], z_low[0], z_high[0], z_mean[1], z_low[1], z_high[1],
        *label_summary,
        label_valid,
    )


spatial_statistics = compute_spatial_stats


def summarize_excess(result: SpatialStatsResult, *, null: str = "z") -> list[tuple[float, float]]:
    """Return contiguous radius intervals where observed g exits a null band."""
    if null not in {"z", "label"}:
        raise ValueError("null must be 'z' or 'label'")
    if null == "label" and not result.label_null_valid:
        raise ValueError("label null is invalid when fewer than two labels are present")
    low, high = (
        (result.z_null_low_g, result.z_null_high_g)
        if null == "z"
        else (result.label_null_low_g, result.label_null_high_g)
    )
    outside = (result.observed_g < low) | (result.observed_g > high)
    intervals: list[tuple[float, float]] = []
    starts = np.flatnonzero(outside & np.r_[True, ~outside[:-1]])
    ends = np.flatnonzero(outside & np.r_[~outside[1:], True])
    intervals.extend((float(result.radii[start]), float(result.radii[end])) for start, end in zip(starts, ends))
    return intervals


def _points(positions: np.ndarray) -> np.ndarray:
    points = np.asarray(positions, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or len(points) < 2:
        raise ValueError("positions must have shape (N, 3), with N >= 2")
    return points


def _labels(labels: np.ndarray | None, size: int) -> np.ndarray:
    return np.zeros(size, dtype=np.int64) if labels is None else np.asarray(labels)


def _curves(
    points: np.ndarray,
    labels: np.ndarray,
    radii: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return same-mark g, K, and L-r using whole-population intensity."""
    distances = cKDTree(points).query_pairs(float(radii[-1]), output_type="ndarray")
    pair_distances = np.linalg.norm(points[distances[:, 0]] - points[distances[:, 1]], axis=1) if len(distances) else np.array([])
    same = labels[distances[:, 0]] == labels[distances[:, 1]] if len(distances) else np.array([], dtype=bool)
    counts = np.array([np.sum((pair_distances <= radius) & same) for radius in radii], dtype=float) * 2.0
    volume = float(np.prod(np.ptp(points, axis=0)))
    intensity = len(points) / volume if volume > 0 else 0.0
    k_values = counts / (intensity * len(points)) if intensity > 0 else np.zeros_like(counts)
    edges = np.r_[0.0, radii]
    shell = 4.0 * np.pi * ((edges[1:] ** 3 - edges[:-1] ** 3) / 3.0)
    g_values = np.diff(np.r_[0.0, k_values]) / shell if intensity > 0 else np.zeros_like(k_values)
    return g_values, k_values, np.cbrt(3.0 * k_values / (4.0 * np.pi)) - radii


def z_stratified_null(
    points: np.ndarray,
    bins: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Randomize lateral coordinates within realized z strata, preserving z."""
    if bins < 1:
        raise ValueError("bins must be positive")
    result = points.copy()
    edges = np.linspace(float(np.min(points[:, 2])), float(np.max(points[:, 2])), bins + 1)
    for index in range(bins):
        members = np.flatnonzero((points[:, 2] >= edges[index]) & ((points[:, 2] < edges[index + 1]) | (index == bins - 1)))
        if len(members):
            lateral = points[members, :2]
            lo = np.min(lateral, axis=0)
            hi = np.max(lateral, axis=0)
            result[members, :2] = rng.uniform(lo, hi, size=(len(members), 2))
    return result


def label_permutation_null(
    positions: np.ndarray,
    labels: np.ndarray,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Return bit-identical positions and a composition-preserving labels permutation."""
    return np.asarray(positions).copy(), rng.permutation(np.asarray(labels)).copy()


def _summaries(curves: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return np.mean(curves, axis=0), np.percentile(curves, 2.5, axis=0), np.percentile(curves, 97.5, axis=0)


def _label_null_summaries(
    curves: np.ndarray | None,
    radii: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Summarize label-null g and K curves, or return invalid NaN arrays."""
    if curves is None:
        nan = np.full_like(radii, np.nan)
        return nan, nan.copy(), nan.copy(), nan.copy(), nan.copy(), nan.copy()
    g_mean, g_low, g_high = _summaries(curves[:, 0, :])
    k_mean, k_low, k_high = _summaries(curves[:, 1, :])
    return g_mean, g_low, g_high, k_mean, k_low, k_high
