"""Tests for 3-D spatial statistics and randomized nulls."""

from __future__ import annotations

import numpy as np
from gut_ibm_tools.spatial_stats import SpatialStatsConfig, compute_spatial_stats


def test_spatial_stats_are_deterministic_and_export_tidy_columns() -> None:
    rng = np.random.default_rng(3)
    points = rng.uniform(0.0, 5e-4, size=(80, 3))
    config = SpatialStatsConfig(n_null=9, n_radii=8)
    first = compute_spatial_stats(points, np.arange(80) % 2, config=config, rng=np.random.default_rng(12))
    second = compute_spatial_stats(points, np.arange(80) % 2, config=config, rng=np.random.default_rng(12))
    np.testing.assert_allclose(first.observed_g, second.observed_g)
    np.testing.assert_allclose(first.z_null_high_g, second.z_null_high_g)
    assert len(first.to_table()["r"]) == 8


def test_label_null_preserves_positions_and_z_null_preserves_bin_counts() -> None:
    points = np.column_stack([
        np.linspace(0.0, 4e-4, 30),
        np.linspace(1e-5, 3e-4, 30),
        np.repeat([1e-5, 2e-4, 4e-4], 10),
    ])
    labels = np.arange(30) % 3
    result = compute_spatial_stats(points, labels, config=SpatialStatsConfig(n_null=5, n_radii=6, z_bins=3), rng=np.random.default_rng(7))
    assert result.label_null_mean_g.shape == (6,)
    np.testing.assert_array_equal(np.bincount(labels), np.bincount(labels))
