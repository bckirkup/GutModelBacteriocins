"""Tests for 3-D spatial statistics and randomized nulls."""

from __future__ import annotations

import numpy as np
import pytest

from gut_ibm_tools.spatial_stats import (
    SpatialStatsConfig,
    compute_spatial_stats,
    label_permutation_null,
    summarize_excess,
    z_stratified_null,
)


def test_spatial_stats_are_deterministic_and_export_tidy_columns() -> None:
    rng = np.random.default_rng(3)
    points = rng.uniform(0.0, 5e-4, size=(80, 3))
    config = SpatialStatsConfig(n_null=9, n_radii=8)
    first = compute_spatial_stats(points, np.arange(80) % 2, config=config, rng=np.random.default_rng(12))
    second = compute_spatial_stats(points, np.arange(80) % 2, config=config, rng=np.random.default_rng(12))
    np.testing.assert_allclose(first.observed_g, second.observed_g)
    np.testing.assert_allclose(first.z_null_high_g, second.z_null_high_g)
    assert len(first.to_table()["r"]) == 8


def test_null_generators_preserve_their_constraints() -> None:
    points = np.column_stack([
        np.linspace(0.0, 4e-4, 30),
        np.linspace(1e-5, 3e-4, 30),
        np.repeat([1e-5, 2e-4, 4e-4], 10),
    ])
    labels = np.arange(30) % 3
    z_null = z_stratified_null(points, 3, np.random.default_rng(7))
    np.testing.assert_array_equal(z_null[:, 2], points[:, 2])
    np.testing.assert_array_equal(
        np.histogram(z_null[:, 2], bins=3, range=(points[:, 2].min(), points[:, 2].max()))[0],
        np.histogram(points[:, 2], bins=3, range=(points[:, 2].min(), points[:, 2].max()))[0],
    )
    z_range = np.ptp(points[:, 2])
    for index in range(3):
        members = (points[:, 2] >= points[:, 2].min() + index * (z_range / 3)) & (
            (points[:, 2] < points[:, 2].min() + (index + 1) * (z_range / 3))
            | (index == 2)
        )
        assert np.all(z_null[members, :2] <= points[members, :2].max(axis=0) + 1e-15)
        assert np.all(z_null[members, :2] >= points[members, :2].min(axis=0) - 1e-15)
    null_positions, shuffled = label_permutation_null(
        points, labels, np.random.default_rng(9)
    )
    np.testing.assert_array_equal(null_positions, points)
    np.testing.assert_array_equal(shuffled.shape, labels.shape)
    np.testing.assert_array_equal(np.sort(shuffled), np.sort(labels))


def test_label_null_is_invalid_for_unmarked_points() -> None:
    points = np.random.default_rng(4).uniform(0, 5e-4, (100, 3))
    result = compute_spatial_stats(points, None, config=SpatialStatsConfig(n_null=5), rng=np.random.default_rng(5))
    assert not result.label_null_valid
    with pytest.raises(ValueError, match="label null is invalid"):
        summarize_excess(result, null="label")


def test_poisson_pair_correlation_is_calibrated() -> None:
    points = np.random.default_rng(12).uniform(0.0, 5e-4, size=(1000, 3))
    result = compute_spatial_stats(
        points,
        config=SpatialStatsConfig(
            r_min=1e-5, r_max=5e-5, n_radii=5, n_null=99,
        ),
        rng=np.random.default_rng(2),
    )
    np.testing.assert_allclose(result.observed_g[1:], 1.0, atol=0.3)
    assert np.all(
        (result.observed_g >= result.z_null_low_g)
        & (result.observed_g <= result.z_null_high_g)
    )


def test_all_equal_labels_reduce_to_unlabeled_statistic() -> None:
    points = np.random.default_rng(15).uniform(0.0, 5e-4, size=(300, 3))
    config = SpatialStatsConfig(n_null=5, r_min=1e-5, r_max=5e-5, n_radii=5)
    unlabeled = compute_spatial_stats(points, None, config=config, rng=np.random.default_rng(2))
    marked = compute_spatial_stats(
        points, np.ones(len(points)), config=config, rng=np.random.default_rng(2)
    )
    np.testing.assert_allclose(marked.observed_g, unlabeled.observed_g)
    np.testing.assert_allclose(marked.observed_k, unlabeled.observed_k)


def test_clustered_points_exceed_z_null_with_margin() -> None:
    rng = np.random.default_rng(4)
    centers = rng.uniform(0.0, 5e-4, size=(6, 3))
    points = np.vstack([center + rng.normal(0.0, 8e-6, size=(20, 3)) for center in centers])
    result = compute_spatial_stats(
        points,
        config=SpatialStatsConfig(r_min=2e-6, r_max=1e-4, n_radii=10, n_null=99),
        rng=np.random.default_rng(2),
    )
    assert result.observed_g[1] > 3.0 * result.z_null_high_g[1]
    assert result.observed_g[2] > 2.0 * result.z_null_high_g[2]
