"""Known-answer tests for gut_ibm_tools.analysis metrics."""

from __future__ import annotations

import numpy as np

from gut_ibm_tools import analysis


def test_nearest_neighbor_distances_two_clusters() -> None:
    positions = np.array([
        [0.0, 0.0, 0.0],
        [10e-6, 0.0, 0.0],
    ])
    types = np.array([1, 2], dtype=np.int32)

    nnd = analysis.nearest_neighbor_distances(positions, types)
    assert len(nnd) == 2
    np.testing.assert_allclose(nnd[1][0], 10e-6, rtol=1e-6)
    np.testing.assert_allclose(nnd[2][0], 10e-6, rtol=1e-6)


def test_monochromatic_patch_score_segregated() -> None:
    # Two tight monochromatic clusters separated in space.
    rng = np.random.default_rng(0)
    cluster_a = rng.uniform(0, 5e-6, size=(30, 3))
    cluster_b = rng.uniform(50e-6, 55e-6, size=(30, 3))
    positions = np.vstack([cluster_a, cluster_b])
    types = np.array([1] * 30 + [2] * 30, dtype=np.int32)

    score = analysis.monochromatic_patch_score(positions, types, radius=15e-6)
    assert score > 0.85


def test_monochromatic_patch_score_mixed_lower() -> None:
    rng = np.random.default_rng(1)
    positions = rng.uniform(0, 50e-6, size=(60, 3))
    types = np.array([1, 2] * 30, dtype=np.int32)

    score = analysis.monochromatic_patch_score(positions, types, radius=10e-6)
    assert score < 0.75


def test_exclusion_radius_known_geometry() -> None:
    positions = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [5e-6, 0.0, 0.0],
        [5e-6, 0.0, 0.0],
    ])
    types = np.array([1, 1, 2, 2], dtype=np.int32)

    radius = analysis.exclusion_radius(positions, types, target_type=1)
    np.testing.assert_allclose(radius, 5e-6, rtol=1e-6)


def test_comet_tail_index_downstream_bias() -> None:
    positions = np.column_stack([
        np.linspace(-1e-6, 1e-6, 10),
        np.zeros(10),
        np.zeros(10),
    ])
    concentrations = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 2.0, 3.0, 4.0, 5.0, 6.0])

    ratio = analysis.comet_tail_index(positions, concentrations)
    assert ratio > 1.5


def test_hopkins_statistic_clustered_vs_random() -> None:
    rng = np.random.default_rng(123)
    clustered = rng.normal(loc=25e-6, scale=2e-6, size=(80, 3))
    clustered = np.clip(clustered, 0, None)

    # Hopkins uses np.random internally — fix seed for reproducibility.
    np.random.seed(42)
    h_clustered = analysis.hopkins_statistic(clustered)

    np.random.seed(42)
    uniform = rng.uniform(0, 50e-6, size=(80, 3))
    np.random.seed(42)
    h_uniform = analysis.hopkins_statistic(uniform)

    assert h_clustered > h_uniform


def test_comet_tail_asymmetry_index() -> None:
    positions = np.column_stack([
        np.linspace(-1e-6, 1e-6, 8),
        np.zeros(8),
        np.zeros(8),
    ])
    concentrations = np.array([0.1, 0.2, 0.3, 0.4, 2.0, 3.0, 4.0, 5.0])

    asym = analysis.comet_tail_asymmetry_index(positions, concentrations, flow_direction=0)
    assert asym > 1.0
