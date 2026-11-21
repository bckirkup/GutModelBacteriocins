"""Tests for Kirkup & Riley (2004) fecal-pellet LOD censoring & cage pooling (issue #474)."""

from __future__ import annotations

import math

import numpy as np
import pytest

from gut_ibm_tools import pellet_observation
from gut_ibm_tools.pellet_observation import (
    PelletAssayConfig,
    analytical_detection_probability,
    censor_pellet_sample,
    pool_cage_samples,
)


def test_analytical_detection_probability_kirkup_riley_36_colonies() -> None:
    # At 1% true fecal abundance and 36 colonies plated, P(>=1 colony) = 1 - 0.99^36 = 0.303587
    probs = analytical_detection_probability([0.0, 0.01, 0.10, 1.0], n_draws=36)
    assert math.isclose(float(probs[0]), 0.0, abs_tol=1e-12)
    assert math.isclose(float(probs[1]), 1.0 - (0.99**36), rel_tol=1e-9)
    assert math.isclose(float(probs[2]), 1.0 - (0.90**36), rel_tol=1e-9)
    assert math.isclose(float(probs[3]), 1.0, abs_tol=1e-12)

    probs_two = analytical_detection_probability([0.05], n_draws=36, min_positive_draws=2)
    expected_two = 1.0 - ((0.95**36) + 36.0 * 0.05 * (0.95**35))
    assert math.isclose(float(probs_two[0]), expected_two, rel_tol=1e-9)


def test_censor_pellet_sample_reproducible_and_non_mutating() -> None:
    raw = np.array([995, 5, 0], dtype=np.int64)
    raw_copy = raw.copy()
    cfg = PelletAssayConfig(n_draws=36, lod_fraction=0.05, hosts_per_cage=3)

    res1 = censor_pellet_sample(raw, cfg, rng=np.random.default_rng(42))
    res2 = censor_pellet_sample(raw, cfg, rng=np.random.default_rng(42))
    np.testing.assert_array_equal(raw, raw_copy)
    np.testing.assert_array_equal(res1["sampled_counts"], res2["sampled_counts"])
    assert math.isclose(float(np.sum(res1["censored_fractions"])), 1.0, rel_tol=1e-12)
    # Strain 2 has 0 true count -> never false extinction
    assert bool(res1["false_extinction_mask"][2]) is False


def test_pool_cage_samples_three_mice_minority_censoring() -> None:
    # 3 mice/cage: mouse 0 carries a 1.5% producer minority, mice 1 & 2 carry 100% sensitive
    # Cage-level true producer fraction is 0.5% -> censored in most 36-colony platings
    cage_counts = np.array([
        [15.0, 985.0, 0.0],
        [0.0, 1000.0, 0.0],
        [0.0, 1000.0, 0.0],
    ])
    cfg = PelletAssayConfig(n_draws=36, lod_fraction=0.02, hosts_per_cage=3)
    rng = np.random.default_rng(2026)
    false_extinctions = 0
    n_trials = 200
    for _ in range(n_trials):
        pooled = pool_cage_samples(cage_counts, cfg, rng=rng)
        if bool(pooled["false_extinction_mask"][0]):
            false_extinctions += 1
    # With true p = 0.005 and 36 draws, P(0 colonies) = 0.995^36 = 83.5%
    assert false_extinctions > 140
    assert pooled["n_hosts_pooled"] == 3
    np.testing.assert_array_equal(pooled["hosts_carrying_strain"], np.array([1, 3, 0]))


def test_pellet_validation_rejects_invalid_inputs() -> None:
    with pytest.raises(ValueError, match="n_draws"):
        PelletAssayConfig(n_draws=0)
    with pytest.raises(ValueError, match="lod_fraction"):
        PelletAssayConfig(lod_fraction=1.0)
    with pytest.raises(ValueError, match="non-negative"):
        pellet_observation.censor_pellet_sample([10, -1])
