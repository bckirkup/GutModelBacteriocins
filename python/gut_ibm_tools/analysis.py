"""
Spatial clustering analysis and nearest-neighbor distance (NND) metrics.
Exclusion-radius and NND clustering metrics for validation.
"""

from __future__ import annotations

import numpy as np
from scipy.spatial import KDTree


def nearest_neighbor_distances(positions: np.ndarray, types: np.ndarray) -> dict[int, np.ndarray]:
    """
    NND between competing clones.

    For each agent, compute the distance to the nearest agent of a
    *different* type.  Returns dict mapping type → array of distances.
    """
    result: dict[int, np.ndarray] = {}
    unique_types = np.unique(types)

    if len(unique_types) < 2:
        for t in unique_types:
            result[int(t)] = np.array([])
        return result

    for t in unique_types:
        mask = types == t
        pts = positions[mask]
        if len(pts) == 0:
            result[int(t)] = np.array([])
            continue

        other_mask = ~mask
        other_pts = positions[other_mask]
        if len(other_pts) == 0:
            result[int(t)] = np.array([])
            continue

        other_tree = KDTree(other_pts)
        dists, _ = other_tree.query(pts, k=1)
        result[int(t)] = dists.flatten()

    return result


def inter_type_distances(positions: np.ndarray, types: np.ndarray) -> dict[tuple[int, int], np.ndarray]:
    """
    Compute nearest-neighbor distances between agents of different types.
    Returns dict mapping (type_a, type_b) → NNDs from a to nearest b.
    """
    result: dict[tuple[int, int], np.ndarray] = {}
    unique_types = sorted(np.unique(types))

    for i, t1 in enumerate(unique_types):
        for t2 in unique_types[i + 1 :]:
            mask1 = types == t1
            mask2 = types == t2
            pts1 = positions[mask1]
            pts2 = positions[mask2]
            if len(pts1) == 0 or len(pts2) == 0:
                continue
            tree2 = KDTree(pts2)
            dists, _ = tree2.query(pts1, k=1)
            result[(int(t1), int(t2))] = dists.flatten()

    return result


def spatial_clustering_index(
    positions: np.ndarray,
    types: np.ndarray,
    rng: np.random.Generator | None = None,
) -> dict[int, float]:
    """
    Compute a clustering index (Hopkins statistic variant) per type.
    Values > 0.5 indicate clustering; ~0.5 = random; < 0.5 = regular.
    """
    result: dict[int, float] = {}
    unique_types = np.unique(types)

    for t in unique_types:
        mask = types == t
        pts = positions[mask]
        n = len(pts)
        if n < 10:
            result[int(t)] = 0.5
            continue

        tree = KDTree(pts)

        # Sample m random points
        _rng = rng if rng is not None else np.random.default_rng()
        m = min(n // 2, 100)
        indices = _rng.choice(n, m, replace=False)
        sample = pts[indices]

        # NND from sample to rest
        dists_data, _ = tree.query(sample, k=2)
        nnd_data = dists_data[:, 1]

        # Random reference points
        lo = pts.min(axis=0)
        hi = pts.max(axis=0)
        random_pts = _rng.uniform(lo, hi, size=(m, pts.shape[1]))
        dists_rand, _ = tree.query(random_pts, k=1)
        nnd_rand = dists_rand.flatten()

        # Hopkins statistic
        sum_data = np.sum(nnd_data)
        sum_rand = np.sum(nnd_rand)
        denom = sum_data + sum_rand
        hopkins = sum_rand / denom if denom > 0 else 0.5

        result[int(t)] = float(hopkins)

    return result


def monochromatic_patch_score(
    positions: np.ndarray,
    types: np.ndarray,
    radius: float = 10e-6,
) -> float:
    """
    Measure degree of monochromatic patchiness.
    For each agent, check fraction of same-type neighbors within `radius`.
    Returns mean fraction (1.0 = perfectly monochromatic, 1/n_types = random).
    """
    if len(positions) == 0:
        return 0.0

    tree = KDTree(positions)
    fractions = []

    for i in range(len(positions)):
        neighbors = tree.query_ball_point(positions[i], radius)
        if len(neighbors) <= 1:
            continue
        same = sum(1 for j in neighbors if types[j] == types[i])
        fractions.append(same / len(neighbors))

    return float(np.mean(fractions)) if fractions else 0.0


def exclusion_radius(
    positions: np.ndarray,
    types: np.ndarray,
    target_type: int,
) -> float:
    """Mean distance from *target_type* agents to nearest different-type agent."""
    mask = types == target_type
    target_pts = positions[mask]
    other_pts = positions[~mask]

    if len(target_pts) == 0 or len(other_pts) == 0:
        return 0.0

    tree = KDTree(other_pts)
    dists, _ = tree.query(target_pts, k=1)
    return float(np.mean(dists))


def hopkins_statistic(
    positions: np.ndarray,
    n_samples: int | None = None,
    rng: np.random.Generator | None = None,
) -> float:
    """
    Hopkins clustering statistic over the full point cloud.

    H > 0.7 → significantly clustered; ≈ 0.5 → random.
    """
    n = len(positions)
    if n < 10:
        return 0.5

    m = n_samples if n_samples is not None else min(n // 2, 100)
    m = max(1, min(m, n))

    tree = KDTree(positions)

    _rng = rng if rng is not None else np.random.default_rng()
    indices = _rng.choice(n, m, replace=False)
    sample = positions[indices]
    dists_data, _ = tree.query(sample, k=2)
    nnd_data = dists_data[:, 1]

    lo = positions.min(axis=0)
    hi = positions.max(axis=0)
    random_pts = _rng.uniform(lo, hi, size=(m, positions.shape[1]))
    dists_rand, _ = tree.query(random_pts, k=1)
    nnd_rand = dists_rand.flatten()

    sum_data = np.sum(nnd_data)
    sum_rand = np.sum(nnd_rand)
    denom = sum_data + sum_rand
    return float(sum_rand / denom) if denom > 0 else 0.5


def comet_tail_asymmetry_index(
    positions: np.ndarray,
    concentrations: np.ndarray,
    flow_direction: int = 0,
) -> float:
    """
    Enhanced comet-tail metric measuring downstream elongation.

    *flow_direction* is the axis index (0 = x, 1 = y, 2 = z).
    Returns the ratio of concentration-weighted mean downstream distance
    to upstream distance; values > 1 indicate advective comet-tail.
    """
    if len(positions) == 0 or len(concentrations) == 0:
        return 1.0

    centroid = np.mean(positions, axis=0)
    projections = positions[:, flow_direction] - centroid[flow_direction]

    downstream = projections > 0
    upstream = ~downstream

    if not np.any(downstream) or not np.any(upstream):
        return 1.0

    c = np.abs(concentrations)
    w_down = np.sum(c[downstream] * np.abs(projections[downstream]))
    w_up = np.sum(c[upstream] * np.abs(projections[upstream]))

    return float(w_down / max(w_up, 1e-30))


def comet_tail_index(
    positions: np.ndarray,
    concentrations: np.ndarray,
    flow_direction: np.ndarray = np.array([1, 0, 0]),
) -> float:
    """
    Measure asymmetry of concentration field along flow direction.
    Returns ratio of downstream/upstream mean concentration.
    Values > 1 indicate comet-tail formation.
    """
    if len(positions) == 0 or len(concentrations) == 0:
        return 1.0

    centroid = np.mean(positions, axis=0)
    projections = np.dot(positions - centroid, flow_direction)

    downstream = concentrations[projections > 0]
    upstream = concentrations[projections <= 0]

    mean_down = np.mean(downstream) if len(downstream) > 0 else 0
    mean_up = np.mean(upstream) if len(upstream) > 0 else 1e-30

    return float(mean_down / max(mean_up, 1e-30))
