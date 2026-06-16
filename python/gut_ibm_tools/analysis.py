"""
Spatial clustering analysis and nearest-neighbor distance (NND) metrics.
Generates simulated HiPR-FISH spatial clustering metrics for validation.
"""

from __future__ import annotations

import numpy as np
from scipy.spatial import KDTree


def nearest_neighbor_distances(positions: np.ndarray, types: np.ndarray) -> dict[int, np.ndarray]:
    """
    Compute nearest-neighbor distances between agents of the same type.
    Returns dict mapping type → array of NNDs.
    """
    result: dict[int, np.ndarray] = {}
    unique_types = np.unique(types)

    for t in unique_types:
        mask = types == t
        pts = positions[mask]
        if len(pts) < 2:
            result[int(t)] = np.array([])
            continue
        tree = KDTree(pts)
        dists, _ = tree.query(pts, k=2)  # k=2: self + nearest
        result[int(t)] = dists[:, 1]  # column 1 = nearest non-self

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


def spatial_clustering_index(positions: np.ndarray, types: np.ndarray) -> dict[int, float]:
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
        m = min(n // 2, 100)
        indices = np.random.choice(n, m, replace=False)
        sample = pts[indices]

        # NND from sample to rest
        dists_data, _ = tree.query(sample, k=2)
        nnd_data = dists_data[:, 1]

        # Random reference points
        lo = pts.min(axis=0)
        hi = pts.max(axis=0)
        random_pts = np.random.uniform(lo, hi, size=(m, pts.shape[1]))
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
