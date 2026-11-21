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
        _rng = rng if rng is not None else np.random.default_rng(42)
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

    _rng = rng if rng is not None else np.random.default_rng(42)
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
    reference: np.ndarray | None = None,
    period: float | None = None,
) -> float:
    """
    Enhanced comet-tail metric measuring downstream elongation.

    *flow_direction* is the axis index (0 = x, 1 = y, 2 = z).
    *reference* is the point from which displacement is measured; when
    omitted, the position centroid is used.  When *period* is supplied,
    along-flow displacement uses the minimum-image convention.
    Returns the ratio of concentration-weighted mean downstream distance
    to upstream distance; values > 1 indicate advective comet-tail.
    """
    if len(positions) == 0 or len(concentrations) == 0:
        return 1.0

    if reference is None:
        reference = np.mean(positions, axis=0)
    projections = positions[:, flow_direction] - reference[flow_direction]
    if period is not None:
        projections -= period * np.round(projections / period)

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
    flow_direction: np.ndarray | None = None,
    reference: np.ndarray | None = None,
    period: float | None = None,
) -> float:
    """
    Measure asymmetry of concentration field along flow direction.
    Returns ratio of downstream/upstream mean concentration.
    Values > 1 indicate comet-tail formation.  *reference* overrides the
    default position centroid.  When *period* is supplied, the projected
    along-flow displacement uses the minimum-image convention.
    """
    if len(positions) == 0 or len(concentrations) == 0:
        return 1.0

    if flow_direction is None:
        flow_direction = np.array([1, 0, 0])

    if reference is None:
        reference = np.mean(positions, axis=0)
    projections = np.dot(positions - reference, flow_direction)
    if period is not None:
        projections -= period * np.round(projections / period)

    downstream = concentrations[projections > 0]
    upstream = concentrations[projections <= 0]

    mean_down = np.mean(downstream) if len(downstream) > 0 else 0
    mean_up = np.mean(upstream) if len(upstream) > 0 else 1e-30

    return float(mean_down / max(mean_up, 1e-30))


def _window_slice_arrays(
    times_s: np.ndarray,
    counts_focal: np.ndarray,
    counts_competitor: np.ndarray,
    *,
    pseudocount: float,
    window_s: float | None,
    t_end_s: float | None,
) -> tuple[np.ndarray, np.ndarray]:
    """Extract sorted time and log10((n1 + eps) / (n2 + eps)) arrays on [t_end - window_s, t_end]."""
    if pseudocount <= 0.0:
        raise ValueError("pseudocount must be positive")
    t = np.asarray(times_s, dtype=np.float64).ravel()
    n1 = np.asarray(counts_focal, dtype=np.float64).ravel()
    n2 = np.asarray(counts_competitor, dtype=np.float64).ravel()
    if not (len(t) == len(n1) == len(n2)):
        raise ValueError("times_s, counts_focal, and counts_competitor must have equal length")
    if len(t) == 0:
        return np.array([], dtype=np.float64), np.array([], dtype=np.float64)
    if np.any(n1 < 0.0) or np.any(n2 < 0.0):
        raise ValueError("strain counts must be non-negative")

    order = np.argsort(t, kind="mergesort")
    t = t[order]
    n1 = n1[order]
    n2 = n2[order]

    t_hi = float(t[-1]) if t_end_s is None else float(t_end_s)
    mask = t <= t_hi + 1e-12
    if window_s is not None:
        if window_s <= 0.0:
            raise ValueError("window_s must be positive when provided")
        mask &= t >= (t_hi - float(window_s) - 1e-12)

    t_sel = t[mask]
    ratios = np.log10((n1[mask] + pseudocount) / (n2[mask] + pseudocount))
    return t_sel, ratios


def log_ratio_ols_slope(
    times_s: np.ndarray | list[float],
    counts_focal: np.ndarray | list[int] | list[float],
    counts_competitor: np.ndarray | list[int] | list[float],
    *,
    pseudocount: float = 0.5,
    window_s: float | None = None,
    t_end_s: float | None = None,
    doubling_time_s: float = 2700.0,
    wave_warning_threshold_decades: float = 0.10,
    min_recommended_cycles: float = 2.0,
) -> dict[str, float | int | bool | None]:
    """
    Division-wave-robust ordinary-least-squares slope of log10((n1 + eps)/(n2 + eps)).

    Uniformly initialized cohorts divide in partially synchronized waves whose
    single-step endpoint ratio can oscillate by +-0.2 decades (overstating
    competitive advantage by >2x). Fitting over multiple division cycles and
    inspecting detrended peak-to-trough residuals guards against phase aliasing.
    """
    if doubling_time_s <= 0.0:
        raise ValueError("doubling_time_s must be positive")
    t_sel, y_sel = _window_slice_arrays(
        np.asarray(times_s),
        np.asarray(counts_focal),
        np.asarray(counts_competitor),
        pseudocount=pseudocount,
        window_s=window_s,
        t_end_s=t_end_s,
    )
    n_samples = len(t_sel)
    if n_samples < 2:
        return {
            "slope_decades_per_h": None,
            "slope_se_decades_per_h": None,
            "intercept_log10_ratio": None,
            "n_samples": n_samples,
            "span_s": 0.0,
            "estimated_cycles": 0.0,
            "residual_peak_to_trough_decades": None,
            "wave_aliasing_warning": True,
        }

    span_s = float(t_sel[-1] - t_sel[0])
    if span_s <= 1e-12:
        return {
            "slope_decades_per_h": None,
            "slope_se_decades_per_h": None,
            "intercept_log10_ratio": None,
            "n_samples": n_samples,
            "span_s": span_s,
            "estimated_cycles": 0.0,
            "residual_peak_to_trough_decades": 0.0,
            "wave_aliasing_warning": True,
        }

    t_centered = t_sel - np.mean(t_sel)
    y_mean = float(np.mean(y_sel))
    s_xx = float(np.sum(t_centered * t_centered))
    slope_per_s = float(np.sum(t_centered * (y_sel - y_mean)) / s_xx)
    intercept = float(y_mean - slope_per_s * float(np.mean(t_sel)))
    residuals = y_sel - (intercept + slope_per_s * t_sel)

    if n_samples > 2:
        s2 = float(np.sum(residuals * residuals) / (n_samples - 2))
        slope_se_per_h = float(np.sqrt(max(s2, 0.0) / s_xx) * 3600.0)
    else:
        slope_se_per_h = None

    resid_ptp = float(np.max(residuals) - np.min(residuals))
    est_cycles = span_s / float(doubling_time_s)
    aliasing_warn = bool(
        (est_cycles < min_recommended_cycles)
        or (resid_ptp > wave_warning_threshold_decades)
    )

    return {
        "slope_decades_per_h": slope_per_s * 3600.0,
        "slope_se_decades_per_h": slope_se_per_h,
        "intercept_log10_ratio": intercept,
        "n_samples": n_samples,
        "span_s": span_s,
        "estimated_cycles": est_cycles,
        "residual_peak_to_trough_decades": resid_ptp,
        "wave_aliasing_warning": aliasing_warn,
    }


def wave_averaged_log_ratio(
    times_s: np.ndarray | list[float],
    counts_focal: np.ndarray | list[int] | list[float],
    counts_competitor: np.ndarray | list[int] | list[float],
    *,
    pseudocount: float = 0.5,
    tail_window_s: float = 7200.0,
    t_end_s: float | None = None,
    wave_warning_threshold_decades: float = 0.10,
) -> dict[str, float | int | bool | None]:
    """
    Compute trapezoid time-averaged and median log10 strain ratios over a trailing window.

    Reports ``endpoint_bias_decades = endpoint_log10_ratio - trapezoid_mean_log10_ratio``
    to quantify single-timepoint synchronized division-wave aliasing directly.
    """
    t_sel, y_sel = _window_slice_arrays(
        np.asarray(times_s),
        np.asarray(counts_focal),
        np.asarray(counts_competitor),
        pseudocount=pseudocount,
        window_s=tail_window_s,
        t_end_s=t_end_s,
    )
    n_samples = len(t_sel)
    if n_samples == 0:
        return {
            "trapezoid_mean_log10_ratio": None,
            "tail_median_log10_ratio": None,
            "endpoint_log10_ratio": None,
            "endpoint_bias_decades": None,
            "peak_to_trough_decades": None,
            "n_samples": 0,
            "span_s": 0.0,
            "wave_aliasing_warning": True,
        }

    endpoint = float(y_sel[-1])
    median_val = float(np.median(y_sel))
    ptp = float(np.max(y_sel) - np.min(y_sel))
    span_s = float(t_sel[-1] - t_sel[0])

    if span_s > 1e-12:
        trapz_fn = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
        trap_mean = float(trapz_fn(y_sel, t_sel) / span_s)
    else:
        trap_mean = float(np.mean(y_sel))

    bias = endpoint - trap_mean
    return {
        "trapezoid_mean_log10_ratio": trap_mean,
        "tail_median_log10_ratio": median_val,
        "endpoint_log10_ratio": endpoint,
        "endpoint_bias_decades": bias,
        "peak_to_trough_decades": ptp,
        "n_samples": n_samples,
        "span_s": span_s,
        "wave_aliasing_warning": bool(ptp > wave_warning_threshold_decades or span_s <= 1e-12),
    }
