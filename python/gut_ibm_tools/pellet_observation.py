"""
Fecal-pellet colony/read multinomial censoring and multi-host cage pooling.

Models culture-plate and amplicon/metagenomic detection limits for mouse and
human stool competitions, including Kirkup & Riley (2004, Nature 428:412-414),
where 3 co-housed mice per cage were sampled by plating ~36 colonies per fecal
homogenate (~1%-2.8% strain-fraction limit of detection).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import binom

DEFAULT_PELLET_RNG_SEED = 2004


@dataclass(frozen=True)
class PelletAssayConfig:
    """Observation model parameters for fecal-pellet plating or read sampling."""

    n_draws: int = 36
    lod_fraction: float = 0.01
    hosts_per_cage: int = 3

    def __post_init__(self) -> None:
        if self.n_draws < 1:
            raise ValueError("n_draws must be >= 1")
        if not 0.0 <= self.lod_fraction < 1.0:
            raise ValueError("lod_fraction must be in [0, 1)")
        if self.hosts_per_cage < 1:
            raise ValueError("hosts_per_cage must be >= 1")


def analytical_detection_probability(
    true_fractions: np.ndarray | list[float],
    *,
    n_draws: int = 36,
    min_positive_draws: int = 1,
) -> np.ndarray:
    """
    Exact binomial probability P(X >= min_positive_draws) in *n_draws* independent picks.

    At ``min_positive_draws == 1``, equals ``1 - (1 - p) ** n_draws`` (e.g. at
    ``p = 0.01`` and ``n_draws = 36``, detection probability is only ``30.4%``,
    so a stable 1% mucosal minority reads as undetected in 69.6% of pellets).
    """
    if n_draws < 1:
        raise ValueError("n_draws must be >= 1")
    if min_positive_draws < 1 or min_positive_draws > n_draws:
        raise ValueError("min_positive_draws must be in [1, n_draws]")

    p = np.asarray(true_fractions, dtype=np.float64)
    if np.any(p < -1e-12) or np.any(p > 1.0 + 1e-12):
        raise ValueError("true_fractions entries must lie in [0, 1]")
    p = np.clip(p, 0.0, 1.0)

    if min_positive_draws == 1:
        return 1.0 - np.power(1.0 - p, n_draws)

    # Exact complement sum_{k=0}^{m-1} C(n, k) p^k (1-p)^(n-k) in log-safe form
    return np.asarray(1.0 - binom.cdf(min_positive_draws - 1, n_draws, p), dtype=np.float64)


def censor_pellet_sample(
    strain_counts: np.ndarray | list[int] | list[float],
    config: PelletAssayConfig | None = None,
    *,
    rng: np.random.Generator | None = None,
) -> dict[str, np.ndarray | int | float]:
    """
    Simulate multinomial colony/read sampling and LOD censoring on a single fecal sample.

    Raw simulation counts are never mutated; returns true fractions alongside
    sampled draw counts, censored observed fractions, and false-extinction flags.
    """
    if config is None:
        config = PelletAssayConfig()
    if rng is None:
        rng = np.random.default_rng(DEFAULT_PELLET_RNG_SEED)

    counts = np.asarray(strain_counts, dtype=np.float64).ravel()
    if np.any(counts < 0.0):
        raise ValueError("strain_counts must be non-negative")

    total = float(np.sum(counts))
    n_strains = len(counts)
    if total <= 0.0 or n_strains == 0:
        zeros = np.zeros(n_strains, dtype=np.float64)
        zeros_int = np.zeros(n_strains, dtype=np.int64)
        zeros_bool = np.zeros(n_strains, dtype=bool)
        return {
            "true_fractions": zeros,
            "sampled_counts": zeros_int,
            "sampled_fractions": zeros,
            "detected_mask": zeros_bool,
            "censored_fractions": zeros,
            "false_extinction_mask": zeros_bool,
            "n_draws": config.n_draws,
            "lod_fraction": config.lod_fraction,
        }

    true_fracs = counts / total
    sampled_counts = rng.multinomial(config.n_draws, true_fracs).astype(np.int64)
    sampled_fracs = sampled_counts.astype(np.float64) / float(config.n_draws)
    detected = (sampled_counts >= 1) & (sampled_fracs >= config.lod_fraction)

    censored = np.where(detected, sampled_fracs, 0.0)
    censored_sum = float(np.sum(censored))
    if censored_sum > 0.0:
        censored = censored / censored_sum

    true_present = counts > 0.0
    false_extinction = true_present & (~detected)

    return {
        "true_fractions": true_fracs,
        "sampled_counts": sampled_counts,
        "sampled_fractions": sampled_fracs,
        "detected_mask": detected,
        "censored_fractions": censored,
        "false_extinction_mask": false_extinction,
        "n_draws": config.n_draws,
        "lod_fraction": config.lod_fraction,
    }


def pool_cage_samples(
    host_strain_counts: np.ndarray | list[list[float]],
    config: PelletAssayConfig | None = None,
    *,
    host_weights: np.ndarray | list[float] | None = None,
    rng: np.random.Generator | None = None,
) -> dict[str, np.ndarray | int | float]:
    """
    Pool fecal shedding across co-housed hosts in a cage before multinomial plating.

    Matches cage-level homogenate protocols (e.g. 3 mice/cage in Kirkup & Riley
    2004). *host_strain_counts* has shape ``(n_hosts, n_strains)``.
    """
    if config is None:
        config = PelletAssayConfig()
    matrix = np.asarray(host_strain_counts, dtype=np.float64)
    if matrix.ndim != 2:
        raise ValueError("host_strain_counts must be a 2D array of shape (n_hosts, n_strains)")
    if np.any(matrix < 0.0):
        raise ValueError("host_strain_counts must be non-negative")

    n_hosts = matrix.shape[0]
    if n_hosts == 0:
        raise ValueError("host_strain_counts must contain at least one host")

    if host_weights is None:
        weights = np.ones(n_hosts, dtype=np.float64)
    else:
        weights = np.asarray(host_weights, dtype=np.float64).ravel()
        if len(weights) != n_hosts or np.any(weights < 0.0) or float(np.sum(weights)) <= 0.0:
            raise ValueError("host_weights must match n_hosts with positive sum")

    pooled_abundance = np.sum(matrix * weights[:, None], axis=0)
    per_host_present = np.sum(matrix > 0.0, axis=0).astype(np.int64)
    sample_res = censor_pellet_sample(pooled_abundance, config=config, rng=rng)

    return {
        **sample_res,
        "n_hosts_pooled": n_hosts,
        "hosts_carrying_strain": per_host_present,
    }
