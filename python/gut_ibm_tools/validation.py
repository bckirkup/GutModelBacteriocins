"""
Validation metrics against empirical datasets.
Tests the model's predictions against:
1. HiPR-FISH spatial clustering (mouse gut)
2. Longitudinal metagenomics (human, 70-80% retention)
"""

from __future__ import annotations

import numpy as np

from .analysis import (
    comet_tail_index,
    monochromatic_patch_score,
    nearest_neighbor_distances,
    spatial_clustering_index,
)
from .hdf5_reader import GutIBMData


def validate_spatial_signatures(data: GutIBMData, step: str) -> dict[str, float]:
    """
    Validate spatial signatures against HiPR-FISH expectations.

    Metrics:
    - monochromatic_score: should be > 0.7 for Enterobacteriaceae
    - clustering_index: should show significant clustering (> 0.6)
    - comet_tail_ratio: should be > 1.5 indicating advective asymmetry
    """
    agents = data.get_agents(step)
    grid = data.get_grid(step)

    positions = np.column_stack([agents["x"], agents["y"], agents["z"]])
    types = agents["type"]

    # Monochromatic patchiness (same-type neighbors)
    mono_score = monochromatic_patch_score(positions, types)

    # Spatial clustering
    clustering = spatial_clustering_index(positions, types)
    mean_clustering = float(np.mean(list(clustering.values()))) if clustering else 0.5

    # Comet-tail asymmetry
    bacteriocin = grid.get("bacteriocin", np.zeros(1))
    # Rough grid-position mapping (1D for now)
    n = len(bacteriocin)
    grid_pos = np.column_stack([
        np.linspace(0, 1e-3, n),
        np.zeros(n),
        np.zeros(n),
    ])
    comet = comet_tail_index(grid_pos, bacteriocin)

    return {
        "monochromatic_score": mono_score,
        "mean_clustering_index": mean_clustering,
        "comet_tail_ratio": comet,
    }


def validate_genomic_signatures(data: GutIBMData) -> dict[str, float]:
    """
    Validate longitudinal genomic predictions.

    Metrics:
    - resident_retention: should be 70-80%
    - dominant_bi_loci: residents should have more BI clusters
    - transient_receptor_downreg: transients should show lower expression
    """
    steps = data.steps
    if len(steps) < 2:
        return {"error": -1.0}

    first_step = steps[0]
    last_step = steps[-1]

    first_agents = data.get_agents(first_step)
    last_agents = data.get_agents(last_step)

    first_lineages = set(np.unique(first_agents.get("lineage", np.array([]))))
    last_lineages = set(np.unique(last_agents.get("lineage", np.array([]))))

    if not first_lineages:
        return {"resident_retention": 0.0}

    retained = first_lineages & last_lineages
    retention = len(retained) / len(first_lineages)

    # BI loci complexity of residents vs transients
    last_lin = data.get_lineage(last_step)
    n_bi = last_lin.get("num_bi_loci", np.array([]))

    if len(n_bi) > 0:
        last_lineage_ids = last_agents.get("lineage", np.array([]))
        resident_mask = np.isin(last_lineage_ids, list(retained))
        transient_mask = ~resident_mask

        resident_bi = float(np.mean(n_bi[resident_mask])) if np.any(resident_mask) else 0
        transient_bi = float(np.mean(n_bi[transient_mask])) if np.any(transient_mask) else 0
    else:
        resident_bi = 0.0
        transient_bi = 0.0

    # Receptor expression of transients
    btuB_expr = last_lin.get("btuB_expression", np.array([]))
    if len(btuB_expr) > 0:
        last_lineage_ids = last_agents.get("lineage", np.array([]))
        transient_mask = ~np.isin(last_lineage_ids, list(retained))
        transient_expr = float(np.mean(btuB_expr[transient_mask])) if np.any(transient_mask) else 1.0
    else:
        transient_expr = 1.0

    return {
        "resident_retention": retention,
        "resident_mean_bi_loci": resident_bi,
        "transient_mean_bi_loci": transient_bi,
        "transient_mean_btuB_expression": transient_expr,
    }
