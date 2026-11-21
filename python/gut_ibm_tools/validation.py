"""
Validation metrics against empirical datasets.
Tests the model's predictions against:
1. Exclusion-radius / NND spatial clustering (VADI §75)
2. Longitudinal metagenomics (human, 70-80% retention)
"""

from __future__ import annotations

import numpy as np

from .analysis import (
    comet_tail_asymmetry_index,
    comet_tail_index,
    exclusion_radius,
    hopkins_statistic,
    monochromatic_patch_score,
    nearest_neighbor_distances,
)
from .hdf5_reader import GutIBMData


def _toxin_source_reference(
    agents: dict[str, np.ndarray],
    grid_positions: np.ndarray,
) -> np.ndarray:
    """Choose a biomass-weighted producer reference with documented fallbacks.

    Producer agents are identified by ``n_bi_loci > 0``.  If that field or a
    producer population is unavailable, all agents are used.  With no agents,
    the domain midpoint (the mean of the uniform cell centers) is used.
    """
    if len(agents.get("x", [])) == 0:
        return np.mean(grid_positions, axis=0)

    agent_positions = np.column_stack([
        agents["x"], agents["y"], agents["z"],
    ])
    biomass = np.asarray(agents.get("biomass", np.ones(len(agent_positions))))
    producer_mask = np.asarray(
        agents.get("n_bi_loci", np.zeros(len(agent_positions))) > 0,
        dtype=bool,
    )
    selected = producer_mask if np.any(producer_mask) else np.ones(
        len(agent_positions), dtype=bool
    )
    weights = biomass[selected]
    if len(weights) == 0 or not np.all(np.isfinite(weights)) or np.sum(weights) <= 0:
        weights = np.ones(np.count_nonzero(selected))
    return np.average(agent_positions[selected], axis=0, weights=weights)


def validate_spatial_signatures(
    data: GutIBMData,
    step: str,
    rng: np.random.Generator | None = None,
    *,
    depth_window_m: tuple[float, float] | None = None,
    depth_axis: int = 2,
) -> dict[str, float]:
    """
    Validate spatial signatures via exclusion-radius clustering (VADI §75).

    Retained metrics:
        monochromatic_score  – should be > 0.7
        comet_tail_ratio     – should be > 1.5
        The comet-tail flow axis is positive x, the model's distal-flow
        convention from ``AdvectionField::velocity``; it is not inferred from
        HDF5 storage order.  Its reference is the biomass-weighted centroid of
        producer agents, then all agents, then the domain midpoint when no
        agents exist.  The x period is the domain extent because x is
        periodic.  These reference and wrapping choices are part of the metric
        definition.
    New exclusion-radius metrics:
        mean_exclusion_radius – mean over all types
        hopkins_statistic     – H > 0.7 → significantly clustered
        nnd_mean              – grand mean of inter-type NND
        comet_tail_asymmetry  – concentration-weighted downstream elongation
    """
    agents = data.get_agents(step)
    grid = data.get_grid(step)
    expected_length = int(np.prod(data.grid_shape()))

    positions = np.column_stack([agents["x"], agents["y"], agents["z"]])
    types = agents["type"]
    if depth_window_m is not None:
        lo_m, hi_m = depth_window_m
        if lo_m > hi_m:
            raise ValueError("depth_window_m lower bound must be <= upper bound")
        if depth_axis not in (0, 1, 2):
            raise ValueError("depth_axis must be 0, 1, or 2")
        slab_mask = (positions[:, depth_axis] >= lo_m) & (positions[:, depth_axis] <= hi_m)
        positions = positions[slab_mask]
        types = types[slab_mask]

    # --- retained metrics ------------------------------------------------
    mono_score = monochromatic_patch_score(positions, types)

    bacteriocin_keys = [
        k for k in grid
        if k.startswith("bacteriocin_") and k != "bacteriocin_lumped"
    ]
    if bacteriocin_keys:
        toxin_fields = []
        for key in bacteriocin_keys:
            field = np.asarray(grid[key]).ravel()
            if len(field) != expected_length:
                raise ValueError(
                    f"flattened toxin field length {len(field)} for {key} "
                    f"does not match grid size {expected_length}"
                )
            toxin_fields.append(field)
        bacteriocin = np.sum(toxin_fields, axis=0)
    else:
        bacteriocin = grid.get(
            "bacteriocin_BtuB",
            grid.get("bacteriocin_lumped", np.zeros(1)),
        )
        bacteriocin = np.asarray(bacteriocin).ravel()
    if len(bacteriocin) != expected_length:
        raise ValueError(
            "flattened toxin field length "
            f"{len(bacteriocin)} does not match grid size {expected_length}"
        )
    grid_pos = data.grid_cell_centers()
    if depth_window_m is not None:
        lo_m, hi_m = depth_window_m
        slab_agents = {
            k: np.asarray(v)[slab_mask]
            for k, v in agents.items()
            if len(v) == len(slab_mask)
        }
        grid_slab = (grid_pos[:, depth_axis] >= lo_m) & (grid_pos[:, depth_axis] <= hi_m)
        if np.any(grid_slab):
            grid_pos = grid_pos[grid_slab]
            bacteriocin = bacteriocin[grid_slab]
        reference = _toxin_source_reference(slab_agents, grid_pos)
    else:
        reference = _toxin_source_reference(agents, grid_pos)
    period = data.grid_periods()[0]
    comet = comet_tail_index(
        grid_pos,
        bacteriocin,
        flow_direction=np.array([1.0, 0.0, 0.0]),
        reference=reference,
        period=period,
    )

    # --- new exclusion-radius metrics ------------------------------------
    unique_types = np.unique(types)

    # Localized exclusion radii per phylogroup
    excl_radii = {
        int(t): exclusion_radius(positions, types, int(t))
        for t in unique_types
    }
    mean_excl = float(np.mean(list(excl_radii.values()))) if excl_radii else 0.0

    # Hopkins clustering over the full point cloud
    hopkins = hopkins_statistic(positions, rng=rng)

    # NND between competing clones
    nnd = nearest_neighbor_distances(positions, types)
    all_nnd = np.concatenate(list(nnd.values())) if nnd else np.array([])
    nnd_mean = float(np.mean(all_nnd)) if len(all_nnd) > 0 else 0.0

    # Enhanced comet-tail asymmetry (concentration-weighted)
    comet_asym = comet_tail_asymmetry_index(
        grid_pos,
        bacteriocin,
        flow_direction=0,
        reference=reference,
        period=period,
    )

    return {
        "monochromatic_score": mono_score,
        "comet_tail_ratio": comet,
        "mean_exclusion_radius": mean_excl,
        "hopkins_statistic": hopkins,
        "nnd_mean": nnd_mean,
        "comet_tail_asymmetry": comet_asym,
    }


def compute_lineage_retention_metrics(
    first_lineage_ids: np.ndarray,
    last_lineage_ids: np.ndarray,
    *,
    lod_relative_abundance: float = 0.0,
) -> dict[str, float]:
    """
    Compute unweighted founder-lineage persistence and abundance-weighted retention.

    Parameters
    ----------
    first_lineage_ids, last_lineage_ids:
        Per-agent founder/lineage integer identifiers at initial and final steps.
    lod_relative_abundance:
        Minimum final relative abundance in ``[0, 1)`` required for a lineage to
        count as detected at the final endpoint (models metagenomic/culture LOD).

    Returns
    -------
    dict with keys:
        - ``resident_retention``: unweighted fraction of initial lineages detected at end
        - ``founder_lineage_persistence``: explicit alias of ``resident_retention``
        - ``abundance_weighted_retention``: initial-abundance-weighted mass of surviving
          founder lineages detected above *lod_relative_abundance*
        - ``final_resident_abundance_fraction``: share of final population belonging to
          detected initial founder lineages
    """
    if not 0.0 <= lod_relative_abundance < 1.0:
        raise ValueError("lod_relative_abundance must be in [0, 1)")

    first_arr = np.asarray(first_lineage_ids).ravel()
    last_arr = np.asarray(last_lineage_ids).ravel()
    if len(first_arr) == 0:
        return {
            "resident_retention": 0.0,
            "founder_lineage_persistence": 0.0,
            "abundance_weighted_retention": 0.0,
            "final_resident_abundance_fraction": 0.0,
        }

    init_ids, init_counts = np.unique(first_arr, return_counts=True)
    init_weights = {
        int(lid): float(cnt) / float(len(first_arr))
        for lid, cnt in zip(init_ids, init_counts)
    }

    if len(last_arr) == 0:
        return {
            "resident_retention": 0.0,
            "founder_lineage_persistence": 0.0,
            "abundance_weighted_retention": 0.0,
            "final_resident_abundance_fraction": 0.0,
        }

    final_ids, final_counts = np.unique(last_arr, return_counts=True)
    n_final = float(len(last_arr))
    detected_final = {
        int(lid): float(cnt) / n_final
        for lid, cnt in zip(final_ids, final_counts)
        if (float(cnt) / n_final) >= lod_relative_abundance and cnt > 0
    }

    retained = set(init_weights) & set(detected_final)
    unweighted = float(len(retained)) / float(len(init_weights))
    init_weighted = float(sum(init_weights[lid] for lid in retained))
    final_share = float(sum(detected_final[lid] for lid in retained))

    return {
        "resident_retention": unweighted,
        "founder_lineage_persistence": unweighted,
        "abundance_weighted_retention": init_weighted,
        "final_resident_abundance_fraction": final_share,
    }


def validate_genomic_signatures(
    data: GutIBMData,
    *,
    lod_relative_abundance: float = 0.0,
) -> dict[str, float]:
    """
    Validate longitudinal genomic predictions.

    Metrics:
    - resident_retention / founder_lineage_persistence: unweighted founder survival
    - abundance_weighted_retention: initial-share-weighted founder persistence
    - final_resident_abundance_fraction: final population share of surviving founders
    - resident_mean_bi_loci: residents should have more BI clusters
    - transient_mean_btuB_expression: transients should show lower expression
    """
    steps = data.steps
    if len(steps) < 2:
        return {"error": -1.0}

    first_step = steps[0]
    last_step = steps[-1]

    first_agents = data.get_agents(first_step)
    last_agents = data.get_agents(last_step)

    first_arr = np.asarray(first_agents.get("lineage", np.array([])))
    last_arr = np.asarray(last_agents.get("lineage", np.array([])))

    if len(first_arr) == 0:
        return {"resident_retention": 0.0}

    retention_metrics = compute_lineage_retention_metrics(
        first_arr,
        last_arr,
        lod_relative_abundance=lod_relative_abundance,
    )
    first_lineages = set(np.unique(first_arr))
    if len(last_arr) > 0:
        final_ids, final_counts = np.unique(last_arr, return_counts=True)
        n_final = float(len(last_arr))
        detected_final_ids = {
            int(lid)
            for lid, cnt in zip(final_ids, final_counts)
            if (float(cnt) / n_final) >= lod_relative_abundance and cnt > 0
        }
    else:
        detected_final_ids = set()
    retained = first_lineages & detected_final_ids

    # BI loci complexity of residents vs transients
    last_lin = data.get_lineage(last_step)
    n_bi = last_lin.get("num_bi_loci", np.array([]))

    if len(n_bi) > 0:
        resident_mask = np.isin(last_arr, list(retained))
        transient_mask = ~resident_mask

        resident_bi = float(np.mean(n_bi[resident_mask])) if np.any(resident_mask) else 0.0
        transient_bi = float(np.mean(n_bi[transient_mask])) if np.any(transient_mask) else 0.0
    else:
        resident_bi = 0.0
        transient_bi = 0.0

    # Receptor expression of transients
    btub_expr = last_lin.get("btuB_expression", np.array([]))
    if len(btub_expr) > 0:
        transient_mask = ~np.isin(last_arr, list(retained))
        transient_expr = float(np.mean(btub_expr[transient_mask])) if np.any(transient_mask) else 1.0
    else:
        transient_expr = 1.0

    return {
        **retention_metrics,
        "resident_mean_bi_loci": resident_bi,
        "transient_mean_bi_loci": transient_bi,
        "transient_mean_btuB_expression": transient_expr,
    }
