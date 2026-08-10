"""Colony-scale spatial observables for GutIBM agent snapshots."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from .hdf5_reader import GutIBMData
from .path_utils import prepare_output_file, validate_input_path

PRODUCER_THRESHOLDS: tuple[int, ...] = (113, 527, 1361)
AGENT_COLUMNS = ("agent_id", "colony_id", "type", "n_bi_loci", "x", "y", "z")


@dataclass(frozen=True)
class ColonyConfig:
    """DBSCAN and diagnostic settings."""

    eps: float | None = None
    min_samples: int = 3
    eps_rule: str = "median_nnd"
    eps_factor: float = 2.0
    knee_k: int = 4
    sensitivity_factors: tuple[float, ...] = (0.5, 1.0, 2.0, 4.0)
    producer_thresholds: tuple[int, ...] = PRODUCER_THRESHOLDS


@dataclass
class ColonyCatalog:
    """Joinable per-agent and per-colony colony observables."""

    agents: dict[str, np.ndarray]
    colonies: dict[str, np.ndarray]
    eps: float
    eps_diagnostics: dict[str, float | str]
    eps_sensitivity: dict[str, np.ndarray]

    def write_csv(self, agent_path: str | Path, colony_path: str | Path) -> None:
        """Write both tidy tables after validating output paths."""
        _write_table(prepare_output_file(agent_path), self.agents)
        _write_table(prepare_output_file(colony_path), self.colonies)


def estimate_eps(
    positions: np.ndarray,
    *,
    rule: str = "median_nnd",
    factor: float = 2.0,
    k: int = 4,
) -> tuple[float, dict[str, float | str]]:
    """Estimate DBSCAN ``eps`` using median NND or a k-distance knee.

    ``median_nnd`` returns ``factor * median(nearest-neighbour distance)``.
    ``knee`` uses the maximum perpendicular distance from the line joining
    the endpoints of the sorted k-th-neighbour curve.  Diagnostics preserve
    the rule and measured distances for auditability.
    """
    points = _positions(positions)
    if len(points) < 2:
        raise ValueError("at least two positions are required")
    distances, _ = cKDTree(points).query(points, k=2)
    nnd = distances[:, 1]
    if rule == "median_nnd":
        median = float(np.median(nnd))
        return factor * median, {
            "rule": rule,
            "factor": factor,
            "median_nnd": median,
        }
    if rule != "knee":
        raise ValueError("rule must be 'median_nnd' or 'knee'")
    kth = cKDTree(points).query(points, k=min(k + 1, len(points)))[0][:, -1]
    curve = np.sort(kth)
    x = np.linspace(0.0, 1.0, len(curve))
    y = (curve - curve[0]) / (curve[-1] - curve[0]) if curve[-1] > curve[0] else np.zeros_like(curve)
    line = y[0] + x * (y[-1] - y[0])
    index = int(np.argmax(np.abs(y - line)))
    return float(curve[index]), {"rule": rule, "k": float(k), "knee_index": float(index)}


def dbscan_colonies(
    positions: np.ndarray,
    *,
    eps: float,
    min_samples: int,
) -> np.ndarray:
    """Return faithful DBSCAN labels (-1 noise) using a sparse core graph."""
    points = _positions(positions)
    if eps <= 0 or min_samples < 1:
        raise ValueError("eps must be positive and min_samples must be positive")
    tree = cKDTree(points)
    pairs = np.asarray(list(tree.query_pairs(eps)), dtype=np.int64)
    adjacency = _symmetric_graph(pairs, len(points))
    degree = np.asarray(adjacency.sum(axis=1)).ravel() + 1
    core = degree >= min_samples
    labels = np.full(len(points), -1, dtype=np.int64)
    if np.any(core):
        core_indices = np.flatnonzero(core)
        core_graph = adjacency.tocsr()[core_indices][:, core_indices]
        count, component = connected_components(core_graph, directed=False)
        labels[core_indices] = component
        for index in np.flatnonzero(~core):
            neighbours = tree.query_ball_point(points[index], eps)
            core_neighbours = [item for item in neighbours if core[item]]
            if core_neighbours:
                labels[index] = int(component[np.flatnonzero(core_indices == min(core_neighbours))[0]])
        if count == 0:
            labels.fill(-1)
    return labels


def build_colony_catalog(
    agents: dict[str, np.ndarray],
    config: ColonyConfig | None = None,
) -> ColonyCatalog:
    """Build per-agent and per-colony tables from ``GutIBMData.get_agents``."""
    cfg = config or ColonyConfig()
    positions = np.column_stack([_array(agents, name, len(agents["id"])) for name in ("x", "y", "z")])
    eps, diagnostics = (
        (float(cfg.eps), {"rule": "configured", "eps": float(cfg.eps)})
        if cfg.eps is not None
        else estimate_eps(positions, rule=cfg.eps_rule, factor=cfg.eps_factor, k=cfg.knee_k)
    )
    labels = dbscan_colonies(positions, eps=eps, min_samples=cfg.min_samples)
    diagnostics["noise_fraction"] = float(np.mean(labels < 0))
    sensitivity_eps = eps * np.asarray(cfg.sensitivity_factors, dtype=float)
    sensitivity_labels = [dbscan_colonies(positions, eps=value, min_samples=cfg.min_samples) for value in sensitivity_eps]
    sensitivity_counts = np.array([np.sum(value >= 0) for value in sensitivity_labels])
    sensitivity_colonies = np.array([len(np.unique(value[value >= 0])) for value in sensitivity_labels])
    sensitivity_noise = 1.0 - sensitivity_counts / max(len(positions), 1)
    agent_table = {
        "agent_id": np.asarray(agents["id"], dtype=np.int64).copy(),
        "colony_id": labels,
        "type": _array(agents, "type", len(positions), dtype=np.int64),
        "n_bi_loci": _array(agents, "n_bi_loci", len(positions), dtype=np.int64),
        "x": positions[:, 0],
        "y": positions[:, 1],
        "z": positions[:, 2],
    }
    colony_table = _colony_table(agent_table, labels, agents, cfg.producer_thresholds)
    return ColonyCatalog(
        agent_table,
        colony_table,
        eps,
        diagnostics,
        {"factor": np.asarray(cfg.sensitivity_factors), "eps": sensitivity_eps,
         "n_colonies": sensitivity_colonies, "noise_fraction": sensitivity_noise},
    )


def colony_catalog_from_hdf5(data: GutIBMData, step: str, config: ColonyConfig | None = None) -> ColonyCatalog:
    """Build a catalog from one open HDF5 agent step."""
    return build_colony_catalog(data.get_agents(step), config)


def producer_threshold_flags(
    n_producers: int,
    thresholds: Iterable[int] = PRODUCER_THRESHOLDS,
) -> dict[str, bool]:
    """Return whether a colony reaches each named producer threshold."""
    return {f"reaches_{threshold}": n_producers >= threshold for threshold in thresholds}


def _colony_table(
    table: dict[str, np.ndarray],
    labels: np.ndarray,
    agents: dict[str, np.ndarray],
    thresholds: tuple[int, ...],
) -> dict[str, np.ndarray]:
    """Compute tidy statistics for non-noise DBSCAN components."""
    lineage = agents.get("lineage_id")
    ids = np.unique(labels[labels >= 0])
    rows: list[dict[str, float | int | bool]] = []
    for colony_id in ids:
        member = labels == colony_id
        xyz = np.column_stack([table[axis][member] for axis in ("x", "y", "z")])
        centroid = xyz.mean(axis=0)
        radius = np.sqrt(np.mean(np.sum((xyz - centroid) ** 2, axis=1)))
        distances = np.sqrt(np.sum((xyz - centroid) ** 2, axis=1))
        genotypes = _genotypes(table["type"][member], table["n_bi_loci"][member], lineage[member] if lineage is not None else None)
        producers = int(np.sum(table["n_bi_loci"][member] > 0))
        row: dict[str, float | int | bool] = {
            "colony_id": int(colony_id), "n_members": int(np.sum(member)),
            "centroid_x": float(centroid[0]), "centroid_y": float(centroid[1]), "centroid_z": float(centroid[2]),
            "radius_of_gyration": float(radius), "max_member_distance_from_centroid": float(np.max(distances)),
            "n_genotypes": len(set(genotypes)), "dominant_genotype_fraction": max(Counter(genotypes).values()) / len(genotypes),
            "is_clonal": len(set(genotypes)) == 1, "n_producers": producers,
            "producer_fraction": producers / int(np.sum(member)), "mean_z": float(np.mean(table["z"][member])),
        }
        row.update(producer_threshold_flags(producers, thresholds))
        rows.append(row)
    result = {key: np.asarray([row[key] for row in rows]) for key in rows[0]} if rows else {}
    if rows:
        centroids = np.column_stack([result[f"centroid_{axis}"] for axis in ("x", "y", "z")])
        if len(rows) > 1:
            distances = cKDTree(centroids).query(centroids, k=2)[0][:, 1]
            neighbours = cKDTree(centroids).query(centroids, k=2)[1][:, 1]
            result["nn_colony_id"] = result["colony_id"][neighbours]
            result["nn_colony_distance"] = distances
        else:
            result["nn_colony_id"] = np.array([-1], dtype=np.int64)
            result["nn_colony_distance"] = np.array([np.nan])
    return result


def _genotypes(types: np.ndarray, loci: np.ndarray, lineage: np.ndarray | None) -> list[tuple[int, ...]]:
    if lineage is None:
        return [(int(kind), int(count)) for kind, count in zip(types, loci)]
    return [(int(kind), int(parent), int(count)) for kind, parent, count in zip(types, lineage, loci)]


def _positions(positions: np.ndarray) -> np.ndarray:
    points = np.asarray(positions, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("positions must have shape (N, 3)")
    return points


def _array(
    agents: dict[str, np.ndarray],
    name: str,
    size: int,
    dtype: np.dtype | type = float,
) -> np.ndarray:
    if name in agents:
        return np.asarray(agents[name], dtype=dtype)
    return np.zeros(size, dtype=dtype)


def _symmetric_graph(pairs: np.ndarray, size: int) -> coo_matrix:
    if pairs.size == 0:
        return coo_matrix((size, size))
    rows = np.concatenate([pairs[:, 0], pairs[:, 1]])
    cols = np.concatenate([pairs[:, 1], pairs[:, 0]])
    return coo_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(size, size))


def _write_table(path: Path, table: dict[str, np.ndarray]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        columns = list(table)
        writer.writerow(columns)
        for row in zip(*(table[column] for column in columns)):
            writer.writerow(row)


def main() -> None:
    """Write agent and colony CSV tables for an HDF5 step."""
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=validate_input_path)
    parser.add_argument("step")
    parser.add_argument("--agent-output", default="agents.csv")
    parser.add_argument("--colony-output", default="colonies.csv")
    parser.add_argument("--eps", type=float)
    parser.add_argument("--min-samples", type=int, default=3)
    args = parser.parse_args()
    config = ColonyConfig(eps=args.eps, min_samples=args.min_samples)
    with GutIBMData(args.input) as data:
        catalog = colony_catalog_from_hdf5(data, args.step, config)
    catalog.write_csv(args.agent_output, args.colony_output)
