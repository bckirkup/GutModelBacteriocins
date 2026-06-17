"""
Visualization utilities for GutIBM output.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    HAS_MPL = True
except ImportError:
    HAS_MPL = False

from .hdf5_reader import GutIBMData


def plot_agent_positions(
    data: GutIBMData,
    step: str,
    projection: str = "xy",
    output: Optional[str] = None,
) -> None:
    """Plot agent positions colored by type."""
    if not HAS_MPL:
        raise ImportError("matplotlib required for visualization")

    agents = data.get_agents(step)
    x, y, z = agents["x"], agents["y"], agents["z"]
    types = agents["type"]
    meta = data.get_metadata(step)

    fig, ax = plt.subplots(figsize=(10, 8))

    if projection == "xy":
        scatter = ax.scatter(x * 1e6, y * 1e6, c=types, cmap="Set1",
                             s=2, alpha=0.6)
        ax.set_xlabel("x (um)")
        ax.set_ylabel("y (um)")
    elif projection == "xz":
        scatter = ax.scatter(x * 1e6, z * 1e6, c=types, cmap="Set1",
                             s=2, alpha=0.6)
        ax.set_xlabel("x (um) [distal]")
        ax.set_ylabel("z (um) [epithelium -> lumen]")

    plt.colorbar(scatter, ax=ax, label="Agent type")
    ax.set_title(f"t = {meta['time']:.0f}s, N = {meta['num_agents']}")
    ax.set_aspect("equal")

    if output:
        fig.savefig(output, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()


def plot_population_timeseries(
    data: GutIBMData,
    output: Optional[str] = None,
) -> None:
    """Plot population size over time."""
    if not HAS_MPL:
        raise ImportError("matplotlib required for visualization")

    times, counts = data.time_series("num_agents")

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(times / 3600, counts, "b-", linewidth=1.5)
    ax.set_xlabel("Time (hours)")
    ax.set_ylabel("Number of agents")
    ax.set_title("Population Dynamics")
    ax.grid(True, alpha=0.3)

    if output:
        fig.savefig(output, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()


def plot_lineage_composition(
    data: GutIBMData,
    output: Optional[str] = None,
) -> None:
    """Plot lineage composition (stacked area) over time."""
    if not HAS_MPL:
        raise ImportError("matplotlib required for visualization")

    steps = data.steps
    times = []
    lineage_data: dict[int, list[int]] = {}

    for step in steps:
        meta = data.get_metadata(step)
        agents = data.get_agents(step)
        times.append(meta["time"])

        lineages = agents.get("lineage", np.array([]))
        unique, counts = np.unique(lineages, return_counts=True)
        for lin, cnt in zip(unique, counts):
            lin_int = int(lin)
            if lin_int not in lineage_data:
                lineage_data[lin_int] = [0] * len(times[:-1])
            lineage_data[lin_int].append(int(cnt))

        # Pad missing lineages
        for lin in lineage_data:
            if len(lineage_data[lin]) < len(times):
                lineage_data[lin].append(0)

    fig, ax = plt.subplots(figsize=(12, 6))
    t_hours = np.array(times) / 3600

    labels = sorted(lineage_data.keys())
    data_arrays = np.array([lineage_data[lb] for lb in labels])

    ax.stackplot(t_hours, data_arrays, labels=[f"L{lb}" for lb in labels],
                 alpha=0.7)
    ax.set_xlabel("Time (hours)")
    ax.set_ylabel("Agent count")
    ax.set_title("Lineage Composition Over Time")
    ax.legend(loc="upper right", fontsize=8, ncol=3)

    if output:
        fig.savefig(output, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()
