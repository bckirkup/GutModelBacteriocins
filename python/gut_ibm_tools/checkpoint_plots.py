"""Snapshot and timeseries plots for GutIBM checkpoint analysis."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .checkpoint_summary import DEAD_STATE
from .hdf5_reader import GutIBMData
from .path_utils import prepare_output_file

try:
    import matplotlib.pyplot as plt

    HAS_MPL = True
except ImportError:
    HAS_MPL = False

_MATPLOTLIB_REQUIRED = "matplotlib required: pip install 'gut-ibm-tools[viz]'"
UM = 1e6
HEATMAP_SPECIES = (
    "oxygen",
    "carbon",
    "acetate",
    "bacteriocin_BtuB",
    "bacteriocin_FepA",
    "bacteriocin_CirA",
    "bacteriocin_FhuA",
    "ai2",
    "iron",
)


def _require_mpl() -> None:
    if not HAS_MPL:
        raise ImportError(_MATPLOTLIB_REQUIRED)


def _save(fig: Any, path: Path) -> None:
    out = prepare_output_file(path)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _mid_slice(volume: np.ndarray) -> np.ndarray:
    """Return a 2D xz-ish mid-y slice; handle (nx,ny,nz) layout."""
    if volume.ndim != 3:
        return np.asarray(volume, dtype=float)
    mid = volume.shape[1] // 2
    return np.asarray(volume[:, mid, :], dtype=float).T


def plot_agents_scatter(
    data: GutIBMData,
    step: str,
    output: Path,
    *,
    projection: str = "xz",
) -> None:
    """2D scatter of alive agents colored by strain type."""
    _require_mpl()
    agents = data.get_agents(step)
    if "x" not in agents:
        return
    alive = agents.get("state", np.zeros(len(agents["x"]))) != DEAD_STATE
    x = agents["x"][alive] * UM
    y = agents["y"][alive] * UM
    z = agents["z"][alive] * UM
    types = agents["type"][alive] if "type" in agents else np.zeros(np.sum(alive))

    fig, ax = plt.subplots(figsize=(10, 6))
    if projection == "xy":
        scatter = ax.scatter(x, y, c=types, cmap="Set1", s=8, alpha=0.7)
        ax.set_xlabel("x (um)")
        ax.set_ylabel("y (um)")
    else:
        scatter = ax.scatter(x, z, c=types, cmap="Set1", s=8, alpha=0.7)
        ax.set_xlabel("x (um)")
        ax.set_ylabel("z (um)")
    plt.colorbar(scatter, ax=ax, label="type")
    summary = data.get_metadata(step)
    ax.set_title(f"Agents t={summary['time']:.0f}s N={int(np.sum(alive))}")
    ax.set_aspect("equal")
    _save(fig, output)


def plot_density_vs_depth(data: GutIBMData, step: str, output: Path) -> None:
    """Histogram of agent z-position (mucus depth)."""
    _require_mpl()
    agents = data.get_agents(step)
    if "z" not in agents:
        return
    alive = agents.get("state", np.zeros(len(agents["z"]))) != DEAD_STATE
    z = agents["z"][alive] * UM
    types = agents["type"][alive] if "type" in agents else np.zeros(len(z))

    fig, ax = plt.subplots(figsize=(8, 5))
    for type_id in sorted(np.unique(types)):
        ax.hist(
            z[types == type_id],
            bins=30,
            alpha=0.5,
            label=f"type {int(type_id)}",
        )
    ax.set_xlabel("z (um)")
    ax.set_ylabel("count")
    ax.set_title("Cell density vs depth")
    ax.legend()
    _save(fig, output)


def plot_grid_heatmaps(data: GutIBMData, step: str, out_dir: Path) -> list[Path]:
    """Write mid-y heatmaps for available chemical species."""
    _require_mpl()
    volumes = data.get_grid_volumes(step)
    written: list[Path] = []
    for name in HEATMAP_SPECIES:
        if name not in volumes:
            continue
        slice2d = _mid_slice(volumes[name])
        fig, ax = plt.subplots(figsize=(8, 4))
        im = ax.imshow(slice2d, origin="lower", aspect="auto", cmap="viridis")
        plt.colorbar(im, ax=ax, label=name)
        ax.set_xlabel("x index")
        ax.set_ylabel("z index")
        ax.set_title(f"{name} (mid-y slice)")
        path = out_dir / f"heatmap_{name}.png"
        _save(fig, path)
        written.append(path)
    return written


def plot_kill_zone(data: GutIBMData, step: str, output: Path) -> None:
    """Overlay agents on max bacteriocin field (BtuB preferred)."""
    _require_mpl()
    volumes = data.get_grid_volumes(step)
    toxin_name = next(
        (n for n in ("bacteriocin_BtuB", "bacteriocin_FepA", "bacteriocin_CirA", "bacteriocin_FhuA")
         if n in volumes),
        None,
    )
    if toxin_name is None:
        return
    toxin = _mid_slice(volumes[toxin_name])
    agents = data.get_agents(step)
    if "x" not in agents:
        return
    alive = agents.get("state", np.zeros(len(agents["x"]))) != DEAD_STATE
    nx, _ny, nz = data.grid_shape
    dx = data.grid_dx if data.grid_dx > 0 else 1.0
    x_idx = np.clip((agents["x"][alive] / dx).astype(int), 0, max(nx - 1, 0))
    z_idx = np.clip((agents["z"][alive] / dx).astype(int), 0, max(nz - 1, 0))
    types = agents["type"][alive] if "type" in agents else np.zeros(np.sum(alive))

    fig, ax = plt.subplots(figsize=(10, 5))
    im = ax.imshow(toxin, origin="lower", aspect="auto", cmap="magma", alpha=0.85)
    plt.colorbar(im, ax=ax, label=toxin_name)
    ax.scatter(x_idx, z_idx, c=types, cmap="Set1", s=10, edgecolors="white", linewidths=0.3)
    ax.set_xlabel("x index")
    ax.set_ylabel("z index")
    ax.set_title(f"Kill zone: {toxin_name} + agents")
    _save(fig, output)


def write_snapshot_figures(checkpoint: Path, out_dir: Path) -> list[Path]:
    """Generate all snapshot figures for one checkpoint into *out_dir*."""
    _require_mpl()
    written: list[Path] = []
    with GutIBMData(checkpoint) as data:
        steps = data.steps
        if not steps:
            raise ValueError(f"no steps in {checkpoint}")
        step = steps[-1]
        agents_path = out_dir / "agents_xz.png"
        plot_agents_scatter(data, step, agents_path, projection="xz")
        if agents_path.exists():
            written.append(agents_path)
        xy_path = out_dir / "agents_xy.png"
        plot_agents_scatter(data, step, xy_path, projection="xy")
        if xy_path.exists():
            written.append(xy_path)
        dens_path = out_dir / "density_vs_depth.png"
        plot_density_vs_depth(data, step, dens_path)
        if dens_path.exists():
            written.append(dens_path)
        written.extend(plot_grid_heatmaps(data, step, out_dir))
        kill_path = out_dir / "kill_zone.png"
        plot_kill_zone(data, step, kill_path)
        if kill_path.exists():
            written.append(kill_path)
    return written


def _plot_columns(
    frame: pd.DataFrame,
    columns: list[str],
    output: Path,
    *,
    title: str,
    ylabel: str,
) -> None:
    _require_mpl()
    present = [c for c in columns if c in frame.columns]
    if not present or "time" not in frame.columns:
        return
    fig, ax = plt.subplots(figsize=(10, 5))
    t_hours = frame["time"].to_numpy(dtype=float) / 3600.0
    for col in present:
        ax.plot(t_hours, frame[col].to_numpy(dtype=float), label=col, linewidth=1.5)
    ax.set_xlabel("time (h)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    _save(fig, output)


def write_timeseries_figures(frame: pd.DataFrame, out_dir: Path) -> list[Path]:
    """Plot standard trend panels from a summary DataFrame."""
    _require_mpl()
    written: list[Path] = []
    type_cols = sorted(c for c in frame.columns if c.startswith("n_type") or c.startswith("n_by_type_"))
    # Prefer agent-derived typeN counts; also plot n_by_type_* from summary.
    pop_cols = sorted(c for c in frame.columns if c.startswith("n_type"))
    if not pop_cols:
        pop_cols = [c for c in type_cols if c.startswith("n_by_type_")]
    if "n_alive" in frame.columns:
        pop_cols = ["n_alive", *pop_cols]
    elif "n_total" in frame.columns:
        pop_cols = ["n_total", *pop_cols]

    specs: list[tuple[str, list[str], str, str]] = [
        ("population.png", pop_cols, "Population", "count"),
        (
            "events_kills.png",
            [
                c
                for c in frame.columns
                if c.startswith("event_") and ("kill" in c or c.endswith("deaths"))
            ],
            "Kill / death events",
            "count",
        ),
        (
            "chem_means.png",
            [
                c
                for c in frame.columns
                if c.startswith("chem_mean_") or c.startswith("grid_mean_")
            ],
            "Chemical means",
            "concentration",
        ),
        (
            "spatial.png",
            [
                c
                for c in ("hopkins", "monochromatic_score", "rg_global", "mean_z")
                if c in frame.columns
            ],
            "Spatial metrics",
            "value",
        ),
    ]
    for filename, columns, title, ylabel in specs:
        path = out_dir / filename
        if not columns:
            continue
        _plot_columns(frame, columns, path, title=title, ylabel=ylabel)
        if path.exists():
            written.append(path)

    # Phase plot: producer-ish count vs bacteriocin
    x_col = next((c for c in ("n_type1", "n_by_type_1", "n_alive") if c in frame.columns), None)
    y_col = next(
        (c for c in ("grid_mean_bacteriocin_BtuB", "chem_max_toxin_BtuB") if c in frame.columns),
        None,
    )
    if x_col and y_col and HAS_MPL:
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.plot(frame[x_col], frame[y_col], "o-", markersize=3)
        ax.set_xlabel(x_col)
        ax.set_ylabel(y_col)
        ax.set_title("Phase: population vs bacteriocin")
        ax.grid(True, alpha=0.3)
        phase_path = out_dir / "phase_pop_vs_toxin.png"
        _save(fig, phase_path)
        written.append(phase_path)

    # Deduplicate while preserving order
    seen: set[Path] = set()
    unique: list[Path] = []
    for path in written:
        if path not in seen:
            seen.add(path)
            unique.append(path)
    return unique
