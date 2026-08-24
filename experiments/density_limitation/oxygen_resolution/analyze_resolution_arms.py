#!/usr/bin/env python3
"""Regenerate oxygen resolution sweep metrics from HDF5 arm outputs."""

import argparse
import csv
import json
import math
import re
from pathlib import Path

import h5py
import numpy as np

ARM_NAMES = ("res2_dt60", "res4_dt60", "res6_dt60", "res2_dt10")
CHECKPOINT_TIMES = (600.0, 1800.0, 3600.0)


def _steps(handle: h5py.File) -> list[str]:
    return sorted(
        handle["summary"].keys(),
        key=lambda key: int(key.rsplit("_", 1)[-1]),
    )


def _species_names(flux: h5py.Group) -> list[str]:
    return [
        bytes(row[row != 0]).decode()
        for row in flux["species_names"][:]
    ]


def _value(flux: h5py.Group, key: str, index: int) -> float:
    return float(flux[key][index]) if key in flux else 0.0


def _decode(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode()
    return str(value)


def _profile_depth(profile: np.ndarray, grid_dx: float) -> float:
    threshold = float(profile[0]) / math.e
    for index in range(1, len(profile)):
        if profile[index] <= threshold:
            before = float(profile[index - 1])
            after = float(profile[index])
            fraction = ((threshold - before) / (after - before)
                        if after != before else 0.0)
            return ((index - 1) + fraction) * grid_dx * 1.0e6
    return float("nan")


def _wall_time(log_path: Path) -> float:
    matches = re.findall(
        r"WALL_TIME_SECONDS=([0-9.eE+-]+)",
        log_path.read_text(),
    )
    return float(matches[-1]) if matches else float("nan")


def _row(group: h5py.Group, previous: dict | None,
         agent_seconds: float) -> tuple[dict, float]:
    step = int(group["step"][0])
    time_s = float(group["time"][0])
    live_agents = float(group["num_agents"][0])
    dt_s = float(group["dt"][0])
    if previous is not None:
        agent_seconds += 0.5 * (
            previous["live_agents"] + live_agents
        ) * (time_s - previous["time_s"])

    flux = group["nutrient_flux"]
    species = _species_names(flux)
    oxygen = species.index("oxygen")
    carbon = species.index("carbon")
    funded_cumulative = _value(
        flux, "agent_uptake_cumulative", oxygen
    )
    demanded_cumulative = _value(
        flux, "uptake_demand_cumulative", oxygen
    )
    funded_interval = _value(flux, "agent_uptake_interval", oxygen)
    demanded_interval = _value(flux, "uptake_demand_interval", oxygen)
    row = {
        "step": step,
        "time_s": time_s,
        "dt_s": dt_s,
        "live_agents": live_agents,
        "funded_growth_o2_cumulative_mol": funded_cumulative,
        "demanded_growth_o2_cumulative_mol": demanded_cumulative,
        "funded_growth_o2_interval_mol": funded_interval,
        "demanded_growth_o2_interval_mol": demanded_interval,
        "funded_fraction_cumulative": (
            funded_cumulative / demanded_cumulative
            if demanded_cumulative > 0.0 else float("nan")
        ),
        "funded_fraction_interval": (
            funded_interval / demanded_interval
            if demanded_interval > 0.0 else float("nan")
        ),
        "funded_growth_o2_per_live_agent_s_mol": (
            funded_cumulative / agent_seconds
            if agent_seconds > 0.0 else float("nan")
        ),
        "maintenance_paid_o2_cumulative_mol": _value(
            flux, "maintenance_cumulative", oxygen
        ),
        "maintenance_paid_o2_interval_mol": _value(
            flux, "maintenance_interval", oxygen
        ),
        "maintenance_shortfall_o2_cumulative_mol": _value(
            flux, "maintenance_shortfall_cumulative", oxygen
        ),
        "maintenance_shortfall_o2_interval_mol": _value(
            flux, "maintenance_shortfall_interval", oxygen
        ),
        "agent_o2_realized_cumulative_mol": funded_cumulative,
        "agent_o2_realized_interval_mol": funded_interval,
        "flora_o2_realized_cumulative_mol": _value(
            flux, "vbf_sink_cumulative", oxygen
        ),
        "flora_o2_realized_interval_mol": _value(
            flux, "vbf_sink_interval", oxygen
        ),
        "oxygen_reaction_clip_cumulative_mol": _value(
            flux, "reaction_clip_cumulative", oxygen
        ),
        "oxygen_reaction_clip_interval_mol": _value(
            flux, "reaction_clip_interval", oxygen
        ),
        "carbon_reaction_clip_cumulative_mol": _value(
            flux, "reaction_clip_cumulative", carbon
        ),
        "carbon_reaction_clip_interval_mol": _value(
            flux, "reaction_clip_interval", carbon
        ),
        "mean_realized_fermentation_fraction": float(
            group["mean_realized_fermentation_fraction"][0]
        ),
        "mean_oxygen_mol_per_m3": float(group["chem/mean_oxygen"][0]),
        "agent_seconds_to_time": agent_seconds,
    }
    return row, agent_seconds


def _extract_arm(arm_dir: Path) -> tuple[list[dict], dict]:
    with h5py.File(arm_dir / "output.h5", "r") as handle:
        steps = _steps(handle)
        rows = []
        previous = None
        agent_seconds = 0.0
        for step_name in steps:
            row, agent_seconds = _row(
                handle["summary"][step_name],
                previous,
                agent_seconds,
            )
            rows.append(row)
            previous = row

        provenance = handle["run_provenance"]
        last_grid = max(
            handle["grid"].keys(),
            key=lambda key: int(key.rsplit("_", 1)[-1]),
        )
        grid = handle["grid"][last_grid]["oxygen"][:]
        profile = grid.mean(axis=(1, 2))
        with (arm_dir / "input.json").open() as config_file:
            config = json.load(config_file)
        grid_dx = float(config["grid_dx"])
        final = {
            "termination_cause": _decode(
                provenance["termination_cause"][()]
            ),
            "termination_step": int(provenance["termination_step"][0]),
            "termination_time_s": float(provenance["termination_time"][0]),
            "git_sha": _decode(provenance["git_sha"][()]),
            "grid_step": int(last_grid.rsplit("_", 1)[-1]),
            "grid_time_s": float(last_grid.rsplit("_", 1)[-1])
            * float(config["bio_dt"]),
            "oxygen_surface_mol_per_m3": float(profile[0]),
            "oxygen_one_over_e_depth_um": _profile_depth(profile, grid_dx),
            "wall_time_s": _wall_time(arm_dir / "full.log"),
        }
    return rows, final


def _nearest(rows: list[dict], target: float) -> dict:
    return min(rows, key=lambda row: abs(row["time_s"] - target))


def _write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main(run_root: Path, output_root: Path) -> None:
    all_rows = []
    final_rows = []
    for arm in ARM_NAMES:
        rows, final = _extract_arm(run_root / arm)
        for row in rows:
            row["arm"] = arm
        all_rows.extend(rows)
        final["arm"] = arm
        for target in CHECKPOINT_TIMES:
            checkpoint = _nearest(rows, target)
            for key, value in checkpoint.items():
                if key not in {"step", "time_s", "dt_s"}:
                    final[f"{key}_at_{int(target)}s"] = value
        final_rows.append(final)

    output_root.mkdir(parents=True, exist_ok=True)
    _write_csv(output_root / "resolution_timeseries.csv", all_rows)
    _write_csv(output_root / "resolution_metrics.csv", final_rows)
    for arm in ARM_NAMES:
        arm_dir = run_root / arm
        with h5py.File(arm_dir / "output.h5", "r") as handle:
            last_grid = max(
                handle["grid"].keys(),
                key=lambda key: int(key.rsplit("_", 1)[-1]),
            )
            profile = handle["grid"][last_grid]["oxygen"][:].mean(axis=(1, 2))
        with (arm_dir / "input.json").open() as config_file:
            grid_dx = float(json.load(config_file)["grid_dx"])
        profile_path = output_root / f"{arm}_oxygen_profile_last.csv"
        with profile_path.open("w", newline="") as output:
            writer = csv.writer(output)
            writer.writerow(["z_index", "z_um", "oxygen_mean_mol_per_m3"])
            writer.writerows(
                (index, index * grid_dx * 1.0e6, float(value))
                for index, value in enumerate(profile)
            )
    print(output_root / "resolution_timeseries.csv")
    print(output_root / "resolution_metrics.csv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    arguments = parser.parse_args()
    main(arguments.run_root.expanduser().resolve(),
         arguments.output_root.expanduser().resolve())
