import argparse
import csv
import itertools
import json
import math
from pathlib import Path

import h5py
import numpy as np

ROOT = Path("/home/ubuntu/gutibm-campaign/ros-counterfactual")


def keys(group):
    return sorted(group.keys(), key=lambda value: int(value.rsplit("_", 1)[1]))


def scalar(dataset):
    value = np.asarray(dataset[()])
    return float(value.reshape(-1)[0])


def integer(dataset):
    value = np.asarray(dataset[()])
    return int(value.reshape(-1)[0])


def text(value):
    if isinstance(value, bytes):
        return value.decode()
    value = np.asarray(value)
    if value.dtype.kind == "S":
        return value.reshape(-1)[0].decode()
    return str(value)


def names(flux):
    result = []
    for row in flux["species_names"][:]:
        result.append(bytes(row[row != 0]).decode())
    return result


def mean_or_nan(values):
    return float(np.mean(values)) if len(values) else math.nan


def field(flux, name, index, default=0.0):
    if name not in flux:
        return default
    return float(flux[name][index])


def event(events, name, default=0):
    if name not in events:
        return default
    return integer(events[name])


def event_rate(events, name, default=0.0):
    if name not in events:
        return default
    return scalar(events[name])


def trajectory_means(file, rows):
    for row in rows:
        step_name = f"step_{row['step']:06d}"
        if step_name not in file["agents"]:
            row["mean_biomass"] = math.nan
            row["mean_mu_realized"] = math.nan
            continue
        agents = file["agents"][step_name]
        row["mean_biomass"] = mean_or_nan(np.asarray(agents["biomass"][:]))
        row["mean_mu_realized"] = mean_or_nan(
            np.asarray(agents["mu_realized"][:])
        )


def delivery_reduction_paths(file):
    paths = []

    def visit(name, value):
        if "delivery_reduction" in name:
            paths.append(name)

    file.visititems(visit)
    return sorted(paths)


def analyze_arm(arm):
    path = ROOT / arm / "output.h5"
    with h5py.File(path, "r") as file:
        summary = file["summary"]
        rows = []
        for key in keys(summary):
            record = summary[key]
            flux = record["nutrient_flux"]
            species = names(flux)
            oxygen = species.index("oxygen")
            carbon = species.index("carbon")
            events = record["events"]
            stocks = record["stocks"]
            chem = record["chem"]
            components = {
                name: event_rate(events, name)
                for name in (
                    "sos_basal_rate",
                    "sos_post_division_rate",
                    "sos_nuclease_cross_induction_rate",
                    "sos_ros_rate",
                )
            }
            rows.append(
                {
                    "step": integer(record["step"]),
                    "time_s": scalar(record["time"]),
                    "N": integer(record["n_total"]),
                    "mean_fermentation_fraction": scalar(
                        record["mean_realized_fermentation_fraction"]
                    ),
                    "mean_mu_by_type": np.asarray(
                        record["mean_mu_by_type"][:]
                    ).tolist(),
                    "funded_growth_oxygen": field(
                        flux, "agent_uptake_cumulative", oxygen
                    ),
                    "demanded_growth_oxygen": field(
                        flux, "uptake_demand_cumulative", oxygen
                    ),
                    "funded_growth_oxygen_interval": field(
                        flux, "agent_uptake_interval", oxygen
                    ),
                    "demanded_growth_oxygen_interval": field(
                        flux, "uptake_demand_interval", oxygen
                    ),
                    "funded_maintenance_oxygen": field(
                        flux, "maintenance_cumulative", oxygen
                    ),
                    "maintenance_oxygen_shortfall": field(
                        flux, "maintenance_shortfall_cumulative", oxygen
                    ),
                    "agent_oxygen_removal": field(
                        flux, "agent_uptake_cumulative", oxygen
                    ),
                    "flora_oxygen_removal": field(
                        flux, "vbf_sink_cumulative", oxygen
                    ),
                    "boundary_oxygen_delivery": field(
                        flux, "boundary_cumulative", oxygen
                    ),
                    "oxygen_reaction_clips": field(
                        flux, "reaction_clip_cumulative", oxygen
                    ),
                    "funded_growth_carbon": field(
                        flux, "agent_uptake_cumulative", carbon
                    ),
                    "demanded_growth_carbon": field(
                        flux, "uptake_demand_cumulative", carbon
                    ),
                    "funded_growth_carbon_interval": field(
                        flux, "agent_uptake_interval", carbon
                    ),
                    "demanded_growth_carbon_interval": field(
                        flux, "uptake_demand_interval", carbon
                    ),
                    "funded_maintenance_carbon": field(
                        flux, "maintenance_cumulative", carbon
                    ),
                    "maintenance_carbon_shortfall": field(
                        flux, "maintenance_shortfall_cumulative", carbon
                    ),
                    "agent_carbon_removal": field(
                        flux, "agent_uptake_cumulative", carbon
                    ),
                    "carbon_reaction_clips": field(
                        flux, "reaction_clip_cumulative", carbon
                    ),
                    "divisions": event(events, "cumulative_divisions"),
                    "sos_inductions": event(events, "cumulative_sos_inductions"),
                    "mortality_colicin": event(
                        events, "cumulative_mortality_colicin"
                    ),
                    "mortality_cdi": event(events, "cumulative_mortality_cdi"),
                    "mortality_lysis": event(
                        events, "cumulative_mortality_lysis"
                    ),
                    "outflow_boundary": event(
                        events, "cumulative_outflow_boundary"
                    ),
                    "outflow_washout": event(
                        events, "cumulative_outflow_washout"
                    ),
                    "bacteriostatic_live_agents": event(
                        stocks, "bacteriostatic_live_agents"
                    ),
                    "washout_trapped_live_agents": event(
                        stocks, "washout_trapped_live_agents"
                    ),
                    "acetate_mean": scalar(chem["mean_acetate"])
                    if "mean_acetate" in chem
                    else math.nan,
                    "mean_oxygen": scalar(chem["mean_oxygen"]),
                    "mean_carbon": scalar(chem["mean_carbon"]),
                    "sos_components_interval": components,
                    "sos_rate_total_interval": event_rate(
                        events, "sos_rate_total"
                    ),
                    "sos_components_cumulative": {
                        name: event_rate(events, f"cumulative_{name}")
                        for name in components
                    },
                    "sos_rate_total_cumulative": event_rate(
                        events, "cumulative_sos_rate_total"
                    ),
                }
            )

        trajectory_means(file, rows)
        first = rows[0]
        final = rows[-1]
        config = json.loads((ROOT / arm / "input.json").read_text())
        agent_seconds = sum(
            0.5 * (left["N"] + right["N"])
            * (right["time_s"] - left["time_s"])
            for left, right in itertools.pairwise(rows)
        )
        live_means = []
        for row in rows:
            components = row["sos_components_interval"]
            live_means.append(
                {
                    name: value / row["N"] if row["N"] else math.nan
                    for name, value in components.items()
                }
            )
            row["sos_components_interval_per_live_agent"] = live_means[-1]

        last_grid = keys(file["grid"])[-1]
        grid = file["grid"][last_grid]
        final_grid_metrics = {}
        for nutrient in ("oxygen", "carbon"):
            values = np.asarray(grid[nutrient][:])
            final_grid_metrics[f"final_{nutrient}_mean"] = float(values.mean())
            final_grid_metrics[f"final_{nutrient}_minimum"] = float(values.min())

        final_agent_step = last_grid if last_grid in file["agents"] else None
        if final_agent_step is not None:
            agents = file["agents"][final_agent_step]
            for nutrient in ("oxygen", "carbon"):
                values = np.asarray(grid[nutrient][:])
                dx = config["grid_dx"]
                cell_values = []
                for x, y, z in zip(
                    agents["x"][:], agents["y"][:], agents["z"][:]
                ):
                    ix = min(max(int(float(x) / dx), 0), values.shape[2] - 1)
                    iy = min(max(int(float(y) / dx), 0), values.shape[1] - 1)
                    iz = min(max(int(float(z) / dx), 0), values.shape[0] - 1)
                    cell_values.append(float(values[iz, iy, ix]))
                final_grid_metrics[f"final_{nutrient}_agent_cells"] = cell_values
                mean = final_grid_metrics[f"final_{nutrient}_mean"]
                final_grid_metrics[f"final_{nutrient}_agent_cell_ratios"] = [
                    value / mean if mean else math.nan for value in cell_values
                ]
        provenance = file["run_provenance"]
        result = {
            "arm": arm,
            "grid_dx_um": config["grid_dx"] * 1e6,
            "oxygen_k_ROS": config["oxygen.k_ROS"],
            "oxygen_vbf_sink": config["oxygen.vbf_sink"],
            "initial_N": first["N"],
            "final_N": final["N"],
            "max_N": max(row["N"] for row in rows),
            "N_decrease": first["N"] - final["N"],
            "termination_cause": text(provenance["termination_cause"][()]),
            "termination_detail": text(provenance["termination_detail"][()]),
            "termination_time_s": scalar(provenance["termination_time"]),
            "termination_wall_seconds": scalar(
                provenance["termination_wall_seconds"]
            ),
            "N_trajectory": [
                {"step": row["step"], "time_s": row["time_s"], "N": row["N"]}
                for row in rows
            ],
            "cumulative_divisions": final["divisions"],
            "population_loss": {
                "mortality_lysis": final["mortality_lysis"],
                "mortality_colicin": final["mortality_colicin"],
                "mortality_cdi": final["mortality_cdi"],
                "outflow_washout_trapped": final["outflow_washout"],
                "outflow_boundary": final["outflow_boundary"],
                "channel_sum": (
                    final["mortality_lysis"]
                    + final["mortality_colicin"]
                    + final["mortality_cdi"]
                    + final["outflow_washout"]
                    + final["outflow_boundary"]
                ),
                "closure_matches_N_decrease": (
                    first["N"] - final["N"]
                    == final["mortality_lysis"]
                    + final["mortality_colicin"]
                    + final["mortality_cdi"]
                    + final["outflow_washout"]
                    + final["outflow_boundary"]
                ),
            },
            "sos": {
                "cumulative_inductions": final["sos_inductions"],
                "interval_and_cumulative": [
                    {
                        "step": row["step"],
                        "time_s": row["time_s"],
                        "live_agents": row["N"],
                        "interval_components": row["sos_components_interval"],
                        "interval_total": row["sos_rate_total_interval"],
                        "interval_components_per_live_agent": row[
                            "sos_components_interval_per_live_agent"
                        ],
                        "cumulative_components": row[
                            "sos_components_cumulative"
                        ],
                        "cumulative_total": row[
                            "sos_rate_total_cumulative"
                        ],
                    }
                    for row in rows
                ],
            },
            "oxygen": {
                "funded_growth": final["funded_growth_oxygen"],
                "demanded_growth": final["demanded_growth_oxygen"],
                "funded_growth_fraction": (
                    final["funded_growth_oxygen"]
                    / final["demanded_growth_oxygen"]
                    if final["demanded_growth_oxygen"]
                    else math.nan
                ),
                "funded_maintenance": final["funded_maintenance_oxygen"],
                "demanded_maintenance": (
                    final["funded_maintenance_oxygen"]
                    + final["maintenance_oxygen_shortfall"]
                ),
                "maintenance_shortfall": final["maintenance_oxygen_shortfall"],
                "agent_removal": final["agent_oxygen_removal"],
                "flora_removal": final["flora_oxygen_removal"],
                "agent_share": (
                    final["agent_oxygen_removal"]
                    / (
                        final["agent_oxygen_removal"]
                        + final["flora_oxygen_removal"]
                    )
                    if final["agent_oxygen_removal"]
                    + final["flora_oxygen_removal"]
                    else math.nan
                ),
                "funded_minus_agent_removal": (
                    final["funded_growth_oxygen"]
                    - final["agent_oxygen_removal"]
                ),
                "boundary_delivery": final["boundary_oxygen_delivery"],
                "reaction_clips": final["oxygen_reaction_clips"],
            },
            "carbon": {
                "funded_growth": final["funded_growth_carbon"],
                "demanded_growth": final["demanded_growth_carbon"],
                "funded_maintenance": final["funded_maintenance_carbon"],
                "demanded_maintenance": (
                    final["funded_maintenance_carbon"]
                    + final["maintenance_carbon_shortfall"]
                ),
                "maintenance_shortfall": final["maintenance_carbon_shortfall"],
                "agent_removal": final["agent_carbon_removal"],
                "funded_minus_agent_removal": (
                    final["funded_growth_carbon"]
                    - final["agent_carbon_removal"]
                ),
                "reaction_clips": final["carbon_reaction_clips"],
            },
            "fermentation_fraction_mean": mean_or_nan(
                [row["mean_fermentation_fraction"] for row in rows[1:]]
            ),
            "fermentation_fraction_final": final[
                "mean_fermentation_fraction"
            ],
            "final_acetate_mean": final["acetate_mean"],
            "mean_mu_realized_trajectory": [
                {
                    "step": row["step"],
                    "time_s": row["time_s"],
                    "mean_mu_realized": row["mean_mu_realized"],
                }
                for row in rows
            ],
            "mean_biomass_trajectory": [
                {
                    "step": row["step"],
                    "time_s": row["time_s"],
                    "mean_biomass": row["mean_biomass"],
                }
                for row in rows
            ],
            "final_grid": final_grid_metrics,
            "delivery_reduction_fields": delivery_reduction_paths(file),
            "agent_seconds": agent_seconds,
            "raw_rows": rows,
        }
        if result["oxygen"]["reaction_clips"] != 0.0:
            raise RuntimeError(f"{arm}: nonzero oxygen reaction clip")
        if result["carbon"]["reaction_clips"] != 0.0:
            raise RuntimeError(f"{arm}: nonzero carbon reaction clip")
        if (
            result["oxygen"]["funded_minus_agent_removal"] > 1.0e-30
            or result["carbon"]["funded_minus_agent_removal"] > 1.0e-30
        ):
            raise RuntimeError(f"{arm}: funded uptake exceeds realized removal")
        return result


def write_outputs(results):
    (ROOT / "metrics.json").write_text(
        json.dumps(results, indent=2, allow_nan=True) + "\n"
    )
    rows = []
    for result in results:
        rows.append(
            {
                "arm": result["arm"],
                "grid_dx_um": result["grid_dx_um"],
                "oxygen_k_ROS": result["oxygen_k_ROS"],
                "oxygen_vbf_sink": result["oxygen_vbf_sink"],
                "initial_N": result["initial_N"],
                "final_N": result["final_N"],
                "N_decrease": result["N_decrease"],
                "cumulative_divisions": result["cumulative_divisions"],
                "termination_cause": result["termination_cause"],
                "termination_time_s": result["termination_time_s"],
                "termination_wall_seconds": result[
                    "termination_wall_seconds"
                ],
                "oxygen_funded_growth": result["oxygen"]["funded_growth"],
                "oxygen_demanded_growth": result["oxygen"]["demanded_growth"],
                "oxygen_funded_growth_fraction": result["oxygen"][
                    "funded_growth_fraction"
                ],
                "oxygen_funded_maintenance": result["oxygen"][
                    "funded_maintenance"
                ],
                "oxygen_demanded_maintenance": result["oxygen"][
                    "demanded_maintenance"
                ],
                "oxygen_maintenance_shortfall": result["oxygen"][
                    "maintenance_shortfall"
                ],
                "oxygen_agent_removal": result["oxygen"]["agent_removal"],
                "oxygen_flora_removal": result["oxygen"]["flora_removal"],
                "oxygen_agent_share": result["oxygen"]["agent_share"],
                "oxygen_boundary_delivery": result["oxygen"][
                    "boundary_delivery"
                ],
                "carbon_funded_growth": result["carbon"]["funded_growth"],
                "carbon_demanded_growth": result["carbon"][
                    "demanded_growth"
                ],
                "carbon_funded_maintenance": result["carbon"][
                    "funded_maintenance"
                ],
                "carbon_demanded_maintenance": result["carbon"][
                    "demanded_maintenance"
                ],
                "carbon_maintenance_shortfall": result["carbon"][
                    "maintenance_shortfall"
                ],
                "carbon_agent_removal": result["carbon"]["agent_removal"],
                "fermentation_fraction_mean": result[
                    "fermentation_fraction_mean"
                ],
                "fermentation_fraction_final": result[
                    "fermentation_fraction_final"
                ],
                "final_acetate_mean": result["final_acetate_mean"],
                "oxygen_reaction_clips": result["oxygen"]["reaction_clips"],
                "carbon_reaction_clips": result["carbon"]["reaction_clips"],
                "mortality_lysis": result["population_loss"][
                    "mortality_lysis"
                ],
                "mortality_colicin": result["population_loss"][
                    "mortality_colicin"
                ],
                "mortality_cdi": result["population_loss"]["mortality_cdi"],
                "outflow_washout_trapped": result["population_loss"][
                    "outflow_washout_trapped"
                ],
                "outflow_boundary": result["population_loss"][
                    "outflow_boundary"
                ],
                "loss_channel_sum": result["population_loss"]["channel_sum"],
                "loss_closure": result["population_loss"][
                    "closure_matches_N_decrease"
                ],
                "delivery_reduction_fields": ",".join(
                    result["delivery_reduction_fields"]
                ),
            }
        )
    with (ROOT / "metrics.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-arm")
    args = parser.parse_args()
    if args.check_arm:
        analyze_arm(args.check_arm)
        return
    arms = [
        "A_ros0_res2",
        "A_ctrl_res2",
        "B_ros0_res2",
        "B_ctrl_res2",
        "B_ros0_res6",
        "B_ctrl_res6",
    ]
    results = [analyze_arm(arm) for arm in arms]
    write_outputs(results)
    print(f"Wrote {ROOT / 'metrics.json'}")
    print(f"Wrote {ROOT / 'metrics.csv'}")


if __name__ == "__main__":
    main()
