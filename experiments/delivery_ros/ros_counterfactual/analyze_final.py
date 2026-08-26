import argparse
import csv
import itertools
import json
import math
import re
from pathlib import Path

import h5py
import numpy as np

ROOT = Path("/home/ubuntu/gutibm-campaign/ros-counterfactual")
ARMS = [
    "A_ros0_res2",
    "A_ctrl_res2",
    "B_ros0_res2",
    "B_ctrl_res2",
    "B_ros0_res6",
    "B_ctrl_res6",
]
INVALID_ARM = "invalid_B_ros0_res2_100_founders"


def ordered_keys(group):
    return sorted(group.keys(), key=lambda value: int(value.rsplit("_", 1)[1]))


def scalar(dataset):
    value = np.asarray(dataset[()])
    return float(value.reshape(-1)[0])


def integer(dataset):
    value = np.asarray(dataset[()])
    return int(value.reshape(-1)[0])


def decode(value):
    if isinstance(value, bytes):
        return value.decode()
    value = np.asarray(value)
    if value.dtype.kind == "S":
        return value.reshape(-1)[0].decode()
    return str(value)


def species_names(flux):
    return [bytes(row[row != 0]).decode() for row in flux["species_names"][:]]


def mean_or_nan(values):
    return float(np.mean(values)) if len(values) else math.nan


def flux_value(flux, name, index):
    return float(flux[name][index]) if name in flux else 0.0


def event_count(events, name):
    return integer(events[name]) if name in events else 0


def event_rate(events, name):
    return scalar(events[name]) if name in events else 0.0


def find_dataset_paths(file, fragment):
    paths = []

    def visitor(name, value):
        if fragment in name:
            paths.append(name)

    file.visititems(visitor)
    return sorted(paths)


def parse_wrapper_wall(path):
    matches = re.findall(
        r"wrapper_wall_seconds=([0-9]+(?:\.[0-9]+)?)", path.read_text()
    )
    return float(matches[-1]) if matches else math.nan


def grid_agent_cells(grid, agents, nutrient, dx):
    values = np.asarray(grid[nutrient][:])
    result = []
    if agents is None:
        return result
    for x, y, z in zip(agents["x"][:], agents["y"][:], agents["z"][:]):
        ix = min(max(int(float(x) / dx), 0), values.shape[2] - 1)
        iy = min(max(int(float(y) / dx), 0), values.shape[1] - 1)
        iz = min(max(int(float(z) / dx), 0), values.shape[0] - 1)
        result.append(float(values[iz, iy, ix]))
    return result


def analyze_arm(arm):
    arm_dir = ROOT / arm
    config = json.loads((arm_dir / "input.json").read_text())
    path = arm_dir / "output.h5"
    with h5py.File(path, "r") as file:
        rows = []
        for key in ordered_keys(file["summary"]):
            record = file["summary"][key]
            flux = record["nutrient_flux"]
            names = species_names(flux)
            oxygen = names.index("oxygen")
            carbon = names.index("carbon")
            events = record["events"]
            stocks = record["stocks"]
            chem = record["chem"]
            component_names = (
                "sos_basal_rate",
                "sos_post_division_rate",
                "sos_nuclease_cross_induction_rate",
                "sos_ros_rate",
            )
            components = {
                name: event_rate(events, name) for name in component_names
            }
            cumulative_components = {
                name: event_rate(events, f"cumulative_{name}")
                for name in component_names
            }
            rows.append(
                {
                    "step": integer(record["step"]),
                    "time_s": scalar(record["time"]),
                    "N": integer(record["n_total"]),
                    "mean_fermentation_fraction": scalar(
                        record["mean_realized_fermentation_fraction"]
                    ),
                    "funded_growth_oxygen": flux_value(
                        flux, "agent_uptake_cumulative", oxygen
                    ),
                    "demanded_growth_oxygen": flux_value(
                        flux, "uptake_demand_cumulative", oxygen
                    ),
                    "funded_growth_oxygen_interval": flux_value(
                        flux, "agent_uptake_interval", oxygen
                    ),
                    "demanded_growth_oxygen_interval": flux_value(
                        flux, "uptake_demand_interval", oxygen
                    ),
                    "funded_maintenance_oxygen": flux_value(
                        flux, "maintenance_cumulative", oxygen
                    ),
                    "demanded_maintenance_oxygen": (
                        flux_value(flux, "maintenance_cumulative", oxygen)
                        + flux_value(
                            flux, "maintenance_shortfall_cumulative", oxygen
                        )
                    ),
                    "maintenance_oxygen_shortfall": flux_value(
                        flux, "maintenance_shortfall_cumulative", oxygen
                    ),
                    "agent_growth_oxygen_removal": flux_value(
                        flux, "agent_uptake_cumulative", oxygen
                    ),
                    "agent_maintenance_oxygen_removal": flux_value(
                        flux, "maintenance_cumulative", oxygen
                    ),
                    "flora_oxygen_removal": flux_value(
                        flux, "vbf_sink_cumulative", oxygen
                    ),
                    "boundary_oxygen_delivery": flux_value(
                        flux, "boundary_cumulative", oxygen
                    ),
                    "oxygen_reaction_clips": flux_value(
                        flux, "reaction_clip_cumulative", oxygen
                    ),
                    "funded_growth_carbon": flux_value(
                        flux, "agent_uptake_cumulative", carbon
                    ),
                    "demanded_growth_carbon": flux_value(
                        flux, "uptake_demand_cumulative", carbon
                    ),
                    "funded_maintenance_carbon": flux_value(
                        flux, "maintenance_cumulative", carbon
                    ),
                    "demanded_maintenance_carbon": (
                        flux_value(flux, "maintenance_cumulative", carbon)
                        + flux_value(
                            flux, "maintenance_shortfall_cumulative", carbon
                        )
                    ),
                    "maintenance_carbon_shortfall": flux_value(
                        flux, "maintenance_shortfall_cumulative", carbon
                    ),
                    "agent_growth_carbon_removal": flux_value(
                        flux, "agent_uptake_cumulative", carbon
                    ),
                    "agent_maintenance_carbon_removal": flux_value(
                        flux, "maintenance_cumulative", carbon
                    ),
                    "carbon_reaction_clips": flux_value(
                        flux, "reaction_clip_cumulative", carbon
                    ),
                    "divisions": event_count(events, "cumulative_divisions"),
                    "sos_inductions": event_count(
                        events, "cumulative_sos_inductions"
                    ),
                    "mortality_lysis": event_count(
                        events, "cumulative_mortality_lysis"
                    ),
                    "mortality_colicin": event_count(
                        events, "cumulative_mortality_colicin"
                    ),
                    "mortality_cdi": event_count(
                        events, "cumulative_mortality_cdi"
                    ),
                    "outflow_washout": event_count(
                        events, "cumulative_outflow_washout"
                    ),
                    "outflow_boundary": event_count(
                        events, "cumulative_outflow_boundary"
                    ),
                    "immigrations": event_count(
                        events, "cumulative_immigrations"
                    ),
                    "bacteriostatic_live_agents": event_count(
                        stocks, "bacteriostatic_live_agents"
                    ),
                    "washout_trapped_live_agents": event_count(
                        stocks, "washout_trapped_live_agents"
                    ),
                    "acetate_mean": (
                        scalar(chem["mean_acetate"])
                        if "mean_acetate" in chem
                        else math.nan
                    ),
                    "mean_oxygen": scalar(chem["mean_oxygen"]),
                    "mean_carbon": scalar(chem["mean_carbon"]),
                    "sos_components_interval": components,
                    "sos_components_cumulative": cumulative_components,
                    "sos_rate_total_interval": event_rate(
                        events, "sos_rate_total"
                    ),
                    "sos_rate_total_cumulative": event_rate(
                        events, "cumulative_sos_rate_total"
                    ),
                }
            )

        first = rows[0]
        final = rows[-1]
        total_live_steps = sum(row["N"] for row in rows)
        for row in rows:
            row["sos_components_interval_per_live_agent"] = {
                name: (
                    value / row["N"] if row["N"] else math.nan
                )
                for name, value in row["sos_components_interval"].items()
            }
            row["sos_components_cumulative_per_live_agent_step"] = {
                name: (
                    value / total_live_steps if total_live_steps else math.nan
                )
                for name, value in row["sos_components_cumulative"].items()
            }

        agent_trajectory = []
        for key in ordered_keys(file["agents"]):
            agents = file["agents"][key]
            agent_step = (
                integer(agents["step"])
                if "step" in agents
                else int(key.rsplit("_", 1)[1])
            )
            agent_trajectory.append(
                {
                    "step": agent_step,
                    "time_s": scalar(agents["time"])
                    if "time" in agents
                    else agent_step * config["bio_dt"],
                    "live_agents": len(agents["id"][:]),
                    "mean_biomass": mean_or_nan(
                        np.asarray(agents["biomass"][:])
                    ),
                    "mean_mu_realized": mean_or_nan(
                        np.asarray(agents["mu_realized"][:])
                    ),
                }
            )

        last_grid_key = ordered_keys(file["grid"])[-1]
        grid = file["grid"][last_grid_key]
        matching_agents = file["agents"].get(last_grid_key, None)
        final_grid = {}
        for nutrient in ("oxygen", "carbon", "acetate"):
            if nutrient not in grid:
                continue
            values = np.asarray(grid[nutrient][:])
            mean = float(values.mean())
            cells = grid_agent_cells(
                grid, matching_agents, nutrient, config["grid_dx"]
            )
            final_grid[nutrient] = {
                "domain_mean": mean,
                "domain_minimum": float(values.min()),
                "agent_cells": cells,
                "agent_cell_ratios": [
                    value / mean if mean else math.nan for value in cells
                ],
            }

        loss_channels = {
            "mortality_lysis": final["mortality_lysis"],
            "mortality_colicin": final["mortality_colicin"],
            "mortality_cdi": final["mortality_cdi"],
            "outflow_washout_trapped": final["outflow_washout"],
            "outflow_boundary": final["outflow_boundary"],
        }
        gross_losses = sum(loss_channels.values())
        net_decrease = first["N"] - final["N"]
        divisions = final["divisions"]
        immigrations = final["immigrations"]
        adjusted_net_decrease = net_decrease + divisions + immigrations
        area = config["domain_x"] * config["domain_y"]
        agent_seconds = sum(
            0.5 * (left["N"] + right["N"])
            * (right["time_s"] - left["time_s"])
            for left, right in itertools.pairwise(rows)
        )

        oxygen_funded = (
            final["funded_growth_oxygen"]
            + final["funded_maintenance_oxygen"]
        )
        oxygen_realized = (
            final["agent_growth_oxygen_removal"]
            + final["agent_maintenance_oxygen_removal"]
        )
        carbon_funded = (
            final["funded_growth_carbon"]
            + final["funded_maintenance_carbon"]
        )
        carbon_realized = (
            final["agent_growth_carbon_removal"]
            + final["agent_maintenance_carbon_removal"]
        )

        result = {
            "arm": arm,
            "grid_dx_um": config["grid_dx"] * 1e6,
            "oxygen_k_ROS": config["oxygen.k_ROS"],
            "oxygen_ros_driver": config.get(
                "oxygen.ros_driver", "ambient (default)"
            ),
            "configured_founders": sum(
                strain["count"] for strain in config["initial_strains"]
            ),
            "initial_N": first["N"],
            "final_N": final["N"],
            "max_N": max(row["N"] for row in rows),
            "termination_cause": decode(
                file["run_provenance/termination_cause"][()]
            ),
            "termination_detail": decode(
                file["run_provenance/termination_detail"][()]
            ),
            "termination_time_s": scalar(
                file["run_provenance/termination_time"]
            ),
            "termination_wall_seconds": scalar(
                file["run_provenance/termination_wall_seconds"]
            ),
            "wrapper_wall_seconds": parse_wrapper_wall(
                arm_dir / "full.log"
            ),
            "N_trajectory": [
                {"step": row["step"], "time_s": row["time_s"], "N": row["N"]}
                for row in rows
            ],
            "cumulative_divisions": divisions,
            "population_loss": {
                **loss_channels,
                "gross_loss_channel_sum": gross_losses,
                "net_N_decrease": net_decrease,
                "adjusted_N_decrease_for_divisions_and_immigrations": (
                    adjusted_net_decrease
                ),
                "gross_channel_sum_matches_adjusted_decrease": (
                    gross_losses == adjusted_net_decrease
                ),
            },
            "sos": {
                "cumulative_inductions": final["sos_inductions"],
                "interval_and_cumulative": [
                    {
                        "step": row["step"],
                        "time_s": row["time_s"],
                        "live_agents": row["N"],
                        "interval_components": row[
                            "sos_components_interval"
                        ],
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
                        "cumulative_components_per_live_agent_step": row[
                            "sos_components_cumulative_per_live_agent_step"
                        ],
                    }
                    for row in rows
                ],
                "cumulative_component_sum_matches_total": all(
                    math.isclose(
                        sum(item["cumulative_components"].values()),
                        item["cumulative_total"],
                        rel_tol=1e-12,
                        abs_tol=1e-30,
                    )
                    for item in (
                        {
                            "cumulative_components": row[
                                "sos_components_cumulative"
                            ],
                            "cumulative_total": row[
                                "sos_rate_total_cumulative"
                            ],
                        }
                        for row in rows
                    )
                ),
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
                "demanded_maintenance": final[
                    "demanded_maintenance_oxygen"
                ],
                "maintenance_shortfall": final[
                    "maintenance_oxygen_shortfall"
                ],
                "agent_growth_removal": final[
                    "agent_growth_oxygen_removal"
                ],
                "agent_maintenance_removal": final[
                    "agent_maintenance_oxygen_removal"
                ],
                "agent_total_removal": oxygen_realized,
                "flora_removal": final["flora_oxygen_removal"],
                "agent_share": (
                    oxygen_realized / (oxygen_realized + final["flora_oxygen_removal"])
                    if oxygen_realized + final["flora_oxygen_removal"]
                    else math.nan
                ),
                "boundary_delivery": final["boundary_oxygen_delivery"],
                "boundary_mean_flux": (
                    final["boundary_oxygen_delivery"]
                    / (area * final["time_s"])
                    if final["time_s"]
                    else math.nan
                ),
                "funded_total": oxygen_funded,
                "realized_total": oxygen_realized,
                "funded_minus_realized": oxygen_funded - oxygen_realized,
                "reaction_clips": final["oxygen_reaction_clips"],
                "reaction_clip_relative_to_growth_uptake": (
                    final["oxygen_reaction_clips"]
                    / final["funded_growth_oxygen"]
                    if final["funded_growth_oxygen"]
                    else math.nan
                ),
            },
            "carbon": {
                "funded_growth": final["funded_growth_carbon"],
                "demanded_growth": final["demanded_growth_carbon"],
                "funded_maintenance": final["funded_maintenance_carbon"],
                "demanded_maintenance": final[
                    "demanded_maintenance_carbon"
                ],
                "maintenance_shortfall": final[
                    "maintenance_carbon_shortfall"
                ],
                "agent_growth_removal": final[
                    "agent_growth_carbon_removal"
                ],
                "agent_maintenance_removal": final[
                    "agent_maintenance_carbon_removal"
                ],
                "agent_total_removal": carbon_realized,
                "funded_total": carbon_funded,
                "realized_total": carbon_realized,
                "funded_minus_realized": carbon_funded - carbon_realized,
                "reaction_clips": final["carbon_reaction_clips"],
                "reaction_clip_relative_to_growth_uptake": (
                    final["carbon_reaction_clips"]
                    / final["funded_growth_carbon"]
                    if final["funded_growth_carbon"]
                    else math.nan
                ),
            },
            "fermentation_fraction_mean": mean_or_nan(
                [row["mean_fermentation_fraction"] for row in rows[1:]]
            ),
            "fermentation_fraction_final": final[
                "mean_fermentation_fraction"
            ],
            "final_acetate_mean": final_grid.get("acetate", {}).get(
                "domain_mean", math.nan
            ),
            "acid_inhibition_enabled": bool(
                config["metabolism.acid_inhibition_enabled"]
            ),
            "mean_mu_realized_trajectory": [
                {
                    "step": item["step"],
                    "time_s": item["time_s"],
                    "live_agents": item["live_agents"],
                    "mean_mu_realized": item["mean_mu_realized"],
                }
                for item in agent_trajectory
            ],
            "mean_biomass_trajectory": [
                {
                    "step": item["step"],
                    "time_s": item["time_s"],
                    "live_agents": item["live_agents"],
                    "mean_biomass": item["mean_biomass"],
                }
                for item in agent_trajectory
            ],
            "live_agent_stocks": [
                {
                    "step": row["step"],
                    "time_s": row["time_s"],
                    "bacteriostatic_live_agents": row[
                        "bacteriostatic_live_agents"
                    ],
                    "washout_trapped_live_agents": row[
                        "washout_trapped_live_agents"
                    ],
                }
                for row in rows
            ],
            "final_grid": final_grid,
            "delivery_reduction_fields": find_dataset_paths(
                file, "delivery_reduction"
            ),
            "mucin_liberation_fields": find_dataset_paths(file, "mucin"),
            "raw_rows": rows,
            "agent_seconds": agent_seconds,
            "nonzero_clip_intervals": [
                {
                    "step": row["step"],
                    "time_s": row["time_s"],
                    "oxygen_clip": row["oxygen_reaction_clips"],
                    "oxygen_clip_relative": (
                        row["oxygen_reaction_clips"]
                        / row["funded_growth_oxygen"]
                        if row["funded_growth_oxygen"]
                        else math.nan
                    ),
                    "carbon_clip": row["carbon_reaction_clips"],
                    "carbon_clip_relative": (
                        row["carbon_reaction_clips"]
                        / row["funded_growth_carbon"]
                        if row["funded_growth_carbon"]
                        else math.nan
                    ),
                }
                for row in rows
                if row["oxygen_reaction_clips"] > 0.0
                or row["carbon_reaction_clips"] > 0.0
            ],
            "max_clip_relative_to_cumulative_growth_uptake": {
                "oxygen": max(
                    (
                        row["oxygen_reaction_clips"]
                        / row["funded_growth_oxygen"]
                        if row["funded_growth_oxygen"]
                        else 0.0
                    )
                    for row in rows
                ),
                "carbon": max(
                    (
                        row["carbon_reaction_clips"]
                        / row["funded_growth_carbon"]
                        if row["funded_growth_carbon"]
                        else 0.0
                    )
                    for row in rows
                ),
            },
            "clip_threshold_exceeded": any(
                (
                    row["oxygen_reaction_clips"]
                    / row["funded_growth_oxygen"]
                    if row["funded_growth_oxygen"]
                    else 0.0
                )
                > 1e-6
                or (
                    row["carbon_reaction_clips"]
                    / row["funded_growth_carbon"]
                    if row["funded_growth_carbon"]
                    else 0.0
                )
                > 1e-6
                for row in rows
            ),
        }
        return result


def flatten(result):
    loss = result["population_loss"]
    oxygen = result["oxygen"]
    carbon = result["carbon"]
    return {
        "arm": result["arm"],
        "grid_dx_um": result["grid_dx_um"],
        "oxygen_k_ROS": result["oxygen_k_ROS"],
        "configured_founders": result["configured_founders"],
        "initial_N": result["initial_N"],
        "final_N": result["final_N"],
        "max_N": result["max_N"],
        "cumulative_divisions": result["cumulative_divisions"],
        "termination_cause": result["termination_cause"],
        "termination_time_s": result["termination_time_s"],
        "termination_wall_seconds": result["termination_wall_seconds"],
        "wrapper_wall_seconds": result["wrapper_wall_seconds"],
        "mortality_lysis": loss["mortality_lysis"],
        "mortality_colicin": loss["mortality_colicin"],
        "mortality_cdi": loss["mortality_cdi"],
        "outflow_washout_trapped": loss["outflow_washout_trapped"],
        "outflow_boundary": loss["outflow_boundary"],
        "gross_loss_channel_sum": loss["gross_loss_channel_sum"],
        "net_N_decrease": loss["net_N_decrease"],
        "loss_closure": loss[
            "gross_channel_sum_matches_adjusted_decrease"
        ],
        "sos_inductions": result["sos"]["cumulative_inductions"],
        "oxygen_funded_growth": oxygen["funded_growth"],
        "oxygen_demanded_growth": oxygen["demanded_growth"],
        "oxygen_funded_growth_fraction": oxygen[
            "funded_growth_fraction"
        ],
        "oxygen_funded_maintenance": oxygen["funded_maintenance"],
        "oxygen_demanded_maintenance": oxygen["demanded_maintenance"],
        "oxygen_maintenance_shortfall": oxygen["maintenance_shortfall"],
        "oxygen_agent_total_removal": oxygen["agent_total_removal"],
        "oxygen_flora_removal": oxygen["flora_removal"],
        "oxygen_agent_share": oxygen["agent_share"],
        "oxygen_boundary_delivery": oxygen["boundary_delivery"],
        "oxygen_boundary_mean_flux": oxygen["boundary_mean_flux"],
        "oxygen_reaction_clips": oxygen["reaction_clips"],
        "oxygen_clip_relative": oxygen[
            "reaction_clip_relative_to_growth_uptake"
        ],
        "carbon_funded_growth": carbon["funded_growth"],
        "carbon_demanded_growth": carbon["demanded_growth"],
        "carbon_funded_maintenance": carbon["funded_maintenance"],
        "carbon_demanded_maintenance": carbon["demanded_maintenance"],
        "carbon_maintenance_shortfall": carbon["maintenance_shortfall"],
        "carbon_agent_total_removal": carbon["agent_total_removal"],
        "carbon_reaction_clips": carbon["reaction_clips"],
        "carbon_clip_relative": carbon[
            "reaction_clip_relative_to_growth_uptake"
        ],
        "fermentation_fraction_mean": result[
            "fermentation_fraction_mean"
        ],
        "fermentation_fraction_final": result[
            "fermentation_fraction_final"
        ],
        "final_acetate_mean": result["final_acetate_mean"],
        "acid_inhibition_enabled": result["acid_inhibition_enabled"],
        "delivery_reduction_fields": ",".join(
            result["delivery_reduction_fields"]
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-arm")
    args = parser.parse_args()
    if args.check_arm:
        result = analyze_arm(args.check_arm)
        print(
            args.check_arm,
            "clip_threshold_exceeded=",
            result["clip_threshold_exceeded"],
        )
        return
    results = [analyze_arm(arm) for arm in ARMS]
    invalid_record = {
        "arm": INVALID_ARM,
        "status": "invalid_superseded_run",
        "configured_founders": 100,
        "source_arm": "B_ros0_res2",
        "reason": (
            "Initial preparation used 80 founders for strain 1 while "
            "retaining 20 founders for strain 2."
        ),
        "preserved_output": str(ROOT / INVALID_ARM),
    }
    (ROOT / "metrics.json").write_text(
        json.dumps(
            {
                "arms": results,
                "invalid_arms": [invalid_record],
                "clip_policy": {
                    "relative_threshold": 1e-6,
                    "denominator": (
                        "cumulative agent growth uptake of the same species"
                    ),
                },
            },
            indent=2,
            allow_nan=True,
        )
        + "\n"
    )
    rows = [flatten(result) for result in results]
    with (ROOT / "metrics.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {ROOT / 'metrics.json'}")
    print(f"Wrote {ROOT / 'metrics.csv'}")


if __name__ == "__main__":
    main()
