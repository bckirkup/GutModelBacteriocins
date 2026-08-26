#!/usr/bin/env python3
"""Q_O2_per_gen analysis: per-cell funded oxygen per generation, and its regime.

Two independent estimators, both reported:
  Q_div = (funded growth O2 + funded maintenance O2) / cumulative divisions
  Q_mu  = funded O2 / cell_seconds * ln(2) / mean_mu_realized

Validity gates printed alongside, because Q is only usable as the divisor for
the lysis coefficient if the arm is genuinely uncrowded and chemically valid:
  rationing factor == 1.0, infeasible == 0, clips == 0, min oxygen >= 0.
"""
import json
import math
from pathlib import Path

import h5py
import numpy as np

ROOT = Path("/home/ubuntu/gutibm-campaign/q-o2-per-gen")
ARMS = ["Q_n2_res2", "Q_n2_res4", "Q_n2_res6",
        "Q_n8_res4", "Q_n20_res4", "Q_n80_res4"]


def sorted_keys(group):
    return sorted(group.keys(), key=lambda v: int(v.rsplit("_", 1)[1]))


def scalar(dataset):
    return float(np.asarray(dataset[()]).reshape(-1)[0])


def species_names(flux):
    return [bytes(r[r != 0]).decode() for r in flux["species_names"][:]]


def field(flux, name, index, default=0.0):
    return float(flux[name][index]) if name in flux else default


def cause_of(fh):
    value = fh["run_provenance"]["termination_cause"][()]
    if isinstance(value, bytes):
        return value.decode()
    return str(np.asarray(value).reshape(-1)[0])


def analyze(arm):
    cfg = json.loads((ROOT / arm / "input.json").read_text())
    with h5py.File(ROOT / arm / "output.h5", "r") as fh:
        keys = sorted_keys(fh["summary"])
        last = fh["summary"][keys[-1]]
        flux = last["nutrient_flux"]
        ox = species_names(flux).index("oxygen")

        # cell-seconds: trapezoidal integral of N over the summary trajectory
        times, counts = [], []
        for key in keys:
            rec = fh["summary"][key]
            times.append(scalar(rec["nutrient_flux"]["interval_end_time"]))
            counts.append(scalar(rec["n_total"]))
        cell_seconds = 0.0
        for (t0, n0), (t1, n1) in zip(zip(times, counts), zip(times[1:], counts[1:])):
            cell_seconds += 0.5 * (n0 + n1) * (t1 - t0)

        events = last["events"]
        divisions = scalar(events["cumulative_divisions"]) \
            if "cumulative_divisions" in events else 0.0
        funded = (field(flux, "agent_uptake_cumulative", ox)
                  + field(flux, "maintenance_cumulative", ox))
        mu = scalar(last["mean_mu_realized"]) \
            if "mean_mu_realized" in last else float("nan")

        # worst rationing factor seen anywhere in the run (minimum semantics)
        factor = min(field(fh["summary"][k]["nutrient_flux"],
                           "delivery_rationing_factor_cumulative", ox, 1.0)
                     for k in keys)
        infeasible = max(field(fh["summary"][k]["nutrient_flux"],
                               "delivery_infeasible_cumulative", ox, 0.0)
                         for k in keys)

        min_oxygen = float("nan")
        if "grid" in fh:
            grid_keys = sorted_keys(fh["grid"])
            mins = []
            for key in grid_keys:
                node = fh["grid"][key]
                if "oxygen" in node:
                    mins.append(float(np.min(node["oxygen"][()])))
            if mins:
                min_oxygen = min(mins)

        t_gen_div = cell_seconds / divisions if divisions > 0 else float("nan")
        t_gen_mu = math.log(2.0) / mu if mu and mu > 0 else float("nan")
        return {
            "arm": arm,
            "grid_dx_um": cfg["grid_dx"] * 1e6,
            "founders": cfg["initial_strains"][0]["count"],
            "final_N": scalar(last["n_total"]),
            "divisions": divisions,
            "cell_seconds": cell_seconds,
            "funded_oxygen": funded,
            "demanded_oxygen": field(flux, "uptake_demand_cumulative", ox),
            "funded_fraction": (funded / field(flux, "uptake_demand_cumulative", ox)
                                if field(flux, "uptake_demand_cumulative", ox) else float("nan")),
            "o2_per_cell_s": funded / cell_seconds if cell_seconds > 0 else float("nan"),
            "T_gen_div_h": t_gen_div / 3600.0,
            "T_gen_mu_h": t_gen_mu / 3600.0,
            "Q_div": funded / divisions if divisions > 0 else float("nan"),
    "Q_mu": (funded / cell_seconds * t_gen_mu
             if cell_seconds > 0 and math.isfinite(t_gen_mu) else float("nan")),
            "rationing_factor": factor,
            "infeasible": infeasible,
            "clips_oxygen": field(flux, "reaction_clip_cumulative", ox),
            "min_oxygen": min_oxygen,
            "termination": cause_of(fh),
        }


rows = [analyze(a) for a in ARMS if (ROOT / a / "output.h5").exists()]
(ROOT / "metrics.json").write_text(json.dumps(rows, indent=1))

cols = ["arm", "grid_dx_um", "founders", "final_N", "divisions", "funded_fraction",
        "o2_per_cell_s", "T_gen_div_h", "Q_div", "Q_mu", "rationing_factor",
        "infeasible", "clips_oxygen", "min_oxygen", "termination"]
print("\t".join(cols))
for r in rows:
    print("\t".join(f"{r[c]:.4g}" if isinstance(r[c], float) else str(r[c])
                    for c in cols))


def spread(subset, key):
    vals = [r[key] for r in subset if r[key] == r[key]]
    return max(vals) / min(vals) if vals and min(vals) > 0 else float("nan")


grid_arm = [r for r in rows if r["founders"] == 2]
print(f"\ngrid invariance at 2 founders: Q_div spread = {spread(grid_arm, 'Q_div'):.3f}"
      f", Q_mu spread = {spread(grid_arm, 'Q_mu'):.3f}")
ladder = sorted((r for r in rows if r["grid_dx_um"] == 4.0),
                key=lambda r: r["founders"])
print("density ladder at 4 um: "
      + ", ".join(f"N0={r['founders']} Q_div={r['Q_div']:.3e} "
                  f"factor={r['rationing_factor']:.4g}" for r in ladder))
