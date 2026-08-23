"""Audit whether agent oxygen respiration is actually removed from the O2 field,
or refused by positivity clipping after the explicit reaction update.

Lead-authored measurement (Spec 13 review seam S11). No new runs are needed:
every supply-ladder arm (PR #314) already ran with `oxygen.enabled = true` and
`oxygen.metabolic_switch_enabled = true`, and per-species reaction clips are
accounted, so the answer is already in those HDF5 files.

For each arm this reports, at the last written step:

  reaction_clip_cumulative[oxygen]   mol of O2 sink refused by conc >= 0
  boundary_cumulative[oxygen]        mol of O2 delivered across boundaries
  reaction_clip_cumulative[carbon]   the same quantity for carbon (delivery mode)
  mean_oxygen first vs last          whether agents deplete the field at all

and an independent estimate of the respiratory demand, integrated over the
population trajectory from the summary series:

  demand = sum_i (q_maintenance + q_consumption * mu_i * (1 - ferm_i)) * N_i * dt_i

with the shipped defaults q_maintenance = 1e-18 mol/s/cell and
q_consumption = 1e-14 mol/cell. `mu_i` is the population-weighted mean of
`mean_mu_by_type`, so the estimate is coarse; it is an order-of-magnitude
reference for the clip, not a closure identity.

Interpretation: clip / demand ~ 1 means the entire respiratory sink is refused,
i.e. the pre-#308 carbon defect, unfixed for oxygen. clip >> boundary means the
refused amount exceeds all oxygen the domain ever received.
"""

import argparse
import pathlib

import h5py
import numpy as np

Q_MAINTENANCE = 1.0e-18  # mol/s/cell, oxygen.q_maintenance default
Q_CONSUMPTION = 1.0e-14  # mol/cell, oxygen.q_consumption default


def species_names(flux_group):
    return [bytes(row[row != 0]).decode() for row in flux_group["species_names"][:]]


def audit(path):
    with h5py.File(path, "r") as handle:
        summary = handle["summary"]
        steps = sorted(summary.keys())
        times, counts, mus, ferms, oxy = [], [], [], [], []
        for step in steps:
            group = summary[step]
            times.append(group["nutrient_flux/interval_end_time"][0])
            counts.append(group["num_agents"][0])
            by_type = group["n_by_type"][:]
            mu_by_type = group["mean_mu_by_type"][:]
            total = max(by_type.sum(), 1.0)
            mus.append(float((mu_by_type * by_type).sum() / total))
            ferms.append(group["mean_realized_fermentation_fraction"][0])
            oxy.append(group["chem/mean_oxygen"][0])

        flux = summary[steps[-1]]["nutrient_flux"]
        names = species_names(flux)
        i_o2 = names.index("oxygen")
        i_c = names.index("carbon")
        clip = flux["reaction_clip_cumulative"][:]
        boundary = flux["boundary_cumulative"][:]

    times = np.asarray(times)
    counts = np.asarray(counts, dtype=float)
    mus = np.asarray(mus)
    ferms = np.asarray(ferms)
    dt = np.diff(np.concatenate([[0.0], times]))
    per_cell = Q_MAINTENANCE + Q_CONSUMPTION * np.maximum(mus, 0.0) * (1.0 - ferms)
    demand = float((per_cell * counts * dt).sum())

    return {
        "n_final": counts[-1],
        "clip_oxygen": float(clip[i_o2]),
        "clip_carbon": float(clip[i_c]),
        "boundary_oxygen": float(boundary[i_o2]),
        "demand_oxygen": demand,
        "oxygen_first": float(oxy[0]),
        "oxygen_last": float(oxy[-1]),
    }


def main(root):
    arms = sorted(p for p in root.glob("*/output.h5"))
    if not arms:
        raise SystemExit(f"no */output.h5 under {root}")
    header = (
        f"{'arm':10s} {'N':>6s} {'clip_O2':>10s} {'demand_O2':>10s} "
        f"{'clip/dem':>8s} {'bdry_O2':>10s} {'clip/bdry':>10s} "
        f"{'clip_C':>8s} {'O2 first->last':>26s}"
    )
    print(header)
    for arm in arms:
        row = audit(arm)
        print(
            f"{arm.parent.name:10s} {row['n_final']:6.0f} "
            f"{row['clip_oxygen']:10.3e} {row['demand_oxygen']:10.3e} "
            f"{row['clip_oxygen'] / row['demand_oxygen']:8.2f} "
            f"{row['boundary_oxygen']:10.3e} "
            f"{row['clip_oxygen'] / row['boundary_oxygen']:10.3e} "
            f"{row['clip_carbon']:8.1e} "
            f"{row['oxygen_first']:12.4e} -> {row['oxygen_last']:10.4e}"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "run_root",
        type=pathlib.Path,
        help="directory containing one subdirectory per arm, each with output.h5",
    )
    main(parser.parse_args().run_root.expanduser().resolve())
