#!/usr/bin/env python3
"""Report carrying capacity against available carbon for the supply ladder.

Reads each arm's ``output.h5`` and reports, per arm:

* ``N_plateau``  -- mean population over the final quarter of the run
* ``slope``      -- (final quarter mean) / (previous quarter mean); ~1.0 is a
                    plateau, <1 still declining, >1 still climbing, so a
                    capacity number from a non-plateaued arm is not a capacity
* ``block``      -- carbon nutrient blocking fraction in the final window: the
                    share of realized carbon removal taken by agents rather
                    than by the flora. This is the observable the campaign
                    should report capacity against, not any input constant.
* ``funded``     -- realized growth uptake / uptake demand in the final window
* ``short``      -- maintenance shortfall rate (unpaid maintenance)
* ``clip``       -- cumulative reaction clip; nonzero means carbon vanished
* ``cause``      -- authoritative /run_provenance termination cause

Requires h5py: use the project interpreter (python/.venv311/bin/python), not
the system python3.
"""

import argparse
import json
import pathlib

import h5py
import numpy as np

SPECIES = "carbon"


def _decode(raw):
    value = raw[0] if not isinstance(raw, bytes) else raw
    return value.decode() if isinstance(value, bytes) else str(value)


def _species_index(group) -> int:
    names = [bytes(row[row != 0]).decode()
             for row in group["nutrient_flux/species_names"][:]]
    return names.index(SPECIES)


def _steps(handle):
    return sorted(handle["summary"].keys(),
                  key=lambda key: int(key.rsplit("_", 1)[-1]))


def read_arm(path: pathlib.Path) -> dict:
    with h5py.File(path, "r") as handle:
        steps = _steps(handle)
        times = np.array([float(handle["summary"][s]["time"][0])
                          for s in steps])
        pops = np.array([float(handle["summary"][s]["n_total"][0])
                         for s in steps])
        last = handle["summary"][steps[-1]]
        index = _species_index(last)

        def flux(key: str) -> float:
            return float(last[f"nutrient_flux/{key}"][index])

        config = json.loads(_decode(handle["run_provenance/resolved_config"][()]))
        # /run_provenance/resolved_config does not emit the epithelial carbon
        # boundary keys, so an output file cannot state the supply it was run
        # at; fall back to the arm's own input.json and report the gap.
        supply = config.get("carbon.epithelial_flux")
        if supply is None:
            arm_config = json.loads((path.parent / "input.json").read_text())
            supply = arm_config["carbon.epithelial_flux"]
        quarter = max(len(pops) // 4, 1)
        late = float(np.mean(pops[-quarter:]))
        prev = float(np.mean(pops[-2 * quarter:-quarter]))
        demand = flux("uptake_demand_interval")
        growth = flux("agent_uptake_interval")
        return {
            "supply": supply,
            "flora": config["vbf_carbon_sink_vmax"],
            "n_plateau": late,
            "slope": late / prev if prev > 0 else float("nan"),
            "n_end": pops[-1],
            "n_max": pops.max(),
            "block": flux("nutrient_blocking_fraction"),
            "funded": growth / demand if demand > 0 else float("nan"),
            "short": flux("maintenance_shortfall_interval"),
            "vbf": flux("vbf_sink_interval"),
            "clip": flux("reaction_clip_cumulative"),
            "cause": _decode(handle["run_provenance/termination_cause"][()]),
            "t_end": times[-1],
            "times": times,
            "pops": pops,
        }


def main(run_root: pathlib.Path) -> None:
    arms = sorted(p.parent for p in run_root.glob("*/output.h5"))
    if not arms:
        raise SystemExit(f"no arms with output.h5 under {run_root}")

    rows = [(arm.name, read_arm(arm / "output.h5")) for arm in arms]

    print("== capacity against available carbon ==")
    header = (f"{'arm':<12} {'flux':>10} {'vbf_vmax':>10} {'N_plateau':>10}"
              f" {'slope':>7} {'N_end':>8} {'N_max':>8} {'block':>8}"
              f" {'clip':>10} {'cause'}")
    print(header)
    for name, row in rows:
        print(f"{name:<12} {row['supply']:>10.3e} {row['flora']:>10.3e}"
              f" {row['n_plateau']:>10.1f} {row['slope']:>7.3f}"
              f" {row['n_end']:>8.0f} {row['n_max']:>8.0f}"
              f" {row['block']:>8.4f} {row['clip']:>10.2e} {row['cause']}")

    print("\n== funding in the final window (mol/s where dimensional) ==")
    print(f"{'arm':<12} {'funded':>8} {'shortfall':>12} {'vbf_sink':>12}"
          f" {'t_end_h':>12}")
    for name, row in rows:
        print(f"{name:<12} {row['funded']:>8.4f} {row['short']:>12.3e}"
              f" {row['vbf']:>12.3e} {row['t_end'] / 3600.0:>12.2f}")

    print("\n== population trajectory (t_h, N, quarter means) ==")
    for name, row in rows:
        stride = max(len(row["times"]) // 8, 1)
        trace = " ".join(
            f"{t / 3600.0:.1f}:{int(n)}"
            for t, n in zip(row["times"][::stride], row["pops"][::stride])
        )
        print(f"{name}: {trace}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=pathlib.Path, required=True,
                        help="directory containing <arm>/output.h5")
    args = parser.parse_args()
    main(args.run_root)
