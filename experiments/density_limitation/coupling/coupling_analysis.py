"""Analyse the agent_carbon_coupling probe (#312 Change 1). Lead-authored.

Reads the four arm outputs and reports, per arm:
  - the observable: nutrient blocking fraction for carbon, in the final
    summary window and its trajectory
  - population trajectory, final and maximum N, and plateau evidence
  - funded fraction = realized growth uptake / demanded uptake, plus the
    maintenance shortfall, since coupling brakes by starving delivery
  - carbon budget as final-window interval rates
  - the authoritative termination record, so a starved arm is reported as
    a closure_violation result rather than read as a completed run

Expected layout under the run root:

    coupling/<arm>/output.h5

The claim under test is that blocking fraction rises monotonically with
coupling and that carrying capacity tracks blocking fraction, not the coupling
constant. Report campaign results against blocking fraction.
"""

import argparse
import os
import pathlib
import sys

import h5py
import numpy as np

ARMS = ["c0", "c1e21", "c1e20", "c1e19"]
COUPLING = {"c0": 0.0, "c1e21": 1.0e-21, "c1e20": 1.0e-20, "c1e19": 1.0e-19}
DOMAIN_VOLUME_M3 = 9.6e-5 * 9.6e-5 * 3.0e-4


def species_names(grp):
    raw = grp["nutrient_flux/species_names"][:]
    return [bytes(row[row != 0]).decode() for row in raw]


def steps(f):
    return sorted(f["summary"].keys())


def scalar(grp, key):
    return float(grp[key][0])


def text(f, key):
    if key not in f:
        return "absent"
    raw = f[key][()]
    if isinstance(raw, np.ndarray):
        raw = raw[0] if raw.size else b""
    return raw.decode() if isinstance(raw, bytes) else str(raw)


def arm_report(root, name):
    path = root / name / "output.h5"
    if not path.exists():
        print(f"{name}: MISSING {path}")
        return None
    with h5py.File(path, "r") as f:
        sk = steps(f)
        last = f["summary"][sk[-1]]
        ic = species_names(last).index("carbon")

        n_t = np.array([scalar(f["summary"][s], "n_total") for s in sk])
        t_h = np.array([scalar(f["summary"][s], "time") for s in sk]) / 3600.0
        block = np.array(
            [f["summary"][s]["nutrient_flux/nutrient_blocking_fraction"][ic]
             if "nutrient_flux/nutrient_blocking_fraction" in f["summary"][s]
             else np.nan for s in sk])

        dt_int = (scalar(last, "nutrient_flux/interval_end_time")
                  - scalar(last, "nutrient_flux/interval_start_time"))
        i_up = last["nutrient_flux/agent_uptake_interval"][ic] / dt_int
        i_maint = last["nutrient_flux/maintenance_interval"][ic] / dt_int
        i_short = (last["nutrient_flux/maintenance_shortfall_interval"][ic]
                   / dt_int)
        i_vbf = last["nutrient_flux/vbf_sink_interval"][ic] / dt_int
        i_bnd = last["nutrient_flux/boundary_interval"][ic] / dt_int
        demand_key = "nutrient_flux/uptake_demand_interval"
        i_demand = (last[demand_key][ic] / dt_int
                    if demand_key in last else float("nan"))
        clip = (last["nutrient_flux/reaction_clip_cumulative"][ic]
                if "nutrient_flux/reaction_clip_cumulative" in last
                else float("nan"))

        # plateau evidence: mean N over the last quarter vs the quarter before
        q = max(len(n_t) // 4, 1)
        n_late = float(np.mean(n_t[-q:]))
        n_prev = float(np.mean(n_t[-2 * q:-q])) if len(n_t) > q else n_late

        return {
            "arm": name, "coupling": COUPLING[name],
            "t_end_h": t_h[-1], "n_end": int(n_t[-1]), "n_max": int(n_t.max()),
            "n_late": n_late, "n_prev": n_prev,
            "dens_end": n_t[-1] / DOMAIN_VOLUME_M3 / 1e6,
            "block_end": float(block[-1]),
            "block_max": float(np.nanmax(block)),
            "i_up": i_up, "i_maint": i_maint, "i_short": i_short,
            "i_vbf": i_vbf, "i_bnd": i_bnd, "i_demand": i_demand,
            "funded": i_up / i_demand if i_demand else float("nan"),
            "clip": clip,
            "cause": text(f, "run_provenance/termination_cause"),
            "detail": text(f, "run_provenance/termination_detail"),
            "traj": list(zip(t_h[::max(len(t_h) // 8, 1)],
                             n_t[::max(len(t_h) // 8, 1)],
                             block[::max(len(t_h) // 8, 1)])),
        }


def main(run_root: pathlib.Path) -> None:
    root = run_root / "coupling"
    rows = [r for r in (arm_report(root, a) for a in ARMS) if r]
    if not rows:
        sys.exit("no arms")

    print("\n== the observable: carbon nutrient blocking fraction ==")
    print(f"{'arm':<8}{'coupling':>11}{'block_end':>11}{'block_max':>11}"
          f"{'N_end':>7}{'N_max':>7}{'N_late':>9}{'N_prev':>9}")
    for r in rows:
        print(f"{r['arm']:<8}{r['coupling']:>11.2e}{r['block_end']:>11.4f}"
              f"{r['block_max']:>11.4f}{r['n_end']:>7d}{r['n_max']:>7d}"
              f"{r['n_late']:>9.1f}{r['n_prev']:>9.1f}")

    print("\n== carbon budget, final-window rates (mol/s) ==")
    print(f"{'arm':<8}{'delivery':>11}{'growth':>11}{'demand':>11}"
          f"{'funded':>9}{'maint':>11}{'shortfall':>11}{'vbf_sink':>11}")
    for r in rows:
        print(f"{r['arm']:<8}{r['i_bnd']:>11.3e}{r['i_up']:>11.3e}"
              f"{r['i_demand']:>11.3e}{r['funded']:>9.4f}"
              f"{r['i_maint']:>11.3e}{r['i_short']:>11.3e}"
              f"{r['i_vbf']:>11.3e}")

    print("\n== termination record and closure ==")
    print(f"{'arm':<8}{'t_end_h':>9}{'clip_cum':>12}  cause / detail")
    for r in rows:
        print(f"{r['arm']:<8}{r['t_end_h']:>9.2f}{r['clip']:>12.3e}  "
              f"{r['cause']} / {r['detail']}")

    print("\n== trajectories (h, N, blocking fraction) ==")
    for r in rows:
        pts = " ".join(f"{t:.0f}h:{int(n)}/{b:.3f}" for t, n, b in r["traj"])
        print(f"{r['arm']:<8}{pts}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--run-root",
        default=os.environ.get("GUTIBM_CAMPAIGN_ROOT", "."),
        help="campaign root containing coupling/<arm>/output.h5 (default: "
        "GUTIBM_CAMPAIGN_ROOT or current directory)",
    )
    main(pathlib.Path(parser.parse_args().run_root).expanduser().resolve())
