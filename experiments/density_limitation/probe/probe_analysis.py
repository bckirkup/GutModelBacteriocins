"""Analyse the #302/#304 mechanism probe. Lead-authored.

Reads the six arm outputs and reports, per arm:
  - population trajectory, final N, guard halt time and halt density
  - mean realized fermentation fraction (is the aerobic branch masked?)
  - carbon: epithelial delivery, agent growth uptake, maintenance draw,
    maintenance shortfall, VBF sink - all as interval rates in the last window
    and as cumulative totals
  - realized carbon cost, measured as carbon uptake per unit biomass gained
  - acetate accumulation, oxygen and carbon mean concentration
  - division / lysis / outflow per agent per hour in the last 6 h

Expected layout under the run root:
```
probe/<arm>/output.h5
```
The six arm directories are named by ``ARMS`` below.
"""

import argparse
import os
import pathlib
import sys

import h5py
import numpy as np

ARMS = ["off_f018", "maint_f018", "full_f018",
        "full_f030", "full_f060", "full_f100"]
DOMAIN_VOLUME_M3 = 9.6e-5 * 9.6e-5 * 3.0e-4


def species_names(grp):
    raw = grp["nutrient_flux/species_names"][:]
    out = []
    for row in raw:
        out.append(bytes(row[row != 0]).decode())
    return out


def steps(f):
    return sorted(f["summary"].keys())


def scalar(grp, key):
    return float(grp[key][0])


def arm_report(root, name):
    path = root / name / "output.h5"
    if not path.exists():
        print(f"{name}: MISSING {path}")
        return None
    with h5py.File(path, "r") as f:
        sk = steps(f)
        last = f["summary"][sk[-1]]
        names = species_names(last)
        ic = names.index("carbon")
        ia = names.index("acetate") if "acetate" in names else None

        n_t = np.array([scalar(f["summary"][s], "n_total") for s in sk])
        t_h = np.array([scalar(f["summary"][s], "time") for s in sk]) / 3600.0
        ferm = np.array([scalar(f["summary"][s],
                                "mean_realized_fermentation_fraction")
                         for s in sk])
        halt = int(scalar(last, "halt_reason_code"))
        halt_dens = scalar(last, "halt_density_cells_per_mL")

        cum = last["nutrient_flux/boundary_cumulative"][ic]
        up = last["nutrient_flux/agent_uptake_cumulative"][ic]
        maint = last["nutrient_flux/maintenance_cumulative"][ic]
        short = last["nutrient_flux/maintenance_shortfall_cumulative"][ic]
        vbf = last["nutrient_flux/vbf_sink_cumulative"][ic]

        # interval rates in the final window
        dt_int = (scalar(last, "nutrient_flux/interval_end_time")
                  - scalar(last, "nutrient_flux/interval_start_time"))
        i_bnd = last["nutrient_flux/boundary_interval"][ic] / dt_int
        i_up = last["nutrient_flux/agent_uptake_interval"][ic] / dt_int
        i_maint = last["nutrient_flux/maintenance_interval"][ic] / dt_int
        i_short = last["nutrient_flux/maintenance_shortfall_interval"][ic] / dt_int
        i_vbf = last["nutrient_flux/vbf_sink_interval"][ic] / dt_int

        # realized carbon cost = cumulative growth uptake / biomass gained
        ag = sorted(f["agents"].keys())
        bm0 = float(np.sum(f["agents"][ag[0]]["biomass"][:]))
        bm1 = float(np.sum(f["agents"][ag[-1]]["biomass"][:]))
        divisions = scalar(last, "events/cumulative_divisions")
        # biomass gained includes biomass lost with dead/outflowed agents, so
        # this is a lower bound on gross production
        cost = up / max(bm1 - bm0, 1e-30)

        acet = (last["nutrient_flux/vbf_source_cumulative"][ia]
                if ia is not None else float("nan"))
        mean_c = scalar(last, "chem/mean_carbon")
        mean_o2 = scalar(last, "chem/mean_oxygen")

        # last 6 h event rates per agent per hour
        win = t_h >= t_h[-1] - 6.0
        idx = np.where(win)[0]
        a, b = idx[0], idx[-1]
        span = max(t_h[b] - t_h[a], 1e-9)
        mean_n = max(float(np.mean(n_t[a:b + 1])), 1e-9)

        def rate(key):
            d = (scalar(f["summary"][sk[b]], key)
                 - scalar(f["summary"][sk[a]], key))
            return d / span / mean_n

        div_r = rate("events/cumulative_divisions")
        lys_r = rate("events/cumulative_mortality_lysis")
        out_r = (rate("events/cumulative_outflow_boundary")
                 + rate("events/cumulative_outflow_washout"))

        stat = {
            "arm": name, "t_end_h": t_h[-1], "n_end": int(n_t[-1]),
            "n_max": int(n_t.max()), "halt": halt, "halt_dens": halt_dens,
            "dens_end": n_t[-1] / DOMAIN_VOLUME_M3 / 1e6,
            "ferm_end": ferm[-1], "ferm_max": ferm.max(),
            "bnd_cum": cum, "up_cum": up, "maint_cum": maint,
            "short_cum": short, "vbf_cum": vbf,
            "i_bnd": i_bnd, "i_up": i_up, "i_maint": i_maint,
            "i_short": i_short, "i_vbf": i_vbf,
            "cost": cost, "acet": acet, "mean_c": mean_c, "mean_o2": mean_o2,
            "div_r": div_r, "lys_r": lys_r, "out_r": out_r,
            "divisions": divisions, "bm0": bm0, "bm1": bm1,
            "traj": list(zip(t_h[::max(len(t_h) // 12, 1)],
                             n_t[::max(len(t_h) // 12, 1)])),
        }
        return stat


def main(run_root: pathlib.Path) -> None:
    root = run_root / "probe"
    rows = [r for r in (arm_report(root, a) for a in ARMS) if r]
    if not rows:
        sys.exit("no arms")

    print("\n== population and mode ==")
    print(f"{'arm':<12}{'t_end':>7}{'N_end':>7}{'N_max':>7}{'halt':>6}"
          f"{'dens_end':>11}{'ferm_end':>10}{'ferm_max':>10}")
    for r in rows:
        print(f"{r['arm']:<12}{r['t_end_h']:>7.1f}{r['n_end']:>7d}"
              f"{r['n_max']:>7d}{r['halt']:>6d}{r['dens_end']:>11.3e}"
              f"{r['ferm_end']:>10.4f}{r['ferm_max']:>10.4f}")

    print("\n== carbon budget, final-window rates (mol/s) ==")
    print(f"{'arm':<12}{'delivery':>11}{'growth':>11}{'maint':>11}"
          f"{'shortfall':>11}{'vbf_sink':>11}{'maint/deliv':>12}")
    for r in rows:
        ratio = r["i_maint"] / r["i_bnd"] if r["i_bnd"] else float("nan")
        print(f"{r['arm']:<12}{r['i_bnd']:>11.3e}{r['i_up']:>11.3e}"
              f"{r['i_maint']:>11.3e}{r['i_short']:>11.3e}"
              f"{r['i_vbf']:>11.3e}{ratio:>12.3f}")

    print("\n== carbon cost, cumulative (mol) and per kg biomass ==")
    print(f"{'arm':<12}{'delivered':>12}{'growth':>12}{'maint':>12}"
          f"{'short':>12}{'cost_mol/kg':>13}{'dbiomass_kg':>13}")
    for r in rows:
        print(f"{r['arm']:<12}{r['bnd_cum']:>12.3e}{r['up_cum']:>12.3e}"
              f"{r['maint_cum']:>12.3e}{r['short_cum']:>12.3e}"
              f"{r['cost']:>13.3e}{r['bm1'] - r['bm0']:>13.3e}")

    print("\n== environment and demography ==")
    print(f"{'arm':<12}{'mean_C':>11}{'mean_O2':>11}{'acetate_src':>13}"
          f"{'div/a/h':>10}{'lys/a/h':>10}{'out/a/h':>10}")
    for r in rows:
        print(f"{r['arm']:<12}{r['mean_c']:>11.3e}{r['mean_o2']:>11.3e}"
              f"{r['acet']:>13.3e}{r['div_r']:>10.4f}{r['lys_r']:>10.4f}"
              f"{r['out_r']:>10.4f}")

    print("\n== trajectories (h, N) ==")
    for r in rows:
        pts = " ".join(f"{t:.0f}h:{int(n)}" for t, n in r["traj"])
        print(f"{r['arm']:<12}{pts}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--run-root",
        default=os.environ.get("GUTIBM_CAMPAIGN_ROOT", "."),
        help="campaign root containing probe/<arm>/output.h5 (default: "
        "GUTIBM_CAMPAIGN_ROOT or current directory)",
    )
    main(pathlib.Path(parser.parse_args().run_root).expanduser().resolve())
