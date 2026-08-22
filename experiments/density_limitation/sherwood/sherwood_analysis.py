"""Sherwood-vs-none analysis: does the diffusive uptake cap create a carrying
capacity below the dysbiosis guard, and by how much does it diverge from `none`?

Reuses the bracket loaders so the two experiments are read identically. Adds the
quantities that only matter once uptake can be capped:

  * `funded / demanded` carbon per arm and over time - exactly 1.0 means the cap
    never bound and the arm is behaviourally identical to `none`;
  * the fraction of agent-steps that were uptake-limited;
  * whether density plateaus: the final-24 h mean density and its trend, versus
    the guard.

Divergence is reported per matched flux as a paired none-vs-sherwood
comparison, since the arms share seed, flux and image.

Expected layout under the run root:
```
sherwood/out/<arm>/output.h5
```
The analysis summary is written to ``sherwood/sherwood_analysis.json``.
"""

import argparse
import importlib.util
import json
import os
import pathlib

import h5py

MULTIPLIERS = [0.14, 0.18, 0.22]
LIMITS = ["none", "sherwood"]
GUARD = 1.0e8
HORIZON_S = 604800.0
WINDOW_S = 86400.0

def carbon_index(nf):
    species = [s.tobytes().split(b"\0")[0].decode()
               for s in nf["species_names"][()]]
    return species.index("carbon")


def uptake_series(h5):
    """Per-snapshot uptake-limit telemetry, ordered by step."""
    out = []
    for key in h5["summary"]:
        g = h5[f"summary/{key}"]
        nf = g["nutrient_flux"]
        ci = carbon_index(nf)
        out.append({
            "step": int(g["step"][()][0]),
            "time": float(g["time"][()][0]),
            "n_total": int(g["n_total"][()][0]),
            "limited_interval": float(
                nf["uptake_limited_agents_interval"][()][ci]),
            "limited_cumulative": float(
                nf["uptake_limited_agents_cumulative"][()][ci]),
            "uptake_interval": float(nf["agent_uptake_interval"][()][ci]),
            "demand_interval": float(nf["uptake_demand_interval"][()][ci]),
        })
    out.sort(key=lambda s: s["step"])
    return out


def load(root, analysis, mult, limit):
    stem = f"sherwood_f{round(mult * 100):03d}_{limit}_s1001"
    path = os.path.join(root, stem, "output.h5")
    if not os.path.exists(path):
        return None
    with h5py.File(path, "r") as h5:
        snaps = analysis.snapshots(h5)
        prov = analysis.provenance(h5)
        scale = analysis.density_scale(h5)
        ups = uptake_series(h5)
    fin = snaps[-1]
    return {
        "stem": stem, "mult": mult, "limit": limit, "snaps": snaps,
        "prov": prov, "scale": scale, "uptake": ups,
        "final": fin,
        "reached_h": fin["time"] / 3600.0,
        "guard_fired": prov["halt_code"] != 0,
        "peak_density": max(s["n_total"] for s in snaps) * scale,
        "window": analysis.window_rates(snaps, scale, fin["time"], WINDOW_S),
    }


def main(run_root):
    root = os.path.join(run_root, "sherwood", "out")
    bracket_path = pathlib.Path(__file__).parents[1] / "bracket" / \
        "bracket_analysis.py"
    spec = importlib.util.spec_from_file_location(
        "bracket_analysis", bracket_path)
    ba = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ba)
    arms = {}
    for mult in MULTIPLIERS:
        for limit in LIMITS:
            arm = load(root, ba, mult, limit)
            if arm is None:
                print(f"MISSING: f{round(mult * 100):03d}/{limit}")
                continue
            arms[(mult, limit)] = arm

    if not arms:
        raise SystemExit("no arms loaded")

    matched = min(a["final"]["time"] for a in arms.values())
    print(f"matched simulated time: {matched / 3600.0:.2f} h "
          f"(horizon {HORIZON_S / 3600.0:.0f} h)\n")

    print("=== outcome per arm ===")
    print("flux   limit      reach_h  guard  peak_dens   "
          "mean_dens(24h)  final_n  n1/n2     funded/demanded  limited_frac")
    rows = []
    for (mult, limit), a in sorted(arms.items()):
        fin = a["final"]
        funded = (fin["carbon_uptake"] / fin["carbon_demand"]
                  if fin["carbon_demand"] else None)
        # agent-steps that were uptake-limited, over agent-steps simulated
        agent_steps = sum(s["n_total"] for s in a["snaps"])
        limited_frac = (a["uptake"][-1]["limited_cumulative"] / agent_steps
                        if agent_steps else None)
        w = a["window"]
        rows.append({
            "mult": mult, "limit": limit, "reached_h": a["reached_h"],
            "guard_fired": a["guard_fired"],
            "halt_code": a["prov"]["halt_code"],
            "peak_density": a["peak_density"],
            "mean_density_24h": w["mean_density"] if w else None,
            "final_n": fin["n_total"], "n1": fin["n1"], "n2": fin["n2"],
            "funded_over_demanded": funded, "limited_frac": limited_frac,
            "div_per_agent_h": w["div_per_agent_h"] if w else None,
            "loss_per_agent_h": (w["outflow_per_agent_h"] + w["lysis_per_agent_h"]
                                 + w["colicin_per_agent_h"]) if w else None,
            "mu1": w["mu1"] if w else None,
        })
        print(f"{mult:.2f}x  {limit:9s}  {a['reached_h']:7.1f}  "
              f"{'YES' if a['guard_fired'] else ' no':5s}  "
              f"{a['peak_density']:.3e}  "
              f"{(w['mean_density'] if w else float('nan')):.3e}      "
              f"{fin['n_total']:5d}  {fin['n1']:4d}/{fin['n2']:<4d}  "
              f"{'n/a' if funded is None else f'{funded:15.6f}'}  "
              f"{'n/a' if limited_frac is None else f'{limited_frac:12.4f}'}")

    print("\n=== paired divergence at matched flux (sherwood vs none) ===")
    for mult in MULTIPLIERS:
        n = arms.get((mult, "none"))
        s = arms.get((mult, "sherwood"))
        if not n or not s:
            continue
        wn, ws = n["window"], s["window"]

        def rel(a_, b_):
            if a_ in (None, 0) or b_ is None:
                return None
            return (b_ - a_) / a_

        print(f"  {mult:.2f}x  guard: none {'YES' if n['guard_fired'] else 'no'}"
              f" at {n['reached_h']:.1f} h vs sherwood "
              f"{'YES' if s['guard_fired'] else 'no'} at {s['reached_h']:.1f} h")
        if wn and ws:
            d_mu = rel(wn["mu1"], ws["mu1"])
            d_div = rel(wn["div_per_agent_h"], ws["div_per_agent_h"])
            d_dens = rel(wn["mean_density"], ws["mean_density"])
            for label, val in (("mu1", d_mu), ("div/agent/h", d_div),
                               ("mean density", d_dens)):
                print(f"      {label:14s} "
                      f"{'n/a' if val is None else f'{val * 100:+8.2f}%'}")

    print("\n=== does sherwood bind? (per-interval, sherwood arms only) ===")
    for mult in MULTIPLIERS:
        s = arms.get((mult, "sherwood"))
        if not s:
            continue
        ups = [u for u in s["uptake"] if u["n_total"] > 0]
        binding = [u for u in ups if u["demand_interval"] > 0
                   and u["uptake_interval"] < u["demand_interval"] * (1 - 1e-9)]
        worst = min((u["uptake_interval"] / u["demand_interval"]
                     for u in binding), default=None)
        print(f"  {mult:.2f}x: {len(binding)}/{len(ups)} intervals with the cap "
              f"binding; tightest funded fraction "
              f"{'n/a' if worst is None else f'{worst:.4f}'}")
        if binding:
            first = binding[0]
            print(f"        first bound at {first['time'] / 3600.0:.1f} h "
                  f"with n={first['n_total']} "
                  f"({first['n_total'] * s['scale']:.2e} cells/mL)")

    out = os.path.join(run_root, "sherwood", "sherwood_analysis.json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump({"matched_time_s": matched, "guard_cells_per_mL": GUARD,
                   "arms": rows}, fh, indent=1, default=str)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--run-root",
        default=os.environ.get("GUTIBM_CAMPAIGN_ROOT", "."),
        help="campaign root containing sherwood/out/ (default: "
        "GUTIBM_CAMPAIGN_ROOT or current directory)",
    )
    main(pathlib.Path(parser.parse_args().run_root).expanduser().resolve())
