"""Delivery-flux bracket analysis: which epithelial flux holds a population
under the dysbiosis guard for a full 168 h without decaying to extinction.

Reads the 12 bracket runs (6 flux multipliers x 2 corrinoid regimes, seed 1001)
and reports per arm:

  * whether the guard fired, and at what simulated time and density;
  * the density trajectory: peak, and the mean over the final 24 h of whatever
    trajectory exists, with its margin below the 1e8 cells/mL guard threshold;
  * survival: final n_total and per-type counts, and whether type 1 or type 2
    went extinct (the 0.10x anchor decayed to 16 agents, so "did not trip the
    guard" is not by itself a pass);
  * the demographic balance in the final 24 h - divisions, colicin deaths,
    lysis deaths and outflow per live agent per hour - because a flat
    population held by high turnover is a different regime from a flat
    population at low turnover, and only the second is a usable baseline;
  * mean realized mu by type, skipping intervals where a lineage is absent
    (mean_mu_by_type reads exactly 0 for an absent type, which silently
    deflates any average taken over it);
  * carbon accounting: epithelial supply, VBF sink, agent uptake and the
    demanded-vs-funded uptake ratio, to confirm the delivery rate the arm was
    configured with is the delivery rate it actually received.

Comparisons across arms are made at matched simulated time (the shortest
trajectory in the set), since guard-halted arms stop early and a final-state
comparison would confound delivery rate with censoring.

Expected layout under the run root:
```
bracket/out/<arm>/output.h5
```
The analysis summary is written to ``bracket/bracket_analysis.json``.
"""
import argparse
import glob
import json
import os
import pathlib

import h5py

J_DIR = 1.0756e-8
MULTIPLIERS = [0.10, 0.14, 0.18, 0.22, 0.26, 0.30]
KD_TAGS = [("kd1e6", 1.0e-6), ("kd1e4", 1.0e-4)]
SEED = 1001
GUARD = 1.0e8  # cells/mL
HORIZON_S = 604800.0
WINDOW_S = 86400.0  # final-24 h window


def arm_stem(mult, kd_tag):
    return f"bracket_f{round(mult * 100):03d}_{kd_tag}_s{SEED}"


def snapshots(h5):
    """Every summary snapshot, ordered by step."""
    out = []
    for key in h5["summary"]:
        g = h5[f"summary/{key}"]
        n_by_type = g["n_by_type"][()].tolist()
        ev = g["events"]
        nf = g["nutrient_flux"]
        species = [
            s.tobytes().split(b"\0")[0].decode()
            for s in nf["species_names"][()]
        ]
        ci = species.index("carbon")
        out.append({
            "step": int(g["step"][()][0]),
            "time": float(g["time"][()][0]),
            "n_total": int(g["n_total"][()][0]),
            "n1": n_by_type[1],
            "n2": n_by_type[2],
            "mu1": float(g["mean_mu_by_type"][()][1]),
            "mu2": float(g["mean_mu_by_type"][()][2]),
            "divisions": int(ev["cumulative_divisions"][()][0]),
            "colicin": int(ev["cumulative_mortality_colicin"][()][0]),
            "lysis": int(ev["cumulative_mortality_lysis"][()][0]),
            "outflow": int(ev["cumulative_outflow_boundary"][()][0])
            + int(ev["cumulative_outflow_washout"][()][0]),
            "bacteriostatic": int(
                g["stocks/bacteriostatic_live_agents"][()][0]),
            "carbon_boundary": float(nf["boundary_cumulative"][()][ci]),
            "carbon_vbf_sink": float(nf["vbf_sink_cumulative"][()][ci]),
            "carbon_uptake": float(nf["agent_uptake_cumulative"][()][ci]),
            "carbon_demand": float(nf["uptake_demand_cumulative"][()][ci]),
        })
    out.sort(key=lambda s: s["step"])
    return out


def provenance(h5):
    p = h5["run_provenance"]
    return {
        "halt_code": int(p["halt_reason_code"][()][0]),
        "termination_time": float(p["termination_time"][()][0]),
        "halt_density": float(p["halt_density_cells_per_mL"][()][0]),
        "git_sha": p["git_sha"][()],
        "image_digest": (p["container_image_digest"][()]
                         if "container_image_digest" in p else None),
    }


def density_scale(h5):
    """cells/mL per agent, from the run's own resolved config.

    Derived rather than assumed: halt_density_cells_per_mL is cross-checked
    against it below, so a domain-size or unit change cannot pass silently.
    """
    cfg = json.loads(h5["run_provenance/resolved_config"][()])
    vol_m3 = float(cfg["domain_x"]) * float(cfg["domain_y"]) \
        * float(cfg["domain_z"])
    return 1.0 / (vol_m3 * 1.0e6)


def window_rates(snaps, scale, end_time, window):
    """Per-live-agent per-hour demographic rates over (end-window, end]."""
    inside = [s for s in snaps if s["time"] > end_time - window
              and s["time"] <= end_time]
    if len(inside) < 2:
        return None
    first, last = inside[0], inside[-1]
    hours = (last["time"] - first["time"]) / 3600.0
    live = [s["n_total"] for s in inside if s["n_total"] > 0]
    mean_live = sum(live) / len(live) if live else 0.0
    if hours <= 0 or mean_live <= 0:
        return None

    def rate(key):
        return (last[key] - first[key]) / hours / mean_live

    # mu averaged only over intervals where the lineage is present: an absent
    # type reports exactly 0 and would otherwise be averaged in as real data.
    def mu(key, nkey):
        vals = [s[key] for s in inside if s[nkey] > 0]
        return sum(vals) / len(vals) if vals else None

    dens = [s["n_total"] * scale for s in inside]
    return {
        "hours": hours,
        "mean_live": mean_live,
        "mean_density": sum(dens) / len(dens),
        "max_density": max(dens),
        "div_per_agent_h": rate("divisions"),
        "colicin_per_agent_h": rate("colicin"),
        "lysis_per_agent_h": rate("lysis"),
        "outflow_per_agent_h": rate("outflow"),
        "mu1": mu("mu1", "n1"),
        "mu2": mu("mu2", "n2"),
        "bacteriostatic_frac": sum(
            s["bacteriostatic"] for s in inside) / max(1.0, sum(
                s["n_total"] for s in inside)),
    }


def load_arm(root, stem):
    path = os.path.join(root, stem, "output.h5")
    if not os.path.exists(path):
        hits = glob.glob(os.path.join(root, "**", f"*{stem}*", "output.h5"),
                         recursive=True)
        if not hits:
            return None
        path = hits[0]
    with h5py.File(path, "r") as h5:
        snaps = snapshots(h5)
        prov = provenance(h5)
        scale = density_scale(h5)
    if not snaps:
        return None
    peak = max(s["n_total"] for s in snaps) * scale
    final = snaps[-1]
    return {
        "stem": stem, "path": path, "prov": prov, "scale": scale,
        "snaps": snaps, "peak_density": peak,
        "final": final,
        "reached_h": final["time"] / 3600.0,
        "guard_fired": prov["halt_code"] != 0,
    }


def main(run_root):
    root = os.path.join(run_root, "bracket", "out")
    arms = {}
    for mult in MULTIPLIERS:
        for kd_tag, _kd in KD_TAGS:
            stem = arm_stem(mult, kd_tag)
            arm = load_arm(root, stem)
            if arm is None:
                print(f"MISSING: {stem}")
                continue
            arm["mult"] = mult
            arm["kd_tag"] = kd_tag
            arms[stem] = arm

    if not arms:
        raise SystemExit("no arms loaded")

    # Consistency check: the guard's own reported density must agree with the
    # scale derived from the resolved config.
    for arm in arms.values():
        if arm["guard_fired"] and arm["prov"]["halt_density"] > 0:
            implied = arm["prov"]["halt_density"] / arm["scale"]
            print(f"{arm['stem']}: halt_density "
                  f"{arm['prov']['halt_density']:.3e} cells/mL implies "
                  f"{implied:.1f} agents (guard threshold {GUARD:.0e} "
                  f"= {GUARD / arm['scale']:.1f} agents)")

    matched = min(a["final"]["time"] for a in arms.values())
    print(f"\nmatched simulated time: {matched / 3600.0:.2f} h "
          f"(full horizon = {HORIZON_S / 3600.0:.0f} h)")

    rows = []
    for stem, arm in arms.items():
        full = window_rates(arm["snaps"], arm["scale"],
                            arm["final"]["time"], WINDOW_S)
        at_matched = window_rates(arm["snaps"], arm["scale"], matched,
                                  WINDOW_S)
        fin = arm["final"]
        carbon_supply = fin["carbon_boundary"]
        row = {
            "stem": stem, "mult": arm["mult"], "kd": arm["kd_tag"],
            "reached_h": arm["reached_h"],
            "guard_fired": arm["guard_fired"],
            "halt_code": arm["prov"]["halt_code"],
            "peak_density": arm["peak_density"],
            "final_n": fin["n_total"], "final_n1": fin["n1"],
            "final_n2": fin["n2"],
            "extinct_n1": fin["n1"] == 0, "extinct_n2": fin["n2"] == 0,
            "survived_full_horizon": (not arm["guard_fired"])
            and fin["time"] >= HORIZON_S - 1.0 and fin["n_total"] > 0,
            "carbon_supply_mol": carbon_supply,
            "carbon_vbf_frac": (fin["carbon_vbf_sink"] / carbon_supply
                                if carbon_supply else None),
            "carbon_agent_frac": (fin["carbon_uptake"] / carbon_supply
                                  if carbon_supply else None),
            "uptake_funded_frac": (fin["carbon_uptake"] / fin["carbon_demand"]
                                   if fin["carbon_demand"] else None),
            "final_window": full,
            "matched_window": at_matched,
        }
        rows.append(row)

    rows.sort(key=lambda r: (r["kd"], r["mult"]))

    print("\n=== at matched time "
          f"({matched / 3600.0:.1f} h), final 24 h window ===")
    hdr = ("kd      flux    reach_h  guard  peak_dens   mean_dens   "
           "n(1/2)      div/a/h  lys/a/h  col/a/h  out/a/h  mu1")
    print(hdr)
    for r in rows:
        w = r["matched_window"]
        if w is None:
            print(f"{r['kd']:6s}  {r['mult']:.2f}x  "
                  f"{r['reached_h']:7.1f}  (too short for a 24 h window)")
            continue
        mu1 = "n/a" if w["mu1"] is None else f"{w['mu1']:.2e}"
        print(f"{r['kd']:6s}  {r['mult']:.2f}x  {r['reached_h']:7.1f}  "
              f"{'YES' if r['guard_fired'] else ' no':5s}  "
              f"{r['peak_density']:.3e}  {w['mean_density']:.3e}  "
              f"{r['final_n1']:4d}/{r['final_n2']:<4d}  "
              f"{w['div_per_agent_h']:7.3f}  {w['lysis_per_agent_h']:7.3f}  "
              f"{w['colicin_per_agent_h']:7.3f}  "
              f"{w['outflow_per_agent_h']:7.3f}  {mu1}")

    print("\n=== carbon accounting (cumulative, end of each run) ===")
    print("kd      flux    supply_mol   VBF/supply  agents/supply  "
          "funded/demanded")
    for r in rows:
        vbf = r["carbon_vbf_frac"]
        ag = r["carbon_agent_frac"]
        fu = r["uptake_funded_frac"]
        print(f"{r['kd']:6s}  {r['mult']:.2f}x  "
              f"{r['carbon_supply_mol']:.4e}  "
              f"{'n/a' if vbf is None else f'{vbf:10.3f}'}  "
              f"{'n/a' if ag is None else f'{ag:13.5f}'}  "
              f"{'n/a' if fu is None else f'{fu:15.4f}'}")

    print("\n=== verdict per arm ===")
    for r in rows:
        if r["guard_fired"]:
            verdict = (f"guard at {r['reached_h']:.1f} h - "
                       "too much delivery for a 168 h horizon")
        elif r["final_n"] == 0:
            verdict = "extinct - too little delivery"
        elif r["extinct_n1"] or r["extinct_n2"]:
            lost = "producers" if r["extinct_n1"] else "susceptibles"
            verdict = (f"survived the horizon but lost the {lost} - "
                       "not a coexistence baseline")
        elif not r["survived_full_horizon"]:
            verdict = (f"stopped at {r['reached_h']:.1f} h without the guard "
                       "- inspect (timeout or application failure)")
        else:
            verdict = (f"held {r['final_n']} agents for 168 h at "
                       f"{r['peak_density'] / GUARD:.2f}x the guard threshold")
        print(f"  {r['kd']:6s} {r['mult']:.2f}x: {verdict}")

    out = os.path.join(run_root, "bracket", "bracket_analysis.json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump({"matched_time_s": matched, "guard_cells_per_mL": GUARD,
                   "arms": rows}, fh, indent=1, default=str)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--run-root",
        default=os.environ.get("GUTIBM_CAMPAIGN_ROOT", "."),
        help="campaign root containing bracket/out/ (default: "
        "GUTIBM_CAMPAIGN_ROOT or current directory)",
    )
    main(pathlib.Path(parser.parse_args().run_root).expanduser().resolve())
