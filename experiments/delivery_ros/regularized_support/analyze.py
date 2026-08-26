#!/usr/bin/env python3
"""Acceptance analysis for PR #335: population-scale resolution invariance."""
import json
from pathlib import Path

import h5py
import numpy as np

ROOT = Path("/home/ubuntu/gutibm-campaign/regularized-support")
ARMS = ["P_ff10_res2", "P_ff10_res4", "P_ff10_res6",
        "P_ff0_res2", "P_ff0_res4", "P_ff0_res6"]


def sk(group):
    return sorted(group.keys(), key=lambda v: int(v.rsplit("_", 1)[1]))


def s(dataset):
    return float(np.asarray(dataset[()]).reshape(-1)[0])


def i(dataset):
    return int(np.asarray(dataset[()]).reshape(-1)[0])


def species(flux):
    return [bytes(r[r != 0]).decode() for r in flux["species_names"][:]]


def f(flux, name, index):
    return float(flux[name][index]) if name in flux else 0.0


def analyze(arm):
    with h5py.File(ROOT / arm / "output.h5", "r") as fh:
        steps = sk(fh["summary"])
        rec = fh["summary"][steps[-1]]
        flux = rec["nutrient_flux"]
        names = species(flux)
        ox, ca = names.index("oxygen"), names.index("carbon")
        ev = rec["events"]
        # delivery_* datasets are per-interval; accumulate over all summaries
        red = {"carbon": 0.0, "oxygen": 0.0}
        ret = {"carbon": 0.0, "oxygen": 0.0}
        for key in steps:
            fl = fh["summary"][key]["nutrient_flux"]
            for sp, idx in (("carbon", ca), ("oxygen", ox)):
                red[sp] = max(red[sp], f(fl, "delivery_reduction_cumulative", idx))
                ret[sp] = max(ret[sp], f(fl, "delivery_retry_events_cumulative", idx))
        prov = fh["run_provenance"]
        cause = prov["termination_cause"][()]
        if isinstance(cause, bytes):
            cause = cause.decode()
        else:
            cause = str(np.asarray(cause).reshape(-1)[0])
        out = {
            "arm": arm,
            "grid_dx": json.loads((ROOT / arm / "input.json").read_text())["grid_dx"],
            "radius": json.loads((ROOT / arm / "input.json").read_text())[
                "metabolism.delivery_far_field_radius"],
            "final_N": i(rec["n_total"]),
            "peak_N": max(i(fh["summary"][k]["n_total"]) for k in steps),
            "divisions": i(ev["cumulative_divisions"]) if "cumulative_divisions" in ev else 0,
            "mortality_lysis": i(ev["cumulative_mortality_lysis"]) if "cumulative_mortality_lysis" in ev else 0,
            "outflow_boundary": i(ev["cumulative_outflow_boundary"]) if "cumulative_outflow_boundary" in ev else 0,
            "funded_oxygen": f(flux, "agent_uptake_cumulative", ox),
            "demanded_oxygen": f(flux, "uptake_demand_cumulative", ox),
            "funded_carbon": f(flux, "agent_uptake_cumulative", ca),
            "demanded_carbon": f(flux, "uptake_demand_cumulative", ca),
            "maint_shortfall_oxygen": f(flux, "maintenance_shortfall_cumulative", ox),
            "maint_shortfall_carbon": f(flux, "maintenance_shortfall_cumulative", ca),
            "clips_oxygen": f(flux, "reaction_clip_cumulative", ox),
            "clips_carbon": f(flux, "reaction_clip_cumulative", ca),
            "reduction_oxygen": red["oxygen"],
            "reduction_carbon": red["carbon"],
            "retries_oxygen": ret["oxygen"],
            "retries_carbon": ret["carbon"],
            "mean_fermentation_fraction": s(rec["mean_realized_fermentation_fraction"]),
            "termination_cause": cause,
        }
        out["funded_fraction_oxygen"] = (
            out["funded_oxygen"] / out["demanded_oxygen"] if out["demanded_oxygen"] else float("nan"))
        out["funded_fraction_carbon"] = (
            out["funded_carbon"] / out["demanded_carbon"] if out["demanded_carbon"] else float("nan"))
        return out


def spread(rows, key):
    vals = [r[key] for r in rows]
    lo, hi = min(vals), max(vals)
    return hi / lo if lo > 0 else float("inf")


rows = [analyze(a) for a in ARMS]
(ROOT / "metrics.json").write_text(json.dumps(rows, indent=1))

hdr = ["arm", "grid_dx", "radius", "final_N", "peak_N", "divisions",
       "funded_fraction_oxygen", "funded_fraction_carbon",
       "reduction_oxygen", "reduction_carbon",
       "retries_oxygen", "retries_carbon",
       "clips_oxygen", "clips_carbon", "termination_cause"]
print("\t".join(hdr))
for r in rows:
    print("\t".join(
        f"{r[k]:.4g}" if isinstance(r[k], float) else str(r[k]) for k in hdr))

on = [r for r in rows if r["radius"] > 0]
off = [r for r in rows if r["radius"] == 0]
print()
for label, group in (("radius 10um", on), ("radius 0", off)):
    print(f"{label}: final_N ratio max/min = {spread(group, 'final_N'):.3f}"
          f", funded_fraction_oxygen ratio = {spread(group, 'funded_fraction_oxygen'):.3f}"
          f", funded_fraction_carbon ratio = {spread(group, 'funded_fraction_carbon'):.3f}")
