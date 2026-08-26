#!/usr/bin/env python3
"""Q_O2_per_gen at low density, on merged main (post-#338 positivity rationing).

Why this campaign exists: the earlier Q_O2_per_gen table (7.86e-14 / 2.44e-14 /
2.38e-14 mol at 2/4/6 um) was measured at bloom density on code that allowed the
oxygen field to go negative, and #338 showed that at bloom the feasible ration is
1.5% of the analytic draw, i.e. the answer there is epithelial supply divided by
the crowd. The lysis coefficient needs the *uncrowded* per-cell quantity, so this
measures it at low founder counts and asks two questions separately:

  A. grid invariance   - Q_n2_res{2,4,6}: does Q stop tracking grid_dx?
  B. density onset     - Q_n{2,8,20,80}_res4: at what density does Q start to
                         bend, i.e. where does crowding replace per-cell physics?

Nothing here is fitted. Differences between arms are grid_dx and founder count
only. ROS hazard is off (oxygen.k_ROS = 0, retired in #332) and founders carry no
plasmids, so colicin kills cannot confound the ledger: the only loss channels are
washout and boundary outflow.

Primary estimator, which needs no mu bookkeeping:
    Q_O2_per_gen = (agent_uptake_cumulative + maintenance_cumulative) / divisions
Cross-check via the growth rate:
    Q = funded_O2 / cell_seconds * ln(2) / mean_mu_realized
Both are reported; disagreement between them is itself a finding.
"""
import json
import os
import shutil
import subprocess
import time

BASE = "/home/ubuntu/gutibm-campaign/far-field-invariance/C_ff10_res2/input.json"
ROOT = "/home/ubuntu/gutibm-campaign/q-o2-per-gen"
BINARY = os.environ.get(
    "GUTIBM_BINARY", "/home/ubuntu/repos/GutModelBacteriocins/build/gut_ibm")

# arm -> (grid_dx, founder count)
ARMS = {
    "Q_n2_res2": (2.0e-6, 2),
    "Q_n2_res4": (4.0e-6, 2),
    "Q_n2_res6": (6.0e-6, 2),
    "Q_n8_res4": (4.0e-6, 8),
    "Q_n20_res4": (4.0e-6, 20),
    "Q_n80_res4": (4.0e-6, 80),
}

TOTAL_TIME = 43200.0  # 12 h: several generations even at T_gen ~ 4 h.

COMMENT = [
    "Q_O2_per_gen measurement on merged main after #338 (positivity rationing).",
    "Nothing here is fitted. Only grid_dx and founder count vary between arms.",
    "oxygen.k_ROS = 0 and founders carry no plasmids, so lysis cannot confound",
    "the oxygen ledger; the only loss channels are washout and boundary outflow.",
    "Read per arm: cumulative funded growth + maintenance oxygen, cumulative",
    "divisions, cell-seconds from the N trajectory, delivery_rationing_factor",
    "(must be 1.0 where the patch is genuinely uncrowded), delivery_infeasible",
    "(must be 0), reaction clips (must be 0), and /run_provenance termination.",
    "Q is only usable as a lysis-coefficient divisor if it agrees across",
    "2/4/6 um AND the rationing factor at that density is 1.0.",
]


def build_configs():
    with open(BASE) as fh:
        base = json.load(fh)
    for arm, (dx, count) in ARMS.items():
        cfg = dict(base)
        cfg["grid_dx"] = dx
        cfg["total_time"] = TOTAL_TIME
        cfg["output_interval"] = 300
        cfg["initial_strains"] = [{
            "type": 1,
            "count": count,
            "mu_max": 0.00055,
            "plasmids": [],
            "conjugative": False,
        }]
        cfg["hdf5"] = {
            "schedule": {
                "summary": 1,
                "agents": 5,
                "grid": 20,
                "grid_species": ["oxygen", "carbon"],
                "lineage": 0,
                "genome": 0,
            }
        }
        cfg["_comment"] = COMMENT + [
            f"arm = {arm}. grid_dx = {dx}. founders = {count}."]
        d = os.path.join(ROOT, arm)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "input.json"), "w") as fh:
            json.dump(cfg, fh, indent=1)


def run():
    for arm in ARMS:
        d = os.path.join(ROOT, arm)
        start = time.monotonic()
        with open(os.path.join(d, "full.log"), "w") as log:
            rc = subprocess.call([BINARY, "input.json"], cwd=d,
                                 stdout=log, stderr=subprocess.STDOUT)
        wall = time.monotonic() - start
        with open(os.path.join(d, "wrapper_wall_seconds.txt"), "w") as fh:
            fh.write(f"{wall}\n")
        print(f"{arm}: rc={rc} wall={wall:.1f}s", flush=True)


if __name__ == "__main__":
    assert os.path.exists(BINARY), BINARY
    assert os.path.exists(BASE), BASE
    os.makedirs(ROOT, exist_ok=True)
    shutil.copy(BINARY, os.path.join(ROOT, "gut_ibm.binary"))
    build_configs()
    run()
