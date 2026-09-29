#!/usr/bin/env python3
"""Generate the Stage D grid-refinement inputs (7 arms x 3 new seeds).

Replicates stage_d_runs() from prepare.py on a refined lysis-target grid so
every non-lysis parameter is identical to the deployed Stage D inputs; only
SEEDS and the D-axis (target, sos_lysis_prob) pairs are overridden. Must be
run from a checkout at the execution source SHA (29fd563c) so the generated
inputs carry the correct provenance; this script itself is a record
committed later.

Arms: plasmid-free null and ColE1 zero-lysis carrier (both controls,
re-run on new seeds), plus producer targets 0.25%, 0.5%, 0.75%, 1.5%, 3%
per generation — probing onset and curvature between the deployed
0/1/2/5% grid.

sos_lysis_prob mapping: each division draws the lysis hazard on both
mother and daughter, so realized per-generation lysis is 1-(1-p)^2;
p = 1 - sqrt(1 - t) with sos_basal_rate = 0 reproduces nominal target t.
"""
import hashlib
import math
import shutil
import sys
from pathlib import Path

PKG = Path(__file__).resolve().parent
sys.path.insert(0, str(PKG))
import prepare

EXEC_SHA = "29fd563ce73b28059cb439694e92466f451ddb4a"
IMAGE_DIGEST = "sha256:fafc77ee802882bee6e47321515057aed20c0326d622bfb287862ef15304f748"
NEW_SEEDS = [20261009, 20261011, 20261013]
TARGETS = [0, 0.0025, 0.005, 0.0075, 0.015, 0.03]

prepare.SEEDS = NEW_SEEDS
prom = dict(prepare.CONTRACT["planning_promotions"])
kd = float(prom["selected_kd_corrinoid_btuB_mol_m3"])
b12 = float(prom["selected_b12_initial_conc_mol_m3"])
amp0 = float(prom["selected_mucin_charge_amplitude"])
fixed = {**prepare.BASE_UPDATES, "kd_corrinoid_btuB": kd, "b12_initial_conc": b12}

runs = []
for target in TARGETS:
    p = 0 if target == 0 else 1.0 - math.sqrt(1.0 - target)
    for seed in NEW_SEEDS:
        rid = f"D_target{target:.4f}_producer_s{seed}"
        updates = {**fixed, prepare.MUCIN_AMPLITUDE_KEY: amp0,
                   "sos_basal_rate": 0.0, "sos_lysis_prob": p}
        axes = {"nominal_target_per_generation": target,
                "sos_lysis_prob": p, "selected_amplitude": amp0}
        runs.append((rid, prepare.make(
            "D", rid, "producer", seed, axes,
            [prepare.strain(1, ["ColE1"]), prepare.strain(2)],
            EXEC_SHA, updates)))
for seed in NEW_SEEDS:
    rid = f"D_plasmid_free_null_s{seed}"
    updates = {**fixed, prepare.MUCIN_AMPLITUDE_KEY: amp0,
               "sos_basal_rate": 0.0, "sos_lysis_prob": 0}
    axes = {"control": "plasmid_free_null", "selected_amplitude": amp0,
            "selected_kd": kd, "selected_b12": b12}
    runs.append((rid, prepare.make(
        "D", rid, "plasmid_free_null", seed, axes,
        [prepare.strain(1), prepare.strain(2)], EXEC_SHA, updates)))

outdir = PKG / "generated" / "stage_D_refine"
if outdir.exists():
    shutil.rmtree(outdir)
entries = []
for idx, (rid, cfg) in enumerate(runs):
    ip = outdir / "jobs" / str(idx) / "input.json"
    prepare.dump(ip, cfg)
    entries.append({
        "array_index": idx,
        "run_id": rid,
        "stage": "D",
        "input_relpath": str(ip.relative_to(PKG)),
        "input_sha256": hashlib.sha256(ip.read_bytes()).hexdigest(),
        "output_relpath": f"generated/stage_D_refine/results/{idx}/output.h5.gz",
    })
manifest = {
    "schema_version": 2,
    "campaign_id": prepare.CONTRACT["campaign_id"],
    "stage": "D",
    "model_baseline_sha": prepare.BASELINE,
    "execution_source_sha": EXEC_SHA,
    "container_image_digest": IMAGE_DIGEST,
    "mpi_ranks": 1,
    "gpu_required": True,
    "attempt_timeout_s": 7200,
    "gate_requires": "C_transport_gate",
    "gate_locked": True,
    "promotion_inputs": {**prom, "selected_kd_corrinoid_btuB_mol_m3": kd,
                         "selected_b12_initial_conc_mol_m3": b12,
                         "selected_mucin_charge_amplitude": amp0},
    "grid_refinement": {
        "extends_stage": "D",
        "extends_arrays": ["47460bfa-1558-47b6-9a48-47b1bd9a7b0d",
                            "734ab0fb-7df4-4023-80d7-e703075d8b23"],
        "deployed_targets": prepare.CONTRACT["axes"]["D"]["realized_lysis_target_per_generation"],
        "new_targets": TARGETS,
        "new_seeds": NEW_SEEDS,
        "note": "Refined lysis-target grid probing onset and curvature between the "
                "deployed 0/1/2/5% settings; both controls re-run on new seeds.",
    },
    "runs": entries,
}
prepare.dump(outdir / "manifest.json", manifest)
print("wrote", len(entries), "jobs to", outdir)
