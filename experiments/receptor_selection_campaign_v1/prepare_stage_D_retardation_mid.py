#!/usr/bin/env python3
"""Generate the Stage D retardation-sweep middle-point fills (amp 15 x 2 targets x 3 seeds).

Companion to prepare_stage_D_retardation.py. On the refine seeds
(20261009/1011/1013) the deployed amp-15 gradient has no measurements at
targets 1% and 5% — the original grid and seeds10 arrays ran those
targets on different seeds. These six runs supply the same-seed amp-15
anchors so every cell of the three-amplitude x three-target comparison
is paired within seed. Identical construction otherwise.
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
SEEDS = [20261009, 20261011, 20261013]
TARGETS = [0.01, 0.05]
AMPLITUDE = 15.0

prepare.SEEDS = SEEDS
prom = dict(prepare.CONTRACT["planning_promotions"])
kd = float(prom["selected_kd_corrinoid_btuB_mol_m3"])
b12 = float(prom["selected_b12_initial_conc_mol_m3"])
fixed = {**prepare.BASE_UPDATES, "kd_corrinoid_btuB": kd, "b12_initial_conc": b12}

runs = []
for target in TARGETS:
    p = 1.0 - math.sqrt(1.0 - target)
    for seed in SEEDS:
        rid = f"D_amp{AMPLITUDE:g}_target{target:.4f}_producer_s{seed}"
        updates = {**fixed, prepare.MUCIN_AMPLITUDE_KEY: AMPLITUDE,
                   "sos_basal_rate": 0.0, "sos_lysis_prob": p}
        axes = {"nominal_target_per_generation": target,
                "sos_lysis_prob": p, "selected_amplitude": AMPLITUDE}
        runs.append((rid, prepare.make(
            "D", rid, "producer", seed, axes,
            [prepare.strain(1, ["ColE1"]), prepare.strain(2)],
            EXEC_SHA, updates)))

outdir = PKG / "generated" / "stage_D_retardation_mid"
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
        "output_relpath": f"generated/stage_D_retardation_mid/results/{idx}/output.h5.gz",
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
                         "selected_b12_initial_conc_mol_m3": b12},
    "grid_refinement": {
        "extends_stage": "D",
        "extends_arrays": ["adfd5af1-3fd1-4ea8-b3d2-6b7d40609d66"],
        "new_targets": TARGETS,
        "amplitude": AMPLITUDE,
        "seeds": SEEDS,
        "note": "Same-seed amp-15 anchors at targets 1% and 5% for the "
                "retardation sweep — the deployed amp-15 grid ran those "
                "targets on other seeds. Pairs each new arm against the "
                "refine array's same-seed plasmid-free null.",
    },
    "runs": entries,
}
prepare.dump(outdir / "manifest.json", manifest)
print("wrote", len(entries), "jobs to", outdir)
