#!/usr/bin/env python3
"""Generate the Stage D retardation-sensitivity inputs (2 amplitudes x 3 targets x 3 seeds).

Replicates stage_d_runs() from prepare.py so every parameter except the
D-axis (target, sos_lysis_prob) pair and the mucin-charge amplitude is
identical to the deployed Stage D inputs. Must be run from a checkout at
the execution source SHA (29fd563c) so the generated inputs carry the
correct provenance; this script itself is a record committed later.

Purpose: the deployed Stage D gradient was measured entirely at the
selected amplitude 15, where retardation_from_pI(9.0) resolves to ~13.5.
SPEC13 requires the response at both transport endpoints — amplitude 0
(ColE1 retardation 1.2, literature D_eff ~3.3e-11 m2/s) and amplitude 60
(retardation ~50.2, the pI-derived shipped value). If the dose-response
moves along the amplitude axis, the gradient measures the transport
parameter, not the lysis prior.

Arms: producer targets 0.25%, 1%, 5% per generation — knee, mid-curve,
and plateau of the amp-15 gradient — at amplitudes 0 and 60, on the same
three seeds as the stage-d-refine array (20261009/1011/1013). The
existing amp-15 refine runs on those seeds at the same targets complete
the three-amplitude comparison, and the existing same-seed plasmid-free
nulls and zero-lysis carriers remain valid shared controls (a null
releases no colicin, and a P=0 carrier emits no toxin, so amplitude does
not act on either). No new seeds, no new controls.

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
SEEDS = [20261009, 20261011, 20261013]
TARGETS = [0.0025, 0.01, 0.05]
AMPLITUDES = [0.0, 60.0]

prepare.SEEDS = SEEDS
prom = dict(prepare.CONTRACT["planning_promotions"])
kd = float(prom["selected_kd_corrinoid_btuB_mol_m3"])
b12 = float(prom["selected_b12_initial_conc_mol_m3"])
fixed = {**prepare.BASE_UPDATES, "kd_corrinoid_btuB": kd, "b12_initial_conc": b12}

runs = []
for amp in AMPLITUDES:
    for target in TARGETS:
        p = 1.0 - math.sqrt(1.0 - target)
        for seed in SEEDS:
            rid = f"D_amp{amp:g}_target{target:.4f}_producer_s{seed}"
            updates = {**fixed, prepare.MUCIN_AMPLITUDE_KEY: amp,
                       "sos_basal_rate": 0.0, "sos_lysis_prob": p}
            axes = {"nominal_target_per_generation": target,
                    "sos_lysis_prob": p, "selected_amplitude": amp}
            runs.append((rid, prepare.make(
                "D", rid, "producer", seed, axes,
                [prepare.strain(1, ["ColE1"]), prepare.strain(2)],
                EXEC_SHA, updates)))

outdir = PKG / "generated" / "stage_D_retardation"
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
        "output_relpath": f"generated/stage_D_retardation/results/{idx}/output.h5.gz",
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
        "new_amplitudes": AMPLITUDES,
        "new_targets": TARGETS,
        "seeds": SEEDS,
        "note": "Retardation-sensitivity sweep: same targets and seeds as "
                "the amp-15 refine arms at the two transport endpoints "
                "(retardation 1.2 and ~50.2 for ColE1). Existing "
                "same-seed nulls/carriers remain valid shared controls.",
    },
    "runs": entries,
}
prepare.dump(outdir / "manifest.json", manifest)
print("wrote", len(entries), "jobs to", outdir)
