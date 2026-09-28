#!/usr/bin/env python3
"""Generate the Stage D seed-extension inputs (5 arms x 10 new seeds).

Reuses stage_d_runs() from prepare.py so every arm parameter is identical to
the original Stage D generation; only SEEDS is overridden. Must be run from a
checkout at the execution source SHA (29fd563c) so the generated inputs carry
the correct provenance; this script itself is a record committed later.

Verification performed before submission (2026-09-28): all 50 generated
inputs are byte-identical to the deployed stage-d S3 inputs except `seed`,
`_campaign.run_id`, and `_campaign.paired_seed`.
"""
import hashlib
import json
import shutil
import sys
from pathlib import Path

PKG = Path(__file__).resolve().parent
sys.path.insert(0, str(PKG))
import prepare

EXEC_SHA = "29fd563ce73b28059cb439694e92466f451ddb4a"
IMAGE_DIGEST = "sha256:fafc77ee802882bee6e47321515057aed20c0326d622bfb287862ef15304f748"
NEW_SEEDS = [20260919, 20260921, 20260923, 20260925, 20260927,
             20260929, 20261001, 20261003, 20261005, 20261007]

prepare.SEEDS = NEW_SEEDS
prom = dict(prepare.CONTRACT["planning_promotions"])
kd = float(prom["selected_kd_corrinoid_btuB_mol_m3"])
b12 = float(prom["selected_b12_initial_conc_mol_m3"])
amp0 = float(prom["selected_mucin_charge_amplitude"])
fixed = {**prepare.BASE_UPDATES, "kd_corrinoid_btuB": kd, "b12_initial_conc": b12}

runs = prepare.stage_d_runs(kd, b12, amp0, fixed, EXEC_SHA)
outdir = PKG / "generated" / "stage_D2"
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
        "output_relpath": f"generated/stage_D2/results/{idx}/output.h5.gz",
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
    "seed_extension": {
        "extends_stage": "D",
        "extends_array": "47460bfa-1558-47b6-9a48-47b1bd9a7b0d",
        "original_seeds": sorted(prepare.CONTRACT["seeds"]),
        "new_seeds": NEW_SEEDS,
        "note": "Identical 5 arms rerun against 10 fresh seeds; all other parameters unchanged.",
    },
    "runs": entries,
}
prepare.dump(outdir / "manifest.json", manifest)
print("wrote", len(entries), "jobs to", outdir)
