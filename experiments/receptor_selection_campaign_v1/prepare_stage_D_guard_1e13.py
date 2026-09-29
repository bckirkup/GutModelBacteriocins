#!/usr/bin/env python3
"""Generate the Stage D dysbiosis-guard diagnostic inputs (4 arms x 2 seeds).

Reuses stage_d_runs() from prepare.py so every arm parameter is identical to
the deployed Stage D generation; the ONLY physics change is
``dysbiosis_threshold``: 1e10 -> 1e13 cells/mL. The diagnostic asks whether all
arms reach the 21600 s horizon (previous runs guard-halted at ~16800-20400 s)
and whether the arm ordering is retained when the guard no longer truncates.

This is a diagnostic, not new independent replicates: the threshold changes the
stopping regime, so outputs are not pooled with the 13-seed Stage D set.

Must be run from a checkout at the execution source SHA (29fd563c) so the
generated inputs carry the correct provenance; this script is a record
committed later.

Arm subset (per the diagnostic spec): plasmid-free null plus producer targets
0.000, 0.010 and 0.050 -- the 0.020 arm is dropped; seeds are the two original
contract seeds 20260911 and 20260913 so each run pairs against an already-run
same-seed Stage D output for ordering comparison.
"""
import hashlib
import shutil
import sys
from pathlib import Path

PKG = Path(__file__).resolve().parent
sys.path.insert(0, str(PKG))
import prepare

EXEC_SHA = "29fd563ce73b28059cb439694e92466f451ddb4a"
IMAGE_DIGEST = "sha256:fafc77ee802882bee6e47321515057aed20c0326d622bfb287862ef15304f748"
SEEDS = [20260911, 20260913]
DYSBIOSIS_THRESHOLD = 1e13
KEEP_TARGETS = {0.0, 0.01, 0.05}

prepare.SEEDS = SEEDS
prom = dict(prepare.CONTRACT["planning_promotions"])
kd = float(prom["selected_kd_corrinoid_btuB_mol_m3"])
b12 = float(prom["selected_b12_initial_conc_mol_m3"])
amp0 = float(prom["selected_mucin_charge_amplitude"])
fixed = {**prepare.BASE_UPDATES, "kd_corrinoid_btuB": kd, "b12_initial_conc": b12,
         "dysbiosis_threshold": DYSBIOSIS_THRESHOLD}

runs = prepare.stage_d_runs(kd, b12, amp0, fixed, EXEC_SHA)
runs = [(rid, cfg) for rid, cfg in runs
        if cfg["_campaign"]["arm"] == "plasmid_free_null"
        or cfg["_campaign"]["axes"]["nominal_target_per_generation"] in KEEP_TARGETS]

outdir = PKG / "generated" / "stage_D_guard_1e13"
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
        "output_relpath": f"generated/stage_D_guard_1e13/results/{idx}/output.h5.gz",
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
    "guard_diagnostic": {
        "extends_stage": "D",
        "changed_parameters": {"dysbiosis_threshold": DYSBIOSIS_THRESHOLD},
        "previous_dysbiosis_threshold": 1e10,
        "seeds": SEEDS,
        "arms": "plasmid_free_null + producer targets {0.000, 0.010, 0.050} (0.020 dropped)",
        "note": ("Diagnostic, not pooled replicates: raising the guard threshold "
                 "changes the stopping regime. Pass condition: every arm reaches "
                 "21600 s (termination_reason_code == 0) and arm ordering is "
                 "retained vs the 13-seed Stage D readout."),
    },
    "runs": entries,
}
prepare.dump(outdir / "manifest.json", manifest)
print("wrote", len(entries), "jobs to", outdir)
