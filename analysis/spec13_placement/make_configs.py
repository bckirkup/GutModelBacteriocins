"""Emit the spec13 placement-decision arm configs.

One representative arm of the C/S/R selection campaign inside a single
`layer2.audit_patch_index` live patch: producer (ColE1), sensitive
(default BtuB), and resistant (receptor_expression {"BtuB": 0.0})
founders co-inoculated in the same patch, lysis lever swept, paired
seeds.  See docs/SPEC13_LYSIS_SELECTION.md for the campaign design and
analysis/spec13_live_patch/ for the live-patch contract.

Arms:
  crypt  x ret{lit,ship} x lysis{0, .02, .10}  (ret=ship only 0/.10)
  prox   x ret{lit}      x lysis{0, .02, .10}
x 5 seeds.  Mutation and conjugation are zeroed so the three-strain
composition is fixed (SPEC13_LYSIS_SELECTION SS4a caveat 3): lineage
identity stays strain identity for the whole run.
"""

import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

SEEDS = [20261027, 20261028, 20261029, 20261030, 20261031]
LYSIS_FULL = [0.0, 0.02, 0.10]
LYSIS_SHIP = [0.0, 0.10]

STRAINS = [
    # C: ColE1 producer — the only strain carrying an SOS-lysis locus,
    # so it is the only strain that pays the lysis cost and the only
    # toxin source.
    {"type": 1, "count": 10, "mu_max": 5e-4, "plasmids": ["ColE1"]},
    # S: sensitive, default receptors.
    {"type": 2, "count": 10, "mu_max": 5e-4},
    # R: btuB-null phenocopy — immune to ColE1/E2 entry (SPEC13
    # _LYSIS_SELECTION SS4a); cost comes from the Km_b12 clamp only.
    {"type": 3, "count": 10, "mu_max": 5e-4,
     "receptor_expression": {"BtuB": 0.0}},
]

# Zero the composition-changing knobs: mutation would mint de-novo
# resistance, conjugation would turn S into C.
ZEROED = {
    "receptor_mutation_rate": 0.0,
    "partial_resistance_rate": 0.0,
    "bi_duplication_rate": 0.0,
    "bi_recombination_rate": 0.0,
    "super_killer_rate": 0.0,
    "compensatory_rate": 0.0,
    "immunity_escape_prob": 0.0,
    "base_transfer_rate": 0.0,  # conjugation
    "sos_basal_rate": 0.0,      # per-division route carries the whole lever
}

PATCH_TYPES = {
    "crypt": {
        "layer2.crypt_fraction": 1.0,
        "layer2.proximal_fraction": 0.0,
        "layer2.distal_fraction": 0.0,
    },
    "prox": {
        "layer2.crypt_fraction": 0.0,
        "layer2.proximal_fraction": 1.0,
        "layer2.distal_fraction": 0.0,
    },
}


def arm_config(patch_type, ret, lysis, seed, run_dir):
    cfg = {
        "total_time": 86400,
        "bio_dt": 60,
        "seed": seed,
        "sos_lysis_prob": lysis,
        "initial_strains": STRAINS,
        "layer2.enabled": True,
        "layer2.n_patches": 10,
        "layer2.initial_patch_fraction": 1.0,
        "layer2.initial_founder_cells": 30,
        "layer2.patch_lateral_m": 60e-6,
        "layer2.patch_depth_m": 60e-6,
        "layer2.contraction_rate_per_min": 0.5,
        "layer2.agent_loss_fraction": 0.3,
        "layer2.transit_half_life_s": 14400,
        "layer2.reattach_prob_per_transit": 0.05,
        "layer2.establish_prob_single": 0.01,
        "layer2.establish_ratio": 10.0,
        "layer2.invalid_patch_fraction_stop": 1.0,
        "layer2.summary_interval_steps": 30,
        "layer2.audit_patch_index": 0,
        "layer2.audit_halt": False,
        "layer2.audit_lineage_file": os.path.join(run_dir, "lineage.csv"),
        "layer2.provenance_file": os.path.join(run_dir, "prov.json"),
        "layer2.timeseries_file": os.path.join(run_dir, "ts.csv"),
    }
    cfg.update(ZEROED)
    cfg.update(PATCH_TYPES[patch_type])
    if ret == "lit":
        # Literature transport analogue (SPEC13 SS4b): D_eff =
        # 7e-11/1.2 = 5.8e-11 m2/s vs shipped 4e-11/50.2 = 8e-13.
        cfg["plasmid_overrides"] = {
            "ColE1": {"diff_coeff": 7e-11, "retardation": 1.2}
        }
    return cfg


def main():
    arms = []
    for patch_type, ret, lyses in (
        ("crypt", "lit", LYSIS_FULL),
        ("crypt", "ship", LYSIS_SHIP),
        ("prox", "lit", LYSIS_FULL),
    ):
        for lysis in lyses:
            for seed in SEEDS:
                arms.append((patch_type, ret, lysis, seed))
    manifest = []
    for patch_type, ret, lysis, seed in arms:
        name = f"{patch_type}_{ret}_ly{lysis:.3g}_s{seed}"
        run_dir = os.path.join(OUT, name)
        os.makedirs(run_dir, exist_ok=True)
        cfg = arm_config(patch_type, ret, lysis, seed, run_dir)
        path = os.path.join(HERE, "configs", f"{name}.json")
        with open(path, "w") as fh:
            json.dump(cfg, fh, indent=2)
            fh.write("\n")
        manifest.append({"name": name, "config": path, "out": run_dir})
    with open(os.path.join(HERE, "configs", "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)
        fh.write("\n")
    print(f"wrote {len(manifest)} configs")


if __name__ == "__main__":
    main()
