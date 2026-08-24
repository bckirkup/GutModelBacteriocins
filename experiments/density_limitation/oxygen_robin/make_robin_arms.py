#!/usr/bin/env python3
"""Can a depletable epithelium fund respiration, and is supply or density the limit?

Why this experiment exists
--------------------------
The measurement in #325 found growth respiration funded at exactly zero in all
eight arms of the funded-O2 x placement x consumption matrix, and located two
reasons that placement was not one of them. First, `calculate_delivery_funding`
pays maintenance before growth, and oxygen maintenance was never fully paid, so
growth respiration was structurally zero at every depth. Second, the epithelial
face was Dirichlet plus a per-step maintained 25 um exponential, which cannot be
depleted by consumption: the O2 z-profile was unchanged between arms whose
background sink differed by 840x.

#326 made the oxygen epithelial boundary configurable, so the mucus-surface
concentration can now emerge from `J = k*(C_tissue - C_surface)` instead of
being imposed. This experiment asks the two questions #325 left open, and it is
built so that either answer is informative:

  Q1. With a depletable, finite-rate epithelium, is growth respiration funded at
      all -- and does it survive past the initial transient?
  Q2. Is the binding constraint supply (the apical flux) or density (the ~300x
      excess founder density #325 identified)?

Arms (4, deliberately not a full cross)
---------------------------------------
All arms use the `funded` driver, delivery-limited carbon and oxygen uptake, and
the anatomy-derived placement of #324. Placement is held fixed because #325
showed it does not affect funding.

`dirichlet`  -- control. Imposed 25 um profile at 5.5e-2 mol/m^3, i.e. the
                configuration of #325's `funded_anatomic_shipped`. It must
                reproduce zero funded growth O2; if it does not, something in
                #326 changed the Dirichlet path and nothing else here is
                interpretable.
`robin`       -- k = 1.2e-6 m/s, C_tissue = 5.0e-2 mol/m^3. Both sourced:
                the transfer coefficient is measured rat mucus+mucosa
                permeability, the reservoir is tissue-side O2. Gradient off, so
                the profile must be built by the boundary.
`robin_k10`   -- k = 1.2e-5 m/s. A deliberate 10x overshoot of the sourced
                permeability. This is a sensitivity probe, not a candidate
                value: it bounds how much of any negative result is attributable
                to the flux magnitude rather than to the funding rule.
`robin_sparse`-- k = 1.2e-6 m/s with 8 founders instead of 80. Same supply, 10x
                lower density, 2.9e6 cells/mL instead of 2.9e7. This separates
                Q2 from Q1: if `robin` cannot fund respiration but
                `robin_sparse` can, the limit is density and the single-patch
                geometry is the problem, which is the Spec 13 metapopulation
                argument. If neither can, it is supply.

Held fixed
----------
Seed, horizon, geometry, carbon supply (1.00x J_dir), flora (1.00x) and guard
match ../supply_ladder and ../oxygen_placement so population numbers stay
comparable. `oxygen.vbf_sink` stays at the shipped 1e-3 1/s: #325 showed the
0.84 1/s "consumption-consistent" value has 98% of its demand clipped because
the oxygen is not there to consume, so it measures nothing until supply changes.
`anaerobic_maintenance_factor` is left at its shipped 15 -- its value is an open
question and changing it in the same experiment would confound the boundary.

Read per arm
------------
Cumulative funded vs demanded growth O2 (the headline: is it still exactly 0);
funded maintenance O2 against maintenance shortfall; the emergent O2 surface
concentration and the depth at which it falls to 1/e of that surface value,
compared with the 25 um the model used to impose; whether any funded respiration
is confined to the first hours (an initialization transient) or persists;
fermentation fraction; N trajectory and termination cause; carbon clip (must
stay 0).

Predictions, fixed before running
---------------------------------
1. `dirichlet` reproduces zero funded growth O2 exactly.
2. `robin` delivers ~5e-8 mol/m^2/s at a ~5 uM emergent surface, roughly 9x the
   6e-9 mol/m^2/s the Dirichlet face actually delivered in #325, which is enough
   to fund basal maintenance for a few hundred cells rather than two. So funded
   growth O2 should become nonzero for the first time.
3. The emergent 1/e depth will not be 25 um. There is no reason it should be:
   25 um was an imposed constant that only ever self-consistent with the
   1000x-low boundary value.
4. `robin_sparse` funds a larger per-cell share than `robin`. If the ratio is
   close to 10x, supply is simply divided among cells and density is the whole
   story.

Nothing here is fitted. No arm is tuned to make respiration appear, and the
`dirichlet` control exists so that a positive result can be attributed to the
boundary change rather than to anything else that has landed since #325.

Run with the project interpreter so ``gut_ibm_tools`` is importable, e.g.
``PYTHONPATH=python /home/ubuntu/.venvs/gutibm/bin/python3 <this script> ...``.
"""

import argparse
import json
import os
import pathlib

from gut_ibm_tools.path_utils import (
    validate_input_path,
    validate_output_path,
    write_json_file,
)

SEED = 1001
HORIZON_S = 86400
J_DIR = 1.0756e-08
VBF_VMAX = 5.5e-05
GUARD_CELLS_PER_ML = 1.0e10

# Imposed apical concentration for the Dirichlet control (~42 mmHg dissolved),
# matching #325 so the control is comparable arm for arm.
DIRICHLET_O2 = 5.5e-2

# Tissue-side reservoir for the Robin arms, mol/m^3. Sourced tissue-side O2.
C_TISSUE = 5.0e-2

# Robin mass-transfer coefficient, m/s. 1.2e-6 is measured rat mucus+mucosa
# permeability; 1.2e-5 is a 10x sensitivity bound, not a candidate value.
K_SOURCED = 1.2e-6
K_OVERSHOOT = 1.2e-5

# Anatomy-derived placement (#324). Held fixed: #325 showed placement does not
# affect funding.
ANATOMIC = {
    "initial_population.placement": "anatomic",
    "initial_population.anatomic_exclusion_floor": 20.0e-6,
    "initial_population.anatomic_exponential_scale": 40.0e-6,
    "initial_population.anatomic_outer_extent": 150.0e-6,
}

# arm -> (boundary mode, k, founder scale)
ARMS = {
    "dirichlet": ("dirichlet", 0.0, 1.0),
    "robin": ("robin", K_SOURCED, 1.0),
    "robin_k10": ("robin", K_OVERSHOOT, 1.0),
    "robin_sparse": ("robin", K_SOURCED, 0.1),
}


def scaled_strains(base: dict, scale: float) -> list:
    strains = json.loads(json.dumps(base["initial_strains"]))
    for strain in strains:
        strain["count"] = max(1, round(strain["count"] * scale))
    return strains


def build_arm(base: dict, tag: str) -> dict:
    mode, k_transfer, founder_scale = ARMS[tag]
    cfg = dict(base)
    cfg["total_time"] = HORIZON_S
    cfg["seed"] = SEED
    cfg["gpu_enabled"] = False
    cfg["uptake_limit"] = "delivery"
    cfg["carbon_z_gradient"] = False
    cfg["carbon.epithelial_flux"] = J_DIR
    cfg["vbf_carbon_sink_vmax"] = VBF_VMAX
    cfg["dysbiosis_threshold"] = GUARD_CELLS_PER_ML

    cfg["oxygen.enabled"] = True
    cfg["oxygen.delivery_uptake_enabled"] = True
    cfg["oxygen.metabolic_switch_enabled"] = True
    cfg["oxygen.respiration_driver"] = "funded"
    cfg["oxygen.vbf_sink"] = 1.0e-3
    cfg["oxygen.epithelial_boundary"] = mode
    if mode == "robin":
        # Under Robin, epithelial_conc is the tissue-side reservoir and the
        # imposed profile must be off: the surface value is predicted, not set.
        cfg["oxygen.epithelial_conc"] = C_TISSUE
        cfg["oxygen.epithelial_transfer_coeff"] = k_transfer
        cfg["oxygen.z_gradient"] = False
    else:
        cfg["oxygen.epithelial_conc"] = DIRICHLET_O2
        cfg["oxygen.z_gradient"] = True

    cfg["acetate.enabled"] = True
    cfg["metabolism.acid_inhibition_enabled"] = True

    for key in ("initial_population.placement",
                "initial_population.z_min",
                "initial_population.z_max"):
        cfg.pop(key, None)
    cfg.pop("initial_population", None)
    cfg.update(ANATOMIC)
    cfg["initial_strains"] = scaled_strains(base, founder_scale)
    founders = sum(s["count"] for s in cfg["initial_strains"])

    cfg["hdf5_file"] = "output.h5"
    cfg["hdf5"] = {
        "schedule": {
            "summary": 1,
            "agents": 5,
            # 60 summary steps = 1 h simulated. Grid output is what makes the
            # emergent O2 profile measurable against the retired 25 um.
            "grid": 60,
            "grid_species": ["oxygen", "carbon", "acetate"],
            "lineage": 0,
            "genome": 0,
        }
    }

    cfg["_comment"] = [
        ("Robin epithelial O2 measurement run. Nothing here is fitted or "
         "calibrated."),
        (f"arm = {tag}. oxygen.epithelial_boundary = {mode}."),
        ("dirichlet is the control: it reproduces the configuration of #325's "
         "funded_anatomic_shipped and must reproduce zero funded growth O2. "
         "If it does not, #326 moved the Dirichlet path and no other arm here "
         "is interpretable."),
        (f"Robin k = {k_transfer:.3g} m/s, C_tissue = {C_TISSUE:.3g} mol/m^3. "
         "1.2e-6 is measured rat mucus+mucosa permeability; 1.2e-5 is a 10x "
         "sensitivity bound on the flux magnitude, not a candidate value. "
         "Under Robin the imposed 25 um exponential is off, so the O2 profile "
         "is built by the boundary rather than restored every step."),
        (f"{founders} founders ({founder_scale:.2g}x). #325 measured that 80 "
         "founders in this domain is 2.9e7 cells/mL, ~300x the culture-based "
         "healthy target of 1e4-1e5 CFU/mL. The sparse arm separates a supply "
         "limit from a density limit; it does not claim to reach the healthy "
         "target, which a single patch of this size cannot represent."),
        ("oxygen.vbf_sink stays at the shipped 1e-3 1/s: #325 showed the "
         "0.84 1/s consumption-consistent value has 98% of its demand clipped "
         "because the oxygen is not there to consume."),
        ("anaerobic_maintenance_factor is left at its shipped 15. Its value is "
         "an open question; changing it here would confound the boundary."),
        ("Read per arm: cumulative funded vs demanded growth O2, funded "
         "maintenance O2 against maintenance shortfall, emergent O2 surface "
         "concentration and 1/e depth against the retired 25 um, whether any "
         "funded respiration persists past the first hours or is a transient, "
         "fermentation fraction, N trajectory, carbon clip (must stay 0), and "
         "/run_provenance termination cause."),
    ]
    return cfg


def main(base_config: pathlib.Path, out_root: pathlib.Path) -> None:
    base = json.loads(validate_input_path(base_config).read_text())
    out_root = out_root.expanduser().resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    out_root = validate_output_path(out_root / "input.json").parent
    previous_cwd = pathlib.Path.cwd()
    os.chdir(out_root)
    try:
        for tag in ARMS:
            cfg = build_arm(base, tag)
            path = pathlib.Path(tag) / "input.json"
            write_json_file(path, cfg, indent=1)
            founders = sum(s["count"] for s in cfg["initial_strains"])
            print(f"{tag}: boundary={cfg['oxygen.epithelial_boundary']} "
                  f"k={ARMS[tag][1]:.3g} founders={founders} "
                  f"-> {out_root / path}")
    finally:
        os.chdir(previous_cwd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-config", type=pathlib.Path, required=True)
    parser.add_argument("--out-root", type=pathlib.Path, required=True)
    args = parser.parse_args()
    main(args.base_config, args.out_root)
