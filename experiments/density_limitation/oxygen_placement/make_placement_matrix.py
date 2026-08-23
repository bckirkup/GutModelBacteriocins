#!/usr/bin/env python3
"""Generate the funded-O2 x placement x O2-consumption matrix.

Why this experiment exists
--------------------------
Three separate things had to be fixed before any aerobic number from this model
could be believed, and all three are now merged: delivery-limited oxygen uptake
(#319), the unmetered epithelial Dirichlet face in the delivery solvers (#320),
and the metabolic switch reading ambient O2 the field never paid for (#323).
Every previously measured oxygen arm predates #320, so its boundary-supply
numbers are the broken ones and none of them can be reused here.

This matrix asks one question: is there a depth at which E. coli can fund
respiration in this model, and if not, which of placement or O2 consumption is
responsible?

Factors, fully crossed (2 x 2 x 2 = 8 arms)
-------------------------------------------
``driver`` -- ``ambient`` sets the respiratory fraction from local O2
concentration (the historical behaviour); ``funded`` sets it from
``funded_growth_O2 / demanded_growth_O2`` after diffusion, VBF and agent
competition. ``ambient`` is the control that shows what the model used to
claim; it is not a scientific alternative, because it reports respiration the
field did not supply.

``placement`` -- ``zslab`` is the 0-50 um uniform band that every arm we have
measured actually used, including the supply ladder of #314. ``anatomic`` is the
anatomy-derived policy: a hard 20 um exclusion floor (Enterobacteriaceae
bacteria-free gap, Swidsinski FISH), exponential with 40 um scale, truncated at
150 um (Duncan / Mondragon-Palomino outer-mucus enrichment), and no crypt
founders. Note that ``anatomic`` is mostly *deeper* than ``zslab``: the point is
not to move founders into the oxic shell but to stop placing 40% of them inside
a zone where they are never observed.

``o2_sink`` -- the first-order background O2 uptake rate constant. ``shipped``
is 1e-3 1/s. ``consistent`` is 0.84 1/s, and that number is derived, not
chosen: oxygen ships with a hard-coded 25 um ``z_gradient_lambda``
(``k_z_lambda`` in ``finalize_config``), while sqrt(D/k) at the shipped sink is
1449 um -- 58x deeper than the profile the model imposes. The imposed 25 um is
only self-consistent with the 1000x-low ``epithelial_conc = 55e-6``: at that
value and 5e8 cells/mL the consumption-derived depth is 21.5 um. Correcting the
boundary value alone therefore yields a ~680 um oxic layer, i.e. oxic mucus
everywhere, which is not what colonic mucus does. ``consistent`` is the k for
which sqrt(D/k) = 50 um at the physiological boundary value, so this axis tests
whether the shallow oxic shell survives being derived from consumption instead
of imposed by a constant.

Held fixed
----------
``oxygen.epithelial_conc = 5.5e-2 mol/m^3`` (~42 mmHg dissolved) in every arm.
The shipped 55e-6 is the same 10^3 unit error as the 25 um lambda and is not
worth another eight arms. Carbon supply is 1.00x J_dir with flora at 1.00x, and
seed, horizon, geometry and guard match ../supply_ladder so population numbers
remain comparable to the audited ladder.

Read per arm
------------
funded vs demanded growth O2; realized fermentation fraction against both mu
and local O2; acetate production; net growth and washout/outflow; final N; and
the z-distribution of survivors against the z-profile of O2. Grid output is
enabled here (it was off in the earlier arms) precisely so the O2 profile can be
compared with the imposed 25 um exponential.

Predictions, fixed before running
---------------------------------
1. ``ambient`` x anything reports a low fermentation fraction that the funded
   ledger does not support: the discrepancy between the two drivers at the same
   placement is the size of the accounting error #323 removed.
2. ``funded`` collapses respiration at both placements under the ``shipped``
   sink, because the depth question is not what limits it -- supply is.
3. The population nevertheless persists as fermenters. E. coli is facultative;
   unfunded respiration is a mode switch, not death. Population collapse in a
   ``funded`` arm would be a finding about ``anaerobic_mu_factor`` or carbon
   supply, not a reason to revert the switch.
4. If ``consistent`` recovers a thin funded shell, the shallow oxic layer is a
   consequence of consumption and the 25 um constant should be retired in favour
   of it. If it does not, apical O2 supply cannot fund respiration at any depth
   in this geometry, and the aerobic axis is inert for a physical reason rather
   than a plumbing one.

Nothing here is fitted. No arm is tuned to make respiration appear.

Run with the project interpreter so ``gut_ibm_tools`` is importable, e.g.
``PYTHONPATH=python /home/ubuntu/.venvs/gutibm/bin/python3 <this script> ...``.
"""

import argparse
import itertools
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

# ~42 mmHg dissolved. Held fixed across the matrix; see module docstring.
APICAL_O2 = 5.5e-2

# First-order background O2 uptake rate constant, 1/s.
#   shipped     -- code default; sqrt(D/k) = 1449 um
#   consistent  -- k for which sqrt(D/k) = 50 um at D = 2.1e-9 m^2/s
O2_SINK = {
    "shipped": 1.0e-3,
    "consistent": 8.4e-1,
}

# Anatomy-derived placement (PR #324). Values are measured, not fitted.
ANATOMIC = {
    "initial_population.placement": "anatomic",
    "initial_population.anatomic_exclusion_floor": 20.0e-6,
    "initial_population.anatomic_exponential_scale": 40.0e-6,
    "initial_population.anatomic_outer_extent": 150.0e-6,
}

# The band every previously measured arm actually used.
ZSLAB = {
    "initial_population.placement": "z_slab",
    "initial_population.z_min": 0.0,
    "initial_population.z_max": 5.0e-5,
}

PLACEMENT = {"zslab": ZSLAB, "anatomic": ANATOMIC}
DRIVERS = ("ambient", "funded")


def arm_tag(driver: str, placement: str, sink: str) -> str:
    return f"{driver}_{placement}_{sink}"


def build_arm(base: dict, driver: str, placement: str, sink: str) -> dict:
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
    cfg["oxygen.respiration_driver"] = driver
    cfg["oxygen.epithelial_conc"] = APICAL_O2
    cfg["oxygen.vbf_sink"] = O2_SINK[sink]

    cfg["acetate.enabled"] = True
    cfg["metabolism.acid_inhibition_enabled"] = True

    for key in ("initial_population.placement",
                "initial_population.z_min",
                "initial_population.z_max"):
        cfg.pop(key, None)
    cfg.pop("initial_population", None)
    cfg.update(PLACEMENT[placement])

    cfg["hdf5_file"] = "output.h5"
    cfg["hdf5"] = {
        "schedule": {
            "summary": 1,
            "agents": 5,
            # Grid output was off in the pre-#320 arms, which is why the O2
            # z-profile could never be compared with the imposed 25 um
            # exponential. 60 summary steps = 1 h of simulated time.
            "grid": 60,
            "grid_species": ["oxygen", "carbon", "acetate"],
            "lineage": 0,
            "genome": 0,
        }
    }

    cfg["_comment"] = [
        ("Funded-O2 x placement x O2-consumption matrix. Measurement run: "
         "nothing here is fitted or calibrated."),
        (f"respiration_driver = {driver}. ambient sets the respiratory "
         "fraction from local O2 concentration; funded sets it from "
         "funded/demanded growth O2 after diffusion, VBF and agent "
         "competition (#323). ambient is a control for the size of the "
         "accounting error, not a scientific alternative."),
        (f"placement = {placement}. zslab is the 0-50 um band every earlier "
         "arm used; anatomic is the measured policy of #324 (20 um "
         "Enterobacteriaceae exclusion floor, 40 um exponential scale, 150 um "
         "outer extent, no crypt founders)."),
        (f"oxygen.vbf_sink = {O2_SINK[sink]:.3g} 1/s ({sink}). Shipped 1e-3 "
         "gives sqrt(D/k) = 1449 um against the hard-coded 25 um oxygen "
         "z_gradient_lambda; 0.84 is the k for which sqrt(D/k) = 50 um. This "
         "axis tests whether a shallow oxic shell survives being derived from "
         "consumption rather than imposed by a constant."),
        (f"oxygen.epithelial_conc = {APICAL_O2:.3e} mol/m^3 (~42 mmHg "
         "dissolved) in every arm; the shipped 55e-6 is the same 10^3 unit "
         "error as the 25 um lambda. The code default is unchanged."),
        ("All earlier oxygen arms predate #320 (unmetered epithelial Dirichlet "
         "face in the delivery solvers) and cannot be compared with these."),
        ("Read per arm: funded vs demanded growth O2, realized fermentation "
         "fraction against mu and local O2, acetate, net growth and outflow, "
         "final N, survivor z-distribution against the O2 z-profile, carbon "
         "clip (must stay 0), and /run_provenance termination cause."),
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
        for driver, placement, sink in itertools.product(
                DRIVERS, PLACEMENT, O2_SINK):
            tag = arm_tag(driver, placement, sink)
            cfg = build_arm(base, driver, placement, sink)
            path = pathlib.Path(tag) / "input.json"
            write_json_file(path, cfg, indent=1)
            print(f"{tag}: driver={driver} placement={placement} "
                  f"vbf_sink={O2_SINK[sink]:.3g} -> {out_root / path}")
    finally:
        os.chdir(previous_cwd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-config", type=pathlib.Path, required=True)
    parser.add_argument("--out-root", type=pathlib.Path, required=True)
    args = parser.parse_args()
    main(args.base_config, args.out_root)
