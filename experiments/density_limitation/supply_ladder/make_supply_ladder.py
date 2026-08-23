#!/usr/bin/env python3
"""Generate the carbon supply x flora-competition ladder for delivery mode.

Why this experiment exists
--------------------------
Every density mechanism probed so far (maintenance #302/#306, delivery #308,
voxel-local ``agent_carbon_coupling`` #312/#313) was measured at a single
supply, ``1.00x J_dir``, where the founding population is maintenance-starved:
funded fraction 0.128, N falls 80 -> 47 in 6 h with every optional mechanism
switched off. A brake added to a population that is already dying cannot be
observed, so those probes returned nulls by construction.

This ladder measures the observable instead of adding a mechanism: carrying
capacity as a function of the carbon actually available to the modelled
population. Two axes, because there are exactly two ways to give the agents
more carbon:

* ``supply`` -- multiplier on ``carbon.epithelial_flux``. 1.00x is the
  literature directly-measured host flux; multipliers above 1 are
  supra-physiological and are included to locate the capacity curve, not as
  claims about the colon.
* ``flora`` -- multiplier on ``vbf_carbon_sink_vmax``. At 1.00x J_dir the
  background flora takes ~98.6% of all realized carbon removal, so lowering it
  raises the agents' share without inventing host flux. 0.0 is the
  no-competitor control, which bounds the achievable population.

Pass criterion for the campaign (decided before running):
capacity must be monotone in available carbon, and the nutrient blocking
fraction must rise with N. If N_plateau is flat across the grid, the delivery
brake is not setting capacity and no coupling constant will change that.

Measurement conventions
-----------------------
``uptake_limit=delivery`` so biomass is funded only by realized field removal
(#308), ``carbon_z_gradient=false`` so the delivery sink is not measured
through the profile decomposition (#310), ``gpu_enabled=false`` because
delivery rejects the GPU path, one fixed seed so arms differ only by the axes.

``dysbiosis_threshold`` is deliberately raised to 1e10 cells/mL. The shipped
guard at 1e8 (~276 agents in this domain) would halt a high-supply arm before
it plateaus, which would report the guard rather than the capacity. Any arm
whose plateau sits above 1e8 cells/mL is therefore a *guard-exceeding* result
and must not be read as a physiological population.

Each arm holds ~2.7 GB RSS; run them strictly one at a time.
"""

import argparse
import importlib.util
import json
import os
import pathlib
import sys
import types

_REPO_PYTHON = pathlib.Path(__file__).resolve().parents[3] / "python"
sys.path.insert(0, str(_REPO_PYTHON))

try:
    from gut_ibm_tools.path_utils import (
        validate_input_path,
        validate_output_path,
        write_json_file,
    )
except ModuleNotFoundError as error:
    if error.name != "h5py":
        raise
    package = types.ModuleType("gut_ibm_tools")
    package.__path__ = [str(_REPO_PYTHON / "gut_ibm_tools")]
    sys.modules["gut_ibm_tools"] = package
    module_path = _REPO_PYTHON / "gut_ibm_tools" / "path_utils.py"
    spec = importlib.util.spec_from_file_location(
        "gut_ibm_tools.path_utils", module_path
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {module_path}") from error
    path_utils = importlib.util.module_from_spec(spec)
    sys.modules["gut_ibm_tools.path_utils"] = path_utils
    spec.loader.exec_module(path_utils)
    validate_input_path = path_utils.validate_input_path
    validate_output_path = path_utils.validate_output_path
    write_json_file = path_utils.write_json_file

SEED = 1001
HORIZON_S = 86400  # 24 h; capacity needs a plateau, not a transient
J_DIR = 1.0756e-08  # 1.00x directly-measured epithelial carbon flux
VBF_VMAX = 5.5e-05  # shipped background flora carbon sink
GUARD_CELLS_PER_ML = 1.0e10  # raised for measurement; see module docstring

SUPPLY = [1.0, 2.0, 4.0, 8.0, 16.0]
FLORA = [1.0, 0.5, 0.0]


def arm_tag(supply: float, flora: float) -> str:
    return f"s{supply:g}_v{flora:g}".replace(".", "p")


def main(base_config: pathlib.Path, out_root: pathlib.Path) -> None:
    base = json.loads(validate_input_path(base_config).read_text())
    out_root = out_root.expanduser().resolve()
    out_root = validate_output_path(out_root / "input.json").parent
    previous_cwd = pathlib.Path.cwd()
    os.chdir(out_root)

    try:
        for supply in SUPPLY:
            for flora in FLORA:
                cfg = dict(base)
                cfg["total_time"] = HORIZON_S
                cfg["seed"] = SEED
                cfg["gpu_enabled"] = False
                cfg["uptake_limit"] = "delivery"
                cfg["carbon_z_gradient"] = False
                cfg["carbon.epithelial_flux"] = J_DIR * supply
                cfg["vbf_carbon_sink_vmax"] = VBF_VMAX * flora
                cfg["dysbiosis_threshold"] = GUARD_CELLS_PER_ML
                cfg["hdf5_file"] = "output.h5"
                cfg["_comment"] = [
                    "Carbon supply x flora-competition ladder, delivery mode.",
                    "Measurement run: nothing here is fitted or calibrated.",
                    (f"supply = {supply:g}x J_dir "
                     f"(carbon.epithelial_flux = {J_DIR * supply:.4e}); "
                     "1.00x is the literature directly-measured flux, higher "
                     "multipliers are supra-physiological probes of the capacity "
                     "curve, not claims about the colon."),
                    (f"flora = {flora:g}x vbf_carbon_sink_vmax "
                     f"({VBF_VMAX * flora:.3e}); at 1.00x the flora takes ~98.6% "
                     "of realized carbon removal, so this axis raises the agents' "
                     "share without inventing host flux. 0.0 bounds capacity."),
                    (f"dysbiosis_threshold raised to {GUARD_CELLS_PER_ML:.0e} "
                     "cells/mL so a high arm plateaus instead of halting on the "
                     "guard; any plateau above 1e8 cells/mL is guard-exceeding "
                     "and is not a physiological population."),
                    ("Read per arm: N_plateau over the final quarter and its "
                     "slope, nutrient blocking fraction, funded fraction, "
                     "maintenance shortfall, reaction clip, and "
                     "/run_provenance termination cause."),
                ]

                output_path = pathlib.Path(arm_tag(supply, flora)) / "input.json"
                write_json_file(output_path, cfg, indent=1)
                print(f"{arm_tag(supply, flora)}: supply={supply:g}x "
                      f"flora={flora:g}x -> {out_root / output_path}")
    finally:
        os.chdir(previous_cwd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-config", type=pathlib.Path, required=True)
    parser.add_argument("--out-root", type=pathlib.Path, required=True)
    args = parser.parse_args()
    main(args.base_config, args.out_root)
