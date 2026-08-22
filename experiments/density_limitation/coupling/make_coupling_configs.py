"""Generate the agent_carbon_coupling probe configs for the amended Spec 12
Change 1 (merged in #312).

Lead-authored. Purpose is measurement, not calibration.

Why this probe exists
---------------------
Spec 12 proposes `agent_carbon_coupling = 1e-16 mol/s/agent` and a sweep
{1e-17, 1e-16, 1e-15}. At `grid_dx = 2 um` the voxel is 8 fL, the entire
background VBF carbon sink in one voxel is 5.5e-5 * 8e-18 = 4.4e-22 mol/s, and
measured per-agent carbon demand is ~5e-20 mol/s, so the proposed values are
2e3-2e5x per-agent demand and every proposed arm sits above total voxel
starvation. The demand-anchored sweep below spans 2%-200% of per-agent demand
instead (see docs/SPEC12_DENSITY_LIMITATION.md).

What is being read, and why it is not the coupling constant
-----------------------------------------------------------
`agent_carbon_coupling` is a model constant, not an observable: the same value
means different competition at different densities, geometries and supply
rates. The observable is the nutrient blocking fraction added in #312,

    (agent_uptake + realized_maintenance)
    / (agent_uptake + realized_maintenance + vbf_sink)

per summary interval and species. The claim this probe tests is that coupling
moves blocking fraction monotonically, and that carrying capacity tracks
blocking fraction rather than the constant. Everything downstream (the RPS /
coexistence campaign) should report against blocking fraction.

Arms (seed 1001, delivery-limited uptake, guard ON, colicin silent):
  c0     - coupling 0. The exact historical uniform sink; the sensitivity
           control without which a coupling effect proves nothing.
  c1e21  - 2% of measured per-agent demand.
  c1e20  - 20%.
  c1e19  - 200%. Expected to be the arm that starves occupied voxels; a
           `closure_violation` halt here is the #311 gate working, not a bug,
           and is a result rather than a failure of the probe.

Delivery mode is required (realized removal must fund growth, otherwise the
coupling changes the field without changing biology) and delivery mode rejects
`gpu_enabled`, so these are CPU arms and must run sequentially: one arm is
~2.7 GB RSS on a ~7.9 GB box.

Base config is the merged probe geometry
(inputs_probe/<sha>/full_f100/input.json) at 1.00x J_dir, so the population
trajectory is comparable to the arms already analysed in ../probe/.
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
    write_json_file = path_utils.write_json_file

SEED = 1001
HORIZON_S = 21600  # 6 h; blocking fraction is a per-interval diagnostic
PER_AGENT_DEMAND = 5.0e-20  # mol/s, measured in the #310 validation run
# (tag, agent_carbon_coupling, fraction of measured per-agent demand)
ARMS = [
    ("c0", 0.0, 0.0),
    ("c1e21", 1.0e-21, 0.02),
    ("c1e20", 1.0e-20, 0.20),
    ("c1e19", 1.0e-19, 2.00),
]


def main(base_config: pathlib.Path, out_root: pathlib.Path) -> None:
    base = json.loads(validate_input_path(base_config).read_text())
    out_root = out_root.expanduser().resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    previous_cwd = pathlib.Path.cwd()
    os.chdir(out_root)

    try:
        for tag, coupling, demand_fraction in ARMS:
            cfg = dict(base)
            cfg["total_time"] = HORIZON_S
            cfg["seed"] = SEED
            cfg["gpu_enabled"] = False
            cfg["uptake_limit"] = "delivery"
            cfg["vbf_agent_carbon_coupling"] = coupling
            cfg["hdf5_file"] = "output.h5"
            cfg["_comment"] = [
                "agent_carbon_coupling probe for merged #312 (Spec 12 Change 1).",
                "Measurement run, not calibration: nothing here is being fitted.",
                (f"vbf.agent_carbon_coupling = {coupling:.3e} mol/s/agent = "
                 f"{demand_fraction:.2f}x the measured per-agent carbon demand "
                 f"{PER_AGENT_DEMAND:.1e} mol/s. Spec 12's own 1e-16 is ~2000x "
                 "full per-agent demand at grid_dx=2um and is not swept here."),
                ("uptake_limit=delivery so growth is funded only by realized "
                 "field removal; delivery rejects gpu_enabled, so this is a CPU "
                 "arm and arms must run sequentially (~2.7 GB RSS each)."),
                ("Read per arm: nutrient blocking fraction (the observable), "
                 "population trajectory and plateau, funded fraction, "
                 "maintenance shortfall, carbon in occupied vs empty voxels, "
                 "and /run_provenance termination cause."),
                ("A closure_violation halt in the highest arm is the expected "
                 "signature of voxel starvation, not a defect."),
            ]

            output_path = pathlib.Path(tag) / "input.json"
            write_json_file(output_path, cfg, indent=1)
            print(f"{tag}: agent_carbon_coupling={coupling:.3e} "
                  f"({demand_fraction:.2f}x demand) -> {out_root/output_path}")
    finally:
        os.chdir(previous_cwd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--base-config",
        required=True,
        help="probe full_f100 input.json to inherit geometry and supply from")
    parser.add_argument(
        "--out-root",
        default=os.environ.get("GUTIBM_COUPLING_ROOT", "."),
        help="directory to write <arm>/input.json into")
    args = parser.parse_args()
    main(pathlib.Path(args.base_config), pathlib.Path(args.out_root))
