#!/usr/bin/env python3
"""Is the flora the thief, and what does fermentative maintenance actually cost?

Why this experiment exists
--------------------------
#327 measured that under a depletable (Robin) epithelium the background flora
takes 99.0-99.1% of everything delivered, independent of the transfer
coefficient and of founder count, while the oxygen field stayed uniform to
within 2.3% across the whole 300 um domain. Two explanations were confounded in
that measurement, and #328 showed why:

  * The VBF oxygen sink was integrated *explicitly* (``conc += reac*dt``) and
    applied before the implicit solve that integrates the agents' sinks, so the
    flora consumed from the undepleted pool and the agents drew on the residual.
  * An explicit first-order sink is stable only for ``k*dt < 1``. At
    ``bio_dt = 60 s`` any sink strong enough to set a penetration depth is far
    outside that bound, so the "consumption-consistent" 0.84 1/s arm in #325
    clipped 98% of its demand. I read that as the field having no oxygen to
    give; it was the integrator.

#328 moves the flora sink into the same implicit solve and splits realized
removal by rate share, so both defects are gone and two questions that were
previously unaskable become measurable.

  Q1. Attribution. Is the flora the binding constraint on agent oxygen, or is
      it voxel-scale resupply? The ``s0`` arms answer this directly: with the
      flora sink switched off entirely, whatever the agents still cannot fund is
      not attributable to competition.
  Q2. Geometry. With a sink that can now be integrated at any magnitude, does a
      1/e penetration depth actually emerge, and where? ``sqrt(D/k)`` with the
      repository's own ``oxygen.D_free = 2.1e-9`` predicts 50 um at 0.84 1/s and
      25 um at 3.36 1/s. That prediction has never been tested, only asserted.

And one that has been open since #325:

  Q3. What does paying maintenance fermentatively cost? ``anaerobic_maintenance
      _factor`` shipped at 15 with no citation, and at 15 every arm of #327 died
      as a fermenter inside 7 h, so nothing measured so far is a survival
      result.

Stages (run ``ladder`` first; it determines which ``oxygen`` arms are worth it)
------------------------------------------------------------------------------
``--stage ladder`` -- the maintenance axis alone, at the oxygen configuration of
    #327's ``robin`` arm (k = 1.2e-6 m/s, C_tissue = 5.0e-2 mol/m^3, shipped
    ``vbf_sink`` = 1e-3 1/s, 80 founders). Four arms:

      amf1p0  -- 1.0. Maintenance carbon costs the same fermentatively as
                 aerobically. This is the value called for on the grounds that
                 NGAM is ATP-invariant (Varma & Palsson 1994).
      amf4p1  -- 4.1. The factor already used for *growth* carbon
                 (``anaerobic_carbon_cost_factor``). If maintenance carbon
                 should scale by the anaerobic ATP-yield penalty, this is the
                 value that keeps one yield ratio in the model instead of two.
      amf8p3  -- 8.3. The endpoint of the same derivation from Palsson's W3110
                 numbers (NGAM/26 aerobic vs NGAM/3 anaerobic). This is the arm
                 that separates "maintenance deserves its own stoichiometry"
                 from "maintenance shares growth's".
      amf15   -- 15, the shipped value. Not a candidate: it is the reproduction
                 control. Every result we have (the supply ladder, the placement
                 matrix, #327) was measured at 15, so without this arm a change
                 in survival cannot be separated from the new sink code in #328.

    The maintenance factor acts monotonically on carbon cost, so what this stage
    measures is a threshold, not four independent points.

``--stage oxygen --amf A [B ...]`` -- the oxygen axis, crossed only with the
    maintenance values that bracket the threshold the ladder finds. Per ``--amf``
    value, four sink arms:

      s0     -- ``vbf_sink`` = 0. No flora oxygen consumption at all. The
                attribution arm for Q1, and the only arm in which the agents'
                funding ceiling is set by supply and transport alone.
      s1e-3  -- the shipped value, i.e. #327's configuration re-run on the #328
                sink. It must reproduce #327's steady-state zero funded growth
                O2 at the same maintenance factor; if it does not, #328 changed
                more than the integration order. Note this arm is *not* expected
                to be bit-identical to #327 for the Dirichlet control: with the
                z-gradient on, the old explicit sink consumed against the
                perturbation only and never saw the profile-borne oxygen.
      s0p84  -- D/(50 um)^2. Predicted 1/e depth 50 um.
      s3p36  -- D/(25 um)^2. Predicted 1/e depth 25 um, i.e. the shell the model
                used to impose as a hard-coded constant.

    Plus one arm independent of ``--amf``:

      ceiling -- a single founder at ``vbf_sink`` = 0 and the sourced boundary.
                One cell competing with nothing establishes the per-cell funding
                ceiling this geometry can support. If even this arm cannot fund
                growth respiration, no density and no flora setting will, and
                the limit is the apical flux or the voxel-scale transport.

    The two high-sink values are *hypotheses about geometry*, not sourced
    parameters, and are labelled as such. They are derived from the depth the
    literature reports for dense mucus, inverted through this model's own
    diffusivity. If a plausible depth requires a flora sink that starves the
    agents even with proportional competition, that is the finding.

Held fixed across every arm
---------------------------
Seed, horizon (24 h), geometry, carbon supply (1.00x J_dir), flora carbon
(1.00x), guard, the anatomy-derived placement of #324, the ``funded``
respiration driver, delivery-limited carbon and oxygen uptake, and the Robin
boundary at the sourced k = 1.2e-6 m/s with C_tissue = 5.0e-2 mol/m^3. The
imposed z-gradient is off everywhere: the profile is predicted, never restored.

Read per arm
------------
Cumulative funded vs demanded growth O2 and the funded fraction; funded
maintenance O2 against maintenance shortfall; the agent share and the flora
share of realized oxygen removal (these must sum to the total; #328 asserts it
per voxel and this checks it per run); the emergent surface concentration and
the depth at which oxygen falls to 1/e of it, against the 25 um the model used
to impose and against ``sqrt(D/k)``; whether funded respiration persists past
the opening hours or is an initialization transient; fermentation fraction;
acetate; N trajectory, termination cause, and z-distribution of survivors;
oxygen and carbon reaction clips, which must stay at 0 -- a nonzero oxygen clip
at high sink means #328 did not actually remove the stability bound.

Predictions, fixed before running
---------------------------------
1. ``s0p84`` and ``s3p36`` produce zero oxygen reaction clip. If they clip,
   the sink is still being integrated explicitly somewhere.
2. A 1/e crossing appears in both high-sink arms, near the ``sqrt(D/k)`` depth.
   This is the first geometry in which an oxic-anoxic transition can exist.
3. Raising the flora sink 300-3400x *reduces* the agents' share of delivery even
   with proportional competition, because the split is by rate and the flora's
   rate is what grew. So the high-sink arms should fund respiration for the
   agents nearest the epithelium and less than the shipped arm elsewhere. A
   niche, not a uniform improvement.
4. ``s0`` funds the largest share of demand of any arm. If it still cannot fund
   growth respiration, competition was never the constraint and Q1 resolves
   against the flora-as-thief reading.
5. The maintenance ladder has a survival threshold somewhere in it, which is the
   point of running four values rather than adopting one.

Nothing here is fitted. No arm is tuned to make respiration appear.

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

# Robin boundary, held fixed. k is measured rat mucus+mucosa permeability;
# C_TISSUE is tissue-side O2. Both sourced (#327).
K_SOURCED = 1.2e-6
C_TISSUE = 5.0e-2

# Oxygen diffusivity as this repository ships it. Used only to invert a target
# penetration depth into a first-order rate constant, k = D / delta^2.
D_OXYGEN = 2.1e-9

# Flora oxygen sink values, 1/s. 0 is the attribution arm; 1e-3 is shipped;
# the two high values are geometry hypotheses, not sourced parameters.
SINK_OFF = 0.0
SINK_SHIPPED = 1.0e-3
SINK_DELTA_50UM = D_OXYGEN / (50.0e-6 ** 2)
SINK_DELTA_25UM = D_OXYGEN / (25.0e-6 ** 2)

SINK_ARMS = {
    "s0": SINK_OFF,
    "s1e-3": SINK_SHIPPED,
    "s0p84": SINK_DELTA_50UM,
    "s3p36": SINK_DELTA_25UM,
}

# Maintenance ladder. 15 is the shipped value and the reproduction control.
LADDER = {
    "amf1p0": 1.0,
    "amf4p1": 4.1,
    "amf8p3": 8.3,
    "amf15": 15.0,
}

# Maintenance factor for the ladder stage's oxygen configuration is the arm
# variable; for the ceiling arm it is pinned to the cheapest value so that a
# failure to fund respiration cannot be blamed on carbon maintenance.
CEILING_AMF = 1.0

# Anatomy-derived placement (#324). Held fixed: #325 showed placement does not
# affect funding.
ANATOMIC = {
    "initial_population.placement": "anatomic",
    "initial_population.anatomic_exclusion_floor": 20.0e-6,
    "initial_population.anatomic_exponential_scale": 40.0e-6,
    "initial_population.anatomic_outer_extent": 150.0e-6,
}


def single_founder(base: dict) -> list:
    """One founder of the first strain, for the funding-ceiling arm."""
    strains = json.loads(json.dumps(base["initial_strains"]))
    first = strains[0]
    first["count"] = 1
    return [first]


def common_config(base: dict, vbf_sink: float, amf: float) -> dict:
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
    cfg["oxygen.epithelial_boundary"] = "robin"
    cfg["oxygen.epithelial_conc"] = C_TISSUE
    cfg["oxygen.epithelial_transfer_coeff"] = K_SOURCED
    cfg["oxygen.z_gradient"] = False
    cfg["oxygen.vbf_sink"] = vbf_sink
    cfg["oxygen.anaerobic_maintenance_factor"] = amf

    cfg["acetate.enabled"] = True
    cfg["metabolism.acid_inhibition_enabled"] = True

    for key in ("initial_population.placement",
                "initial_population.z_min",
                "initial_population.z_max"):
        cfg.pop(key, None)
    cfg.pop("initial_population", None)
    cfg.update(ANATOMIC)

    cfg["hdf5_file"] = "output.h5"
    cfg["hdf5"] = {
        "schedule": {
            "summary": 1,
            "agents": 5,
            # 60 summary steps = 1 h simulated. Grid output is what makes the
            # emergent O2 profile and its 1/e depth measurable.
            "grid": 60,
            "grid_species": ["oxygen", "carbon", "acetate"],
            "lineage": 0,
            "genome": 0,
        }
    }
    return cfg


def provenance(tag: str, vbf_sink: float, amf: float, founders: int) -> list:
    notes = [
        ("Oxygen competition measurement run on the #328 implicit, "
         "rate-proportional flora sink. Nothing here is fitted or calibrated."),
        (f"arm = {tag}. oxygen.vbf_sink = {vbf_sink:.4g} 1/s. "
         f"oxygen.anaerobic_maintenance_factor = {amf:.3g}. "
         f"{founders} founder(s)."),
        (f"Robin boundary held fixed at k = {K_SOURCED:.3g} m/s, "
         f"C_tissue = {C_TISSUE:.3g} mol/m^3, both sourced (#327). The imposed "
         "z-gradient is off: the O2 profile is predicted, not restored."),
    ]
    if vbf_sink == SINK_OFF:
        notes.append(
            "vbf_sink = 0 is the attribution arm: with no flora oxygen "
            "consumption, anything the agents still cannot fund is not "
            "attributable to competition.")
    elif vbf_sink == SINK_SHIPPED:
        notes.append(
            "vbf_sink = 1e-3 is the shipped value and the reproduction control "
            "for #327 on the new sink code. It is not expected to be "
            "bit-identical where a z-gradient was active: the old explicit "
            "sink consumed against the perturbation only, never the "
            "profile-borne oxygen (the same error #310 fixed for carbon).")
    else:
        depth_um = (D_OXYGEN / vbf_sink) ** 0.5 * 1.0e6
        notes.append(
            f"vbf_sink = {vbf_sink:.4g} 1/s is a geometry hypothesis, not a "
            f"sourced parameter: it is D/delta^2 with D = {D_OXYGEN:.3g} m^2/s "
            f"for a target 1/e depth of {depth_um:.0f} um. Before #328 this "
            "value could not be integrated at bio_dt = 60 s at all "
            "(k*dt >> 1), which is why #325 measured 98% of its demand "
            "clipped. The oxygen reaction clip must now be 0; if it is not, "
            "the sink is still explicit somewhere.")
    notes.append(
        "anaerobic_maintenance_factor: 15 is the shipped, uncited value and "
        "the reproduction control; 4.1 matches the factor already applied to "
        "growth carbon; 8.3 is the endpoint of the Palsson W3110 ATP-yield "
        "derivation; 1.0 asserts maintenance carbon is yield-invariant. The "
        "factor is monotone in carbon cost, so the ladder measures a survival "
        "threshold, not four independent points.")
    notes.append(
        "Read per arm: cumulative funded vs demanded growth O2, funded "
        "maintenance O2 against shortfall, agent share vs flora share of "
        "realized O2 removal (must sum to the total), emergent surface "
        "concentration and 1/e depth against sqrt(D/k) and against the "
        "retired 25 um constant, persistence vs transient of any funded "
        "respiration, fermentation fraction, acetate, N trajectory, survivor "
        "z-distribution, oxygen and carbon reaction clips (must stay 0), and "
        "/run_provenance termination cause.")
    return notes


def build_ladder_arm(base: dict, tag: str) -> dict:
    amf = LADDER[tag]
    cfg = common_config(base, SINK_SHIPPED, amf)
    cfg["initial_strains"] = json.loads(json.dumps(base["initial_strains"]))
    founders = sum(s["count"] for s in cfg["initial_strains"])
    cfg["_comment"] = provenance(tag, SINK_SHIPPED, amf, founders)
    return cfg


def build_oxygen_arm(base: dict, sink_tag: str, amf: float) -> dict:
    cfg = common_config(base, SINK_ARMS[sink_tag], amf)
    cfg["initial_strains"] = json.loads(json.dumps(base["initial_strains"]))
    founders = sum(s["count"] for s in cfg["initial_strains"])
    tag = f"{sink_tag}_amf{amf:g}"
    cfg["_comment"] = provenance(tag, SINK_ARMS[sink_tag], amf, founders)
    return cfg


def build_ceiling_arm(base: dict) -> dict:
    cfg = common_config(base, SINK_OFF, CEILING_AMF)
    cfg["initial_strains"] = single_founder(base)
    cfg["_comment"] = provenance("ceiling", SINK_OFF, CEILING_AMF, 1) + [
        ("One founder, no flora oxygen sink, cheapest maintenance: the "
         "per-cell funding ceiling this geometry can support. If growth "
         "respiration is unfunded here, no density or flora setting will fund "
         "it and the limit is the apical flux or voxel-scale transport."),
    ]
    return cfg


def arm_configs(base: dict, stage: str, amfs: list) -> dict:
    if stage == "ladder":
        return {tag: build_ladder_arm(base, tag) for tag in LADDER}
    arms = {}
    for amf in amfs:
        for sink_tag in SINK_ARMS:
            arms[f"{sink_tag}_amf{amf:g}"] = build_oxygen_arm(
                base, sink_tag, amf)
    arms["ceiling"] = build_ceiling_arm(base)
    return arms


def main(base_config: pathlib.Path, out_root: pathlib.Path, stage: str,
         amfs: list) -> None:
    if stage == "oxygen" and not amfs:
        raise SystemExit(
            "--stage oxygen requires --amf, i.e. the maintenance values "
            "bracketing the survival threshold the ladder stage measured")
    base = json.loads(validate_input_path(base_config).read_text())
    out_root = out_root.expanduser().resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    out_root = validate_output_path(out_root / "input.json").parent
    previous_cwd = pathlib.Path.cwd()
    os.chdir(out_root)
    try:
        for tag, cfg in arm_configs(base, stage, amfs).items():
            path = pathlib.Path(tag) / "input.json"
            write_json_file(path, cfg, indent=1)
            founders = sum(s["count"] for s in cfg["initial_strains"])
            print(f"{tag}: vbf_sink={cfg['oxygen.vbf_sink']:.4g} "
                  f"amf={cfg['oxygen.anaerobic_maintenance_factor']:g} "
                  f"founders={founders} -> {out_root / path}")
    finally:
        os.chdir(previous_cwd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-config", type=pathlib.Path, required=True)
    parser.add_argument("--out-root", type=pathlib.Path, required=True)
    parser.add_argument("--stage", choices=("ladder", "oxygen"),
                        required=True)
    parser.add_argument("--amf", type=float, nargs="*", default=[],
                        help="maintenance factors for --stage oxygen")
    args = parser.parse_args()
    main(args.base_config, args.out_root, args.stage, args.amf)
