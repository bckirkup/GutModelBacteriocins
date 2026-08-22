"""Generate the mechanism probe configs for #302 (Pirt carbon maintenance) and
#304 (O2 metabolic mode + acid inhibition), both now merged.

Lead-authored. Purpose is measurement, not calibration: read what the new
mechanisms actually do at campaign geometry before re-bracketing delivery.

Anchors (see ../spec12_audit.md in the campaign checkout):
  Pirt aerobic maintenance  ~0.04 g glucose/gDW/h -> 2.1e-5 mol/(s kg wet)
  oxygen.anaerobic_maintenance_factor default 15  -> 3.15e-4 mol/(s kg),
    matching the independent anaerobic estimate 3.2e-4.
  Measured mean agent biomass 8.14e-16 kg (Sherwood final agent dump), so the
  per-agent anaerobic draw is ~2.56e-19 mol/s.

Predicted capacity, N* = supply / (rate_eff * biomass), with the epithelial
area 9.216e-9 m^2:
  0.18x J_dir -> supply 1.784e-17 mol/s -> N* ~ 70 agents
  0.30x       -> 2.97e-17              -> N* ~ 116
  0.60x       -> 5.95e-17              -> N* ~ 232
  1.00x       -> 9.91e-17              -> N* ~ 387
The dysbiosis guard at 1e8 cells/mL is ~276 agents in this 2.765e-12 m^3
domain, so the prediction is that the brake binds *below* the guard up to
~0.7x and above it thereafter. That is the graded claim this probe tests: a
few different delivery rates should give a few different plateau densities,
and the low arms should plateau rather than halt.

Arms (seed 1001, 24 h, guard ON, kd_b12_btuB at the 1e-6 default so these are
comparable to every arm already analysed):
  off   - both mechanisms disabled, 0.18x. Reproduces the prior bloom; the
          sensitivity control without which a plateau proves nothing.
  maint - carbon maintenance only (aerobic anchor, no mode switch), 0.18x.
          Isolates the linear-in-N carbon sink from the O2/acid machinery.
  full  - maintenance + O2 metabolic mode + acetate, at 0.18/0.30/0.60/1.00x.

Read per arm: mean realized fermentation fraction (is the aerobic branch
masked, as predicted, or inert?), realized carbon cost, maintenance draw and
maintenance shortfall, agent uptake vs epithelial delivery, acetate
concentration, population trajectory and guard crossing time.

The run root contains corrinoid/corr_kd1e4_s1001.json as the source
configuration; generated configs are written under probe/.
"""

import argparse
import json
import os
import pathlib

J_DIR = 1.0756e-8          # mol/m^2/s, measured Dirichlet face flux
PIRT_AEROBIC = 2.1e-5      # mol C/(s kg biomass)
SEED = 1001
HORIZON_S = 86400          # 24 h
# (tag, flux multiplier, maintenance on, metabolic mode on)
ARMS = [
    ("off_f018", 0.18, False, False),
    ("maint_f018", 0.18, True, False),
    ("full_f018", 0.18, True, True),
    ("full_f030", 0.30, True, True),
    ("full_f060", 0.60, True, True),
    ("full_f100", 1.00, True, True),
]


def main(run_root: pathlib.Path) -> None:
    base = json.loads(
        (run_root / "corrinoid" / "corr_kd1e4_s1001.json").read_text())
    out = run_root / "probe"
    out.mkdir(parents=True, exist_ok=True)

    for tag, mult, maintenance, mode in ARMS:
        cfg = dict(base)
        cfg.pop("notes", None)
        cfg["total_time"] = HORIZON_S
        cfg["seed"] = SEED
        cfg["carbon.epithelial_boundary"] = "flux"
        cfg["carbon.epithelial_flux"] = J_DIR * mult
        cfg["kd_b12_btuB"] = 1.0e-6

        cfg["metabolism.carbon_maintenance_rate"] = (
            PIRT_AEROBIC if maintenance else 0.0)

        cfg["oxygen.enabled"] = bool(mode)
        cfg["oxygen.metabolic_switch_enabled"] = bool(mode)
        cfg["acetate.enabled"] = bool(mode)
        cfg["metabolism.acid_inhibition_enabled"] = bool(mode)

        cfg["_comment"] = [
            "Mechanism probe for merged #302 (Pirt carbon maintenance) and",
            "#304 (O2 metabolic mode + undissociated-acetate inhibition).",
            "Measurement run, not calibration: nothing here is being fitted.",
            f"carbon.epithelial_flux = {mult:.2f}x J_dir ({J_DIR:.4e}).",
            (f"metabolism.carbon_maintenance_rate = "
             f"{cfg['metabolism.carbon_maintenance_rate']:.3e} mol C/(s kg), "
             "the Pirt aerobic anchor (0.04 g glucose/gDW/h); the anaerobic "
             "draw follows from oxygen.anaerobic_maintenance_factor = 15."),
            "Metabolic mode " + ("ON" if mode else "OFF")
            + ": aerobic/anaerobic carbon cost factors 1.0/4.1 (cost is "
              "substrate-per-biomass here, so fermentation is the expensive "
              "mode - the spec's 0.25 had it inverted).",
            ("Dysbiosis guard stays ON at 1e8 cells/mL (~276 agents in this "
             "2.765e-12 m^3 domain)."),
            ("Predicted plateau N* = supply/(rate_eff*biomass) with rate_eff "
             "3.15e-4 anaerobic and biomass 8.14e-16 kg."),
            ("kd_b12_btuB at the 1e-6 default (colicin silent) so the density "
             "trajectory is comparable to the bracket and Sherwood arms."),
        ]

        name = f"probe_{tag}_s{SEED}.json"
        (out / name).write_text(json.dumps(cfg, indent=1) + "\n")
        print(f"{name}: flux={cfg['carbon.epithelial_flux']:.5e} "
              f"maint={cfg['metabolism.carbon_maintenance_rate']:.3e} "
              f"mode={mode}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--run-root",
        default=os.environ.get("GUTIBM_CAMPAIGN_ROOT", "."),
        help="campaign root containing corrinoid/ and probe/ (default: "
        "GUTIBM_CAMPAIGN_ROOT or current directory)",
    )
    main(pathlib.Path(parser.parse_args().run_root).expanduser().resolve())
