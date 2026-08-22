"""Generate the epithelial-delivery bracket configs for a guard-safe 168 h run.

Lead-authored. The only quantities that vary across arms are
`carbon.epithelial_flux` and the seed. Everything else is copied verbatim from
the corrinoid-sweep arm (itself the flux-0.3x VBF-sweep arm) so the bracket is
comparable to every run already analysed.

Measured Dirichlet face flux J_dir = 1.0756e-8 mol/m^2/s (stage-1 controls);
the 0.3x arm used 3.2268e-9. Known endpoints on pre-#300 code, seed-matched:
  0.30x -> blooms, dysbiosis guard halts at 23.4 h
  0.10x -> never blooms, decays to 16 agents (lysis 0.125 > division 0.047 /h)
so the guard-safe steady band, if it exists, is interior to [0.10, 0.30].

The run root contains ``corrinoid/corr_kd1e4_s1001.json`` as the source
configuration; generated configs are written under ``bracket/``.
"""

import argparse
import json
import os
import pathlib

J_DIR = 1.0756e-8
MULTIPLIERS = [0.10, 0.14, 0.18, 0.22, 0.26, 0.30]
# Two colicin regimes, because the corrinoid affinity is unresolved: 1e-6 is
# the shipped default (competitive factor 1001, target occupancy 5.4e-3, no
# realized kills - colicin silent) and 1e-4 is factor 11 (occupancy 0.14,
# susceptibles 20 -> 39, colicin biting). If the guard-safe flux is the same in
# both, the corrinoid choice is separable from delivery calibration.
KD_B12 = [1.0e-6, 1.0e-4]
SEED = 1001
HORIZON_S = 604800  # 168 h

def main(run_root: pathlib.Path) -> None:
    base = run_root / "corrinoid" / "corr_kd1e4_s1001.json"
    out = run_root / "bracket"
    source = json.loads(base.read_text())
    out.mkdir(parents=True, exist_ok=True)

    for mult, kd in ((m, k) for m in MULTIPLIERS for k in KD_B12):
        cfg = dict(source)
        cfg["total_time"] = HORIZON_S
        cfg["seed"] = SEED
        cfg["carbon.epithelial_boundary"] = "flux"
        cfg["carbon.epithelial_flux"] = J_DIR * mult

        cfg["kd_b12_btuB"] = kd

        cfg["_comment"] = [
            "Epithelial delivery bracket for a guard-safe 168 h horizon.",
            (f"carbon.epithelial_flux = {mult:.2f}x J_dir (J_dir = "
             f"{J_DIR:.4e} mol/m^2/s measured in the stage-1 controls)."),
            "Dysbiosis guard stays ON at 1e8 cells/mL: the question is which",
            "delivery rate holds a steady population under the guard for the",
            "full week, not how to evade the guard.",
            "Outcome read per arm: guard crossing time (or none), final",
            "population, mean realized mu, and division vs lysis per agent",
            "per hour in the last 24 h.",
            "Runs post-#300/#301, so retardation is pI-derived (ColE1 50.2,",
            "ColB 1.27, ColE2 2.04) rather than the old hardcoded values.",
            f"kd_b12_btuB = {kd:.0e}: colicin "
            + ("silent (competitive factor 1001)" if kd <= 1e-6
               else "biting (competitive factor 11)") + ".",
        ]
        cfg.pop("notes", None)

        kd_tag = f"kd{kd:.0e}".replace("e-0", "e").replace("-", "")
        name = f"bracket_f{round(mult * 100):03d}_{kd_tag}_s{SEED}.json"
        (out / name).write_text(json.dumps(cfg, indent=1) + "\n")
        print(f"{name}: flux={cfg['carbon.epithelial_flux']:.5e} "
              f"kd={cfg['kd_b12_btuB']:.0e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--run-root",
        default=os.environ.get("GUTIBM_CAMPAIGN_ROOT", "."),
        help="campaign root containing corrinoid/ and bracket/ (default: "
        "GUTIBM_CAMPAIGN_ROOT or current directory)",
    )
    main(pathlib.Path(parser.parse_args().run_root).expanduser().resolve())
