#!/usr/bin/env python3
"""Section 8 follow-up: is the in-vitro exclusion recoverable at any corrinoid level?

The first pass FAILED V1 and did so in a specific, mechanistic way.  In a
well-mixed 100 um cube at 1e10 cells/mL, the ColE2 producer lysed 298-345 cells
to buy 35-40 receptor-mediated kills -- 0.12 kills per producer death -- so the
producer strain finished *below* its own paired null.  Colicin is not merely
weak here; the strategy is net-negative by about 8x.

That is not a diffusion problem.  ColE2's mucin retardation is
`retardation_from_pI(6.5)` = 1.2 + 60/10^1.85 = 2.04, not the 50 that ColE1 (pI
9.0) carries, so the nuclease colicin is almost unretarded.  A point burst of
5e4 molecules over `burst_release_tau = 300 s` gives a steady-state
`Q/(4*pi*D_eff*r)` of ~1.3e-6 mol/m^3 at 1 um and ~2.6e-7 at 5 um, i.e. at or
above `kd_colicinE_btuB = 5e-7` across a plume containing hundreds of cells.
Two orders of magnitude of kills are missing against that estimate.

They are missing into corrinoid competition.  `FixReceptor::toxin_occupancy`
is competitive:

    apparent_kd = kd_colicinE_btuB * (1 + [corrinoid] / kd_b12_btuB)

and the shipped field is `b12_initial_conc = 1e-3 mol/m^3` (1 uM) against
`kd_b12_btuB = 1e-6` (1 nM).  BtuB is therefore corrinoid-saturated 1000-fold
and the apparent Kd for colicin E is 5e-4 mol/m^3, not 5e-7.  The header
comment on that constant already names it "the key unknown governing colicin-E
competition against the ~1 uM corrinoid pool", and the Spec 14 brief lists
corrinoid availability as one of the four things to inspect if the reversal is
absent.  It is the whole of the absence.

This ladder asks the question that follows: at what corrinoid concentration,
if any, does the producer win in the well-mixed arm?  That matters because the
in-vitro half of the reversal is run in medium, and medium is not 1 uM free
corrinoid -- so a producer that wins in vitro and coexists in vivo may be a
corrinoid contrast rather than a spatial one.

Confound handled by design, not by argument: corrinoid is also a growth
substrate (`km_b12 = 1e-6`), so lowering it slows *both* strains.  Every level
carries its own no-producer null with identical founders and seeds, and the
reported statistic is always the paired producer-minus-null difference, so a
shared growth penalty cancels.  Nothing is tuned to make the producer win; the
levels bracket the shipped value from 1000x above kd_b12 down to below it.
"""
from __future__ import annotations

import json
import pathlib

from prepare import COMMON, PRODUCER, SEEDS, SENSITIVE, VITRO

ROOT = pathlib.Path(__file__).resolve().parent

# mol/m^3.  1e-3 is the shipped value (1 uM, 1000x kd_b12_btuB); 1e-6 is kd
# itself; 1e-9 is a corrinoid-poor medium.
B12_LEVELS = [1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-9]


def tag(level: float) -> str:
    return f"{level:.0e}".replace("-", "m").replace("+", "")


def main() -> None:
    for level in B12_LEVELS:
        for treatment in ("producer", "null"):
            producer = (dict(PRODUCER) if treatment == "producer"
                        else dict(PRODUCER, plasmids=[]))
            for seed in SEEDS:
                cfg = dict(COMMON)
                cfg.update(VITRO)
                cfg["seed"] = seed
                cfg["b12_initial_conc"] = level
                cfg["initial_strains"] = [producer, dict(SENSITIVE)]
                name = f"b12{tag(level)}_{treatment}_s{seed}"
                out = ROOT / "arms_corrinoid" / name
                out.mkdir(parents=True, exist_ok=True)
                cfg["_comment"] = [
                    "Section 8 follow-up: corrinoid competition ladder.",
                    "Well-mixed in-vitro analogue; only b12_initial_conc and",
                    "the presence of the ColE2 plasmid differ across arms.",
                    "apparent_kd(colicin E) = 5e-7 * (1 + b12/1e-6).",
                    f"arm = {name}.",
                ]
                (out / "input.json").write_text(json.dumps(cfg, indent=1) + "\n")


if __name__ == "__main__":
    main()
