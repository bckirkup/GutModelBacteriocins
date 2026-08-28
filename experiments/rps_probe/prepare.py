#!/usr/bin/env python3
"""Generate the btuB-null resistant-strain probe arms.

Two questions, both authored against the shipped code paths:

  A (efficacy) Does a founder strain declared with `receptor_expression:
      {"BtuB": 0.0}` actually survive a ColE1 producer that kills the
      otherwise-identical sensitive strain?  Kill hazard is
      `kill_rate * occupancy * immunity`, and occupancy is linear in
      `agent.receptor_expr[BtuB]` (src/fixes/fix_receptor.cpp), so the
      prediction is zero receptor-mediated deaths for type 3.

  B (cost) What does that knockout cost?  Growth enters through
      `Km_b12 = km_b12 / (expr_btuB * lig_aff)` with `expr_btuB` clamped
      at 0.01 (src/fixes/fix_metabolism.cpp), i.e. a 100x Km inflation.
      At the shipped `b12_initial_conc = 1e-3` against `km_b12 = 1e-6`
      the knockout is still B12-saturated, so the predicted cost is
      small; the sweep asks at what corrinoid level the cost becomes a
      real fitness difference.  No producer in these arms.

Nothing here is fitted.  Only the keys named in each arm label differ.
Carbon amplitude is x2 the shipped default because at x1 the patch
capacity measured in docs/CARBON_LADDER_CAMPAIGN.md is 122 agents from
100 founders -- there is no growth headroom in which a fitness
difference could express itself.
"""
from __future__ import annotations

import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent

BASE = {
    "total_time": 43200.0,
    "bio_dt": 60,
    "output_interval": 600,
    "domain_x": 0.0002,
    "domain_y": 0.0002,
    "domain_z": 0.0001,
    "grid_dx": 4e-06,
    "dysbiosis_threshold": 1.0e10,
    "mucus_thickness": 0.0001,
    "radial_turnover": 5400,
    "distal_transit": 43200,
    "vbf_density": 1.0e11,
    "vbf_viscosity": 0.01,
    "vbf_mucin_liberation": 5e-05,
    "vbf_carbon_sink_vmax": 5.5e-05,
    "carbon.boundary_conc": 0.01,
    "carbon_z_gradient": True,
    "carbon_z_lambda": 2.5e-05,
    "metabolism.uptake_limit": "delivery",
    "oxygen.k_ROS": 0.0,
    "initial_population.placement": "z_slab",
    "initial_population.z_min": 0.0,
    "initial_population.z_max": 2e-05,
    "hdf5_file": "output.h5",
    "hdf5": {
        "schedule": {
            "summary": 1,
            "agents": 10,
            "grid": 0,
            "lineage": 0,
            "genome": 0,
        }
    },
}

# type 1 = ColE1 producer (C), type 2 = sensitive (S), type 3 = btuB-null (R).
PRODUCER = {"type": 1, "count": 40, "mu_max": 5.5e-4,
            "plasmids": ["ColE1"], "conjugative": False}
SENSITIVE = {"type": 2, "count": 40, "mu_max": 5.5e-4,
             "plasmids": [], "conjugative": False}
RESISTANT = {"type": 3, "count": 40, "mu_max": 5.5e-4,
             "plasmids": [], "conjugative": False,
             "receptor_expression": {"BtuB": 0.0}}

SEEDS_A = [20260901, 20260902, 20260903]
SEEDS_B = [20260901, 20260902]
B12_LEVELS = [1e-3, 1e-4, 3e-5, 1e-5]


def emit(name: str, cfg: dict, note: str) -> None:
    out = ROOT / name
    out.mkdir(parents=True, exist_ok=True)
    cfg = dict(cfg)
    cfg["_comment"] = [
        "btuB-null resistant-strain probe; see docs/SPEC13_LYSIS_SELECTION.md.",
        "type 1 = ColE1 producer, type 2 = sensitive, type 3 = btuB-null.",
        "Read per arm: per-type agent counts from /agents/<step>/type, the",
        "population ledger, and /run_provenance termination cause.",
        f"arm = {name}. {note}",
    ]
    (out / "input.json").write_text(json.dumps(cfg, indent=1) + "\n")


def main() -> None:
    for seed in SEEDS_A:
        cfg = dict(BASE)
        cfg["seed"] = seed
        cfg["initial_strains"] = [PRODUCER, SENSITIVE, RESISTANT]
        emit(f"A_three_strain_s{seed}", cfg, "three-strain efficacy arm")

    # Null arm: same three founders, no producer plasmid anywhere.  If S and R
    # still diverge here, the divergence is not colicin.
    for seed in SEEDS_A[:1]:
        cfg = dict(BASE)
        cfg["seed"] = seed
        cfg["initial_strains"] = [
            dict(PRODUCER, plasmids=[]), SENSITIVE, RESISTANT]
        emit(f"A_null_no_producer_s{seed}", cfg, "null arm, no ColE1")

    for level in B12_LEVELS:
        tag = f"{level:.0e}".replace("-", "m").replace("+", "")
        for seed in SEEDS_B:
            cfg = dict(BASE)
            cfg["seed"] = seed
            cfg["b12_initial_conc"] = level
            cfg["initial_strains"] = [dict(SENSITIVE, count=60),
                                      dict(RESISTANT, count=60)]
            emit(f"B_b12_{tag}_s{seed}", cfg,
                 f"resistance-cost arm, b12_initial_conc={level}")


if __name__ == "__main__":
    main()
