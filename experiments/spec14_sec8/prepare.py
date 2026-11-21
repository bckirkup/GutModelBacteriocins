#!/usr/bin/env python3
"""Spec 14 Section 8: the in-vitro / in-vivo reversal, on current main, no code changes.

Target (Henrot et al. preprint, unreviewed): an E-type endonuclease colicin
producer excludes a lambda-lysogen within ~6 h in vitro, yet both persist at
comparable levels over 10 d in dixenic mice.  Per the spec's own sequencing
note, strain 2 here is a plain sensitive non-producer: the claim under test is
that *spatial structure plus density limitation* produce the reversal, which
does not require phage.

Strains (identical except for the plasmid):
  type 1  producer   plasmids ["ColE2"]  -- the E-type nuclease in
                                            src/genome/plasmid.cpp
                                            (is_nuclease = true, BtuB target,
                                             SOS_LYSIS release)
  type 2  sensitive  plasmids []

Conditions:
  vitro  100 um cube, no mucus turnover, no peristalsis, no crypts, no
         motility, uniform rich carbon.  ~1.2e8 cells/mL at 120 founders, and
         every agent is within `toxin_cutoff` (200 um default) of every other,
         so no agent can be spared by the QSSA support.

         Note what this arm cannot be: a closed batch culture.  Carbon enters
         through Dirichlet walls, so the medium is never exhausted and the
         population is limited only by mechanical crowding -- a first pass at
         a 50 um cube ran to 1.3e5 agents, 1e12 cells/mL, roughly half the
         domain volume as cells.  The in-vitro analogue is therefore run to a
         stationary-analogue density of 1e10 cells/mL, at which the dysbiosis
         guard halts it, and the strain ratio is read at that endpoint.  If
         that endpoint arrives before 6 h, say so rather than quoting a 6 h
         number the run never reached.
  vivo   200x200x100 um mucus patch, 5400 s radial turnover, peristalsis on,
         crypts on, motility on, shipped carbon z-gradient.  Same 120
         founders, i.e. ~3e4 cells/mL: four orders below the vitro arm.  That
         density gap is the intended "density limitation" contrast and is also
         the reason the two arms are not matched on anything else.

Null arms carry the same two types with the producer's plasmid removed.  They
are not optional: the earlier RPS probe (docs/SPEC13_LYSIS_SELECTION.md)
measured a *larger* type-ratio separation in the no-producer null than with a
producer, because ~30 survivors sweep a patch stochastically.  Any separation
in the producer arm has to be read against its own null, within seed.

Five seeds per arm, identical across arms, so every comparison is paired.
Nothing here is fitted; the only keys that differ between conditions are the
ones named in the condition dicts below.
"""
from __future__ import annotations

import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent

SEEDS = [20260901, 20260902, 20260903, 20260904, 20260905]

COMMON = {
    "bio_dt": 60,
    "grid_dx": 4e-06,
    "vbf_density": 1.0e11,
    "vbf_viscosity": 0.01,
    "vbf_mucin_liberation": 5e-05,
    "vbf_carbon_sink_vmax": 5.5e-05,
    "metabolism.uptake_limit": "delivery",
    "oxygen.k_ROS": 0.0,
    # Raised off the 1e8 cells/mL default, which fires before either arm has
    # a measurable ratio, but not disabled: 1e10 cells/mL is the density at
    # which these domains are physically jammed, and past it the run is
    # measuring contact mechanics.  Every arm reports its termination cause.
    "dysbiosis_threshold": 1.0e10,
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

VITRO = {
    "total_time": 21600.0,      # 6 h, the window in which exclusion is claimed
    "output_interval": 300.0,
    "domain_x": 0.0001,
    "domain_y": 0.0001,
    "domain_z": 0.0001,
    "grid_dx": 2e-06,
    "mucus_thickness": 0.0001,
    "radial_turnover": 1.0e9,   # turnover effectively off
    "distal_transit": 1.0e9,
    "peristaltic_enabled": False,
    "crypts_enabled": False,
    "motility.enabled": False,
    "carbon_z_gradient": False,
    "carbon.boundary_conc": 0.05,   # rich medium, 10x the shipped mucus value
    "initial_population.placement": "z_slab",
    "initial_population.z_min": 0.0,
    "initial_population.z_max": 0.0001,
}

VIVO = {
    "total_time": 864000.0,     # 10 d
    "output_interval": 3600.0,
    "domain_x": 0.0002,
    "domain_y": 0.0002,
    "domain_z": 0.0001,
    "mucus_thickness": 0.0001,
    "radial_turnover": 5400,
    "distal_transit": 43200,
    "peristaltic_enabled": True,
    "peristaltic_period": 20.0,
    "peristaltic_amplitude": 0.5,
    "peristaltic_wavelength": 0.001,
    "crypts_enabled": True,
    "motility.enabled": True,
    "carbon_z_gradient": True,
    "carbon.boundary_conc": 0.005,   # shipped mucus value
    "carbon_z_lambda": 2.5e-05,
    "initial_population.placement": "z_slab",
    "initial_population.z_min": 0.0,
    "initial_population.z_max": 2e-05,
}

CONDITIONS = {"vitro": VITRO, "vivo": VIVO}

PRODUCER = {"type": 1, "count": 60, "mu_max": 5.5e-4,
            "plasmids": ["ColE2"], "conjugative": False}
SENSITIVE = {"type": 2, "count": 60, "mu_max": 5.5e-4,
             "plasmids": [], "conjugative": False}

KILL_NULL_OVERRIDES = {
    "kill_rate_colicin": 0.0,
    "kill_rate_microcin": 0.0,
    "sos_cross_induction_rate": 0.0,
}

# Non-BtuB receptor comparators must NEVER use ColE1 (which also binds BtuB).
NON_BTUB_PLASMIDS = ("ColIa", "ColB")
assert "ColE1" not in NON_BTUB_PLASMIDS

TREATMENTS: dict[str, dict] = {
    "producer": {
        "producer": dict(PRODUCER),
        "sensitive": dict(SENSITIVE),
        "overrides": {},
        "note": "ColE2 (BtuB nuclease) producer vs plasmid-free sensitive",
    },
    "null": {
        "producer": dict(PRODUCER, plasmids=[]),
        "sensitive": dict(SENSITIVE),
        "overrides": {},
        "note": "plasmid-free reference null (no plasmid carriage or lysis)",
    },
    "kill_null": {
        "producer": dict(PRODUCER),
        "sensitive": dict(SENSITIVE),
        "overrides": dict(KILL_NULL_OVERRIDES),
        "note": (
            "carriage-matched ColE2 kill-null control (same plasmid carriage "
            "and SOS lysis burden, zero receptor killing / cross-induction)"
        ),
    },
    "immunity_null": {
        "producer": dict(PRODUCER),
        "sensitive": dict(SENSITIVE, plasmids=["ColE2"]),
        "overrides": {},
        "note": (
            "carriage- and immunity-matched ColE2 control (both strains carry "
            "ColE2 plasmid and immunity)"
        ),
    },
    "colIa_producer": {
        "producer": dict(PRODUCER, plasmids=["ColIa"], conjugative=False),
        "sensitive": dict(SENSITIVE),
        "overrides": {},
        "note": "non-BtuB comparator: ColIa (CirA receptor) producer vs sensitive",
    },
    "colIa_kill_null": {
        "producer": dict(PRODUCER, plasmids=["ColIa"], conjugative=False),
        "sensitive": dict(SENSITIVE),
        "overrides": dict(KILL_NULL_OVERRIDES),
        "note": "carriage-matched ColIa (CirA receptor) kill-null control",
    },
    "colB_producer": {
        "producer": dict(PRODUCER, plasmids=["ColB"], conjugative=False),
        "sensitive": dict(SENSITIVE),
        "overrides": {},
        "note": "non-BtuB comparator: ColB (FepA receptor) producer vs sensitive",
    },
    "colB_kill_null": {
        "producer": dict(PRODUCER, plasmids=["ColB"], conjugative=False),
        "sensitive": dict(SENSITIVE),
        "overrides": dict(KILL_NULL_OVERRIDES),
        "note": "carriage-matched ColB (FepA receptor) kill-null control",
    },
}


def build_arm_config(
    cond_name: str,
    treatment: str,
    seed: int,
    *,
    b12_growth_essential: bool | None = None,
    btub_competitor_conc: float | None = None,
    b12_free_fraction: float | None = None,
    extra_overrides: dict | None = None,
) -> dict:
    if cond_name not in CONDITIONS:
        raise ValueError(f"Unknown condition: {cond_name}")
    if treatment not in TREATMENTS:
        raise ValueError(f"Unknown treatment: {treatment}")
    spec = TREATMENTS[treatment]
    for strain_key in ("producer", "sensitive"):
        plasmids = spec[strain_key].get("plasmids", [])
        if treatment.startswith(("colIa", "colB")) and "ColE1" in plasmids:
            raise ValueError("ColE1 targets BtuB and cannot be used as a non-BtuB control")

    cfg = dict(COMMON)
    cfg.update(CONDITIONS[cond_name])
    cfg.update(spec["overrides"])
    if b12_growth_essential is not None:
        cfg["metabolism.b12_growth_essential"] = bool(b12_growth_essential)
    if btub_competitor_conc is not None:
        cfg["receptor.btub_competitor_conc"] = float(btub_competitor_conc)
    if b12_free_fraction is not None:
        cfg["receptor.b12_free_fraction"] = float(b12_free_fraction)
    if extra_overrides:
        cfg.update(extra_overrides)
    cfg["seed"] = seed
    cfg["initial_strains"] = [dict(spec["producer"]), dict(spec["sensitive"])]
    return cfg


def emit(
    name: str,
    cfg: dict,
    note: str,
    *,
    output_root: pathlib.Path | None = None,
) -> pathlib.Path:
    arms_root = (ROOT / "arms") if output_root is None else pathlib.Path(output_root)
    out = arms_root / name
    out.mkdir(parents=True, exist_ok=True)
    cfg = dict(cfg)
    cfg["_comment"] = [
        "Spec 14 Section 8 reversal validation & receptor/carriage controls.",
        "type 1 = producer (or null/control counterpart), type 2 = sensitive.",
        "Read per arm: per-type counts from /agents/<step>/type, the",
        "population ledger (cumulative_mortality_colicin, divisions), and",
        "/run_provenance termination cause.",
        f"arm = {name}. {note}",
    ]
    out_file = out / "input.json"
    out_file.write_text(json.dumps(cfg, indent=1) + "\n")
    return out_file


def main(argv: list[str] | None = None) -> None:
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate Spec 14 §8 V1/V2 competition and control arms."
    )
    parser.add_argument(
        "--output-dir",
        type=pathlib.Path,
        default=None,
        help="Target directory for generated arm folders (default: arms/).",
    )
    parser.add_argument(
        "--conditions",
        nargs="+",
        choices=tuple(CONDITIONS.keys()),
        default=list(CONDITIONS.keys()),
        help="Conditions to generate (default: vitro vivo).",
    )
    parser.add_argument(
        "--treatments",
        nargs="+",
        choices=tuple(TREATMENTS.keys()),
        default=list(TREATMENTS.keys()),
        help="Treatments to generate (default: all producer/null/control arms).",
    )
    parser.add_argument(
        "--b12-growth-essential",
        dest="b12_growth_essential",
        action="store_true",
        default=None,
        help="Explicitly require exogenous B12 for Monod growth.",
    )
    parser.add_argument(
        "--no-b12-growth-essential",
        dest="b12_growth_essential",
        action="store_false",
        help="Decouple Monod growth from exogenous B12 (wild-type MetE prototroph).",
    )
    parser.add_argument(
        "--btub-competitor-conc",
        type=float,
        default=None,
        help="Explicit free BtuB competing ligand concentration (mol/m^3).",
    )
    parser.add_argument(
        "--b12-free-fraction",
        type=float,
        default=None,
        help="Fraction of bulk B12 field available as free BtuB competitor in [0, 1].",
    )
    args = parser.parse_args(argv)

    for cond_name in args.conditions:
        for treatment in args.treatments:
            spec = TREATMENTS[treatment]
            for seed in SEEDS:
                cfg = build_arm_config(
                    cond_name,
                    treatment,
                    seed,
                    b12_growth_essential=args.b12_growth_essential,
                    btub_competitor_conc=args.btub_competitor_conc,
                    b12_free_fraction=args.b12_free_fraction,
                )
                emit(
                    f"{cond_name}_{treatment}_s{seed}",
                    cfg,
                    f"{cond_name} condition, {treatment} treatment ({spec['note']})",
                    output_root=args.output_dir,
                )


if __name__ == "__main__":
    main()
