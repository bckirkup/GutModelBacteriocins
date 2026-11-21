#!/usr/bin/env python3
"""Spec 13 Layer-3 mucosal occupancy 2 x 3 factorial sensitivity generator (#467 / GAP-05).

Generates deterministic multi-seed JSON configurations crossing:
  - luminal axial profile regime: {"gradient", "uniform"}
  - luminal-to-mucosal reattachment multiplier on baseline reattach_p0 = 0.05:
      0.1x (0.005), 0.3x (0.015), 1.0x (0.050)
  - seeds: (1001, 1002, 1003) -> 18 configs total.

No arbitrary <100% binary occupancy pass/fail gate is imposed; downstream
readouts evaluate post-purge (>48 h) N >= 1, 10, 50 patch density metrics,
72 h terminal provenance endpoints, stool CFU/g, and regional flatness.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
BASELINE_REATTACH_P0 = 0.05
REATTACH_SCALES: tuple[tuple[str, float], ...] = (
    ("p0_0p1x", 0.1),
    ("p0_0p3x", 0.3),
    ("p0_1p0x", 1.0),
)
ARMS: tuple[str, ...] = ("gradient", "uniform")
DEFAULT_SEEDS: tuple[int, ...] = (1001, 1002, 1003)


def _load_spec13_readout():
    mod_path = REPO / "analysis" / "spec13_gradient" / "run_readout.py"
    spec = importlib.util.spec_from_file_location("spec13_gradient_readout", mod_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Unable to load {mod_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build_sensitivity_config(
    arm: str,
    scale_label: str,
    scale_factor: float,
    seed: int,
    *,
    out_dir_rel: str = "experiments/spec13_occupancy_sensitivity/out",
) -> dict:
    if arm not in ARMS:
        raise ValueError(f"Unknown arm: {arm}")
    if scale_factor <= 0.0:
        raise ValueError(f"scale_factor must be positive, got {scale_factor}")

    readout = _load_spec13_readout()
    reattach_p0 = round(BASELINE_REATTACH_P0 * scale_factor, 10)
    stem = f"{arm}_{scale_label}_seed{seed}"

    layer3 = dict(readout.LAYER3)
    layer3["enabled"] = True
    layer3["uniform_profile"] = arm == "uniform"
    layer3["reattach_p0"] = reattach_p0
    layer3["timeseries_file"] = f"{out_dir_rel}/{stem}_timeseries.csv"
    layer3["provenance_file"] = f"{out_dir_rel}/{stem}_prov.json"

    layer2 = dict(readout.LAYER2)
    layer2["n_patches"] = layer3["n_patches_per_region"]

    return {
        "_comment": [
            "Spec 13 Layer-3 mucosal occupancy 2x3 sensitivity arm (Issue #467 / GAP-05).",
            f"arm={arm}, reattach_p0={reattach_p0} ({scale_factor}x baseline {BASELINE_REATTACH_P0}), seed={seed}.",
        ],
        "seed": int(seed),
        "total_time": readout.HORIZON_S,
        "bio_dt": readout.DT,
        "hdf5": {"enabled": False},
        "initial_strains": [
            {
                "type": 1,
                "count": 100,
                "mu_max": 5.0e-4,
                "plasmids": ["ColE1"],
                "conjugative": False,
            }
        ],
        "layer2": layer2,
        "layer3": layer3,
    }


def generate_all(
    output_dir: Path | None = None,
    seeds: tuple[int, ...] = DEFAULT_SEEDS,
) -> list[Path]:
    target_dir = (HERE / "configs") if output_dir is None else Path(output_dir)
    target_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for arm in ARMS:
        for scale_label, scale_factor in REATTACH_SCALES:
            for seed in seeds:
                cfg = build_sensitivity_config(
                    arm, scale_label, scale_factor, seed
                )
                out_path = target_dir / f"{arm}_{scale_label}_seed{seed}.json"
                out_path.write_text(json.dumps(cfg, indent=2) + "\n")
                written.append(out_path)
    return written


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Generate 2x3 Spec 13 Layer-3 mucosal occupancy sensitivity configs."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory for generated JSON configs (default: configs/).",
    )
    parser.add_argument(
        "--seeds",
        type=int,
        nargs="+",
        default=list(DEFAULT_SEEDS),
        help="Random seeds per factorial cell (default: 1001 1002 1003).",
    )
    args = parser.parse_args(argv)
    written = generate_all(output_dir=args.output_dir, seeds=tuple(args.seeds))
    print(f"Generated {len(written)} sensitivity config(s) in {(args.output_dir or (HERE / 'configs'))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
