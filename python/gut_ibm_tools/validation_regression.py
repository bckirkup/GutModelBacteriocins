"""
EARI/VADI validation regression checks for CI and release workflows.

Compares simulation HDF5 output against documented biological targets
(README, EARI.md, VADI.md) and optional golden-file baselines for
deterministic short-run regression (issue #56).
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .fish_observation import flatten_fish_metrics, validate_fish_observability
from .hdf5_reader import GutIBMData
from .validation import validate_genomic_signatures, validate_spatial_signatures

# Documented validation targets (full-length simulations).
# Short CI runs use golden-file comparison instead of these thresholds.
VALIDATION_TARGETS: dict[str, dict[str, Any]] = {
    "monochromatic_score": {
        "min": 0.7,
        "description": "HiPR-FISH monochromatic patchiness",
        "references": ["VADI §75", "README validation targets"],
    },
    "comet_tail_ratio": {
        "min": 1.5,
        "description": "Advective comet-tail asymmetry (downstream/upstream toxin)",
        "references": ["EARI comet-tail", "VADI §35", "README validation targets"],
    },
    "comet_tail_asymmetry": {
        "min": 1.0,
        "description": "Concentration-weighted downstream elongation",
        "references": ["VADI §75", "docs/MECHANISMS.md"],
    },
    "hopkins_statistic": {
        "min": 0.7,
        "description": "Significant spatial clustering",
        "references": ["VADI §75"],
    },
    "resident_retention": {
        "min": 0.70,
        "max": 0.80,
        "description": "Longitudinal resident lineage retention",
        "references": ["EARI longitudinal metagenomics", "README validation targets"],
    },
    "resident_mean_bi_loci": {
        "min": 0.0,
        "description": "Residents retain BI locus complexity",
        "references": ["EARI BI locus evolution"],
    },
}

# Metrics compared against golden baselines (stochastic metrics need fixed seed).
GOLDEN_METRICS = frozenset({
    "monochromatic_score",
    "comet_tail_ratio",
    "comet_tail_asymmetry",
    "mean_exclusion_radius",
    "hopkins_statistic",
    "nnd_mean",
    "resident_retention",
    "resident_mean_bi_loci",
    "transient_mean_bi_loci",
    "transient_mean_btuB_expression",
})

# FISH observability metrics (issue #25) — deterministic given fixed RNG seed.
FISH_GOLDEN_METRICS = frozenset({
    "plasmid_dna_fish_detection_fraction",
    "plasmid_dna_fish_mean_snr",
    "immunity_hipr_detectable",
    "immunity_hcr_detectable",
    "n_detected",
    "detection_fraction",
    "monochromatic_score_detected",
    "immunity_mrna_detection_fraction",
    "immunity_mrna_hcr_detection_fraction",
    "colicin_plasmid_detection_fraction",
    "rrna_phylogroup_detection_fraction",
})

# Biological invariants enforced on every CI FISH check (VADI §75).
FISH_TARGETS: dict[str, dict[str, Any]] = {
    "immunity_hipr_detectable": {
        "max": 0.5,
        "description": "Immunity mRNA below HiPR-FISH detection (single-digit copies)",
        "references": ["VADI §75", "issue #25"],
    },
    "immunity_hcr_detectable": {
        "min": 0.5,
        "description": "HCR-FISH amplification rescues low-copy immunity mRNA",
        "references": ["VADI §75", "issue #25"],
    },
    "plasmid_dna_fish_detection_fraction": {
        "min": 0.5,
        "description": "Multicopy plasmid DNA-FISH detects colicin carriers",
        "references": ["docs/MECHANISMS.md", "issue #25"],
    },
    "rrna_phylogroup_detection_fraction": {
        "min": 0.9,
        "description": "Multicopy rRNA operons detectable by DNA-FISH",
        "references": ["VADI §75", "issue #25"],
    },
}


@dataclass(frozen=True)
class ValidationFailure:
    metric: str
    message: str

    def __str__(self) -> str:
        return f"{self.metric}: {self.message}"


def evaluate_metrics(
    h5_path: str | Path,
    *,
    step: str | None = None,
    random_seed: int = 4092,
) -> dict[str, float]:
    """Compute spatial and genomic validation metrics from an HDF5 file."""
    np.random.seed(random_seed)
    with GutIBMData(h5_path) as data:
        if not data.steps:
            raise ValueError(f"No step groups found in {h5_path}")
        target_step = step if step is not None else data.steps[-1]
        spatial = validate_spatial_signatures(data, target_step)
        genomic = validate_genomic_signatures(data)
        if "error" in genomic:
            raise ValueError("Genomic validation requires at least two HDF5 steps")
        return {**spatial, **genomic}


def check_thresholds(metrics: dict[str, float]) -> list[ValidationFailure]:
    """Check metrics against documented EARI/VADI validation targets."""
    failures: list[ValidationFailure] = []

    for name, spec in VALIDATION_TARGETS.items():
        if name not in metrics:
            continue
        value = metrics[name]
        if "min" in spec and value < spec["min"]:
            refs = ", ".join(spec.get("references", []))
            failures.append(ValidationFailure(
                name,
                f"{value:.4g} < min {spec['min']} ({spec['description']}; {refs})",
            ))
        if "max" in spec and value > spec["max"]:
            refs = ", ".join(spec.get("references", []))
            failures.append(ValidationFailure(
                name,
                f"{value:.4g} > max {spec['max']} ({spec['description']}; {refs})",
            ))

    return failures


def evaluate_fish_metrics(
    h5_path: str | Path,
    *,
    step: str | None = None,
    fish_seed: int = 25,
) -> dict[str, float]:
    """Compute FISH observability metrics from an HDF5 file (issue #25)."""
    with GutIBMData(h5_path) as data:
        if not data.steps:
            raise ValueError(f"No step groups found in {h5_path}")
        target_step = step if step is not None else data.steps[-1]
        rng = np.random.default_rng(fish_seed)
        raw = validate_fish_observability(data, target_step, rng=rng)
        return flatten_fish_metrics(raw)


def check_fish_targets(metrics: dict[str, float]) -> list[ValidationFailure]:
    """Check FISH metrics against documented detection-limit invariants."""
    failures: list[ValidationFailure] = []

    for name, spec in FISH_TARGETS.items():
        if name not in metrics:
            continue
        value = metrics[name]
        if "min" in spec and value < spec["min"]:
            refs = ", ".join(spec.get("references", []))
            failures.append(ValidationFailure(
                name,
                f"{value:.4g} < min {spec['min']} ({spec['description']}; {refs})",
            ))
        if "max" in spec and value > spec["max"]:
            refs = ", ".join(spec.get("references", []))
            failures.append(ValidationFailure(
                name,
                f"{value:.4g} > max {spec['max']} ({spec['description']}; {refs})",
            ))

    return failures


def compare_fish_golden(
    metrics: dict[str, float],
    golden: dict[str, Any],
    *,
    rtol: float = 1e-4,
    atol: float = 1e-6,
) -> list[ValidationFailure]:
    """Compare computed FISH metrics to a golden baseline within tolerance."""
    failures: list[ValidationFailure] = []
    expected = golden.get("metrics", golden)

    for name in FISH_GOLDEN_METRICS:
        if name not in expected:
            continue
        if name not in metrics:
            failures.append(ValidationFailure(f"fish.{name}", "missing from computed metrics"))
            continue
        exp = float(expected[name])
        got = float(metrics[name])
        if not np.isclose(got, exp, rtol=rtol, atol=atol):
            failures.append(ValidationFailure(
                f"fish.{name}",
                f"got {got:.6g}, expected {exp:.6g} (rtol={rtol}, atol={atol})",
            ))

    return failures


def write_fish_golden(metrics: dict[str, float], path: str | Path, *, scenario: str) -> None:
    payload = {
        "scenario": scenario,
        "metrics": {k: float(v) for k, v in metrics.items() if k in FISH_GOLDEN_METRICS},
    }
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
        f.write("\n")


def compare_golden(
    metrics: dict[str, float],
    golden: dict[str, Any],
    *,
    rtol: float = 1e-4,
    atol: float = 1e-6,
) -> list[ValidationFailure]:
    """Compare computed metrics to a golden baseline within tolerance."""
    failures: list[ValidationFailure] = []
    expected = golden.get("metrics", golden)

    for name in GOLDEN_METRICS:
        if name not in expected:
            continue
        if name not in metrics:
            failures.append(ValidationFailure(name, "missing from computed metrics"))
            continue
        exp = float(expected[name])
        got = float(metrics[name])
        if not np.isclose(got, exp, rtol=rtol, atol=atol):
            failures.append(ValidationFailure(
                name,
                f"got {got:.6g}, expected {exp:.6g} (rtol={rtol}, atol={atol})",
            ))

    return failures


def load_golden(path: str | Path) -> dict[str, Any]:
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def write_golden(metrics: dict[str, float], path: str | Path, *, scenario: str) -> None:
    payload = {
        "scenario": scenario,
        "metrics": {k: float(v) for k, v in metrics.items() if k in GOLDEN_METRICS},
    }
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
        f.write("\n")


def run_validation(
    h5_path: str | Path,
    *,
    golden_path: str | Path | None = None,
    fish_golden_path: str | Path | None = None,
    check_targets: bool = False,
    enforce_fish_targets: bool = False,
    rtol: float = 1e-4,
    atol: float = 1e-6,
) -> tuple[dict[str, float], dict[str, float] | None, list[ValidationFailure]]:
    metrics = evaluate_metrics(h5_path)
    failures: list[ValidationFailure] = []
    fish_metrics: dict[str, float] | None = None

    if golden_path is not None:
        golden = load_golden(golden_path)
        failures.extend(compare_golden(metrics, golden, rtol=rtol, atol=atol))

    if check_targets:
        failures.extend(check_thresholds(metrics))

    if fish_golden_path is not None or enforce_fish_targets:
        fish_metrics = evaluate_fish_metrics(h5_path)
        if fish_golden_path is not None:
            fish_golden = load_golden(fish_golden_path)
            failures.extend(compare_fish_golden(
                fish_metrics, fish_golden, rtol=rtol, atol=atol,
            ))
        if enforce_fish_targets:
            failures.extend(check_fish_targets(fish_metrics))

    return metrics, fish_metrics, failures


def _print_metrics(metrics: dict[str, float], *, prefix: str = "") -> None:
    label = f"{prefix} " if prefix else ""
    for key in sorted(metrics):
        print(f"  {label}{key}: {metrics[key]:.6g}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Validate GutIBM HDF5 output against EARI/VADI targets or golden baselines.",
    )
    parser.add_argument("h5_file", type=Path, help="Path to GutIBM HDF5 output")
    parser.add_argument(
        "--golden",
        type=Path,
        help="Golden JSON baseline for regression comparison",
    )
    parser.add_argument(
        "--check-targets",
        action="store_true",
        help="Enforce full EARI/VADI validation thresholds (full-length runs)",
    )
    parser.add_argument("--rtol", type=float, default=1e-4, help="Relative tolerance for golden compare")
    parser.add_argument("--atol", type=float, default=1e-6, help="Absolute tolerance for golden compare")
    parser.add_argument(
        "--write-golden",
        type=Path,
        metavar="PATH",
        help="Write computed metrics to a golden JSON file and exit",
    )
    parser.add_argument(
        "--scenario",
        default="eari_vadi_ci",
        help="Scenario label stored in --write-golden output",
    )
    parser.add_argument(
        "--fish-golden",
        type=Path,
        help="Golden JSON baseline for FISH observability regression (issue #25)",
    )
    parser.add_argument(
        "--check-fish-targets",
        action="store_true",
        help="Enforce FISH detection-limit invariants (HiPR vs HCR vs plasmid DNA-FISH)",
    )
    parser.add_argument(
        "--write-fish-golden",
        type=Path,
        metavar="PATH",
        help="Write computed FISH metrics to a golden JSON file and exit",
    )
    parser.add_argument(
        "--fish-scenario",
        default="eari_vadi_ci_fish",
        help="Scenario label stored in --write-fish-golden output",
    )
    args = parser.parse_args(argv)

    if not args.h5_file.is_file():
        print(f"ERROR: HDF5 file not found: {args.h5_file}", file=sys.stderr)
        return 2

    try:
        metrics = evaluate_metrics(args.h5_file)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    if args.write_fish_golden is not None:
        fish_metrics = evaluate_fish_metrics(args.h5_file)
        write_fish_golden(fish_metrics, args.write_fish_golden, scenario=args.fish_scenario)
        print(f"Wrote FISH golden metrics to {args.write_fish_golden}")
        _print_metrics(
            {k: fish_metrics[k] for k in sorted(fish_metrics) if k in FISH_GOLDEN_METRICS},
            prefix="fish",
        )
        return 0

    if args.write_golden is not None:
        write_golden(metrics, args.write_golden, scenario=args.scenario)
        print(f"Wrote golden metrics to {args.write_golden}")
        _print_metrics({k: metrics[k] for k in sorted(metrics) if k in GOLDEN_METRICS})
        return 0

    failures: list[ValidationFailure] = []
    fish_metrics: dict[str, float] | None = None
    if args.golden is not None:
        golden = load_golden(args.golden)
        failures.extend(compare_golden(metrics, golden, rtol=args.rtol, atol=args.atol))
    if args.check_targets:
        failures.extend(check_thresholds(metrics))

    if args.fish_golden is not None or args.check_fish_targets:
        fish_metrics = evaluate_fish_metrics(args.h5_file)
        if args.fish_golden is not None:
            fish_golden = load_golden(args.fish_golden)
            failures.extend(compare_fish_golden(
                fish_metrics, fish_golden, rtol=args.rtol, atol=args.atol,
            ))
        if args.check_fish_targets:
            failures.extend(check_fish_targets(fish_metrics))

    print("Validation metrics:")
    _print_metrics(metrics)
    if fish_metrics is not None:
        print("FISH observability metrics:")
        _print_metrics(fish_metrics, prefix="fish")

    if failures:
        print("\nValidation failures:", file=sys.stderr)
        for failure in failures:
            print(f"  - {failure}", file=sys.stderr)
        return 1

    modes = []
    if args.golden is not None:
        modes.append(f"golden ({args.golden})")
    if args.fish_golden is not None:
        modes.append(f"fish golden ({args.fish_golden})")
    if args.check_targets:
        modes.append("EARI/VADI targets")
    if args.check_fish_targets:
        modes.append("FISH targets")
    label = " and ".join(modes) if modes else "metrics only"
    print(f"\nValidation passed ({label}).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
