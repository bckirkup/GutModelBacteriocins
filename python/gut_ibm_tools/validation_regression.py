"""Invariant checks for EARI/VADI validation outputs."""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .fish_observation import flatten_fish_metrics, validate_fish_observability
from .hdf5_reader import GutIBMData
from .path_utils import validate_input_path
from .validation import validate_genomic_signatures, validate_spatial_signatures

# Documented validation targets for full-length simulations.
REF_VADI_75 = "VADI §75"
REF_README_VALIDATION = "README validation targets"
REF_ISSUE_25 = "issue #25"

VALIDATION_TARGETS: dict[str, dict[str, Any]] = {
    "monochromatic_score": {
        "min": 0.7,
        "description": "HiPR-FISH monochromatic patchiness",
        "references": [REF_VADI_75, REF_README_VALIDATION],
    },
    "comet_tail_ratio": {
        "min": 1.5,
        "description": "Advective comet-tail asymmetry (downstream/upstream toxin)",
        "references": ["EARI comet-tail", "VADI §35", REF_README_VALIDATION],
    },
    "comet_tail_asymmetry": {
        "min": 1.0,
        "description": "Concentration-weighted downstream elongation",
        "references": [REF_VADI_75, "docs/MECHANISMS.md"],
    },
    "hopkins_statistic": {
        "min": 0.7,
        "description": "Significant spatial clustering",
        "references": [REF_VADI_75],
    },
    "resident_retention": {
        "min": 0.70,
        "max": 0.80,
        "description": "Longitudinal resident lineage retention",
        "references": ["EARI longitudinal metagenomics", REF_README_VALIDATION],
    },
    "resident_mean_bi_loci": {
        "min": 0.0,
        "description": "Residents retain BI locus complexity",
        "references": ["EARI BI locus evolution"],
    },
}

VALIDATION_METRICS = frozenset({
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

FISH_METRICS = frozenset({
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

UNIT_INTERVAL_METRICS = frozenset({
    "monochromatic_score",
    "immunity_hipr_detectable",
    "immunity_hcr_detectable",
    "plasmid_dna_fish_detection_fraction",
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
        "references": [REF_VADI_75, REF_ISSUE_25],
    },
    "immunity_hcr_detectable": {
        "min": 0.5,
        "description": "HCR-FISH amplification rescues low-copy immunity mRNA",
        "references": [REF_VADI_75, REF_ISSUE_25],
    },
    "plasmid_dna_fish_detection_fraction": {
        "min": 0.5,
        "description": "Multicopy plasmid DNA-FISH detects colicin carriers",
        "references": ["docs/MECHANISMS.md", REF_ISSUE_25],
    },
    "rrna_phylogroup_detection_fraction": {
        "min": 0.9,
        "description": "Multicopy rRNA operons detectable by DNA-FISH",
        "references": [REF_VADI_75, REF_ISSUE_25],
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
    rng = np.random.default_rng(random_seed)
    with GutIBMData(h5_path) as data:
        if not data.steps:
            raise ValueError(f"No step groups found in {h5_path}")
        target_step = step if step is not None else data.steps[-1]
        spatial = validate_spatial_signatures(data, target_step, rng=rng)
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


def check_invariants(
    data: GutIBMData,
    metrics: dict[str, float],
    fish_metrics: dict[str, float],
) -> list[ValidationFailure]:
    """Check output liveness, finite metrics, bounds, and relationships."""
    failures: list[ValidationFailure] = []
    all_metrics = {**metrics, **fish_metrics}
    for name in VALIDATION_METRICS | FISH_METRICS:
        if name not in all_metrics:
            failures.append(ValidationFailure(name, "missing from computed metrics"))
        elif not np.isfinite(all_metrics[name]):
            failures.append(ValidationFailure(name, "must be finite"))

    for name in UNIT_INTERVAL_METRICS:
        if name in all_metrics and not 0.0 <= all_metrics[name] <= 1.0:
            failures.append(ValidationFailure(
                name, f"{all_metrics[name]:.6g} is outside [0, 1]",
            ))

    for name in ("n_detected",):
        if name in all_metrics and all_metrics[name] < 0.0:
            failures.append(ValidationFailure(name, "must be non-negative"))

    summary_steps = data.steps_for("summary")
    if not summary_steps:
        failures.append(ValidationFailure("summary", "at least one summary interval is required"))
        return failures

    last_step = data.steps[-1]
    last_summary = data.get_summary(last_step)
    if last_summary["num_agents"] <= 0:
        failures.append(ValidationFailure("num_agents", "final population must be non-zero"))
    agents = data.get_agents(last_step)
    n_agents = len(agents.get("id", agents.get("x", [])))
    if "n_detected" in fish_metrics and fish_metrics["n_detected"] > n_agents:
        failures.append(ValidationFailure(
            "n_detected",
            f"{fish_metrics['n_detected']:.6g} exceeds analyzed agents {n_agents}",
        ))

    diagonal = data.domain_diagonal
    for name in ("nnd_mean", "mean_exclusion_radius"):
        if name in metrics and not 0.0 < metrics[name] < diagonal:
            failures.append(ValidationFailure(
                name,
                f"{metrics[name]:.6g} must be strictly between zero and "
                f"domain diagonal {diagonal:.6g}",
            ))

    retention = metrics.get("resident_retention")
    if retention is not None and not 0.0 < retention < 1.0:
        failures.append(ValidationFailure(
            "resident_retention", "must be strictly inside (0, 1)",
        ))

    if (
        "immunity_hcr_detectable" in fish_metrics
        and "immunity_hipr_detectable" in fish_metrics
        and fish_metrics["immunity_hcr_detectable"]
        < fish_metrics["immunity_hipr_detectable"]
    ):
        failures.append(ValidationFailure(
            "immunity_detection",
            "HCR detection cannot be lower than HiPR detection",
        ))
    if (
        "rrna_phylogroup_detection_fraction" in fish_metrics
        and "immunity_mrna_detection_fraction" in fish_metrics
        and fish_metrics["rrna_phylogroup_detection_fraction"]
        < fish_metrics["immunity_mrna_detection_fraction"]
    ):
        failures.append(ValidationFailure(
            "rRNA_detection",
            "multicopy rRNA detection cannot be lower than single-copy mRNA detection",
        ))

    return failures


def run_validation(
    h5_path: str | Path,
    *,
    check_targets: bool = False,
    enforce_fish_targets: bool = False,
) -> tuple[dict[str, float], dict[str, float] | None, list[ValidationFailure]]:
    safe_h5_path = validate_input_path(h5_path)
    with GutIBMData(safe_h5_path) as data:
        metrics = evaluate_metrics(safe_h5_path)
        fish_metrics = evaluate_fish_metrics(safe_h5_path)
        failures = check_invariants(data, metrics, fish_metrics)
    if check_targets:
        failures.extend(check_thresholds(metrics))
    if enforce_fish_targets:
        failures.extend(check_fish_targets(fish_metrics))

    return metrics, fish_metrics, failures


def _print_metrics(metrics: dict[str, float], *, prefix: str = "") -> None:
    label = f"{prefix} " if prefix else ""
    for key in sorted(metrics):
        print(f"  {label}{key}: {metrics[key]:.6g}")


def _build_validation_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate GutIBM HDF5 output against EARI/VADI and FISH invariants.",
    )
    parser.add_argument("h5_file", type=Path, help="Path to GutIBM HDF5 output")
    parser.add_argument(
        "--check-targets",
        action="store_true",
        help="Enforce full EARI/VADI validation thresholds (full-length runs)",
    )
    parser.add_argument(
        "--check-fish-targets",
        action="store_true",
        help="Enforce FISH detection-limit invariants (HiPR vs HCR vs plasmid DNA-FISH)",
    )
    return parser


def _load_metrics_or_exit(h5_file: Path) -> dict[str, float] | None:
    if not h5_file.is_file():
        print(f"ERROR: HDF5 file not found: {h5_file}", file=sys.stderr)
        return None
    try:
        return evaluate_metrics(h5_file)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return None


def _collect_validation_failures(
    args: argparse.Namespace,
    metrics: dict[str, float],
) -> tuple[list[ValidationFailure], dict[str, float] | None]:
    failures: list[ValidationFailure] = []
    fish_metrics = evaluate_fish_metrics(args.h5_file)

    with GutIBMData(args.h5_file) as data:
        failures.extend(check_invariants(data, metrics, fish_metrics))
    if args.check_targets:
        failures.extend(check_thresholds(metrics))
    if args.check_fish_targets:
        failures.extend(check_fish_targets(fish_metrics))

    return failures, fish_metrics


def _validation_mode_label(args: argparse.Namespace) -> str:
    modes = []
    if args.check_targets:
        modes.append("EARI/VADI targets")
    if args.check_fish_targets:
        modes.append("FISH targets")
    return "invariants" + (f" and {' and '.join(modes)}" if modes else "")


def _handle_validate(args: argparse.Namespace, metrics: dict[str, float]) -> int:
    failures, fish_metrics = _collect_validation_failures(args, metrics)

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

    print(f"\nValidation passed ({_validation_mode_label(args)}).")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = _build_validation_parser().parse_args(argv)

    metrics = _load_metrics_or_exit(args.h5_file)
    if metrics is None:
        return 2

    return _handle_validate(args, metrics)


if __name__ == "__main__":
    raise SystemExit(main())
