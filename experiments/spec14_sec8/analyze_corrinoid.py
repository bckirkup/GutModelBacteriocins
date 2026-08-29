#!/usr/bin/env python3
"""Corrinoid ladder: at what corrinoid level does the producer stop losing?

Per level, paired within seed against its own no-producer null:

  L1  selection coefficient: the slope of log10(n1/n2) in decades per hour,
      fitted over the second half of the run, paired against the arm's own
      null.  Negative means the producer strain is worse off for carrying the
      plasmid.

      Not the final-sample ratio.  The two strains divide in alternating
      synchronized waves once the culture is dense, so log10(n1/n2) at the
      last output step oscillates by +-0.2 decades with the phase of the
      wave: at b12=1e-5 the last eight samples run 1.10, 1.12, 0.68, 1.18,
      1.27, 1.30, 0.78, 1.44, and reading the final one alone reports +0.16
      decades where the wave-averaged value is +0.065.  A statistic sampled
      at one instant of a synchronized population is measuring the phase.
      The tail median is reported alongside the slope for the same reason.
  L2  kills per producer lysis.  This is the efficacy of the strategy itself
      and it is the number that should move with corrinoid: each lysis is one
      certain producer death, so the strategy cannot pay below 1.0 unless the
      killed cells' descendants matter more than the lysed cell's would have.
  L3  apparent Kd = kd_colicinE_btuB * (1 + b12 / kd_b12_btuB), reported so the
      measured efficacy can be read against the binding model that predicts it.

Ratios against extinct or near-extinct denominators are reported INSUFFICIENT,
never quoted.
"""
from __future__ import annotations

import json
import math
import statistics
from pathlib import Path

from analyze import ARMS, MIN_DENOM, analyze, fmt, log_ratio  # noqa: F401

ROOT = Path("/home/ubuntu/gutibm-campaign/spec14-sec8")
LADDER = ROOT / "arms_corrinoid"

KD_COLICIN_E_BTUB = 5.0e-7
KD_B12_BTUB = 1.0e-6


def apparent_kd(b12: float) -> float:
    return KD_COLICIN_E_BTUB * (1.0 + b12 / KD_B12_BTUB)


MIN_COUNT = 20      # below this a log ratio is founder noise, not selection
TAIL_SAMPLES = 8    # ~40 min at the 300 s output interval


def wave_robust(series: list[dict]) -> tuple[float | None, float | None]:
    """Return (tail median, slope in decades/hour) of log10(n1/n2).

    Both are taken over multiple division waves so neither depends on which
    phase of the wave the run happened to stop in.
    """
    ratios = [math.log10(s["n1"] / s["n2"]) for s in series
              if s["n1"] >= MIN_COUNT and s["n2"] >= MIN_COUNT]
    if len(ratios) < TAIL_SAMPLES:
        return None, None
    tail = ratios[-TAIL_SAMPLES:]
    half = ratios[len(ratios) // 2:]
    hours = [i * 300.0 / 3600.0 for i in range(len(half))]
    mean_t = sum(hours) / len(hours)
    mean_y = sum(half) / len(half)
    denom = sum((t - mean_t) ** 2 for t in hours)
    slope = (sum((t - mean_t) * (y - mean_y) for t, y in zip(hours, half))
             / denom) if denom > 0 else None
    return statistics.median(tail), slope


def load(arm_dir: Path) -> dict:
    """Reuse analyze.analyze() by pointing it at the ladder directory."""
    import analyze as base
    base.ARMS = LADDER
    record = base.analyze(arm_dir.name)
    cfg = json.loads((arm_dir / "input.json").read_text())
    record["b12"] = float(cfg["b12_initial_conc"])
    record["treatment"] = arm_dir.name.split("_")[1]
    record["kills_per_lysis"] = (
        record["mortality_colicin"] / record["mortality_lysis"]
        if record["mortality_lysis"] else None)
    record["tail_median"], record["slope_per_h"] = wave_robust(record["series"])
    return record


def main() -> int:
    dirs = sorted(p for p in LADDER.iterdir() if (p / "output.h5").exists())
    missing = sorted(p.name for p in LADDER.iterdir()
                     if not (p / "output.h5").exists())
    records = [load(p) for p in dirs]
    print(f"arms analyzed: {len(records)}; missing output: {missing}")
    print()

    header = ("b12 mol/m3", "apparent_kd", "paired slope decades/h",
              "paired tail median", "kills/lysis (producer)",
              "kills/div (producer)", "lysis (producer)",
              "final agents (producer)", "termination", "usable seeds")
    print("| " + " | ".join(header) + " |")
    print("|" + "---|" * len(header))

    for b12 in sorted({r["b12"] for r in records}, reverse=True):
        level = [r for r in records if r["b12"] == b12]
        by_key = {(r["treatment"], r["seed"]): r for r in level}
        slopes, tails = [], []
        for seed in sorted({r["seed"] for r in level}):
            prod, null = by_key.get(("producer", seed)), by_key.get(("null", seed))
            if not prod or not null:
                continue
            if prod["slope_per_h"] is None or null["slope_per_h"] is None:
                continue
            slopes.append(prod["slope_per_h"] - null["slope_per_h"])
            tails.append(prod["tail_median"] - null["tail_median"])
        prods = [r for r in level if r["treatment"] == "producer"]
        kpl = [r["kills_per_lysis"] for r in prods
               if r["kills_per_lysis"] is not None]
        kpd = [r["kills_per_division"] for r in prods
               if r["kills_per_division"] is not None]
        lys = [r["mortality_lysis"] for r in prods]
        finals = [r["n1"] + r["n2"] for r in prods]
        causes = sorted({r["termination_cause"] for r in prods})
        enough = len(slopes) >= 3
        print("| " + " | ".join([
            f"{b12:.0e}", f"{apparent_kd(b12):.3g}",
            f"{statistics.median(slopes):+.4f}" if enough else "INSUFFICIENT",
            f"{statistics.median(tails):+.3f}" if enough else "INSUFFICIENT",
            f"{statistics.median(kpl):.3g}" if kpl else "INSUFFICIENT",
            f"{statistics.median(kpd):.3g}" if kpd else "INSUFFICIENT",
            f"{statistics.median(lys):.0f}" if lys else "-",
            f"{statistics.median(finals):.0f}" if finals else "-",
            ",".join(causes),
            str(len(slopes)),
        ]) + " |")

    print()
    print("Termination causes:")
    for cause in sorted({r["termination_cause"] for r in records}):
        n = sum(1 for r in records if r["termination_cause"] == cause)
        print(f"  {cause}: {n} arms")

    (ROOT / "results_corrinoid.json").write_text(
        json.dumps(records, indent=1, default=str) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
