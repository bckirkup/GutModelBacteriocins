# Review: Spec 12 Addendum A (Damköhler scale-gating) and Spec 14 §8

Status: review of an externally supplied package, recomputed against the
shipped parameter set. No code changed. This document decides what is
implementable as written, what must be re-derived first, and in what order.

Package under review:

| Supplied document | Claim |
|---|---|
| Spec 12 Addendum A | `agent_carbon_coupling` is not broken; it is scale-gated below the `Da = 1` aggregate radius |
| `spec12_da_report.md` + `spec12_da_sweep.csv` | Numerical spherical Michaelis–Menten reaction–diffusion backing for the above |
| Spec 14 (prophage induction) | Accumulated DNA damage, independent lysogen state, Hill induction, deferred free-phage field |
| Devin brief | Ordering: Spec 14 §8 validation first, then A1, A2, then Spec 14 Changes 1–2 |

## 1. Verdict

The reframing is accepted in substance and rejected in its numbers.

"`agent_carbon_coupling` is a no-op below the `Da = 1` aggregate radius, and
the VBF continuum sink is the correct sub-crossover representation" is a
better-supported claim than "the coupling does not work", and it is consistent
with what Spec 12 measured. It survives recomputation on this model's own
parameters — but the crossover radius moves by about 2x, the proposed
regression thresholds do not hold at either parameter set, and the diagnostic
as designed cannot measure the radii it is being asked to validate.

The supplied `Da` values are literature-anchored, not GutIBM-anchored. Every
number in Addendum A §3–§4 uses `q_cell = 8.33e-19 mol/cell/s`,
`C_inf(O2) = 20 µM` and `C_inf(carbon) = 100 µM`. None of the three is what
this repository ships.

## 2. Recomputation on shipped parameters

Same formula as the Addendum, `R_50 = sqrt(D * C_inf / (rho * q_cell))`, with
the values actually compiled in:

| Quantity | Shipped value | Source | Addendum value |
|---|---|---|---|
| `D` carbon | `5.0e-10 m²/s`, retardation 1.0 | `src/io/input_parser.cpp` default chemical table | `6.0e-10` |
| `C_inf` carbon | `5.0e-3 mol/m³` (5 µM) | same table; `km_carbon` is the same 5e-3 | `0.1 mol/m³` (100 µM) |
| `q_cell` carbon | `1.44e-19 mol/cell/s` | `mu_max 5.0e-4 /s` x `biomass 5.76e-16 kg` x `yield_carbon 0.5 mol/kg` | `8.33e-19` |
| `D` O2 | `2.1e-9 m²/s` | `OxygenConfig::D_free` | `2.0e-9` |
| `C_inf` O2 | `5.5e-5 mol/m³` (0.055 µM) | `OxygenConfig::epithelial_conc` | `0.02 mol/m³` (20 µM) |
| `q_cell` O2 | `6.0e-18 mol/cell/s` | `q_consumption 1e-14 x mu_max` + `q_maintenance 1e-18` | `8.33e-19` |

Resulting half-depletion radii:

| Local density | Carbon `R_50`, shipped | Carbon `R_50`, Addendum | O2 `R_50`, shipped | O2 `R_50`, Addendum |
|---|---|---|---|---|
| `1e10` cells/mL | 41.7 µm | 87.1 µm | 1.39 µm | 70.8 µm |
| `1e11` cells/mL | 13.2 µm | 27.6 µm | 0.44 µm | 22.4 µm |
| `1e12` cells/mL | 4.2 µm | 8.7 µm | 0.14 µm | 7.1 µm |

Two conclusions, and they point in opposite directions.

**Carbon: the Addendum is right, and the crossover is nearer than it says.**
Shipped carbon `R_50` is about 2x smaller than the Addendum's at every density.
The brief's own geometric argument — one agent in a 2 µm voxel is
`1.25e11 cells/mL`, which is arithmetically correct — puts the shipped carbon
crossover at **11.8 µm**, not 19.6 µm. That is ~6 voxels of radius, about 900
agents at one per voxel, not ~10 voxels and ~3,900 agents. The scale-gating
story gets *stronger*: the crossover sits inside the Spec 13 patch dimension
and inside the 10 µm delivery support that has been the default since #336,
rather than just above them.

**Oxygen: the Addendum's O2 column cannot be applied to any shipped config.**
At `epithelial_conc = 5.5e-5 mol/m³` the O2 crossover is *sub-voxel*: a single
agent alone in one 2 µm voxel is already at `Da ≈ 5`, the opposite of the
Addendum's `Da(single agent) = 0.003`. This is not a new finding — it is the
known 10³ unit offset recorded in `docs/PARAMETERS.md` §Oxygen and
`docs/SPEC13_IMPLEMENTATION_REVIEW.md` S11: 42 mmHg of dissolved O2 is 55 µM,
i.e. `5.5e-2 mol/m³`, and the shipped default is `5.5e-5`, with `Km = 1e-6`
carrying the same offset so the Monod *ratio* is preserved while the absolute
inventory is not. Recomputed at the sourced `5.5e-2 mol/m³`, O2 `R_50` is
13.9 µm at `1e11` cells/mL and the Addendum's band is recovered.

So an O2 Damköhler diagnostic run against shipped defaults would report
`Da >> 1` in every occupied voxel and would be measuring the unit offset, not
biology. Either the diagnostic ships carbon-only, or it ships with that caveat
attached to its output.

## 3. Defects in Change A1 as specified

**(a) The regression thresholds contradict the design's own `R_50`.** The test
asks for `da_max < 0.05` on a 4-voxel-radius blob at one agent per voxel. Since
`Da = (R_coh / R_50)²` and that blob's `R_coh` is its own radius (8 µm),
`da_max` is 0.17 at the Addendum's `R_50 = 19.6 µm` and **0.46** at the shipped
carbon `R_50 = 11.8 µm`. A correct implementation fails this assertion at
either parameter set. The thresholds must be derived from the `R_50` implied by
the config under test, not carried over from the literature table.

**(b) `R_coh` saturates, so `da_max` cannot resolve large aggregates.** `R_coh`
is the equivalent-sphere radius of the *occupied* volume within a cubic
neighbourhood of half-width `W`. At the default `W = 5`, a fully occupied
neighbourhood is `11³` voxels, so `R_coh` is capped at
`cbrt(3 * 1331 * dx³ / 4π) = 13.6 µm` and `Da_loc` saturates at
`(13.6/R_50)²` — 0.48 with the Addendum's `R_50`, i.e. the 15-voxel blob
*fails* its own `da_max > 0.5` assertion, and no blob larger than ~14 µm can be
distinguished from any other. The filter window is the measurement ceiling.
Either `W` scales with the aggregate being probed, or the diagnostic reports
saturation explicitly, or `R_coh` comes from a connected-component extent
rather than a fixed box.

**(c) The A2 crossing prediction is anchored to the wrong `R_50`.** A2's
acceptance criterion is "crossing near 22 µm at 1e11, near 7 µm at 1e12, report
if off by more than ~2x". Shipped-parameter carbon predicts 13.2 µm and 4.2 µm
— already a factor 2.1 from the stated targets before the solver is consulted.
Run against the supplied targets, A2 would report a spurious solver
discrepancy. Its prediction must be recomputed from the config it runs.

**(d) Two smaller items.** The realized-`q` reference cites
`qssa_solver.cpp:496-497`; the expression is at ~937–940 on current main. And
A2 requires growth, division and motility off, but not mechanics: at
`1e12 cells/mL` a 2 µm voxel holds 8 agents of 0.5 µm radius (~52% volume
fraction), so contact forces will disperse the clump during the run unless
mechanics is also disabled or the clump is held by an explicit constraint.

None of these are reasons to reject A1/A2. They are reasons not to implement
them from the supplied constants.

## 4. Spec 14 §8: run it, but it is not yet a clean test

The reversal target — an E-type nuclease colicin producer excludes a λ-lysogen
in ~6 h in vitro, both persist to 10 d in dixenic mice — is the right shape of
validation for this model, and §8 needs no code. Two measured facts from this
session's earlier work mean it cannot be interpreted at face value:

- **The producer treatment is nearly inert at shipped potency.** The RPS probe
  measured 11 receptor-mediated kills against 7,629 divisions with
  `retardation = 50` hard-coded for ColE1/E2. A "coexistence" result in the
  *in vivo* arm is therefore predicted by colicin weakness alone, with no
  spatial explanation required. The in-vitro arm is the discriminating one: if
  the producer cannot exclude the sensitive strain in a well-mixed voxel
  either, §8 has measured colicin potency and said nothing about spatial
  structure. That is the third branch of the spec's own interpretation guide,
  and it is the outcome I expect.
- **A single patch cannot report which strain wins.** Three identical-treatment
  seeds ended R-only, S-only and C+S; founders collapse to ~30 live in ~2 h and
  the survivors sweep. Any §8 in-vivo arm therefore needs >=5 seeds with paired
  within-seed contrasts, and its "coexistence" must be separated from
  stochastic non-resolution before it is called coexistence.

Run order below reflects this: §8 first, but scored as a colicin-potency
measurement with a spatial contrast, not as a validation pass/fail.

## 5. Recommended order

1. **Spec 14 §8, no code.** Well-mixed vs full-mucus, producer vs sensitive
   non-producer, >=5 seeds per arm, paired. Report abundance ratio vs time and
   time to 100-fold exclusion, plus realized receptor-mediated kills per
   division so a null is attributable. Check `toxin_cutoff` (200 µm default)
   and the QSSA support in the well-mixed arm before concluding anything.
2. **Re-derive the Da constants for this repository** (done above) and fold the
   corrected table into Addendum A before implementing it.
3. **A1**, carbon-first, with thresholds computed from the config under test
   and the `R_coh` saturation resolved. Default off, byte-identical when
   disabled.
4. **A2**, with predictions recomputed per config and mechanics disabled.
5. **Spec 14 Changes 1–2** only if §8 shows spatial structure alone does not
   account for the reversal.
6. **Change 3 (free phage) stays deferred.** A ~12-virion burst as a continuum
   field is a questionable representation at that copy number, independent of
   whether the ecology needs it.

## 6. Unchanged prohibitions

The Addendum's own §5.3 list stands and is consistent with this repository's
existing decisions: do not raise `agent_carbon_coupling` to make an effect
visible, do not raise per-cell uptake beyond literature support, do not remove
the VBF continuum sink — it is the correct sub-crossover representation, which
is precisely what the scale-gating argument establishes.

Two preprint caveats carry into anything built on this package: both source
papers (bioRxiv 10.64898/2026.05.02.721930 and 10.64898/2026.05.26.727859) are
unreviewed, and the Spec 14 Hill parameters (`hill_n`, `damage_half_max`) are
calibration assumptions, not measurements, as Spec 14 itself states.
