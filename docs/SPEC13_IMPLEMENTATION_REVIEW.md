# Spec 13 — Implementation Review: Gaps, Seams, and Phased Plan

Companion to [`SPEC13_MULTISCALE.md`](SPEC13_MULTISCALE.md), which holds the
spec as received. This document is repository analysis, not the lead's spec: it
records what Spec 13 collides with in the existing code, what it leaves
undefined, and the order in which it can be built so that each phase produces a
falsifiable result. **No part of Spec 13 is implemented.**

## 1. The reframing that motivates the architecture

The measurement that ended the Spec 12 line of work (PR #314) is not discarded
by Spec 13; it is re-interpreted as a Layer 1 input.

| Quantity | Value | Source |
|---|---|---|
| Within-patch carbon-limited capacity at 1.00x epithelial flux | ~30 agents, `1.1e7` cells/mL | PR #314, 48×48×150 patch, full flora |
| Agent granularity at that geometry | 1 agent = `3.6e5` cells/mL | model's own `global_density_cells_per_mL` |
| Healthy mucosal *E. coli* | `1e4`–`1e5` CFU/mL | Elliott et al. 2013, via Spec 13 |

Absolute comparison of rows 1 and 3 says the model is ~100x too dense and
cannot be fixed by supply (the ladder is monotone in supply, so no rung is
lower). Spec 13's claim is that the comparison is category-wrong: row 1 is a
*colonized* patch and row 3 is a *segment mean*, related by occupancy:

```
segment_mean_density ≈ occupied_fraction × within_patch_density
```

`1e4`–`1e5` CFU/mL at segment scale from a `1e7` cells/mL colonized patch
requires **0.1–1% occupancy**. That is the quantitative content of "small
transient clonal clusters", and it is the first thing Phase 1 must be able to
produce or fail to produce. It also means occupancy — not density — is the
model's primary state variable at Layer 2, and that the healthy state is
*mostly empty patches*, which is why the per-patch guard cannot be an absolute
density (Spec 13's own guard rows are 0.14–1.4 agents per patch; see the
annotation in the spec document).

## 2. Precompute vs. live IBM — the affordability split

The lead's direction is that most micropatch simulation should become
precomputed lookup rather than live IBM. Agreed, with an explicit boundary,
because the wrong side of the split silently bakes in the answer.

**Tabulate (fast, mechanism-neutral, monotone):** carbon/O₂-limited growth rate
and patch capacity as a function of (patch type, epithelial carbon flux, O₂
flux, VBF vmax, N). These equilibrate within hours of simulated time, are
monotone in each argument, and interpolate. PR #314's ladder is already the
first 15 rows of this table.

**Do not tabulate (spatially structured):** the bacteriocin producer advantage,
and the fragmentation response to a contraction. A colicin kill radius is
50–200 µm against a 100–300 µm patch, so the outcome depends on where clones
sit relative to each other, not on composition alone — tabulating advantage
against composition assumes the result the campaign is meant to predict.
Fragmentation is spatial by construction (which clones are in the sheared
layer).

**Audit gate.** A fixed fraction of lookup calls (e.g. 1%) is re-run as a live
patch and the discrepancy reported into `/run_provenance`, with a configured
tolerance that halts on exceedance. This is the same "does the mechanism
actually bind" check that caught the earlier wiring defects, applied to the
emulator; without it a stale or out-of-range table is indistinguishable from
biology.

**Cost.** Measured on this repository's own CPU path, a low-occupancy
48×48×150 patch costs ~2 min per simulated day and cost is agent-dominated
(a bloom arm is ~15x that). A 200-patch, 90-day Layer 2 run with live patches
throughout is therefore order `10^2`–`10^3` core-hours — feasible on one
machine, and well under it once healthy (near-empty) patches are served from
the table. **No GPU spend is required for Phase 1**, which matters because the
GPU path is where non-reproducibility enters.

## 3. Seams against existing mechanisms

These are the places where Spec 13 adds a mechanism the repository already has
in another form. Each must be *resolved*, not added alongside.

**S1 — Loss double-counting (highest risk).** The repository already loses
agents through peristaltic `advection`, through washout with
`washout.trap=emergent` (actual transport decides `outflow_boundary`), and
through `crypt_exit_rate`. Spec 13's contraction event is a fourth loss channel
over the same physics. Resolution: at Layer 2, a patch's *internal* washout must
be switched off and contraction becomes the sole export term, or contraction is
defined as a modulation of the existing advection field rather than an agent
sink. The population ledger must show one debit per lost agent; the closure test
is what proves it.

**S2 — Crypt refuge represented twice.** Today a crypt is a per-agent flag with
entry/exit rates and a carrying capacity (`flags.in_crypt`, persisted through
MPI transfer and checkpoints). Spec 13 wants a crypt to be a *patch type* with
its own O₂ and shear. These are alternative models of one refuge. Recommendation:
keep the agent-level crypt for within-patch structure, and let the patch type
set only boundary conditions and disruption probability — with a test asserting
that a crypt-type patch and an in-crypt agent do not both discount the same
disruption.

**S3 — VBF heterogeneity.** Layer 3 wants regional VBF density; today VBF is one
global config. Also note that `vbf.density` cannot serve as the total-bacteria
denominator the spec's "E. coli > 5% of total bacteria" criterion needs: it is
live, but only as a rate coefficient in dynamic mucin liberation (CPU
`src/fields/vbf.cpp` and CUDA `src/gpu/chemistry_kernel.cu`), never as a flora
population, and its default (`1e5` cells/mL) is ~3 orders below the FISH
total-mucus figure the spec cites. Repurposing it would change mucin liberation
rates as a side effect. Either introduce a separate flora-population parameter
or drop that criterion; do not silently reinterpret the existing key.

**S4 — Guard semantics.** The current guard is one global density threshold with
a 7-sample/300 s window. Under Spec 13 the per-patch absolute density guard is
sub-single-agent and meaningless. Proposed replacement: (i) bloom-relative
(sustained ≥10x over initialized baseline) and (ii) spatial (sustained inner-mucus
or epithelial contact) as first-class halts at patch level; (iii) the absolute
ladder applied only to the **segment mean**, where the volume is large enough
for it to mean anything; and (iv) a startup admissibility check that refuses to
run when one agent already exceeds the configured guard — which would have
flagged every configuration in this session on day one.

**S5 — Termination scope.** `TerminationCause` is run-level and authoritative.
A single blooming patch must not terminate a 200-patch run. Layer 2 needs a
per-patch status distinct from the run-level cause, with a run-level cause
reserved for "too many patches invalid".

## 4. Undefined interfaces that must be decided before code

**G1 — What is a patch instance?** Three options, with the costs that decide it:

- *(a) N `Simulation` objects in one process.* Closest to the existing code and
  the natural home for the precompute/lookup split. Blocker: GPU dispatch state
  (`src/gpu/dispatch.cpp`) and mechanics scratch (`src/gpu/mechanics_gpu.cpp`)
  are function-local statics, i.e. one instance per process — verified by
  inspection, not by running two simulations concurrently — so concurrent GPU
  patches would share device/stream and scratch. Viable **CPU-only** without
  changes; GPU would require per-instance device state.
- *(b) Patches as disjoint sub-blocks of one domain.* Reuses the slab MPI
  machinery, but chemistry couples across patch borders unless barriers are
  added, and the chemical grid is allocated globally per rank (measured ~16 GiB
  for two ranks at the large scaling grid), so it bloats badly.
- *(c) One patch per MPI rank.* Conflicts with the existing slab decomposition
  of a single global domain.

Recommendation: **(a), CPU-only, for Phase 1.**

**G2 — Contraction timing vs `bio_dt`.** `bio_dt` is 60 s; contractions are
0.2–1.0/min. Contraction should be a per-step Bernoulli event with
`p = rate·dt`, not a separately resolved timer, and it should fire as a
`post_step` hook **after** physics/migration and **before** washout/cleanup, so
it uses settled positions and books losses through the existing cleanup path.
No timestep module reordering.

**G3 — RNG streams.** Per-patch contraction and the luminal pool need
independent, patch-indexed streams so results are invariant to patch iteration
order, rank count, and CPU/GPU. Otherwise the existing parity tests
(`openmp_parity`, `compare_gpu_parity.sh`) become unfixable.

**G4 — Checkpoint/restart.** Luminal pool contents, contraction phase, per-patch
occupancy and lineage identity must all persist. If they do not, a Spot resume
silently resets the metapopulation — the same failure shape as the `mu_max`
resume wipe. Agent serialization itself is reusable as-is
(`agent_transfer.cpp` already packs crypt state, affinities, immunity escape and
genome), and the luminal pool should use it rather than a second format.

**G5 — Statistical power.** With 9–25 patches and 0.1–1% occupancy, the
expected number of occupied patches is `<<1`: the default patch count cannot
estimate the statistic it exists to produce. Layer 2 needs either
`O(10^2–10^3)` patches or a many-seed occupancy-*probability* formulation, and
the choice must be made before any result is quoted.

**G6 — Mucus volume vs surface density.** Reported observables should include a
surface density (CFU/cm²) alongside volumetric, since the biopsy-derived numbers
depend on an assumed 100–500 µm mucus thickness — a factor of 5 that must be
stated per comparison rather than absorbed into a fit.

## 5. Phased plan

Each phase ends in a measurement that can fail. No phase calibrates a parameter
to match an observation.

**Phase 0 — this document.** Spec committed, seams named, interfaces decided.
Deliverable: agreement on G1/G2 and on the S1 loss resolution. No code.

**Phase 1 — Layer 2 with lookup patches (CPU).** Patch registry, contraction as
a `post_step` Bernoulli event, luminal pool with reattachment and distal loss,
occupancy as a first-class observable, tabulated growth/capacity from the #314
ladder, and the 1% live-patch audit gate. Off by default.
*Falsifiable result:* does any physiological contraction rate / disruption
fraction combination produce a stationary occupancy in the 0.1–1% band with
segment mean `1e4`–`1e5` CFU/mL? If nothing does, Spec 13's central claim is
wrong and we will know inside one phase.
*Gates:* population-ledger closure across patches and pool; loss booked once
(S1); reproducibility under fixed seed and patch-order permutation (G3);
checkpoint round-trip of pool and occupancy (G4).

**Phase 2 — live patches where structure matters.** Promote occupied patches
from lookup to live `Simulation` instances, with the audit gate measuring the
substitution error. *Result:* the discrepancy between tabulated and live patch
trajectories, which bounds every Phase 1 conclusion.

**Phase 3 — Layer 3 chain.** Per-segment parameter sets, unidirectional luminal
transit, distal shedding as an emergent output. *Result:* retention time and
shedding rate against independent data not used to set any parameter.

**Phase 4 — perturbation and competition.** Antibiotic PK/PD (Spec 13 §
antibiotic module) and only then the bacteriocin producer/sensitive/resistant
campaigns, run live in patches per §2.

## 6. Explicit non-goals for Phase 1

Not implementing: antibiotics, Layer 3, GPU patch execution, a rewrite of the
existing washout/advection physics, or any change to frozen HDF5 field
meanings. Not adopting: Spec 13's absolute per-patch density guard rows
(unrepresentable, see S4), or any literature value as a fitted parameter — the
Elliott and Swidsinski figures are validation targets and remain so.
