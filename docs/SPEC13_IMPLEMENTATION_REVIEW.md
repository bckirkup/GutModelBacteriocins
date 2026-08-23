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

**S1 — Two loss channels with different reseeding kernels (revised).** An
earlier draft of this review called Spec 13's contraction a redundant fourth
loss channel and proposed merging it with washout. That was wrong, and the
project lead's objection is the correct physics: washout and matrix failure are
different events with different fates for the cell.

- *Washout* is single-cell detachment into flow. The cell tumbles in the lumen
  and its probability of founding a new patch is low.
- *Contraction* is failure of the mucin gel. A gel fragment departs with its
  clonal cluster intact, so what arrives downstream is a multi-cell seed with a
  substantially higher establishment probability, and with a composition (which
  clones travelled together) that matters for bacteriocin interference.

Both therefore coexist at Layer 2, and the modelling requirement is not one
merged channel but **two channels feeding one luminal pool with distinct
reseeding kernels** — `p_establish(single)` << `p_establish(fragment)`, with
fragment size and composition carried along.

What is actually in the code today, which is narrower than that earlier draft
implied (`src/core/simulation.cpp`, washout stage):

- In the default `emergent` mode there is exactly **one** departure event: an
  agent that transport has carried to `z >= z_max` is set `DEAD`, booked as
  `outflow_boundary`, and recorded through `lineage_.record_washout`. Advection
  is what moves it there; it is not a second, separate sink. `outflow_washout`
  only exists in the non-default `imposed` mode, where the `mu < gamma`
  comparison removes cells before transport reaches the lumen — the two modes
  are alternatives, never both.
- Agents with `flags.in_crypt` are skipped entirely, which is the existing
  refuge behaviour (see S2).
- Departure is **terminal**: the agent is deleted at the boundary. There is no
  export pool and no reattachment, i.e. the model currently hard-codes
  `p_reattach = 0` for the one channel it has.

So the real gap is the opposite of double-counting: the model has a single
terminal single-cell loss and no representation of cluster-preserving
fragmentation at all. Layer 2 must add (i) a luminal pool that receives
departures instead of deleting them, (ii) contraction as a second, fragment-wise
departure, and (iii) two reattachment probabilities. Only one thing must not be
double-counted — an agent already carried past `z_max` by advection in a step
must not also be taken by that step's contraction event; the population ledger
must show one debit per departing agent, and the closure test is what proves it.

Usefully, `outflow_boundary` is already the flux the Layer 3 shedding
observable needs — today it is a death counter, and it becomes a stool-export
rate once the pool exists to receive it.

**S2 — Crypt refuge represented twice.** Today a crypt is a per-agent flag with
entry/exit rates and a carrying capacity (`flags.in_crypt`, persisted through
MPI transfer and checkpoints). Spec 13 wants a crypt to be a *patch type* with
its own O₂ and shear. These are alternative models of one refuge.

Resolved by S6: patch type is a property of the **host location**, so the patch
type carries the location's boundary conditions and disruption probability,
while the agent-level `in_crypt` flag remains within-patch structure. A test
must assert that a crypt-*type* patch and an `in_crypt` agent do not both
discount the same disruption event.

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

**S6 — Patch type is host anatomy, not a stochastic label.** A location's
character is fixed by the host: a crypt is a crypt at every contraction, with
its own depth, O₂ at base, shear protection and mucus thickness, and proximal
vs distal position is likewise fixed. Consequences:

- Patch type is **persistent per-patch identity**, assigned once at
  initialization and checkpointed (G4), never re-drawn per contraction event.
- The *fractions* of each type, and the per-type boundary conditions, are
  anatomical inputs that can be sourced independently (crypt density per unit
  area, mucus thickness, epithelial O₂). They are therefore **not free
  parameters available for fitting** to an occupancy target.
- Occupancy is conditional on type, which changes the observable: the segment
  mean is a type-weighted sum `Σ_t f_t · occ_t · density_t`, not a single
  occupancy number. The interesting and falsifiable version of Spec 13's claim
  is that persistent occupancy concentrates in crypt patches acting as refuges
  while exposed patches are mostly empty and transiently reseeded — which is a
  much sharper prediction than an aggregate 0.1–1%.

**S7 — The anatomy is axially non-uniform, so Layer 3 is a gradient, not a
chain of interchangeable segments.** Crypt density, mucus thickness, epithelial
O₂, VBF density, pH and contraction frequency all vary proximal→distal.
Therefore the S6 quantities are **functions of axial position**: `f_t(x)` and
the per-type boundary conditions both vary, and the type-weighted mean becomes
`Σ_t f_t(x) · occ_t(x) · density_t(x)`. Consequences:

- Comparisons must state their axial position. A biopsy-derived density is a
  measurement *somewhere*, and the model's answer at the wrong `x` is not a
  disagreement.
- Shedding is dominated by conditions at the distal end rather than by an
  average patch, so retention time and stool flux are sensitive to the distal
  end of the profile specifically — the least well-characterised end.
- The gradient must be specified as a **profile with sourced endpoints and a
  stated interpolation** (e.g. monotone in mucus thickness and crypt density
  between proximal and distal anchors). A per-segment table of independent
  values would make every segment a free parameter, which is the fitting route
  this project has explicitly refused.
- A same-parameters sanity arm (gradient collapsed to uniform) belongs in the
  Layer 3 gates, so the gradient's contribution is separable from the rest of
  the mechanism.

**S8 — The luminal growth compartment is a population, not a pool, and cannot
be agents.** The revised Layer 3 has detached cells *replicating* in transit,
which is what supplies the ~10x amplification between mucosal detachment and
stool. Consequences:

- Magnitude forbids agents: `1e8` CFU/g over ~100 g of content is `1e10` cells
  per region. The compartment must be **per-lineage counts** with an anaerobic
  growth term, so the S1 luminal pool splits into a transport buffer (agents,
  briefly) and a residence compartment (counts).
- Composition still has to survive, because the bacteriocin question needs to
  know which clones travelled and grew together. Counts must therefore be
  indexed by lineage/BI-locus genotype, not aggregated, and reattachment must
  instantiate an agent from a count with its genotype intact. The
  `agent_transfer.cpp` serialization covers that instantiation; it does not
  cover residence.
- Growth in the compartment competes for luminal carbon and produces acetate,
  so the compartment needs its own scalar chemistry with the same yield and
  acid-inhibition parameters as Layer 1 — implemented once and shared, or the
  two layers will disagree about physiology.
- The population ledger must extend across the mucosa/lumen boundary: a
  departure is a transfer, not a death, and luminal division and washout to
  stool are new ledger events. The existing closure machinery is the right home
  for this and must be extended rather than bypassed.

**S9 — The wall-to-stool equation is degenerate against its only target.**
`C_stool` is set by `k_mucus` (0.3–1.0/h), `k_0` and `α_HAPC`, `f_edge`
(explicitly 0–1, i.e. unconstrained), `r_luminal`, `T_transit` and `p_0` —
seven quantities, several with an order of range, against one validation number
(`1e6`–`1e8` CFU/g, itself two orders wide). Many combinations reproduce it, so
agreement with stool CFU/g is **not evidence** for any of them. Before Phase 3
is run, additional independent observables must be nominated, e.g.:

- *Burstiness.* HAPCs are ~6/day, 95% antegrade, 80% daytime. Contraction-driven
  detachment therefore predicts a pulsed shedding time series with a diurnal
  envelope, while mucus turnover alone predicts a smooth one. The temporal
  structure of shedding discriminates the mechanisms; the mean does not.
- *Growth state in transit.* A luminal amplification of ~10x over 25–40 h
  implies a specific actively-dividing fraction on arrival, which is measurable
  independently of CFU/g.
- *Recovery kinematics.* The time to re-establish stool counts after clearance
  is set by reattachment and mucosal regrowth, not by the steady-state balance,
  so it constrains a different combination of the same parameters.

Report Phase 3 against a *set* of observables and state which parameter
combination each one pins. A single-number match must not be presented as
validation.

**S10 — The axial abundance gradient is not supported, which makes flatness a
prediction rather than an input.** The revision records that the
Enterobacteriaceae fraction does **not** vary significantly proximal→distal
(P=0.09) and instructs uniform initialization at ~1–2% (0.1–5% individual
range). This narrows S7 rather than contradicting it: the *physiology* gradients
(mucus, pH, transit, motility, O₂ prior) remain, but E. coli *abundance* must
not be tuned axially. The corresponding test is sharp and worth stating as a
Phase 3 gate: with those gradients imposed, the model must produce a
**mucosally flat** Enterobacteriaceae fraction while stool counts are set by the
distal luminal compartment. A model that produces a strong mucosal axial
gradient is wrong even if its stool number is right.

**S11 — The revision's apical O₂ prior is 4–20x below the repository default,
and the O₂ sink is the one uptake path that never got the delivery fix.** The
revision distinguishes tissue-side electrode values (30 mmHg serosal, 39 mmHg
sigmoid) from the apical mucus surface, for which it states there are **no
segment data** and offers a broad prior of 2–10 mmHg. The repository default is
`oxygen.epithelial_conc = 55.0e-6` annotated `~42 mmHg` — i.e. the tissue-side
number used as the apical boundary condition, which is the boundary that drives
the entire Spec 12 aerobic/overflow phase plane.

Two further things need settling before that boundary is trusted, both grounded
in the code rather than in the spec:

- *The value and its own annotation disagree by 10³.* The field's units are SI
  (`mol/m³`; `grid_dx = 2e-6 m`, `cell_volume = 8e-18 m³`), and 42 mmHg of
  dissolved O₂ is ~55 µM = `5.5e-2 mol/m³`, not `5.5e-5`. `Km = 1.0e-6 mol/m³`
  (1 nM, against a literature O₂ affinity of ~0.1–1 µM) carries the same 10³
  offset, so the Monod *ratio* is preserved while the absolute O₂ inventory is
  not — and respiration is booked in absolute moles per cell.
- *Agent respiration is an explicit zero-order sink, clipped before diffusion.*
  `apply_agent_reactions` does `reac(oxygen, cell) -= o2_use / cell_vol` with
  `o2_use = q_consumption·µ·(1−ferm) + q_maintenance`; there is no `uptake_limit`
  path for oxygen. That is exactly the defect class fixed for carbon in #308 and
  #310, unfixed here. At `q_maintenance = 1e-18 mol/s/cell` over an 8 fL voxel
  this demands `~7.5 mol/m³` per 60 s step against a boundary value of
  `5.5e-5`, so an occupied voxel can only be funded by within-step resupply —
  which is precisely what an explicit pre-diffusion sink cannot do.

**This has now been measured, and the sink is refused in full.** Per-species
reaction clips are already accounted (`add_reaction_clip`), and all fifteen
supply-ladder arms (PR #314) ran with `oxygen.enabled = true` and the metabolic
switch on, so the measurement needed no new runs — see
`experiments/density_limitation/oxygen_audit/`. At 24 h, cumulative oxygen
reaction clip is `5.7e-12` mol at N=28 rising to `8.2e-10` mol at N=6317, which
is **1.2–1.4x an independent estimate of total respiratory demand** (integrating
`q_maintenance·N + q_consumption·µ·(1−ferm)·N` over the trajectory) and
`5e3`–`5e5` times the **total oxygen ever delivered across the epithelial
boundary** (`1.1e-15`–`1.6e-15` mol). Carbon clip is exactly `0` in the same
arms, so this is oxygen-specific and not a general closure failure.

Consequences, all confirmed in the same data:

- *The oxygen field is static.* Mean O₂ is `4.586e-6 mol/m³` at N=28 and
  `4.584e-6` at N=6317 — a 0.04% difference across a 200x population range, and
  unchanged from step 0. Agents cannot deplete oxygen at all, so there is no
  oxygen gradient, no anoxic core, and no density feedback through O₂.
- *The Spec 12 phase plane is one-dimensional in practice.* With a static field,
  `f_O2 = C/(Km+C) = 0.821` everywhere and always, so the only live term in
  `resp_capacity = f_O2 · mu_crit / max(µ, mu_crit)` is µ. The fermentation
  fractions reported in #313 and #314 (0.33 → 0.79 with density) are therefore
  **growth-rate driven, not O₂ driven**; the earlier description of them as
  "purely O₂-driven" is withdrawn.
- *Correcting the 10³ unit offset does not fix it.* At `5.5e-2 mol/m³` a voxel
  holds `4.4e-19` mol while one agent demands `6e-17` mol per 60 s step — still
  ~100x the voxel inventory. Oxygen needs the same delivery-limited implicit
  treatment as carbon (#308/#310); no boundary value rescues an explicit
  pre-diffusion sink.

So S11 is not a spec question but an implementation prerequisite: **no apical O₂
value should be adopted, and no O₂-dependent Layer 2/3 result should be
believed, until oxygen has an `uptake_limit`-style delivery path and the units
are reconciled.** The oxygen-dependent parts of Spec 12 (`mu_crit`, the
overflow axis, `anaerobic_mu_factor`) remain untested by any run to date.

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

**G5 — Statistical power, stratified by type.** With 9–25 patches and 0.1–1%
occupancy, the expected number of occupied patches is `<<1`: the spec's default
patch count cannot estimate the statistic it exists to produce. Given S6 this
resolves to a **type-stratified** formulation: occupancy is estimated per patch
type, with enough patches (or enough seeds) per type for the rarest type that
carries signal, and the type fractions taken from anatomy rather than chosen for
convenience. Aggregate occupancy is then a derived quantity, and no aggregate
number should be quoted without its per-type breakdown.

**G6 — Mucus volume vs surface density.** Reported observables should include a
surface density (CFU/cm²) alongside volumetric, since the biopsy-derived numbers
depend on an assumed 100–500 µm mucus thickness — a factor of 5 that must be
stated per comparison rather than absorbed into a fit.

## 5. Phased plan

Each phase ends in a measurement that can fail. No phase calibrates a parameter
to match an observation.

**Phase 0 — this document.** Spec committed, seams named, interfaces decided.
Deliverable: agreement on G1/G2, on the S1 two-channel loss model, and on the
Layer 2 formulation question (patch count vs occupancy probability, G5). No
code.

**Phase 1 — Layer 2 with lookup patches (CPU).** Patch registry, contraction as
a `post_step` Bernoulli event, a luminal pool fed by both departure channels of
S1 with separate single-cell and fragment reattachment probabilities and distal
loss, occupancy as a first-class observable, tabulated growth/capacity from
the #314 ladder, and the 1% live-patch audit gate. Off by default.
*Falsifiable result:* does any physiological contraction rate / disruption
fraction combination produce a stationary occupancy in the 0.1–1% band with
segment mean `1e4`–`1e5` CFU/mL? If nothing does, Spec 13's central claim is
wrong and we will know inside one phase. The discriminating parameter is the
fragment-vs-single establishment ratio of S1, since it is what decides whether
reseeding can balance loss at low occupancy; report the result against it.
*Gates:* population-ledger closure across patches and pool; each departing
agent debited exactly once and attributed to one channel (S1); reproducibility
under fixed seed and patch-order permutation (G3); checkpoint round-trip of
pool and occupancy (G4).

**Phase 2 — live patches where structure matters.** Promote occupied patches
from lookup to live `Simulation` instances, with the audit gate measuring the
substitution error. *Result:* the discrepancy between tabulated and live patch
trajectories, which bounds every Phase 1 conclusion.

**Phase 3 — Layer 3 gradient.** Axial profiles for the S6/S7 anatomical
quantities (sourced endpoints plus a stated interpolation, not per-segment free
values), unidirectional luminal transit, distal shedding as an emergent output.
*Result:* retention time and shedding rate against independent data not used to
set any parameter, reported at a stated axial position, with a
gradient-collapsed-to-uniform arm to separate the gradient's contribution.

**Phase 4 — perturbation and competition.** Antibiotic PK/PD (Spec 13 §
antibiotic module) and only then the bacteriocin producer/sensitive/resistant
campaigns, run live in patches per §2.

## 6. Explicit non-goals for Phase 1

Not implementing: antibiotics, Layer 3, GPU patch execution, a rewrite of the
existing washout/advection physics, or any change to frozen HDF5 field
meanings. Not adopting: Spec 13's absolute per-patch density guard rows
(unrepresentable, see S4), or any literature value as a fitted parameter — the
Elliott and Swidsinski figures are validation targets and remain so.
