# Coherence Diagnosis: density, spatial structure, and bacteriocin efficacy

Why genomic diversity, spatial structure, and bacteriocin efficacy currently fail
to cohere in GutIBM at realized cell densities — and the cheap experimental
program that resolves it without week-long runs.

Companion to [`docs/MULTI_SCALE_EXPERIMENTATION.md`](MULTI_SCALE_EXPERIMENTATION.md),
which covers *how* to compose runs. This document covers *what is wrong* and
*what to measure*.

---

## 0. Summary

The incoherence is not a deep conceptual impasse about sparse populations. It is
three concrete, verified defects, each of which independently destroys one leg of
the tripod:

| # | Defect | Consequence |
|---|--------|-------------|
| **D1** | The bacteriocin Green's function has **no decay length**. Degradation enters only as a time-dependent amplitude on the whole burst, never as `exp(-r/l)`. | The pI → retardation → *Lethal Core* vs *Lethal Halo* mechanism — the model's central biophysical claim — cannot be produced by this kernel. Retardation rescales amplitude only, so the *more* retarded "core" toxin has the *larger* field at every radius. Backwards. |
| **D2** | Receptor occupancy is **saturated by 4–8 orders of magnitude** everywhere inside the evaluation cutoff. | The lethal radius equals `toxin_cutoff` (200 µm), a numerical truncation parameter. The toxin-scape is a binary step function with no spatial gradient, so geometry, density and clustering are mathematically incapable of mattering. |
| **D3** | Spatial structure is **not measured** and there is **no exogenous immigration**. | The C++ spatial scalars are hardcoded placeholders; the Python patch metric is estimated from ~1% of agents at realized density. Endogenous diversity generation runs at ~0.8 events/day. |

Once D1 and D2 are fixed, the model's *own* plasmid library already places it in
the one interesting regime — screening length comparable to inter-colony spacing —
which is precisely where all three empirical observations can co-occur (§2).

The tractability answer (§3) follows from the same diagnosis: the unit of ecology
is the **microcolony**, not the cell; the week-long run is a *compound* of one
elementary encounter event; and the QSSA toxin grid — which dominates Stage 3
cost and memory — is an analytic function of agent positions and never needs to
be stored or, for most questions, simulated at all.

All numbers below are derived from source, with citations, and were independently
re-derived a second time before publication. Where a quantity is a conditional
upper bound or an illustrative scale rather than a measured model output, it is
labelled as such — that distinction matters here, because several of the claims
below assert that committed results are artifacts.

---

## 1. The three defects, with numbers

### D1 — The kernel has no length scale

`GreensFunction::single_kernel()` implements
([`src/diffusion/greens_function.cpp:182-204`](../src/diffusion/greens_function.cpp)):

```
C(r) = Q / (4 pi D_eff r) * exp( (U.r - |U| r) / (2 D_eff) ),    D_eff = D_free / R
```

There is no `exp(-r/l)` term. Degradation exists, but only as a scalar amplitude
on the whole burst, applied uniformly in space
([`src/diffusion/qssa_solver.cpp:288-305`](../src/diffusion/qssa_solver.cpp)):

```cpp
factor = std::exp(-burst.decay_rate * age);   // function of age, not of r
```

A point source releasing at rate `Q` into a medium with first-order degradation
`lambda` has the **screened (Yukawa)** steady state, not the Coulomb one:

```
C(r) = Q exp(-r/l) / (4 pi D_eff r),     l = sqrt(D_eff / lambda)
```

`l` is the length that separates a core from a halo. (It is not the model's only
length scale — the implemented system also has advective lengths, the 200 µm
toxin cutoff, the grid spacing, and the domain dimensions — but it is the only
one that is *biophysical*, and it is absent.) Dropping it has two consequences:

1. **The core/halo dichotomy cannot emerge.** With a pure `1/r` kernel the only
   effect of retardation is the `1/D_eff` prefactor, so raising retardation
   raises the concentration *at every radius equally*. ColE1 (`R = 50`) therefore
   projects a **larger** field than ColB (`R = 1.5`) at all distances — the exact
   inverse of the intended "basic toxins bind mucin and stay local" biology
   ([`src/fixes/fix_bacteriocin.h:14-15,31-33`](../src/fixes/fix_bacteriocin.h)).
2. **The information needed to fix it is already in the repo.** The plasmid
   library already carries per-toxin `D_free`, `retardation`, and
   `protease_half_life` (ColE1 [`src/genome/plasmid.cpp:29-43`](../src/genome/plasmid.cpp),
   ColB [`:65-79`](../src/genome/plasmid.cpp), MccV [`:115-128`](../src/genome/plasmid.cpp)).
   Combining each toxin's effective diffusivity with its own protease half-life
   gives exactly the intended separation for free. **These are illustrative
   scales, not implemented screening lengths** — nothing in the code computes
   them today:

| Toxin | `D_free` (m²/s) | `R` | `D_eff` (m²/s) | half-life (s) | **`l = sqrt(D_eff/lambda)`** | Character |
|-------|------|-----|--------|-----------|----------|-----------|
| ColE1 | 4e-11 | 50 | 8.0e-13 | 1800 | **46 µm** | Lethal **core** |
| ColB | 4e-11 | 1.5 | 2.7e-11 | 900 | **186 µm** | Lethal **halo** |
| MccV | 1e-10 | 1.2 | 8.3e-11 | 7200 | **930 µm** | Field-scale |

That is the core/halo dichotomy, quantitatively, from parameters already
committed. It is simply not wired into the kernel.

**Comet tails have the same problem.** Advective distortion requires
`Pe = |U| r / (2 D_eff) >~ 1`. With `v_max ~ 2.3e-8 m/s`
([`src/fields/advection.cpp:21-29`](../src/fields/advection.cpp)):

| Toxin | `r` at `Pe = 1` (distal `v_max`) | using radial `v_max` (1.85e-8 m/s) | Comet tail visible in a 2 mm domain? |
|-------|------|------|--------------------------------------|
| ColE1 | 69 µm | 86 µm | Marginal — comparable to the NND |
| ColB | 2.3 mm | 2.9 mm | No — larger than the domain |
| MccV | 7.2 mm | 9.0 mm | No |

The kernel actually uses the full local velocity vector `(v_x, 0, v_z)`, whose
magnitude varies with `z` under the `z^1.5` profile
([`src/fields/advection.cpp:43-50`](../src/fields/advection.cpp),
[`src/diffusion/greens_function.cpp:190-204`](../src/diffusion/greens_function.cpp)),
so both columns are bounding estimates rather than a single exact radius. The
conclusion is insensitive to the choice.

So comet tails are a *core*-toxin phenomenon at these flow rates, and are
essentially absent for halo toxins. The committed golden
`comet_tail_ratio = 0.94` ([`python/tests/fixtures/eari_vadi_ci_golden.json`](../python/tests/fixtures/eari_vadi_ci_golden.json))
is consistent with this, and is currently being locked in as "correct".

### D2 — Occupancy is saturated, so the kill radius is a numerical parameter

The kill law has no threshold and no Hill term; it is hazard accumulation on
receptor occupancy ([`src/fixes/fix_receptor.cpp:91-204`](../src/fixes/fix_receptor.cpp)):

```
occupancy = expr * [T] / (Kd_app + [T]),   H += kill_rate * occupancy * immunity * dt,   P_kill = 1 - exp(-H)
```

Now put the actual source strength in. A lysing ColE1 cell contributes
`Q = colicin_release_rate * burst_size / burst_molecules = 1e-18 * 1e5/1e4 = 1e-17 mol/s`
([`src/fixes/fix_bacteriocin.cpp:118-142`](../src/fixes/fix_bacteriocin.cpp),
[`src/diffusion/qssa_solver.h:40-42`](../src/diffusion/qssa_solver.h),
[`src/fixes/fix_bacteriocin.h:39`](../src/fixes/fix_bacteriocin.h)), against
`kd_colicinE_btuB = 5e-10 mol/m³` ([`src/fixes/fix_receptor.h:26-51`](../src/fixes/fix_receptor.h)):

Evaluating the kernel with zero flow (equivalently, exactly downstream-aligned,
where the advection exponent vanishes):

| Distance from **one** producer | `C` (mol/m³) | `C / Kd` | Occupancy |
|---|---|---|---|
| 10 µm | 9.9e-2 | 2.0e8 | 1.000000 |
| 50 µm (≈ realized NND) | 2.0e-2 | 4.0e7 | 1.000000 |
| 200 µm (cutoff) | 5.0e-3 | 1.0e7 | 1.000000 |

The radius at which `C` would fall to `Kd` is **1989 metres**. Every toxin is
saturating at the realized NND: ColB reaches `3.0e4 * Kd`
(`kd_colicinB_fepA = 2e-9`) and even continuously secreted MccV reaches
`64 * Kd` (`kd_colicinIa_cirA = 3e-9`, its CirA target,
[`src/fixes/fix_receptor.cpp:137-159`](../src/fixes/fix_receptor.cpp)).
The upstream advection exponent modulates this by at most ~3× for ColE1 and
negligibly for the others — invisible against a 10⁴–10⁸ margin.

Therefore, for a non-immune cell anywhere within the 200 µm cutoff, with unit
receptor expression and affinities and negligible ligand competition:

```
H = kill_rate_colicin * 1 * 1 * dt = 1e-3 * 60 = 0.06 per step  ->  P_kill = 5.8% per minute
```

A saturated-hazard median lifetime of ~12 minutes, **independent of distance**,
and abruptly zero hazard 200.001 µm away. (Occupancy is not unconditionally 1 —
it still depends on receptor expression, toxin/ligand affinities and ligand
competition, [`src/fixes/fix_receptor.cpp:91-204`](../src/fixes/fix_receptor.cpp) —
but the toxin term itself is saturated by 4–8 orders of magnitude, so occupancy
is no longer a function of *distance*.) The toxin-scape is a hard sphere whose
radius is the QSSA evaluation cutoff. This is why spatial structure cannot
matter: there is no gradient for geometry to interact with. It is also why the
Kd sweeps in the Stage 3 campaign cannot move anything — a 10× change in Kd moves
occupancy from 0.99999997 to 0.9999997.

Two independent errors stack here:

1. **A finite burst is modelled as a sustained source.** `1e5` molecules is
   `1.66e-19 mol`. At `Q = 1e-17 mol/s` the cell emits the entire stated burst
   every **17 ms**. The source is not literally permanent — bursts are pruned
   after five decay constants, ~3.6 h for ColE1
   ([`src/core/simulation.cpp:1422-1440`](../src/core/simulation.cpp)) — but over
   just one protease half-life it emits `1.1e10` molecules, **~10⁵× the stated
   burst size**. Molecule inventory is never conserved. A one-shot release should
   be a decaying transient whose *time-integrated dose* is
   `Q_total/(4 pi D_eff r)`, not a steady state fed by a near-inexhaustible
   source. Order-of-magnitude check: `1.66e-19 mol` diluted into a 50 µm-radius
   sphere is `3.2e-7 mol/m³`, versus the `2.0e-2` the code produces — a ~10⁵
   overestimate, matching.
2. **A possible unit slip in `Kd`.** Receptor Kd values are declared in `mol/m³`
   ([`src/fixes/fix_receptor.h:26-40`](../src/fixes/fix_receptor.h)), consistent
   with the chemical fields. But `5e-10`, `2e-9`, `3e-9` read like *molar*
   receptor affinities; as `mol/m³` they are `0.5 pM`–`3 pM`, ~10³ tighter than
   published TonB-dependent-transporter affinities. This is a **hypothesis about
   provenance, not a confirmed defect** — it needs checking against whatever
   source the values were taken from. Note that defect (1) alone is sufficient to
   produce the saturation.

**Side effect worth noting:** with saturated occupancy and the default
`immunity_binding_affinity = 1.0` ([`src/core/types.h:83-91`](../src/core/types.h)),
`immunity_factor = 1e-3` still yields `H = 6e-5` per step for an *immune*
resident — a **45% cumulative probability of self-death over a 7-day run**
(10080 steps at `bio_dt = 60`). This is a conditional upper bound on the same
saturation assumption, but it means the resident population may be substantially
killed by its own toxin as an artifact of D2.

### D3 — Structure is unmeasured; diversity has no exogenous source

**Spatial observables are non-functional on both sides.**

The C++ summary writes literal constants
([`src/io/hdf5_writer.cpp:551-562`](../src/io/hdf5_writer.cpp)):

```cpp
if (live >= 10) {
  mean_nnd = 5.0e-6;      // placeholder
  hopkins  = 0.5;         // placeholder
  mono     = 0.5;         // placeholder
}
```

Every `spatial/*` scalar in every HDF5 file is a placeholder, not a measurement.

The Python metric is scale-mismatched to the realized density.
`monochromatic_patch_score` probes a **fixed 10 µm radius**
([`python/gut_ibm_tools/analysis.py:120-143`](../python/gut_ibm_tools/analysis.py)),
but the nearest-neighbour distance is 20–40 µm. Agents are placed uniformly in
`x`/`y` but only in the **lower half of `z`**
([`src/core/simulation.cpp:344-350`](../src/core/simulation.cpp)), so the occupied
slab is 50 µm, not 100 µm. For a Poisson cloud at the two campaign densities:

| Config | cells/mm³ (occupied slab) | Poisson NND | Expected agents with ≥1 neighbour within 10 µm |
|---|---|---|---|
| `eari_vadi_validation` (160 cells, 0.5 mm box) | 12,800 | 21.9 µm | ~8 / 160 (5%) |
| Stage 3 (600 cells, 2 mm box) | 3,000 | 38.4 µm | ~8 / 600 (1.3%) |

So at Stage 3 the patch metric is estimated from roughly **1% of the population**
and ignores the other 99%. That alone disqualifies it as a structure measure.

The committed golden is worse than that: it is exactly
`monochromatic_score: 0.0`
([`python/tests/fixtures/eari_vadi_ci_golden.json`](../python/tests/fixtures/eari_vadi_ci_golden.json)).
That value is *only reachable through the empty-sample paths*: the neighbour
query includes the focal agent itself, so any agent that survives the
`len(neighbors) <= 1` guard contributes a strictly positive fraction, and a
non-empty `fractions` list cannot average to zero
([`python/gut_ibm_tools/analysis.py:130-143`](../python/gut_ibm_tools/analysis.py)).
Since the companion FISH golden records 69 detected agents, the population was
not empty. Therefore **no agent in that run had a neighbour within 10 µm**, and
the empty-sample sentinel has been regression-locked as a scientific result.

Separately, the metric's null is not 0.5 and is composition-sensitive: because
the focal agent is counted in its own neighbourhood and isolated agents are
dropped, a purely random cloud returns something near the majority-type fraction
(0.76–0.79 in 75/25 test realizations; Stage 3 is 500/600 ≈ 83/17). It therefore
conflates community composition with spatial segregation and cannot isolate
structure.

**Diversity generation is far too slow, and there is no immigration.** Mutation
fires only on division ([`src/fixes/fix_mutation.cpp:15-60`](../src/fixes/fix_mutation.cpp)).
At Stage 3, taking the optimistic upper bound in which every cell divides
continuously at its configured `mu_max` (**40,511 divisions/day**) and using the
campaign's own `bi_duplication_rate = 2e-5`
([`experiments/diversity_campaign/stage3_campaign/3a_baseline.json:15-17`](../experiments/diversity_campaign/stage3_campaign/3a_baseline.json)):

| Event | Rate | Per day | Per 7-day run |
|---|---|---|---|
| BI duplication | 2e-5 | **0.810** | ~5.7 |
| BI recombination | 5e-6 | **0.203** | ~1.4 |
| Receptor mutation | 1e-7 | **0.0041** | ~0.03 |
| Super-killer emergence | 1e-8 | **4.05e-4** | one per **6.8 years** |

Realized rates will be *lower*, since growth limitation, death and washout all
reduce division below `mu_max`.

And there is **no exogenous immigration anywhere**: "immigrant" is only an
initial strain type; nothing arrives after `t = 0`, and a checkpoint fork cannot
inject agents (§3.4). Endogenous diversity *can* increase — mutation
([`src/fixes/fix_mutation.cpp:29-60`](../src/fixes/fix_mutation.cpp)) and
conjugation ([`src/fixes/fix_conjugation.cpp:61-81`](../src/fixes/fix_conjugation.cpp))
both generate new BI states — but at under one event per day it cannot shape a
7-day outcome.

So standing diversity is effectively a fixed initial condition. "Resident
retention" is not a diversity dynamic; it is a survival statistic on a closed
founder set. Any claim about within-host diversity *maintenance* or about
across-host diversity is currently unsupported by construction.

---

## 2. What coherence actually requires

Strip away the defects and the science question is a single dimensionless
comparison between two lengths:

- **`l`** — the toxin screening length, `sqrt(D_eff / lambda)` (D1), or more
  precisely the dose-based lethal radius `r_L(N)` around a colony of `N` producers
- **`d`** — the inter-**colony** spacing, set by density and clustering

| Regime | Ecology | Diversity | Spatial signature |
|--------|---------|-----------|-------------------|
| `l << d` | Colonies never interact | High, neutral drift | Random (Hopkins ≈ 0.5) |
| `l >> d` | Effectively well-mixed antagonism | Collapses to the best killer | None — no structure to see |
| **`l ~ d`** | **Percolation-like; interaction is patchy and history-dependent** | **Low but non-zero; strong priority effects** | **Monochromatic patches, exclusion zones** |

Only the middle regime reproduces all three empirical observations at once. And
with the repo's own parameters, ColE1 sits at `l = 46 µm` against a Stage 3
Poisson NND of 38 µm: **`l/d ~ 1`**. The model is already in the interesting
regime (for single cells; §2.1 shows the colony-scale comparison is the one that
matters, and it is more favourable still). It
cannot express it because D1 removed `l` and D2 replaced `r_L` with the cutoff.

### 2.1 The unit of ecology is the microcolony, not the cell

This is the resolution of the "low density can't work" objection, and it should
be stated explicitly because it inverts the intuition.

A single lysing cell cannot deliver a lethal dose to a neighbour 50 µm away —
that is the correct physical conclusion from the corrected burst model, and it is
why the low-density objection *feels* fatal. But time-integrated dose from a
co-located clone scales linearly in `N`:

```
Dose(r) = N * Q_total * exp(-r/l) / (4 pi D_eff r)
```

so a lethal radius comparable to `l` requires `N ~ 10²–10³` co-located producers.
Cells reach that within hours at a 21-minute doubling time.

**Clustering is therefore not an obstacle to bacteriocin efficacy — it is the
precondition for it.** The apparent paradox comes from evaluating efficacy at the
mean-field single-cell density (where it is indeed impossible) rather than at the
colony scale (where it is comfortable). The correct state variable is the
microcolony: position, size, BI genotype, immunity set. Everything downstream
follows from treating it that way.

---

## 3. Making it tractable

The reframing pays for itself computationally. Once the colony is the unit, a
7-day host-scale run is a **compound of one elementary event**: a propagule
arrives at distance `d` from an established colony and is either intercepted
(killed), washed out, or establishes. Measure the per-encounter outcome on short
forks; compose it analytically over a week; spend the rare expensive run on
*validating the composition*, not on exploring parameters.

Four instruments, in increasing cost. Note that the first two involve **no
simulation at all**.

### 3.1 T0 — Toxin-scape probe (zero simulation, seconds)

A Python tool (`gut_ibm_tools/toxin_scape.py`) that reads any existing agent dump
or Tier-2 checkpoint and evaluates the analytic kernel to report:

- **`phi_L`** — the lethal-dose volume fraction of the domain
- the **dose distribution** at every non-producer's position
- **`p_int = phi_L`** — the interception probability for a uniformly-arriving immigrant
- the realized **inter-colony NND distribution**, from a DBSCAN colony catalog

This is the single number that connects density, structure and efficacy, and it
is computable on states already on S3, in seconds, for free. It is also the
diagnostic that would have caught D2 on day one: today it returns `phi_L = 1.0`.

Run this **first**, on the existing `burnin_seed1001` checkpoint, before changing
any physics.

### 3.2 T1 — Kernel calibration (zero simulation, seconds)

After wiring the Yukawa term and the finite-burst source, sweep
`(D_free, retardation, protease_half_life, burst_size, Kd)` purely analytically
and report `l` and `r_L(N)` for `N = 1, 10, 100, 1000`. Accept a parameter set
only if:

- ColE1 is a core and ColB a halo (`l_ColE1 < l_ColB`, separated by ≳3×)
- `r_L(N ~ 100)` lands within a factor of ~3 of the realized inter-colony NND
- a single cell (`N = 1`) is **not** lethal at the NND

No IbM run is needed to satisfy any of these. This is a `docs/PARAMETERS.md`
calibration exercise plus a unit test, and it replaces the entire Stage 3 Kd
sweep as a tuning device.

### 3.3 T2 — Colony challenge assay (minutes each, hundreds of them)

The elementary event, as a designed experiment. Restore a burn-in, **inject** `k`
immigrants at controlled distances from established colonies, run 1–6 simulated
hours, and record each immigrant's fate — killed / washed out / dividing — against
its distance to the nearest producer colony and that colony's size.

Output: the **interception kernel** `p_kill(d, N_colony)` and the washout rate,
measured rather than assumed. Replicate across burn-in seeds and injection
placements. Each run is minutes; the ensemble is hundreds of runs for less than
one L3 run's cost.

This requires the one genuine code change: **post-restore agent injection at fork
time** (§3.4).

### 3.4 The one blocking capability gap

`apply_checkpoint_snapshot()` reconstructs the population strictly from the
checkpoint arrays ([`src/core/simulation.cpp:428-513`](../src/core/simulation.cpp),
[`src/io/hdf5_reader.cpp:116-243`](../src/io/hdf5_reader.cpp)). Changing
`initial_strains` in a fork overlay does **not** add agents. There is no
immigration mechanism anywhere in the stepping loop.

Needed: an `immigration` config block applied after checkpoint restore —
`{count, strain, placement: uniform|at_distance|z_slab, schedule: pulse|continuous}`.
This is small, and it unlocks both T2 and any claim about diversity maintenance
(§1 D3). Without it, forks can only explore parameters, never perturbations, and
the only diversity source is endogenous mutation/conjugation at under one event
per day.

### 3.5 T3 — Compound model (instant)

A colony-level birth / interception / washout process parameterized entirely by
T1 and T2 outputs plus the burn-in's colony spatial statistics. It predicts 7-day
retention and across-seed diversity spread by cheap Monte Carlo. The expensive L3
run is then spent testing this prediction on 2–3 seeds.

This is the actual answer to "we can't run it full length in large numbers": you
don't. You run it a few times to *falsify a cheap surrogate*, and if the surrogate
survives, you explore in the surrogate.

### 3.6 Observable fixes required for any of the above

None of T0–T3 can be evaluated with today's observables:

1. **Delete or implement the C++ placeholder spatial scalars.** Writing `0.5`
   into an HDF5 dataset named `hopkins_statistic` is worse than writing nothing.
2. **Make the patch metric scale-aware.** Replace the fixed 10 µm probe with the
   pair-correlation function `g(r)` / Ripley's `K` over 10–500 µm, and always
   report the composition-matched randomized null alongside it. A single scalar
   at a hardcoded radius cannot work across a 4× density range.
3. **Emit a colony catalog** (DBSCAN at `eps ≈ NND/2`): colony count, sizes,
   genotype purity, inter-colony NND. This is the state variable of §2 and it is
   currently not recorded anywhere.
4. **Log per-kill provenance** — killer colony ID, distance, accumulated dose — so
   interception is directly observed rather than inferred from population curves.

### 3.7 Storage: stop saving the grid

Stage 3 is dominated by ~50M grid cells and ~8 GB of chemistry per rank
([`experiments/diversity_campaign/README.md:157-160`](../experiments/diversity_campaign/README.md)).

For the QSSA species this is **pure waste**: the bacteriocin field is, by
construction, an analytic superposition of kernels over agent positions. It is a
deterministic function of the agent state. It never needs to be written to disk,
and for most questions it never needs to be evaluated on a grid at all — only at
the ~10³ agent positions where occupancy is actually read.

Retention policy that follows:

| Artifact | Cadence | Size |
|----------|---------|------|
| Colony catalog + spatial summary | High (every summary step) | O(10²) rows |
| Event stream (kills w/ distance + dose, divisions, washouts, lineage) | Every step | O(10³)/day |
| Agent dump (positions, genomes) | Moderate | O(10³) rows |
| Full grid | **Fork points only** | 8 GB |

That is the difference between "we cannot save everything" and "we save
everything that carries information".

---

## 4. Ordered plan

| Step | Work | Cost | Unblocks |
|------|------|------|----------|
| 1 | **T0 probe**, run on the existing `burnin_seed1001` checkpoint | Hours | Confirms D2 in your own data before any physics change |
| 2 | Fix the **burst source term** (finite dose, decaying transient) and add **Yukawa screening**; audit the `Kd` units | Hours | D1, D2. Goldens will move — re-anchor them deliberately |
| 3 | **Observable fixes** (§3.6): kill the placeholders, `g(r)`, colony catalog, kill provenance | Hours | Makes anything measurable |
| 4 | **Immigration / fork-time injection** (§3.4) | ~1 session | T2, and any diversity claim at all |
| 5 | **T1 sweep** + **T2 assay ensemble** across existing burn-in seeds | ~1 session | Calibrated `l`, `r_L(N)`, measured `p_kill(d, N)` |
| 6 | **T3 compound model**, validated against 2–3 existing 7-day runs | ~1 session | Retention and across-host diversity without new week-long runs |

No new burn-ins are required until step 6, and step 5's validation reuses the
long runs already completed. Steps 1–3 change no science, only what is computed
and recorded.

---

## 5. What this changes about existing results

Stated plainly, because it affects anything already reported:

- Any Stage 3 **Kd sweep** result is uninformative: occupancy is saturated across
  the entire swept range (§1 D2).
- Any **spatial structure** result is unmeasured: the C++ scalars are constants,
  and the Python patch score is estimated from ~1% of agents — and in the
  committed golden, provably from none of them (§1 D3).
- Any **core vs halo** or **comet tail** result is unproduced by the physics: the
  kernel has no length scale and halo-toxin Péclet numbers are ≪1 (§1 D1).
- Any **diversity maintenance** claim is unsupported: no exogenous immigration,
  and under ~6 BI events per 7-day run even on the optimistic bound (§1 D3).
- The committed EARI/VADI goldens encode several of these artifacts as expected
  values and must be re-anchored, not preserved, when steps 2–3 land.

Resident **retention** and **washout** results are not implicated by D1/D2 in the
same way, but they are inflated by the 45% self-death artifact noted in §1 D2 and
should be re-run after step 2.
