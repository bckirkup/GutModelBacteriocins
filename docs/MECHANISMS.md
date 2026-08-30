# Biological Mechanism Reference

Detailed descriptions of each Fix module, their biological basis, and implementation.

---

## Adaptive Timestep Selection

When `adaptive_dt_enabled = true`, `compute_adaptive_dt()` selects the biological timestep each iteration using CFL-like constraints:

```
dt = dt_max
dt = min(dt, dt_growth_limit / max(|mu_realized|))   // growth constraint
if (sos_count > 5)  dt = min(dt, 10)                  // lysis cascade
if (sos_count > 20) dt = min(dt, 2)
if (density > 1e15) dt = min(dt, 10)                  // overcrowding
dt = clamp(dt * dt_safety, dt_min, dt_max)
```

The adaptive dt is recomputed before every call to `step()`. During quiescent periods, `dt` rises toward `dt_max` (up to 300 s), improving throughput. During mass lysis or rapid growth, `dt` drops to resolve the fast dynamics (down to `dt_min` = 1 s). The safety factor (default 0.8) provides a margin below the stability limit.

The `run()` loop uses `while (time < total_time)` with the adaptive dt, clamping the final step to avoid overshooting `total_time`.

## Immigration

Immigration runs in the pre-step phase immediately after ghost removal and
before ghost exchange. This leaves the pool containing only real local agents,
so ghost indices remain valid and nearest-biomass distances count each resident
once. The new cell is then exchanged, coupled to the grid, and processed by
biology, chemistry, physics, post-step, migration, and washout in its first
step; injecting later would grant it a free step.

Pulse events fire when `current_step - run_start_step == immigration.step`;
therefore a checkpoint fork can use a step relative to the fork, independent
of the checkpoint's absolute step. Continuous events draw a Poisson number
with mean `rate * dt`.

`at_distance` proposes candidates on shells around replicated live-biomass
anchors, using isotropic directions and the requested distance. The global
minimum-image distance to the nearest live pre-existing agent is still
computed for every proposal, so a proposal near an unrelated cluster is
rejected unless its true nearest-biomass distance meets the tolerance. This
targets the colony surface because toxin exposure is dominated by nearby
producing cells, and remains well-defined for multiple colonies.
`distance_reference: centroid` instead anchors proposals on the globally
reduced centroid and is intended for a single compact cluster; its centroid is
not min-image-corrected across periodic wraps. `z_slab` samples uniformly
within its configured z bounds. If neither candidate batch meets
`distance_tolerance`, the cell is skipped with a rank-0 warning rather than
placed at the wrong distance.

Immigration uses a dedicated RNG seeded from `seed` with a fixed mixing
constant. Every rank draws identical event counts and candidates. The owning
rank alone constructs each cell, using its stride-allocated `AgentPool` tag;
this guarantees globally unique IDs without migration duplicates.

## Founder Placement

The default `initial_population.placement: legacy` policy is retained for
backward compatibility: it samples founder z positions uniformly from
`[domain.lo.z, 0.5 * domain.hi.z)`. The `z_slab` policy samples uniformly from
its explicit `initial_population.z_min` and `z_max` bounds.

The anatomy-derived `anatomic` policy is available for E. coli founders. It
keeps x and y uniform over the domain, and samples depth as

```
z = anatomic_exclusion_floor + Exponential(anatomic_exponential_scale)
```

with the draw rejected and resampled while `z >= anatomic_outer_extent`.
Resampling, rather than clamping, avoids creating an artificial probability
spike at the outer truncation point. The defaults are a 20 µm hard epithelial
exclusion floor, a 40 µm exponential scale, and a 150 µm outer extent. No
founder is marked as `in_crypt` under this policy, even when the configured
crypt geometry includes its sampled position.

These values are anatomy-derived rather than fitted parameters: the 20 µm
floor represents the Enterobacteriaceae-free epithelial gap observed by
Swidsinski et al. using FISH, while the 150 µm outer-mucus extent reflects the
enrichment described by Duncan and Mondragón-Palomino. The absence of crypt
founders is a data-driven constraint from healthy FISH observations, not an
implementation limitation. Immigration currently has separate placement
handling; applying the anatomy policy to reseeding is a future consistency
decision.

---

## Fix Architecture (NUFEB-inspired)

Each biological rule is encapsulated as a **Fix** — a modular computation unit called once per biological timestep. Default registration order:

1. `fix_metabolism` — Monod growth, Fur iron regulation, division, death
2. `fix_quorum_sensing` — AI-2 production / Lsr import (Spec 11)
3. `fix_bacteriocin` — SOS/phage/microcin release (Spec 2)
4. `fix_receptor` — Competitive binding at TBDTs, toxin-mediated killing
5. `fix_motility` — Run-and-reverse swimming, chemotaxis (Spec 3 / 10v2 / 11)
6. `fix_conjugation` — Horizontal gene transfer via F-pili
7. `fix_cdi` — Contact-dependent inhibition (Spec 3)
8. `fix_mutation` — BI locus evolution, receptor downregulation, super-killers
9. `fix_mechanics` — Hertzian repulsion, optional EPS adhesion

The QSSA solver deposits bacteriocin fields immediately before `fix_receptor` (so killing sees the current step's toxin sources), then runs the full chemistry pass (nutrients, VBF, reactions) after biology. Advection operates in the physics module. Motility displacement is applied in `module_physics()` after advection.

---

## 1. fix_metabolism — Triple Monod Growth

### Delivery-limited carbon uptake

With `metabolism.uptake_limit: "delivery"`, requested growth and maintenance
carbon are held until chemistry completes. This mode is CPU-only:
`gpu_enabled=true` together with `metabolism.uptake_limit="delivery"` is
rejected during config finalization because CUDA parity is not implemented
yet. Each owned agent contributes an analytic prescribed mass
`min(demand_total, Sherwood)` to its cell. The aggregate is applied as a
zero-order right-hand-side draw, `-F_cell/(3*V_cell)`, in each directional
backward-Euler solve. This is not the retired legacy zero-order sink: that sink
was applied and clipped to voxel content before diffusion, while this sink is
inside implicit diffusion and is capped at analytic diffusive delivery rather
than unbounded demand. The implicit solve therefore supplies the diffusive
neighbourhood rather than restricting the draw to the agent's voxel, with no
conductance or `k/3` splitting bias. Realized removal equals funded uptake by
construction. If a genuinely starved neighbourhood would become negative,
local prescribed draws are reduced first and the solve is retried. Each local
retry halves prescribed mass in the physical delivery-radius neighbourhood of
every negative owned cell, using periodic x/y and clipped z support semantics.
The local pass runs for up to four attempts and stops early when no affected
prescribed value can be changed. When more than one quarter of owned cells are
negative, the local pass is skipped: dilation of a domain-wide deficit covers
the domain and would reproduce the global rationing at higher cost. If local
reductions do not restore positivity, the original prescribed field is
restored and one global scalar factor is bisected for 12 iterations; the
largest feasible factor is applied uniformly as the final guarantee. MPI
feasibility decisions are collective, including the negative-cell fraction,
so every rank performs the same solves. The delivery-reduction ledger is exactly
original prescribed mass minus final prescribed mass, and retry events include
local and bisection solves. The emitted rationing factor is the minimum
per-cell factor over owned cells in the interval, not a mean. If non-delivery
sources still leave a negative owned cell at factor zero, the run continues,
emits a species/step/minimum/count diagnostic, and increments the per-species
infeasible-step counter. Ghost agents do not contribute prescribed draws or
accounting writes.

### Delivery closure auditing

Delivery runs enable a closure guard by default. After each chemistry step,
globally reduced carbon demand and realized removal are compared; consecutive
steps with positive demand and zero realized removal beyond
`closure.zero_realization_grace_steps` terminate with
`closure_violation`. Set `closure.enforce_delivery_realization` to `false` to
disable this guard. The optional
`closure.enforce_reaction_clip` guard compares interval reaction clipping with
`closure.reaction_clip_tolerance_fraction`; it defaults to `false` because
legacy `uptake_limit: "none"` intentionally clips substantial requested
carbon. Closure decisions are synchronized across MPI ranks and a closure
violation returns a nonzero process status.

Every run records one authoritative termination cause under
`/run_provenance/termination_cause_code` and
`/run_provenance/termination_cause`:

| Code | Cause |
|---:|---|
| 0 | `horizon_reached` |
| 1 | `dysbiosis_guard` |
| 2 | `population_stop` |
| 3 | `stop_requested` |
| 4 | `closure_violation` |
| 5 | `incomplete_unknown` |

The pessimistic `incomplete_unknown` marker is flushed when provenance is
created and is overwritten only when `Simulation::run()` reaches a recorded
termination. Thus a crash or OOM kill remains distinguishable from a completed
run. `termination_detail` carries a short diagnostic string and
`termination_wall_seconds` records elapsed wall time; existing termination
step/time datasets retain their historical meanings.

**Biological basis:** *E. coli* growth in the gut requires carbon (mucin-derived monosaccharides), iron (via siderophores through multiple TBDTs), and vitamin B12 (via BtuB). Growth rate follows multiplicative Monod kinetics. When `oxygen.enabled`, the legacy aerobic boost applies unless `oxygen.metabolic_switch_enabled` is true.

**Equation:**
```
mu = mu_max * [C]/(Km_C + [C]) * monod_iron * [B12]/(Km_B12 + [B12])
```

**Aerobic boost (Spec 1, when `oxygen.enabled`):**
```
monod_O2_boost = 1 + boost_max * [O2] / (Km_O2 + [O2])
mu *= monod_O2_boost
```
Agents consume O₂ proportional to `mu_realized`. `Simulation::local_O2()` and `ros_induction_rate()` expose local O₂ for Spec 2 SOS coupling.

When the metabolic switch is enabled, it supersedes the aerobic boost. Oxygen
availability limits respiratory capacity, and the complementary fraction is
fermentative; high growth at high oxygen can therefore enter aerobic overflow.
The realized fermentation fraction relaxes over `oxygen.tau_metabolic_switch`
and is inherited at division. It controls growth multiplier, carbon cost,
growth-associated oxygen consumption, fermentative acetate, and the carbon
maintenance factor. Basal oxygen maintenance remains unconditional. The
fermentative acetate term replaces the legacy overflow term while scavenging
is unchanged.

When `metabolism.uptake_limit` is `delivery` and
`oxygen.delivery_uptake_enabled` is true, growth-associated and maintenance
respiration use the same per-agent implicit delivery sink as carbon. Realized
oxygen removal is distributed in proportion to per-agent respiratory demand;
the explicit QSSA oxygen sink remains the default when the flag is false.

`oxygen.respiration_driver` selects how that realized fermentation fraction is
updated. The default `ambient` driver preserves the concentration-based
behavior above. With `funded`, the fraction is driven by the funded growth-O₂
capacity:

```
capacity = clamp01(funded_growth_o2 / demanded_growth_o2)
instantaneous_fermentation = 1 - capacity
```

Funded mode requires `oxygen.delivery_uptake_enabled=true` and
`metabolism.uptake_limit="delivery"`; it is rejected otherwise because no
funded quantity exists. The mode is **one step lagged**: this step's growth-O₂
demand is computed from the previous step's fermentation fraction, and this
step's realized delivery updates the fraction for the next step. If demanded
growth O₂ is zero, the ratio is undefined and the fraction is left unchanged.
Funded respiration is CPU-only and is rejected with GPU execution.
This switch does not change `oxygen.epithelial_conc`.

The oxygen epithelial boundary is configured with
`oxygen.epithelial_boundary`, `oxygen.epithelial_transfer_coeff`, and
`oxygen.epithelial_flux` (or their underscore aliases). The default
`dirichlet` mode preserves the historical imposed epithelial concentration.
In `robin` mode, `oxygen.epithelial_conc` is the tissue-side reservoir
concentration, not the imposed mucus-surface concentration; the surface value
emerges from vascular delivery minus epithelial and bacterial consumption.
The sourced experiment values are `k = 1.2e-6 m/s` for transfer and
`C_tissue = 5.0e-2 mol/m³` for the tissue-side reservoir, and are not new
defaults. The imposed 25 µm exponential oxygen profile must be disabled with
`oxygen.z_gradient=false` (or `oxygen_z_gradient=false`) in Robin mode.

Acid inhibition uses undissociated acetate via Henderson–Hasselbalch. Its pH
comes from `bacteriocin.mucin_charge.ph`, the environment pH lever introduced
for pI-derived retardation; no duplicate pH parameter is used. At pH 6 the
undissociated fraction is approximately 5.4%, so physiological colonic
acetate is generally below half-inhibition and the mechanism mainly acts in
acidified microenvironments.

See [`SPEC12_DENSITY_LIMITATION.md`](SPEC12_DENSITY_LIMITATION.md) for the
amended Spec 12 mechanism and implementation contract.

See the proposed multi-scale architecture (not yet implemented) in
[`SPEC13_MULTISCALE.md`](SPEC13_MULTISCALE.md) and [`SPEC13_IMPLEMENTATION_REVIEW.md`](SPEC13_IMPLEMENTATION_REVIEW.md).

**Graded iron uptake (Issue #10):** Rather than relying solely on FepA, iron acquisition uses four receptor systems in parallel with different affinities:

| Receptor | Siderophore | Km (nM) | Role |
|----------|-------------|---------|------|
| FepA | Enterobactin | 10 | Primary (highest affinity) |
| IroN | Salmochelin | 50 | Secondary (glycosylated enterobactin) |
| IutA | Aerobactin | 100 | Secondary (hydroxamate) |
| Fiu | Catecholates | 200 | Tertiary (broad specificity) |

```
iron_uptake = Σ expr_i * [Fe]/(Km_i + [Fe])   for i ∈ {FepA, IroN, IutA, Fiu}
monod_iron  = iron_uptake / (1 + expr_IroN + expr_IutA + expr_Fiu)
```

This replaces the previous binary FepA-dependent penalty (`Km_Fe / expr_FepA`). When FepA is downregulated to resist colicin B/D, cells switch to secondary systems rather than complete iron starvation. The normalization ensures wild-type cells (all receptors at 1.0) maintain equivalent growth.

**Note:** FhuA (ferrichrome) is NOT included as a secondary iron fallback because it transports fungal ferrichrome, which is not an endogenous enterobactin pathway (corrects EARI §70). Ferrichrome is represented only as an optional ambient/boundary ligand field and is disabled by default.

### Apo-enterobactin → ferric enterobactin → iron

Siderophore chemistry is enabled by default because enterobactin production is
near-universal among iron-limited *E. coli*. `siderophore` represents
apo-enterobactin.
Chelation consumes apo-enterobactin and free iron and produces the distinct
`ferric_enterobactin` species. FepA-mediated reimport consumes ferric
enterobactin and returns iron to the extracellular field, preserving the
tracked reaction chain. Secretion and reimport are specific rates in
mol/(s·kg) multiplied by biomass density. The secretion rate is a constrained
estimate; the reimport capacity is grounded in approximately 35,000 FepA
copies per iron-starved cell and approximately five TonB-limited transport
cycles per minute per FepA (Smallwood et al. 2016; Newton et al. 2010),
converted using the default cell mass. The exact conversion is
`8.33e-6 mol/(s·kg)` for the code-derived default cell mass, rounded to the
configured `1e-5 mol/(s·kg)`. The former secretion-rate-proportional recapture
term has been removed.

This chemistry does not make an isolated cell a meaningful FepA competitor.
With diffusion enabled, the measured local source-cell FeEnt concentration
immediately after the reaction update in a maintained single-cell assay was
`1e-11 mol/m³` at `1e-4 mol/m³` iron and `6e-15 mol/m³` at `1e-8 mol/m³` iron,
versus `kd_enterobactin = 1e-6 mol/m³`. The earlier `8e-16 mol/m³` value was
the post-diffusion field average, not the reaction-stage source-cell pulse.

The corrected chemistry was measured at 1, 4, 16, and 64 agents co-located in
one grid cell with a maintained apo-enterobactin background:

| Iron condition | N=1 | N=4 | N=16 | N=64 |
|---|---:|---:|---:|---:|
| `1e-4 mol/m³`, apo `4e-9 mol/m³` | `5e-12` | `1e-12` | `3e-13` | `7e-14` |
| `1e-8 mol/m³`, apo `3e-8 mol/m³` | `1e-15` | `3e-16` | `8e-17` | `2e-17` |

Values are reaction-stage FeEnt in `mol/m³`. Chelation is a solution-phase
reaction and is applied in every grid cell containing apo-enterobactin and
iron, independent of occupancy. Secretion and FepA reimport remain
biomass-weighted and confined to occupied cells. Thus increasing local
biomass increases the aggregate FepA sink without multiplying the
solution-phase chelation term, so local FeEnt decreases with co-location.
Diffusive escape further lowers the source-cell concentration. In the
EARI/VADI run, the final saved domain-wide FeEnt field had mean
`6.37e-8 mol/m³` and maximum `1.07e-7 mol/m³`, corresponding to competition
factors of approximately `1.064` and `1.107` relative to
`kd_enterobactin = 1e-6 mol/m³`. This is a spatially flat background effect,
not local co-location-driven protection: it cannot generate spatial structure
in colicin B susceptibility. It is a global parameter shift wearing the
costume of a spatial mechanism.

Enabling this chemistry by default adds two full implicit diffusion solves per
biological step, measured at `+15`–`17 ms/step` of chemistry time
(`+44`–`50%`) against the same benchmark with siderophore chemistry disabled.

This is different from bacteriocin dose. Lysing producers release toxin but do
not reimport it, so toxin dose does scale with co-located producers; the
microcolony threshold from #213 is consequently a genuine density effect,
whereas local FeEnt competition remains unreachable even though the
domain-wide background produces a modest uniform detuning.

The FeEnt reimport step uses a positivity-preserving backward-Euler
Michaelis–Menten update. With the literature-grounded `Vmax` and the model's
60 s biological timestep, the low-concentration pseudo-first-order sink can
be tens of inverse seconds; an explicit `reac * dt` update would overshoot
and clamp FeEnt to zero. The implicit solve remains nonnegative and bounded
by the biomass-scaled `Vmax` for all tested timesteps from `1e-6` through
`600 s`.

**Penalties applied:**
- **BtuB loss** (expr < 0.5): Activates MetE pathway for B12-independent methionine synthesis. Base cost = `metE_penalty` (default 5%). The MetE cost is further amplified by local acetate concentration (see below). Additionally, concentration-dependent ethanolamine utilization loss applies:
  ```
  eut_effect = eut_max_penalty * [EA] / (eut_km + [EA])
  mu *= (1 - metE_eff - eut_effect)
  ```
  At homeostatic [EA] (~0.5 mM, `eut_km` = 0.1 mM) the eut penalty is ~8.3%. In inflamed gut ([EA] >> Km) it approaches `eut_max_penalty` (10%). With zero ethanolamine the penalty vanishes.
- **Acetate inhibition of MetE** (VADI §87): Colonic acetate (60–100 mM) inhibits the MetE enzyme. The effective MetE penalty is scaled via Michaelis-Menten kinetics:
  ```
  acetate_factor = 1 + (max_factor - 1) * [acetate] / (Km + [acetate])
  metE_eff = metE_penalty * acetate_factor
  ```
  At 80 mM acetate (Km = 40 mol/m³, max_factor = 2.5), the penalty doubles from 5% to 10%. This strengthens the Combinatorial Washout Trap: BtuB-downregulated cells face a larger proteome burden in the acetate-rich colon.
- **Plasmid maintenance**: 2% per BI locus (reduced by compensatory mutations, see fix_mutation). Capped at 10%.
- **Maintenance energy**: Subtracted from mu_realized after all growth terms.

**Division:** When biomass >= `division_threshold * initial_mass` (default 2x), cell divides. Daughter is offset by one cell diameter in a random direction.

**Death:** Biomass below minimum threshold → cell dies.

### Non-growth-associated carbon maintenance

When `carbon_maintenance_rate` is nonzero, each owner agent removes
`carbon_maintenance_rate × biomass × dt` mol of carbon from its voxel every
biological step. This sink is independent of the growth-rate
`maintenance_rate` tax and is applied to bacteriostatic or otherwise
non-growing agents as well as growing agents. Its budget follows the
configured `uptake_limit` model: `none` charges the full request, `sherwood`
applies the same quasi-steady diffusive delivery cap as growth uptake, and
`voxel` limits the request to instantaneous voxel stock with a per-voxel
reservation. The unfunded amount is recorded in mol as
`maintenance_shortfall_*`; the separately reported
`maintenance_limited_agents_*` fields count agent-steps that were capped.

---

## 2. fix_bacteriocin — Toxin Release (Spec 2)

**Three release modes** (`ReleaseMode` on each `BICluster`):

| Mode | Plasmids | Mechanism |
|------|----------|-----------|
| `SOS_LYSIS` | ColE1, ColE2, ColM | LexA-regulated suicide lysis |
| `PHAGE_LYSIS` | ColB, ColIa | Temperate prophage induction |
| `CONTINUOUS` | MccV | Secretion without lysis |

### a) Multi-trigger SOS induction (Group A colicins)

SOS probability per biological timestep:

```
rate_total = sos_basal_rate
           + (just_divided ? sos_lysis_prob / bio_dt : 0)
           + sos_cross_induction_rate * [nuclease_bacteriocin]
           + ROS_rate
p_sos = 1 - exp(-rate_total * dt)
```

- `just_divided` is set on parent and daughter during division in `fix_metabolism`, cleared at the start of each `Simulation::step()`
- Cross-induction uses a per-agent QSSA exposure sampled from only
  `is_nuclease` sources after MPI source exchange. It is receptor-agnostic and
  therefore works identically for grid and agent toxin evaluation without
  adding a chemical species. BI-locus immunity attenuates exposure from
  nuclease toxins using the minimum matching immunity factor.
- After SOS induction: 5-minute delay (`sos_timer = 300 s`), then lysis with per-colicin `burst_size`.

The ROS driver is selected by `oxygen.ros_driver`. The `ambient` driver
retains the historical term
`oxygen.k_ROS * local_O2 * mu_realized`, where `local_O2` is the ambient
voxel concentration and is not a funded uptake amount. Its shipped
coefficient default is `oxygen.k_ROS = 0.0`, so ambient ROS mortality is
disabled unless a positive coefficient is explicitly configured. This default
retires an uncited coefficient whose `k_ROS = 1e2` product with oxygen no cell
paid for accounted for 99.5% of cumulative SOS hazard and every population
loss in the oxygen series (#332). With
`oxygen.ros_driver = funded`, the ROS term is

```
ROS_rate = oxygen.k_ROS_respiratory
           * respired_oxygen_rate / biomass
```

This is the specific-flux normalization: `k_ROS_respiratory` has units of
kg/mol, and `respired_oxygen_rate` is the funded growth plus maintenance
respiratory oxygen flux in mol/s, committed by delivery-limited metabolism and
divided by the agent biomass amount. Alternatively, when
`oxygen.k_ROS_funded > 0`, funded mode uses the absolute-flux normalization:

```
ROS_rate = oxygen.k_ROS_funded * respired_oxygen_rate
```

Here `k_ROS_funded` has units of mol⁻¹ and the coefficient is calibrated
directly against funded oxygen per generation, without biomass or
`mu_realized` normalization. Funded ROS uses the value committed by the
previous step: metabolism commits funded oxygen in `post_chemistry`, while
`FixBacteriocin` runs during biology. This one-step lag is intentional and
accepted. Exactly one of `oxygen.k_ROS_respiratory` and
`oxygen.k_ROS_funded` must be positive in funded mode. The absolute form is
the normalization used for the per-generation calibration described in
`docs/PARAMETERS.md`.

Respiratory ROS is only one SOS/colicin induction route. Bile salts, nutrient
stress, and SOS-independent operon regulation are not represented explicitly,
so a fitted coefficient stands in for several routes. The intended validation
target is in-vivo colicin-mediated displacement kinetics at Spec 13 scale,
not a scalar lysis fraction.

The event ledger records the resulting death as `mortality_lysis`; this is
separate from the `sos_inductions` and `phage_inductions` counters.
The same death is also recorded in kill provenance with cause `LYSIS`.

### b) Phage-mediated lysis (Group B colicins)

Founder strains may configure receptor genotype with a per-strain
`receptor_expression` object mapping receptor names such as `BtuB` and `FepA`
to expression levels in `[0, 1]`. These values initialize active and genomic
receptor expression; values below `0.2` use the same resistant-state rule as
mutation. B12/corrinoid initial and boundary concentrations are set together
with `b12.initial_conc` (or its flat and `corrinoid` aliases).

BI transport is controlled by each plasmid's library entry. The optional
`plasmid_overrides` object can override `retardation`, `diff_coeff`, and
`burst_size` per plasmid. Canonical names and legacy aliases resolve through
`PlasmidLibrary::find()`, and transferred BI clusters retain the donor's
configured values. The former global pI retardation keys are removed:
per-plasmid retardation is authoritative. Unknown keys are hard configuration
failures by default; `GUTIBM_STRICT_CONFIG=0` explicitly restores the
warning-and-ignore behavior for exploratory configurations.

ColB and ColIa use prophage induction, not SOS suicide:

```
rate_per_s = phage_induction_rate / (ln(2) / mu_realized)
```

On induction: `sos_timer = 60 s` (shorter lytic cycle), then burst release.
The resulting death is included in `mortality_lysis`, while induction remains
counted separately as `phage_inductions`.
It is likewise recorded in kill provenance with cause `LYSIS`.

### c) Continuous microcin secretion

MccV (`CONTINUOUS` mode) exports peptide without lysis:
- Applies a one-per-step growth penalty (`microcin_mu_penalty`, default 3%)
- Contributes steady-state QSSA point sources

**Per-colicin burst sizes** (default library): ColE1 1e5, ColE2 5e4, ColB/ColIa 1e4, ColM 2e5.

Lysis releases a finite inventory of `burst_size / AVOGADRO` mol with an
exponential `burst_release_tau` timescale. The source integral therefore equals
the configured burst inventory.

### Population stocks

Summary output includes two instantaneous live-population stocks:

- `bacteriostatic_live_agents`: live agents with
  `mu_realized < fixes.metabolism.bacteriostasis_threshold`;
- `washout_trapped_live_agents`: live agents with
  `mu_realized < washout_rate(z)`, evaluated at each agent's own z.

These are instantaneous stocks, not cumulative counters. They can fall as well
as rise, are not part of population closure, and must never be summed over
time. Bacteriostasis is a viable, non-growing state: falling below
`bacteriostasis_threshold` classifies a viable, non-reproducing state; it does
not kill a cell. Such a cell can leave only through outflow or an explicit
mortality mechanism.

Washout has two labelled modes. In the default `emergent` mode, the
washout-trapped predicate is an observation only; transport must carry the
agent to the luminal boundary before `outflow_boundary` is recorded. The
resulting `washout_trapped_live_agents` stock is therefore meaningful. In
`imposed` mode, the predicate removes non-crypt agents immediately and records
`outflow_washout`, so that stock is approximately zero for those agents.
Neither stock is a closure term.

**pI-dependent mucin retardation:**
| Class | pI range | Retardation | Behavior |
|-------|----------|-------------|----------|
| Lethal Core | > 8.5 | High R | Binds mucin glycoproteins, concentrates near producer |
| Lethal Halo | < 7.0 | Low R | Repelled by anionic mucin, spreads widely |
| Neutral | 7.0–8.5 | Intermediate R | Intermediate diffusion |

Retardation is derived from net charge using the configured local pH:

`R(pI) = r_min + amplitude / (1 + 10^((dz_half - (pI - pH)) / width))`

With the defaults (`r_min=1.2`, `amplitude=60`, `dz_half=1.35`,
`width=1`, `pH=7`), higher pI at fixed pH produces higher retardation,
whereas higher pH at fixed pI produces lower retardation. The
`classify_by_pI()` function remains the single source of truth for the
Lethal Core/Halo class thresholds; class labels and transport retardation are
resolved together when a BI cluster is created or mutated.

*Note:* The HALO threshold was raised from pI < 6.0 to pI < 7.0 (VADI §74) to correctly classify bacteriocins secreted as acidic complexes (e.g. Colicin E2 with Im2, net complex pI ≈ 6.5).

The effective diffusion coefficient: `D_eff = D_free / R`

---

## 3. fix_receptor — Competitive Binding at TBDTs

**The Double-Bind hypothesis:** TonB-dependent transporters serve dual roles as nutrient importers AND bacteriocin receptors. Cells cannot selectively block toxins without also losing nutrient uptake.

**Four receptor systems modeled** (each reads its own QSSA toxin field by
default):

| Receptor | Nutrient ligand | Toxin field | Toxin ligand | Kd_nutrient | Kd_toxin |
|----------|----------------|-------------|--------------|-------------|----------|
| BtuB | Vitamin B12 | `bacteriocin_BtuB` | Colicin E1/E2/E5/E7/E8/E9/K | 1 nM (`1e-6 mol/m³`) | 0.5 nM (`5e-7 mol/m³`) |
| FepA | Enterobactin-Fe | `bacteriocin_FepA` | Colicin B/D | 1 nM (`1e-6 mol/m³`) | 2 nM (`2e-6 mol/m³`) |
| CirA | Linearized siderophore | `bacteriocin_CirA` | Colicin Ia, Microcin V | 50 nM (`5e-5 mol/m³`) | 3 nM (`3e-6 mol/m³`) |
| FhuA | Ferrichrome | `bacteriocin_FhuA` | Colicin M | 10 nM (`1e-5 mol/m³`) | 2.5 nM (`2.5e-6 mol/m³`) |

CirA ligand is `cirA_linearized_fraction × [ferric_enterobactin]` when
siderophore chemistry is enabled; otherwise zero. FepA competition likewise
uses ferric enterobactin, but diffusion-limited local FeEnt remains far below
its Kd and falls further with co-location; the only competition present is the
uniform domain-wide background described above.
Ferrichrome remains disabled by default because
*E. coli* does not synthesize it and no defensible gut concentration is
available; consequently FhuA has no ambient ferrichrome competition by default.

**Competitive binding model (Michaelis-Menten with competitive inhibition):**
```
Apparent_Kd = Kd_tox * (1 + [Ligand]/Kd_ligand)
Occupancy = receptor_expr * [Tox] / (Apparent_Kd + [Tox])
P(kill) = 1 - exp(-kill_rate * Occupancy * dt)
```

When nutrients are scarce, receptors are unoccupied and "unlocked" — maximum toxin sensitivity. When nutrients are abundant, competitive inhibition protects cells.

**Immunity and affinity-neutralization (VADI §57):** Agents carrying a BI cluster targeting the same receptor get protection scaled by `immunity_factor * immunity_binding_affinity`. For wild-type clusters (`immunity_binding_affinity = 1.0`) this yields the full 1000× protection (`immunity_factor = 0.001`). Super-killer immunity-escape variants have reduced `immunity_binding_affinity` (typically 0.01–0.3), so the effective protection drops to 3–100× — making the agent partially vulnerable to the novel toxin even though it carries the ancestral immunity gene. The best (lowest `eff`) among all matching BI clusters is used.

---

## 4. fix_conjugation — Horizontal Gene Transfer

**Biological basis:** F-pili mediated conjugation transfers BI clusters between cells in physical contact. The mating-pair stabilization (MPS) is sensitive to local shear rate.

**Transfer probability:**
```
P_transfer = P_base * exp(-gamma / gamma_crit) * dt
```
where `gamma` is the local shear rate (from advection field) and `gamma_crit` is the critical shear for pilus retraction.

**Requirements:**
- Donor must have `has_conjugative_plasmid = true`
- Cells must be within contact distance (`contact_radius`, default 2 um)
- Recipient must not already carry the same toxin_id

**F-pili length heterogeneity (VADI §55):** When `pili_heterogeneity = true`, each conjugation attempt samples its effective contact radius from `uniform(pili_length_min, pili_length_max)` (default 1–4 μm). This replaces the fixed `contact_radius` for that attempt. In vivo F-pili exhibit significant length heterogeneity; shorter pili restrict transfer to immediate neighbors while longer pili can bridge larger gaps. The stochastic sampling means each donor–recipient pair evaluation uses an independently drawn pilus length per timestep.

**Effect:** Recipient gains a copy of one randomly-selected BI cluster from the donor. This enables the spatial spreading of bacteriocin-immunity cassettes through the population.

---

## 5. fix_mutation — Stochastic Genome Evolution

**Six mutation classes:**

### a) BI locus duplication (rate: 10^-5 per division)
Tandem duplication of an existing BI cluster. Expands the agent's toxin range. Capped at `max_bi_loci` (default 8).

### b) BI locus recombination (rate: 5 × 10^-6 per division)
Swaps immunity proteins between two BI clusters. Creates novel toxin/immunity combinations that may escape existing immunity.

### c) Receptor downregulation (rate: 10^-7 per division)
Reduces expression of a randomly-selected receptor by `receptor_reduction` (default 0.1). When expression drops below 0.2, agent transitions to RESISTANT phenotype. Provides toxin resistance at the cost of nutrient uptake.

### d) Partial resistance — extracellular loop missense (rate: 5 × 10^-7 per division)
Missense mutations in receptor extracellular loops can abrogate colicin binding while retaining sufficient affinity for native ligands (VADI §78, §101). Unlike full receptor downregulation:
- **Receptor expression is unchanged** — no MetE pathway penalty, no loss of nutrient uptake capacity
- **Toxin Kd increased 10–100×** — `toxin_affinity` set to [0.01, 0.1], effectively raising the concentration required for colicin binding
- **Ligand Kd within 2× of wild-type** — `ligand_affinity` set to [0.5, 1.0], preserving nutrient transport

This allows strains to potentially bypass the Combinatorial Washout Trap by resisting toxins without incurring the metabolic double-bind penalty. The `toxin_occupancy()` function in fix_receptor scales both Kd values:
```
eff_Kd_tox    = Kd_tox    / toxin_affinity     (larger → less toxin binding)
eff_Kd_ligand = Kd_ligand / ligand_affinity     (slightly larger → slightly less nutrient binding)
```

### e) Super-killer emergence (rate: 10^-8 per division)
Generates a novel toxin variant from an existing BI cluster. The new toxin has:
- Unique toxin_id
- Gaussian-perturbed pI (std = 0.5)
- 30% chance of switching target receptor
- **Immunity-escape mutation** (probability `immunity_escape_prob`, default 50%): reduces `immunity_binding_affinity` to a value drawn uniformly from [`escape_affinity_lo`, `escape_affinity_hi`] (default [0.01, 0.3]). This models the VADI §57 affinity-neutralization matrix where structural changes in the novel toxin reduce the binding affinity of the cognate immunity protein, creating variants that partially evade immunity.

The immunity-escape mechanism operates through the agent's own BI cluster: `immunity_binding_affinity` on a given cluster determines how well that cluster's immunity protein cross-neutralizes incoming toxins targeting the same receptor. When a neighbouring cell produces an escape-variant toxin, the target agent's existing immunity provides only partial protection (effective factor = `immunity_factor × immunity_binding_affinity`).

### f) Compensatory mutation (rate: 10^-6 per division)
Chromosomal mutations that ameliorate plasmid metabolic cost. Reduces per-locus cost by `compensatory_reduction` (default 0.005), capped so that at most 75% of the original cost can be eliminated. This models the well-documented phenomenon where bacteria rapidly evolve to reduce the fitness cost of newly acquired plasmids (VADI Section 79).

---

---

## 7. fix_motility — Active Cell Swimming (Spec 3 / Spec 10v2)

**Biological basis:** E. coli in colonic mucus swims at ~7.76 µm/s using a mucus-adapted run-and-reverse pattern (not classic run-and-tumble). Lazova et al. 2016 showed Aer-mediated aerotaxis is the primary directional cue for colonization through mucus-like hydrogels with vertical O₂ gradients. Carbon chemotaxis (Tar/Tsr), energy taxis, surface sensing, and mucin viscosity further shape spatial distribution.

**Model:** `FixMotility::pre_step()` advances every run/stop transition that
occurs within the biological timestep and accumulates the corresponding swim
displacement. `module_physics()` applies that displacement after advection:
```
update_chemotaxis (Weber–Fechner run-timer bias)
while remaining_dt > 0:
    event_dt = min(remaining_dt, current_run_or_stop_timer)
    if running:
        speed = effective_swim_speed(agent)   # energy × surface × mucin
        displacement += swim_direction * speed * event_dt
    advance timer and transition state when the event ends
x += displacement
```

**Directional taxis (run-length modulation, additive):**
- **Aerotaxis** (`aerotaxis_enabled`, default on): Weber–Fechner on the O₂ field; sensitivity `aerotaxis_sensitivity` (default 4.0). Primary cue when oxygen is present.
- **Carbon chemotaxis** (`chemotaxis_enabled`): Weber–Fechner on carbon; sensitivity `chi_carbon` (default 2.0). Floor `chemotaxis_threshold` avoids division by zero.
- **AI-2 chemotaxis** (`quorum_sensing.ai2_chemotaxis`, Spec 11): Weber–Fechner on the AI-2 field; sensitivity `chi_ai2` (default 3.0). Requires `quorum_sensing.enabled`.

**Speed modulation (multiplicative):**
- **Energy taxis**: `speed *= floor + (1-floor) * mu_realized/mu_max`
- **Surface sensing** (opt-in): linear ramp from `surface_sensing_floor` at the epithelium to full speed beyond `surface_sensing_depth`
- **Mucin drag** (opt-in): `speed *= ref / (ref + [mucin])`

Cluster suppression reduces reorientation rate when neighbor density exceeds threshold.

---

## 7b. fix_quorum_sensing — AI-2 Quorum Sensing (Spec 11)

**Biological basis:** *E. coli* produces autoinducer-2 (AI-2) via LuxS and imports it via the Lsr ABC transporter. LsrB–Tsr chemotaxis toward AI-2 promotes conspecific clustering (Hegde et al. 2011; Laganenka et al. 2025), raising local toxin, CDI, and conjugation encounter rates.

**Model:** When `quorum_sensing.enabled`, the `ai2` chemical species is registered (diffusing small molecule, no z-gradient, zero IC/BC). `FixQuorumSensing::compute()` deposits per-agent reactions:

```
production = ai2_basal_rate + ai2_growth_coupled * max(mu_realized, 0)
import     = lsr_vmax * [AI-2] / (lsr_km + [AI-2])
reac      += production / cell_vol − import / cell_vol − decay_rate * [AI-2]
```

Chemotaxis is handled by `fix_motility` (Weber–Fechner on AI-2). Characteristic diffusion length √(D/decay) ≈ 70 µm sets the attraction range.

---

## 8. fix_cdi — Contact-Dependent Inhibition (Spec 3)

**Biological basis:** CDI (CdiA/CdiB) systems deliver toxic effectors upon direct cell-cell contact. Strain-specific immunity protects cognate partners. Corpse barriers at colony interfaces self-limit killing.

**Logic:** For each CDI+ agent, neighbors within `contact_radius` without matching immunity are killed with `P = 1 - exp(-kill_rate * dt)`. CDI kills record `death_time`; corpses persist as mechanical obstacles for `corpse_persistence` seconds and block line-of-sight CDI to cells behind them.

---

## 9. Fur-Regulated Receptor Expression (Spec 3)

**Biological basis:** The Ferric Uptake Regulator (Fur) represses iron-uptake genes when iron is abundant and derepresses them under starvation — the "Vulnerability Paradox" where iron-depleted zones have higher colicin susceptibility.

**Implementation:** In `fix_metabolism`, before graded iron Monod:
```
fur_factor = 1 + upregulation_max * Km / (Km + [iron])
receptor_expr[r] = min(receptor_expr_base[r] * fur_factor, receptor_max)
```
for iron-related receptors only. Mutations modify `receptor_expr_base`.

---

## 10. fix_mechanics — Soft-Sphere Mechanical Repulsion

**Biological basis:** Bacterial cells are deformable spheres with an elastic modulus of ~0.1–1 MPa (measured by AFM). When two cells overlap due to growth or flow, they exert a repulsive contact force following Hertzian contact mechanics.

**Hertzian contact model:**
```
overlap = r_i + r_j - d       (positive when cells overlap)
F = k * overlap^(3/2)         (Hertzian, non-linear)
```

This is physically more accurate than a linear spring (`F ∝ overlap`) because:
1. Contact area grows with indentation → force rises faster than linearly
2. Matches AFM force-displacement curves on bacterial cells
3. Prevents excessive interpenetration at large overlaps

**Overdamped force application:** Mucus mechanics use Stokes mobility rather
than inertial displacement. For each agent in a pair:
```
mobility_i = 1 / (6 * pi * mu * r_i)
delta_x_i = F * dt * mobility_i
```
Here `mu` is the existing `vbf_viscosity` parameter, so the displacement
relationship is `delta_x = F * dt / (6 * pi * mu * r)`. Pair movement remains
equal and opposite, with each agent's share partitioned by its mobility rather
than by mass. The total mechanics displacement accumulated by one agent in a
biological step is capped at `0.1 * r` to prevent numerical explosions while
retaining the clamp count in HDF5 output.

**EPS-mediated adhesion (optional):**

When `adhesion_enabled = true`, cells within `adhesion_range` of contact (but not overlapping) experience a short-range attractive force:
```
gap = d - (r_i + r_j)
F_adhesion = adhesion_strength * (1 - gap / adhesion_range)
```
This models extracellular polymeric substance (EPS) bridging, relevant for:
- Biofilm microcolony formation
- Kin-recognition clustering
- Protection from advective washout

The adhesion force decays linearly to zero at `adhesion_range`, preventing long-range artifacts.

**Calibration:** The default `hertz_k = 1e-6 N/m^1.5` corresponds to the effective spring constant for two 0.5 µm radius *E. coli* cells with Young's modulus ~0.5 MPa and Poisson ratio ~0.5.

---

## Chemical Transport — Implicit Nutrients and QSSA Toxins

### Nutrient and small-molecule diffusion

An explicit 3-D stencil is unusable at the biological timestep: for O₂ at `D = 2.1e-9 m²/s`, `dt = 60 s`, and `dx = 5 µm`, the diffusion number is `D·dt/dx² = 5040`, versus the explicit stability limit `1/6`. GutIBM therefore uses backward-Euler directional splitting for enabled nutrient fields.

Each directional pass solves a tridiagonal system in O(cells): x and y are periodic and the luminal z face has zero flux. The epithelial z=0 face defaults to a fixed concentration (`dirichlet`); configured species may instead use `robin`, `J = k(C_epi - c0)`, or `flux`, a fixed delivery rate in mol/m²/s. Delivery modes solve all `nz` cells with the epithelial cell as an unknown and record the realized post-solve exchange in the nutrient ledger. The method is L-stable; non-delivery reactions are clamped nonnegative after the solve, while prescribed delivery uses the local-then-global positivity rationing described above rather than clipping a draw. For species with a configured exponential z-gradient, diffusion acts on departures from that prescribed background profile rather than erasing it; the gradient is rejected with Robin or flux because its pinned reference assumes a Dirichlet boundary. The chemistry order is rank-local agent reactions → MPI sum → global VBF coupling → concentration update → implicit diffusion → boundary enforcement. On GPU-active steps, host-written reactions are uploaded before any device reaction kernel runs, so the device reaction buffer is the single accumulated source for the step; resulting concentrations are synchronized back to the host after device integration. Delivery uptake remains explicitly CPU-only because the GPU diffusion path has no equivalent prescribed-sink retry/rationing loop; configuration finalization rejects that unsupported combination rather than silently diverging.

### Bacteriocin QSSA solver

Bacteriocins remain on the analytical QSSA path because they are point-source bursts or continuous producer sources with receptor-specific Green's-function fields. The default `chemistry.toxin_lumping = "per_receptor"` keeps separate BtuB, FepA, CirA, and FhuA fields. The scientific `lumped` variant instead superposes every toxin source into one `bacteriocin_lumped` field and makes all receptors read it, while retaining receptor-specific binding, toxin affinity, and immunity in the kill calculation. Its explicit modelling costs are that a colicin to which an agent is immune can still contribute to exposure through another receptor, and that the nuclease SOS/ROS cross-induction path sees total toxin burden rather than only nuclease-colicin/BtuB toxin.

The `chemistry.species_subset` modelling position controls which chemical
species and mechanisms exist. `full` preserves this complete model.
`nutrient_only` removes bacteriocin fields and disables bacteriocin release,
QSSA toxin chemistry, and receptor killing to isolate nutrient competition.
`carbon_only` retains carbon alone and also disables oxygen, acetate,
ethanolamine, mucin, siderophore, ferrichrome, quorum sensing, and motility
taxis terms that read removed fields, Fur regulation, iron/B12/eut uptake
terms, the VBF iron sink, and dynamic mucin; it asks whether carbon
competition alone is sufficient for retention and clustering. Zeroing the VBF
iron sink and disabling dynamic mucin are intentional parts of this modelling
position, not hidden parameter overrides. Initialization validates required
species from mechanism enablement. An enabled mechanism with a missing species
throws `ConfigError` naming the mechanism, species, and enabling configuration
key instead of silently turning that term off. In lumped mode, the nuclease
SOS/ROS cross-induction reader also sees total toxin burden, not a nuclease-only
field; this is an explicit second approximation cost.

**Method:** At each biological timestep:
1. Collect all active toxin sources (SOS-lysed cells as bursts, microcin producers as steady sources)
2. In MPI runs, allgather the compact source records in rank-major order so every rank solves from the same global source list. This exchanges source records, not grid-sized field data; ranks with no local sources still participate.
3. Evaluate the analytical Green's function at each grid cell:
   ```
   C(r) = Q / (4 * pi * D_eff * r)    (3D point source, steady state)
   ```
4. Apply Method of Images for bounded mucus domain (z = 0 epithelium, z = h lumen)
5. Superpose contributions from all sources within cutoff radius

**Taylor-Aris dispersion:** The model does not include this enhancement. At
shipped parameters its effect is negligible for the species where the
long-time Taylor-Aris limit is valid, while the regime where it is large does
not satisfy that limit. Adding it honestly would require an anisotropic
transport kernel, plus re-derivation of the image series and Robin correction
table. See the decision and measurements in
`docs/EXTERNAL_AUDIT_2026-08.md`.

### FMM Far-Field Acceleration

When `QSSAConfig::use_fmm = true`, the QSSA solver uses a kernel-independent Fast Multipole Method to accelerate the Green's function superposition:

1. **Build**: Construct an octree over the globally exchanged active toxin sources with Cartesian multipole moments up to `fmm_expansion_order` (default 2).
2. **Upward pass (P2M + M2M)**: Accumulate particle charges into leaf moments, then shift and combine into parent nodes.
3. **M2L + L2L**: Translate well-separated remote multipoles into local expansions and propagate downward for O(N+M) grid evaluation.
4. **Near-field** (distance ≤ `toxin_cutoff`): Evaluate exact advection-diffusion Green's function per-source, preserving accuracy for nearby lethal core/halo interactions.
5. **Far-field** (distance > `toxin_cutoff`): Deposit `total_fmm − near_exact` per grid cell using the precomputed local expansions.

The opening angle `theta` controls the well-separated criterion. Expansion order `p` gives error ~ theta^p for smooth fields. Order 1 recovers monopole Barnes-Hut; orders 2–3 add dipole/quadrupole/octupole terms via source-position Taylor coefficients of the advection-diffusion kernel (not assuming a 1/r Laplace kernel).

---

## Advection — Dual-Vector Mucus Flow

**Two flow components:**
1. **Radial** (z-axis, epithelium → lumen): Driven by mucus shedding, turnover ~1.5 h
2. **Distal** (x-axis, proximal → distal): Peristaltic transit, ~12 h

**Velocity profile:** Power-law with exponent alpha (default 1.5):
```
v(z) = v_max * (z/h)^alpha
```
Near-zero at epithelium (z = 0), maximum at lumen (z = h). This creates a spatial refuge near the epithelium where flow is minimal, allowing resident colonies to persist.

**Washout criterion:** When `mu_realized < gamma_flow = v_radial(z) / h`, the agent cannot grow fast enough to resist radial flow and is expelled.

### Peristaltic Mixing (VADI §77)

Real colonic flow is not steady-state — periodic slow-wave contractions modulate the local flow velocity at ~15–30 s intervals. GutIBM models this via an oscillatory perturbation:

```
v_effective(pos, t) = v_base(pos) * peristaltic_factor(pos, t)
peristaltic_factor  = 1 + A * sin(2π t/T − 2π x/λ)
```

where:
- **A** = amplitude (default 0.5, i.e. ±50% modulation)
- **T** = period (default 20 s)
- **λ** = wavelength (0 = spatially uniform oscillation; >0 = propagating contractile wave along x-axis)

**Effect on agents:** Bacteria near the lumen experience periodically enhanced and reduced advective drag. This creates transient "washout surges" that can dislodge cells with marginal growth rates (`mu_realized ≈ gamma_flow`), and quiescent intervals that allow brief recovery. The net effect increases the selective pressure for fast-growing cells near the luminal surface while leaving epithelium-proximal residents relatively unaffected (due to the power-law velocity profile).

**Propagating wave mode:** When `wavelength > 0`, the oscillation has a spatial phase offset along the x-axis (distal direction), mimicking a contractile wave propagating aborally. Cells at different x-positions experience the peak flow at different times, which enhances longitudinal mixing and dispersal.

---

## MPI Domain Decomposition

The simulation domain is partitioned across MPI ranks using 1D slab decomposition along the x-axis (distal flow direction). This enables multi-node scaling for large agent populations.

### Slab Decomposition
The domain is divided into cell-aligned slabs along the x-axis. Each rank owns agents within its contiguous slab bounds `[local_lo_x, local_hi_x)`, with boundaries at chemical-grid cell faces and widths balanced to within one cell.

### Ghost Layers
Before the biology module, each rank exchanges **ghost agents** — copies of agents within `ghost_width` of slab boundaries — with its neighbors. This ensures correct neighbor queries (spatial hash lookups) for cross-boundary interactions (conjugation, mechanical repulsion, receptor binding).

Ghost agents are read-only and discarded before the physics module to avoid double-counting position updates.

### Agent Migration
After the physics module (advection + mechanics), agents that have moved past their slab boundary are serialized and sent to the appropriate neighbor rank via `MPI_Sendrecv`. The full agent state — position, velocity, metabolism, genome, BI clusters — is transferred.

### Global Statistics
`MPI_Allreduce` aggregates per-rank counts and growth rate sums to produce global agent count and mean growth rate. These are used for output and lineage tracking.
HDF5 summary interval and cumulative event counters are likewise globally reduced once per summary, so they describe the same global population as `n_total`.
Death-channel counters include `mortality_lysis` for actual SOS/phage lysis deaths;
induction counters are not death counts.
The summary also contains the instantaneous `bacteriostatic_live_agents` and
`washout_trapped_live_agents` stocks; these are not event counters and are not
part of population closure.

### HDF5 Parallel I/O
When `hdf5.parallel = true`, the file is opened with `H5Pset_fapl_mpio` and agent data is written using collective hyperslab operations — each rank writes its local agents at a computed offset. In slab chemistry mode, all ranks pack owned x-cells for grid output and rank 0 assembles the global datasets; halo planes are never written. Checkpoint/restart restores owned global cells and refreshes concentration halos. The GPU mirror remains unsupported with slab chemistry.

---

## Crypt Refugia — Zero-Flow Zones (VADI §80, §98–99)

**Biological basis:** The model’s original assumption that advective flow is inescapable throughout the domain is unrealistic. Intestinal crypts and the firmly adherent inner mucus layer have γ_flow ≈ 0, providing physical refugia where the Washout Trap threshold (μ_realized < γ_flow) can be bypassed entirely.

**Implementation:**

When `crypts_enabled = true`, a crypt zone is defined as `z < lo_z + crypt_depth`. Within this zone:

1. **Zero flow:** `velocity()`, `radial_velocity()`, `distal_velocity()`, and `washout_rate()` all return 0. Agents are not advected.
2. **Washout bypass:** `check_washout()` skips agents with `in_crypt = true`, regardless of their μ_realized.
3. **Stochastic migration:** Each timestep, `crypt_migration(dt)` processes:
   - **Exit:** Crypt agents exit with probability `1 - exp(-crypt_exit_rate * dt)`, placed at `z = crypt_depth + ε`.
   - **Entry:** Flow-zone agents near the crypt boundary enter with probability `1 - exp(-crypt_entry_rate * dt)`, placed at random z within the crypt zone. Entry is blocked when the crypt population reaches `crypt_carrying_capacity`.
4. **Agent flag:** Each agent carries an `in_crypt` boolean, set during initialization for agents spawned inside the crypt zone and updated dynamically by migration.

**Derivative issue:** The current implementation uses a simple z-threshold model. More complex crypt geometry (discrete invaginations, variable crypt density) may be warranted for spatially heterogeneous mucosa.

---
## Viscoelastic Background Field (VBF)

The 99% obligate anaerobic microbiota is modeled as a continuum rather than discrete agents:
- **Physical drag**: Stokes-like force opposing agent velocity
- **Nutrient sink**: Background consumption at volumetric rate
- **Mucin liberation**: Monosaccharide release from mucin glycoproteins (carbon source)
- **Carrying capacity**: Local density limit for the simulation domain

### Implicit carbon competition sink

The VBF carbon sink follows Monod kinetics,
`vmax * c / (Km + c)`, integrated implicitly over each biological timestep.
The backward-Euler solve reports the realized concentration removal rather than
the unconstrained instantaneous demand, so the sink cannot remove more carbon
than the local field contains. That realized amount is written both into the
carbon reaction rate and `vbf_sink` accounting. Any residual reaction
positivity clip is recorded separately in the nutrient-flux summary. Agent-side
uptake remains an independent pathway and can still overdraw a cell.

With `vbf.agent_carbon_coupling` nonzero, the VBF `vmax` in each voxel gains
`agent_carbon_coupling * n_owned_live_agents / V_cell`. The default zero keeps
the historical sink unchanged; ghost and dead agents are excluded. The
nutrient-flux summary reports the derived blocking fraction. Its numerator is
realized agent carbon removal, `agent_uptake + maintenance`; unpaid
maintenance shortfall is excluded. The denominator is that realized agent
removal plus realized VBF sink removal.

---

## z-Dependent Nutrient Gradient

**Biological basis:** Mucin-derived monosaccharides (the primary carbon source for *Enterobacteriaceae* in the gut) are liberated by goblet cells at the epithelial surface (z=0). The concentration of free monosaccharides is highest near the epithelium and decays exponentially into the lumen.

**Implementation:**

### Chemical field initialization
When `z_gradient_enabled` is set for a chemical species, the initial concentration follows:
```
C(z) = C_max * exp(-z_rel / lambda_mucin)
```
where `z_rel` is the distance from the epithelium and `lambda_mucin` is the characteristic decay length (~25 μm by default). The default Dirichlet boundary at z=0 maintains `boundary_conc` as the peak value. Carbon can opt into finite-rate delivery with `carbon.epithelial_boundary` (or its underscore alias): Robin uses `carbon.epithelial_transfer_coeff` and the existing `boundary_conc` as `C_epi`; flux uses `carbon.epithelial_flux`. The post-solve exchange is recorded in mol as `beta(C_epi - c0_after)V` for Robin or `J·A·dt` for flux. Oxygen exposes the same boundary keys; in Robin mode its `epithelial_conc` is the tissue-side reservoir rather than the imposed mucus-surface concentration, and its imposed oxygen profile must be disabled because the surface value is predicted by delivery and consumption.


Agent uptake can be limited to what diffusion can deliver. With
`uptake_limit=sherwood`, an agent's per-step uptake is capped at
`N_max = 4*pi*D_eff*r_agent*C_local*dt` (Sherwood number 2 for an isolated
sphere at rest), with `D_eff = D_free / retardation` taken from the same
species record the field solver uses. Delivery mode uses that same ceiling as
a prescribed mass inside the implicit nutrient solve. The realized fraction is
the Liebig minimum over the capped species (carbon and iron; the corrinoid
pool is not depleted, Spec 6 §3), clamped to `[0, 1]`, and it scales the
realized uptake, the biomass increment, and `mu_realized`, so division,
bacteriostasis classification, and the washout/VADI `mu` versus `gamma_flow`
comparison all see the funded rate. `uptake_limit=voxel` caps at the voxel
content `C_local*V_cell` and exists only to measure the grid artifact by
divergence from `sherwood`; it is not a biological model. `uptake_limit=none`
is the default and leaves growth unfunded as before. The Sherwood relation is
defined using a far-field concentration. For compatibility,
`delivery_far_field_radius=0` evaluates it at the agent's own voxel, which
double-counts local depletion and therefore introduces a resolution-dependent
ceiling. In campaign measurements the own-voxel concentration was
approximately 250-fold below the far field at 2 µm grid spacing, but only
approximately 1.4-fold below it at 6 µm. A positive
`delivery_far_field_radius` instead uses a volume-weighted mean over cell
centres inside the radius, including the agent cell with its actual volume
weight. Delivery-mode runs now use a 10 µm support by default. At that
radius, single-cell physical depletion is approximately 0.008% of the
far-field concentration, so the support regularizes the unresolved sink
rather than doing biological work. The same spherical support is used for prescribed delivery
deposition: the funded mass is divided uniformly across included cell
centres, with the complete mass renormalized over cells remaining inside the
domain when the support reaches a nonperiodic z boundary. Periodic x/y
coordinates wrap as in the concentration read. The grid cannot resolve the
sub-micron depletion boundary layer around a cell, so this regularizes the
point sink over a physical, grid-independent support; the approximately
500 µm diffusion length over 60 s makes a support of approximately 10 µm
well inside the diffusive reach. Radius zero remains intentionally compatible
but grid-dependent: it reads the own voxel and deposits the full prescribed
mass in that single voxel. Positive radii are refused in slab chemistry
because distributed deposition into halo cells cannot preserve ownership and
mass accounting safely; slab configurations must explicitly set
`metabolism.delivery_far_field_radius=0.0` to opt into the grid-dependent
single-voxel model. The parser and simulation initialization both refuse a
positive radius in slab chemistry. This changes shipped behaviour: delivery
mode now deposits over a 10 µm support by default. Previously reported
supply-limited carrying capacities, including the carbon ladder in #314, were
measured under single-voxel deposition and therefore move under this default.
In the settled 80-founder, ROS-off, six-hour population campaign, support-on
final populations were 201/207/204 at 2/4/6 µm with funded oxygen fraction
0.988 at all three resolutions; support-off final populations were 76/123/171
with funded fractions 0.620/0.901/0.948.
### VBF mucin liberation coupling
When `mucin_z_gradient_enabled`, the monosaccharide release rate applied by the VBF also varies with z:
```
rate(z) = mucin_liberation * exp(-z_rel / mucin_z_gradient_lambda)
```
This couples the ongoing carbon source term to the same spatial gradient, maintaining the profile during simulation.

**Spec 1 extensions:** When `mucin.enabled`, liberation follows `k_liberation * vbf_density * [mucin] / (Km + [mucin])` with goblet secretion at z=0. When `acetate.enabled`, VBF adds z-weighted acetate production and epithelial uptake. When `oxygen.enabled`, VBF applies a constant O₂ sink.

### Biological consequences
- Agents near the epithelium (z ~ 0) have access to more carbon, leading to faster Monod growth.
- Agents further from the epithelium face carbon limitation, reducing their competitive fitness.
- Combined with the advective washout trap (flow increases with z), this creates a strong spatial advantage for epithelium-attached colonies: high nutrients, low flow.

---

## Spatial Validation — Exclusion-Radius Clustering (VADI §75)

**Background:** The original plan to validate strain-specific spatial patterns via HiPR-FISH targeting immunity mRNA was abandoned because immunity transcripts exist in single-digit copy numbers per cell — below the detection threshold of standard HiPR-FISH probes.

**Replacement approach:** Exclusion-radius and NND clustering metrics that can be validated with DNA-FISH phylogroup probes or HCR-FISH amplification, both of which target multicopy sequences.

**Metrics computed by `validate_spatial_signatures()`:**

| Metric | What it measures | Empirical target |
|--------|-----------------|------------------|
| `monochromatic_score` | Same-type neighbor fraction | > 0.7 |
| `comet_tail_ratio` | Downstream/upstream toxin concentration | > 1.5 |
| `mean_exclusion_radius` | Mean distance to nearest competing-type boundary | Phylogroup-dependent |
| `hopkins_statistic` | Global spatial clustering (Hopkins H) | > 0.7 |
| `nnd_mean` | Grand mean NND between competing clones | Phylogroup-dependent |
| `comet_tail_asymmetry` | Concentration-weighted downstream elongation | > 1.0 |

The exclusion radius captures the characteristic "no-go zone" each bacteriocin-producing strain creates around itself; NND between competing clones quantifies how far apart rival phylogroups settle; the Hopkins statistic confirms non-random spatial arrangement.
