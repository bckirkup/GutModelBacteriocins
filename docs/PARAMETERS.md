# Parameter Reference

All configurable parameters for the GutIBM simulation, with defaults and biological justification.

For **batch runner** manifests (multi-run sweeps), see [BATCH_RUNNER.md](BATCH_RUNNER.md).

---

## Time Control

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `total_time` | 86400 | s | Simulation duration (24 h) |
| `bio_dt` | 60 | s | Biological timestep |
| `output_interval` | 3600 | s | Console progress + in-memory lineage snapshot interval (not HDF5) |
| `seed` | 42 | — | Random number generator seed |
| `dysbiosis_threshold` | 1e8 | cells/mL | Spec 5 §4 bloom boundary: with the sampling controls below, halt only after a persistent above-boundary rise whose net window trend is positive and whose increment means do not decelerate. See [OPERATING_ENVELOPE.md](OPERATING_ENVELOPE.md). `0` disables the check |
| `dysbiosis_sampling_interval` | 300 | s | Simulated-time interval between global-density samples for the bloom guardrail |
| `dysbiosis_sample_count` | 7 | samples | Consecutive samples required: 7 samples span 30 minutes from first to last; all must be above the boundary, the net rise must be positive, and the mean increment in the second half must be at least the first half (both positive) |

**Guidance:** `bio_dt` should be ≤ 60 s for accurate growth dynamics. Larger values speed up simulation but may miss fast-timescale events (SOS induction, toxin killing).

### Adaptive Timestep

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `adaptive_dt_enabled` | false | — | Enable CFL-like adaptive timestep |
| `dt_min` | 1 | s | Minimum allowed timestep |
| `dt_max` | 300 | s | Maximum allowed timestep |
| `dt_safety` | 0.8 | — | Safety factor applied after constraint computation |
| `dt_growth_limit` | 0.1 | — | Maximum mu × dt product allowed |

When `adaptive_dt_enabled = true`, the simulation replaces the fixed `bio_dt` loop with a `while (time < total_time)` loop that recomputes `dt` each step based on three constraints:

1. **Growth rate**: `dt ≤ dt_growth_limit / max(|mu_realized|)` — prevents large fractional biomass changes per step.
2. **SOS cascade**: `dt ≤ 10 s` when >5 SOS-induced agents; `dt ≤ 2 s` when >20 — resolves lysis burst dynamics.
3. **Agent density**: `dt ≤ 10 s` when density exceeds 10^15 /m^3 — avoids numerical overlap in crowded regions.

After constraint evaluation, `dt` is multiplied by `dt_safety` and clamped to `[dt_min, dt_max]`.

When disabled (`adaptive_dt_enabled = false`), the fixed `bio_dt` is used as before.

## Initial population

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `initial_population.placement` | `legacy` | — | `legacy` preserves the historical `[domain.lo.z, 0.5 * domain.hi.z]` founder band; `z_slab` uses explicit bounds; `anatomic` uses anatomy-derived truncated-exponential depth |
| `initial_population.z_min`, `initial_population.z_max` | — | m | Inclusive lower and exclusive upper bounds for the founder `z_slab`; both must lie inside the domain and `z_min < z_max` |
| `initial_population.anatomic_exclusion_floor` | 20e-6 | m | Hard minimum founder depth for `anatomic`; no founder is placed below the epithelial exclusion zone |
| `initial_population.anatomic_exponential_scale` | 40e-6 | m | Exponential depth scale above the exclusion floor for `anatomic`; must be positive |
| `initial_population.anatomic_outer_extent` | 150e-6 | m | Exclusive outer truncation for `anatomic`; draws at or above this depth are rejected and resampled, and the extent must exceed the exclusion floor |

The placement policy is scenario-wide and applies to every entry in
`initial_strains`. With the default `legacy` policy, founder placement remains
derived from the domain height for backward compatibility. Use `z_slab` when
changing `domain_z` without moving founders, for example to keep them in the
mucus layer while extending the lumen. Use `anatomic` to apply the hard 20 µm
epithelial exclusion, a 40 µm exponential depth distribution, and 150 µm
outer truncation. The latter two values describe anatomy reported by
Swidsinski FISH observations and Duncan/Mondragón-Palomino outer-mucus
enrichment; they are not fitted to simulation trajectories. Anatomic founders
never enter crypts because healthy FISH observations constrain crypt absence.
x and y remain uniform. Immigration has separate placement handling and is not
changed by this policy.

## Immigration

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `immigration.enabled` | false | — | Enable propagule injection |
| `immigration.count` | 1 | cells/event | Cells injected per event |
| `immigration.strain_index` | 0 | index | Entry in `initial_strains` used to construct immigrants |
| `immigration.placement` | `uniform` | — | `uniform`, `at_distance`, or `z_slab` |
| `immigration.distance` | 0 | m | Target nearest-biomass/centroid distance |
| `immigration.distance_tolerance` | 0 | m | Maximum placement error |
| `immigration.distance_reference` | `nearest_agent` | — | Nearest live cell or global centroid |
| `immigration.z_min`, `immigration.z_max` | 0 | m | Bounds for `z_slab` placement |
| `immigration.schedule` | `pulse` | — | One pulse or Poisson `continuous` schedule |
| `immigration.step` | 0 | relative step | Pulse step relative to run start |
| `immigration.rate` | 0 | events/s | Continuous event rate |

`at_distance` proposes positions on shells around sampled live-biomass
anchors, then verifies each proposal using the true global minimum-image
distance to the nearest live, pre-existing agent. The centroid convention
anchors on the global centroid and is intended only for a single compact
cluster; it is not min-image-corrected across periodic wraps. If no proposal
meets `distance_tolerance` after the retry, the cell is skipped with a rank-0
warning rather than placed at the wrong distance.

---

## Domain

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `domain.lo` | {0,0,0} | m | Lower corner |
| `domain.hi` | {1e-3, 1e-3, 100e-6} | m | Upper corner (1mm × 1mm × 100um) |
| `domain.grid_dx` | 2e-6 | m | Fine chemistry grid spacing |
| `domain.chemistry_stride.x/y/z` | 1 | — | Integer chemistry-cell coarsening per axis; x/y probe lateral microgradient sensitivity while z remains fine for epithelial O2/mucin stratification |
| `domain.grid_halo_width` | 1 | cells | Slab-chemistry halo width; must cover `ceil(domain.ghost_width / domain.grid_dx)` and is not a tuning knob |
| `domain.hash_cell_size` | 10e-6 | m | Spatial hash bucket size |
| `domain.periodic` | {true, true, false} | — | Periodicity per axis |
| `domain.mpi_decomp_axis` | 0 | — | Axis for 1D slab decomposition (0=x) |
| `domain.ghost_width` | 10e-6 | m | Ghost layer thickness for cross-rank neighbor queries |

**Biological context:** The domain represents a patch of colonic mucus layer. x,y are periodic (infinite mucosa plane). z spans from epithelium (z=0) to luminal surface (z=h).

**MPI decomposition:** The domain is partitioned into cell-aligned slabs along `mpi_decomp_axis` (default: x, the distal flow direction). Each rank owns a contiguous half-open range of grid cells; physical slab widths can differ by at most one cell. `ghost_width` should be ≥ `hash_cell_size` to ensure correct neighbor queries across slab boundaries.

For `chemistry.decomposition=slab`, `domain.grid_halo_width` is the number of chemistry cells retained around each owned slab. It must be large enough to cover the physical `domain.ghost_width` at the chemistry spacing; configure it from those two values rather than treating it as a numerical tuning parameter.

---

## Advection

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `advection.radial_turnover` | 5400 | s | Mucus radial turnover (1.5 h) |
| `advection.mucus_thickness` | 100e-6 | m | Mucus layer depth |
| `washout.trap` | `emergent` | — | `emergent` removes only at the luminal boundary; `imposed` retains the explicit `mu_realized < washout_rate(z)` comparison variant |
| `advection.distal_transit_time` | 43200 | s | Peristaltic transit (12 h) |
| `advection.distal_length` | 1e-3 | m | Domain length for transit calc |
| `advection.profile_alpha` | 1.5 | — | Flow profile exponent |
| `advection.taylor_aris_enabled` | true | — | Enable Taylor-Aris dispersion |

**Velocity profile:** `v(z) = v_max * (z/h)^alpha`

- `alpha = 1.0`: linear shear (Couette-like)
- `alpha = 1.5`: intermediate (default)
- `alpha = 2.0`: parabolic (Poiseuille-like)

**Taylor-Aris:** Enhances effective longitudinal diffusion via:
`D_eff = D_mol + U(z)^2 * h^2 / (210 * D_mol)`

This captures shear-enhanced spreading of toxins in the mucus flow. The
constant `210` is the classical result for fully-developed parabolic
(Poiseuille) flow (`profile_alpha = 2`); for other profile exponents the
prefactor differs, so with the default `profile_alpha = 1.5` this is an
order-of-magnitude approximation rather than an exact coefficient.

### Peristaltic Mixing (VADI §77)

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `advection.peristaltic_enabled` | false | — | Enable oscillatory flow modulation |
| `advection.peristaltic_period` | 20 | s | Slow-wave period (colonic: 15–30 s) |
| `advection.peristaltic_amplitude` | 0.5 | — | ±50% modulation of flow velocity |
| `advection.peristaltic_wavelength` | 0.0 | m | 0 = uniform, >0 = propagating contractile wave |

**Modulation:** When enabled, both radial and distal velocities are multiplied by:
```
peristaltic_factor = 1 + amplitude * sin(2π t/period − 2π x/wavelength)
```
With `wavelength = 0`, the spatial phase offset is omitted (uniform oscillation everywhere).

---

## Crypt Refugia

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `advection.crypts_enabled` | false | — | Enable crypt zero-flow zones |
| `advection.crypt_depth` | 10e-6 | m | Depth of crypt zone below epithelium (z < lo_z + depth) |
| `advection.crypt_exit_rate` | 1e-4 | 1/s | Per-second probability of agent exiting crypt |
| `advection.crypt_entry_rate` | 5e-5 | 1/s | Per-second probability of agent entering crypt |
| `advection.crypt_carrying_capacity` | 50 | — | Maximum agents in crypt region |

**Biological context:** Intestinal crypts and the firmly adherent inner mucus layer provide physical refugia where advective flow is effectively zero (γ_flow ≈ 0). Agents in this zone bypass the Washout Trap criterion (μ_realized < γ_flow), enabling colonization persistence even under metabolic stress (VADI §80, §98–99).

**Config file keys:** `crypts_enabled`, `crypt_depth`, `crypt_exit_rate`, `crypt_entry_rate`, `crypt_carrying_capacity`.

---

## Chemical Species – z-Dependent Gradient

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `ChemicalSpec.z_gradient_enabled` | false (true for carbon) | — | Enable exponential z-decay from epithelium |
| `ChemicalSpec.z_gradient_lambda` | 25e-6 | m | Characteristic decay length |
| `carbon_z_gradient` | true | — | Config file key for carbon gradient |
| `carbon_z_lambda` | 25e-6 | m | Config file key for carbon decay length |
| `carbon.epithelial_boundary` / `carbon_epithelial_boundary` | `dirichlet` | — | Epithelial z=0 mode: fixed concentration, Robin delivery, or fixed flux |
| `carbon.epithelial_transfer_coeff` / `carbon_epithelial_transfer_coeff` | 0 | m/s | Robin mass-transfer coefficient `k` |
| `carbon.epithelial_flux` / `carbon_epithelial_flux` | 0 | mol/m²/s | Fixed epithelial delivery flux `J` |
| `uptake_limit` / `metabolism.uptake_limit` | `none` | — | Agent-side uptake limitation model: `none` (unfunded growth, default), `sherwood` (per-agent diffusive delivery cap using the configured delivery concentration), or `voxel` (diagnostic-only voxel-content cap, not biology) |
| `delivery_far_field_radius` / `metabolism.delivery_far_field_radius` | `1.0e-5` | m | Shared delivery neighborhood radius, defaulting to a 10 µm physical support for both the volume-weighted Sherwood concentration read and prescribed delivery deposition. A single-cell depletion at 10 µm is approximately 0.008% of far-field concentration, so this radius is not doing biological work. Supports clipped by nonperiodic z boundaries are renormalized over included cells. Positive values are refused with slab chemistry; slab configurations must explicitly set this key to `0.0` to opt into the grid-dependent single-voxel model. |

When delivery-mode prescribed sinks would make an owned concentration
negative, rationing first acts locally: for up to eight retries, prescribed
mass is halved in the physical `delivery_far_field_radius` neighbourhood of
every negative owned cell. The neighbourhood uses the delivery support
geometry (periodic x/y and clipped z), rather than a voxel-count radius. If
the local pass cannot change any affected prescribed value or does not restore
positivity, a uniform global factor is used only as a final guarantee. Its
largest feasible value is found by 12 bisection iterations and applied
uniformly. Feasibility decisions are collective under MPI.

`delivery_reduction` is the original owned prescribed mass minus the final
owned prescribed mass. Retry events include both local and bisection solves.
The emitted rationing factor is the minimum factor among owned cells in the
interval, not a mean. If non-delivery sources remain negative at global factor
zero, the run continues, emits a diagnostic, and increments the infeasible
counter. In replicated chemistry, delivery reduction is combined with
`MPI_MAX` because every rank holds the same global field; in slab chemistry,
owned-cell reductions use `MPI_SUM`. Collective retry and infeasible counters
use `MPI_MAX`, while rationing factors use `MPI_MIN`. These delivery counters
are reduced once per chemistry step.

This default changes shipped delivery-mode behavior: prescribed mass is now
deposited over a 10 µm support unless radius zero is selected explicitly.
Settled population-scale measurements (80 founders, ROS off, six hours) found
support-on final populations of 201/207/204 at 2/4/6 µm, with funded oxygen
fraction 0.988 at all three resolutions. The corresponding support-off values
were final populations 76/123/171 and funded fractions 0.620/0.901/0.948.
Previously reported supply-limited carrying capacities, including the carbon
ladder in #314, were measured under single-voxel deposition and move under
this default.

**Profile:** `C(z) = C_max * exp(-z_rel / lambda)` where `z_rel` is the distance from the epithelium (z=0).

**Biological basis:** Mucin-derived monosaccharides are liberated primarily at the epithelial surface where host goblet cells secrete mucin glycoproteins. The anaerobic background degrades mucin locally, so the concentration of free monosaccharides is highest near z=0 and decays exponentially into the lumen. A characteristic length of ~25 μm places most carbon within the inner mucus layer.

### Nutrient Diffusion (Spec 7)

Nutrients and small molecules use backward-Euler directional splitting in `ChemicalField::apply_diffusion`. Each x/y line is solved with periodic boundaries and the luminal z face has zero flux. The epithelial z=0 face defaults to fixed concentration (`dirichlet`); `robin` uses `J = k(C_epi - c0)` and `flux` uses a fixed `J`, with the bottom cell included as an unknown in both delivery modes. The solve is L-stable and concentrations are clamped nonnegative, so it runs once per biological timestep without explicit CFL substeps. Boundary exchange is recorded after the solve in mol. `z_gradient_enabled` is incompatible with Robin and flux modes because its pinned reference profile assumes a Dirichlet boundary.

| Species | Effective configuration | Diffusion enabled |
|---------|-------------------------|-------------------|
| Carbon | `D_free = 5e-10 m²/s` | yes |
| Iron | `D_free = 7e-10 m²/s` | yes |
| Corrinoid (B12) | `D_free = 5e-10 m²/s` | yes |
| Oxygen | `oxygen.D_free` (default `2.1e-9 m²/s`) | when oxygen is enabled |
| Acetate | `acetate.D_free` (default `1.2e-9 m²/s`) | yes |
| Ethanolamine | `D_free = 1e-9 m²/s` | yes |
| Siderophore | `siderophore.D_free` (default `1e-10 m²/s`) | when siderophore is enabled |
| Mucin | immobile polymer field | no |
| Bacteriocins | analytical QSSA Green's functions | no grid diffusion |

`ChemicalSpec.diffusion_enabled` is currently a programmatic species property rather than an input-file key. Reactions from rank-local agents are summed with `MPI_Allreduce` before the shared VBF coupling and diffusion solve. CUDA runs apply diffusion on the host and synchronize the resulting concentrations back to the device.

---

## Nutrient Cycle (Spec 6)

Spec 6 makes the **metabolism Fix the single canonical site for per-agent
nutrient uptake**, resolving a double-count in which carbon/iron/B12 were
consumed both by `FixMetabolism::grow_agent` (yield-based, `reac -= d_biomass *
yield_X / (cell_vol * dt)`) and again by `QSSASolver::solve_nutrient_depletion`
(stoichiometry-based). The QSSA carbon/iron/B12 terms and the GPU
`nutrient_depletion_kernel` have been removed; `solve_nutrient_depletion` now
applies **only** aerobic O₂ respiration (which has no counterpart in the
metabolism Fix).

| Nutrient | Uptake site | Field behavior |
|----------|-------------|----------------|
| Carbon | metabolism Fix (`yield_carbon`) | Sourced by VBF mucin liberation, bounded by the VBF Monod carbon sink (now active by default, ~1 mM equilibrium) |
| Iron | metabolism Fix (`yield_iron`) + siderophore coupling | VBF first-order sink; Fur-regulated receptor uptake |
| O₂ | `solve_nutrient_depletion` — Pirt respiration `(q_consumption × μ_realized + q_maintenance) / cell_vol` | Supplied at the epithelial boundary; per-agent respiration (growth-associated **+ density-coupled maintenance**) + first-order VBF background sink |
| Corrinoid (B12) | **not consumed** | Constant field pinned at 1 µM = `1e-3 mol/m³` (see below) |

**Corrinoid (B12) is a constant pool, not a depletable field (Spec 6 §3).** The
B12 species now represents the *total bioavailable corrinoid pool* (~1 µM,
`initial_conc = boundary_conc = 1e-3`), the great majority of which are
non-cobalamin analogs produced by the anaerobic majority at rates far exceeding
E. coli demand. It is neither produced nor consumed in the model, so the field
stays pinned at 1 µM. (This replaces the Spec 5 `vbf_b12_production` source,
which has been removed; `yield_b12` is retained in config for compatibility but
no longer removes corrinoid from the field.)

**Competitive binding & colicin E (Receptor Ligand Parameterization).** BtuB is
both the corrinoid importer and the colicin-E receptor, so ambient corrinoid
competitively blocks colicin E: `apparent_Kd = kd_colicinE_btuB × (1 + [corrinoid] / kd_corrinoid_btuB)`.
With genuine concentration units, the corrinoid field is `1e-3 mol/m³` and
`kd_corrinoid_btuB` is `1e-6 mol/m³`, giving a factor of 1001. Raising the
corrinoid field from 1 nM to 1 µM increases this competitive factor
by ~1000×, making colicin E markedly less potent — the single most consequential
downstream effect of the nutrient-cycle rework. `kd_corrinoid_btuB`
(alias of `kd_b12_btuB`, default `1e-6 mol/m³` = 1 nM) is flagged as the **key unknown**; a
sweep over `{1e-6 … 1e-3}` mol/m³ is recommended future work (out of scope for this
change).

---

## VBF (Viscoelastic Background Field)

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `vbf.density` | 1e11 | #/m^3 | Anaerobic background density |
| `vbf.drag_coeff` | 1e-9 | N·s/m | Stokes drag coefficient |
| `vbf.nutrient_sink` | 1e-4 | 1/s | First-order iron uptake rate constant (sink is `-nutrient_sink · [iron]`, concentration-dependent — **not** a zero-order mol/m³/s removal) |
| `vbf.mucin_liberation` | 5e-5 | mol/m^3/s | Peak monosaccharide release (at z=0) |
| `vbf.carrying_cap` | 1e12 | #/m^3 | Local carrying capacity |
| `vbf.viscosity` | 0.01 | Pa·s | Effective viscosity (~10× water) |
| `vbf.mucin_z_gradient_enabled` | true | — | z-dependent mucin liberation rate |
| `vbf.mucin_z_gradient_lambda` | 25e-6 | m | Liberation decay length from epithelium |
| `vbf_carbon_sink_vmax` | 5.5e-5 | mol/m³/s | Monod carbon consumption by the anaerobic majority (Spec 5 §1 / Spec 6 §1). Activated by default: set just above the mucin liberation rate (5e-5) so bulk carbon settles to a ~1 mM equilibrium instead of accumulating without bound. `0` restores pre-Spec-6 unbounded accumulation |
| `vbf_carbon_sink_km` | 1e-4 | mol/m³ | Half-saturation for the VBF carbon sink |
| `vbf.agent_carbon_coupling` | 0.0 | mol/s/agent | Additional VBF carbon competition proportional to owned live-agent density; zero preserves the historical sink. For 2 µm cells (8 fL), use demand-anchored values such as 0, 1e-21, 1e-20, and 1e-19 mol/s/agent; values large enough to starve a voxel are expected to trigger the delivery closure gate. |

**Mucin liberation profile:** When `mucin_z_gradient_enabled`, the liberation rate varies as:
`rate(z) = mucin_liberation * exp(-z_rel / mucin_z_gradient_lambda)`

This ensures the carbon source term fed into the chemical field is strongest near the epithelium and decays toward the lumen, consistent with the chemical concentration gradient.

---

## QSSA Solver

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `qssa.toxin_cutoff` | 200e-6 | m | Green's function evaluation radius for toxins |
| `chemistry.toxin_evaluation` | `grid` | — | Toxin exposure evaluation: `grid` samples each agent's chemistry cell center; `agents` evaluates the same Green's-function superposition at each agent position |
| `chemistry.toxin_lumping` | `per_receptor` | — | Scientific toxin spatial model: `per_receptor` keeps four target-specific fields; `lumped` puts all sources in one field read by every receptor |
| `qssa.nutrient_cutoff` | 50e-6 | m | Cutoff for nutrient depletion zones |
| `qssa.colicin_release_rate` | 1e-18 | mol/s | Burst release per lysed cell |
| `qssa.microcin_secretion` | 1e-20 | mol/s | Continuous secretion rate |

### Barnes-Hut Acceleration

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `qssa.use_fmm` | false | — | Enable FMM octree for far-field aggregation |
| `qssa.fmm_theta` | 0.5 | — | Opening angle parameter (0→exact, 1→fast/approximate) |
| `qssa.fmm_expansion_order` | 2 | — | Multipole order: 1=monopole, 2=dipole+quadrupole, 3=octupole |

**Scaling note:** At 10^6 agents, naive O(N × M) evaluation is expensive. The cutoff radius limits each grid cell to nearby sources only, giving effective O(N) via spatial hashing. When `use_fmm` is true, distant sources beyond the cutoff are aggregated via a kernel-independent FMM with Cartesian multipole expansions (M2M, M2L, L2L), giving O(N+M) far-field cost after preprocessing. The opening angle `fmm_theta` controls accuracy: smaller values are more accurate but slower. Increase `fmm_expansion_order` for tighter error bounds (~theta^p) without changing theta. Typical values: `fmm_theta` 0.3 (conservative), 0.5 (balanced), 0.7 (fast); `fmm_expansion_order` 2 (default).

`chemistry.toxin_evaluation` is a modelling variant, not a performance toggle.
`grid` is the default and preserves the existing behaviour and fingerprints.
`agents` uses the actual agent position, so it can change killing when a cell is
offset from its chemistry-cell centre. Scheduled grid output is still
materialised in agent mode, but `bacteriocin_max_*` summary values become maxima
over sampled agent positions rather than over every grid cell. GPU execution
currently rejects `agents` mode because sampling remains CPU-only.

`chemistry.toxin_lumping` is also a modelling position, not a performance flag.
`per_receptor` is the default and preserves the existing four-field behaviour:
each toxin threatens an agent through the receptor it targets. `lumped` carries
the superposition of all toxin sources in one `bacteriocin_lumped` field and
lets every receptor read that local burden, while receptor-specific binding
parameters, toxin affinities, and immunity remain unchanged in the killing
calculation. This deliberately asks whether receptor-specific spatial targeting
matters, rather than changing the per-colicin biology. Its modelling cost is
that a colicin to which an agent is immune can contribute to the exposure it
experiences through a different receptor. In lumped mode the four
`bacteriocin_max_*` summary keys remain present and are identical by
construction; grid output contains one bacteriocin dataset.

`chemistry.species_subset` is likewise a modelling position, not a
performance flag. `full` is today's complete species set. `nutrient_only`
removes all bacteriocin fields and turns off bacteriocin/QSSA and receptor
killing, isolating nutrient competition. `carbon_only` retains carbon alone
and turns off toxin/receptor chemistry, oxygen, acetate, ethanolamine, mucin,
siderophore, ferrichrome, quorum sensing, motility taxis terms that read
removed fields, Fur regulation, iron/B12/eut uptake terms, the VBF iron sink,
and dynamic mucin; it asks whether spatial carbon competition alone
reproduces retention and clustering. Setting the VBF iron sink to zero and
disabling dynamic mucin are intentional parts of this modelling position, not
hidden parameter overrides. The chosen subset and exhaustive
disabled-mechanism list are printed in the startup audit line. Required-species
validation is always on: if an enabled mechanism lacks a required species,
initialization throws a `ConfigError` naming the mechanism, species, and
specific configuration key to change. A missing species is not a silent
mechanism disable.

---

## Protease Degradation (Spec 1)

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `protease.enabled` | true | — | Apply first-order protease decay to toxin QSSA sources |
| `protease.default_half_life` | 1800 | s | Default colicin half-life when not set on BI cluster |
| `protease.dilution_rate` | 1e-4 | 1/s | Fallback dilution for steady-state microcin decay |

Per-colicin `protease_half_life` is set on each `BICluster` in the plasmid library (ColE1/E2: 1800 s, ColB/M: 900 s, ColIa: 2400 s, MccV: 7200 s).

---

## Oxygen (Spec 1)

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `oxygen.enabled` | false | — | Enable oxygen chemical species and aerobic growth boost |
| `oxygen.epithelial_conc` | 55e-6 | mol/m³ | Dirichlet O₂ at epithelium (~0.042 mmHg); under Robin, the tissue-side reservoir concentration |
| `oxygen.epithelial_boundary` / `oxygen_epithelial_boundary` | `dirichlet` | — | Epithelial z=0 mode: fixed concentration, Robin delivery, or fixed flux |
| `oxygen.epithelial_transfer_coeff` / `oxygen_epithelial_transfer_coeff` | 0 | m/s | Robin mass-transfer coefficient `k`; must be positive for Robin. The sourced experiment value is `1.2e-6 m/s`, not a new default |
| `oxygen.epithelial_flux` / `oxygen_epithelial_flux` | 0 | mol/m²/s | Fixed epithelial delivery flux `J` for flux mode |
| `oxygen.z_gradient` / `oxygen_z_gradient` | true | — | Imposed exponential z-profile toggle; must be false for Robin or flux |
| `oxygen.respiration_driver` | `ambient` | — | Fermentation-mode driver: `ambient` preserves concentration-based switching; `funded` uses realized growth-O₂ delivery and requires oxygen delivery uptake with `metabolism.uptake_limit=delivery` |
| `oxygen.ros_driver` | `ambient` | — | SOS ROS driver: `ambient` uses ambient voxel oxygen (the compatibility default); `funded` uses funded respiratory oxygen flux and requires oxygen delivery uptake with `metabolism.uptake_limit=delivery` |
| `oxygen.D_free` | 2.1e-9 | m²/s | O₂ diffusion coefficient |
| `oxygen.Km` | 1e-6 | mol/m³ | Monod half-saturation for aerobic boost |
| `oxygen.boost_max` | 2.0 | — | Max growth multiplier above fermentation baseline |
| `oxygen.q_consumption` | 1e-14 | mol/cell | Growth-associated agent O₂ consumption (× μ_realized) |
| `oxygen.q_maintenance` | 1e-18 | mol/s/cell | Basal (density-coupled) agent O₂ respiration, applied per living cell regardless of growth so the field tracks agent density |
| `oxygen.vbf_sink` | 1e-3 | 1/s | VBF background O₂ uptake — **first-order** rate constant (`reac −= vbf_sink × [O₂]`), not a zero-order removal |
| `oxygen.k_ROS` | 1e2 | — | Uncited ambient-path ROS induction coefficient. Its units are not established in the existing documentation; its product with ambient oxygen concentration set mortality in every run to date. |
| `oxygen.k_ROS_respiratory` | 0.0 | — | Dimensionless funded ROS-induction yield per unit specific respiratory oxygen flux; required to be positive with `oxygen.ros_driver=funded` |

The default `oxygen.epithelial_conc = 55e-6 mol/m³` is approximately
`0.042 mmHg` dissolved oxygen, not the approximately `42 mmHg` stated by the
former comment. A dissolved concentration corresponding to approximately
`42 mmHg` would be about `5.5e-2 mol/m³`. The numerical default remains
unchanged; experiment configurations may set the sourced tissue-side value
`C_tissue = 5.0e-2 mol/m³` explicitly. The sourced Robin transfer coefficient
is `k = 1.2e-6 m/s` (rat mucus+mucosa permeability); neither value is a new
default.

---

## Dynamic Acetate (Spec 1)

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `acetate.enabled` | false | — | Enable dynamic acetate production/consumption (static 80 mM when false) |
| `acetate.D_free` | 1.2e-9 | m²/s | Acetate diffusion coefficient |
| `acetate.vbf_production` | 1e-3 | mol/m³/s | VBF fermentation acetate source |
| `acetate.vbf_consumption` | 2e-4 | mol/m³/s | VBF cross-feeding sink |
| `acetate.overflow_threshold` | 3e-4 | 1/s | Growth rate above which agents overflow acetate |
| `acetate.overflow_rate` | 1e-15 | mol/s/cell | Overflow secretion per agent |
| `acetate.scavenge_rate` | 1e-15 | mol/s/cell | Max acetate scavenging per agent |
| `acetate.scavenge_Km` | 5.0 | mol/m³ | Half-saturation for scavenging |
| `acetate.epithelial_uptake` | 5e-4 | mol/m³/s | Colonocyte uptake at z=0 |

---

## Dynamic Mucin (Spec 1)

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `mucin.enabled` | false | — | Track mucin polymer field and dynamic liberation |
| `mucin.initial_conc` | 1e-2 | mol/m³ | Initial mucin concentration |
| `mucin.secretion_rate` | 1e-4 | mol/m³/s | Goblet cell secretion at epithelium |
| `mucin.Km_degradation` | 1e-3 | mol/m³ | Half-saturation for VBF mucin degradation |
| `mucin.k_liberation` | 1e-4 | 1/s | Rate constant for mucin → monosaccharide conversion |

---

## Metabolism

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `metabolism.division_threshold` | 2.0 | — | Biomass ratio for division |
| `metabolism.metE_penalty` | 0.05 | — | MetE pathway base cost (BtuB loss) |
| `metabolism.metE_acetate_km` | 40.0 | mol/m³ | Half-saturation for acetate inhibition of MetE |
| `metabolism.metE_acetate_max_factor` | 2.5 | — | Max scaling factor at saturating acetate |
| `metabolism.eut_km` | 0.1e-3 | mol/m³ | Ethanolamine half-saturation for eut utilization penalty |
| `metabolism.eut_max_penalty` | 0.10 | — | Max eut penalty when ethanolamine abundant |
| `metabolism.maintenance_rate` | 1e-5 | 1/s | Maintenance energy |
| `metabolism.carbon_maintenance_rate` | 0.0 | mol C/(s·kg biomass) | Non-growth-associated carbon substrate consumption; independent of `maintenance_rate` |
| `oxygen.metabolic_switch_enabled` | false | bool | Replace legacy oxygen boost with respiratory, overflow, and fermentative metabolism |
| `oxygen.mu_crit` | 9.7e-5 | 1/s | Respiratory capacity threshold; with `resp_capacity = f_O2 * mu_crit / max(mu, mu_crit)`, this 0.35 h⁻¹ onset keeps the growth-rate/overflow axis active for example strains with `mu_max = 5.5e-4 /s` |
| `oxygen.aerobic_mu_factor` | 1.0 | multiplier | Growth factor at zero fermentation |
| `oxygen.anaerobic_mu_factor` | 0.55 | multiplier | Growth factor at full fermentation |
| `oxygen.aerobic_carbon_cost_factor` | 1.0 | multiplier | Multiplier on substrate-per-biomass carbon cost |
| `oxygen.anaerobic_carbon_cost_factor` | 4.1 | multiplier | Fermentative multiplier on substrate-per-biomass carbon cost |
| `oxygen.tau_metabolic_switch` | 3600 | s | Realized metabolic-state time constant |
| `oxygen.ferm_acid_yield` | 1.0 | acetate/carbon-equivalent | Literal acetate yield from fermentation |
| `oxygen.anaerobic_maintenance_factor` | 15.0 | multiplier | Maintenance factor at full fermentation |
| `acid_inhibition_enabled` | false | bool | Enable undissociated-acetate growth inhibition |
| `acid_inhibition_max` | 0.8 | fraction | Maximum acid inhibition |
| `acid_inhibition_Ki` | 50 | mol/m³ | Half-inhibition concentration of undissociated acetate |
| `acetate_pKa` | 4.76 | pH units | Acetate dissociation pKa |
| `metabolism.km_iron_primary` | 10e-6 | mol/m³ | FepA iron Km (10 nM) |
| `metabolism.km_iron_iroN` | 50e-6 | mol/m³ | IroN salmochelin Km (50 nM) |
| `metabolism.km_iron_iutA` | 100e-6 | mol/m³ | IutA aerobactin Km (100 nM) |
| `metabolism.km_iron_fiu` | 200e-6 | mol/m³ | Fiu catecholate Km (200 nM) |

**Graded iron uptake:** Iron acquisition uses multiple receptor systems in parallel. FepA is primary (highest affinity), with IroN, IutA, and Fiu as secondary fallbacks. When FepA is downregulated (e.g. to resist colicin B), cells switch to these secondary pathways rather than experiencing complete iron starvation. The effective iron Monod term sums contributions from all receptors weighted by expression level, then normalizes to preserve wild-type growth at full expression.

**MetE pathway:** When BtuB expression < 0.5, cells must synthesize methionine via the MetE pathway instead of the B12-dependent MetH pathway. MetE requires ~5% of the proteome.

**Ethanolamine utilization loss:** The eut operon is B12-dependent. When BtuB is downregulated the penalty is proportional to local [ethanolamine] via Monod kinetics: `eut_effect = eut_max_penalty * [EA] / (eut_km + [EA])`. In inflamed gut (high ethanolamine) the cost is much larger than at homeostatic levels.

**Acetate inhibition of MetE (VADI §87):** Colonic acetate (60–100 mM) severely inhibits the MetE enzyme. The effective penalty is scaled by local acetate concentration via Michaelis-Menten kinetics:
```
metE_eff = metE_penalty * (1 + (max_factor - 1) * [acetate] / (Km + [acetate]))
```
At physiological colonic acetate (80 mM, Km = 40 mol/m³), the effective penalty rises from 5% to ~10%. At saturating acetate the penalty approaches 12.5% (`base × max_factor`). This strengthens the Combinatorial Washout Trap by increasing the metabolic cost of BtuB downregulation in acetate-rich environments.

`maintenance_rate` is a growth-rate tax. It is distinct from
`carbon_maintenance_rate`, which removes carbon substrate at
`rate × biomass × dt` even when realized growth is zero or negative. The
default carbon rate is exactly zero for backward compatibility.

The non-growth-associated maintenance anchors are:

- Pirt's maintenance coefficient is approximately `0.04 g glucose/gDW/h`,
  or `0.22 mmol glucose/gDW/h`.
- Varma & Palsson's independent W3110 NGAM cross-check is
  `7.6 mmol ATP/gDW/h`; at approximately `30 ATP/glucose` aerobically this is
  approximately `0.25 mmol glucose/gDW/h`. The two estimates agree within
  approximately `15%`.
- With `30%` dry fraction, `CELL_DENSITY_DEFAULT = 1100 kg/m³`, and glucose
  molecular weight `180 g/mol`, the aerobic anchor is approximately
  `2.1e-5 mol/(s·kg biomass)`, or `1.7e-20 mol/s` for the `8.14e-16 kg`
  Sherwood-campaign agent.
- At approximately `2 ATP/glucose` fermentatively, the anaerobic requirement is
  approximately `15×` higher: `3.2e-4 mol/(s·kg)`, or
  `2.6e-19 mol/s/cell`. This is the same order as the agent's full-growth
  carbon demand of `2.0e-19 mol/s`.

When `oxygen.metabolic_switch_enabled` is enabled, carbon maintenance
interpolates between these aerobic and anaerobic anchors using the realized
fermentation fraction. With the switch disabled, the configured maintenance
rate is unchanged.

Fermentative secretion currently writes literal acetate. Canonical mixed-acid
fermentation is approximately `1 acetate + 1 ethanol + 2 formate per glucose`;
a lumped total-acid model would therefore be near three acid equivalents and
would require a lower effective pKa. This implementation deliberately models
acetate only. The implemented `agent_carbon_coupling` parameter adds a
per-owned-agent VBF carbon sink; see
[`SPEC12_DENSITY_LIMITATION.md`](SPEC12_DENSITY_LIMITATION.md) for sizing and
closure-gate interactions.
the measured carbon maintenance sink is the model's density-coupled substrate
draw.
Finally, the default `yield_carbon = 0.5 mol/kg` is approximately `6.4×`
cheaper than the literature comparison. Varma's `0.524 gDW/g glucose`, with
`30%` dry fraction and `180 g/mol`, implies approximately `3.2 mol glucose per
kg wet biomass`. The default yield remains unchanged.

The proposed multi-scale architecture, which would make these patch-level
parameters per-patch and per-segment, is not implemented; see
[`SPEC13_MULTISCALE.md`](SPEC13_MULTISCALE.md) and
[`SPEC13_IMPLEMENTATION_REVIEW.md`](SPEC13_IMPLEMENTATION_REVIEW.md).

The acid-inhibition `Ki` above is expressed in undissociated acetate units.
The cited Russell & Diez-Gonzalez threshold is approximately 50 mM
undissociated acid, matching the repository's 50 mol/m³ default. The amended
specification's 20 mol/m³ total-acetate value would be about 1.1 mol/m³
undissociated at pH 6, roughly 45 times more sensitive, so it is not adopted.
Likewise, `oxygen.ferm_acid_yield` is acetate per carbon-equivalent consumed,
not acetate per mole of glucose; the specification's 0.67 mol/mol-glucose
quantity is therefore different rather than contradictory.

---

## Fur-Regulated Receptors (Spec 3)

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `fur.enabled` | true | — | Enable Fur-regulated dynamic receptor expression |
| `fur.Km` | 1e-5 | mol/m³ | Iron concentration for half-max Fur repression |
| `fur.upregulation_max` | 10.0 | — | Max fold-upregulation under iron starvation (Spec 6 §4.2; raised 4→10, still conservative vs measured 35–56× Fur-regulon induction; capped by `receptor_max`) |
| `fur.receptor_max` | 5.0 | — | Cap on effective receptor expression |

When enabled, iron-uptake receptors (FepA, FhuA, IroN, IutA, Fiu, CirA) are upregulated under low local iron, increasing colicin susceptibility (Vulnerability Paradox). Mutations modify `receptor_expr_base`; Fur scales effective `receptor_expr` each metabolism step. Fur receptor modulation and per-agent acetate coupling run in the device metabolism path when GPU metabolism is active; per-cell siderophore field chemistry remains host-side.

---

## Contact-Dependent Inhibition (Spec 3)

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `cdi.enabled` | true | — | Enable CDI contact killing |
| `cdi.kill_rate` | 5e-4 | 1/s | Killing rate per contact pair |
| `cdi.contact_radius` | 1e-6 | m | Max CDI delivery distance |
| `cdi.corpse_persistence` | 300 | s | Dead-cell obstacle lifetime |

Per-strain JSON keys: `cdi_type`, `cdi_immunity` on `initial_strains` entries. CDI kills set `death_time` for delayed corpse removal; other death paths remove agents immediately.

---

## Active Motility (Spec 3 / Spec 10v2)

Run-and-reverse swimming with modular behavioral modes. Directional taxis
(aerotaxis, carbon chemotaxis) modulate **run duration** via Weber–Fechner
fractional sensing (`ΔC/(C·dt)`). Energy taxis, surface sensing, and mucin
drag modulate **swim speed** multiplicatively.

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `motility.enabled` | true | — | Enable active swimming |
| `motility.swim_speed` | 7.76e-6 | m/s | Mean swimming speed in mucus |
| `motility.run_mean_duration` | 1.0 | s | Mean run duration |
| `motility.stop_probability` | 0.3 | — | P(stop) per reorientation |
| `motility.stop_duration` | 0.5 | s | Mean stop duration |
| `motility.chemotaxis_enabled` | false | — | Enable carbon (Tar/Tsr) chemotaxis |
| `motility.chi_carbon` | 2.0 | — | Weber–Fechner carbon sensitivity |
| `motility.chemotaxis_threshold` | 1e-6 | mol/m³ | Floor for fractional sensing |
| `motility.aerotaxis_enabled` | true | — | Enable Aer-mediated aerotaxis on O₂ |
| `motility.aerotaxis_sensitivity` | 4.0 | — | Weber–Fechner O₂ sensitivity (primary cue) |
| `motility.energy_taxis_enabled` | true | — | Reduce speed under metabolic stress |
| `motility.energy_taxis_floor` | 0.1 | — | Speed fraction at `mu_realized = 0` |
| `motility.surface_sensing_enabled` | false | — | Reduce speed near epithelium |
| `motility.surface_sensing_depth` | 10e-6 | m | Depth of surface-sensing zone |
| `motility.surface_sensing_floor` | 0.3 | — | Speed fraction at z = epithelium |
| `motility.mucin_drag_enabled` | false | — | Viscosity drag from mucin field |
| `motility.mucin_drag_reference` | 1e-2 | mol/m³ | Mucin conc at half-speed |
| `motility.cluster_suppress_radius` | 10e-6 | m | Cluster detection radius |
| `motility.cluster_suppress_threshold` | 5 | — | Neighbors to suppress tumbling |
| `motility.cluster_tumble_factor` | 0.2 | — | Tumble rate multiplier in cluster center |

`motility.chi_oxygen` was removed in Spec 10v2; use `motility.aerotaxis_sensitivity`.

---

## AI-2 Quorum Sensing (Spec 11)

Opt-in autoinducer-2 field with LuxS production, Lsr import, and optional
Weber–Fechner chemotaxis (via `fix_motility`).

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `quorum_sensing.enabled` | false | — | Master switch; registers `ai2` species |
| `quorum_sensing.ai2_basal_rate` | 1e-20 | mol/s/cell | Constitutive LuxS production |
| `quorum_sensing.ai2_growth_coupled` | 1e-16 | mol/cell | Growth-associated production (`× mu_realized`) |
| `quorum_sensing.lsr_vmax` | 1e-18 | mol/s/cell | Max Lsr import rate |
| `quorum_sensing.lsr_km` | 1e-7 | mol/m³ | Lsr half-saturation (~100 nM) |
| `quorum_sensing.ai2_D_free` | 5e-10 | m²/s | AI-2 free diffusion coefficient |
| `quorum_sensing.ai2_decay_rate` | 1e-4 | 1/s | Background first-order decay |
| `quorum_sensing.ai2_chemotaxis` | false | — | Enable AI-2 run-length bias in motility |
| `quorum_sensing.chi_ai2` | 3.0 | — | Weber–Fechner AI-2 sensitivity |

---

## Receptor Binding

| Parameter | Config key | Default | Units | Description |
|-----------|------------|---------|-------|-------------|
| `receptor.kd_b12_btuB` | `kd_b12_btuB` / `kd_corrinoid_btuB` | 1e-6 | mol/m^3 | BtuB affinity for the dominant corrinoid analog (1 nM; Spec 6 / Receptor Ligand Parameterization). **Key unknown**: with the corrinoid pool at ~1 µM this Kd governs how strongly corrinoid competitively blocks colicin E at BtuB (see note below). `kd_corrinoid_btuB` is an alias for the same field |
| `receptor.kd_colicinE_btuB` | `kd_colicinE_btuB` | 5e-7 | mol/m^3 | Colicin E affinity for BtuB (0.5 nM) |
| `receptor.kd_enterobactin` | `kd_enterobactin` | 1e-6 | mol/m^3 | FepA affinity for ferric enterobactin (1 nM; Spec 6 §4.1) |
| `receptor.kd_colicinB_fepA` | `kd_colicinB_fepA` | 2e-6 | mol/m^3 | Colicin B affinity for FepA (2 nM) |
| `receptor.kd_lin_enterobactin` | `kd_lin_enterobactin` | 5e-5 | mol/m^3 | Linearized enterobactin for CirA (50 nM) |
| `receptor.kd_colicinIa_cirA` | `kd_colicinIa_cirA` | 3e-6 | mol/m^3 | Colicin Ia affinity for CirA (3 nM) |
| `receptor.kd_colicinM_fhuA` | `kd_colicinM_fhuA` | 2.5e-6 | mol/m^3 | Colicin M affinity for FhuA (2.5 nM) |
| `receptor.kd_ferrichrome` | `kd_ferrichrome` | 1e-5 | mol/m^3 | Ferrichrome affinity for FhuA (10 nM) |
| `receptor.kill_rate_colicin` | `kill_rate_colicin` | 1e-3 | 1/s | Single-hit colicin kill rate |
| `receptor.kill_rate_microcin` | `kill_rate_microcin` | 5e-4 | 1/s | Microcin kill rate (slower) |
| `receptor.cirA_linearized_fraction` | `cirA_linearized_fraction` | 0.3 | — | Fraction of ferric enterobactin represented as the CirA linearized ligand |
| `chemical.b12_initial_conc` | `b12.initial_conc`, `b12_initial_conc`, `corrinoid.initial_conc`, `corrinoid_initial_conc` | 1e-3 | mol/m³ | Uniform initial and boundary concentration of the corrinoid/B12 pool |

---

## Bacteriocin

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `bacteriocin.sos_lysis_prob` | 0.01 | — | SOS induction probability per division (active when `just_divided`) |
| `bacteriocin.sos_basal_rate` | 1e-6 | 1/s | Spontaneous SOS rate |
| `bacteriocin.sos_cross_induction_rate` | 1e3 | 1/s per mol/m³ | Nuclease provoker rate, applied to per-agent exposure sampled from nuclease BI sources |
| `bacteriocin.D_free_colicin` | 4e-11 | m^2/s | Free diffusion (~50kDa protein) |
| `bacteriocin.burst_release_tau` | 300 | s | Exponential release timescale; total delivered dose is tau-invariant, while tau sets peak versus duration |
| `bacteriocin.microcin_mu_penalty` | 0.03 | — | Growth cost of microcin secretion |
| `bacteriocin.mucin_charge.r_min` | 1.2 | — | Minimum bacteriocin retardation |
| `bacteriocin.mucin_charge.amplitude` | 60.0 | — | Maximum charge-driven retardation excess |
| `bacteriocin.mucin_charge.dz_half` | 1.35 | pI − pH | Half-transition net-charge offset |
| `bacteriocin.mucin_charge.width` | 1.0 | pI − pH | Sigmoid width in net-charge units |
| `bacteriocin.mucin_charge.ph` | 7.0 | pH | Reference local mucus pH |

For each library cluster, retardation is resolved from pI and local pH as:

`R(pI) = r_min + amplitude / (1 + 10^((dz_half - (pI - pH)) / width))`

The resolved value is serialized with the BI cluster and is therefore preserved
through MPI transfer and checkpoint restore. Per-plasmid retardation overrides
remain authoritative. Per-plasmid defaults in `PlasmidLibrary` also include
`release_mode`, `is_nuclease`, `diff_coeff`, `burst_size`, and
`phage_induction_rate` (ColB/ColIa: 1e-4 /generation). Optional overrides use:

```json
"plasmid_overrides": {
  "ColE1": {
    "retardation": 25.0,
    "diff_coeff": 8.0e-11,
    "burst_size": 200000.0
  }
}
```

Names resolve through `PlasmidLibrary::find()`, so existing aliases such as
`colicin_E1` remain accepted and are emitted canonically. The removed global
`retardation_basic`, `retardation_acidic`, and `retardation_neutral` keys are
not used; per-plasmid retardation is authoritative. Unknown keys warn by
default, while `GUTIBM_STRICT_CONFIG=1` makes them hard configuration errors.

QSSA maintains four receptor-specific toxin fields (`bacteriocin_BtuB`, `bacteriocin_FepA`, `bacteriocin_CirA`, `bacteriocin_FhuA`) by default. With `chemistry.toxin_lumping = "lumped"`, it instead maintains one `bacteriocin_lumped` field containing all sources, and all receptors read that field. Nuclease colicin bursts continue to use the BtuB cross-induction path in the default mode; lumped mode makes that path read total toxin burden as well. The lumped approximation means a colicin an agent is immune to can contribute to exposure through another receptor, and can also contribute to SOS/ROS cross-induction; these are intentional scientific costs of removing receptor-specific spatial targeting.

---

## Siderophore (Spec 4)

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `siderophore.enabled` | true | — | Register `siderophore` species and secretion/chelation in metabolism; enabled by default because iron-limited *E. coli* commonly produces enterobactin |
| `siderophore.secretion_rate` | 1e-5 | mol/(s·kg) | Constrained estimate for enterobactin secretion, scaled by Fur activity |
| `siderophore.D_free` | 1e-10 | m²/s | Free siderophore diffusion |
| `siderophore.chelation_rate` | 1e3 | m³/mol/s | Iron–siderophore chelation sink |
| `siderophore.Km_reimport` | 1e-6 | mol/m³ | FepA-mediated ferric-enterobactin reimport; converted from the erroneous `1e-9` default |
| `siderophore.Vmax_reimport` | 1e-5 | mol/(s·kg) | FepA ferric-enterobactin transport capacity, grounded in FepA copy number and TonB-limited turnover |

Siderophore chemistry is enabled by default because enterobactin production is
near-universal among iron-limited *E. coli*. Enabling it does not restore
meaningful FepA competition for an isolated cell: diffusion carries apo
enterobactin away faster than local chelation, so the realized
`ferric_enterobactin` concentration remains far below
`receptor.kd_enterobactin`. In a maintained-background single-cell assay
(`grid_dx = 5 µm`), the local source-cell FeEnt concentration immediately
after the reaction update was `1e-11 mol/m³` at `1e-4 mol/m³` iron and
`6e-15 mol/m³` at `1e-8 mol/m³` iron. The previously reported `8e-16 mol/m³`
was the post-diffusion field average; diffusion spreads the reaction-stage
source pulse across the grid before the next biological step.

The corrected chemistry was measured with 1, 4, 16, and 64 agents co-located
in one grid cell while maintaining the local apo-enterobactin background.
Reaction-stage FeEnt decreased with occupancy:

| Iron condition | N=1 | N=4 | N=16 | N=64 |
|---|---:|---:|---:|---:|
| `1e-4 mol/m³`, apo `4e-9 mol/m³` | `5e-12` | `1e-12` | `3e-13` | `7e-14` |
| `1e-8 mol/m³`, apo `3e-8 mol/m³` | `1e-15` | `3e-16` | `8e-17` | `2e-17` |

Values are in `mol/m³`. Chelation is a solution-phase reaction and is now
applied in every grid cell containing apo-enterobactin and iron, independent
of occupancy. Secretion and FepA reimport remain biomass-weighted and are
confined to occupied cells. Increasing local biomass therefore increases the
aggregate FepA sink without multiplying the solution-phase chelation term,
driving local FeEnt down rather than up. Diffusive escape further removes apo
and FeEnt from the source cell. In the EARI/VADI run, diffusion and
solution-phase chelation produce a domain-wide background with mean FeEnt
`6.37e-8 mol/m³` and maximum `1.07e-7 mol/m³` at the final saved step.
Relative to `kd_enterobactin = 1e-6 mol/m³`, that is a mean competition factor
of approximately `1.064` and a maximum of approximately `1.107`.
The competition is therefore a spatially flat background effect rather than
local co-location-driven protection: increasing occupancy decreases the local
source-cell FeEnt, so it cannot create spatial structure in colicin B
susceptibility. It is a global parameter shift wearing the costume of a
spatial mechanism.

This differs from the colicin result: bacteriocin dose scales with co-located
producers because lysing producers are not also a sink for the toxin. The
microcolony threshold identified for colicin in #213 is therefore a real
density effect, whereas FeEnt competition here is a domain-wide detuning and
local co-location-driven competition remains unreachable.

Chelation consumes apo-enterobactin (`siderophore`) and free iron and produces
`ferric_enterobactin`. FepA reimport consumes ferric enterobactin and returns
iron; both secretion and reimport are specific rates in mol/(s·kg), multiplied
by biomass density. The secretion value is a constrained estimate rather than
a direct measurement. The reimport value uses approximately 35,000 FepA per
iron-starved cell and approximately five transport cycles per minute per FepA
(Smallwood et al. 2016; Newton et al. 2010), converted using the default cell
mass from `CELL_RADIUS_DEFAULT` and `CELL_DENSITY_DEFAULT`. The selected
`1e-5 mol/(s·kg)` is a rounded capacity; the exact default-cell-mass
conversion of the cited estimate is `8.33e-6 mol/(s·kg)`.

## Ferrichrome ambient field

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `ferrichrome.enabled` | false | — | Register the ambient ferrichrome field |
| `ferrichrome.initial_conc` | 0 | mol/m³ | Initial bulk ferrichrome concentration |
| `ferrichrome.boundary_conc` | 0 | mol/m³ | Boundary ferrichrome concentration |

Ferrichrome is an exogenous ligand: *E. coli* does not synthesize it, so it is
not a secretion product. The default is deliberately disabled because no
defensible gut ferrichrome concentration is currently available.

`receptor.cirA_linearized_fraction` (default 0.3) scales ferric-enterobactin
concentration as the CirA linearized ligand when siderophore chemistry is enabled.

---

## Mechanics (Cell-Cell Contact)

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `mechanics.hertz_k` | 1e-6 | N/m^1.5 | Hertzian spring constant |
| `mechanics.hertzian_enabled` | true | — | Use Hertzian (F∝δ^1.5) vs linear (F∝δ) |
| `mechanics.adhesion_enabled` | false | — | Enable EPS-mediated adhesion |
| `mechanics.adhesion_strength` | 1e-12 | N | Maximum adhesion force |
| `mechanics.adhesion_range` | 0.5e-6 | m | Range beyond contact for adhesion |

**Hertzian contact model:** `F = hertz_k * overlap^(3/2)` where `overlap = r_i + r_j - d`. Calibrated from AFM measurements of bacterial elastic modulus (~0.1–1 MPa). Only applies when cells physically overlap (`overlap > 0`).

Mechanical mobility uses the existing `vbf_viscosity` value, so that
parameter also sets the Stokes drag scale for contact and adhesion relaxation.
The mechanics summary writes `mechanics/displacement_clamps` and
`mechanics/cumulative_displacement_clamps` to expose use of the per-step
displacement safety cap.

`fixes.metabolism.bacteriostasis_threshold` is the threshold used to classify
the instantaneous `bacteriostatic_live_agents` stock. Its default of `1e-6 /s`
identifies a viable non-reproducing state on the model's simulation horizons.
It does not kill cells or cause starvation mortality.
Each agent's accumulated mechanics displacement is capped at `0.1 * radius`
per biological step.

**EPS adhesion:** When enabled, cells within `adhesion_range` of contact experience an attractive force that decays linearly with gap distance. Models extracellular polymeric substance (EPS) bridging for biofilm-like clustering.

---

## Conjugation

| Parameter | Config key | Default | Units | Description |
|-----------|------------|---------|-------|-------------|
| `conjugation.pili_length` | `pili_length` | 4e-6 | m | Max F-pilus reach (heterogeneity off) |
| `conjugation.base_transfer_rate` | `base_transfer_rate` | 1e-4 | 1/s | Conjugation events per s per pair |
| `conjugation.shear_critical` | `shear_critical` | 10.0 | 1/s | Critical shear for MPS |
| `conjugation.plasmid_copy_cost` | `plasmid_copy_cost` | 0.02 | — | Metabolic cost per transferred plasmid |
| `conjugation.pili_heterogeneity` | `pili_heterogeneity` | false | — | Enable per-event F-pilus length sampling |
| `conjugation.pili_length_min` | `pili_length_min` | 1e-6 | m | Min F-pilus length (uniform lower bound) |
| `conjugation.pili_length_max` | `pili_length_max` | 4e-6 | m | Max F-pilus length (uniform upper bound) |

**Pili heterogeneity (VADI §55):** In vivo F-pili are 1–4 μm with significant length heterogeneity. When `pili_heterogeneity = true`, each conjugation attempt samples its effective contact radius from `uniform(pili_length_min, pili_length_max)` instead of using the fixed `pili_length`. The expected mean reach is 2.5 μm.

---

## Mutation

| Parameter | Config key | Default | Units | Description |
|-----------|------------|---------|-------|-------------|
| `mutation.bi_duplication_rate` | `bi_duplication_rate` | 1e-5 | per division | BI locus duplication |
| `mutation.bi_recombination_rate` | `bi_recombination_rate` | 5e-6 | per division | BI locus recombination |
| `mutation.receptor_mutation_rate` | `receptor_mutation_rate` | 1e-7 | per division | Receptor downregulation |
| `mutation.super_killer_rate` | `super_killer_rate` | 1e-8 | per division | Novel toxin variant |
| `mutation.compensatory_rate` | `compensatory_rate` | 1e-6 | per division | Plasmid cost amelioration |
| `mutation.receptor_reduction` | `receptor_reduction` | 0.1 | — | Expression drop per mutation |
| `mutation.partial_resistance_rate` | `partial_resistance_rate` | 5e-7 | per division | Extracellular loop missense mutation |
| `mutation.compensatory_reduction` | `compensatory_reduction` | 0.005 | — | Per-locus cost reduction |
| `mutation.max_bi_loci` | `max_bi_loci` | 8 | — | Maximum BI clusters per genome |
| `mutation.immunity_escape_prob` | `immunity_escape_prob` | 0.5 | — | Fraction of super-killers with immunity escape |
| `mutation.escape_affinity_lo` | `escape_affinity_lo` | 0.01 | — | Lower bound of reduced binding affinity |
| `mutation.escape_affinity_hi` | `escape_affinity_hi` | 0.3 | — | Upper bound of reduced binding affinity |

---

## OpenMP Threading

| Parameter | Default | Description |
|-----------|---------|-------------|
| `GUTIBM_USE_OPENMP` | `OFF` | CMake option to enable OpenMP parallelism |

**Build with OpenMP:**
```bash
cmake -B build -DGUTIBM_USE_OPENMP=ON
cmake --build build -j$(nproc)
```

**Thread control (runtime):**
```bash
export OMP_NUM_THREADS=8     # limit thread count
export OMP_SCHEDULE=dynamic  # override schedule policy
```

**Parallelized regions:**
- Green's function superposition (`superpose_to_grid`): thread-local accumulation with dynamic scheduling
- QSSA per-receptor bacteriocin field deposition (`bacteriocin_BtuB/FepA/CirA/FhuA`): static schedule over grid cells
- QSSA nutrient depletion: dynamic schedule over agents with atomic grid updates
- Agent metabolism (growth rate + biomass): static schedule, atomic nutrient consumption
- Receptor kill probability: static schedule (precomputed in parallel, applied serially)
- Agent advection and physics: static schedule over agents
- Chemical field reaction application: static schedule over grid cells
- Grid coupling update: static schedule over agents

**Thread safety:** Shared chemical field updates use `#pragma omp atomic` to prevent
race conditions when multiple agents occupy the same grid cell. The mechanical repulsion
loop (pairwise neighbor interactions) remains serial due to cross-agent writes.

---

## GPU Acceleration

| Parameter | Default | Description |
|-----------|---------|-------------|
| `GUTIBM_USE_CUDA` | `OFF` | CMake option to compile CUDA kernels |
| `gpu_enabled` | `false` | Runtime toggle (input file); requires CUDA build |
| `gpu_device_id` | `-1` | CUDA device index; `-1` = `MPI_rank % num_devices` |
| `profile_steps` | `false` | Per-step wall-clock profiling on rank 0 (see `docs/SCALING.md`) |

**Build with CUDA:**
```bash
cmake -B build -DGUTIBM_USE_MPI=ON -DGUTIBM_USE_HDF5=ON -DGUTIBM_USE_CUDA=ON
cmake --build build -j$(nproc)
```

**Example input file:**
```
gpu_enabled true
gpu_device_id 0
```

**Memory:** Default domain (~12.5M cells × 6 species × 2 arrays) uses ~1.2 GB VRAM for the chemical field mirror alone. Plan for ≥4 GB VRAM for typical agent counts.

**Parity:** GPU results use relaxed floating-point tolerance vs CPU (see `tests/test_gpu_smoke.cpp`). Bit-identical reproducibility is not guaranteed.

---

## HDF5 Output (Spec 4 layered schema)

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `hdf5.filename` / `hdf5_file` | `gut_ibm_output.h5` | — | Output file path |
| `hdf5.enabled` | true | — | Master switch (also off when all schedule intervals are 0) |
| `hdf5.schedule.summary` | 1 | steps | Per-step summary stats + globally reduced interval and cumulative event counters, including `mortality_lysis`, plus instantaneous `stocks/` and `mechanics/` groups; the latter contains `displacement_clamps` and `cumulative_displacement_clamps` |
| `hdf5.schedule.agents` | 5 | steps | Lightweight agent arrays |
| `hdf5.schedule.grid` | 0 | steps | 3D chemical grids (0 = disabled) |
| `hdf5.schedule.lineage` | 100 | steps | Lineage tracker arrays |
| `hdf5.schedule.genome` | 100 | steps | Full genome / BI locus tables |
| `hdf5.schedule.provenance` | 0 | steps | Per-kill provenance records (0 = disabled) |
| `hdf5.compression` | `gzip` | — | Grid compression: `none` or `gzip` |
| `hdf5.compression_level` | 4 | — | gzip level (0–9) when compression is `gzip` |
| `hdf5.parallel` | false | — | MPI-parallel agent gather on rank 0 |

### Closed midstream restarts (Tier 2)

Separate from the live analysis trail. When enabled, Simulation writes
`{directory}/step_NNNNNN.h5` (create → agents/lineage/genome/grid/summary → close)
every `interval_steps`, plus a final file on exit. AWS `entry.sh` uploads these
immutably and points `latest.json` at the newest.

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `restart.enabled` | false | — | Master switch for closed restart writes |
| `restart.directory` | `restart` | — | Local directory for `step_*.h5` files |
| `restart.interval_steps` | 0 | steps | Cadence (0 = off even if enabled) |
| `checkpoint_file` | — | — | Path to an existing restart/analysis HDF5 to resume |
| `checkpoint_step` | — | — | Optional step group; empty = latest |

Nested JSON example:

```json
"restart": {
  "enabled": true,
  "directory": "restart",
  "interval_steps": 60
}
```

HDF5 nested JSON example:

```json
"hdf5": {
  "file": "output.h5",
  "compression": "gzip",
  "schedule": {
    "summary": 1,
    "agents": 5,
    "grid": 0,
    "lineage": 100,
    "genome": 100,
    "provenance": 0
  }
}
```

File layout: `/summary/step_NNNNNN/`, `/agents/step_NNNNNN/`, `/grid/step_NNNNNN/` (3D, optional gzip), `/lineage/`, `/genome/`. File attributes include `gutibm_version=4`, `nx`, `ny`, `nz`, and `grid_dx_x`, `grid_dx_y`, `grid_dx_z`; legacy `grid_dx` remains readable.

Every HDF5 file also contains a root-level `/run_provenance/` record written
once when the file is first populated. `resolved_config` is the complete
post-default, post-override configuration in parser-compatible JSON, so it
can be saved and fed back to `gut_ibm` to reproduce or fork the run.
`git_sha`, `version`, compile-time feature flags, and `mpi_rank_count` identify
the binary and execution shape. `container_image_digest` and `job_id` are
included only when supplied through `GUTIBM_IMAGE_DIGEST` and
`AWS_BATCH_JOB_ID`; absent optional datasets are not fabricated.
After termination, the same group records `termination_reason_code`,
`termination_step`, and `termination_time`. A dysbiosis halt has
`termination_reason_code=1` and additionally records
`halt_reason_code=1`, `halt_density_cells_per_mL`, `halt_step`, and
`halt_time`. A run that reaches `total_time` records
`completed_total_time=1` and `termination_reason_code=0`; this makes a
guard-censored artifact distinct from a complete-horizon artifact.
The complete termination-code enumeration is `0` for reaching `total_time`,
`1` for a dysbiosis guard halt, and `2` for another early exit such as the
population-stop path.

`/run_provenance/` is deliberately separate from `/provenance/`. The latter
remains the per-kill mechanism record controlled by
`hdf5.schedule.provenance`; it is not run identity or configuration metadata.
The same run-provenance record is present in normal output files and closed
midstream restart/checkpoint files.

---

## Fix Plugins

Optional `fixes` array in the input JSON selects which Fix modules run and in what order.
When omitted, all registered defaults are used in registry order:

`metabolism` → `bacteriocin` → `receptor` → `conjugation` → `mutation` → `mechanics`

| Config key | Maps to | Description |
|------------|---------|-------------|
| `fixes` | `SimulationConfig::enabled_fixes` | JSON string array of Fix names |

Example (growth + mechanics only):

```json
"fixes": ["metabolism", "mechanics"]
```

Unknown Fix names log a warning and are skipped. Register new Fix modules in
`src/fixes/fix_registry.cpp` without editing `simulation.cpp`.

---

## Input Config Format

Simulation configs are strict JSON. Use `"_comment"` (string or array) for
human-readable notes — see [CONFIG_FORMAT.md](CONFIG_FORMAT.md).

## Initial Population

Each strain in `initial_strains` has:

| Field | Default | Description |
|-------|---------|-------------|
| `type` | — | Integer strain identifier |
| `count` | — | Number of initial agents |
| `mu_max` | `5e-4` | Maximum specific growth rate (1/s) for the strain's agents. This is the **only** place the max growth rate is configured — it is a per-strain property (`Agent::mu_max`), scaled each step by the Monod terms in `FixMetabolism` (`mu = mu_max · monod_carbon · monod_iron · monod_b12`). There is no global `metabolism` default growth rate. |
| `plasmids` | `[]` | List of plasmid names (from `PlasmidLibrary`) |
| `conjugative` | `false` | Whether the strain can initiate conjugation (HGT) |
| `cdi_type` | `0` | CDI system identifier delivered by this strain (`0` = none); see [Contact-Dependent Inhibition](#contact-dependent-inhibition-spec-3) |
| `cdi_immunity` | `0` | CDI immunity identifier this strain carries (`0` = none) |
| `receptor_expression` | `{}` | Map of receptor names (`BtuB`, `FepA`, `Tsx`, `FhuA`, `IroN`, `Fiu`, `CirA`, `IutA`) to expression levels in `[0, 1]`; values below `0.2` make the founder resistant |

Example:
```json
{
  "type": 1,
  "count": 1000,
  "mu_max": 5e-4,
  "plasmids": ["ColE1"],
  "conjugative": false
}
```

Delivery-mode chemistry first retries a negative solve with local reductions in
the affected cells. If local reductions cannot restore positivity, it bisects
a single global prescribed-mass factor for 12 iterations and applies the
largest feasible factor. MPI feasibility decisions are collective.
`delivery_reduction_*` remains the original owned prescribed mass minus the
final owned prescribed mass; retry counters include local and bisection solves.
The nutrient-flux summary also emits the minimum rationing factor over the
interval and run, plus an infeasible-step counter for negativity that remains
at factor zero. Delivery uptake is CPU-only; GPU plus delivery uptake is
rejected during configuration because GPU chemistry does not implement the
same rationing loop.
