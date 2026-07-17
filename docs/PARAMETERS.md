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
| `dysbiosis_threshold` | 0 | cells/mL | Spec 5 §4 safety net: halt the run if global agent density exceeds this value. `0` disables the check |

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

---

## Domain

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `domain.lo` | {0,0,0} | m | Lower corner |
| `domain.hi` | {1e-3, 1e-3, 100e-6} | m | Upper corner (1mm × 1mm × 100um) |
| `domain.grid_dx` | 2e-6 | m | Grid cell size |
| `domain.hash_cell_size` | 10e-6 | m | Spatial hash bucket size |
| `domain.periodic` | {true, true, false} | — | Periodicity per axis |
| `domain.mpi_decomp_axis` | 0 | — | Axis for 1D slab decomposition (0=x) |
| `domain.ghost_width` | 10e-6 | m | Ghost layer thickness for cross-rank neighbor queries |

**Biological context:** The domain represents a patch of colonic mucus layer. x,y are periodic (infinite mucosa plane). z spans from epithelium (z=0) to luminal surface (z=h).

**MPI decomposition:** The domain is partitioned into equal-width slabs along `mpi_decomp_axis` (default: x, the distal flow direction). Each rank owns one slab. `ghost_width` should be ≥ `hash_cell_size` to ensure correct neighbor queries across slab boundaries.

---

## Advection

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `advection.radial_turnover` | 5400 | s | Mucus radial turnover (1.5 h) |
| `advection.mucus_thickness` | 100e-6 | m | Mucus layer depth |
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

**Profile:** `C(z) = C_max * exp(-z_rel / lambda)` where `z_rel` is the distance from the epithelium (z=0).

**Biological basis:** Mucin-derived monosaccharides are liberated primarily at the epithelial surface where host goblet cells secrete mucin glycoproteins. The anaerobic background degrades mucin locally, so the concentration of free monosaccharides is highest near z=0 and decays exponentially into the lumen. A characteristic length of ~25 μm places most carbon within the inner mucus layer.

### Nutrient Diffusion (Spec 7)

Nutrients and small molecules use backward-Euler directional splitting in `ChemicalField::apply_diffusion`. Each x/y line is solved with periodic boundaries; z uses a fixed epithelial concentration at z=0 and zero flux at the luminal face. The solve is L-stable and concentrations are clamped nonnegative, so it runs once per biological timestep without explicit CFL substeps. When `z_gradient_enabled` is true, the configured exponential profile is treated as a prescribed environmental background and diffusion smooths departures from that profile.

| Species | Effective configuration | Diffusion enabled |
|---------|-------------------------|-------------------|
| Carbon | `D_free = 5e-10 m²/s` | yes |
| Iron | `D_free = 7e-10 m²/s` | yes |
| Corrinoid (B12) | `D_free = 5e-10 m²/s` | yes |
| Oxygen | `oxygen.D_free` (default `2.1e-9 m²/s`) | when oxygen is enabled |
| Acetate | `acetate.D_free` (default `1.2e-9 m²/s`) | yes |
| Ethanolamine | `D_free = 1e-9 m²/s` | yes |
| Siderophore | `siderophore.D_free` (default `5e-10 m²/s`) | when siderophore is enabled |
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
| Corrinoid (B12) | **not consumed** | Constant field pinned at 1 µM (see below) |

**Corrinoid (B12) is a constant pool, not a depletable field (Spec 6 §3).** The
B12 species now represents the *total bioavailable corrinoid pool* (~1 µM,
`initial_conc = boundary_conc = 1e-6`), the great majority of which are
non-cobalamin analogs produced by the anaerobic majority at rates far exceeding
E. coli demand. It is neither produced nor consumed in the model, so the field
stays pinned at 1 µM. (This replaces the Spec 5 `vbf_b12_production` source,
which has been removed; `yield_b12` is retained in config for compatibility but
no longer removes corrinoid from the field.)

**Competitive binding & colicin E (Receptor Ligand Parameterization).** BtuB is
both the corrinoid importer and the colicin-E receptor, so ambient corrinoid
competitively blocks colicin E: `apparent_Kd = kd_colicinE_btuB × (1 + [corrinoid] / kd_corrinoid_btuB)`.
Raising the corrinoid field from 1 nM to 1 µM increases this competitive factor
by ~1000×, making colicin E markedly less potent — the single most consequential
downstream effect of the nutrient-cycle rework. `kd_corrinoid_btuB`
(alias of `kd_b12_btuB`, default 1 nM) is flagged as the **key unknown**; a
sweep over `{1e-9 … 1e-6}` is recommended future work (out of scope for this
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

**Mucin liberation profile:** When `mucin_z_gradient_enabled`, the liberation rate varies as:
`rate(z) = mucin_liberation * exp(-z_rel / mucin_z_gradient_lambda)`

This ensures the carbon source term fed into the chemical field is strongest near the epithelium and decays toward the lumen, consistent with the chemical concentration gradient.

---

## QSSA Solver

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `qssa.toxin_cutoff` | 200e-6 | m | Green's function evaluation radius for toxins |
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
| `oxygen.epithelial_conc` | 55e-6 | mol/m³ | Dirichlet O₂ at epithelium (~42 mmHg) |
| `oxygen.D_free` | 2.1e-9 | m²/s | O₂ diffusion coefficient |
| `oxygen.Km` | 1e-6 | mol/m³ | Monod half-saturation for aerobic boost |
| `oxygen.boost_max` | 2.0 | — | Max growth multiplier above fermentation baseline |
| `oxygen.q_consumption` | 1e-14 | mol/cell | Growth-associated agent O₂ consumption (× μ_realized) |
| `oxygen.q_maintenance` | 1e-18 | mol/s/cell | Basal (density-coupled) agent O₂ respiration, applied per living cell regardless of growth so the field tracks agent density |
| `oxygen.vbf_sink` | 1e-3 | 1/s | VBF background O₂ uptake — **first-order** rate constant (`reac −= vbf_sink × [O₂]`), not a zero-order removal |
| `oxygen.k_ROS` | 1e2 | — | ROS induction rate coefficient (Spec 2 hook) |

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

---

## Fur-Regulated Receptors (Spec 3)

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `fur.enabled` | true | — | Enable Fur-regulated dynamic receptor expression |
| `fur.Km` | 1e-5 | mol/m³ | Iron concentration for half-max Fur repression |
| `fur.upregulation_max` | 10.0 | — | Max fold-upregulation under iron starvation (Spec 6 §4.2; raised 4→10, still conservative vs measured 35–56× Fur-regulon induction; capped by `receptor_max`) |
| `fur.receptor_max` | 5.0 | — | Cap on effective receptor expression |

When enabled, iron-uptake receptors (FepA, FhuA, IroN, IutA, Fiu, CirA) are upregulated under low local iron, increasing colicin susceptibility (Vulnerability Paradox). Mutations modify `receptor_expr_base`; Fur scales effective `receptor_expr` each metabolism step. GPU metabolism fast-path is disabled when Fur is enabled.

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

## Receptor Binding

| Parameter | Config key | Default | Units | Description |
|-----------|------------|---------|-------|-------------|
| `receptor.kd_b12_btuB` | `kd_b12_btuB` / `kd_corrinoid_btuB` | 1e-9 | mol/m^3 | BtuB affinity for the dominant corrinoid analog (Spec 6 / Receptor Ligand Parameterization). **Key unknown**: with the corrinoid pool at ~1 µM this Kd governs how strongly corrinoid competitively blocks colicin E at BtuB (see note below). `kd_corrinoid_btuB` is an alias for the same field |
| `receptor.kd_colicinE_btuB` | `kd_colicinE_btuB` | 5e-10 | mol/m^3 | Colicin E affinity for BtuB |
| `receptor.kd_enterobactin` | `kd_enterobactin` | 1e-9 | mol/m^3 | FepA affinity for ferric enterobactin (Spec 6 §4.1; tightened 10 nM → 1 nM to match measured high-affinity siderophore uptake) |
| `receptor.kd_colicinB_fepA` | `kd_colicinB_fepA` | 2e-9 | mol/m^3 | Colicin B affinity for FepA |
| `receptor.kd_lin_enterobactin` | `kd_lin_enterobactin` | 5e-8 | mol/m^3 | Linearized enterobactin for CirA |
| `receptor.kd_colicinIa_cirA` | `kd_colicinIa_cirA` | 3e-9 | mol/m^3 | Colicin Ia affinity for CirA |
| `receptor.kill_rate_colicin` | `kill_rate_colicin` | 1e-3 | 1/s | Single-hit colicin kill rate |
| `receptor.kill_rate_microcin` | `kill_rate_microcin` | 5e-4 | 1/s | Microcin kill rate (slower) |
| `receptor.cirA_linearized_fraction` | `cirA_linearized_fraction` | 0.3 | — | Fraction of siderophore as CirA ligand |

---

## Bacteriocin

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `bacteriocin.sos_lysis_prob` | 0.01 | — | SOS induction probability per division (active when `just_divided`) |
| `bacteriocin.sos_basal_rate` | 1e-6 | 1/s | Spontaneous SOS rate |
| `bacteriocin.sos_cross_induction_rate` | 1e3 | 1/s per mol/m³ | Nuclease provoker rate (reads `bacteriocin_BtuB` field) |
| `bacteriocin.retardation_basic` | 50.0 | — | R for pI > 8.5 (Lethal Core) |
| `bacteriocin.retardation_acidic` | 1.5 | — | R for pI < 7.0 (Lethal Halo) |
| `bacteriocin.retardation_neutral` | 5.0 | — | R for 7.0–8.5 |
| `bacteriocin.D_free_colicin` | 4e-11 | m^2/s | Free diffusion (~50kDa protein) |
| `bacteriocin.burst_molecules` | 1e4 | — | Reference burst size for scaling per-BI `burst_size` |
| `bacteriocin.microcin_mu_penalty` | 0.03 | — | Growth cost of microcin secretion |

Per-plasmid defaults in `PlasmidLibrary`: `release_mode`, `is_nuclease`, `burst_size`, `phage_induction_rate` (ColB/ColIa: 1e-4 /generation).

QSSA maintains four receptor-specific toxin fields (`bacteriocin_BtuB`, `bacteriocin_FepA`, `bacteriocin_CirA`, `bacteriocin_FhuA`). Nuclease colicin bursts deposit into `bacteriocin_BtuB` for cross-induction. `fix_receptor` reads only the field matching each receptor target.

---

## Siderophore (Spec 4)

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `siderophore.enabled` | false | — | Register `siderophore` species and secretion/chelation in metabolism |
| `siderophore.secretion_rate` | 1e-15 | mol/s per biomass | Enterobactin secretion scaled by Fur activity |
| `siderophore.D_free` | 4e-11 | m²/s | Free siderophore diffusion |
| `siderophore.chelation_rate` | 1e3 | m³/mol/s | Iron–siderophore chelation sink |
| `siderophore.Km_reimport` | 1e-6 | mol/m³ | FepA-mediated siderophore–iron reimport |
| `siderophore.recapture_fraction` | 0.1 | — | Local iron recapture near producer |

`receptor.cirA_linearized_fraction` (default 0.3) scales siderophore concentration as the CirA nutrient ligand when siderophore is enabled.

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
| `hdf5.schedule.summary` | 1 | steps | Per-step summary stats + event counters |
| `hdf5.schedule.agents` | 5 | steps | Lightweight agent arrays |
| `hdf5.schedule.grid` | 0 | steps | 3D chemical grids (0 = disabled) |
| `hdf5.schedule.lineage` | 100 | steps | Lineage tracker arrays |
| `hdf5.schedule.genome` | 100 | steps | Full genome / BI locus tables |
| `hdf5.compression` | `none` | — | Grid compression: `none` or `gzip` |
| `hdf5.compression_level` | 4 | — | gzip level (0–9) when compression is `gzip` |
| `hdf5.parallel` | false | — | MPI-parallel agent gather on rank 0 |

Nested JSON example:

```json
"hdf5": {
  "file": "output.h5",
  "compression": "gzip",
  "schedule": {
    "summary": 1,
    "agents": 5,
    "grid": 0,
    "lineage": 100,
    "genome": 100
  }
}
```

File layout: `/summary/step_NNNNNN/`, `/agents/step_NNNNNN/`, `/grid/step_NNNNNN/` (3D, optional gzip), `/lineage/`, `/genome/`. File attributes include `gutibm_version=4`, `nx`, `ny`, `nz`, `grid_dx`.

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
