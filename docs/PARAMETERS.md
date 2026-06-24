# Parameter Reference

All configurable parameters for the GutIBM simulation, with defaults and biological justification.

---

## Time Control

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `total_time` | 86400 | s | Simulation duration (24 h) |
| `bio_dt` | 60 | s | Biological timestep |
| `output_interval` | 3600 | s | HDF5 dump interval |
| `seed` | 42 | — | Random number generator seed |

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

This captures shear-enhanced spreading of toxins in the mucus flow.

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

---

## VBF (Viscoelastic Background Field)

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `vbf.density` | 1e11 | #/m^3 | Anaerobic background density |
| `vbf.drag_coeff` | 1e-9 | N·s/m | Stokes drag coefficient |
| `vbf.nutrient_sink` | 1e-4 | mol/m^3/s | Background nutrient consumption |
| `vbf.mucin_liberation` | 5e-5 | mol/m^3/s | Peak monosaccharide release (at z=0) |
| `vbf.carrying_cap` | 1e12 | #/m^3 | Local carrying capacity |
| `vbf.viscosity` | 0.01 | Pa·s | Effective viscosity (~10× water) |
| `vbf.mucin_z_gradient_enabled` | true | — | z-dependent mucin liberation rate |
| `vbf.mucin_z_gradient_lambda` | 25e-6 | m | Liberation decay length from epithelium |

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
| `qssa.use_fmm` | false | — | Enable Barnes-Hut octree for far-field aggregation |
| `qssa.fmm_theta` | 0.5 | — | Opening angle parameter (0→exact, 1→fast/approximate) |

**Scaling note:** At 10^6 agents, naive O(N × M) evaluation is expensive. The cutoff radius limits each grid cell to nearby sources only, giving effective O(N) via spatial hashing. When `use_fmm` is true, distant sources beyond the cutoff are aggregated via a Barnes-Hut octree monopole approximation, reducing total cost to O(N log N). The opening angle `fmm_theta` controls accuracy: smaller values are more accurate but slower. Typical values: 0.3 (conservative), 0.5 (balanced), 0.7 (fast).

---

## Metabolism

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `metabolism.mu_max_default` | 5e-4 | 1/s | Default max growth rate |
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

## Receptor Binding

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `receptor.kd_b12_btuB` | 1e-9 | mol/m^3 | B12 affinity for BtuB |
| `receptor.kd_colicinE_btuB` | 5e-10 | mol/m^3 | Colicin E affinity for BtuB |
| `receptor.kd_enterobactin` | 1e-8 | mol/m^3 | Enterobactin affinity for FepA |
| `receptor.kd_colicinB_fepA` | 2e-9 | mol/m^3 | Colicin B affinity for FepA |
| `receptor.kd_lin_enterobactin` | 5e-8 | mol/m^3 | Linearized enterobactin for CirA |
| `receptor.kd_colicinIa_cirA` | 3e-9 | mol/m^3 | Colicin Ia affinity for CirA |
| `receptor.kill_rate_colicin` | 1e-3 | 1/s | Single-hit colicin kill rate |
| `receptor.kill_rate_microcin` | 5e-4 | 1/s | Microcin kill rate (slower) |
| `receptor.immunity_factor` | 0.001 | — | 1000× immunity protection |

---

## Bacteriocin

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `bacteriocin.sos_lysis_prob` | 0.01 | — | SOS induction per division |
| `bacteriocin.sos_basal_rate` | 1e-6 | 1/s | Spontaneous SOS rate |
| `bacteriocin.retardation_basic` | 50.0 | — | R for pI > 8.5 (Lethal Core) |
| `bacteriocin.retardation_acidic` | 1.5 | — | R for pI < 7.0 (Lethal Halo) |
| `bacteriocin.retardation_neutral` | 5.0 | — | R for 7.0–8.5 |
| `bacteriocin.D_free_colicin` | 4e-11 | m^2/s | Free diffusion (~50kDa protein) |
| `bacteriocin.burst_molecules` | 1e4 | — | Molecules per lysis burst |
| `bacteriocin.microcin_mu_penalty` | 0.03 | — | Growth cost of microcin secretion |

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

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `conjugation.base_transfer_prob` | 1e-3 | — | Base HGT probability |
| `conjugation.contact_radius` | 2e-6 | m | Max pilus reach (used when heterogeneity off) |
| `conjugation.shear_crit` | 10.0 | 1/s | Critical shear for MPS |
| `conjugation.pili_heterogeneity` | false | — | Enable per-event F-pilus length sampling |
| `conjugation.pili_length_min` | 1e-6 | m | Min F-pilus length (uniform lower bound) |
| `conjugation.pili_length_max` | 4e-6 | m | Max F-pilus length (uniform upper bound) |

**Pili heterogeneity (VADI §55):** In vivo F-pili are 1–4 μm with significant length heterogeneity. When `pili_heterogeneity = true`, each conjugation attempt samples its effective contact radius from `uniform(pili_length_min, pili_length_max)` instead of using the fixed `contact_radius`. The expected mean reach is 2.5 μm.

---

## Mutation

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `mutation.bi_duplication_rate` | 1e-5 | per division | BI locus duplication |
| `mutation.bi_recombination_rate` | 5e-6 | per division | BI locus recombination |
| `mutation.receptor_mutation_rate` | 1e-7 | per division | Receptor downregulation |
| `mutation.super_killer_rate` | 1e-8 | per division | Novel toxin variant |
| `mutation.compensatory_rate` | 1e-6 | per division | Plasmid cost amelioration |
| `mutation.receptor_reduction` | 0.1 | — | Expression drop per mutation |
| `mutation.partial_resistance_rate` | 5e-7 | per division | Extracellular loop missense mutation |
| `mutation.compensatory_reduction` | 0.005 | — | Per-locus cost reduction |
| `mutation.max_bi_loci` | 8 | — | Maximum BI clusters per genome |
| `mutation.immunity_escape_prob` | 0.5 | — | Fraction of super-killers with immunity escape |
| `mutation.escape_affinity_lo` | 0.01 | — | Lower bound of reduced binding affinity |
| `mutation.escape_affinity_hi` | 0.3 | — | Upper bound of reduced binding affinity |

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
- QSSA bacteriocin field deposition: static schedule over grid cells
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

## HDF5 Output

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `hdf5.filename` | `gut_ibm_output.h5` | — | Output file path |
| `hdf5.dump_every` | 100 | steps | Steps between dumps |
| `hdf5.parallel` | false | — | MPI-parallel I/O (each rank writes its agents via collective hyperslab) |
| `hdf5.enabled` | true | — | Set false to suppress output |

---

## Initial Population

Each strain in `initial_strains` has:

| Field | Description |
|-------|-------------|
| `type` | Integer strain identifier |
| `count` | Number of initial agents |
| `mu_max` | Maximum growth rate |
| `plasmids` | List of plasmid names (from PlasmidLibrary) |
| `conjugative` | Whether the strain can conjugate |

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
