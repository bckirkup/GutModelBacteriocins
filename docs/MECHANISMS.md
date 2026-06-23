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

---

## Fix Architecture (NUFEB-inspired)

Each biological rule is encapsulated as a **Fix** — a modular computation unit called once per biological timestep. Fixes are executed in order:

1. `fix_metabolism` — Monod growth, division, death
2. `fix_bacteriocin` — SOS lysis, microcin secretion, toxin release
3. `fix_receptor` — Competitive binding at TBDTs, toxin-mediated killing
4. `fix_conjugation` — Horizontal gene transfer via F-pili
5. `fix_mutation` — Stochastic BI locus evolution, receptor downregulation

The QSSA solver and advection field operate between fix passes (chemistry and physics modules respectively).

---

## 1. fix_metabolism — Triple Monod Growth

**Biological basis:** *E. coli* growth in the gut requires carbon (mucin-derived monosaccharides), iron (via siderophores through multiple TBDTs), and vitamin B12 (via BtuB). Growth rate follows multiplicative Monod kinetics.

**Equation:**
```
mu = mu_max * [C]/(Km_C + [C]) * monod_iron * [B12]/(Km_B12 + [B12])
```

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

**Note:** FhuA (ferrichrome) is NOT included as a secondary iron fallback because it transports fungal ferrichrome, which is not an endogenous enterobactin pathway (corrects EARI §70).

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

---

## 2. fix_bacteriocin — Toxin Release

**Two pathways implemented:**

### a) SOS-mediated lysis (colicins)
Large-protein colicins (30–80 kDa) are released only upon cell death:
- Spontaneous SOS induction at basal rate (`sos_basal_rate`, default 10^-6 /s)
- After induction, 5-minute delay (`sos_timer`), then lysis
- Burst release of ~10^4 toxin molecules as instantaneous point source
- The QSSA solver treats lysed cells as burst sources for Green's function superposition

### b) Continuous microcin secretion (new)
Small peptide microcins (<10 kDa, e.g. MccV) are exported without lysis:
- Applies a growth rate penalty (`microcin_mu_penalty`, default 3%) to producing cells
- Contributes steady-state continuous sources to the toxin field
- Per EARI model Section 51: microcins impose 2–5% mu_max cost

**pI-dependent diffusion classification:**
| Class | pI range | Retardation | Behavior |
|-------|----------|-------------|----------|
| Lethal Core | > 8.5 | R = 50 | Binds mucin glycoproteins, concentrates near producer |
| Lethal Halo | < 7.0 | R = 1.5 | Repelled by anionic mucin, spreads widely |
| Neutral | 7.0–8.5 | R = 5.0 | Intermediate diffusion |

*Note:* The HALO threshold was raised from pI < 6.0 to pI < 7.0 (VADI §74) to correctly classify bacteriocins secreted as acidic complexes (e.g. Colicin E2 with Im2, net complex pI ≈ 6.5). The `classify_by_pI()` function in `src/genome/plasmid.h` is the single source of truth for these thresholds.

The effective diffusion coefficient: `D_eff = D_free / R`

---

## 3. fix_receptor — Competitive Binding at TBDTs

**The Double-Bind hypothesis:** TonB-dependent transporters serve dual roles as nutrient importers AND bacteriocin receptors. Cells cannot selectively block toxins without also losing nutrient uptake.

**Three receptor systems modeled:**

| Receptor | Nutrient ligand | Toxin ligand | Kd_nutrient | Kd_toxin |
|----------|----------------|--------------|-------------|----------|
| BtuB | Vitamin B12 | Colicin E1/E2/E5/E7/E8/E9/K | 1 nM | 0.5 nM |
| FepA | Enterobactin-Fe | Colicin B/D | 10 nM | 2 nM |
| CirA | Linearized enterobactin | Colicin Ia, Microcin V | 50 nM | 3 nM |

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

## 6. fix_mechanics — Soft-Sphere Mechanical Repulsion

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

**Force application:** Equal-and-opposite impulses using reduced mass weighting:
```
push_i = F * dt / m_i * (m_j / (m_i + m_j))
push_j = F * dt / m_j * (m_i / (m_i + m_j))
```
This ensures momentum conservation and correct size-dependent response.

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

## QSSA Solver — Quasi-Steady-State Diffusion

**Why not FTCS?** At 1 um grid resolution, the CFL stability criterion forces sub-millisecond timesteps for explicit diffusion. Since small-molecule diffusion equilibrates in microseconds-to-minutes while cell biology operates on hour timescales, we exploit this separation via QSSA.

**Method:** At each biological timestep:
1. Collect all active toxin sources (SOS-lysed cells as bursts, microcin producers as steady sources)
2. Evaluate the analytical Green's function at each grid cell:
   ```
   C(r) = Q / (4 * pi * D_eff * r)    (3D point source, steady state)
   ```
3. Apply Method of Images for bounded mucus domain (z = 0 epithelium, z = h lumen)
4. Superpose contributions from all sources within cutoff radius

**Taylor-Aris dispersion:** The shear profile in the mucus layer enhances longitudinal spreading:
```
D_eff = D_mol + U(z)^2 * h^2 / (210 * D_mol)
```
This is computed via `AdvectionField::taylor_aris_D_eff()`.

### Barnes-Hut Far-Field Acceleration

When `QSSAConfig::use_fmm = true`, the QSSA solver uses a Barnes-Hut octree to accelerate the Green's function superposition from O(N × M) to O(N log N):

1. **Build**: Construct an octree over all active toxin sources, computing monopole moments (total source strength, source-weighted centroid) at each internal node.
2. **Near-field** (distance ≤ `toxin_cutoff`): Evaluate exact advection-diffusion Green's function per-source, preserving accuracy for nearby lethal core/halo interactions.
3. **Far-field** (distance > `toxin_cutoff`): Traverse the octree; accept a node's monopole approximation when `cell_size / distance < theta` (opening-angle criterion). The monopole uses the node's aggregate strength at its source-weighted centroid with source-rate-weighted average diffusion parameters.

The opening angle `theta` controls the accuracy/speed trade-off. At `theta = 0`, all interactions are computed exactly (equivalent to the original method). At `theta = 0.5`, the monopole error is bounded by ~theta² for smooth fields, giving a good balance for bacteriocin concentration fields.

**Note:** This is a simplified Barnes-Hut (monopole-only) implementation. A full Fast Multipole Method (FMM) with higher-order multipole expansions would provide O(N) scaling with tighter error bounds — see derivative issue for future work.

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
The domain is divided into `nprocs` equal-width slabs along the x-axis. Each rank owns agents within its slab bounds `[local_lo_x, local_hi_x)`.

### Ghost Layers
Before the biology module, each rank exchanges **ghost agents** — copies of agents within `ghost_width` of slab boundaries — with its neighbors. This ensures correct neighbor queries (spatial hash lookups) for cross-boundary interactions (conjugation, mechanical repulsion, receptor binding).

Ghost agents are read-only and discarded before the physics module to avoid double-counting position updates.

### Agent Migration
After the physics module (advection + mechanics), agents that have moved past their slab boundary are serialized and sent to the appropriate neighbor rank via `MPI_Sendrecv`. The full agent state — position, velocity, metabolism, genome, BI clusters — is transferred.

### Global Statistics
`MPI_Allreduce` aggregates per-rank counts and growth rate sums to produce global agent count and mean growth rate. These are used for output and lineage tracking.

### HDF5 Parallel I/O
When `hdf5.parallel = true`, the file is opened with `H5Pset_fapl_mpio` and agent data is written using collective hyperslab operations — each rank writes its local agents at a computed offset. Grid, metadata, and lineage data are written by rank 0 only.

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

---

## z-Dependent Nutrient Gradient

**Biological basis:** Mucin-derived monosaccharides (the primary carbon source for *Enterobacteriaceae* in the gut) are liberated by goblet cells at the epithelial surface (z=0). The concentration of free monosaccharides is highest near the epithelium and decays exponentially into the lumen.

**Implementation:**

### Chemical field initialization
When `z_gradient_enabled` is set for a chemical species, the initial concentration follows:
```
C(z) = C_max * exp(-z_rel / lambda_mucin)
```
where `z_rel` is the distance from the epithelium and `lambda_mucin` is the characteristic decay length (~25 μm by default). The Dirichlet boundary at z=0 maintains `boundary_conc` as the peak value.

### VBF mucin liberation coupling
When `mucin_z_gradient_enabled`, the monosaccharide release rate applied by the VBF also varies with z:
```
rate(z) = mucin_liberation * exp(-z_rel / mucin_z_gradient_lambda)
```
This couples the ongoing carbon source term to the same spatial gradient, maintaining the profile during simulation.

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
