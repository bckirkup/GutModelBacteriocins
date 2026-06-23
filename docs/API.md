# API Reference

Class and function reference for the GutIBM simulation engine.

---

## Core Classes

### `Simulation`
**Header:** `src/core/simulation.h`

Main orchestration engine. Manages the timestep loop inspired by NUFEB's `nufeb_run`.

| Method | Description |
|--------|-------------|
| `init(SimulationConfig&)` | Initialize all subsystems from config |
| `run()` | Execute the full simulation loop |
| `step(Real dt)` | Single biological timestep |
| `agents()` | Access the `AgentPool` |
| `domain()` | Access the `Domain` |
| `chemical_field()` | Access the `ChemicalField` |
| `advection()` | Access the `AdvectionField` |
| `vbf()` | Access the VBF continuum |
| `qssa()` | Access the QSSA solver |
| `lineage_tracker()` | Access the `LineageTracker` |
| `rng()` | Access the random number generator |
| `time()` | Current simulation time (s) |
| `step_count()` | Number of completed steps |

**Timestep structure:**
```
step(dt):
  1. Pre-step: zero reactions, rebuild spatial hash, update grid coupling
  2. module_biology(dt)   — execute all Fix modules
  3. module_chemistry()   — QSSA steady-state fields
  4. module_physics(dt)   — advection, drag, mechanical repulsion
  5. Post-step: fix post-processing (division, lysis completion)
  6. Cleanup: check washout, remove dead agents
```

---

### `Agent`
**Header:** `src/core/agent.h`

Single bacterial cell. Each agent has spatial, metabolic, receptor, phenotype, and genomic state.

| Field | Type | Description |
|-------|------|-------------|
| `tag` | `TagID` | Globally unique identifier |
| `type` | `Int` | Species/phylogroup index |
| `x` | `Vec3` | Position in meters |
| `v` | `Vec3` | Velocity in m/s |
| `radius` | `Real` | Cell radius (m) |
| `biomass` | `Real` | Dry biomass (kg) |
| `mu_max` | `Real` | Intrinsic max growth rate (1/s) |
| `mu_realized` | `Real` | Current growth rate after penalties |
| `receptor_expr` | `array<Real, 7>` | Expression level per receptor (0–1) |
| `state` | `PhenoState` | NORMAL, RESISTANT, SOS_INDUCED, DEAD |
| `genome` | `Genome` | BI clusters, lineage, mutations |
| `age` | `Real` | Time since last division (s) |

**Factory:** `Agent::create_default(tag, type, pos, mu_max)` — creates a wild-type cell with default parameters.

---

### `AgentPool`
**Header:** `src/core/agent.h`

Vector-backed container for agents with O(1) removal (swap-and-pop).

| Method | Description |
|--------|-------------|
| `size()` | Number of agents |
| `push_back(Agent)` | Add an agent |
| `remove(Int idx)` | Remove by index (swap with last) |
| `next_tag()` | Generate next unique TagID |
| `operator[](Int)` | Index access |

---

### `Domain`
**Header:** `src/core/domain.h`

3D simulation domain with periodic boundary conditions (x,y) and bounded z (epithelium→lumen).

| Method | Description |
|--------|-------------|
| `init(DomainConfig&)` | Set up grid and spatial hash |
| `lo()`, `hi()` | Domain bounds (Vec3) |
| `nx()`, `ny()`, `nz()` | Grid dimensions |
| `ncells()` | Total grid cells |
| `apply_pbc(Vec3&)` | Apply periodic boundary conditions |
| `pos_to_grid(Vec3, ix, iy, iz)` | Position → grid indices |
| `cell_index(ix, iy, iz)` | Grid indices → flat index |
| `min_image_delta(a, b)` | Minimum-image displacement |
| `spatial_hash()` | Access `SpatialHash` |

---

### `SpatialHash`
**Header:** `src/core/spatial_hash.h`

Uniform grid spatial index for O(N) neighbor queries.

| Method | Description |
|--------|-------------|
| `clear()` | Reset all buckets |
| `insert(Int idx, Vec3 pos)` | Insert agent index at position |
| `query_neighbors(Vec3 pos)` | Return agent indices in adjacent cells |

Cell size default: 10 um (configurable via `hash_cell_size`).

---

## Field Classes

### `ChemicalField`
**Header:** `src/fields/chemical_field.h`

Grid-based storage for chemical species concentrations and reaction rates.

| Method | Description |
|--------|-------------|
| `init(Domain&, vector<ChemicalSpec>&)` | Initialize species grid |
| `find(string name)` | Species index by name (-1 if absent) |
| `conc(species, cell)` | Concentration (read/write) |
| `reac(species, cell)` | Reaction rate (read/write) |
| `zero_reactions()` | Reset all reaction rates to zero |
| `apply_boundaries(Domain&)` | Apply boundary conditions |

Default species: `carbon`, `iron`, `b12`, `bacteriocin`.

---

### `AdvectionField`
**Header:** `src/fields/advection.h`

Dual-vector mucus flow field.

| Method | Description |
|--------|-------------|
| `init(AdvectionConfig&, Domain&)` | Compute max velocities |
| `velocity(Vec3 pos)` | Flow velocity at position |
| `radial_velocity(Real z)` | Radial component at height z |
| `distal_velocity(Real z)` | Distal component at height z |
| `shear_rate(Vec3 pos)` | Local shear magnitude |
| `advect(Vec3& pos, Real dt)` | Move a position by flow |
| `washout_rate(Real z)` | Dilution rate at height z |
| `taylor_aris_D_eff(z, D_mol)` | Shear-enhanced dispersion coefficient |

---

### `VBF`
**Header:** `src/fields/vbf.h`

Viscoelastic Background Field (anaerobic microbiota continuum).

| Method | Description |
|--------|-------------|
| `init(VBFConfig&, Domain&)` | Initialize background field |
| `apply_nutrient_coupling(chem, dt)` | Sink/source on chemical field |
| `drag_force(Vec3 vel)` | Stokes-like drag on agent velocity |
| `viscosity()` | Effective viscosity (Pa·s) |

---

## Diffusion

### `GreensFunction`
**Header:** `src/diffusion/greens_function.h`

Analytical Green's function kernels for point sources in bounded domain.

| Method | Description |
|--------|-------------|
| `steady_state(r, params)` | C = Q/(4πD_eff r) for point source |
| `with_advection(r, pos, params)` | Includes flow-advected asymmetry |
| `with_images(r, z_src, z_tgt, params)` | Method of Images for bounded z |

### `QSSASolver`
**Header:** `src/diffusion/qssa_solver.h`

Quasi-steady-state toxin field computation.

| Method | Description |
|--------|-------------|
| `init(QSSAConfig&, Domain&, AdvectionField&)` | Setup |
| `solve_bacteriocin_field(agents, chem, idx)` | Compute toxin field |
| `solve_nutrient_depletion(agents, chem)` | Compute depletion zones |
| `point_concentration(target, sources, params)` | Single-point query |

---

## Fix Modules

All fixes inherit from `Fix` (header: `src/fixes/fix.h`):

```cpp
class Fix {
  virtual void init() {}
  virtual void pre_step(Real dt) {}
  virtual void compute(Real dt) = 0;
  virtual void post_step(Real dt) {}
};
```

| Fix | Config struct | Key parameters |
|-----|--------------|----------------|
| `FixMetabolism` | `MetabolismConfig` | `mu_max_default`, `division_threshold`, `metE_penalty`, `eut_penalty`, `maintenance_rate` |
| `FixBacteriocin` | `BacteriocinConfig` | `sos_basal_rate`, `retardation_basic/acidic/neutral`, `D_free_colicin`, `burst_molecules`, `microcin_mu_penalty` |
| `FixReceptor` | `ReceptorConfig` | `kd_*` binding affinities, `kill_rate_colicin/microcin`, `immunity_factor` |
| `FixConjugation` | `ConjugationConfig` | `base_transfer_prob`, `contact_radius`, `shear_crit` |
| `FixMutation` | `MutationConfig` | Per-division rates for duplication, recombination, receptor downreg, super-killer, compensatory |

---

## Genome

### `Genome`
**Header:** `src/core/types.h`

| Field | Type | Description |
|-------|------|-------------|
| `lineage_id` | `TagID` | Ancestral lineage tag |
| `parent_id` | `TagID` | Immediate parent |
| `generation` | `uint32_t` | Division count from ancestor |
| `bi_loci` | `vector<BICluster>` | Bacteriocin-immunity clusters |
| `receptor_expression` | `array<Real, 7>` | Per-receptor expression |
| `has_conjugative_plasmid` | `bool` | Can initiate HGT |
| `mutations` | `uint32_t` | Accumulated mutation count |
| `plasmid_cost_amelioration` | `Real` | Compensatory reduction |

### `BICluster`
**Header:** `src/core/types.h`

| Field | Description |
|-------|-------------|
| `toxin_id` | Bacteriocin identity |
| `immunity_id` | Cognate immunity protein |
| `target` | Which receptor the toxin hijacks |
| `bclass` | LETHAL_CORE, LETHAL_HALO, or NEUTRAL |
| `pI` | Isoelectric point |
| `diff_coeff` | Free diffusion coefficient (m²/s) |
| `retardation` | Mucin retardation factor |
| `molecular_weight` | Daltons |

### `PlasmidLibrary`
**Header:** `src/genome/plasmid.h`

Pre-characterized colicin/microcin entries:

| Entry | pI | Target | Class | MW (kDa) |
|-------|-----|--------|-------|----------|
| ColE1 | 9.2 | BtuB | Core | 57 |
| ColE2 | 6.5 | BtuB | Halo | 62 (complex) |
| ColB | 5.5 | FepA | Halo | 55 |
| ColIa | 5.8 | CirA | Halo | 69 |
| ColM | 4.8 | FhuA | Halo | 29 |
| MccV | 8.8 | CirA | Core | 8.8 |

---

## I/O

### `HDF5Writer`
**Header:** `src/io/hdf5_writer.h`

Writes simulation state in HDF5 format compatible with nufeb_tools.

```
/step_NNNN/
  atoms/        → id, type, x, y, z, radius, biomass, mu, state, lineage
  grid/         → carbon, iron, b12, bacteriocin
  metadata/     → time, num_agents, num_lineages
```

Set `HDF5Config.enabled = false` to suppress I/O (useful for testing).

### `LineageTracker`
**Header:** `src/genome/lineage_tracker.h`

Tracks genealogy and allele time-series for validation.

| Method | Description |
|--------|-------------|
| `record_mutation(tag, type, lineage)` | Log mutation event |
| `record_lysis(tag, pos, lineage)` | Log SOS lysis |
| `record_washout(tag, lineage, pos)` | Log washout event |
| `take_snapshot(time, lineages)` | Periodic population census |
| `resident_retention(since_time)` | Fraction of original lineages surviving |
| `dominant_lineage()` | Most abundant lineage tag |

---

## Python Toolkit

Package: `python/gut_ibm_tools/`

### `hdf5_reader`
Read GutIBM HDF5 output files into pandas DataFrames.

### `analysis`
Spatial statistics and exclusion-radius clustering metrics.

| Function | Description |
|----------|-------------|
| `nearest_neighbor_distances(positions, types)` | NND between competing clones (nearest different-type agent). Returns `dict[int, ndarray]`. |
| `inter_type_distances(positions, types)` | Pairwise inter-type NND. Returns `dict[tuple[int,int], ndarray]`. |
| `spatial_clustering_index(positions, types)` | Hopkins statistic variant per type. |
| `monochromatic_patch_score(positions, types, radius)` | Fraction of same-type neighbors within radius. |
| `exclusion_radius(positions, types, target_type)` | Mean distance from target-type agents to nearest different-type agent. |
| `hopkins_statistic(positions, n_samples=None)` | Hopkins clustering statistic over the full point cloud. H > 0.7 → clustered. |
| `comet_tail_asymmetry_index(positions, concentrations, flow_direction=0)` | Concentration-weighted downstream elongation ratio. |
| `comet_tail_index(positions, concentrations, flow_direction)` | Downstream/upstream mean concentration ratio. |

### `validation`
Compare simulation output to empirical targets using exclusion-radius clustering (VADI §75).

`validate_spatial_signatures(data, step)` returns:
- `monochromatic_score` – same-type neighbor fraction (target > 0.7)
- `comet_tail_ratio` – advective asymmetry (target > 1.5)
- `mean_exclusion_radius` – mean inter-type boundary distance
- `hopkins_statistic` – global clustering (target > 0.7)
- `nnd_mean` – grand mean of competing-clone NND
- `comet_tail_asymmetry` – concentration-weighted elongation

`validate_genomic_signatures(data)` returns:
- `resident_retention` – fraction of original lineages surviving (target 70–80%)
- `resident_mean_bi_loci` / `transient_mean_bi_loci` – BI cluster counts
- `transient_mean_btuB_expression` – receptor downregulation in transients

### `visualization`
Plot agent distributions, chemical fields, lineage time-series.
