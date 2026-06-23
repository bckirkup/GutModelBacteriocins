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
| `global_agent_count()` | Total agents across all MPI ranks |
| `global_mu_avg()` | Global mean growth rate (via MPI_Allreduce) |

**Timestep structure:**
```
step(dt):
  1. Pre-step: clear ghosts, zero reactions
  2. Exchange ghost agents (MPI boundary layers)
  3. Rebuild spatial hash, update grid coupling
  4. module_biology(dt)   — execute all Fix modules (uses ghosts for neighbor queries)
  5. Clear ghost agents
  6. module_chemistry()   — QSSA steady-state fields
  7. module_physics(dt)   — crypt migration, advection, drag, repulsion
  8. Post-step: fix post-processing (division, lysis completion)
  9. Migrate agents that crossed slab boundaries (MPI_Sendrecv)
 10. Cleanup: check washout (crypt agents exempt), remove dead agents
 11. MPI_Allreduce for global statistics
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
| `receptor_expr` | `array<Real, 8>` | Expression level per receptor (0–1) |
| `state` | `PhenoState` | NORMAL, RESISTANT, SOS_INDUCED, DEAD |
| `genome` | `Genome` | BI clusters, lineage, mutations |
| `age` | `Real` | Time since last division (s) |
| `in_crypt` | `bool` | True when agent resides in a crypt refugium |

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
| `init(DomainConfig&)` | Set up grid, spatial hash, and slab decomposition |
| `lo()`, `hi()` | Global domain bounds (Vec3) |
| `nx()`, `ny()`, `nz()` | Grid dimensions |
| `ncells()` | Total grid cells |
| `apply_pbc(Vec3&)` | Apply periodic boundary conditions |
| `pos_to_grid(Vec3, ix, iy, iz)` | Position → grid indices |
| `cell_index(ix, iy, iz)` | Grid indices → flat index |
| `min_image_delta(a, b)` | Minimum-image displacement |
| `spatial_hash()` | Access `SpatialHash` |
| `local_lo_x()`, `local_hi_x()` | This rank's slab bounds along decomposition axis |
| `ghost_width()` | Ghost layer thickness (m) |
| `is_local(Vec3)` | Check if position is within this rank's slab |
| `owner_rank(Vec3)` | Determine which rank owns a position |
| `rank_lo()`, `rank_hi()` | Neighbor ranks in ±x direction (-1 if none) |

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
| `washout_rate(Real z)` | Dilution rate at height z (0 in crypt zone) |
| `taylor_aris_D_eff(z, D_mol)` | Shear-enhanced dispersion coefficient |
| `in_crypt_zone(Real z)` | True if z falls within the crypt zero-flow zone |

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

### `Octree`
**Header:** `src/diffusion/octree.h`

Barnes-Hut octree for O(N log N) far-field source aggregation. Distant source clusters are approximated by a single monopole at the source-weighted centroid. Controlled by the opening-angle parameter `theta`.

| Method | Description |
|--------|-------------|
| `build(positions, strengths, Domain&)` | Construct octree from source positions and strengths |
| `evaluate_far_field(target, theta, near_cutoff, gf, avg_params)` | Far-field monopole contribution at a target point |
| `evaluate_field(target, theta, near_cutoff, gf, positions, params, avg_params)` | Full near+far field at target |
| `num_nodes()` | Number of tree nodes |
| `empty()` | Whether the tree has zero nodes |

**`OctreeNode` fields:**

| Field | Type | Description |
|-------|------|-------------|
| `center` | `Vec3` | Geometric center of the cell |
| `half_size` | `Real` | Half the side length |
| `total_source_strength` | `Real` | Sum of source rates in subtree |
| `center_of_source` | `Vec3` | Source-weighted centroid |
| `children[8]` | `int` | Child node indices (-1 = no child) |
| `sources` | `vector<int>` | Leaf source indices |
| `is_leaf` | `bool` | Whether this is a leaf node |

### `QSSASolver`
**Header:** `src/diffusion/qssa_solver.h`

Quasi-steady-state toxin field computation. Supports optional Barnes-Hut acceleration via `QSSAConfig::use_fmm`.

| Method | Description |
|--------|-------------|
| `init(QSSAConfig&, Domain&, AdvectionField&)` | Setup |
| `solve_bacteriocin_field(agents, chem, idx)` | Compute toxin field (exact or Barnes-Hut depending on `use_fmm`) |
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
| `FixMetabolism` | `MetabolismConfig` | `mu_max_default`, `division_threshold`, `metE_penalty`, `eut_penalty`, `maintenance_rate`, `km_iron_primary`, `km_iron_iroN`, `km_iron_iutA`, `km_iron_fiu` |
| `FixBacteriocin` | `BacteriocinConfig` | `sos_basal_rate`, `retardation_basic/acidic/neutral`, `D_free_colicin`, `burst_molecules`, `microcin_mu_penalty` |
| `FixReceptor` | `ReceptorConfig` | `kd_*` binding affinities, `kill_rate_colicin/microcin`, `immunity_factor` |
| `FixConjugation` | `ConjugationConfig` | `base_transfer_prob`, `contact_radius`, `shear_crit` |
| `FixMutation` | `MutationConfig` | Per-division rates for duplication, recombination, receptor downreg, super-killer, compensatory; `immunity_escape_prob`, `escape_affinity_lo/hi` |
| `FixMechanics` | `MechanicsConfig` | `hertz_k`, `hertzian_enabled`, `adhesion_enabled`, `adhesion_strength`, `adhesion_range` |

---

### `FixMechanics`
**Header:** `src/fixes/fix_mechanics.h`

Computes Hertzian contact forces between overlapping agents and optional EPS-mediated adhesion. Replaces the previous inline linear repulsion in `module_physics()`.

```cpp
FixMechanics(Simulation& sim, const MechanicsConfig& cfg);
void compute(Real dt) override;
```

**Hertzian repulsion:** For each neighbor pair with overlap > 0, applies `F = hertz_k * overlap^1.5` as equal-and-opposite impulses weighted by reduced mass.

**EPS adhesion:** When enabled, cells separated by gap < `adhesion_range` receive a linearly decaying attractive force `F = adhesion_strength * (1 - gap/range)`.

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
| `receptor_expression` | `array<Real, 8>` | Per-receptor expression |
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
| `immunity_binding_affinity` | Cross-protection efficacy (1.0 = full, 0.0 = none) |

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

## OpenMP Parallelism

GutIBM supports shared-memory parallelism via OpenMP. Enable with `-DGUTIBM_USE_OPENMP=ON`.

**Compile-time guard:** All OpenMP pragmas are wrapped in `#ifdef GUTIBM_OPENMP`, so the code compiles cleanly with or without OpenMP.

**Parallelized functions:**

| Function | File | Strategy |
|----------|------|----------|
| `GreensFunction::superpose_to_grid()` | `greens_function.cpp` | Thread-local grid vectors, merged via `omp critical` |
| `QSSASolver::solve_bacteriocin_field()` | `qssa_solver.cpp` | `omp parallel for` over grid cells (deposition) |
| `QSSASolver::solve_nutrient_depletion()` | `qssa_solver.cpp` | `omp parallel for` with `omp atomic` on shared grid |
| `FixMetabolism::compute()` | `fix_metabolism.cpp` | `omp parallel for` over agents; atomic nutrient writes |
| `FixReceptor::compute()` | `fix_receptor.cpp` | Parallel kill-prob precomputation; serial RNG application |
| `Simulation::module_physics()` | `simulation.cpp` | `omp parallel for` over advection pass |
| `Simulation::module_chemistry()` | `simulation.cpp` | `omp parallel for` over reaction application |
| `Simulation::update_grid_coupling()` | `simulation.cpp` | `omp parallel for` over agents |

**Thread safety notes:**
- Chemical field reaction writes use `#pragma omp atomic` for concurrent cell access
- `GreensFunction::superpose_to_grid()` uses thread-local accumulation buffers to avoid false sharing
- RNG-dependent operations (kill decisions, mutation) remain serial
- Mechanical repulsion (pairwise neighbor push) remains serial due to symmetric position updates

---

## GPU Acceleration (Planned)

A stub for future CUDA/OpenCL kernels is located in `src/gpu/`. See `src/gpu/README.md` for the planned architecture covering Green's function evaluation, spatial hash operations, agent updates, and field operations.

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
