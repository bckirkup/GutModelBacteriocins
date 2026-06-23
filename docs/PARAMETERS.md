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

---

## VBF (Viscoelastic Background Field)

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `vbf.density` | 1e11 | #/m^3 | Anaerobic background density |
| `vbf.drag_coeff` | 1e-9 | N·s/m | Stokes drag coefficient |
| `vbf.nutrient_sink` | 1e-4 | mol/m^3/s | Background nutrient consumption |
| `vbf.mucin_liberation` | 5e-5 | mol/m^3/s | Monosaccharide release |
| `vbf.carrying_cap` | 1e12 | #/m^3 | Local carrying capacity |
| `vbf.viscosity` | 0.01 | Pa·s | Effective viscosity (~10× water) |

---

## QSSA Solver

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `qssa.toxin_cutoff` | 200e-6 | m | Green's function evaluation radius for toxins |
| `qssa.nutrient_cutoff` | 50e-6 | m | Cutoff for nutrient depletion zones |
| `qssa.colicin_release_rate` | 1e-18 | mol/s | Burst release per lysed cell |
| `qssa.microcin_secretion` | 1e-20 | mol/s | Continuous secretion rate |

**Scaling note:** At 10^6 agents, naive O(N × M) evaluation is expensive. The cutoff radius limits each grid cell to nearby sources only, giving effective O(N) via spatial hashing.

---

## Metabolism

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `metabolism.mu_max_default` | 5e-4 | 1/s | Default max growth rate |
| `metabolism.division_threshold` | 2.0 | — | Biomass ratio for division |
| `metabolism.metE_penalty` | 0.05 | — | MetE pathway cost (BtuB loss) |
| `metabolism.eut_penalty` | 0.03 | — | Ethanolamine loss (BtuB loss) |
| `metabolism.maintenance_rate` | 1e-5 | 1/s | Maintenance energy |

**MetE pathway:** When BtuB expression < 0.5, cells must synthesize methionine via the MetE pathway instead of the B12-dependent MetH pathway. MetE requires ~5% of the proteome, and cells also lose ethanolamine utilization (normally B12-dependent).

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
| `bacteriocin.retardation_acidic` | 1.5 | — | R for pI < 6.0 (Lethal Halo) |
| `bacteriocin.retardation_neutral` | 5.0 | — | R for 6.0–8.5 |
| `bacteriocin.D_free_colicin` | 4e-11 | m^2/s | Free diffusion (~50kDa protein) |
| `bacteriocin.burst_molecules` | 1e4 | — | Molecules per lysis burst |
| `bacteriocin.microcin_mu_penalty` | 0.03 | — | Growth cost of microcin secretion |

---

## Conjugation

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `conjugation.base_transfer_prob` | 1e-3 | — | Base HGT probability |
| `conjugation.contact_radius` | 2e-6 | m | Max pilus reach |
| `conjugation.shear_crit` | 10.0 | 1/s | Critical shear for MPS |

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
| `mutation.compensatory_reduction` | 0.005 | — | Per-locus cost reduction |
| `mutation.max_bi_loci` | 8 | — | Maximum BI clusters per genome |

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
  "plasmids": ["colicin_E1"],
  "conjugative": false
}
```
