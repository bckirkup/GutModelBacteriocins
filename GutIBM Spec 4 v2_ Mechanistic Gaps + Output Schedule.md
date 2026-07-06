# Devin Task: Spec 4 — Mechanistic Gaps + HDF5 Output Controls

## Context

A code audit identified three mechanistic gaps that compromise the biological
fidelity of the simulation, plus an output sizing problem that makes long runs
impractical. This spec addresses all four.

---

## 1. Per-Receptor Toxin Fields (CRITICAL)

### Problem
All colicins and microcins contribute to a single "bacteriocin" chemical
species. The kill probability code in `fix_receptor.cpp` discriminates by
receptor type (BtuB, FepA, CirA), but reads the same `tox_conc` for all three.
When a cell carrying ColE1 (BtuB-targeted) and ColB (FepA-targeted) lyses,
both toxins go into one field. A susceptible cell then sees the combined
concentration at ALL receptor types, overpredicting kill probability.

This is mechanistically wrong: a BtuB-targeted colicin can't enter through
FepA. The kill probability at each receptor should depend on the concentration
of toxins that actually target that receptor.

### Fix
Replace the single "bacteriocin" species with per-receptor toxin fields:

```
bacteriocin_BtuB   — colicins E1, E2, E5, E7, E8, E9, K
bacteriocin_FepA   — colicins B, D
bacteriocin_CirA   — colicin Ia, microcin V
bacteriocin_FhuA   — colicin M
```

Register these as separate chemical species in `ChemicalField`. Use
`species_names.h` constants:
```cpp
constexpr const char* BACTERIOCIN_BTUB = "bacteriocin_BtuB";
constexpr const char* BACTERIOCIN_FEPA = "bacteriocin_FepA";
constexpr const char* BACTERIOCIN_CIRA = "bacteriocin_CirA";
constexpr const char* BACTERIOCIN_FHUA = "bacteriocin_FhuA";
```

#### Changes to fix_bacteriocin.cpp
In `lyse_agent()`, deposit each BI cluster's burst into the field matching
its `target` receptor:
```cpp
Int species_idx;
switch (bi.target) {
    case ReceptorType::BtuB: species_idx = chem.find(species::BACTERIOCIN_BTUB); break;
    case ReceptorType::FepA: species_idx = chem.find(species::BACTERIOCIN_FEPA); break;
    case ReceptorType::CirA: species_idx = chem.find(species::BACTERIOCIN_CIRA); break;
    case ReceptorType::FhuA: species_idx = chem.find(species::BACTERIOCIN_FHUA); break;
    default: continue;
}
```

Microcin continuous secretion similarly deposits into the receptor-specific field.

#### Changes to fix_receptor.cpp
In `compute_kill_prob()`, each receptor block reads its own toxin field:
```cpp
// BtuB block:
Int i_tox_btuB = chem.find(species::BACTERIOCIN_BTUB);
Real tox_conc_btuB = (i_tox_btuB >= 0) ? chem.conc(i_tox_btuB, cell) : 0.0;
Real occ = toxin_occupancy(tox_conc_btuB, ligand, ...);
```

Similarly for FepA and CirA blocks.

Remove the single `bacteriocin` species lookup at the top of the function.
If all four per-receptor fields are zero, return 0.0 early.

#### Changes to qssa_solver.cpp
The QSSA solver must handle multiple toxin species. Each burst source carries
its target receptor type. The Green's function evaluation is per-species:
```
for each burst source:
    deposit into grid for species matching source.target_receptor
```

The pI-dependent retardation, protease decay, and advection-diffusion kernel
all operate per-source and are already parameterized by the BICluster fields,
so no change to the kernel math is needed — just routing to the correct
grid species.

#### Changes to cross-induction (fix_bacteriocin.cpp)
`local_nuclease_toxin()` in simulation.cpp should read from `bacteriocin_BtuB`
specifically (nuclease colicins E2/E7/E8/E9 all target BtuB), not the old
combined field.

#### Memory impact
4 toxin species × 800,000 grid cells × 8 bytes × 2 arrays = 51 MB additional
VRAM. Acceptable on the 3070 Ti.

---

## 2. Siderophore Production and Iron Depletion Zones

### Biology
E. coli secretes enterobactin (a catecholate siderophore) to scavenge Fe³⁺.
Enterobactin has an extremely high affinity for iron (Kf ~10⁴⁹). The
siderophore-iron complex is reimported via FepA.

This creates two important ecological effects:
1. **Local iron enrichment** around the producer (siderophore cloud)
2. **Iron depletion** in surrounding regions (siderophore strips iron from
   the environment faster than abiotic diffusion replenishes it)

The depletion zone is the "nutritional shell" from the ISIM framework: a
region around a resident colony where iron is depleted, forcing immigrants
to upregulate FepA/IroN/FhuA (via Fur derepression), which increases their
vulnerability to colicins.

Currently, iron is consumed by Monod uptake but there is no siderophore
secretion/chelation/reimport cycle.

### Implementation

Add `siderophore` as a new chemical species.

**Source:** Each agent secretes siderophore at a rate proportional to Fur
activation (low iron → high siderophore production):
```cpp
Real fur_activity = 1.0 - S_iron / (fur_cfg.Km + S_iron);  // 0 at high iron, 1 at low
Real sid_rate = cfg_.siderophore_secretion * fur_activity * biomass;
chem.reac(i_sid, cell) += sid_rate;
```

**Iron chelation:** Siderophore + free iron → siderophore-iron complex.
Model as a second-order reaction consuming both species:
```cpp
Real chelation = k_chelation * [siderophore] * [iron];
chem.reac(i_iron, cell) -= chelation;
chem.reac(i_sid, cell)  -= chelation;
```

**Reimport:** Agent consumes siderophore-iron complex via FepA:
```cpp
Real reimport = expr_FepA * [siderophore-iron] / (Km_reimport + [siderophore-iron]);
```

For simplicity, we can skip tracking the siderophore-iron complex as a
separate species. Instead, model the net effect: siderophore secretion
creates a zone of enhanced iron uptake near the producer and iron depletion
further away.

**Simplified implementation:**
```cpp
// In fix_metabolism, after iron Monod calculation:
if (cfg_.siderophore_enabled) {
    // Producer enhances local iron: effective [Fe] is higher near siderophore sources
    // Immigrant in a depletion zone: effective [Fe] is lower

    // Siderophore depletion radius: D_sid / (k_chelation * [Fe_background])
    // Use QSSA-style analytical depletion zone around each agent with bi_loci

    // Simple version: add iron sink proportional to agent density with BI loci
    // and iron source proportional to siderophore secretion near producer
    Real sid_secretion = fur_activity * cfg_.siderophore_rate;
    chem.reac(i_iron, cell) += sid_secretion;  // producer recaptures iron locally
    // Iron depletion in surrounding cells handled by diffusion + consumption
}
```

### Config Parameters

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `siderophore.enabled` | false | — | Enable siderophore dynamics |
| `siderophore.secretion_rate` | 1e-18 | mol/s/cell | Base enterobactin secretion |
| `siderophore.D_free` | 5e-10 | m²/s | Siderophore diffusion (~700 Da) |
| `siderophore.chelation_rate` | 1e6 | 1/(mol/m³)/s | Second-order chelation rate |
| `siderophore.Km_reimport` | 1e-7 | mol/m³ | FepA reimport half-saturation |

### Files to Modify
- `src/core/species_names.h` — add SIDEROPHORE constant
- `src/fields/chemical_field.h/.cpp` — register siderophore species
- `src/fixes/fix_metabolism.cpp` — siderophore secretion and iron reimport
- `src/io/input_parser.cpp` — parse siderophore.* parameters
- `src/io/chem_environment_config.h` — add SiderophoreConfig struct

---

## 3. CirA Ligand Correction

### Problem
Line 118 of `fix_receptor.cpp` uses `iron * 0.1` as a proxy for linearized
enterobactin concentration at CirA:
```cpp
Real ligand = (i_iron >= 0) ? chem.conc(i_iron, cell) * 0.1 : 0.0;
```

This is arbitrary. Linearized enterobactin is a degradation product of
enterobactin (produced by esterase-mediated hydrolysis). Its concentration
depends on siderophore turnover, not iron directly.

### Fix

**If siderophore species is implemented (§2):** Use siderophore concentration
as the CirA ligand, with a conversion factor for the linearized fraction:
```cpp
Int i_sid = chem.find(species::SIDEROPHORE);
Real ligand = (i_sid >= 0) ? chem.conc(i_sid, cell) * 0.3 : 0.0;
// ~30% of enterobactin pool is linearized in steady state
```

**If siderophore is not yet implemented:** At minimum, replace the arbitrary
`0.1` with a configurable parameter:
```cpp
Real ligand = (i_iron >= 0) ? chem.conc(i_iron, cell) * cfg_.cirA_iron_proxy_factor : 0.0;
```

Add `cirA_iron_proxy_factor` to `ReceptorConfig` (default 0.1).

---

## 4. HDF5 Output Schedule Matrix

### Problem
A 7-day simulation at 5μm resolution produces ~45 MB per full snapshot.
Fixed-interval dumping of everything wastes disk on data that's rarely needed
at full temporal resolution, while missing the fine-grained population dynamics
that summary stats can capture cheaply.

### Design: Configurable Output Schedule

The output schedule is a matrix of (data layer, interval) pairs. Each layer
has its own dump cadence. The user configures this as a JSON object; anything
not specified falls back to defaults.

```json
"hdf5": {
    "file": "output.h5",
    "compression": "gzip",
    "compression_level": 4,
    "schedule": {
        "summary":    1,
        "agents":     5,
        "grid":       0,
        "grid_species": [],
        "lineage":    100,
        "genome":     100
    }
}
```

| Layer | Default interval (steps) | Description | Size per dump |
|-------|------------------------:|-------------|---------------|
| `summary` | **1** (every step) | On-the-fly aggregate statistics + event counters | ~500 bytes |
| `agents` | **5** | Per-agent position, type, mu_realized, receptor_expr, state | ~80 KB at 200 agents |
| `grid` | **0** (disabled) | Full chemical field grids for selected species | ~45 MB (all) or ~6 MB (1 species) |
| `grid_species` | `[]` (none) | Which species to include in grid dumps. Empty = none. `["all"]` = everything. | — |
| `lineage` | **100** | Lineage tracker snapshot | ~10 KB |
| `genome` | **100** | Full genome/BI-loci dump per agent | ~50 KB |

**Interval = 0 means disabled.** Interval = 1 means every step. Interval = N
means every N-th step.

The `hdf5.every` legacy key is still supported: if present and `schedule` is
absent, it sets all layers to that interval (backward compatible).

### Summary Layer (default: every step)

The cheapest and most important output. Computed on-the-fly, never requires
reading grid data back:

```
summary/step_NNNNNN:
  time: Real
  dt: Real
  n_total: Int
  n_by_type: Int[max_types]
  n_in_crypt: Int
  n_by_state: Int[4]    // NORMAL, RESISTANT, SOS_INDUCED, DEAD_THIS_STEP
  mean_z_by_type: Real[max_types]
  std_z_by_type: Real[max_types]
  mean_mu_by_type: Real[max_types]
  mean_receptor_expr: Real[NUM_RECEPTORS]

  // Event counters (accumulated since last summary dump)
  events.sos_inductions: Int
  events.phage_inductions: Int
  events.colicin_kills: Int
  events.cdi_kills: Int
  events.washout_deaths: Int
  events.boundary_deaths: Int     // pushed past z_max
  events.starvation_deaths: Int   // biomass below threshold
  events.divisions: Int
  events.conjugation_transfers: Int
  events.mutations: Int

  // Chemical field summaries (no grid needed)
  chem.mean_carbon: Real
  chem.mean_iron: Real
  chem.mean_oxygen: Real
  chem.max_toxin_BtuB: Real
  chem.max_toxin_FepA: Real
  chem.max_toxin_CirA: Real
  chem.max_toxin_FhuA: Real

  // Spatial statistics (computed only if n_total >= 10)
  spatial.hopkins_statistic: Real
  spatial.mean_nnd: Real          // nearest-neighbor distance
  spatial.monochromatic_score: Real
```

### Event Counters

Add to `Simulation`:
```cpp
struct StepEvents {
    Int sos_inductions = 0;
    Int phage_inductions = 0;
    Int colicin_kills = 0;
    Int cdi_kills = 0;
    Int washout_deaths = 0;
    Int boundary_deaths = 0;
    Int starvation_deaths = 0;
    Int divisions = 0;
    Int conjugation_transfers = 0;
    Int mutations = 0;
    void reset() { *this = {}; }
};
StepEvents step_events_;
```

Each Fix increments the relevant counter during `compute()`:
- `fix_bacteriocin`: `step_events_.sos_inductions++` on SOS trigger
- `fix_receptor`: `step_events_.colicin_kills++` on kill
- `fix_cdi`: `step_events_.cdi_kills++` on kill
- `fix_metabolism`: `step_events_.divisions++` on division,
  `step_events_.starvation_deaths++` on biomass death
- `fix_conjugation`: `step_events_.conjugation_transfers++` on transfer
- `fix_mutation`: `step_events_.mutations++` on any mutation
- `check_washout()`: `step_events_.washout_deaths++` on mu < gamma death,
  `step_events_.boundary_deaths++` on z >= z_max death

The summary writer reads the accumulated counters and resets them after writing.
If `summary` interval > 1, counters accumulate across multiple steps before dump.

### Agents Layer (default: every 5 steps)

Per-agent arrays. Written as 1D datasets under `agents/step_NNNNNN/`:
- `id`, `type`, `x`, `y`, `z`, `mu_realized`, `state`, `biomass`
- `receptor_expr[NUM_RECEPTORS]` (flattened)
- `in_crypt` (bool)
- `n_bi_loci` (Int)

Does NOT include full genome (BI cluster details, affinity arrays) — that's
in the `genome` layer at lower cadence.

### Grid Layer (default: disabled)

Full 3D chemical field arrays. Only written for species listed in
`grid_species`. Each species is a 3D dataset `grid/step_NNNNNN/species_name`.

If `compression` is "gzip", apply chunk-based deflate:
```cpp
hsize_t chunk_dims[] = {static_cast<hsize_t>(nx),
                        static_cast<hsize_t>(ny),
                        static_cast<hsize_t>(nz)};
H5Pset_chunk(plist, 3, chunk_dims);
H5Pset_deflate(plist, compression_level);
```

Toxin fields are mostly zero (sparse sources) → expect 10-100× compression.

### Lineage + Genome Layers (default: every 100 steps)

Written at low cadence for post-hoc phylogenetic reconstruction:
- `lineage/step_NNNNNN`: parent-child tag pairs
- `genome/step_NNNNNN`: full BI cluster arrays, toxin/ligand affinity,
  receptor expression base, mutations count, plasmid cost amelioration

### Implementation

In `hdf5_writer.cpp`, the `dump()` method checks which layers are due:
```cpp
void HDF5Writer::dump(Int step, ...) {
    if (schedule_.summary > 0 && step % schedule_.summary == 0)
        write_summary(step, ...);
    if (schedule_.agents > 0 && step % schedule_.agents == 0)
        write_agents(step, ...);
    if (schedule_.grid > 0 && step % schedule_.grid == 0)
        write_grid(step, ...);
    if (schedule_.lineage > 0 && step % schedule_.lineage == 0)
        write_lineage(step, ...);
    if (schedule_.genome > 0 && step % schedule_.genome == 0)
        write_genome(step, ...);
}
```

### Parsing

In `input_parser.cpp`, parse the `hdf5` sub-object:
```cpp
if (key == "hdf5") {
    // Parse as JSON object with schedule sub-object
    auto hdf5_obj = parse_json_object(val);
    cfg.hdf5.filename = hdf5_obj.get("file", cfg.hdf5.filename);
    cfg.hdf5.compression = hdf5_obj.get("compression", "none");
    cfg.hdf5.compression_level = hdf5_obj.get("compression_level", 4);
    if (hdf5_obj.has("schedule")) {
        auto sched = hdf5_obj["schedule"];
        cfg.hdf5.schedule.summary = sched.get("summary", 1);
        cfg.hdf5.schedule.agents = sched.get("agents", 5);
        cfg.hdf5.schedule.grid = sched.get("grid", 0);
        cfg.hdf5.schedule.lineage = sched.get("lineage", 100);
        cfg.hdf5.schedule.genome = sched.get("genome", 100);
        // grid_species parsing
    }
}
```

### Files to Modify
- `src/io/hdf5_writer.h/.cpp` — schedule-based dump, compression, summary writer
- `src/io/input_parser.h/.cpp` — parse hdf5 sub-object with schedule
- `src/core/simulation.h/.cpp` — add StepEvents, pass to writer
- All fix modules — increment event counters
- `docs/PARAMETERS.md` — document hdf5 schedule matrix

## Priority Order

1. **Per-receptor toxin fields** (§1) — mechanistically critical, affects
   every kill calculation
2. **HDF5 output controls** (§4) — practically critical, blocks long runs
3. **Siderophore production** (§2) — completes the Vulnerability Paradox
4. **CirA ligand correction** (§3) — small fix, depends on §2

## Tests to Add

- `test_per_receptor_toxin`: ColE1 burst only appears in bacteriocin_BtuB
  field, not in bacteriocin_FepA
- `test_receptor_specific_kill`: Agent with FepA=0, BtuB=1 is killed by
  BtuB-targeted toxin but not FepA-targeted toxin at same concentration
- `test_siderophore_depletion`: Iron concentration decreases near
  siderophore-producing agent
- `test_hdf5_selective_output`: Config with agents_every=10, grid_every=100
  produces more agent snapshots than grid snapshots
- `test_hdf5_compression`: Compressed file is smaller than uncompressed
- `test_summary_events`: Event counters increment correctly for each event type

## Biological References

- Enterobactin Kf ~10⁴⁹ for Fe³⁺ (Raymond et al., PNAS 2003)
- Linearized enterobactin is the CirA ligand (Buchanan et al., 2007)
- Siderophore-mediated iron depletion zones create the "nutritional shell"
  (ISIM framework, knowledge_base §Metabolic Niche Partitioning)
- Host NGAL sequesters enterobactin but not salmochelin (Flo et al., Nature 2004)
- Colicin receptor specificity: E→BtuB, B/D→FepA, Ia/Ib→CirA, M→FhuA
  (knowledge_base §Mechanisms of Bacteriocin Translocation)
