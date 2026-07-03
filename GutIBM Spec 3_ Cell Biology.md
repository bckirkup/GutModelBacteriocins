# Devin Task: Cell Biology Expansion (CDI, Fur-Regulated Receptors, Active Motility)

## Context

GutIBM models individual E. coli cells with metabolic trade-offs, receptor
binding, and passive advection. Three important cell-level behaviors are
missing: contact-dependent killing, dynamic receptor regulation, and
active motility.

## 1. Contact-Dependent Inhibition (CDI) — new Fix module

### Biology
CDI systems (CdiA/CdiB) deliver toxic effector domains directly into
neighboring cells upon physical contact. Two classes:
- **Class I (PFT):** Pore-forming toxins — moderate killing
- **Class II (tRNase):** tRNA nucleases — produce significantly higher
  "neighbor index" (stronger spatial clustering effect)

CDI requires direct cell-cell contact — it's the close-quarters complement
to diffusible colicins. T6SS is similar but more relevant for interspecies
competition (Klebsiella vs. E. coli); CDI is intraspecific.

Key dynamics from the knowledge base:
- T6SS/CDI efficacy is self-limiting via "corpse barriers" (dead cells
  accumulate at the interface and block further contact)
- Survival against CDI assault depends on the invader:resident ratio at
  the microcolony boundary
- Internal growth outpaces contact-interface death if invaders form
  microcolonies first

CDI immunity is strain-specific (like colicin immunity).

### Implementation

Create `src/fixes/fix_cdi.h` and `src/fixes/fix_cdi.cpp`.

**Agent state:** Add a `cdi_type` (uint16_t) and `cdi_immunity` (uint16_t)
to the Agent genome. CDI-positive agents deliver their effector to neighbors;
CDI-immune agents with matching immunity are protected.

**Logic (per timestep):**
```
for each agent A with cdi_type > 0:
    for each neighbor B within contact_radius:
        if B.cdi_immunity != A.cdi_type:
            P(kill) = 1 - exp(-k_cdi * dt)
            if kill:
                B.state = DEAD
                // Optionally: B becomes a "corpse" that blocks further CDI
                // by occupying space (handled by mechanics — dead cells
                // aren't removed immediately)
```

**Corpse barrier:** Dead cells from CDI are not immediately removed. They
persist as inert obstacles for `corpse_persistence` seconds (default 300 s),
during which they occupy space and block CDI contact with cells behind them.
After persistence time, they are removed. This self-limits CDI killing at
colony interfaces.

### Config Parameters

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `cdi.enabled` | false | — | Enable CDI fix |
| `cdi.kill_rate` | 5e-4 | 1/s | CDI killing rate per contact pair |
| `cdi.contact_radius` | 1e-6 | m | Max distance for CDI delivery |
| `cdi.corpse_persistence` | 300 | s | How long dead cells block further CDI |

**Initial strains:** Add `cdi_type` and `cdi_immunity` to strain definition:
```json
{"type": 1, "count": 50, "mu_max": 5.5e-4, "plasmids": ["ColE1"],
 "cdi_type": 1, "cdi_immunity": 1, "conjugative": false}
```

### Files to Create/Modify
- `src/fixes/fix_cdi.h` — new
- `src/fixes/fix_cdi.cpp` — new
- `src/fixes/fix_registry.cpp` — register "cdi" fix
- `src/core/types.h` — add cdi_type, cdi_immunity to Agent genome
- `src/io/input_parser.cpp` — parse cdi.* and per-strain cdi fields

---

## 2. Fur-Regulated Dynamic Receptor Expression

### Biology
The Ferric Uptake Regulator (Fur) is the master regulator of iron homeostasis
in E. coli. When iron is abundant, Fur represses iron-uptake genes (fepA,
fhuA, iroN, iutA, fiu, cirA) and iron-responsive microcin genes (MccH47,
MccM, MccV). When iron is scarce, Fur derepression causes massive
upregulation of these receptors.

This is the "Vulnerability Paradox": iron limitation forces cells to
upregulate the very receptors that colicins hijack. A resident colony that
depletes local iron creates a zone where immigrants are maximally vulnerable
to colicins.

Currently, receptor expression levels (`receptor_expr[]`) are fixed per agent
(initialized at 1.0, only changed by mutation). They should dynamically
respond to local iron concentration.

### Implementation

In `fix_metabolism.cpp`, after computing the iron Monod term, add
Fur-regulated receptor scaling:

```cpp
// Fur regulation: low iron → high receptor expression
// High iron → Fur represses receptors to basal level
Real iron_local = chemical_field.conc("iron", agent.grid_cell);

for (int r = 0; r < NUM_RECEPTORS; ++r) {
    if (!is_iron_receptor(r)) continue;  // only iron-related receptors

    Real base_expr = agent.receptor_expr_base[r];  // genetic baseline (mutations affect this)
    Real fur_factor = 1.0 + fur_upregulation_max * Km_fur / (Km_fur + iron_local);
    // When iron >> Km_fur: fur_factor → 1.0 (basal expression)
    // When iron << Km_fur: fur_factor → 1.0 + fur_upregulation_max (~5×)

    agent.receptor_expr[r] = std::min(base_expr * fur_factor, receptor_max);
}
```

This means:
- In iron-replete zones (near siderophore-producing residents): receptors
  are at basal level, low colicin susceptibility
- In iron-depleted zones (the "nutritional shell"): receptors are upregulated
  5×, dramatically increasing colicin susceptibility

**New agent fields:**
- `receptor_expr_base[NUM_RECEPTORS]` — genetic baseline (affected by mutations)
- `receptor_expr[NUM_RECEPTORS]` — effective expression (base × Fur regulation)

The existing `receptor_expr` becomes the effective (dynamic) value. The
mutation system in `fix_mutation.cpp` should modify `receptor_expr_base`
instead of `receptor_expr`.

### Config Parameters

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `fur.enabled` | false | — | Enable Fur-regulated receptor dynamics |
| `fur.Km` | 1e-5 | mol/m³ | Iron concentration for half-max Fur repression |
| `fur.upregulation_max` | 4.0 | — | Max fold-upregulation under iron starvation |
| `fur.receptor_max` | 5.0 | — | Cap on effective receptor expression |

### Iron field requirements
This assumes iron is tracked as a chemical species with local depletion by
siderophore-producing agents. The existing iron Monod term provides the
concentration. If iron is not yet a dynamic field (just a parameter), it
may need to be promoted to a full chemical species — or at minimum, local
iron depletion zones around siderophore producers should be computed by
the QSSA solver (analogous to nutrient depletion zones, `nutrient_cutoff`).

### Files to Modify
- `src/core/types.h` — add `receptor_expr_base[]` to Agent
- `src/fixes/fix_metabolism.cpp` — Fur regulation logic
- `src/fixes/fix_mutation.cpp` — modify `receptor_expr_base` not `receptor_expr`
- `src/io/input_parser.cpp` — parse fur.* parameters
- Agent initialization — set `receptor_expr_base` = 1.0 for all receptors

---

## 3. Active Cell Motility

### Biology
E. coli in the gut is not passively carried by mucus flow. It actively swims
using flagella in a distinctive pattern:
- Average speed in mucus: 7.76 μm/s (reduced to 5.63 with lipids)
- "Alternating back-and-forth run phases interrupted by short stops without
  reorientation or tumbling" (NOT the classic run-and-tumble of liquid culture)
- Track linearity: 18.3% (highly tortuous)
- Flagellar activity is activated 9.3× when associated with mucin-producing cells
- In chemotactic clusters: center cells suppress tumbling, boundary cells
  tumble more → stable macroscopic cluster

Motility affects:
- Colonization: motile cells can swim toward favorable niches (O₂, carbon)
- Escape from colicin zones: cells could potentially swim out of comet-tails
- Cluster formation: chemotactic aggregation creates microcolonies
- Washout resistance: swimming against radial flow (if fast enough)

### Implementation

Add active velocity to agent physics. Each timestep, agent velocity is:
```
v_total = v_advection(z) + v_swim
```

**Swimming model (mucus-adapted run-and-reverse):**
```cpp
struct MotilityState {
    Vec3 swim_direction;    // unit vector
    Real swim_speed;        // current speed (μm/s)
    Real run_timer;         // time until next reorientation
    bool is_stopped;        // in a stop phase
    Real stop_timer;        // time until resume
};
```

Each timestep:
1. Decrement `run_timer` by dt
2. If `run_timer <= 0`: reverse direction (180° flip) or stop
   - P(stop) = 0.3 per reorientation event
   - P(reverse) = 0.7
   - New run duration: exponential(mean = run_mean_duration)
3. If stopped: decrement `stop_timer`, resume when expired
4. Apply `v_swim = swim_speed * swim_direction` to position update

**Chemotaxis (optional, higher complexity):**
Bias run duration based on temporal gradient of attractant (O₂ or carbon):
```
run_duration_modifier = 1 + chi * d[attractant]/dt
```
Where `chi` is the chemotactic sensitivity. Positive gradients (swimming
toward attractant) extend runs; negative gradients shorten them.

**Cluster behavior:**
When agent density within `cluster_radius` exceeds a threshold:
- Reduce tumble/reversal rate (suppress reorientation)
- This creates the observed "suppressed tumbling center" behavior
- Effectively traps cells in dense clusters

### Config Parameters

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `motility.enabled` | false | — | Enable active swimming |
| `motility.swim_speed` | 7.76e-6 | m/s | Mean swimming speed in mucus |
| `motility.run_mean_duration` | 1.0 | s | Mean run duration |
| `motility.stop_probability` | 0.3 | — | P(stop) per reorientation |
| `motility.stop_duration` | 0.5 | s | Mean stop duration |
| `motility.chemotaxis_enabled` | false | — | Enable chemotactic bias |
| `motility.chi_carbon` | 0.1 | — | Chemotactic sensitivity to carbon |
| `motility.chi_oxygen` | 0.1 | — | Chemotactic sensitivity to O₂ |
| `motility.cluster_suppress_radius` | 10e-6 | m | Radius for cluster detection |
| `motility.cluster_suppress_threshold` | 5 | — | Neighbor count to suppress tumbling |
| `motility.cluster_tumble_factor` | 0.2 | — | Tumble rate multiplier in cluster center |

### Files to Create/Modify
- `src/core/types.h` — add MotilityState to Agent
- `src/fixes/fix_motility.h` — new
- `src/fixes/fix_motility.cpp` — new (run-and-reverse + chemotaxis)
- `src/fixes/fix_registry.cpp` — register "motility" fix
- `src/core/simulation.cpp` — initialize motility state
- `src/io/input_parser.cpp` — parse motility.* parameters

---

## Priority Order

1. **Fur-regulated receptors** — directly implements the Vulnerability Paradox,
   modifies existing code only, high biological impact
2. **CDI** — new Fix module, relatively self-contained, important for
   microcolony-scale competition
3. **Active motility** — largest change, affects all spatial dynamics, but
   can be enabled incrementally (swimming first, chemotaxis later, cluster
   behavior last)

## Tests to Add

- `test_fur_upregulation`: Low iron → receptor expression increases
- `test_fur_repression`: High iron → receptor at baseline
- `test_fur_increases_susceptibility`: Low iron + colicin → higher kill rate
  than high iron + colicin
- `test_cdi_kills_neighbor`: CDI+ agent kills adjacent CDI-susceptible agent
- `test_cdi_immunity_protects`: CDI+ agent does NOT kill immune neighbor
- `test_cdi_corpse_barrier`: Accumulated dead cells reduce CDI efficiency
- `test_motility_displacement`: Motile agent moves further than non-motile
  over N timesteps
- `test_chemotaxis_bias`: Agent with chemotaxis moves preferentially toward
  attractant source

## Biological References

- CDI neighbor index: Class II-tRNase > Class I-PFT (knowledge_base §Bacteriocin Mediated Interference)
- T6SS corpse barriers, ratio-dependent survival (knowledge_base §Bacteriocin Mediated Interference)
- Fur regulation of microcins MccH47, MccM, MccV (knowledge_base §Regulatory Triggers)
- E. coli swim speed 7.76 μm/s in mucus (knowledge_base §Spatial Organization)
- Run-and-reverse (not run-and-tumble) in mucus (knowledge_base §Spatial Organization)
- Track linearity 18.3% (knowledge_base §Spatial Organization)
- Flagellar activation 9.3× near mucin cells (knowledge_base §Micro-Scale Biogeography)
- Chemotactic cluster: center suppresses tumbling, boundary tumbles high
  (knowledge_base §Micro-Scale Biogeography)
