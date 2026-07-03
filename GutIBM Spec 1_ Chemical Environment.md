# Devin Task: Chemical Environment Expansion (O₂, VFAs, Protease Decay, Mucin Degradation)

## Context

GutIBM models a patch of colonic mucus with discrete E. coli agents and a
continuum anaerobic background (VBF). The chemical field currently tracks
carbon (mucin-derived monosaccharides), iron, B12, and bacteriocin
concentrations. Several critical chemical species are missing.

The existing infrastructure supports z-dependent gradients and multiple
chemical species via the `ChemicalField` class. This spec adds four new
chemical dynamics to the existing framework.

## 1. Oxygen Gradient

### Biology
Oxygen diffuses from the submucosal vasculature through the epithelium into
the mucus layer. Concentration: ~42 mmHg at the subepithelial mucosa, dropping
to 3-11 mmHg in the colonic lumen. Enterobacteriaceae are facultative
anaerobes that concentrate near the epithelium specifically because they can
use aerobic respiration there (faster growth than fermentation).

E. coli consumes O₂ via cytochrome oxidases (cydAB, cyoABCDE). The anaerobic
background (VBF) also consumes some O₂ but at much lower rates. The steep
gradient is maintained by the balance of epithelial supply and bacterial
consumption.

### Implementation

Add `oxygen` as a new chemical species in `ChemicalField`.

**Source term (epithelial supply):**
```
dO2/dt|source = D_O2 * (O2_epithelial - O2(z)) / mucus_thickness²
```
Implemented as a Dirichlet boundary at z=0: `O2(z=0) = O2_epithelial`.
The QSSA solver handles steady-state diffusion from this boundary.

**Sink terms:**
- E. coli agents: `consumption = q_O2 * mu_realized * biomass` where
  `q_O2` is the specific O₂ consumption rate (mol O₂ per mol biomass per s)
- VBF background: `vbf_O2_sink` (constant volumetric rate, small)

**Effect on growth:**
Add O₂ as a fourth Monod term in `fix_metabolism`:
```
mu = mu_max * monod_carbon * monod_iron * monod_b12 * monod_O2_boost
```
Where `monod_O2_boost` provides a growth rate multiplier (1.0-3.0×) when O₂
is available, reflecting the ~3× faster growth from aerobic vs. fermentative
metabolism. When O₂ = 0, boost = 1.0 (fermentation baseline). This is NOT a
strict requirement — E. coli grows anaerobically — but aerobic growth is faster.

```
monod_O2_boost = 1.0 + O2_boost_max * [O2] / (Km_O2 + [O2])
```

**Effect on SOS (feeds into Spec 2):**
Aerobic growth generates ROS. Provide a method `local_O2(agent)` that
fix_bacteriocin can call to scale the SOS induction rate:
```
rate_ROS = k_ROS * [O2] * mu_realized
```
This creates the feedback: fast aerobic growth near the epithelium →
more ROS → more SOS → more colicin release. But the cell is also in the
low-flow zone near the epithelium, so it's protected from washout.

### Config Parameters

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `oxygen.enabled` | false | — | Enable oxygen dynamics |
| `oxygen.epithelial_conc` | 55e-6 | mol/m³ | O₂ at epithelium (~42 mmHg) |
| `oxygen.D_free` | 2.1e-9 | m²/s | O₂ diffusion in water |
| `oxygen.Km` | 1e-6 | mol/m³ | Monod half-saturation for O₂ boost |
| `oxygen.boost_max` | 2.0 | — | Max growth boost from aerobic respiration |
| `oxygen.q_consumption` | 1e-14 | mol/s/cell | Specific O₂ consumption per agent |
| `oxygen.vbf_sink` | 1e-6 | mol/m³/s | VBF background O₂ consumption |
| `oxygen.k_ROS` | 1e2 | 1/s per mol/m³ | ROS-driven SOS rate coefficient |

### Files to Modify
- `src/fields/chemical_field.h/.cpp` — register "oxygen" species
- `src/fixes/fix_metabolism.cpp` — add O₂ boost to growth, add consumption
- `src/core/simulation.cpp` — initialize O₂ field with z-gradient
- `src/io/input_parser.cpp` — parse oxygen.* parameters
- `docs/PARAMETERS.md` — document oxygen section

---

## 2. Dynamic Volatile Fatty Acid (Acetate) Field

### Biology
The anaerobic majority (Bacteroidetes, Firmicutes) produces SCFAs as
fermentation end-products. Colonic acetate is 60-100 mM, propionate 15-25 mM,
butyrate 10-25 mM. E. coli both produces acetate (overflow metabolism when
carbon is abundant) and consumes it (acetate scavenging via acs when carbon
is scarce).

Acetate already affects the model via MetE inhibition (VADI §87), but the
concentration is currently implicit. Making it a dynamic field creates the
metabolic cross-feeding loop.

### Implementation

Add `acetate` as a new chemical species.

**Source terms:**
- VBF: `vbf_acetate_production` (constant volumetric rate, represents
  anaerobic fermentation). This should be z-dependent — higher near
  epithelium where mucin degradation releases more carbon.
- E. coli agents: `overflow_rate = k_overflow * max(0, mu_realized - mu_threshold) * biomass`
  When growing fast (high carbon), E. coli overflows acetate.

**Sink terms:**
- E. coli agents: `scavenge_rate = k_scavenge * [acetate] / (Km_scavenge + [acetate]) * biomass`
  When carbon is scarce, E. coli switches to acetate consumption.
- VBF: `vbf_acetate_consumption` (cross-feeding by other fermenters)
- Epithelial uptake: colonocytes absorb SCFAs (small sink at z=0)

**Effect on metabolism:**
Replace the static acetate concentration in the MetE penalty calculation
with the dynamic local value from the chemical field:
```
[acetate] = chemical_field.conc("acetate", agent.grid_cell)
```

**Butyrate** — lower priority. Could be added as a separate species that
feeds into epithelial O₂ consumption (butyrate β-oxidation uses O₂),
closing the O₂-SCFA loop. But acetate alone captures the main dynamics.

### Config Parameters

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `acetate.enabled` | false | — | Enable dynamic acetate |
| `acetate.D_free` | 1.2e-9 | m²/s | Acetate diffusion coefficient |
| `acetate.vbf_production` | 1e-3 | mol/m³/s | VBF acetate production |
| `acetate.vbf_consumption` | 2e-4 | mol/m³/s | VBF acetate consumption |
| `acetate.overflow_threshold` | 3e-4 | 1/s | Growth rate above which E. coli overflows |
| `acetate.overflow_rate` | 1e-15 | mol/s/cell | Acetate overflow per agent |
| `acetate.scavenge_Km` | 5.0 | mol/m³ | Half-saturation for scavenging |
| `acetate.epithelial_uptake` | 5e-4 | mol/m³/s | Epithelial absorption rate at z=0 |

---

## 3. Protease Degradation of Bacteriocins

### Biology
Colicins (30-80 kDa proteins) are degraded by intestinal proteases (trypsin,
chymotrypsin, pepsin). Class II bacteriocins and Colicins B and M are
explicitly sensitive to these proteases. The knowledge base states: "uniquely
susceptible to intestinal proteases, sharply limiting in vivo half-lives."

MccJ25 half-life in simulated intestinal fluid: ~2 hours. Full-length colicins
are expected to be shorter — 30 min to 1 hour. Microcins (lasso peptides,
post-translationally modified) are more protease-resistant.

Currently colicins in the QSSA field persist indefinitely. This overestimates
the killing radius.

### Implementation

Add a first-order decay term to the bacteriocin concentration in the QSSA solver:
```
C(r, t) = C_steady_state(r) * exp(-k_decay * t_since_release)
```

Since the QSSA computes steady-state concentrations from active sources, the
simplest implementation is to multiply each source's contribution by an
exponential decay based on the time since that source was created (lysis event):

In `qssa_solver.cpp`, when evaluating Green's function superposition:
```cpp
Real age = current_time - source.creation_time;
Real decay = std::exp(-cfg_.protease_decay_rate * age);
contribution *= decay;
```

This requires tracking `creation_time` on each toxin source (burst event).

For continuously secreted microcins, the decay applies to the steady-state
concentration directly — the QSSA already gives a steady-state; multiply by
`1 / (1 + k_decay / k_dilution)` to account for degradation in the
steady-state balance.

**Per-colicin decay rates:** Add `protease_half_life` to `BICluster`:
- Colicin E1 (pore-forming, 57 kDa): half_life = 1800 s (30 min)
- Colicin E2 (nuclease, 61.5 kDa): half_life = 1800 s
- Colicin B (pore-forming, 54.8 kDa): half_life = 900 s (explicitly protease-sensitive)
- Colicin M (murein inhibitor, 29.5 kDa): half_life = 900 s (explicitly sensitive)
- Colicin Ia (pore-forming, 69.4 kDa): half_life = 2400 s
- Microcin V (8.9 kDa, peptide): half_life = 7200 s (2 hr, lasso-like protection)

### Config Parameters

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `protease.enabled` | true | — | Enable protease degradation |
| `protease.default_half_life` | 1800 | s | Default colicin half-life (30 min) |

Per-colicin values are set in the PlasmidLibrary, not as global config.

### Files to Modify
- `src/core/types.h` — add `protease_half_life` to BICluster
- `src/genome/plasmid.cpp` — set per-colicin half-lives
- `src/diffusion/qssa_solver.cpp` — apply decay to source contributions
- `src/fixes/fix_bacteriocin.cpp` — timestamp burst sources

---

## 4. Dynamic Mucin Degradation

### Biology
E. coli cannot degrade complex mucin oligosaccharides — it relies entirely on
other species (mainly Bacteroides) for monosaccharide liberation. B.
thetaiotaomicron desulfates mucin, releasing monosaccharides and altering
mucin structure. The current VBF has a static `mucin_liberation` rate.

In reality, mucin degradation depends on:
- Local anaerobe density (VBF density)
- Available mucin substrate (which is replenished by goblet cell secretion)
- Enzyme activity (glycosidases, sulfatases, sialidases)

Under high mucus viscosity, E. coli populations are 20× lower because glycans
aren't released efficiently.

### Implementation

Make `mucin_liberation` a dynamic rate that depends on local VBF density
and a substrate availability term:

```
liberation_rate(z) = k_lib * vbf_density(z) * mucin(z) / (Km_mucin + mucin(z))
```

This requires tracking `mucin` as a chemical species:
- **Source:** Goblet cell secretion at z=0 (Dirichlet boundary or flux term)
- **Sink:** VBF-mediated degradation (produces monosaccharides = carbon source)
- **Diffusion:** Very slow (mucin polymers are ~MDa, essentially immobile
  after secretion; transport is by mucus flow, not diffusion)

The carbon source term in the chemical field then becomes:
```
dC/dt|mucin = liberation_rate(z)   // replaces static vbf_mucin_liberation
```

This closes the loop: goblet cells → mucin → anaerobe degradation → 
monosaccharides → E. coli growth.

### Config Parameters

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `mucin.enabled` | false | — | Enable dynamic mucin |
| `mucin.secretion_rate` | 1e-4 | mol/m³/s | Goblet cell secretion at z=0 |
| `mucin.Km_degradation` | 1e-3 | mol/m³ | Half-saturation for anaerobe degradation |
| `mucin.k_liberation` | 1e-4 | 1/s | Rate constant for VBF mucin degradation |

### Files to Modify
- `src/fields/chemical_field.h/.cpp` — register "mucin" species
- `src/fields/vbf.cpp` — dynamic liberation rate
- `src/core/simulation.cpp` — initialize mucin field
- `src/io/input_parser.cpp` — parse mucin.* parameters

---

## Priority Order

1. **Protease degradation** — straightforward, high impact on killing radius
2. **Oxygen gradient** — fundamental to spatial niche, uses existing infrastructure
3. **Dynamic acetate** — closes the metabolic cross-feeding loop
4. **Dynamic mucin** — lower priority, more complex, can be a follow-up

## Tests to Add

- `test_oxygen_gradient`: Verify z-dependent O₂ profile, agent consumption
- `test_O2_growth_boost`: Agents near epithelium grow faster than luminal agents
- `test_protease_decay`: Colicin concentration decreases with time after burst
- `test_acetate_mete_dynamic`: MetE penalty responds to local acetate field
  (extends existing test_acetate_mete)
- `test_mucin_liberation`: Carbon source scales with VBF density

## Biological References

- O₂: 42 mmHg subepithelial → 3-11 mmHg lumen (knowledge_base §Spatial Organization)
- Acetate: 60-100 mM colonic (knowledge_base §Metabolic Niche Partitioning)
- Protease sensitivity: Class II bacteriocins and ColB/ColM explicitly sensitive
  to trypsin, chymotrypsin (knowledge_base §Bacteriocin Mediated Interference)
- MccJ25 half-life ~2 hr in SIF (Yu et al. 2022, Frontiers in Immunology)
- Mucin degradation: E. coli lacks extracellular enzymes, relies on Bacteroides
  (knowledge_base §Metabolic Niche Partitioning)
- High mucus viscosity → 20× lower E. coli populations (knowledge_base §Host Inflammatory States)
