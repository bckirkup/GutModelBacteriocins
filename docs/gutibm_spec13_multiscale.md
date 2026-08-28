# GutIBM Spec 13: Multi-Scale Architecture

## Motivation

The single-patch IBM cannot produce stable E. coli population dynamics because
it asks one mucus volume to exhibit the emergent behavior of a coupled
metapopulation. The model has been trapped on a knife-edge: either net growth
exceeds washout everywhere (overgrowth → guard) or falls below it everywhere
(extinction). No parameter choice resolves this because the real system's
stability comes from spatial coupling across scales the model does not
represent.

The literature (Spec 12 reviews + Greter et al. 2026, Li et al. 2015) reveals
that E. coli in healthy mucus exists as small transient clonal clusters
(~3–4 µm radius) in a gel that is periodically disrupted by peristaltic
contractions (~46 Pa yield stress). Persistence is a metapopulation property:
some patches survive each disturbance and reseed cleared areas. The
bacteriocin competition question — whether colicin-mediated interference
enables coexistence — operates across this patchy landscape.

This spec defines a three-layer hybrid model that couples the existing IBM to
a patch-dynamics layer and a colonic-section layer.

---

## Architecture

### Layer 1: IBM Patch (~100–300 µm lateral, ~50–150 µm deep)

**This is the existing model.** Individual agents, reaction-diffusion of
carbon, O₂, bacteriocin, B12, acetate. Monod growth with the Spec 12
metabolic switch. Cell mechanics, motility. GPU-accelerated diffusion.

Each patch instance resolves:
- Clonal microcolony growth (doubling time ~hours)
- Bacteriocin diffusion and killing (~50–200 µm range)
- O₂ consumption and metabolic switching
- Cell-cell competition within and between nearby microcolonies
- Timescale: seconds to ~10 minutes (inter-contraction interval)

**Patch types** (distinguished by boundary conditions, not code):

| Type | O₂ at base | Carbon flux | Washout rate | Shear protection | Fraction |
|---|---|---|---|---|---|
| **Crypt refuge** | Higher (near epithelium) | Mucin secretion point | Low (sheltered geometry) | High | 5–15% |
| **Exposed proximal** | Moderate | Cross-fed from VBF | Moderate | Low | 40–50% |
| **Exposed distal** | Low (thicker mucus, lower tissue O₂) | Lower | Higher (more flow) | Low | 35–50% |

Patch type sets the Dirichlet/Neumann boundary conditions for O₂ and carbon,
the local VBF density, and the susceptibility to contraction-driven loss.
These are derived from biophysics, not from E. coli spatial data.

### Layer 2: Mucus Segment (~1–5 mm, a grid of patches)

A **network of N patches** (default N = 9–25) representing one ~mm² region of
colonic mucosa. Patches are connected by contraction-driven exchange.

Layer 2 manages:
- **Contraction events**: Periodic (every 1–5 minutes, from colonic motility
  measurements). A contraction:
  1. Selects patches to disrupt (probability based on shear exposure)
  2. Removes a fraction of agents from disrupted patches (fragmentation)
  3. Redistributes removed agents: some go to neighboring patches (local
     reseeding), some go to the luminal pool (washout)
  4. Crypt patches have reduced disruption probability
- **Inter-patch migration**: Between contractions, a low rate of agent exchange
  between adjacent patches represents diffusive spread through intact mucus
  (slow — cells are immobile in gel, but daughter cells at colony edges can
  push into adjacent territory)
- **Luminal pool**: Agents washed out of patches enter a transit compartment.
  They can reattach to patches with probability proportional to patch
  receptivity, or be lost distally. This is not a full lumen simulation —
  it's a holding compartment for the washout/recolonization cycle.

**Contraction parameterization** (from biophysics):

| Parameter | Value | Source |
|---|---|---|
| Contraction frequency | 0.2–1.0 /min | Human/mouse colonic motility |
| Yield stress of mucus | 46 ± 9 Pa | Greter et al. 2026, cecal rheology |
| Fraction of exposed patches disrupted per contraction | 0.3–0.7 | Geometry-dependent; sweep |
| Fraction of agents lost per disruption | 0.1–0.5 | Shear intensity-dependent; sweep |
| Crypt disruption probability | 0.05–0.15 | Reduced by crypt geometry |
| Luminal transit (washout) half-life | 2–8 hours | Colonic transit time |
| Reattachment probability per transit | 0.01–0.10 | Mucus adhesion; sweep |

### Layer 3: Colonic Section (~cm to whole colon)

Three **colonic regions** (not six segments — the data only distinguish three
evidence classes; Ahmed et al. 2007, doi:10.1128/aem.01143-07), each
containing one Layer 2 mucus-segment simulation plus a **luminal growth
compartment**.

The luminal compartment is essential: mucosal detachment alone produces only
~10⁷ CFU/g stool, but observed stool E. coli is ~10⁸ CFU/g. Luminal
replication during transit accounts for the ~10× amplification.

#### Regional parameterization

All values from direct human measurements except where noted.

| Parameter | Region 1: Cecum/Ascending | Region 2: Transverse | Region 3: Desc/Sigmoid/Rectum | Source |
|---|---|---|---|---|
| Content state | Liquid/slurry | Semisolid | Formed | Anatomy |
| Water content | ~90–95% | ~85–92% | ~70–80% | Constrained by stool 74.6% |
| pH | 6.0–6.5 | 6.4–6.8 | 6.6–7.0 | Diet-dependent; ±0.3–0.5 |
| Mucosal bacteria (log₁₀ 16S/mg) | 7.4 (6.8–7.9) | 7.3 (6.9–8.0) | 7.5 (6.9–8.1) | Ahmed et al. 2007 |
| Enterobacteria fraction | ~1.3% | ~2.5% | ~0.8% | Ahmed et al. 2007 (P=0.09, NS) |
| Crypt occupancy by E. coli | ~10–16% | ~5% | <1% | Swidsinski 2005 (mouse) |
| Epithelial O₂ (tissue side) | 30 mmHg serosal | Not measured | 39 mmHg (sigmoid) | Electrode |
| Epithelial O₂ (apical surface) | Broad prior: 2–10 mmHg | 2–10 mmHg | 2–10 mmHg | No segment data |
| Luminal O₂ | 1–10 mmHg | ~1 mmHg | ~1 mmHg | EPR oximetry |
| Inner mucus (bacteria-free) | 40–100 µm | 40–100 µm | 40–100 µm | Swidsinski 2007 |
| Outer mucus (colonized) | 0–200 µm, patchy | 0–200 µm | 0–200 µm | Heterogeneous |
| Dominant motility | Mixing + retrograde | Bidirectional | Rhythmic 2–6/min + storage | Manometry |
| Segment transit time | ~8–15 h | ~5–10 h | ~5–15 h | Total 25–40 h |

HAPCs: ~6/day total, 95% antegrade, 80% daytime, 120 mmHg, propagate >30 cm.
Cyclic motor patterns: 2–6/min, local mixing over 3–5 cm haustral segments.
Retrograde patterns: ~6→35 per 2h postprandially.

Note: Enterobacteria fraction does NOT show a significant proximal-distal
gradient (P=0.09). Initialize E. coli at ~1–2% of mucosal community
uniformly, with wide individual-level uncertainty (0.1–5%).

#### Luminal growth compartment

Each region contains a luminal compartment representing bacteria detached
from mucus and growing in transit:

- **Growth rate**: µ_anaerobic × Monod(carbon_lumen) × (1 - acid_inhibition).
  Luminal E. coli grows anaerobically (no O₂ in lumen) at the lower yield.
  Growth in liquid proximal content is faster than in formed distal content.
- **Inflow**: mucosal detachment + upstream region outflow.
- **Outflow**: downstream region inflow (or stool for Region 3).
- **Mixing**: well-mixed within each region (haustral churning at 2–6/min).

#### Mucosal-to-luminal transfer

Three mechanisms deliver mucosal bacteria to the lumen:

1. **Mucus turnover**: Outer mucus residence ~1–3 h. First-order:
   k_mucus = 0.3–1.0 /h.

2. **Contraction-driven detachment**: Superimposed on mucus turnover.
   k_contract = k₀ + α_HAPC × I_HAPC(t), where I_HAPC is stochastic ~6/day.

3. **Growth-edge shedding**: Daughter cells at colony-lumen interface displaced
   into flow. J_growth = f_edge × µ × N_mucosal. f_edge uncertain (0–1).

#### Wall-to-stool equation

Stool CFU/g is an emergent prediction:

```
C_stool = N_lumen_region3(t_defecation) / M_stool

N_lumen,i+1 = (N_lumen,i × survival + k_detach,i+1 × N_mucosal,i+1 × A_i+1)
              × exp(r_luminal,i+1 × T_transit,i+1)
```

Stool output constraints (validation targets, NOT inputs):
- Wet stool: 100–200 g/day (median 128 g/day)
- Water content: ~75%
- Defecation frequency: ~1.2/day
- Total fecal bacteria: ~0.35–3.2×10¹¹ cells/g
- E. coli: ~10⁶–10⁸ CFU/g (healthy human range)

#### Inter-region coupling

- **Continuous drift**: Luminal content moves distally at ~cm/h.
- **HAPCs** (~6/day): Mass antegrade transfer >30 cm per event.
- **Retrograde patterns**: Proximal retention, increase postprandially.
- **Luminal reattachment**: Bacteria can reattach to mucus downstream.
  p_attach = p₀ × (1 - occupancy_fraction).

## Biophysical Parameterization (Inputs)

All inputs from first-principles measurements, NOT from E. coli spatial
microscopy. The microscopy data (Greter, Swidsinski, Li) is reserved for
validation.

### E. coli cell biology (from Spec 12 + Palsson)
- µ_max aerobic: 5e-4 /s (~1.8 /h in model units)
- µ_max anaerobic: 0.55 × µ_max aerobic (Varma & Palsson 1994)
- Monod K_carbon: 5e-3 mol/m³ (existing)
- Monod K_O2: 3e-7 mol/m³ (Spec 12)
- O₂ consumption: Pirt model, q_maintenance = 1e-18 mol/s/cell (existing)
- Yield aerobic/anaerobic: 4.1× ratio (Varma & Palsson 1994)

### Bacteriocin biophysics (existing model)
- Diffusion coefficient in mucus: from QSSA Green's function
- Killing rate: receptor-dependent (existing)
- Metabolic cost of production: existing
- Receptor trade-offs: BtuB, FepA (existing)

### Mucus rheology (Greter et al. 2026)
- Yield stress: 46 ± 9 Pa
- Below yield: solid-like (cells immobile)
- Above yield: shear thinning, flow

### O₂ supply (Spec 12 literature)
- Epithelial flux: 1–5 µM at mucus surface (not tissue pO₂)
- Luminal: ~0 in distal colon
- D_O2 in mucus: 8e-10 m²/s

### Carbon supply (existing model)
- Mucin liberation: 5e-5 mol/m³/s (existing)
- VBF consumption: Monod, vmax 5.5e-5 mol/m³/s (existing)

### Colonic motility (literature)
- Contraction frequency: 0.2–1.0 /min (human colonic manometry)
- Amplitude: sufficient to exceed mucus yield stress locally
- Regional variation: higher frequency proximally

### Mucus turnover (histology)
- Goblet cell secretion rate: measured
- Mucus thickness by region: 50–400 µm (Johansson, Swidsinski)
- Turnover time: 1–6 hours (region-dependent)

---

## Emergent Predictions (Outputs to Compare with Independent Data)

These are quantities the model predicts from the biophysical inputs above,
testable against independently measured data that was NOT used for
parameterization.

### 1. Strain retention time
- **Prediction**: How long a resident E. coli strain persists before
  replacement. Emerges from patch survival probability × reseeding rate ×
  competitive exclusion within patches.
- **Test data**: Longitudinal metagenomics — Zhao et al. (resident strains
  persist months-years, transients wash out in days), Faith et al. (stability
  of human gut strains over years), Garud et al. (strain replacement dynamics).
- **Diagnostic**: track individual lineage survival curves across the colonic
  section.

### 2. Invasion resistance (minimum viable inoculum)
- **Prediction**: How many cells of a new strain must arrive orally to
  establish in the colon. Emerges from: luminal transit survival × probability
  of landing on a receptive patch × growth rate vs local washout rate ×
  bacteriocin survival.
- **Test data**: Streptomycin mouse model — colonization resistance measured
  as minimum gavage dose for establishment (~10⁴–10⁸ CFU depending on
  resident community).
- **Diagnostic**: sweep inoculum size, report establishment frequency.

### 3. Diversity ceiling
- **Prediction**: Maximum number of coexisting E. coli strains. Emerges from:
  number of distinct patch types × competitive exclusion rate within patches ×
  migration/exchange rate between patches.
- **Test data**: Longitudinal metagenomics shows 1–4 dominant E. coli strains
  per host (Schloissnig et al., Yassour et al., Chen-Liaw/Yilmaz et al. 2024).
- **Diagnostic**: initialize with 10 strains, track how many persist at 30, 90,
  365 days.

### 4. Antibiotic perturbation and recovery
- **Prediction**: Population decline during antibiotic exposure, recovery
  kinetics after cessation. Emerges from: antibiotic kill rate (input) ×
  crypt refuge survival × recolonization rate from surviving patches.
- **Test data**: Streptomycin/ciprofloxacin mouse experiments with time-course
  CFU/g stool. Mondragón-Palomino et al. (2022) showed ciprofloxacin reduced
  fecal load >1000× but crypt populations partially survived and recovery was
  altered.
- **Requires**: Antibiotic module (see §Implementation below).

### 5. Bacteriocin producer advantage
- **Prediction**: Competitive advantage of a colicin-producing strain over a
  sensitive strain, measured as relative abundance over time. Emerges from:
  bacteriocin diffusion range relative to patch size × metabolic cost ×
  patch exchange rate.
- **Test data**: Kirkup & Riley (2004) — colicin producers displaced
  sensitive strains in mouse colonization over 2–4 weeks. The magnitude and
  timescale of displacement are independently measured.
- **Diagnostic**: report producer/sensitive ratio over time, compare to
  experimental time series.

### 6. Stool shedding rate
- **Prediction**: CFU/g stool as an emergent output of mucosal population ×
  washout rate × luminal transit. Not an input.
- **Test data**: Healthy human stool E. coli: 10⁶–10⁸ CFU/g.
- **Diagnostic**: count agents exiting the distal segment per unit time,
  normalize by stool volume.

---

## Antibiotic Perturbation Module

Required for prediction #4. Minimal implementation:

```
For each agent in each patch, each timestep:
  if antibiotic_concentration > MIC:
    kill_probability = 1 - exp(-k_kill * (C_abx / MIC - 1) * dt)
    if rng.random() < kill_probability:
      kill agent
```

Antibiotic concentration decays in the lumen after dosing (first-order,
half-life from pharmacokinetics). Mucus-associated concentration is lower
than luminal (penetration factor, ~0.1–0.5 for aminoglycosides/fluoroquinolones
in mucus). Crypt concentration may be lower still.

Parameters:
| Antibiotic | MIC (E. coli) | k_kill | Luminal half-life | Mucus penetration |
|---|---|---|---|---|
| Streptomycin | ~8 µg/mL | 0.5 /h | ~4 h | 0.2 |
| Ciprofloxacin | ~0.015 µg/mL | 2.0 /h | ~6 h | 0.3 |

These are defensible starting values from pharmacokinetic literature. The
model predicts the *spatial pattern* of survival (crypt vs exposed) and the
recovery kinetics — both independently measurable.

---

## Implementation Strategy

### Phase 1: Layer 2 (patch coupling) on top of existing IBM

1. **Define patch types** as config presets (boundary conditions only)
2. **Implement contraction event** as a periodic callback in the simulation
   loop:
   - Count agents per patch
   - For disrupted patches: remove fraction, redistribute to neighbors/lumen
   - For crypt patches: reduced removal probability
3. **Implement luminal pool** as a simple holding compartment with
   first-order decay (washout) and stochastic reattachment
4. **Run validation**: Does the coupled system produce stable populations
   in the 10⁴–10⁶ CFU/mL range without the dysbiosis guard?

### Phase 2: Layer 3 (colonic chain)

1. **Chain M segments** with unidirectional luminal flow
2. **Set regional gradients** from literature (O₂, mucus thickness, pH)
3. **Run predictions 1–3**: strain retention, invasion resistance, diversity
4. **Compare to longitudinal metagenomics data**

### Phase 3: Antibiotic module + prediction #4

1. **Implement antibiotic concentration decay in lumen/mucus/crypts**
2. **Run perturbation experiments** matching streptomycin mouse protocol
3. **Compare recovery kinetics to published time courses**

### Phase 4: Bacteriocin competition campaigns (prediction #5)

1. **Initialize with producer + sensitive + resistant strains**
2. **Run for 30–90 simulated days**
3. **Compare displacement dynamics to Kirkup & Riley (2004)**

---

## Dysbiosis Guard (Revised)

Based on Spec 12 dysbiosis review (task:4e6a870f-561f-4065-a1a4-bb0bd6d781f5):

**Replace the single guard at 1e8 cells/mL with a multi-level system:**

| Level | Threshold | Action |
|---|---|---|
| Healthy operating range | 10⁴–10⁵ CFU/mL per patch | Normal |
| Alert | 10⁵–10⁶ CFU/mL in any patch | Log warning |
| Dysbiosis guard | 10⁶ CFU/mL sustained in >50% of patches | Stop or activate host feedback |
| Model-invalid | ≥10⁷ CFU/mL anywhere, or E. coli >5% of total bacteria | Hard stop |

Plus:
- **Bloom criterion**: ≥10-fold increase over initial baseline, sustained >1 hour
- **Spatial guard**: persistent penetration into inner mucus or epithelial contact

---

## What This Spec Does NOT Model

Explicit non-goals (complexity we're deliberately excluding):

1. **Host immune response** — no neutrophils, no cytokines, no epithelial
   signaling. The dysbiosis guard is where the model stops instead.
2. **Pathotype selection** — no AIEC adherence advantages, no virulence
   factor expression. All strains are commensal E. coli K-12 variants.
3. **VBF dynamics** — the anaerobic majority remains a static metabolic field.
   No Bacteroides competition, no cross-feeding dynamics.
4. **Full fluid dynamics** — contractions are stochastic events, not CFD.
5. **Diet variation** — carbon supply is constant.
6. **Phage** — not modeled.

Each is a potential future extension, but the model's scientific value comes
from explaining observed population dynamics (strain retention, diversity
ceiling, bacteriocin advantage) from biophysical mechanisms alone, without
invoking these processes.

---

## References

Greter et al. 2026. Emergent spatial structure in the gut microbiota is driven
by bacterial growth and gut contractions. PLOS Biology.
doi:10.1371/journal.pbio.3003772

Li et al. 2015. The outer mucus layer hosts a distinct intestinal microbial
niche. Nature Communications. doi:10.1038/ncomms9292

Mondragón-Palomino et al. 2022. Three-dimensional imaging for the
quantification of spatial patterns in microbiota of the intestinal mucosa.
PNAS. doi:10.1073/pnas.2118483119

Swidsinski et al. 2005. Spatial organization of bacterial flora in normal and
inflamed intestine. J Clin Microbiol. doi:10.1128/jcm.43.7.3380-3389.2005

Elliott et al. 2013. Quantification and characterization of mucosa-associated
and intracellular E. coli in IBD. Inflamm Bowel Dis.
doi:10.1097/mib.0b013e3182a38a92

Varma & Palsson 1994. Stoichiometric flux balance models quantitatively predict
growth and metabolic by-product secretion in E. coli W3110. Appl Environ
Microbiol. doi:10.1128/aem.60.10.3724-3731.1994

Kirkup & Riley 2004. Antibiotic-mediated antagonism leads to a bacterial game
of rock-paper-scissors in vivo. Nature. doi:10.1038/nature02429

Bakkeren, Foster et al. 2025. Strain displacement in microbiomes via ecological
competition. Nature Microbiology. doi:10.1038/s41564-025-02162-w

Baertschi, Jordi, Yilmaz et al. 2026. Strain-level ecological filtering governs
microbial colonization of the human gut. Cell Reports.
doi:10.1016/j.celrep.2026.115629

Ahmed et al. 2007. Mucosa-associated bacterial diversity in relation to human
terminal ileum and colonic biopsy samples. Appl Environ Microbiol.
doi:10.1128/aem.01143-07

Ahmed et al. 2007. Mucosa-associated bacterial diversity in relation to human
terminal ileum and colonic biopsy samples. Appl Environ Microbiol.
doi:10.1128/aem.01143-07

Greter et al. 2026. Emergent spatial structure in the gut microbiota is driven
by bacterial growth and gut contractions. PLOS Biology.
doi:10.1371/journal.pbio.3003772
