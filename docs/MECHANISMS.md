# Biological Mechanism Reference

Detailed descriptions of each Fix module, their biological basis, and implementation.

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
- **BtuB loss** (expr < 0.5): Activates MetE pathway for B12-independent methionine synthesis. Cost = `metE_penalty` (default 5%) + ethanolamine utilization loss `eut_penalty` (default 3%).
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
| Lethal Halo | < 6.0 | R = 1.5 | Repelled by anionic mucin, spreads widely |
| Neutral | 6.0–8.5 | R = 5.0 | Intermediate diffusion |

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

**Immunity:** Agents carrying the cognate immunity protein (from their BI cluster) get 1000x protection (`immunity_factor = 0.001`).

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

**Effect:** Recipient gains a copy of one randomly-selected BI cluster from the donor. This enables the spatial spreading of bacteriocin-immunity cassettes through the population.

---

## 5. fix_mutation — Stochastic Genome Evolution

**Five mutation classes:**

### a) BI locus duplication (rate: 10^-5 per division)
Tandem duplication of an existing BI cluster. Expands the agent's toxin range. Capped at `max_bi_loci` (default 8).

### b) BI locus recombination (rate: 5 × 10^-6 per division)
Swaps immunity proteins between two BI clusters. Creates novel toxin/immunity combinations that may escape existing immunity.

### c) Receptor downregulation (rate: 10^-7 per division)
Reduces expression of a randomly-selected receptor by `receptor_reduction` (default 0.1). When expression drops below 0.2, agent transitions to RESISTANT phenotype. Provides toxin resistance at the cost of nutrient uptake.

### d) Super-killer emergence (rate: 10^-8 per division)
Generates a novel toxin variant from an existing BI cluster. The new toxin has:
- Unique toxin_id
- Gaussian-perturbed pI (std = 0.5)
- 30% chance of switching target receptor
This models the evolution of "super-killers" that bypass cognate immunity.

### e) Compensatory mutation (rate: 10^-6 per division)
Chromosomal mutations that ameliorate plasmid metabolic cost. Reduces per-locus cost by `compensatory_reduction` (default 0.005), capped so that at most 75% of the original cost can be eliminated. This models the well-documented phenomenon where bacteria rapidly evolve to reduce the fitness cost of newly acquired plasmids (VADI Section 79).

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

---

## Viscoelastic Background Field (VBF)

The 99% obligate anaerobic microbiota is modeled as a continuum rather than discrete agents:
- **Physical drag**: Stokes-like force opposing agent velocity
- **Nutrient sink**: Background consumption at volumetric rate
- **Mucin liberation**: Monosaccharide release from mucin glycoproteins (carbon source)
- **Carrying capacity**: Local density limit for the simulation domain
