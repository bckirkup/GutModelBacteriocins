# Spec 12 — Density-Dependent Growth Limitation (as-received, annotated)

This document is the amended Spec 12 as supplied by the project lead, reproduced
below without edits so that the literature-derived reasoning stays auditable.
Everything above the horizontal rule is repository annotation: what is already
implemented, what differs deliberately, which of the spec's numbers and config
keys do not survive checking against the repository's own units and parser, and
which mechanisms the delivery/ROS campaign (#319–#339) added or bypassed after
the spec was written.

## Implementation status

| Spec change | Status in this repository |
|---|---|
| **Change 1** — agent-density-coupled VBF carbon sink | Implemented as `vbf.agent_carbon_coupling` (default `0.0`, so existing configs are unchanged). Counts owned agents only; ghost agents are excluded per the MPI ghost double-count rule in `AGENTS.md`. |
| **Change 2** — O₂-dependent metabolic mode | Implemented before this document, in amended form: `oxygen.metabolic_switch_enabled`, `oxygen.mu_crit`, `oxygen.aerobic_mu_factor` / `anaerobic_mu_factor`, `oxygen.tau_metabolic_switch`, `oxygen.ferm_acid_yield`, and per-agent `realized_fermentation_fraction` (persisted and inherited at division). |
| **Change 3** — fermentation-product growth inhibition | Implemented before this document, in amended form: `metabolism.acid_inhibition_enabled`, `acid_inhibition_max`, `acid_inhibition_Ki`, acting on *undissociated* acetate via Henderson–Hasselbalch at the existing environment pH. |
| Nutrient blocking fraction diagnostic | Implemented as a derived per-interval quantity over the existing agent and VBF carbon accounting. |

## Per-mechanism reconciliation (spec-says / code-does / gap)

The spec text below was written before PRs #319–#339. This table is the
authoritative statement of what the shipped code actually evaluates. Config keys
are quoted exactly as the parser accepts them; unknown keys are warned on
`stderr` and ignored, or abort the run under `GUTIBM_STRICT_CONFIG=1`.

### Change 1 — density-coupled VBF carbon sink

| Item | Spec says | Code does | Gap |
|---|---|---|---|
| Sink law | `vmax_local = vmax_base + coupling · n_agents / V_cell`, implicit positivity-preserving solve unchanged | Identical, in `vbf.cpp :: apply_carbon_sink` (host) and `chemistry_kernel.cu` (device), from the same shared closed-form helper | None |
| Occupancy histogram | "counted via histogram over `grid_cell[]`", wired in `simulation.cpp` | Built in `chemistry_pipeline.cpp :: run_chemistry_pipeline` (host) or `count_agents_per_cell_kernel` (device); dead and ghost agents excluded | Location differs from the spec's file list; behaviour matches. `simulation.cpp` is not involved |
| MPI correctness | Not addressed | Was rank-local in `replicated` chemical decomposition — every rank applies VBF over the whole global grid, so with `nprocs > 1` and nonzero coupling the ranks' carbon fields diverged from each other and from serial. Now globally summed before the sink is applied; the device path defers to the reduced host path in that configuration | **Real defect, found during this reconciliation and fixed in code.** Inert at the default `coupling = 0.0`, but it would have silently corrupted any multi-rank carbon-competition sweep |
| Starting value / sweep | `1e-16 mol/s/agent`; sweep `{0, 1e-17, 1e-16, 1e-15}` | No default other than `0.0` | Spec magnitudes are ~2000× per-agent demand at this discretization; see "Sizing `agent_carbon_coupling`" below. Use `{0, 1e-21, 1e-20, 1e-19}` |

### Change 2 — O₂-dependent metabolic mode

| Item | Spec says | Code does | Gap |
|---|---|---|---|
| Enable flag | `oxygen.metabolic_switch_enabled`, legacy boost preserved when false | Same; the legacy `1 + boost_max · f_O2` path is the `else` branch and is still the default | None |
| Respiratory capacity | `resp_capacity = f_O2 · mu_crit / max(µ, mu_crit)` | Same, in `metabolic_mode::respiratory_capacity` (shared host/device header) | None |
| Growth multiplier | Interpolate `aerobic_mu_factor` → `anaerobic_mu_factor` on the fermentation fraction | Same, via `metabolic_mode::interpolate` | None |
| Yield | `Y_eff = yield_carbon · ((1−ff)·aerobic_yield_factor + ff·anaerobic_yield_factor)`, `anaerobic_yield_factor = 0.25` | Substrate-cost convention: `oxygen.aerobic_carbon_cost_factor = 1.0`, `oxygen.anaerobic_carbon_cost_factor = 4.1` | Reciprocal formulation of the same Varma–Palsson ratio (`1/0.244 ≈ 4.1`). The spec's key names `oxygen.aerobic_yield_factor` / `oxygen.anaerobic_yield_factor` **do not exist** |
| Fermentative acetate | `acetate_prod = ff · ferm_acetate_yield · carbon_consumed`, `ferm_acetate_yield = 0.67` mol acetate/mol glucose | `oxygen.ferm_acid_yield = 1.0` acetate per carbon-equivalent consumed, multiplied by the realized carbon cost | Different quantity, not a different value; the spec key `oxygen.ferm_acetate_yield` does not exist |
| Inertia | Exponential relaxation over `tau_metabolic_switch = 3600 s` | Same, `metabolic_mode::relax` | None |
| Per-agent state | `f_ferm_realized`, initialized `0.0` | `Agent::realized_fermentation_fraction`, persisted in checkpoints, inherited at division, mirrored to the GPU agent pool, and emitted per agent and as `mean_realized_fermentation_fraction` | Name only |
| Driver of the fermentation fraction | Ambient concentration only | `oxygen.respiration_driver` selects `ambient` (spec behaviour, default) or `funded` (#339), where the fraction comes from `funded_growth_O₂ / demanded_growth_O₂` with a one-step lag | **Post-spec mechanism.** In `funded` mode `oxygen.Km` and `oxygen.mu_crit` do not enter the fermentation fraction at all, so the spec's two-axis phase plane — including the aerobic-overflow axis the `mu_crit` correction was made for — is inactive. `funded` requires `metabolism.uptake_limit=delivery` and is CPU-only |
| Maintenance | Not addressed by Change 2 | `oxygen.anaerobic_maintenance_factor = 15.0` interpolates the carbon maintenance rate on the same fermentation fraction | **Code term with no spec counterpart and no citation in this document.** It is a large multiplier and should be sourced or swept before it is relied on |
| Oxygen boundary for Spec 12 runs | `oxygen.epithelial_conc: 5.0e-6` annotated "mucus-surface O₂ (~4 mmHg)" | Parsed as `mol/m³` | **Unit error in the spec.** `5.0e-6 mol/m³` is ≈`0.004 mmHg`; ≈4 mmHg is ≈`5e-3 mol/m³`. The same 1000× error is in the spec's "55 µM (~42 mmHg)" line. See `docs/PARAMETERS.md` §Oxygen, which carries the corrected conversions and the sourced tissue-side value |

### Change 3 — acid-product growth inhibition

| Item | Spec says | Code does | Gap |
|---|---|---|---|
| Inhibition law | `mu *= 1 − acid_inhibition_max · [Ac]/(Ki + [Ac])` on **total** acetate, before maintenance subtraction | Same functional form and same position in the growth expression, but on the **undissociated** fraction via Henderson–Hasselbalch | Deliberate; see "Acid inhibition acts on undissociated acetate" below |
| Half-inhibition constant | `Ki_acetate = 20 mol/m³` total | `acid_inhibition_Ki = 50 mol/m³` undissociated | Deliberate; the spec value contradicts the spec's own citation |
| Config location | `OxygenConfig`, keys `oxygen.acid_inhibition_enabled` / `oxygen.Ki_acetate` / `oxygen.acid_inhibition_max` | `MetabolismConfig`, keys `metabolism.acid_inhibition_enabled`, `metabolism.acid_inhibition_Ki`, `metabolism.acid_inhibition_max` (bare names also accepted) | The spec's `oxygen.*` acid keys **do not exist** |
| pH and pKa | Implicit "gut pH 6.0" | pH is read from the existing environment value `bacteriocin.mucin_charge.ph`; pKa from `metabolism.acetate_pKa` | Acid inhibition is coupled to the bacteriocin module's pH, not to a Spec-12-local pH |

### Diagnostics (spec §Diagnostics, items 1–8)

| # | Spec diagnostic | Status |
|---|---|---|
| 1 | Metabolic state histogram (respiratory / mixed / fermentative) | Partial — per-agent `realized_fermentation_fraction` and `mean_realized_fermentation_fraction` are written; no binned histogram |
| 2 | Overflow acetate split aerobic vs anaerobic | Not implemented |
| 3 | Mean effective yield | Not implemented |
| 4 | Mean acid inhibition | Not implemented |
| 5 | Carbon in top- vs bottom-quartile-density voxels | Not implemented (derivable offline from grid output) |
| 6 | Thiele modulus per z-layer | Not implemented |
| 7 | Metabolic inertia lag | Not implemented |
| 8 | Nutrient blocking fraction | Implemented, but **per species and per output interval, domain-wide** (`agent uptake + maintenance` over `that + VBF sink`), not per voxel as the spec asks |

### Validation config (spec §"Validation config")

Five of its fourteen keys are not parser keys and are silently ignored (or abort
under `GUTIBM_STRICT_CONFIG=1`): `oxygen.aerobic_yield_factor`,
`oxygen.anaerobic_yield_factor`, `oxygen.ferm_acetate_yield`,
`oxygen.acid_inhibition_enabled`, `oxygen.Ki_acetate`, plus
`oxygen.acid_inhibition_max` (the bare `acid_inhibition_max` is accepted). Do
not copy that block. The equivalent working configuration is:

```yaml
oxygen.enabled: true
oxygen.epithelial_conc: 5.0e-6        # see the unit caveat above
oxygen.Km: 3.0e-7
oxygen.metabolic_switch_enabled: true
oxygen.aerobic_mu_factor: 1.0
oxygen.anaerobic_mu_factor: 0.55
oxygen.aerobic_carbon_cost_factor: 1.0
oxygen.anaerobic_carbon_cost_factor: 4.1
oxygen.mu_crit: 9.7e-5
oxygen.ferm_acid_yield: 1.0
oxygen.tau_metabolic_switch: 3600.0
metabolism.acid_inhibition_enabled: true
metabolism.acid_inhibition_Ki: 50.0
metabolism.acid_inhibition_max: 0.8
vbf.agent_carbon_coupling: 1.0e-20     # not 1e-16; see sizing section below
```

Backward compatibility is as the spec claims: `agent_carbon_coupling = 0.0`,
`metabolic_switch_enabled = false`, `acid_inhibition_enabled = false`, and
`epithelial_conc = 55e-6` are the shipped defaults. The separate open question of
whether the *oxygen mortality and uptake* defaults (`oxygen.ros_driver`,
`oxygen.k_ROS`, `metabolism.uptake_limit`) should still be the pre-campaign
values is tracked in `AGENTS.md` and `docs/DELIVERY_ROS_CAMPAIGN.md`, not here.

## Deliberate deviations from the text below

**Yield is expressed as a substrate cost, not a biomass yield.** The spec's
`anaerobic_yield_factor = 0.25` multiplies biomass-per-substrate. The repository
carries the reciprocal, `oxygen.anaerobic_carbon_cost_factor = 4.1`
(substrate-per-biomass), because the uptake path is written in substrate terms.
`1 / 0.244 ≈ 4.1`, so the two express the same Varma–Palsson ratio.

**Acid inhibition acts on undissociated acetate.** The spec inhibits on total
acetate with `Ki_acetate = 20 mol/m³`. The repository instead applies
Henderson–Hasselbalch and inhibits on the undissociated fraction with
`acid_inhibition_Ki = 50 mol/m³ undissociated`, which is precisely the threshold
the spec itself cites (Russell & Diez-Gonzalez 1998: >50% inhibition at ~50 mM
undissociated acetic acid; 50 mM = 50 mol/m³). The spec's `20 mol/m³` total is
≈`1.1 mol/m³` undissociated at pH 6, roughly 45× more sensitive than its own
citation, and its accompanying annotation ("≈66 mM undissociated") does not
follow from `20 mol/m³ × 5.4%`. The repository default is unchanged; the spec
value is not adopted.

**`ferm_acid_yield` is per carbon-equivalent, not per mole of glucose.** The
spec's `ferm_acetate_yield = 0.67 mol acetate / mol glucose` is a different
quantity from the repository's `oxygen.ferm_acid_yield = 1.0`, which is acetate
per carbon-equivalent consumed. This is a units difference rather than a
disagreement.

**`oxygen.mu_crit` is `9.7e-5 /s`, matching the spec.** The former repository
default of `3.0e-4 /s` (1.08 h⁻¹) sat above achievable µ, and since
`resp_capacity = f_O2 · mu_crit / max(µ, mu_crit)`, any such value collapses the
phase plane onto the O₂ axis and the aerobic-overflow axis can never fire — the
same structurally-inert-cap failure the `sherwood` uptake mode had. `9.7e-5 /s`
(0.35 h⁻¹) is the measured overflow onset and is below the example strains'
`mu_max = 5.5e-4 /s`, so the growth-rate axis is live. A sensitivity test pins
this: with O₂ held saturating, agents above `mu_crit` must ferment and secrete
acetate while agents below it must not.

## Sizing `agent_carbon_coupling` (correction to §"Parameter guidance")

The spec's starting value and sweep are several orders of magnitude too large for
this discretization, and every arm of the proposed sweep would zero the carbon in
any occupied voxel:

- At `grid_dx = 2 µm` a voxel is `8e-18 m³` (8 fL).
- The entire background VBF carbon sink in one voxel is
  `5.5e-5 mol/m³/s × 8e-18 m³ = 4.4e-22 mol/s`.
- Measured per-agent carbon demand (delivery-mode validation run, 1800 s /
  30 steps) is `3.6e-14 mol / (≈400 agents × 30 steps) ≈ 5e-20 mol/s`.

So the spec's `1e-16 mol/s/agent` is ≈2000× *full* per-agent demand and ≈2.3e5×
the background sink of the voxel it occupies — not the intended ~5% of agent
demand. Demand-anchored, 5% is ≈`2.5e-21` and 100% is ≈`5e-20`, giving a
defensible sweep of `{0, 1e-21, 1e-20, 1e-19}` (≈2%–200% of per-agent demand).
Values large enough to starve a voxel are expected to halt a `delivery`-mode run
through the closure gate (`demand > 0` with realized removal `== 0`); that is the
gate doing its job, not a defect to be worked around.

---

# GutIBM Spec 12: Density-Dependent Growth Limitation

## Motivation

The model has no density-dependent brake within reach of the dysbiosis guard
(1e8 cells/mL). The Sherwood uptake cap is structurally slack (Monod starvation
arrives at 16× higher concentration than diffusive limitation), and the VBF
carbon sink is spatially uniform and agent-density-independent. Net growth
crosses zero exactly once; there is no basin, so no parameter choice holds a
population for 168 h.

This spec adds three coupled mechanisms that the literature identifies as
physically real at 1e8 cells/mL and absent or inert in the current model:

1. **Agent-density-coupled VBF carbon competition** — VBF consumption increases
   where agents are present, so local carbon drops faster in dense regions.
2. **O₂-dependent metabolic mode with yield reduction and overflow acetate** —
   O₂ depletion forces a shift from respiration to fermentation (lower growth
   rate, lower carbon yield, fermentation product secretion), and high aerobic
   growth rates trigger overflow acetate even with O₂ present.
3. **Fermentation product inhibition** — accumulated acetate inhibits growth,
   closing the feedback loop.

Together these create a consumption–diffusion–product feedback that can produce
a stable carrying capacity without introducing phenomenological density
dependence.

## Literature basis

### Spatial nutrient competition (Change 1)

iDynoMiCS, NUFEB, MetaBiome, and MICRODIMS all couple background consumption
to the same reaction-diffusion grid as agents, with local Monod uptake
proportional to local biomass. Nutrient depletion alone creates starved
interiors but does not produce a finite whole-domain carrying capacity without
an additional loss term. Every published model that achieves steady state adds
bounded space, detachment, death, or washout
(task:19dd2a60-2f30-4b43-9628-7ccd65e5e2e4).

### O₂ depletion and metabolic switching (Change 2)

At microaerobic boundary conditions (1–5 µM at the mucus surface), critical
density for complete O₂ depletion across 100 µm mucus is ~1–5 × 10⁸ cells/mL.
E. coli O₂ consumption: q_O2,max ≈ 15–20 mmol/gDW/h → ~1.5e-18 mol/cell/s;
K_O ≈ 0.1–3 µM (task:91d151e3-a5c3-4a6a-825e-8e8da5d90bba).

No published model achieves stable E. coli carrying capacity from O₂ alone —
E. coli is facultatively anaerobic. The MICRODIMS E. coli IBM generates
no-growth zones from local acidification via the fermentation product cascade
(Tack et al. 2017, doi:10.3389/fmicb.2017.02509). Henson & Phalak (2017) gut
biofilm model required biomass erosion for steady states
(doi:10.1186/s12918-017-0522-1).

### Palsson metabolic framework (Changes 2–3)

The quantitative parameterization of the aerobic–fermentative transition derives
from the Varma–Palsson stoichiometric flux balance studies of E. coli W3110.
Their strain-calibrated model established the oxygen–carbon phase plane that
governs metabolic mode selection (Varma & Palsson 1994,
doi:10.1128/aem.60.10.3724-3731.1994; Varma, Boesch & Palsson 1993,
doi:10.1128/aem.59.8.2465-2473.1993):

**Measured capacities (W3110):**

| Parameter | Value | Source |
|---|---:|---|
| Maximum O₂ uptake | 15 mmol gDW⁻¹ h⁻¹ | Varma & Palsson 1994 |
| Max aerobic glucose uptake | 10.5 mmol gDW⁻¹ h⁻¹ | Varma & Palsson 1994 |
| Max anaerobic glucose uptake | 18.5 mmol gDW⁻¹ h⁻¹ | Varma & Palsson 1994 |
| Aerobic biomass yield | 0.029 gDW/mmol glucose | Varma, Boesch & Palsson 1993 |
| Anaerobic biomass yield | 0.023 gDW/mmol glucose | Varma & Palsson 1994 |
| Aerobic yield (mass basis) | 0.524 gDW/g glucose | Varma, Boesch & Palsson 1993 |
| Anaerobic yield (mass basis) | 0.128 gDW/g glucose | Varma & Palsson 1994 |
| Anaerobic max growth | 0.43 h⁻¹ | Varma & Palsson 1994 |
| NGAM | 7.6 mmol ATP gDW⁻¹ h⁻¹ | Varma & Palsson 1994 (W3110) |
| NGAM | 8.39 mmol ATP gDW⁻¹ h⁻¹ | iAF1260 (Feist et al. 2007) |

The yield ratio is 0.524/0.128 ≈ 4.1×, meaning anaerobic growth consumes ~4×
more carbon per unit biomass than aerobic growth. This is the primary mechanism
by which O₂ depletion accelerates carbon exhaustion.

**Acetate overflow** does not require anaerobiosis. Under aerobic carbon-limited
growth, acetate appears when respiratory O₂ demand hits the 15 mmol/gDW/h
ceiling, at a critical growth rate µ_crit ≈ 0.3–0.4 h⁻¹ (Varma & Palsson
1994). Sensitivity analysis showed acetate predictions were particularly
sensitive to maximum O₂ and glucose uptake capacities (Varma & Palsson 1995,
doi:10.1002/bit.260450110). The Palsson group's later ME-model work (Chen et
al. 2021, doi:10.1371/journal.pcbi.1008596) showed this overflow is a proteome
allocation phenomenon: respiratory machinery has higher per-ATP proteomic cost,
so at high growth rates the cell preferentially routes carbon through
fermentation even when O₂ is available. This is not a pathology — it is the
optimal metabolic strategy under translational resource limitation.

**The oxygen phase plane** predicts a sequence with declining O₂ availability:
(1) fully respiratory, (2) acetate secretion, (3) acetate + formate,
(4) acetate + formate + ethanol near anaerobiosis. Below ~0.88 mmol O₂/gDW/h,
all three major anaerobic products coexist (Varma, Boesch & Palsson 1993).
The canonical mixed-acid stoichiometry is:
1 glucose → 1 acetate + 1 ethanol + 2 formate
(before biomass allocation; PFL-dominated, measured for W3110).

**Dynamic behavior**: Mahadevan, Edwards & Doyle (2002,
doi:10.1016/s0006-3495(02)73903-9) formulated dFBA for E. coli diauxic growth
using the Palsson metabolic network, with kLa = 7.5 h⁻¹ and O₂ reference =
0.21 mM. DynamicME (Yang et al. 2019, doi:10.1186/s12918-018-0675-6) showed
that proteome reallocation between respiratory and fermentative machinery
operates on a timescale of ~1–2 h. Instantaneous metabolic switching is
therefore unrealistic; the model should include inertia.

(task:d5b38f7a-36ab-41e1-b47a-8e0a2a91a9e1)

---

## Change 1: Agent-density-coupled VBF carbon sink

### Current behavior

The VBF carbon sink applies `implicit_carbon_sink(C, vmax, Km, dt)` uniformly
in every voxel with `vmax = 5.5e-5 mol/m³/s`. Agent presence does not modulate
VBF consumption.

### New behavior

Replace `vmax` with a density-modulated value:

```
vmax_local = vmax_base + agent_carbon_coupling * n_agents_in_voxel / V_cell
```

- `vmax_base` = unchanged background (5.5e-5 mol/m³/s)
- `agent_carbon_coupling` = per-cell additional competition rate (mol/s/cell)
- `n_agents_in_voxel` = counted via histogram over `grid_cell[]` assignments

The implicit carbon sink solver (positivity-preserving quadratic) is unchanged;
only `vmax` varies spatially.

### Implementation

#### Config (`VBFConfig`)

```cpp
Real agent_carbon_coupling = 0.0;  // mol/s per agent; 0 = off (backward compat)
```

#### GPU

Add a histogram kernel before VBF coupling:

```cpp
__global__ void count_agents_per_cell_kernel(
    const int* grid_cell, const int* state,
    int num_agents, int* counts, int ncells) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= num_agents) return;
  if (state[i] == 3) return;
  const int cell = grid_cell[i];
  if (cell >= 0 && cell < ncells) atomicAdd(&counts[cell], 1);
}
```

Modify `apply_vbf_at_cell` to read `agent_counts[cell]` and compute
`vmax_eff`. Pass `int* agent_counts` through `VbfLaunchParams`.

#### CPU

Same logic in `vbf.cpp :: apply_carbon_sink`. Compute agent-per-voxel histogram
before calling `apply_nutrient_coupling`.

#### Parameter guidance

Starting value: `agent_carbon_coupling = 1e-16 mol/s/agent` (~5% of direct
agent carbon demand at full growth). Sensitivity range: 1e-17 to 1e-15.

---

## Change 2: O₂-dependent metabolic mode

### Current behavior

O₂ gives a multiplicative bonus:
```cpp
mu *= (1.0 + boost_max * O2 / (Km + O2))    // boost_max = 2.0
```
Agents grow at full baseline rate without O₂; O₂ is purely beneficial.
Epithelial boundary is 55 µM (~42 mmHg tissue pO₂).

### New behavior

When `metabolic_switch_enabled = true`, replace the boost with a two-axis
metabolic mode selector based on the Varma–Palsson phase plane:

**Axis 1: O₂-dependent respiratory fraction**

```cpp
f_O2 = O2 / (Km_O2 + O2);   // 0 = fully anaerobic, 1 = fully aerobic
```

This determines the maximum fraction of carbon flux that can be respired.

**Axis 2: Growth-rate-dependent overflow**

Even with O₂ available, E. coli produces acetate above µ_crit because
respiratory O₂ demand exceeds the uptake ceiling (Varma & Palsson 1994) and
because respiratory proteome becomes too expensive at high translation rates
(Chen et al. 2021). The overflow fraction is:

```cpp
// Respiratory capacity fraction: how much of current growth the respiratory
// machinery can support, given available O₂
resp_capacity = f_O2 * mu_crit / max(mu, mu_crit);
// resp_capacity = 1.0 when mu <= mu_crit and O₂ is saturating
// resp_capacity < 1.0 when mu > mu_crit OR O₂ is limiting

// Fermentation fraction: what can't be respired must be fermented
ferm_fraction = 1.0 - resp_capacity;
```

**Growth rate and yield**

```cpp
// Effective growth factor: weighted average of aerobic and anaerobic rates
mu_factor = (1.0 - ferm_fraction) * aerobic_mu_factor
          + ferm_fraction * anaerobic_mu_factor;
mu *= mu_factor;

// Effective carbon yield: weighted by metabolic mode
Y_eff = yield_carbon * ((1.0 - ferm_fraction) * aerobic_yield_factor
                       + ferm_fraction * anaerobic_yield_factor);
```

**Fermentation product secretion**

```cpp
// Acetate production proportional to fermentative carbon flux
if (d_biomass > 0.0 && ferm_fraction > 0.0) {
  carbon_consumed = d_biomass * Y_eff;
  acetate_prod = ferm_fraction * ferm_acetate_yield * carbon_consumed;
  // → atomicAdd to reac_acetate[cell]
}
```

**Metabolic inertia**

DynamicME shows ~1–2 h proteome reallocation timescale (Yang et al. 2019).
The realized fermentation fraction should not track f_O2 instantaneously:

```cpp
// Exponential smoothing with time constant tau_switch
// f_ferm_realized tracks ferm_fraction with inertia
const double alpha = 1.0 - exp(-dt / tau_metabolic_switch);
f_ferm_realized[i] += alpha * (ferm_fraction - f_ferm_realized[i]);
```

Use `f_ferm_realized` instead of `ferm_fraction` in all yield and secretion
calculations. `tau_metabolic_switch` default = 3600 s (1 h), range 1800–7200 s.

### Config additions (`OxygenConfig`)

```cpp
// Spec 12 §2 — metabolic mode (Palsson phase-plane model)
bool metabolic_switch_enabled = false;  // false = legacy boost (backward compat)

// Growth rate factors (relative to agent mu_max)
// Aerobic max growth ~0.7–0.8 h⁻¹; anaerobic max 0.43 h⁻¹ (W3110)
// → anaerobic_mu_factor ≈ 0.43/0.78 ≈ 0.55
Real aerobic_mu_factor    = 1.0;
Real anaerobic_mu_factor  = 0.55;  // Varma & Palsson 1994

// Carbon yield factors
// Aerobic: 0.524 gDW/g glucose; anaerobic: 0.128 gDW/g glucose
// → ratio = 0.128/0.524 ≈ 0.244
Real aerobic_yield_factor   = 1.0;
Real anaerobic_yield_factor = 0.25;  // Varma & Palsson 1994 / Varma et al. 1993

// Overflow critical growth rate (1/s)
// Acetate overflow begins at µ_crit ≈ 0.35 h⁻¹ = 9.7e-5 /s
// even under saturating O₂ (Varma & Palsson 1994)
Real mu_crit = 9.7e-5;  // 1/s; range 7.5e-5 to 1.25e-4 (0.27–0.45 h⁻¹)

// Fermentation product stoichiometry (mol product / mol carbon consumed)
// PFL-dominated mixed-acid: 1 glucose → 1 acetate + 1 ethanol + 2 formate
// Simplified: acetate as dominant growth-inhibiting product
Real ferm_acetate_yield = 0.67;  // mol acetate / mol glucose catabolized

// Metabolic inertia: proteome reallocation timescale
// DynamicME: ~1–2 h (Yang et al. 2019)
Real tau_metabolic_switch = 3600.0;  // s; range 1800–7200
```

### Boundary condition

For Spec 12 runs:
```yaml
oxygen.epithelial_conc: 5.0e-6    # mucus-surface O₂ (~4 mmHg), not tissue pO₂
oxygen.Km: 3.0e-7                  # 0.3 µM
oxygen.metabolic_switch_enabled: true
```

### Implementation

#### Per-agent state

Add to agent struct / GPU buffers:
- `f_ferm_realized` (double): smoothed fermentation fraction, initialized to 0.0

#### Agent update kernel

```cpp
if (o2_enabled && o2_metabolic_switch) {
  const double s_o2 = conc_oxygen ? conc_oxygen[cell] : 0.0;
  const double f_o2 = s_o2 / (o2_Km + s_o2);

  // Instantaneous fermentation fraction from phase-plane logic
  const double resp_cap = f_o2 * fmin(1.0, mu_crit / fmax(mu, 1e-30));
  const double ferm_inst = 1.0 - fmin(resp_cap, 1.0);

  // Apply metabolic inertia
  const double alpha = 1.0 - exp(-dt / tau_metabolic_switch);
  f_ferm_realized[i] += alpha * (ferm_inst - f_ferm_realized[i]);
  const double ff = f_ferm_realized[i];

  // Scale growth rate
  const double mu_factor = (1.0 - ff) * o2_aerobic_mu_factor
                         + ff * o2_anaerobic_mu_factor;
  mu *= mu_factor;
} else if (o2_enabled && conc_oxygen) {
  // Legacy boost
  const double s_o2 = conc_oxygen[cell];
  mu *= (1.0 + o2_boost_max * s_o2 / (o2_Km + s_o2));
}
```

After biomass update:

```cpp
// Yield-adjusted carbon uptake
double eff_yield = yield_carbon;
if (o2_metabolic_switch) {
  const double ff = f_ferm_realized[i];
  eff_yield = yield_carbon * ((1.0 - ff) * o2_aerobic_yield_factor
                             + ff * o2_anaerobic_yield_factor);
}
const double uptake = d_biomass * eff_yield;

// Fermentation product secretion
if (o2_metabolic_switch && d_biomass > 0.0) {
  const double ff = f_ferm_realized[i];
  if (ff > 0.0 && reac_acetate) {
    const double carbon_consumed = d_biomass * eff_yield;
    const double acetate_prod = ff * ferm_acetate_yield
        * carbon_consumed / (cell_volume * dt);
    atomicAdd(&reac_acetate[cell], acetate_prod);
  }
}
```

---

## Change 3: Acid-product growth inhibition

### Current behavior

Acetate inhibits growth only indirectly via MetE pathway penalty when BtuB
expression is low (`expr_btuB < 0.5`). No general acid stress.

### New behavior

Add a direct acid-stress growth penalty based on local acetate:

```cpp
if (acid_inhibition_enabled && conc_acetate) {
  const double acetate = conc_acetate[cell];
  const double inhibition = acetate / (Ki_acetate + acetate);
  mu *= (1.0 - acid_inhibition_max * inhibition);
}

mu -= maintenance_rate;
```

This applies to ALL agents regardless of BtuB status. The existing MetE penalty
remains as an additional, genotype-specific effect.

### Config

```cpp
bool acid_inhibition_enabled = false;
Real Ki_acetate = 20.0;           // mol/m³ (~1.2 g/L total) half-inhibition
Real acid_inhibition_max = 0.8;   // max fractional reduction at saturating acetate
// At typical gut pH 6.0, ~5.5% of total acetate is undissociated (pKa 4.76).
// Ki = 20 mol/m³ total ≈ 1.1 mol/m³ undissociated ≈ 66 mM undissociated.
// Literature: >50% growth inhibition at ~50 mM undissociated acetic acid
// (Russell & Diez-Gonzalez 1998).
// Sensitivity range: Ki = 5–50 mol/m³.
```

---

## Feedback cascades (complete)

### Cascade A: O₂ → yield → carbon (primary density brake)

```
↑ agent density
  → ↑ O₂ consumption (Pirt: q_maintenance per cell)
    → ↓ local O₂ → ↓ f_O2 → ↑ ferm_fraction
      → ↓ mu_factor (slower growth)                        [Change 2]
      → ↓ Y_eff (4× more carbon per biomass)               [Change 2]
        → ↑ carbon consumption per unit growth
          → ↓ local carbon (amplified by Change 1)          [Change 1]
            → ↓ Monod carbon → ↓ mu
```

### Cascade B: overflow → acetate → inhibition (high-growth brake)

```
↑ growth rate (even with O₂)
  → mu > mu_crit → overflow acetate                         [Change 2]
    → ↑ local [acetate]
      → acid inhibition → ↓ mu                              [Change 3]
        → mu drops below mu_crit → overflow stops
          → acetate consumed/diffuses → inhibition relaxes
            → growth resumes → cycle
```

### Cascade C: anaerobic fermentation → acetate → inhibition (density brake)

```
↑ density → ↓ O₂ → ↑ ferm_fraction
  → ↑ anaerobic acetate secretion                           [Change 2]
    → ↑ local [acetate]
      → acid inhibition → ↓ mu                              [Change 3]
```

Cascades A and C operate in parallel. At equilibrium, net growth
(mu × mu_factor × (1 - acid_inhibition) - maintenance) = 0. Three independent
braking forces create a basin.

---

## Diagnostics

Per output interval, report:

1. **Metabolic state histogram**: fraction of agents with ferm_fraction
   < 0.1 (respiratory), 0.1–0.5 (mixed), > 0.5 (fermentative)
2. **Overflow acetate**: total acetate secretion from agents with f_O2 > 0.5
   (aerobic overflow) vs f_O2 < 0.5 (anaerobic fermentation)
3. **Effective yield**: mean Y_eff across agents
4. **Acid inhibition**: mean (1 - acid_inhibition_max × [Ac]/(Ki + [Ac]))
5. **Carbon competition**: mean [carbon] in top-quartile-density voxels vs
   bottom-quartile
6. **Thiele modulus**: Φ² = n·q_cell·L²/(D_O2·(K_O + C_ref)) per z-layer
7. **Metabolic inertia diagnostic**: mean |ferm_inst - f_ferm_realized| (how
   far the smoothed fraction lags the instantaneous value)
8. **Nutrient blocking fraction**: ratio of agent carbon uptake to total carbon
   uptake (agent + VBF) per voxel. This is the resource-competition intensity
   in the Bakkeren/Foster displacement framework — high blocking means VBF
   dominates the carbon budget and bacteriocin effects are suppressed; low
   blocking means agents control their own resource environment and interference
   competition (bacteriocins) determines outcomes

---

## Validation plan

### Stage 1: Unit verification (no campaign)

a. **Change 1 alone**: `agent_carbon_coupling > 0`, other changes off. Verify
   carbon drops in agent-dense voxels. Guard time should increase.

b. **Change 2 alone**: `metabolic_switch_enabled = true`, `epithelial_conc =
   5e-6`. Verify agents near lumen have higher ferm_fraction and produce
   acetate. Verify aerobic overflow at high mu.

c. **Change 3 alone**: Inject high acetate, verify growth reduction.

d. **All three**: Target equilibrium at < 1e8 cells/mL for 168 h.

### Stage 2: Operating-envelope sweep

3D sweep: `epithelial_conc × agent_carbon_coupling × Ki_acetate`
Target: steady-state density in 500–5000 cells/mm².

### Stage 3: RPS campaign with resource-competition gradient

The `agent_carbon_coupling` parameter is not merely a carrying-capacity knob —
it controls the resource-competition intensity that determines whether
bacteriocin-mediated interference competition can operate (Bakkeren/Foster
2025). The RPS campaign must therefore sweep it as a first-class experimental
axis:

- `agent_carbon_coupling ∈ {0, 1e-17, 1e-16, 1e-15}` — spanning from no
  community resource competition to strong nutrient blocking
- At each coupling level, run the full bacteriocin competition matrix
  (producer × sensitive × resistant)
- Report: bacteriocin effectiveness (displacement rate, coexistence time) as a
  function of nutrient blocking fraction
- Hypothesis: RPS coexistence requires an intermediate coupling value where
  resource competition is low enough for bacteriocins to matter but high enough
  to prevent runaway growth. Too low → fastest grower wins regardless. Too
  high → everyone starves and bacteriocins are irrelevant.
- 168 h horizon at the operating points identified in Stage 2.

---

## Backward compatibility

All off by default:
- `agent_carbon_coupling = 0.0`
- `metabolic_switch_enabled = false`
- `acid_inhibition_enabled = false`
- `epithelial_conc = 55e-6`

Existing configs identical.

---

## Files to modify

| File | Change |
|------|--------|
| `src/fields/vbf.h` | Add `agent_carbon_coupling` to `VBFConfig` |
| `src/fields/vbf.cpp` | Compute `vmax_eff` from agent count |
| `src/gpu/chemistry_kernel.cu` | Agent histogram kernel; modify VBF coupling |
| `src/gpu/gpu_kernels.h` | Declare histogram launcher |
| `src/gpu/gpu_types.h` | `agent_carbon_coupling` + `int* agent_counts` in `VbfLaunchParams` |
| `src/gpu/chemistry_pipeline.h` | Call histogram before VBF |
| `src/io/chem_environment_config.h` | Metabolic switch + acid inhibition fields |
| `src/io/input_parser.cpp` | Parse new keys |
| `src/fixes/fix_metabolism.cpp` | Phase-plane metabolic switch; yield adjustment; acetate secretion; acid inhibition |
| `src/gpu/agent_update_kernel.cu` | Same on GPU; `f_ferm_realized` buffer |
| `src/gpu/agent_pool_gpu.h/.cpp` | `d_f_ferm_realized_` device buffer |
| `src/core/agent.h` | Add `f_ferm_realized` to agent struct |
| `src/core/simulation.cpp` | Wire histogram; pass to VBF |
| `tests/` | Unit tests per change + combined |

---

## Appendix A — implementation brief as received

Reproduced verbatim. Where it conflicts with the annotations above (mu_crit,
acid inhibition units, coupling magnitude, and the yield/cost convention), the
annotations are authoritative and the code follows them.

# GutIBM Spec 12 — Devin Implementation Brief

## Summary

Three changes to add density-dependent growth limitation to the GutIBM
individual-based model. All changes are backward-compatible (off by default).
Branch from `main`.

## Spec document

Full spec with literature citations, code snippets, and parameter tables:
`/workspace/gutibm_spec12_density_limitation.md`

The spec was developed from four literature reviews (spatial nutrient
competition in biofilm IBMs, O₂ depletion in gut mucus, contact inhibition,
and Palsson's E. coli metabolic framework) plus the Yilmaz et al. 2026
strain-level ecological filtering paper.

## Changes in priority order

### Change 1: Agent-density-coupled VBF carbon sink

**What**: The VBF carbon consumption rate in each voxel increases proportionally
to the number of agents in that voxel.

**Key equation**:
```
vmax_local = vmax_base + agent_carbon_coupling * n_agents_in_voxel / V_cell
```

**Files to modify**:
- `src/fields/vbf.h` — Add `agent_carbon_coupling` (Real, default 0.0) to `VBFConfig`
- `src/fields/vbf.cpp` — Pass agent count to `apply_carbon_sink`, compute `vmax_eff`
- `src/gpu/chemistry_kernel.cu` — Add `count_agents_per_cell_kernel` histogram kernel; modify `apply_vbf_at_cell` to use `agent_counts[cell]`
- `src/gpu/gpu_kernels.h` — Declare histogram kernel launcher
- `src/gpu/gpu_types.h` — Add `agent_carbon_coupling` and `int* agent_counts` to `VbfLaunchParams`
- `src/gpu/chemistry_pipeline.h` — Call histogram kernel before VBF coupling
- `src/core/simulation.cpp` — Allocate histogram buffer, wire to VBF coupling
- `src/io/input_parser.cpp` — Parse `vbf.agent_carbon_coupling`

**Default value**: `agent_carbon_coupling = 0.0` (no effect, backward compatible)
**Starting value for runs**: `1e-16 mol/s/agent`
**Campaign sweep**: `{0, 1e-17, 1e-16, 1e-15}` — this parameter controls
resource-competition intensity and determines whether bacteriocin-mediated
displacement can operate (Bakkeren/Foster, Nature Microbiology 2025). It is a
first-class experimental axis, not just a calibration knob.

**Test**: With coupling on, carbon concentration should drop in agent-dense
voxels compared to empty voxels. Guard-crossing time should increase.

**New diagnostic**: Report the nutrient blocking fraction per voxel —
ratio of agent carbon uptake to total (agent + VBF) carbon uptake. Add to
the per-interval output alongside existing flux accounting.

### Change 2: O₂-dependent metabolic mode (Palsson phase-plane)

**What**: Replace the additive O₂ boost with a metabolic mode selector that
has two axes (O₂ availability AND growth rate), matching the Varma-Palsson
oxygen-carbon phase plane. Includes yield reduction, fermentation acetate
secretion, and metabolic inertia.

**Key equations**:
```
f_O2 = O2 / (Km_O2 + O2)
resp_capacity = f_O2 * mu_crit / max(mu, mu_crit)
ferm_fraction = 1.0 - min(resp_capacity, 1.0)

// Exponential smoothing for metabolic inertia
f_ferm_realized += alpha * (ferm_fraction - f_ferm_realized)
// where alpha = 1 - exp(-dt / tau_metabolic_switch)

mu_factor = (1.0 - ff) * aerobic_mu_factor + ff * anaerobic_mu_factor
Y_eff = yield_carbon * ((1-ff) * aerobic_yield_factor + ff * anaerobic_yield_factor)

// Fermentation acetate secretion
acetate_prod = ff * ferm_acetate_yield * carbon_consumed
```

**Files to modify**:
- `src/io/chem_environment_config.h` — Add to `OxygenConfig`:
  - `bool metabolic_switch_enabled = false`
  - `Real aerobic_mu_factor = 1.0`
  - `Real anaerobic_mu_factor = 0.55` (Varma & Palsson 1994: 0.43/0.78)
  - `Real aerobic_yield_factor = 1.0`
  - `Real anaerobic_yield_factor = 0.25` (0.128/0.524 gDW/g glucose)
  - `Real mu_crit = 9.7e-5` (0.35 h⁻¹ in 1/s)
  - `Real ferm_acetate_yield = 0.67`
  - `Real tau_metabolic_switch = 3600.0` (1 hour)
- `src/fixes/fix_metabolism.cpp` — Replace O₂ boost block (lines 433-440) with phase-plane logic; add yield adjustment in carbon uptake; add fermentation acetate secretion
- `src/gpu/agent_update_kernel.cu` — Same replacement on GPU; add `f_ferm_realized` per-agent buffer
- `src/gpu/agent_pool_gpu.h/.cpp` — Add `DeviceBuffer<double> d_f_ferm_realized_`
- `src/core/agent.h` — Add `Real f_ferm_realized = 0.0` to agent struct
- `src/io/input_parser.cpp` — Parse all new oxygen config keys

**Critical**: The existing O₂ boost code path MUST be preserved when
`metabolic_switch_enabled = false`. Use an if/else to select legacy vs new.

**Test**: With switch on and `epithelial_conc = 5e-6`, agents near the lumen
should have higher `f_ferm_realized` and produce acetate. Agents at high growth
rate should produce overflow acetate even with O₂ present.

### Change 3: Acid-product growth inhibition

**What**: Monod-type growth inhibition from accumulated acetate.

**Key equation**:
```
inhibition = acetate / (Ki_acetate + acetate)
mu *= (1.0 - acid_inhibition_max * inhibition)
```

**Files to modify**:
- `src/io/chem_environment_config.h` — Add to `OxygenConfig` or new struct:
  - `bool acid_inhibition_enabled = false`
  - `Real Ki_acetate = 20.0` (mol/m³)
  - `Real acid_inhibition_max = 0.8`
- `src/fixes/fix_metabolism.cpp` — Insert inhibition multiplier after O₂ switch, before maintenance subtraction
- `src/gpu/agent_update_kernel.cu` — Same on GPU
- `src/io/input_parser.cpp` — Parse new keys

**Test**: Inject high acetate concentration, verify growth reduction follows
inhibition curve.

## Parameter summary

All Palsson-anchored parameters (Varma & Palsson 1994,
doi:10.1128/aem.60.10.3724-3731.1994):

| Parameter | Value | Origin |
|---|---|---|
| anaerobic_mu_factor | 0.55 | W3110: 0.43/0.78 h⁻¹ |
| anaerobic_yield_factor | 0.25 | W3110: 0.128/0.524 gDW/g |
| mu_crit | 9.7e-5 /s (0.35 h⁻¹) | Acetate overflow onset |
| tau_metabolic_switch | 3600 s | DynamicME proteome reallocation |
| ferm_acetate_yield | 0.67 | PFL mixed-acid stoichiometry |
| epithelial_conc (for runs) | 5e-6 mol/m³ | Mucus-surface O₂, not tissue |

## Validation config

```yaml
oxygen.enabled: true
oxygen.epithelial_conc: 5.0e-6
oxygen.Km: 3.0e-7
oxygen.metabolic_switch_enabled: true
oxygen.aerobic_mu_factor: 1.0
oxygen.anaerobic_mu_factor: 0.55
oxygen.aerobic_yield_factor: 1.0
oxygen.anaerobic_yield_factor: 0.25
oxygen.mu_crit: 9.7e-5
oxygen.ferm_acetate_yield: 0.67
oxygen.tau_metabolic_switch: 3600.0
oxygen.acid_inhibition_enabled: true
oxygen.Ki_acetate: 20.0
oxygen.acid_inhibition_max: 0.8
vbf.agent_carbon_coupling: 1.0e-16
```
