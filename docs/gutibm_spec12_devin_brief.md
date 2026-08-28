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
