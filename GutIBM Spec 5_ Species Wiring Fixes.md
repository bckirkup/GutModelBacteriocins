# Devin Task: Spec 5 — Fix Chemical Species Wiring Gaps

## Context

A systematic audit of every chemical species in GutIBM found three wiring
gaps that break the model's carrying capacity dynamics. Without these fixes,
carbon accumulates without bound, oxygen never depletes despite agent
consumption, and B12 drains to zero with no replenishment.

These are small, targeted fixes — not new features.

---

## 1. Add VBF Carbon Sink

### Problem
The VBF produces carbon via mucin_liberation but never consumes it. The
`nutrient_sink` config parameter only acts on iron (vbf.cpp line 62:
`reac(iron) -= nutrient_sink * conc(iron)`). Carbon is only consumed by
E. coli agents, which is negligible at small population sizes. Result:
carbon concentration increases monotonically from 0.001 to 0.85 mol/m3
over 24 hours in the null control run.

In reality, the anaerobic majority (10^11 cells/mL, represented by the VBF)
consumes the vast majority of mucin-derived monosaccharides. E. coli
scavenges only the leftovers ("Restaurant Hypothesis").

### Fix
In `vbf.cpp`, add a carbon consumption term in `apply_carbon_source()`:

```cpp
void apply_carbon_source(ChemicalField& chem, Int cell, const VbfCellContext& ctx) {
  // ... existing mucin liberation code (lines 42-53) ...

  // VBF carbon consumption: anaerobes eat most of the liberated carbon
  if (ctx.idx.carbon >= 0) {
    const Real carbon_conc = chem.conc(ctx.idx.carbon, cell);
    chem.reac(ctx.idx.carbon, cell) -= ctx.cfg.carbon_sink * carbon_conc;
  }
}
```

Add `carbon_sink` to `VBFConfig`:
```cpp
Real carbon_sink = 5e-3;  // 1/s first-order consumption rate
```

Calibration: at steady state, carbon source = carbon sink.
`mucin_liberation * z_weight = carbon_sink * [C_ss]`
`[C_ss] = mucin_liberation * z_weight / carbon_sink`
At z=0: `5e-5 / 5e-3 = 0.01 mol/m3` — a thin residual carbon pool.
At z=50um (z_weight ~0.14): `5e-5 * 0.14 / 5e-3 = 0.0014 mol/m3`.

This gives Monod carbon limitation (`Km_carbon` is typically ~1e-4 mol/m3
for glucose): `monod_carbon = 0.01 / (1e-4 + 0.01) = 0.99` near wall,
`0.0014 / (1e-4 + 0.0014) = 0.93` at z=50um. Still mostly carbon-replete.

For stronger limitation, use `carbon_sink = 5e-2`:
`[C_ss] = 5e-5 / 5e-2 = 0.001` → `monod = 0.001/(1e-4 + 0.001) = 0.91`.
Still not very limiting. The Km_carbon may need to be higher, or the
carbon_sink rate needs to be much larger to create real carbon scarcity.

Alternative: use a Monod-style saturating sink (more realistic):
```cpp
chem.reac(ctx.idx.carbon, cell) -=
    ctx.cfg.carbon_sink_vmax * carbon_conc / (ctx.cfg.carbon_sink_km + carbon_conc);
```
With `carbon_sink_vmax = 4.9e-5` (just under liberation) and
`carbon_sink_km = 1e-3`, this ensures the VBF consumes nearly all carbon at
high concentration but leaves a thin residual at low concentration.

Add to `VBFConfig`:
```cpp
Real carbon_sink_vmax = 4.9e-5;   // mol/m3/s max VBF carbon consumption
Real carbon_sink_km = 1e-3;        // mol/m3 half-saturation
```

Add parser entries for `vbf_carbon_sink_vmax` and `vbf_carbon_sink_km`.

### Files to modify
- `src/fields/vbf.h` — add carbon_sink_vmax, carbon_sink_km to VBFConfig
- `src/fields/vbf.cpp` — add carbon consumption in apply_carbon_source()
- `src/io/input_parser.cpp` — parse vbf_carbon_sink_vmax, vbf_carbon_sink_km

---

## 2. Add O2 Consumption by E. coli Agents

### Problem
O2 is read for the aerobic growth boost (fix_metabolism.cpp lines 217-222)
but agents never consume it. No `chem.reac(i_o2, cell) -= ...` exists.
The O2 field stays constant at 2.75e-6 mol/m3 regardless of agent density.

In reality, E. coli consumes O2 via cytochrome oxidases (cydAB, cyoABCDE).
At high cell density near the epithelium, O2 depletion is the primary
mechanism that limits the aerobic niche and forces cells further from the
wall to fermentative growth.

### Fix
In `fix_metabolism.cpp`, in the nutrient consumption section (after line 338),
add O2 consumption:

```cpp
// O2 consumption by agent
if (const auto& o2cfg = sim_.config().chem_env.oxygen; o2cfg.enabled) {
    Int i_o2 = chem.find(species::OXYGEN);
    if (i_o2 >= 0 && cell >= 0) {
        // Consumption proportional to growth rate and biomass
        // (aerobic respiration scales with metabolic activity)
        Real o2_consumed = o2cfg.q_consumption * agent.mu_realized * agent.biomass;
        chem.reac(i_o2, cell) -= o2_consumed / cell_vol;
    }
}
```

Note: `q_consumption` has units mol_O2/(growth_rate * kg_biomass * s).
The existing default is 1e-14 mol/s/cell, but the formula above scales
consumption with mu_realized (faster growth = more O2 used) and biomass.

For a simpler implementation that matches the existing parameter semantics:
```cpp
Real o2_consumed = o2cfg.q_consumption;  // mol/s per cell, flat rate
chem.reac(i_o2, cell) -= o2_consumed / cell_vol;
```

### Self-limiting feedback
With O2 consumption wired in:
- Few cells near epithelium: O2 stays high → aerobic boost → fast growth
- Many cells: O2 depleted → boost drops → growth slows → population stabilizes
- Cells further from wall never see O2 → always fermentative rate

This is the key density-dependent feedback for carrying capacity.

### Files to modify
- `src/fixes/fix_metabolism.cpp` — add O2 consumption reaction term

---

## 3. Add B12 Replenishment Source

### Problem
B12 has an initial concentration (1e-9 mol/m3) but no source term. It is
only consumed by E. coli agents (line 338). Over time, B12 depletes to zero,
forcing all cells into the MetE pathway regardless of BtuB expression.
This is a drift artifact, not biology.

In reality, B12 (cobalamin) in the gut comes from:
- Dietary intake (absorbed primarily in the ileum, but some reaches colon)
- Bacterial synthesis (Propionibacterium, Lactobacillus reuteri, some 
  Firmicutes produce B12)
- The concentrations are low but maintained at steady state

### Fix
Add a B12 source term to the VBF (representing production by the anaerobic
community):

In `vbf.cpp`, add to `apply_vbf_at_cell()`:
```cpp
void apply_b12_source(ChemicalField& chem, Int cell, const VbfCellContext& ctx) {
    Int i_b12 = chem.find(species::B12);
    if (i_b12 < 0) return;
    // VBF produces B12 at a low constant rate, maintaining steady-state
    chem.reac(i_b12, cell) += ctx.cfg.b12_production;
}
```

Add to `VBFConfig`:
```cpp
Real b12_production = 1e-13;  // mol/m3/s — maintains ~1 nM steady state
```

Calibration: at steady state, B12 production = B12 consumption by agents.
B12 consumption per cell ≈ yield_b12 * mu_realized * biomass / cell_vol.
With 100 cells in 25 nL, consumption ≈ 1e-13 mol/m3/s. Matching production
to typical consumption maintains homeostatic levels.

Add parser entry for `vbf_b12_production`.

### Files to modify
- `src/fields/vbf.h` — add b12_production to VBFConfig
- `src/fields/vbf.cpp` — add apply_b12_source(), call from apply_vbf_at_cell()
- `src/io/input_parser.cpp` — parse vbf_b12_production

---

## 4. Add Dysbiosis Threshold Check (bonus, low priority)

Add a configurable density threshold that halts the simulation if exceeded:

```cpp
// In Simulation::step(), after agent count update:
if (cfg_.dysbiosis_threshold > 0) {
    Real density = global_agent_count / domain_volume_mL;
    if (density > cfg_.dysbiosis_threshold) {
        std::cerr << "DYSBIOSIS THRESHOLD EXCEEDED: "
                  << density << " cells/mL > " << cfg_.dysbiosis_threshold
                  << " — halting simulation\n";
        break;  // or set a flag
    }
}
```

Config:
```cpp
Real dysbiosis_threshold = 1e8;  // cells/mL, 0 = disabled
```

This ensures the model stays in the homeostatic regime it's designed for.

---

## Tests to Add

- `test_carbon_sink`: With VBF carbon sink enabled, carbon reaches steady 
  state instead of accumulating. Verify [C] < 0.1 mol/m3 after 1000 steps
  with liberation enabled.
- `test_o2_consumption`: With 100 agents in a small domain, O2 at their 
  z-position should be lower than the initial/boundary value. More agents
  → lower O2.
- `test_b12_steady_state`: With VBF b12_production enabled, B12 should 
  stabilize near the initial value rather than draining to zero.
- `test_dysbiosis_halt`: Simulation with explosive growth should halt when
  density exceeds threshold.

## Priority
1. Carbon sink (biggest impact on carrying capacity)
2. O2 consumption (density-dependent feedback)
3. B12 replenishment (prevents drift artifact)
4. Dysbiosis threshold (safety net)
