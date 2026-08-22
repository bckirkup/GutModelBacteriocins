# Spec 12 (density-dependent growth limitation) — audit against the code

Reference: `gutibm_spec12_density_limitation.md` (Benjamin, 2026-08-13).
Code state: main @ f1a7aee (post #301). Geometry/parameters taken from the
Sherwood campaign arm `sherwood_f018_*_s1001.json`, i.e. the configuration the
"no brake" result was measured in.

Domain constants used throughout:

| quantity | value |
|---|---|
| domain | 96 x 96 x 300 um, `grid_dx = 2 um` |
| voxels | 345,600, each 8.0e-18 m^3 |
| domain volume | 2.765e-12 m^3 |
| epithelial area | 9.216e-9 m^2 |
| delivery (0.18x J_dir) | 1.936e-9 mol/m^2/s = **1.784e-17 mol/s** total |
| guard (1e8 cells/mL) | **276 agents** in this domain |
| bio step | 60 s |

## Verdict

The diagnosis in the Motivation section is right, and two of the three
mechanisms are the right *kind* of mechanism. But as written, Spec 12 would
not produce a basin, and one of its changes has a sign error that would make
the bloom faster. Four numeric conflicts, then the amendment I recommend.

---

## Conflict 1 (blocking): the anaerobic yield factor has the wrong sign in this codebase

Spec 12 Change 2:

```cpp
Y_eff = yield_carbon * (anaerobic_yield_factor + (aerobic - anaerobic) * f_O2);
// anaerobic_yield_factor = 0.25
```

with the stated intent "lower yield -> more carbon per biomass -> faster local
carbon depletion". That intent is correct for a *biomass-per-substrate* yield
(gDW/g glucose). But `yield_carbon` in this code is the **inverse**: substrate
per biomass. From <ref>src/fixes/fix_metabolism.cpp</ref>:

```cpp
const Real demanded = d_biomass * yield;                  // kg * yield -> mol
Real delta_c = d_biomass * cfg_.yield_carbon / (cell_vol * dt);
```

`d_biomass` is kg, the result is mol, so `yield_carbon` is **mol carbon per kg
biomass**. Multiplying it by 0.25 makes a fermenting cell **4x cheaper**, not
4x more expensive: Change 2 as literally specified *removes* carbon pressure
and accelerates the bloom, on top of the 0.35x growth-rate reduction.

The factor must be inverted (anaerobic ~= 4.0 relative to aerobic 1.0), or the
config renamed to `anaerobic_carbon_cost_factor` so the direction is
unambiguous. This is worth fixing in the spec regardless of what we implement,
because the same inversion is easy to reintroduce.

## Conflict 2 (largest lever, not in the spec): `yield_carbon = 0.5` is ~7x too cheap

`yield_carbon = 0.5 mol/kg` (<ref>src/fixes/fix_metabolism.h:43</ref>) against
literature: E. coli aerobic yield ~0.5 gDW/g glucose, dry fraction ~0.3 of the
wet mass the model tracks (`CELL_DENSITY_DEFAULT = 1100 kg/m^3`, so `biomass`
is wet), glucose 180 g/mol:

```
1000 g wet x 0.3 gDW/g / 0.5 gDW/gGlc / 180 g/mol = 3.33 mol glucose / kg wet
```

So the model's cells are **~6.7x cheaper per unit biomass than E. coli**. That
single number is most of the reason agents booked only 4-15% of epithelial
carbon while the VBF ate 139-217% of it, and therefore most of the reason the
mean-field carbon brake sits above the guard rather than below it:

| yield | supply / full-growth demand |
|---|---|
| 0.5 (current) | 88 agents |
| 3.33 (literature, aerobic) | **13 agents** |

At 0.18x J_dir the guard is 276 agents. With the current yield the carbon
brake engages *above* the guard; with the literature yield it engages an order
of magnitude *below* it. This is a calibration correction, not a new mechanism,
and it is a bigger effect on density limitation than any of Spec 12's three
changes. It also has to land *before* the O2 yield factor of Change 2, or that
factor is re-anchoring to a wrong baseline.

## Conflict 3: `agent_carbon_coupling = 1e-16 mol/s/agent` is ~500x its stated value

The spec calls 1e-16 "~5% of direct agent carbon demand at full growth". Direct
demand at full growth is `m x mu_max x Y = 8.14e-16 kg x 5e-4 /s x 0.5 mol/kg =`
**2.04e-19 mol/s** (agent mass measured from the Sherwood final dump). So:

- 1e-16 is **491x** full-growth demand, not 0.05x;
- 1e-16 x 276 agents = 2.8e-14 mol/s, i.e. **1550x the entire delivery budget**
  (1.78e-17 mol/s) — carbon would hard-zero within a step or two of the first
  few divisions, giving extinction rather than capacity;
- the whole recommended sensitivity range 1e-17 ... 1e-15 sits above per-agent
  demand. 5% of demand is ~1e-20.

A range that brackets the guard is ~1e-20 ... 1e-19 mol/s/agent (capacity
`N* = supply / coupling` = 1780 ... 178 agents, guard at 276).

## Conflict 4: O2 saturates ~2 orders of magnitude below the guard, so Change 2 is three constants

Per-cell respiration already exists and is already density-coupled
(<ref>src/io/chem_environment_config.h</ref>): `o2_use = q_consumption * mu +
q_maintenance` = 1e-14 x 1e-4 + 1e-18 ~= **2e-18 mol/s/cell**, which matches
Alexeeva's 1.5e-18. But diffusive supply through 80 um of mucus is
`D A C / L`:

| boundary [O2] | total O2 supply | cells supported at full respiration |
|---|---|---|
| 5 uM (Spec 12) | 1.21e-18 mol/s | **0.6** |
| 55 uM (default) | 1.33e-17 mol/s | **6.7** |

(The background first-order sink is irrelevant by comparison: 1.4e-20 mol/s at
5 uM, and its depletion length `sqrt(D/k)` = 1.45 mm >> the 300 um domain.)

So `f_O2 -> 0` above a handful of agents at *either* boundary value, and stays
there for the whole run. `mu_factor` is then pinned at 0.35, `Y_eff` at its
anaerobic value and `ferm_fraction` at 1.0 — three constants. That is the same
binary failure mode as the VBF sink sweep (0x decisive, 0.5x/1x/2x
indistinguishable): it moves where net growth crosses zero, it does not bend
the curve back down. Spec 12's own literature note anticipates this ("no
published model achieves stable E. coli carrying capacity from O2 alone").

Consequence: the O2 switch is worth having for realism (it makes the mucus
anaerobic, which it should be, and it is the physically honest way to get the
anaerobic yield and fermentation flux), but it must not be counted as one of
the "three independent braking forces". It is a constant multiplier.

## Conflict 5: acid inhibition would be driven by the VBF, not by density; and no local acid peak can exist

Two independent problems.

**(a) The background dominates the field by ~500x.** `AcetateConfig` has
`vbf_production = 1e-3 mol/m^3/s` against `vbf_consumption = 2e-4`, i.e. a net
background source of **8e-4 mol/m^3/s** everywhere. Agent fermentation at
N = 325 and the measured per-agent carbon draw contributes **1.5e-6
mol/m^3/s** — 0.2% of it. So `[acetate]` reaching `Ki_acetate = 20 mol/m^3` is
a property of the VBF and the washout rate, on a clock (~7 h at 8e-4 with no
removal), not a property of agent density. Change 3 as specified is a
time-dependent uniform tax.

**(b) There are no local gradients to exploit anyway.** Diffusion length over
one 60 s bio step: acetate 268 um, O2 355 um, carbon 200 um — all comparable to
or larger than the 300 um domain depth. Every soluble species is well-mixed
within a step. Combined with occupancy (276 agents over 345,600 voxels = 0.08%,
and one agent alone in an 8 fL voxel is a *local* 1.25e11 cells/mL, 1250x the
guard), the per-voxel histogram `n_agents_in_voxel / V_cell` is not a density:
it is 0 or 1 almost everywhere and 1250x the guard where it is 1. Every
per-voxel coupling in Spec 12 therefore acts mean-field, through domain-mean
concentration vs total N. That is a real carrying-capacity mechanism (chemostat
style), but the spatial story in the spec is not what will be executing.

Also, Change 2's fermentation secretion is not a like-for-like replacement of
the existing `acetate_overflow_threshold` path: that path is
`overflow_rate * biomass / cell_vol` = ~1e-13 mol/m^3/s per agent, i.e. **10
orders of magnitude** smaller than the specified fermentative flux (~1.6e-3
mol/m^3/s per occupied voxel). Also worth noting `overflow_rate` and
`scavenge_rate` are documented "mol/s/cell" but used as mol/(s.kg).

---

## What the code actually lacks: maintenance consumes no carbon

<ref>src/fixes/fix_metabolism.cpp</ref> ends the growth-rate calculation with

```cpp
mu -= cfg_.maintenance_rate;   // 1e-5 /s
```

and `grow_agent` removes carbon strictly in proportion to `d_biomass`. So
maintenance is a *growth-rate* tax that costs the cell **nothing** in substrate:
a cell at mu = 0 eats zero carbon and persists for free. Oxygen already has the
correct Pirt form (`q_consumption * mu + q_maintenance`, with the comment
explicitly noting the maintenance term is what makes the field track density).
Carbon does not.

This is precisely the linear-in-N sink that Change 1 is trying to synthesise
with an invented coefficient — and it has a measured one. Pirt
non-growth-associated glucose maintenance m_s ~= 0.04 g glucose/gDW/h gives

```
0.04 x (0.3 x 8.14e-16 kg x 1e3 g/kg) / 180 / 3600 = 1.5e-20 mol/s/cell
```

and therefore a carrying capacity, with no phenomenology, of

```
N* = supply / m_s = 1.78e-17 / 1.5e-20 = ~1180 agents = 4.3e8 cells/mL   (0.18x J_dir)
```

That is 4.3x above the guard — i.e. *in reach*, and graded: `N*` is linear in
delivery flux, so the flux bracket that found no basin becomes a bracket over
capacity, and a guard-safe capacity sits near 0.04-0.05x J_dir (or at the upper
end of the literature m_s range, which anaerobic maintenance justifies). Compare
Sherwood, which was 86x out and could not be brought into reach without moving a
measured constant. Anaerobic maintenance is *higher* than aerobic, so Change 2's
O2 switch feeds this term rather than competing with it.

## Recommended amendment

Keep the spec's structure, re-anchor its numbers, and add the term it is
missing:

1. **Change 0 (new, do first): Pirt maintenance carbon.**
   `metabolism.maintenance_carbon` (mol/s/kg biomass, default 0 = current
   behaviour), removed from the carbon field independently of `d_biomass`, and
   raised under fermentation. Anchor 1.5e-20 mol/s/cell-equivalent. This is the
   brake; it is measured; it is the same Pirt form the O2 field already uses.
2. **Yield re-anchor (Conflict 2).** `yield_carbon` 0.5 -> ~3.3 mol/kg as an
   explicit factor, reported both ways, because it moves the mean-field capacity
   across the guard on its own and every other coefficient here is defined
   relative to it.
3. **Change 2 as specified, with the yield factor inverted** (Conflict 1) and
   *not* counted as a brake: it buys an anaerobic mucus, the higher anaerobic
   carbon cost, and the fermentative flux that Change 3 needs.
4. **Change 3 with `Ki_acetate` re-anchored to what the field can reach**, and
   with the VBF's own acetate production either subtracted from the inhibiting
   pool or acknowledged as a constant tax. Lump formate into acetate as total
   acid load (stoichiometry ~1.0, not 0.67); do not add a formate field.
5. **Change 1 (`agent_carbon_coupling`) demoted to optional**, at 1e-20 ...
   1e-19 mol/s/agent if kept at all — with Change 0 in place it is the same
   mathematical object with a fitted coefficient instead of a measured one, so
   my preference is to drop it and keep the VBF sink density-independent.
6. Keep the diagnostics list as written; drop the Thiele modulus per z-layer in
   favour of a single well-mixedness check, since the diffusion lengths above
   say the answer in advance.

Everything stays off by default, as the spec requires.
