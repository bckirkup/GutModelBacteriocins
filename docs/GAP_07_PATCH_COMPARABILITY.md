# GAP-07: What do 57 and 488 count?

**They are paired cell counts for one patch at the 12 h audit, not two estimates of carrying capacity.** In the [placement provenance][prov], patch 0 of `incub_ship_ly0.1_s20261027` has 488 live agents and 57 shadow colonists at 43,200 s. The recorded relative difference is 7.5614, measured against the shadow. Neither number is a 24 h unconstrained endpoint.

Source baseline: [main at 28242fd](https://github.com/bckirkup/GutModelBacteriocins/tree/28242fd01d164e66f5752b7436d1ec749fbcf2fd). The [placement report][placement] identifies build `778c2af`; the [carbon-ladder report][ladder] identifies build `a115447`. The former has committed per-run records; the latter's raw configs, analyzer and HDF5 are [not in the repository][ladder-provenance].

| What to compare | Ladder used to build PatchTable | 57: tabulated shadow | 488: live placement slot |
|---|---|---|---|
| Population counted | Reported **peak N** over each amplitude arm, not final N or births. One strain, 100 plasmid-free founders; arm A is single-seed. The report does not give the exact peak-sampling cadence or dead-cell filter. | **Snapshot:** length of the shadow colonist vector, initialized from live cells. Births follow a stochastic logistic rule based on count, patch type and volume, not strain-specific growth. | **Snapshot:** agents in the one live simulation whose state is not DEAD, across all lineages. Initially 30 founders (10 each C/S/R); not the segment stock, cumulative divisions, or a colony count. |
| Time window | 12 h requested; all five amplitude arms reached the horizon. The ladder's one early extinction was in [a different arm][ladder-stop]. A finite-window peak need not be an equilibrium. | At each 6 h audit the engine records both counts, then **resets the shadow from live**. At 6 h it recorded 56 live / 49 shadow. The 12 h shadow thus started its window with **56**, not 49 or 30. | Same audit clock and starting population at 6 h; the live slot retains its fields, agents and history. Its 12 h value is not the peak of a separately sampled ladder arm. |
| Geometry and carbon | 200 × 200 × 100 µm (4.0 × 10⁻⁶ mL), 4 µm grid; founders in the first 20 µm of depth. Oxygen off, uptake **delivery**, 10 µm support, not zero radius. | Same nominal 60 µm cube as live (2.16 × 10⁻⁷ mL); no resolved chemical field. At default crypt supply 2×, the table uses the ladder's **1,063** peak: modeled K = 1,063 × 2.16 × 10⁻⁷ / (4.0 × 10⁻⁶) = **57.402** cells. K limits the birth rate; 57 is a realized integer census, not a rounded lookup. | One of ten crypt-type patches. The live slot applies the 2× carbon boundary, disables crypt flow geometry and oxygen, and inherits the incubator's damped losses. Config + build defaults imply a 2 µm grid and uptake **none** (no per-agent delivery cap), not the ladder's delivery mode. |
| Termination | Dysbiosis threshold raised to 10¹⁰ cells/mL; the five peaks were not censored by that guard. | Steps with the active live patch; after the 12 h reset the 18 h record is 488 / 488. That equality is not evidence of predictive recovery. | First sampled halt at 12 h. The segment logged 1,440 steps / 86,400 s, the requested duration, with `termination_cause=incomplete_unknown`. The [outer loop][outer] can exit at the time limit without changing that label. `audit_halt=false`: the discrepancy did **not** terminate the run. |

Table sources: [ladder setup and peaks][ladder], [run config][config], [crypt defaults and audit cadence][defaults], [stored peaks and growth rates][peaks], [volume scaling][capacity], [tabulated births][births], [live-slot overrides][overrides], [grid][grid] and [uptake defaults][uptake], [audit count and reset][audit]. The ladder's founder depth band is not slab-decomposed chemistry: the [parser guard][slab] rejects positive delivery radius with the latter, not the former.

**Halt evidence.** At 30 min, `patch_cells=302` [sums all ten patches][stocks]; the [live lineage census][lineages] sums to 11. At 10.5 and 11 h the live census is 285 and 321, respectively. The [segment series][series] first records one halted crypt at 12 h (zero at 11.5 h). This is consistent with the [per-patch guard][guard] and its [defaults][guard-defaults]: above 10 × 30 founders for one hour, not a threshold crossed by this patch at 30 min. The series records neither the triggering step nor the individual halt reason. Inactive patches [skip both live and shadow dynamics][skip].

**Inference and smallest next count check.** The shadow underpredicts the stock in this particular reset window. Geometry scaling, population composition, uptake mode and halt censoring prevent attribution to any one cause; the live equilibrium is unknown. For a second observation, inspect the same arm's [next archived seed][next-seed] at its 12 h audit: use the 6 h live count as the starting stock and check the [halt flag][next-series] before comparing the two stocks. No new table or run is needed for that question.

## What the carbon numbers establish

| Question | Evidence on main | Limit |
|---|---|---|
| What does `yield_carbon = 0.5` multiply? | A positive biomass increment (kg) times the [realized carbon cost][cost] gives mol modeled substrate. The [non-delivery reaction][reaction] subtracts that amount divided by voxel volume and step duration (mol/m³/s) from the reaction field. [Delivery][delivery] prices demand first and [funds growth from removal][funding]. | Units are **mol substrate per kg modeled biomass**, the reciprocal direction to g dry biomass per g substrate. The multiplication does not identify a glucose molecule or a mole of carbon atoms. |
| Which biomass? | The [member][biomass] says dry kg, but [initialization][mass] equates biomass to sphere volume × [1,100 kg/m³][density]; [growth][growth] sets mechanical mass equal to biomass again. No dry fraction enters these assignments. | The [documented comparison][comparison] assumes 0.524 g dry biomass/g glucose and 180 g/mol. That implies 1/(0.524 × 0.180) ≈ 10.6 mol glucose/kg **dry**, or 3.18 mol glucose/kg **wet** only at 30% dry mass. One mol glucose is six C-mol. Neither a wet/dry convention nor this substrate mapping is established by a numerically balanced field; the reported 6.4× is conditional, not a measured yield error. |
| What does `enforce_reaction_clip` test? | The [clamp][clamp] adds the deficit needed to lift a negative concentration to zero, in mol. The [optional gate][gate] checks recorded carbon clip per step against tolerance fraction × max(growth demand + paid maintenance + maintenance shortfall, machine epsilon). [Defaults][closure-defaults]: off, fraction zero. A strict excess in [Simulation::run][run] stops with `ClosureViolation`. | The flag tests recorded clipping; it does not enable positivity, test every species, require demand to be met, or compare initial with final carbon. The ladder [reports zero clips][ladder-clip]; the placement slot [calls step directly][live-step], not this run-loop gate. Its [closed ledger][ledger] balances cells, not mol. |

**Smallest next balance check.** For one existing carbon-field interval with before/after stocks and fluxes, evaluate the [implemented residual][residual] in mol: initial + signed boundary + gradient source + VBF source − agent uptake − paid maintenance − VBF sink + reaction clip − final. Report the clip separately from the residual; a zero residual can include mass supplied by the clamp. The placement census lacks these terms. Even field closure would not demonstrate elemental balance into biomass and products: that requires a declared substrate and biomass basis first.

[prov]: ../analysis/spec13_placement/out/incub_ship_ly0.1_s20261027_prov.json
[config]: ../analysis/spec13_placement/out/cfg_incub_ship_ly0.1_s20261027.json
[placement]: ../analysis/spec13_placement/REPORT.md#arms
[ladder]: CARBON_LADDER_CAMPAIGN.md#measured-ladder-28-arms-serial-local
[ladder-provenance]: CARBON_LADDER_CAMPAIGN.md#provenance
[ladder-stop]: CARBON_LADDER_CAMPAIGN.md#L232-L234
[ladder-clip]: CARBON_LADDER_CAMPAIGN.md#L199-L206
[defaults]: ../src/layer2/layer2_config.h#L38-L67
[peaks]: ../src/layer2/patch_table.h#L35-L42
[capacity]: ../src/layer2/patch_table.cpp#L33-L40
[births]: ../src/layer2/segment.cpp#L405-L425
[overrides]: ../src/layer2/segment.cpp#L507-L545
[grid]: ../src/core/domain.h#L33-L37
[uptake]: ../src/fixes/fix_metabolism.h#L70-L76
[slab]: ../src/io/input_parser.cpp#L622-L631
[audit]: ../src/layer2/segment.cpp#L874-L904
[outer]: ../src/layer2/segment_io.cpp#L105-L133
[stocks]: ../src/layer2/segment.cpp#L923-L940
[lineages]: ../analysis/spec13_placement/out/incub_ship_ly0.1_s20261027_lineage.csv
[series]: ../analysis/spec13_placement/out/incub_ship_ly0.1_s20261027_ts.csv
[guard]: ../src/layer2/segment.cpp#L793-L813
[guard-defaults]: ../src/layer2/layer2_config.h#L93-L98
[skip]: ../src/layer2/segment.cpp#L659-L696
[next-seed]: ../analysis/spec13_placement/out/incub_ship_ly0.1_s20261028_prov.json
[next-series]: ../analysis/spec13_placement/out/incub_ship_ly0.1_s20261028_ts.csv
[cost]: ../src/fixes/fix_metabolism.cpp#L915-L922
[reaction]: ../src/fixes/fix_metabolism.cpp#L1149-L1167
[delivery]: ../src/fixes/fix_metabolism.cpp#L409-L432
[funding]: ../src/fixes/fix_metabolism.cpp#L536-L564
[biomass]: ../src/core/agent.h#L46-L55
[mass]: ../src/core/agent.cpp#L33-L45
[density]: ../src/core/types.h#L30-L31
[growth]: ../src/fixes/fix_metabolism.cpp#L1124-L1139
[comparison]: PARAMETERS.md#L558-L561
[clamp]: ../src/fields/chemical_field.cpp#L2654-L2665
[gate]: ../src/core/simulation.cpp#L1422-L1468
[closure-defaults]: ../src/io/input_parser.h#L92-L97
[run]: ../src/core/simulation.cpp#L1531-L1541
[live-step]: ../src/layer2/segment.cpp#L624-L636
[ledger]: ../src/layer2/segment.cpp#L962-L971
[residual]: ../src/fields/chemical_field.cpp#L3694-L3723
