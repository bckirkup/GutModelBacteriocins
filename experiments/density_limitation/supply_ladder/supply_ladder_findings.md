# Supply ladder: the delivery brake sets a carrying capacity, and it is supply-set

Fifteen CPU arms, 24 h horizon, delivery-limited uptake, one seed (1001), probe
geometry `48x48x150` at `grid_dx = 2 um` (`2.765e-6 mL`, so one agent is
`3.6e5 cells/mL` and the shipped dysbiosis guard at `1e8 cells/mL` is ~276
agents). Axes: `supply` = multiplier on `carbon.epithelial_flux` (1.00x =
literature directly-measured host flux), `flora` = multiplier on
`vbf_carbon_sink_vmax`. Every arm reached the horizon with **zero cumulative
reaction clip**, so closure holds throughout.

## Capacity is monotone in available carbon, on both axes

`N_plateau` is the mean population over the final quarter; `block` is the
carbon nutrient blocking fraction (the agents' share of realized carbon
removal).

| supply | flora 1.0x | flora 0.5x | flora 0.0x | block (1.0x / 0.5x / 0.0x) |
|---:|---:|---:|---:|---|
| 1x | 30 | 202 | 588 | 0.010 / 0.315 / 1.000 |
| 2x | 305 | 685 | 1012 | 0.279 / 0.641 / 1.000 |
| 4x | 1197 | 1470 | 1788 | 0.628 / 0.814 / 1.000 |
| 8x | 2738 | 3052 | 3437 | 0.811 / 0.905 / 1.000 |
| 16x | 5914 | 6494 | 6638 | 0.904 / 0.952 / 1.000 |

This is the first monotone dose-response the model has produced for density.
Three things follow.

**Above ~2x supply, capacity is linear in supply.** At full flora, 4x -> 16x
gives 1197 -> 5914, i.e. 4.9x population for 4x carbon. The brake is therefore
behaving as a supply-limited carrying capacity rather than a grid artifact: a
fixed carbon budget divided by a fixed per-cell cost.

**Below ~2x, capacity collapses faster than supply.** 2x -> 1x at full flora
drops 305 -> 30, a 10x fall for a 2x cut. That is the maintenance floor: at
1.00x J_dir maintenance consumes essentially the whole delivered budget, so
capacity is set by what is left over after maintenance, not by supply itself.
This is why every mechanism probed at 1.00x returned a null — including
`agent_carbon_coupling` (#313). That regime was the wrong place to look.

**Flora competition only matters near the maintenance floor.** Removing the
flora entirely raises capacity 20x at 1x supply (30 -> 588) but only 12% at 16x
(5914 -> 6638). The background flora is a decisive competitor exactly when host
supply is marginal, and nearly irrelevant when it is abundant.

The blocking fraction tracks all of this monotonically (0.010 at N=30 to 0.904
at N=5914), which is what makes it usable as the campaign's reported observable:
one number that says how much of the carbon the modelled population is actually
winning, comparable across geometry and supply in a way no input constant is.

## Two caveats that bound the result

**The high-supply arms have not converged.** Final-quarter / previous-quarter
slope is 0.90-1.03 for 1x-4x (plateaued) but 0.70-0.76 at 8x-16x, which
overshoot hard and are still falling: 16x peaks at N=14739 near t=12 h and is
down to 6317 by 24 h. Those two rows are upper bounds on capacity, not
capacities; they need a 48-72 h horizon before they can be quoted.

**Every arm at or above 2x supply plateaus above the dysbiosis guard.** The
shipped guard is `1e8 cells/mL` ~ 276 agents here; only `1x` at full flora
(N=30, `1.1e7 cells/mL`) sits below it. These runs raised
`dysbiosis_threshold` to `1e10` deliberately so capacity could be measured
instead of the guard, but it means the guard-safe band and the density-limited
regime barely overlap at this geometry: the only guard-compliant capacity is
the near-starvation one. That has to be settled before the RPS campaign picks
an operating point, and it is a question about geometry and flux, not about the
guard's implementation.

## Funding is maintenance-limited at capacity everywhere

Funded fraction at the plateau is 0.15-0.29 across all fifteen arms and the
maintenance shortfall is positive in every one. So capacity is reached with the
population chronically unable to fund its full demand — consistent with the
delivery brake, and confirming that the shortfall is a physiological state
rather than the wiring fault the old §6 pass criterion treated it as.

## Provenance gap found while analysing

`/run_provenance/resolved_config` does not emit the epithelial carbon boundary
keys (`carbon.epithelial_flux`, `carbon.epithelial_boundary`), so an output
file cannot state the carbon supply it was run at — for a campaign whose main
axis is supply, that is an auditing hole. The analysis script falls back to the
arm's `input.json` and flags it; the emitter needs the missing keys.
