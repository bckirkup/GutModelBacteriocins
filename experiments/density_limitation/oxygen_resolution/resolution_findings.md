# Funded oxygen respiration in delivery mode is a discretization artifact

Measurement only. Nothing here is fitted or calibrated, and no production parameter
is changed by this experiment.

## Question

Every aerobic result in this line of work (#319, #323, #325, #327, #329) was run at
`grid_dx = 2 um` and `bio_dt = 60 s`, and all of them reported the same headline:
growth respiration is funded at approximately zero. The open question was whether
that is a property of the model's uptake physics or of the grid and step size it was
measured on.

The delivery-mode agent sink is a first-order voxel sink whose rate is capped by
`demand / (C_pre * V_cell * dt)`. Both the cap and the voxel content depend
explicitly on `V_cell` and `dt`, so if funding is transport-limited in the physical
sense the funded fraction should be insensitive to them, and if funding is set by the
cap it should move with them.

## Arms

Base config is the flora-off oxygen arm (`oxygen.vbf_sink = 0`,
`anaerobic_maintenance_factor = 1.0`, `uptake_limit = "delivery"`,
`oxygen.delivery_uptake_enabled = true`, `oxygen.respiration_driver = "funded"`,
Robin epithelial boundary at the sourced `k = 1.2e-6 m/s`, imposed z-gradient off).
The flora sink is off so that nothing competes with the agents for oxygen: whatever
the agents cannot fund here is not attributable to competition.

Only `grid_dx`, `bio_dt`, and the output cadence differ between arms. Domain is
96 x 96 x 300 um throughout, so every `grid_dx` divides all three extents exactly.
Horizon is 3600 s, which is long enough for the emergent oxygen profile to
equilibrate and short enough to run four arms on a CPU.

| Arm | `grid_dx` | `bio_dt` | Grid | Wall |
|---|---:|---:|---:|---:|
| `res2_dt60` | 2 um | 60 s | 48 x 48 x 150 | 66 s |
| `res4_dt60` | 4 um | 60 s | 24 x 24 x 75 | 14 s |
| `res6_dt60` | 6 um | 60 s | 16 x 16 x 50 | 6 s |
| `res2_dt10` | 2 um | 10 s | 48 x 48 x 150 | 426 s |

A `res2_dt1` arm was started and abandoned at t = 505 s on wall-time grounds; its
partial output is preserved but is not used here.

## Result

Cumulative values at t = 3600 s. "Normalized funded" is cumulative funded growth
oxygen divided by integrated live-agent-seconds, which is the comparison that
survives a change of `bio_dt`.

| Arm | Funded fraction | Normalized funded (mol/agent/s) | Steps with any funding | Maintenance shortfall (mol) | Mean fermentation fraction | N at 3600 s |
|---|---:|---:|---:|---:|---:|---:|
| `res2_dt60` | 0.0063 | 4.57e-21 | 23 of 60 | 4.50e-14 | 0.630 | 40 |
| `res4_dt60` | 0.4693 | 4.91e-19 | 59 of 60 | 9.12e-16 | 0.334 | 39 |
| `res6_dt60` | 0.6635 | 8.07e-19 | 60 of 60 | 3.24e-16 | 0.193 | 28 |
| `res2_dt10` | 0.6341 | 4.68e-19 | 357 of 360 | 1.23e-15 | 0.207 | 34 |

Oxygen and carbon reaction clips are exactly zero in every arm, and every arm
terminated at `horizon_reached`. Mean oxygen concentration is 4.46e-2 to 4.69e-2
mol/m^3 in all four arms, i.e. the field itself is essentially identical across the
sweep: the oxygen is equally present, and only the agents' access to it changes.

Three statements follow directly.

1. **Funded respiration spans two orders of magnitude across arms whose physics is
   identical.** Normalized funded oxygen rises 107x from 2 um to 4 um and 177x from
   2 um to 6 um at fixed `bio_dt`, and 102x from `bio_dt = 60 s` to `10 s` at fixed
   `grid_dx`. The "growth respiration is funded at zero" result is a property of the
   2 um / 60 s discretization, not of the model's oxygen supply.

2. **The scheme does not converge under refinement, and it moves in opposite
   directions on the two axes.** Refining the grid drives funding down; refining the
   step drives it up. A transport-limited quantity would be insensitive to both. This
   is the signature of the rate cap `demand / (C_pre * V_cell * dt)`, which scales as
   `1 / (V_cell * dt)`.

3. **The metabolic switch inherits the artifact.** The mean realized fermentation
   fraction at 3600 s is 0.63 at 2 um / 60 s and 0.19 at 6 um / 60 s, from the same
   founders in the same oxygen field. Maintenance shortfall differs by up to 140x on
   the same comparison. So the facultative-metabolism trajectory reported in earlier
   arms is also resolution-dependent, in the direction of spuriously excess
   fermentation.

## What this invalidates and what it leaves standing

Invalidated as physical statements, though the runs themselves are unchanged:

- "Growth respiration is funded at exactly zero" (#325, #327, #329). It is funded at
  47-66 percent of demand in the same configuration at 4-6 um or at a 10 s step.
- Any conclusion that the agents cannot reach the oxygen standing in their voxel,
  including the seam flagged in #327 and Edison's reading of it as a funding-rule
  strictness problem. The funding rule is not too strict; the cap is resolution-bound.

Left standing:

- The flora attribution from #329: with the flora off, funding is still incomplete at
  every resolution, so the flora was never the limiting consumer. This measurement
  strengthens that conclusion rather than weakening it.
- The sqrt(3) operator-split property of the emergent oxygen profile (#329), which is
  a separate discretization statement about the sink, not about uptake.
- Carbon-limited population decline: N falls from 80 to 28-40 in every arm within an
  hour, and that ordering does not track the funded-oxygen ordering.

## What it implies for the uptake model

The comparison Edison proposed (analytic Sherwood versus delivery, flora off) cannot
be run as a configuration change: `oxygen.respiration_driver = "funded"` is only
accepted with `uptake_limit = "delivery"`, so there is no analytic-Sherwood oxygen
funding path to compare against. This sweep is the substitute test, and it answers the
same question: at 2 um and 60 s a single cell's analytic Sherwood ceiling
(4*pi*D_eff*r*C*dt = 8.1e-18 mol) is 69 times its voxel's entire oxygen content
(1.2e-19 mol), so the per-voxel cap cannot express the analytic rate at this
resolution no matter what the field holds.

The fix therefore has to make funded uptake independent of `V_cell` and `dt` up to the
transport the field can actually supply, rather than tuning the cap. Removing the cap
alone is not sufficient: the uncapped conductance sink draws 2.6 times the analytic
ceiling through the implicit solve, and returning that excess to the agent's own voxel
enriches it 620x above background, which the gradient regression test in
`tests/test_uptake_limit.cpp` catches. Any candidate fix must keep removal and funding
equal without violating the locality of the removal.

## Reproduce

    python3 experiments/density_limitation/oxygen_resolution/make_resolution_arms.py
    python3 experiments/density_limitation/oxygen_resolution/analyze_resolution_arms.py

Outputs read for this document:

- `o2fund-res/resolution_timeseries.csv` (per-summary-record series, all arms)
- `o2fund-res/resolution_metrics.csv` (checkpoint table)
- `o2fund-res/<arm>/output.h5`, `<arm>/input.json`, `<arm>/full.log`
- `o2fund-res/<arm>_oxygen_profile_last.csv`
