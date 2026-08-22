# Mechanism probe: #302 Pirt carbon maintenance and #304 O2 metabolic mode

Six arms, seed 1001, 24 h, guard on, `kd_b12_btuB` 1e-6 (colicin silent).
Image `gutibm@sha256:541a8444c0bf14fc6516972861b57586fd5bf3cb1dc04f62cc3a4d772084260a`,
source commit `a22e4e7191fd8037d6208d821805ce30314daff6` (post-#302/#304/#305).
Configs and analysis: `make_probe_configs.py`, `probe_analysis.py`.

## Headline

| arm | t_end (h) | N_end | halt | density (cells/mL) | mean f_ferm | carbon cost (mol/kg) |
|---|---:|---:|---:|---:|---:|---:|
| off_f018 | 24.0 | 160 | none | 5.79e7 | 0.000 | 1.61 |
| maint_f018 | 24.0 | 160 | none | 5.79e7 | 0.000 | 1.61 |
| full_f018 | 24.0 | 124 | none | 4.49e7 | 0.156 | 2.37 |
| full_f030 | 18.5 | 368 | dysbiosis | 1.33e8 | 0.271 | 1.20 |
| full_f060 | 8.7 | 325 | dysbiosis | 1.18e8 | 0.221 | 1.33 |
| full_f100 | 4.5 | 378 | dysbiosis | 1.37e8 | 0.550 | 0.85 |

## 1. Maintenance carbon is inert, and not because the coefficient is small

`maint_f018` is identical to `off_f018` in every reported quantity, including
the population trajectory hour by hour and cumulative growth uptake
(1.245e-13 mol in both). The maintenance ledger says why: over 1440 steps the
mechanism drew **6.999e-16 mol** total, and the shortfall counter fired
**1.218e5** times, i.e. on essentially every agent-step
(1.218e5 / 1440 = 85 = the mean population).

The cause is the denominator, not the rate. `charge_carbon_maintenance` clamps
each agent's draw to `carbon_maintenance::available(conc, cell_volume)` — the
*instantaneous* carbon standing in the agent's 2 µm voxel:

```
voxel stock       = 6.615e-4 mol/m^3 * 8.0e-18 m^3 = 5.29e-21 mol
per-agent demand  = 2.1e-5 * f_maint * 7.7e-16 kg * 60 s
                  = 9.7e-19 mol   (maint arm, f_maint = 1)
                  = 3.2e-18 mol   (full arm, f_maint = 3.2 at f_ferm 0.156)
ratio             = 184x  ... 614x  more demanded than the voxel holds
```

So maintenance drains its voxel to zero every step and collects ~0.2–0.5% of
its Pirt demand. Meanwhile growth uptake in the same step is **not** clamped
(`uptake_limit=none`, `uptake_limited_agents_cumulative = 0`), so growth books
carbon the voxel does not contain while maintenance may not. The two carbon
sinks are budgeted against different denominators.

The instantaneous voxel stock is the wrong budget at `bio_dt = 60 s`: carbon
diffuses ~77 µm in 60 s against a 2 µm cell, so the voxel is replenished
thousands of times within a single step. The diffusive delivery to one cell
over the same step is `4*pi*D_eff*r*C*dt` = 2.7e-17 mol, which is **8–28x
larger** than the maintenance demand — i.e. with the physically correct
denominator, maintenance is fully funded and the brake runs at its measured
Pirt coefficient rather than at 0.3% of it.

Also a reporting defect: `maintenance_shortfall_*` increments by 1.0 per
clamped agent (a count of agent-steps) while sitting in a ledger whose other
entries are mol. It reads as 1.218e5 mol of carbon.

## 2. The metabolic mode works, and my masking prediction was wrong

I predicted O2 supply would support ~1 agent's respiration, pinning
`f_ferm` at 1.0 and making the aerobic branch unobservable. Measured:
mean realized fermentation fraction **0.156** at 0.18x rising to **0.550** at
1.00x, spanning 0.018–0.959 across individual agents, with mean O2
4.586e-6 mol/m^3 — above the 1e-6 Monod Km. The phase plane is genuinely
occupied and both branches are live; density and growth rate, not O2 supply
alone, set the mode.

The mode has a real cost: at matched flux and matched time,
`full_f018` vs `off_f018` is 124 vs 160 agents with realized carbon cost
2.37 vs 1.61 mol/kg (+47%, the anaerobic factor 4.1 acting on a partly
fermentative population). That is a drag, not a brake.

## 3. Still no plateau, and delivery is still the driver

The 0.18x arms are still rising at 24 h (5.8e7 cells/mL, below the 1e8 guard,
division 0.15–0.17 vs outflow 0.086 per agent per hour). Every arm at 0.30x
and above trips the guard, monotonically sooner with flux: 18.5 h, 8.7 h,
4.5 h. So the ordering is unchanged from the pre-mechanism bracket, and no arm
exhibits a density-limited steady state.

Acetate is being produced and drained (net -8.5e-12 mol across the epithelial
boundary at 0.18x, -1.6e-12 at 1.00x), but grid output was off in this probe,
so acetate *concentration* — the quantity that acid inhibition reads — was not
measured. It needs a grid-enabled arm before acid inhibition can be assessed.

## Consequence

The brake was implemented correctly and then starved by the voxel clamp. Fix
the denominator (route the maintenance draw through the same `uptake_limit`
model growth uses) before re-bracketing delivery; the re-bracket done on the
current code would measure a mechanism running at 0.3% strength.
