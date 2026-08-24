# Oxygen competition: the flora is not the thief, and the shell is 1.7x too deep

Thirteen CPU arms on the implicit rate-proportional VBF oxygen sink (PR #328,
`465e3d1`), Robin epithelium at the sourced permeability
(`k_transfer = 1.2e-6 m/s`, `C_tissue = 5.0e-2 mol/m3`), imposed z-profile off,
anatomy-derived placement, 80 founders, `bio_dt = 60 s`, 24 h horizon.

Generator: `make_competition_arms.py`. Configs, logs, HDF5 and the metric CSVs
live outside the repository (campaign workspace), not in git.

## What the sink fix bought

Zero oxygen reaction clip in every arm, including `vbf_sink = 0.84 /s` and
`3.36 /s` at a 60 s step, where the old explicit `reac*dt` path had
`k*dt = 50` and `204` and clipped nearly all of the demanded removal. The
"98% of the sink's demand was clipped" result reported in the placement matrix
was an integration artifact, not a statement about the field. High-rate flora
oxygen consumption is measurable for the first time.

## An oxic shell exists, and it is 1.7x deeper than the parameter implies

With the imposed profile off, the field develops a real exponential z-profile
for the first time:

| `vbf_sink` (/s) | surface (mol/m3) | measured 1/e depth | `sqrt(D/k)` | ratio |
|---:|---:|---:|---:|---:|
| 0 | 4.96e-2 | none (uniform) | - | - |
| 1e-3 | 3.99e-2 | none (uniform) | 1449 um | - |
| 0.84 | 2.27e-3 | 86 um | 50 um | 1.72 |
| 3.36 | 1.17e-3 | 44 um | 25 um | 1.76 |

The ratio is the same in both arms and it is `sqrt(3)`. The cause is in the
transport splitting: the delivery sink is applied at `dt/3` in each of the three
directional sweeps (`transport_slab_z` and its replicated twin pass
`context.dt / 3.0`). Total removal over the three sweeps is correct, and the
per-step ledger closes, but the x and y sweeps see a field that is uniform in x
and y, so they contribute only a spatially flat decay; the *shape* is set by the
z sweep balancing full `D` against `k/3`, giving `sqrt(3D/k)`. At `k*dt` of 50
to 200 the operators do not commute and the splitting biases the profile, not
the mass.

Consequence: the sink needed for a 25 um shell in this solver is about
`10 /s`, not `3.36 /s`. Targeting a depth by `sqrt(D/k)` will overshoot by 73%
until the reaction split is made consistent (symmetric splitting, or a smaller
transport sub-step for the sink term).

## The flora is not the thief

The attribution arm settles Edison's Call 3 against the flora hypothesis. With
`vbf_sink = 0` — no flora oxygen consumption at all, agents alone in a 49.6 uM
column:

- cumulative funded growth respiration / demand: **0.0054**
- funded growth respiration per step: **exactly zero at every step**
- oxygen maintenance is never fully paid; the shortfall is 30-40% of demand

Agent intake tracks the boundary influx step for step and decays with it
(2.6e-14 mol at step 1, 2.0e-15 at step 60, 7.6e-17 at the last step) while the
field holds 49.6 uM at every depth, unchanged. So the ledger inference from the
Robin measurement was right: in delivery mode an agent can draw the marginal
influx into its voxel, not the oxygen standing in it.

Raising the sink to make a gradient makes the agent share worse, as predicted:
cumulative agent oxygen 4.08e-13 mol at `1e-3`, 1.23e-13 at `0.84`, 6.09e-14 at
`3.36`, while the flora takes 1.32e-11 to 1.35e-11 — essentially the whole
boundary delivery, which saturates at its ceiling
`k*(C_tissue - C_surface)*A*dt = 3.24e-14 mol/step`. That ceiling is 6.7x the
population's total oxygen maintenance demand, so the boundary is not the binding
constraint either.

## The maintenance ladder has no survival threshold to bracket

`1.0 / 4.1 / 8.3 / 15.0` at the shipped flora sink, and `1.0 / 4.1` crossed with
the four sink values. The factor does what it is supposed to: carbon maintenance
paid and shortfall scale with it monotonically (at step 60, shortfall 0 at 1.0,
1.3e-17 at 4.1, 8.0e-17 at 8.3, 2.1e-16 at 15). But:

- `N_peak = N_initial = 80` in **every** arm: not one net division, ever
- fermentation fraction reaches 0.9986-0.9990 in every arm
- every arm ends `population_stop` between step 395 and 420

So `anaerobic_maintenance_factor = 15` is not what killed the population in the
earlier campaigns, and setting it to 1.0 does not rescue anything. Growth was
never funded — neither respiratory (0 at every step) nor fermentative (no
divisions at 1.0, where fermentative maintenance carbon costs the same as
aerobic). The maintenance factor is a real modelling choice but it is not the
lever, and there is no threshold in `1.0-15` to cross against the oxygen axis.

## Two void arms and a provenance defect

The `ceiling` arm (one founder, flora off, maintenance 1.0) terminated at step 0
with `population_stop`: the population-stop guard fires on a single-agent run
before any biology executes. The per-cell funding ceiling is therefore still
unmeasured; it needs the guard disabled, not a different founder count.

The first ladder batch was stamped `b4328727` — the merge of PR #326 — because
`build-serial/CMakeCache.txt` carried a pinned `GUTIBM_GIT_SHA` override from an
earlier session, and `cmake/generate_git_sha.cmake` prefers the override over
`git rev-parse HEAD`. Those four arms were discarded and re-run after clearing
the override. The stamp is now verified against `HEAD` before a campaign runs.
No other build directory in the tree carried a non-empty override.

## The open question, and the measurement that settles it

Growth respiration is zero with the flora off, with the boundary unsaturated,
and with 49.6 uM standing in the field. What remains is the uptake model itself.
For one cell at that concentration the quasi-steady diffusive delivery is

    4*pi*D*r*C*dt = 4*pi*2.1e-9*0.5e-6*4.96e-2*60 = 3.9e-14 mol/step

against a per-cell oxygen maintenance demand of 6e-17 mol/step — a factor of
650 of headroom. Delivery mode gives that same cell about 2.3e-17 mol/step, some
three orders below the analytic delivery, because its per-voxel first-order sink
can only remove what the 2 um voxel holds (4e-19 mol at 49.6 uM) plus what the
discrete Laplacian resupplies within the step.

The discriminator is therefore `uptake_limit`, not oxygen: run the flora-off arm
under `sherwood` (analytic `4*pi*D_eff*r*C_local`) against `delivery`, and a
`grid_dx` sweep at 2/4/8 um under `delivery`. If `sherwood` funds respiration
and `delivery` does not, the aerobic result to date is a discretisation
statement about voxel-scale resupply rather than a biological one, and the
choice between the two uptake models becomes a scientific claim that has to be
made explicitly.
