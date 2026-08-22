# agent_carbon_coupling probe: the voxel-local density sink is inert at campaign geometry

Measurement run against the merged Spec 12 Change 1 (#312), CPU, delivery-limited
uptake, probe geometry `48x48x150` cells at `grid_dx = 2 um`, `carbon.epithelial_flux
= 1.00x J_dir`, seed 1001, 6 h horizon, 80 founders. Configs from
`make_coupling_configs.py`, analysis from `coupling_analysis.py`.

## Result

Coupling was swept over **seven orders of magnitude**, spanning the
demand-anchored band `{1e-21, 1e-20, 1e-19}` (2%-200% of measured per-agent
demand) and, as bounding arms, `1e-23` and Spec 12's own `1e-16`. Carbon
totals are cumulative over the run; `block` is the carbon nutrient blocking
fraction in the final window.

| arm | coupling (mol/s/agent) | vbf_sink (mol) | agent uptake (mol) | maintenance (mol) | block | N_end | termination |
|---|---:|---:|---:|---:|---:|---:|---|
| c0 | 0 | 2.33742999e-12 | 7.61534103e-15 | 5.08394990e-14 | 0.014211 | 47 | horizon_reached |
| x1e-23 | 1e-23 | 2.33743212e-12 | 7.61487260e-15 | 5.08378432e-14 | 0.014211 | 47 | horizon_reached |
| c1e21 | 1e-21 | 2.33740425e-12 | 7.65768524e-15 | 5.08233042e-14 | 0.014227 | 47 | horizon_reached |
| c1e20 | 1e-20 | 2.33740844e-12 | 7.65666191e-15 | 5.08201621e-14 | 0.014226 | 47 | horizon_reached |
| c1e19 | 1e-19 | 2.33740891e-12 | 7.65654741e-15 | 5.08198130e-14 | 0.014226 | 47 | horizon_reached |
| x1e-18 | 1e-18 | 2.33740896e-12 | 7.65653581e-15 | 5.08197777e-14 | 0.014226 | 47 | horizon_reached |
| x1e-16 | 1e-16 | 2.33740896e-12 | 7.65653454e-15 | 5.08197738e-14 | 0.014226 | 47 | horizon_reached |

Every arm gives the same population trajectory (80 -> 84 -> 47) and the same
final N. Blocking fraction moves by 0.1% relative between zero and nonzero
coupling and is then **flat from 1e-21 upward**: the mechanism saturates just
above zero and a further 1e5x has no effect. No arm halted on the closure gate,
and cumulative reaction clip is 0 everywhere.

This is not a wiring defect. `tests/test_vbf_accounting.cpp` shows a measurable
occupied-versus-empty carbon gap with live agents and exactly zero with dead or
ghost agents, and the resolved config in `/run_provenance/resolved_config`
carries the intended `vbf.agent_carbon_coupling` in every arm. The mechanism is
wired, active, and negligible.

## Why: the ceiling is voxel resupply, not the constant

The coupled sink raises `vmax` in the agent's own 8 fL voxel. Two scales defeat
it at this geometry:

- **Occupancy.** 80 agents in 345,600 voxels is 0.023% of the domain. The entire
  carbon content of every occupied voxel is `2.84e-3 mol/m^3 * 8e-18 m^3 * 80
  = 1.8e-18 mol`, against `6.5e-15 mol` of flora removal per 60 s step. Even
  total starvation of every occupied voxel is ~3e-4 of the flora's own draw.
- **Within-step resupply.** Carbon diffuses ~77 um in 60 s, i.e. ~38 voxels, so
  the voxel is refilled from its neighbourhood far faster than any local sink
  can drain it. Once local `vmax` exceeds the resupply rate the sink is
  supply-limited and the constant stops mattering. That threshold is
  `~1e-21 mol/s/agent` (local `vmax = 1.25e-4` versus background `5.5e-5`),
  which is where the flat response begins.

So `agent_carbon_coupling` is a saturating knob with a ceiling set by
discretization and diffusion. This is the same shape as the retired `sherwood`
cap and the pre-#308 clipped agent sink - a mechanism that is installed, booked
and incapable of binding - except that here the ceiling is physical rather than
a coefficient, so no value of the constant repairs it.

## What would make density competition bind (inference, not yet tested)

The sink needs a support comparable to the carbon diffusion length rather than
one voxel: agents within ~77 um raise `vmax` across that neighbourhood, so the
flora's extra draw competes for the same pool the agent is drawing from. Sizing
that version: adding one background `vmax` over `(77 um)^3 = 4.6e-13 m^3`
requires `2.5e-17 mol/s` of coupled strength in the region, so at the
demand-anchored `5e-20 mol/s/agent` it takes ~500 agents per diffusion volume,
about 1.8e9 cells/mL - an order of magnitude above the dysbiosis guard at 1e8
cells/mL (~276 agents in this domain). At guard-relevant densities a non-local
sink would need `~1e-18` to `1e-17 mol/s/agent`.

That partially reconciles Spec 12's `1e-16`: the dimensional objection to it
stands for a voxel-local sink at `grid_dx = 2 um`, but the remedy is a larger
support, not a smaller constant. A voxel-local sink small enough to be
dimensionally defensible is necessarily too small to matter.

## Two other things this run establishes

- **The population here is maintenance-dominated, not growth-limited by
  coupling.** Cumulative maintenance `5.1e-14 mol` versus growth uptake
  `7.6e-15 mol`, funded fraction 0.128, and N declines 80 -> 47 with coupling
  switched off entirely. Any density brake added at this supply is being added
  to a population that is already starving.
- **The aerobic-overflow axis is available but unused at this supply.** With
  `oxygen.mu_crit = 9.7e-5 /s` (#312) the population runs at `mu_avg ~ 6e-6 /s`
  throughout, so `mu_crit / max(mu, mu_crit) = 1` and the realized fermentation
  fraction that rises to 0.66 is entirely O2-driven. The mu axis is no longer
  structurally inert, but exercising it needs a supply regime where growth
  approaches 0.35 h^-1; the 1.00x J_dir probe is not that regime.

## Consequence for the campaign

Do not spend a coupling sweep as a first-class RPS axis at this geometry: it is
measurably a no-op from 1e-21 to 1e-16. Either implement the diffusion-length
support and re-probe, or drop density coupling and treat the delivery brake plus
maintenance as the density mechanism, reporting outcomes against nutrient
blocking fraction rather than against any coupling constant.
