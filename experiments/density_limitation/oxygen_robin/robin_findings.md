# Robin epithelium: respiration is funded for the first time, and only as a transient

Revision measured: `b4328727d2d3a98a712a989bd56006c55d2efa69` (main, after #326).
Four arms, CPU only, 24 h horizon requested; all four terminated early on
`population_stop`. Generator: `make_robin_arms.py`. Raw artifacts:
`/home/ubuntu/gutibm-campaign/o2robin-1787530388` (`robin_metrics.csv`,
`probe_validation.csv`, per-arm `oxygen_z_profile_last.csv`, run logs).

Every arm uses the `funded` respiration driver, delivery-limited carbon and
oxygen uptake, anatomy-derived placement (#324), and the shipped
`oxygen.vbf_sink = 1e-3 1/s`. `anaerobic_maintenance_factor` was left at its
shipped `15`, deliberately, so that the boundary change is not confounded by an
open parameter question.

## Result 1: the Dirichlet control reproduces #325 exactly

| arm | funded growth O2 (mol) | demanded (mol) | surface (mol/m^3) | 1/e depth |
|---|---|---|---|---|
| `dirichlet` | **0.0** | 2.02e-13 | 0.055 (imposed) | 24.0 um |

Zero funded growth oxygen, and a 24.0 um `1/e` depth against the hard-coded
25 um `z_gradient_lambda`. The control does what it was built to do: nothing in
#326 moved the Dirichlet path, so the other three arms are attributable to the
boundary condition.

## Result 2: Robin funds growth respiration, and it decays to zero

| arm | k (m/s) | funded growth O2 (mol) | / demanded | funded at step 1 / 5 / 60 / 300 / last |
|---|---|---|---|---|
| `robin` | 1.2e-6 | 6.01e-17 | 3.5e-4 | 0 / 0 / 0 / 0 / 0 |
| `robin_k10` | 1.2e-5 | 1.85e-14 | 0.112 | 4.18e-15 / 1.45e-15 / 0 / 0 / 0 |
| `robin_sparse` | 1.2e-6 | 2.79e-17 | 1.8e-3 | 0 / 0 / 0 / 0 / 0 |

This is the first nonzero funded growth respiration this model has produced, so
the prediction that a depletable epithelium would change the funding was right.
But the per-step columns say what the cumulative column hides: at 10x the
sourced permeability, 11% of respiratory demand is funded **in the opening
minutes and then nothing**, and at the sourced permeability the cumulative
3.5e-4 is likewise assembled from a handful of early steps. This is exactly the
transient signature the PR-#326 initial condition was fixed to avoid faking —
the mucus fills from an empty state, and while it is filling there is oxygen no
one has claimed yet. Steady-state funded growth respiration is zero at both
permeabilities.

## Result 3: the flora takes 99% of epithelial delivery, and density does not change that

At steady state, per step:

| arm | boundary in | VBF sink | remainder | maintenance paid |
|---|---|---|---|---|
| `robin` | 6.612e-15 | 6.548e-15 | 6.5e-17 | 6.559e-17 |
| `robin_k10` | 8.050e-15 | 7.979e-15 | 7.1e-17 | 7.248e-17 |
| `robin_sparse` | 6.612e-15 | 6.548e-15 | 6.4e-17 | 6.580e-17 |

The remainder after the flora equals the maintenance the agents were paid, to
within a percent, in every arm. The first-order background sink absorbs 99.0-99.1%
of everything the epithelium delivers, and it does so *regardless of k*: raising
the transfer coefficient tenfold raised steady-state agent-available oxygen from
6.5e-17 to 7.2e-17 mol/step, about 11%, because the flora consumed the rest.

Then the density axis, which was the load-bearing hypothesis going in:
`robin_sparse` has a tenth of the founders and the **same** absolute remainder,
so its per-cell oxygen is identical to `robin`'s, not tenfold better. At step 60
maintenance paid is 1.559e-15 for 80 cells and 1.578e-16 for 8 cells — the same
per cell to three figures. So for oxygen specifically, agent density is not the
binding constraint: the flora sink is. That is a correction to the framing I
carried out of #325, where the two-cell arithmetic made density look like the
whole story. Density remains the correct account of the *carbon* capacity and of
what a single patch can represent; it is not the account of why respiration is
unfunded.

## Result 4: the oxic shell did not become emergent, it disappeared

`oxygen_z_profile_last.csv`, mol/m^3:

| z (um) | `dirichlet` | `robin` | `robin_k10` |
|---|---|---|---|
| 0 | 0.0550 | 0.03974 | 0.04875 |
| 50 | 0.00714 | 0.03947 | 0.04842 |
| 150 | 1.17e-4 | 0.03907 | 0.04793 |
| 298 | 0.0 | 0.03883 | 0.04764 |

The Robin arms have no `1/e` crossing anywhere in the domain: oxygen falls 2.3%
across 300 um of mucus. The emergent surface concentration is 79% of the
tissue-side reservoir at the sourced k and 97% at 10x, so the epithelium is
barely resisting delivery. Retiring the imposed 25 um profile did not produce a
shallower or a deeper shell; it produced a **uniformly oxic mucus column at
~40 uM**, which is not what colonic mucus does.

This is not a surprise on the arithmetic, and it is the same number from #325
seen from the other side: at the shipped sink, `sqrt(D/k) = 1449 um` against a
300 um domain. The domain is five times thinner than the oxygen diffusion length,
so it cannot support an oxygen gradient under *any* boundary condition. The 25 um
shell was never emergent and cannot be made emergent by fixing the boundary; a
gradient of that depth needs the O2 consumption rate raised by roughly three
orders (~0.02 1/s to decay across the domain, ~0.84 1/s for a 50 um shell), which
is a parameter question about colonic mucus O2 consumption and not a plumbing one.

## Open seam: 40 uM standing everywhere, and respiration still unfunded

These two facts are in the same run: the field holds ~40 uM at every depth, and
funded growth respiration is zero at steady state. The identity in Result 3 —
agent-available oxygen tracks `boundary_in - VBF_sink` per step, independent of
agent count and of the standing content — suggests that under delivery-limited
uptake agents are funded from the marginal epithelial influx net of the flora
rather than from the oxygen standing in their own voxel. That is an inference
from the accounting, not something I have read out of the funding code, and it is
the next thing to check: if it is right, the funding rule is stricter than the
physics (a cell sitting in 40 uM cannot respire), and if it is wrong then the
standing 40 uM is unavailable for some other reason that also needs finding.

## All four arms died as fermenters

Every arm: N 80 (or 8) -> 1, `population_stop` at step 343-410, i.e. under 7 h,
with fermentation fraction rising to 0.999 and cumulative carbon clip exactly 0.
No arm reached the 720-step mark, which is why those columns are `nan`. This is
the `anaerobic_maintenance_factor = 15` behaviour identified in #325, held fixed
here on purpose: fermentation is a death sentence in the current parameterization,
so no arm can be read as a statement about long-run persistence. The maintenance
factor has to be settled before any of these arms is re-run for survival.

## What this changes

- A depletable epithelium was necessary and is not sufficient. #326 was worth
  doing — the Dirichlet face could not be depleted and its shell was imposed —
  but the boundary was not what kept respiration at zero.
- The binding constraint on agent oxygen is the background flora sink, not
  agent density and not the boundary condition. `oxygen.vbf_sink` is now the
  parameter to source, and it is doing two contradictory jobs at once: at 1e-3
  1/s it consumes 99% of delivery while being 58x too weak to shape a gradient.
- The 25 um oxic shell should be retired as a claim, not relocated. This
  geometry cannot produce one at the shipped consumption rate.
- Nothing here is a survival result. `anaerobic_maintenance_factor` is unsettled
  and every arm terminated on it.
