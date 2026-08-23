# Funded O2 x placement x O2 consumption: respiration is refused by density, not by depth

Eight arms, `make_placement_matrix.py`, revision `c01d6e6` (PR #324 branch),
CPU only, seed 1001, 86400 s horizon, 1.00x carbon supply, flora 1.00x,
`oxygen.epithelial_conc = 5.5e-2 mol/m^3` (~42 mmHg) in every arm.
Outputs and per-arm ledgers:
`/home/ubuntu/gutibm-campaign/o2placement-1787521946-citable/` (`matrix_metrics.csv`,
per-arm `output.h5`, `oxygen_z_profile_last.csv`).

All eight arms began with 80 agents; the anatomic arms respected the 20 um
floor and the 150 um outer extent (initial range 20.55-148.58 um), the z-slab
arms sat in 0.70-49.46 um, and no arm placed a founder in a crypt.

## The result

**Growth respiration is funded at exactly zero in all eight arms.** Not small:
`agent_uptake` for oxygen is 0.0 at every step, at both placements, under both
drivers, across an 840x range of background O2 consumption. The mechanism is
not depth. `calculate_delivery_funding` pays maintenance first and funds growth
only from the remainder, and oxygen maintenance is never fully paid:

| arm | step | N | O2 boundary in | O2 VBF sink | O2 maint paid | O2 maint short | O2 growth funded |
|---|---:|---:|---:|---:|---:|---:|---:|
| funded_zslab_shipped | 1 | 80 | 3.323e-15 | 7.601e-16 | 2.678e-15 | 2.122e-15 | 0 |
| funded_zslab_shipped | 300 | 46 | 9.080e-16 | 7.607e-16 | 1.605e-16 | 2.720e-15 | 0 |
| funded_anatomic_shipped | 1 | 80 | 2.725e-15 | 7.601e-16 | 2.097e-15 | 2.703e-15 | 0 |
| ambient_zslab_shipped | 1 | 80 | 3.310e-15 | 7.601e-16 | 2.665e-15 | 2.135e-15 | 0 |

(mol per 60 s step.) Since maintenance is unpaid, `growth_funded` is
structurally 0 and `mu_realized *= growth_fraction` sets realized growth to
exactly 0 for every agent in the funded arms.

**The binding quantity is density, not depth.** Per-cell oxygen maintenance is
`q_maintenance * dt = 6.0e-17` mol/step, so 80 cells demand 4.8e-15 mol/step
against a delivered 3.3e-15 mol/step at the transient peak and 9.1e-16 mol/step
once the field equilibrates -- of which the background flora takes 7.6e-16.
The steady-state remainder funds the maintenance of roughly **two cells**. The
flora, not the epithelium, is what E. coli is competing with for oxygen here.

80 agents in this 96x96x300 um domain is 2.9e7 cells/mL, ~300x the
culture-based healthy target of 1e4-1e5 CFU/mL adopted in
`docs/SPEC13_MULTISCALE.md`. At the healthy target this domain would contain
0.03-0.3 cells: it is smaller than the volume per cell. Every result we have
produced from this geometry, including the supply ladder, has therefore been a
bloom-density result, and the aerobic axis was being asked to supply a
population two to three orders too dense to be fed.

## Placement is exonerated

Anatomic placement is *deeper* than the band we had been running, and it
changes nothing about funding: growth respiration is zero at 0-50 um and at
20-150 um alike. "Founders sit outside the oxic shell" is dead as an
explanation -- it was never depth. What placement does change is survival:
both anatomic arms halt on the population floor at 24600 s, while
`ambient_zslab` reaches the horizon at N=21-36, consistent with the ladder's
~28-agent carbon capacity at 1.00x supply. Depth costs carbon, because carbon
also enters at z=0.

## The size of the accounting error #323 removed

Same placement, same field, driver swapped:

| placement | ambient mean f | funded mean f |
|---|---:|---:|
| z_slab 0-50 um | 0.163 | 0.910 |
| anatomic 20-150 um | 0.476 | 0.855 |

The ambient driver reported 84% of growth as respiratory at the shallow band
while the field funded none of it. `ambient_zslab_shipped` survives to the
horizon on exactly that unfunded respiration; under the funded driver the same
configuration collapses. The population that reached the horizon in every
earlier oxygen arm was living on oxygen the model never delivered.

## The 25 um penetration depth is imposed, and consumption cannot reproduce it

The oxygen z-profile at step 60 is indistinguishable across arms whose
background O2 sink differs by 840x:

```
z (um)      0        10        20        30        50       100       150      200
shipped   5.500e-02 3.541e-02 2.373e-02 1.589e-02 7.119e-03 9.300e-04 1.072e-04 0
```

That is the hard-coded `k_z_lambda = 25e-6` exponential, essentially unmodified
by population or by consumption. Raising the sink to the value for which
sqrt(D/k) = 50 um does not deepen or shallow anything: it produces a 6.3e-13
mol/step reaction clip against a 6.4e-13 mol/step demand, i.e. **98% of the
consumption is refused because the oxygen is not there to consume**. A
consumption-derived shallow oxic shell is unreachable at this boundary supply.

So the answer to "is 25 um correct?" is that it is not a prediction at all. The
shallow shell is an imposed initial and reference profile, and the free
parameter that actually sets the aerobic axis is the epithelial oxygen **flux**
(measured here at 6e-9 mol/m^2/s), which nothing in the model or the docs
currently constrains. That, not lambda and not placement, is the next thing to
source.

## Facultative survival was not observed, and the fermentative penalty is why

Your prediction was that funded O2 should make the population fermentative
rather than dead. It made them fermentative (f -> 0.86-0.91) and then dead. The
chain is: unfunded respiration -> f -> 1 -> `anaerobic_maintenance_factor = 15`
and `anaerobic_carbon_cost_factor = 4.1` -> carbon maintenance shortfall ->
carbon growth funding also drops to exactly 0 under the same maintenance-first
rule -> population floor. A 15x maintenance penalty for fermentative growth is
the number to question: it means the model treats the fermentative mode as
nearly unviable, which is not what a facultative anaerobe in the colon does.

## What this does not license

No claim about aerobic or fermentative biology should be made from these arms
either. They establish where the model breaks, not what E. coli does:

1. the density is ~300x too high for the target organism and compartment;
2. the epithelial O2 flux is unsourced;
3. `anaerobic_maintenance_factor` is doing more work than any evidence behind
   it supports;
4. `q_maintenance = 1e-18` mol/s/cell is ~15x above a back-of-envelope
   E. coli basal O2 demand (0.5-1 mmol/gDW/h at ~3e-13 g dry weight) and enters
   every number above linearly.

Nothing here was tuned, and no arm was re-run to obtain a preferred outcome.
