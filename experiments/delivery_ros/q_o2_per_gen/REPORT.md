# Q_O2_per_gen on merged main (post-#338), low density

Binary built from `ffab55b` (merge of #338). Six arms, 12 h requested horizon,
one strain, no plasmids, `oxygen.k_ROS = 0`, 10 µm delivery support (shipped
default). Only `grid_dx` and founder count vary. Scripts: `prepare.py`,
`analyze.py`; raw metrics in `metrics.json`.

| arm | dx (µm) | N0 | divisions | funded O2/cell/s | T_gen (div) | **Q_O2_per_gen** | rationing factor | infeasible | clips | min O2 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Q_n2_res2 | 2 | 2 | 30 | 2.00e-18 | 2.27 h | **1.631e-14** | 1.0 | 0 | 0 | 0 |
| Q_n2_res4 | 4 | 2 | 30 | 2.00e-18 | 2.27 h | **1.633e-14** | 1.0 | 0 | 0 | 0 |
| Q_n2_res6 | 6 | 2 | 30 | 2.00e-18 | 2.27 h | **1.637e-14** | 1.0 | 0 | 0 | 0 |
| Q_n8_res4 | 4 | 8 | 65 | 1.99e-18 | 2.48 h | 1.780e-14 | 1.0 | 0 | 0 | 0 |
| Q_n20_res4 | 4 | 20 | 104 | 1.96e-18 | 3.27 h | 2.313e-14 | 1.0 | 0 | 0 | 0 |
| Q_n80_res4 | 4 | 80 | 383 | 1.78e-18 | 3.33 h | 2.127e-14 | **0.5** | 0 | 0 | 0 |

`Q_O2_per_gen = (agent_uptake_cumulative + maintenance_cumulative) / divisions`,
i.e. funded growth plus funded maintenance oxygen per completed generation. The
demand ledger records growth demand only (`add_uptake_demand(oxygen,
growth_demand)` in `fix_metabolism.cpp`), which is why the reported funded
fraction is ~2: maintenance is funded at about the same magnitude as growth and
is not in the denominator. It is not a double count.

## Findings

1. **Grid invariance achieved.** Q spread across 2/4/6 µm is **1.004×**, against
   **3.3×** on the pre-#335/#336/#338 code. The rationing factor is exactly 1.0
   at N ≤ 20, no infeasible steps, no reaction clips, and the oxygen field never
   goes negative — so this is a per-cell physical quantity, not a discretisation
   or accounting residue.

2. **The divisor is only defined uncrowded.** Q rises 1.63 → 1.78 → 2.31e-14 from
   N0 = 2 → 20, and at N0 = 80 the delivery rationing factor drops to 0.5: supply
   is being shared. Take **Q = 1.63e-14 mol/generation** (≈ 16 fmol) as the
   uncrowded value. Crowded cells respire less per unit time, so a hazard priced
   per funded mol automatically lyses them less — which is the coupling we want.

3. **Resulting lysis coefficients** (`k = -ln(1-P)/Q`, absolute funded-flux form):

   | target lysis/generation | k_ROS_funded (mol⁻¹) |
   |---:|---:|
   | 1% | 6.2e11 |
   | 2% (central) | 1.24e12 |
   | 5% | 3.1e12 |

   Edison's illustrative 2e13 assumed 1 fmol/generation; the model says 16 fmol.
   These apply to the **absolute** funded-flux driver; #332 currently implements a
   *specific* (per-biomass) flux, so that normalisation change has to land before
   the numbers can be used as-is.

4. **Every arm emptied by boundary export, not mortality.** All losses are
   `outflow_boundary` (32 of 32 at N0=2; 462 at N0=80); zero lysis, zero
   starvation, zero washout-trap deaths, with ROS off. N0=80 peaks at 281 cells
   and decays to 22. With T_gen ≈ 2.3 h against a 5400 s radial mucus turnover, a
   single patch is not self-sustaining — so there is no patch-scale carrying
   capacity to quote without reseeding. This is the washout axis Spec 13 exists to
   close, and it now dominates in the absence of the retired ROS term.

5. **Open item before the carbon ladder:** mean `mu_realized` decays monotonically
   from 5.5e-4 to 7.5e-5 /s (7×) over 12 h even at 2–16 cells, with the oxygen
   rationing factor pinned at 1.0. So the decay is not oxygen delivery; the
   candidate is local carbon/mucin supply. T_gen is therefore not stationary and
   the tabulated Q averages over a declining growth rate — worth explaining before
   the #314 carbon capacity re-run rather than after.

## Correction

The `5.5e-4` value above is the step-0 `mu_max` placeholder, not a realized
growth rate. The realized decline is approximately 2.0×, from `1.62e-4` to
`8.2e-5 /s`, and is spatial; see §3 of `docs/DELIVERY_ROS_CAMPAIGN.md`.
