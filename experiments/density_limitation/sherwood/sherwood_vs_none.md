# Sherwood uptake limit vs none — does a diffusive uptake cap create a carrying capacity below the dysbiosis guard?

Six arms, `uptake_limit ∈ {none, sherwood}` × `carbon.epithelial_flux ∈ {0.14, 0.18, 0.22} × J_dir`
(J_dir = 1.0756e-8 mol/m²/s measured), seed 1001, 168 h requested, dysbiosis guard on at
1e8 cells/mL, `kd_b12_btuB` = 1e-6, image
`sha256:9868a520…fd531` (main `f1a7aee`), job definition `gutibm-cuda-bracket:1`,
queue `gutibm-gpu-practice`. All six jobs succeeded on one attempt with no Spot
reclaims; durable checkpoints were validated against `latest.json` by URI, size and
SHA-256 (`checkpoint_validation.jsonl`).

## Result: zero divergence. The cap never binds, in any arm, at any step.

| flux | limit | reached | halt | peak density | final n | n1/n2 | funded/demanded | intervals with cap binding |
|---|---|---:|---|---:|---:|---|---:|---:|
| 0.14x | none | 43.9 h | guard | 1.175e8 | 325 | 325/0 | 1.000000 | — |
| 0.14x | sherwood | 43.9 h | guard | 1.175e8 | 325 | 325/0 | 1.000000 | 0 / 2636 |
| 0.18x | none | 32.9 h | guard | 1.454e8 | 402 | 212/190 | 1.000000 | — |
| 0.18x | sherwood | 32.9 h | guard | 1.454e8 | 402 | 212/190 | 1.000000 | 0 / 1976 |
| 0.22x | none | 25.1 h | guard | 1.042e8 | 288 | 162/126 | 1.000000 | — |
| 0.22x | sherwood | 25.1 h | guard | 1.042e8 | 288 | 162/126 | 1.000000 | 0 / 1506 |

Each `sherwood` arm reproduces its paired `none` arm exactly: identical guard time,
identical final population, identical per-type composition, and 0.00% divergence in
μ, divisions per agent-hour and mean density over the final 24 h. The `uptake_limit`
key is confirmed ingested — each run's own resolved config records `none` or
`sherwood` as authored — so this is the model's answer, not a configuration no-op:
`limit_fraction` returned 1.0 on every agent on every one of ~2600 steps.

Guard times reproduce the delivery bracket's kd1e6 column run-for-run (43.9 / 32.9 /
25.1 h), which is the expected consequence of a cap that never fires.

## Why it cannot bind: the slack is concentration-independent, not an ordering

Correcting an arithmetic error in the first version of this note, which compared the
cap at a lowered concentration against demand measured at the *run's* concentration
and concluded the cap would bind below ~6e-6 mol/m³. It would not: demand falls with
C as well, so the comparison has to be made as a ratio.

Cap per agent per step (`uptake_limit.h`):

```
allowed(C) = 4π·D_eff·r·C·dt
```

Demand per agent per step (`FixMetabolism::grow_agent` → `limit_species`), where
`demanded = Δbiomass · yield_carbon` and `Δbiomass = m·μ·dt` with
`μ = μ_max · C/(km_carbon + C) · (other Monod terms ≤ 1)`:

```
demanded(C) = m·μ_max·Y·dt · C/(km_carbon + C)
```

Both are proportional to C at low concentration, so C cancels:

```
allowed/demanded = 4π·D_eff·r·(km_carbon + C) / (m·μ_max·Y)
```

This is *increasing* in C, so its floor is at C → 0, and that floor is a pure
parameter ratio with no concentration in it:

```
floor = 4π·D_eff·r·km_carbon / (m·μ_max·Y)
```

With the run's own values — `D_eff` = 5.0e-10 m²/s (carbon, retardation 1.0), mean
`r` = 5.58e-7 m and `m` = 8.14e-16 kg read from the final agent dump, `μ_max` =
5.0e-4 /s, `Y` = `yield_carbon` = 0.5, `km_carbon` = 5.0e-3 mol/m³ (`agent.cpp`),
`dt` = 60 s:

```
floor ≈ 86
```

**The diffusive cap is at least 86x slack at every carbon concentration, including
zero.** Checked against the run: μ_realized/μ_max = 0.096 implies C ≈ 5.3e-4 mol/m³,
at which the formula gives a ratio of 95, and the measured per-agent demand of
1.16e-18 mol/step against a cap of 1.11e-16 gives 96. That is why the binding count
is exactly 0 rather than small.

The floor is not an accident of these settings — it is the Berg–Purcell result. A
micron-scale cell is a very efficient absorber; diffusive supply to a single isolated
cell exceeds what its own kinetics can process by ~2 orders of magnitude, so
*single-cell* transport limitation is genuinely not the brake, and the model is right
about that. Real diffusive limitation at high density is a *collective* effect —
neighbouring depletion shells overlapping — which a per-cell formula evaluated at the
local voxel concentration can only express if that concentration actually falls. It
does not fall here, because agents take 4–15% of epithelial delivery while the VBF
consumes 139–217% of it as a homogeneous field they cannot deplete locally.

So there is no ordering to fix. Making the cap bind would require moving one of the
parameters in the floor by ~86x — `km_carbon` down, `μ_max` up, or `D_eff` down —
each of which is a measured quantity, and none of which is justified by wanting the
cap to matter.

This is consistent with the estimate that motivated #297 (~200x headroom at
10⁶–10⁸ cells/mL) and it settles the open question that came with it: the headroom is
not merely large at healthy densities, it cannot be consumed by lowering the
concentration at all, because demand and supply both scale with it.

## Consequences

1. **#297's Sherwood limit is correct physics with no dynamical effect in this
   regime.** It should stay (it is the honest funding model, and it removes the
   unfunded-uptake overdraw as a matter of accounting), but it is not a carrying
   capacity and must not be presented as one. `none` and `sherwood` are
   behaviourally interchangeable for campaign purposes at these densities.

2. **The 99% voxel write-off is confirmed as a discretization artifact, not a carbon
   constraint** — the third independent confirmation. With uptake funded from
   diffusive delivery instead of voxel content, nothing changes, because voxel
   content was never the physical limit.

3. **The model has no density-dependent brake within reach of the guard.** Combined
   with the delivery bracket — division linear in delivery at 0.72·mult per agent
   per hour, losses flat at 0.11/agent/h and delivery-independent — net growth still
   crosses zero exactly once, with nothing bending it back down. There is no basin,
   so no delivery rate holds a population for 168 h and a finer bracket buys a
   stochastic knife edge rather than a steady state.

4. **Therefore the 7-day coexistence horizon is not reachable by parameter choice.**
   Two honest options, and they are different experiments rather than different
   settings:

   - **Matched short horizon.** Run the RPS campaign at 24–36 h, inside which every
     arm is pre-guard and comparable, and treat guard crossing time as a reported
     outcome rather than a censoring nuisance. This measures competition and
     invasion, not long-run coexistence.
   - **Add a density-dependent mechanism deliberately.** The candidates that are
     physically real at 1e8 cells/mL and absent or inert here: volume exclusion /
     contact inhibition of growth (mechanics currently pushes cells apart but does
     not slow them), mucin-layer carrying capacity via local biomass fraction, and
     spatially resolved competition with the VBF (which today consumes 139–217% of
     epithelial delivery as a homogeneous field the agents cannot deplete locally).
     This is a mechanism addition, and it should be specified and reviewed before
     it is fitted.

My recommendation is to do both, in that order: run the campaign at a matched
24–36 h horizon now, since it answers the RPS questions that are answerable, and
raise density limitation as its own design task rather than as a campaign knob.
