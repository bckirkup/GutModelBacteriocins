# Density limitation: state of play

Working notes for the line of work that ran from "why does every run trip the
dysbiosis guard?" to the Pirt maintenance brake. Written as a handoff: what is
settled, what is measured, what is still open, and the exact next experiment.

Everything here is reproducible from the generators and analysis scripts in
this directory. Paths to run outputs are local scratch and are *not* committed;
the configs and the analysis code are.

---

## 1. The question

Every campaign-scale run either tripped the dysbiosis guard (sustained
1e8 cells/mL) within 4–44 h or went extinct. The 7-day RPS horizon was
therefore unreachable, and no delivery rate fixed it. The reason was structural
rather than parametric: **the model had no density-dependent brake.**

Division rate is linear in carbon delivery (0.72 x multiplier per agent per
hour) while losses — advective outflow, washout — are flat at ~0.11 per agent
per hour and delivery-independent. Net growth therefore crosses zero exactly
once, with nothing bending it back down, so the population either grows without
bound or decays to nothing depending on which side of the crossing it starts.

## 2. What was ruled out, with measurements

**Per-cell diffusive transport limitation cannot be the brake (#297's Sherwood
cap).** Demand is `m*mu_max*Y*dt*C/(km+C)` and the cap is `4*pi*D_eff*r*C*dt`,
so the local concentration `C` cancels and the ratio has a floor at `C -> 0`:

```
4*pi*D_eff*r*km_carbon/(m*mu_max*Y) ~= 86
```

The cap is >= 86x slack at *every* concentration including zero — the
Berg-Purcell result that a micron-scale cell is a ~100x-overpowered absorber.
Measured: 0 of ~2600 steps capped in every arm, and the `sherwood` arms were
identical to `none` down to the final cell count and guard time
(`sherwood/sherwood_vs_none.md`). Making it bind would require moving a
measured constant (`km_carbon`, `mu_max` or `D_eff`) by ~86x.

**Real diffusive limitation at high density is collective** — overlapping
depletion shells — and a per-cell formula reading its own voxel can only
express that if the voxel actually drains, which it cannot while the VBF
consumes 139–217% of delivery as a field agents cannot deplete locally.

**Delivery rate is not a free parameter that avoids the guard.** Bracketing
`carbon.epithelial_flux` at 0.10–0.30 x J_dir (J_dir = 1.0756e-8 mol/m^2/s,
measured): every arm >= 0.14x tripped the guard within 44 h, both 0.10x arms
went extinct, and successive 12 h windows of the 0.10x arm swing between
-0.080 and +0.033 per agent per hour — at the crossover the *sign* of growth is
noise, so "coexistence for 168 h" would be a per-seed coin flip
(`bracket/delivery_bracket.md`).

**Well-mixedness kills every per-voxel coupling.** Over one 60 s step the
diffusion lengths are 268 µm (acetate), 355 (O2), 200 (carbon) against a 300 µm
domain, while a few hundred agents occupy ~0.08% of 345,600 voxels. Any
mechanism specified per-voxel executes mean-field at this timestep.

## 3. What was implemented, and why those mechanisms

| PR | Mechanism | Why |
|---|---|---|
| #302 | Pirt carbon maintenance | A cell at mu=0 previously ate nothing: maintenance taxed the growth *rate* but consumed no carbon. O2 already had the correct Pirt form. Measured coefficient (~0.04 g glucose/gDW/h aerobic -> 2.1e-5 mol C/(s kg)), linear-in-N sink, so `N* = supply/(rate*biomass)` is a real carrying capacity rather than a fitted coupling. |
| #304 | O2 selects a metabolic mode | O2 previously only ever *helped*. Now respiratory capacity vs mu_crit overflow (Varma & Palsson) sets a fermentation fraction with 1 h inertia (DynamicME), inherited by daughters; fermentation costs its real carbon and secretes acetate; acid inhibition acts on the *undissociated* fraction via the pH lever from #301. |
| #305 | GPU negative-growth parity | #304 clamped `mu` at zero on the device only, so a starving cell had `mu_realized = 0` on the GPU and negative on the CPU: no biomass shrinkage, no washout-trap signal. |
| #306 | Maintenance budget | See below — #302 was inert as merged. |

Two spec values were deliberately inverted or dropped:

- **`anaerobic_yield_factor = 0.25` -> 4.1.** The literature ratio is
  biomass-per-substrate; this codebase's `yield_carbon` is the reciprocal
  (`demanded = d_biomass * yield`). As literally specified it would have made
  fermentation 4x *cheaper* and accelerated the bloom.
- **`agent_carbon_coupling` not implemented.** 1e-16 is 491x full-growth
  per-agent demand (2.04e-19 mol/s measured) and 1550x the entire delivery
  budget: it yields extinction, not capacity. With Pirt maintenance in place it
  is the same mathematical object with a fitted coefficient instead of a
  measured one.

## 4. The probe, and the bug it found (#306)

Six arms, seed 1001, 24 h, guard on, `kd_b12_btuB = 1e-6` (colicin silent),
mechanisms off / maintenance only / full, and full mode at 0.18/0.30/0.60/1.00x
delivery. See `probe/probe_findings.md` for the tables.

**#302 was inert in every real configuration.** The maintenance-only arm was
identical to the mechanisms-off arm hour by hour, because the draw was clamped
to the carbon *standing in the agent's 2 µm voxel*:

```
voxel stock      = 6.615e-4 mol/m^3 * 8.0e-18 m^3 = 5.29e-21 mol
per-agent demand = rate * f_maint * biomass * 60 s = 9.7e-19 .. 3.2e-18 mol
                 = 184x .. 614x the voxel stock
```

so it collected ~0.3% of its Pirt demand, while growth uptake in the same step
was not clamped at all. The instantaneous voxel stock is the wrong budget at
`bio_dt = 60 s` (carbon diffuses ~77 µm in 60 s against a 2 µm cell, refilling
the voxel thousands of times per step); diffusive delivery to one cell over the
step is ~2.7e-17 mol, 8–28x the demand. #306 routes maintenance through the
same `uptake_limit` model growth uses, and fixes the shortfall field, which was
counting agent-steps in a mol-valued ledger.

**The metabolic mode is not masked, contrary to my prediction.** I expected O2
supply through the mucus to support ~1 agent's respiration, pinning the
fermentation fraction at 1.0 and making the aerobic branch unobservable.
Measured: mean realized fermentation fraction 0.156 at 0.18x rising to 0.550 at
1.00x, individual agents spanning 0.018–0.959, with mean O2 4.586e-6 mol/m^3 —
above its Monod Km. Both branches are live, and density and growth rate rather
than O2 supply alone select the mode. It costs a real 47% more carbon per unit
biomass (124 vs 160 agents at matched flux and matched time).

## 5. Correction to an earlier prediction

I twice told Benjamin to expect the guard-safe delivery band to sit *above* the
0.10–0.30x range that all blew through the guard. **That is wrong, and the
arithmetic points the other way.** With the brake funded, capacity is
proportional to delivery:

```
supply at 0.18x   = 1.936e-9 mol/m^2/s * 9.216e-9 m^2 = 1.784e-17 mol/s
per-agent draw    = 2.1e-5 * f_maint * 8.08e-16 kg
                  = 1.70e-20 mol/s          (f_maint = 1, fully respiratory)
                  = 5.42e-20 mol/s          (f_maint = 3.19 at f_ferm 0.156)
N* = supply/draw  = 1050 agents (3.8e8 cells/mL)   fully respiratory
                  =  329 agents (1.2e8 cells/mL)   at the measured f_ferm
guard             =  276 agents (1e8 cells/mL at the 2.765e-12 m^3 domain)
```

So at the measured fermentation fraction the plateau at 0.18x sits *just above*
the guard, and a guard-safe plateau needs roughly **<= 0.15x**, with extinction
somewhere below ~0.06x. The band is narrow — but unlike the pre-brake case it
is a **stable fixed point rather than a knife edge**: per-capita supply falls as
N rises, so trajectories converge to N* from either side instead of diverging.
The brake is also self-sharpening, since f_ferm rises with density and
`f_maint = 1 + f_ferm*(anaerobic_maintenance_factor - 1)` rises with it.

`anaerobic_maintenance_factor` (default 15) scales N* inversely and is the one
number in this mechanism with no direct anchor — only the aerobic and anaerobic
Pirt endpoints. It should be a campaign axis, not a pinned constant.

## 6. Next experiment

Re-bracket delivery with the brake funded, before anything else:

- flux multipliers **{0.06, 0.09, 0.12, 0.15, 0.20, 0.30} x J_dir**, centred on
  the predicted <= 0.15x band rather than the old range;
- 168 h, guard **on** (it is a model-validity boundary, not a nuisance),
  seed 1001 to locate, then 3 seeds to confirm the winner;
- full mechanisms as in the probe's `full_*` arms;
- `kd_b12_btuB = 1e-6` for comparability, knowing colicin is silent there;
- grid output **enabled** in at least one arm — the probe had it off, so acetate
  *concentration*, which acid inhibition reads, was never measured;
- second axis if cheap: `anaerobic_maintenance_factor` in {5, 15}.

Pass/fail is explicit: a guard-safe arm must show a **plateau** (net per-capita
growth crossing zero with N settling, not merely a slow rise), with
`maintenance_limited_agents` at or near zero — if maintenance is still being
clamped, the brake is again running below strength and the numbers mean nothing.

## 7. Still open

- **Acid inhibition is unevaluated.** Acetate is produced and drained
  (net -8.5e-12 mol across the epithelial boundary at 0.18x), but the field is
  background-dominated (VBF net source 8e-4 mol/m^3/s vs agents' 1.5e-6 at
  N=325), and `Ki = 20 mol/m^3` is ~60x more inhibitory than Russell's ~50 mM
  undissociated supports (the spec's "1.1 mol/m^3 ~= 66 mM undissociated" is a
  unit error; 1.1 mol/m^3 is 1.1 mM, and matching Russell needs ~900 mol/m^3
  total). Re-anchor Ki to what the field can actually reach.
- **`yield_carbon = 0.5` is ~6.7x cheaper than E. coli** (literature implies
  ~3.2–3.33 mol glucose/kg wet). Left alone deliberately: changing it moves
  every prior run, and it is most of why agents booked 4–15% of epithelial
  carbon while the VBF ate 139–217%. It belongs in the campaign as a factor.
- **`kd_b12_btuB` cannot be inherited silently.** Across 1e-6 -> 1e-3
  (competitive factor 1001 -> 2) the susceptible population ends at 73/68/39/3
  cells from the same 20 founders, i.e. from "colicin is irrelevant" to
  "colicin eradicates the susceptibles". The two-ligand analogue argument
  (98% non-Cbl analogue at ~10 nM) lands at ~5e-6. It does **not** move guard
  timing (44/44, 20/19, 18/17 h), so it is separable from delivery calibration
  (`../rps_campaign/corrinoid_sweep.md`).
- **Resistance cost must stay mechanistic.** Receptor loss pays its own
  emergent cost through `Km_b12 = km_b12/expr` with the expression floor
  (~9% at 1 µM corrinoid, ~25% at 0.3 µM), and it moves with the same corrinoid
  knob that sets colicin competition. A fixed `mu_max` penalty on top charges
  the cell twice.
- **Explicit formate is deferred.** Formate's inhibition enters through the same
  undissociated-acid term, so the acetate species carries it as total
  fermentation acid load. It becomes necessary only for formate-specific biology
  (FHL, H2/CO2, cross-feeding) or its own pKa (3.75 vs acetate's 4.76).
- **The RPS 7-day coexistence outcome is still unproven** and depends entirely
  on the re-bracket finding a plateau.

## 8. Operational notes worth not rediscovering

- `gpu-device-tests.yml` triggers on `pull_request` only: **a PR merged before
  its checks land never gets T4 coverage**, and nothing re-runs it on `main`.
  This is how #304 merged with a broken CUDA translation unit.
- A completed AWS Batch job is not a completed simulation horizon. Always read
  `/run_provenance/halt_reason_code` (1 = dysbiosis threshold) and the final
  `time` before comparing arms; compare at **matched simulated time**.
- Local CUDA is unavailable in the dev environment, so `gpu_*` ctest targets
  skip locally. The hosted CUDA compile plus the physical-T4 gate are the only
  authoritative GPU evidence.
