# Epithelial carbon delivery: the arms, what they show, and why the carbon is inaccessible

Written 2026-08-13. Six GPU runs on a Tesla T4, one seed, one image.
Companion machine-generated tables: `epi_delivery_stage1.md` (controls),
`epi_delivery_stage2.md` (arms). This document is the interpretation, the
mechanism, and what I think the decision is.

---

## 1. Provenance — what these runs actually are

Everything below comes from six runs of one image built from `main` at
`57cc8d9` (post #290/#292/#293/#295, i.e. the first code in which host
reactions reach the field, device metabolism runs with Fur, the device carbon
sink is implicit, and the epithelium has selectable delivery modes):

| | |
|---|---|
| image | `994254241749.dkr.ecr.us-east-1.amazonaws.com/gutibm:epi-57cc8d9` |
| digest | `sha256:2baa22af7f442dad911d90b2fcbb0a43a011c18295cc50cd7296da44de1285bd` |
| job definition | `gutibm-cuda-epi:1` (clone of `gutibm-cuda-campaign:4`, image swapped) |
| queue | `gutibm-gpu-practice`, g4dn.xlarge, Tesla T4, `REQUIRE_GPU=1` |
| seed | 1001 in all six |
| horizon | 1440 steps × 60 s = 24 h |

The protected campaign image and job definition were not touched.

| run | job id | ended |
|---|---|---|
| `epi_ctrl_grad_1001` | `2fd6e996-6f15-4e7d-aaba-a52fe305f844` | dysbiosis guard |
| `epi_ctrl_flat_1001` | `fc910417-2f23-44b3-8f8f-56a4dedfa7a3` | dysbiosis guard |
| `epi_flux100_1001` | `bc93d7ae-fc39-4be1-9e8f-23d2c55e8831` | dysbiosis guard |
| `epi_flux030_1001` | `f763af0d-8a79-47fc-a3b8-78ea4e23a661` | full 24 h |
| `epi_flux010_1001` | `e6c33a4c-a68e-4353-9e38-cd2afc5c36fb` | full 24 h |
| `epi_robin100_1001` | `b7673be9-3b2b-42d7-8bb6-0ecdc987823e` | full 24 h (dysbiosis guard at 950) |

All six report `GPU: ON (device 0)` and pass the `REQUIRE_GPU=1` activation
check, so none of them silently fell back to the CPU path.

---

## 2. The arms

The epithelium at `z = 0` used to be one thing: a Dirichlet cell pinned at
`C_epi = 5.0e-3 mol/m^3`. #295 made the delivery law selectable per species,
so the same epithelium can now be an infinite reservoir, a pump with a fixed
rate, or a membrane with a finite transfer coefficient:

| mode | boundary law | what it represents |
|---|---|---|
| `dirichlet` | `c(z=0) = C_epi` | infinite pantry: whatever is eaten is instantly replaced |
| `flux` | `J = J_fixed` | a pump: fixed mol/m^2/s regardless of what is there |
| `robin` | `J = k·(C_epi − c_0)` | a membrane: delivery falls as the near-wall concentration rises |

Two controls first, because the 12 runs behind the earlier analysis predate
the fixes above and could not be used for anything involving nutrients:

- **`epi_ctrl_grad_1001`** — the historical config: Dirichlet plus the
  exponential carbon z-gradient (`carbon_z_gradient=true`, λ = 25 µm).
- **`epi_ctrl_flat_1001`** — the same thing with the gradient off. This is the
  matched control, because the delivery modes *reject* the z-gradient by
  design (a fixed rate and an imposed profile are contradictory boundary
  statements, and the parser now says so rather than approximating).

Then the delivery arms, calibrated from the control rather than guessed. I
measured what the Dirichlet boundary actually delivers per unit area from the
control's own HDF5 flux ledger:

```
J_dir = boundary_interval / (FACE_AREA · Δt)  with FACE_AREA = (9.6e-5 m)^2
      = 1.0756e-08 mol/m^2/s   (peak, and essentially constant)
```

| arm | setting | rationale |
|---|---|---|
| `flux_1.0x` | `carbon.epithelial_flux = 1.0756e-08` | same mean delivery as Dirichlet, but *unpinned*: the near-wall concentration is now free to fall |
| `flux_0.3x` | `3.2268e-09` | a third of it |
| `flux_0.1x` | `1.0756e-09` | a tenth of it |
| `robin_1.0x` | `carbon.epithelial_transfer_coeff = 2.1512e-06 m/s` | `k = J_dir / C_epi`, i.e. delivers exactly `J_dir` when the wall cell is empty and less as it fills |

`flux_1.0x` is the discriminating arm: it holds the *rate* fixed at the
control's value and removes only the pinning. Robin then adds the one thing
a fixed pump cannot do — delivery that responds to local demand.

---

## 3. What the arms show

### 3.1 Delivery rate is a strong lever on the biology

Realized growth rate of the ColE1/ColB lineage, and what the population did
with it:

| arm | realized J (mol/m^2/s) | mean domain carbon | mean μ (1/s) | peak density (×1e8/mL) | h to reach 1× | 24 h outcome |
|---|---:|---:|---:|---:|---:|---|
| `dirichlet` | 1.075e-08 | 2.89e-03 | 1.678e-04 | 1.63 | 8.5 | guard-halted |
| `flux_1.0x` | 1.076e-08 | 2.89e-03 | 1.673e-04 | 1.32 | 8.5 | guard-halted |
| `robin_1.0x` | 5.411e-09 | 1.41e-03 | 1.053e-04 | 1.27 | 15.0 | guard-halted at 15.8 h |
| `flux_0.3x` | 3.227e-09 | 9.35e-04 | 7.677e-05 | 1.40 | 23.4 | reached horizon |
| `flux_0.1x` | 1.076e-09 | 4.82e-04 | 1.232e-05 | 0.29 | never | reached horizon, decayed to 16 agents |

Three things in that table matter.

**The lever works, and it is monotone.** Mean domain carbon is nearly linear
in delivery (1.00 : 0.49 : 0.32 : 0.17 against J of 1.00 : 0.50 : 0.30 :
0.10), and μ follows Monod on that concentration. Cutting epithelial delivery
tenfold cuts realized growth fourteenfold and converts a bloom that hits the
dysbiosis boundary in 8.5 h into a population that never reaches it at all
and decays to 16 agents. These differences are far outside the same-input GPU
spread (13% on final population) and outside seed-to-seed spread on crossing
time (CV 0.16), so they are effects, not realizations.

**`flux_1.0x` reproduces the control almost exactly** — μ 1.673e-4 vs
1.678e-4, same 8.5 h crossing, same guard halt. That is the clean result from
the discriminating arm: *unpinning the wall cell changes nothing while the
rate is unchanged.* The Dirichlet boundary was not smuggling in extra carbon
beyond its realized flux; my "infinite pantry" framing was right about the
mechanism (the pinned cell sat at `km_carbon` forever) but wrong about the
magnitude — at this population size the pantry was never being emptied fast
enough for pinning to matter.

**Robin lands between 0.3× and 1×, not at 1×.** Calibrated to deliver `J_dir`
at an empty wall, it realizes 5.4e-9, i.e. 50%, because the wall cell is
never empty: `c_0` equilibrates around `C_epi/2`. That is the finite-rate
membrane behaving as intended, and it is the only arm whose delivery is set
by the population rather than by me.

### 3.2 But it is a parameter, not a feedback

Within each arm, μ is essentially independent of how many cells are feeding:

| arm | ρ(density, μ) | μ range across density bands |
|---|---:|---|
| `dirichlet` | +0.187 | 1.665e-4 → 1.703e-4 (rises) |
| `flux_1.0x` | −0.050 | 1.666e-4 → 1.703e-4 (flat) |
| `flux_0.3x` | −0.283 | 7.59e-5 → 7.89e-5 (flat) |
| `robin_1.0x` | −0.441 | 1.059e-4 → 1.097e-4 (flat) |
| `flux_0.1x` | +0.789 | 9.54e-6 → 3.19e-5 (rises) |

Robin has the most negative correlation, which is the sign you would want, but
the per-band μ column is still flat to 4% across a 10× density range — the
correlation is picking up a slow trend, not crowding. And the division-vs-
density relationship that started this whole thread survives *unchanged* in
every arm that blooms (ρ = +0.60, +0.63, +0.61, +0.45 — the control's value),
which confirms the stage-1 correction: it is the founder cohort relaxing
toward the rate μ already implies (`ln2/μ` ⇒ 0.88 /agent/h), with density
collinear with time in a single growing run. No arm bends it down at high
density.

`bacteriostatic_live_agents` is 0 in every interval of every arm, including
the arm that decays to 16 agents. Nothing in these runs is ever classified as
viable-but-not-growing; the 0.1× population dies by lysis outpacing division
(0.125 vs 0.047 /agent/h), not by starvation arrest.

### 3.3 The clip is untouched

`reaction_clip / agent_uptake` per arm: **0.9967, 0.9968, 0.9952, 0.9881,
0.9968**. Every arm, including the starving one, books ~100× the carbon it
can be given and writes off ~99% of it. Delivery rate moves μ; it does not
move this at all. Which is the question you asked.

---

## 4. Why the carbon is inaccessible

Four candidate mechanisms. Three of them are not it, and the numbers say so
cleanly.

### 4.1 It is not polymer hindrance

The mucin retardation factor is an explicit per-species parameter
(`ChemicalSpec::retardation`, effective `D = D_free / retardation`). For
carbon it is **1.0** — no hindrance at all. The species that *are* retarded
are the four bacteriocins, at 10.0. So in the current parameterization the
polymer slows toxins by an order of magnitude and does not slow carbon at all.
Whatever the mucus gel does to sugar diffusion in reality, this model is not
currently representing it, and it is therefore not the explanation for the
clip. (It is a knob worth sweeping — see §5.)

### 4.2 It is not slow diffusion

With `D_carbon = 5.0e-10 m^2/s` and `Δt = 60 s`, the diffusion length per
biological step is

```
sqrt(2·D·Δt) = 245 µm
```

against a **2 µm** grid cell and a **300 µm** domain depth. Carbon crosses the
entire mucus layer within one step; the field is close to a
delivery-versus-consumption steady state at all times, which is exactly what
the near-constant `mean_carbon` trace in every arm shows.

The stronger version of the same point, at the scale of a single cell. The
diffusive supply to a sphere of radius `r` sitting in a medium at
concentration `c` is `4πDrc` (Smoluchowski). For `r ≈ 0.6 µm`:

| arm | demand per agent per step | Smoluchowski supply per step | headroom |
|---|---:|---:|---:|
| `dirichlet` | 2.49e-18 mol | 6.53e-16 mol | **262×** |
| `flux_1.0x` | 2.96e-18 | 6.54e-16 | 221× |
| `robin_1.0x` | 1.75e-18 | 3.19e-16 | 182× |
| `flux_0.3x` | 1.67e-18 | 2.11e-16 | 127× |
| `flux_0.1x` | 6.92e-19 | 1.09e-16 | 157× |

Physically, at these concentrations, an *E. coli*-sized cell can be supplied
by diffusion at 100–260× the rate it wants to eat. Even the arm that starves
to 16 agents has 157× headroom on transport. Diffusion is not the limit
anywhere in these runs.

### 4.3 It *is* the background flora, at the domain scale

`VBFConfig::carbon_sink_vmax = 5.5e-5 mol/m^3/s`, Monod with
`km = 1.0e-4 mol/m^3`, applied in **every cell of the domain**, unweighted by
depth (only the mucin *source* carries the z-gradient). Over one 60 s step
that sink can remove

```
vmax·Δt = 3.3e-3 mol/m^3
```

against a standing concentration of **2.89e-3 mol/m^3**. The continuum
anaerobe field can consume more than the entire standing carbon stock of every
cell, every step. It is implicitly integrated (#258) so it cannot overdraw,
but it takes essentially everything available.

Measured over the whole run (mol of carbon):

| arm | epithelial supply | mucin liberation | VBF sink | agent uptake booked | VBF sink / agent uptake |
|---|---:|---:|---:|---:|---:|
| `dirichlet` | 3.535e-12 | 4.11e-13 | 3.951e-12 | 2.727e-13 | 14.5 |
| `flux_1.0x` | 3.242e-12 | 3.77e-13 | 3.623e-12 | 2.208e-13 | 16.4 |
| `flux_0.3x` | 2.569e-12 | 9.95e-13 | 3.575e-12 | 2.090e-13 | 17.1 |
| `robin_1.0x` | 2.841e-12 | 6.56e-13 | 3.507e-12 | 2.064e-13 | 17.0 |
| `flux_0.1x` | 8.565e-13 | 9.95e-13 | 1.864e-12 | 2.007e-14 | 92.9 |

The budget closes on the flora, to three digits: epithelial supply plus mucin
liberation equals the VBF sink (3.535 + 0.411 = 3.946 vs 3.951e-12). **The
background flora eats the epithelium's entire carbon output.** The focal *E.
coli* population's booked demand is 7% of the throughput in the blooming arms
and 1% in the starved one, and the ~99% of it that is clipped is, in mass
terms, what is left after the continuum has been served.

This is the "Restaurant Hypothesis" sink doing exactly what its comment says
it does ("near-complete uptake at high [C]") — but its `vmax` was chosen to
represent "moderate carbon limitation for *E. coli*", and what it actually
produces is a flora that consumes 100% of epithelial delivery. That parameter,
not the boundary condition, is what sets how much carbon our agents ever see.

### 4.4 And it is the voxel rule, at the agent scale

Here is the number I think matters most. Per 60 s step, one agent books

```
uptake = Δbiomass · yield_carbon        (yield_carbon = 0.5 mol/kg)
```

and that is compared against the carbon in the single 2 µm grid cell it
occupies:

| arm | demand / agent / step | carbon in its voxel | ratio |
|---|---:|---:|---:|
| `dirichlet` | 2.49e-18 mol | 2.31e-20 mol | **108×** |
| `flux_1.0x` | 2.96e-18 | 2.31e-20 | 128× |
| `robin_1.0x` | 1.75e-18 | 1.13e-20 | 155× |
| `flux_0.3x` | 1.67e-18 | 7.48e-21 | 223× |
| `flux_0.1x` | 6.92e-19 | 3.85e-21 | 180× |

A single cell demands 100–200× the carbon its own voxel holds, per step, in
every arm. Since the timestep is operator-split — reactions are applied and
clipped at zero, *then* backward-Euler diffusion runs — an agent cannot draw
on the carbon that diffusion will bring into its voxel during the same step,
even though §4.2 says that carbon arrives 200× faster than it is needed. So
the realizable fraction is `voxel stock / demand ≈ 1/108`, and the measured
clip fraction is 0.99. **The 99% write-off is set by `Δx^3` and `Δt`, not by
biology.** It would be ~0.9 at 4 µm cells and ~0 at a 1 s step; it is
invariant to the boundary condition, which is precisely what the five arms
show.

Two consequences worth stating plainly:

- The clip is a **discretization artifact**, not a scarcity signal. It has
  been sitting in the ledger looking like starvation.
- Therefore μ, which is computed from the pre-step concentration and ignores
  realized uptake, is *not* obviously wrong: it is closer to the physically
  correct answer (transport-limited supply, 200× headroom) than a
  voxel-limited uptake would be.

The one place local depletion does appear to become real is the starved arm.
A mean-field Monod prediction from that arm's domain-mean carbon
(`5.5e-4 · 4.82e-4/(5.0e-3 + 4.82e-4)`, minus maintenance 1e-5) gives
≈ 3.9e-5 /s; the arm actually realizes 1.23e-5. Agents there are sitting in
cells materially poorer than the domain mean, i.e. genuine local drawdown at
the wall. In the blooming arms the same comparison agrees within ~20%. I would
not push that further without grid dumps — see §6.

---

## 5. What I think this means, and the decision

**On the science.** Epithelial delivery rate is now a real, calibrated,
monotone control on whether this system blooms into dysbiosis or stays inside
a defensible envelope, and that is a genuinely useful dial: `flux_0.3x` and
`flux_0.1x` are the first configurations that run a full 24 h *inside* the
population bound you asked for, without the guard cutting them off. But
delivery is a parameter I set, not a feedback the population generates. The
only arm where the population sets its own supply is Robin, and even there μ
is flat across a 10× density range.

**On the mechanism.** The reason nothing self-limits is not the boundary and
not diffusion. It is that the flora eats the entire carbon supply at the
domain scale, and that at the agent scale demand is compared against a
single-voxel stock that is 100–200× too small, one step before diffusion
refills it.

**So the throttle question changes shape.** The naive version — throttle
uptake to what the voxel holds and set μ from realized uptake — would impose
that 1/108 discretization factor as biology and collapse growth by two orders
of magnitude. That would be a grid artifact promoted to a finding, and I do
not recommend it. Three defensible options, and I would build them as variants
rather than pick one:

1. **Mass-transfer-limited uptake.** Cap each agent's uptake at the
   Smoluchowski/Sherwood diffusive flux to its own surface, `4πDr·c`, rather
   than at voxel content, and set μ from realized uptake. This is the standard
   IbM treatment and it is physically the right budget. Per §4.2 it would
   rarely bind at current concentrations — which is itself the honest answer
   to "should crowding limit growth here": at 1e8 cells/mL in a 300 µm mucus
   layer, transport is not the limit. It starts to bind as density rises
   toward the top of the dysbiosis window, so it gives us density dependence
   with a defensible mechanism instead of a grid coefficient.
2. **Voxel-limited throttle**, as the literal reading. Cheap to implement,
   and worth running *precisely so the divergence between 1 and 2 is
   measured* — that difference is the size of the discretization artifact,
   expressed in biology.
3. **Sub-step or finer grid**, i.e. remove the artifact rather than model
   around it: chemistry substeps within the 60 s biological step, or a finer z
   near the wall. Most expensive, and it makes 1 and 2 converge.

And one parameter I would put in front of you regardless of the throttle
decision: **`vbf_carbon_sink_vmax = 5.5e-5` is doing more to determine our
results than any boundary condition.** It currently consumes 100% of
epithelial delivery. It is documented as "moderate carbon limitation", it has
never been swept, and whether colicin producers can grow at the wall at all
depends on it. That is a one-line variant axis (`0`, `0.5×`, `1×`, `2×`) and I
think it belongs in any campaign matrix ahead of the delivery arms.

## 6. Caveats

- **One realization per arm.** Same-input GPU spread is ~13% on final
  population, so I am reading only the large trajectory and μ differences as
  effects; I would not defend the 0.3×-vs-0.1× peak-density ordering (1.40 vs
  0.29) as a mechanism without replicates, nor any of the ρ values to more
  than a sign.
- **Guard censoring.** Three arms were halted at the dysbiosis boundary and
  two ran the full 24 h, so all comparisons above are per density band or at
  matched simulated time, never final-state.
- **No grid dumps.** These runs wrote summaries, not carbon fields, so the
  local concentration *at agent positions* is inferred from μ rather than
  measured. A 200-step probe per arm with carbon grid output (~7 min each on
  the practice queue) would measure the wall drawdown directly and settle
  §4.4 without argument. Say the word and I will run it.
- The 0.1× arm's decay to 16 agents is one seed's realization of a population
  whose lysis rate exceeds its division rate; the *sign* is robust (μ is 14×
  down), the endpoint is not.
