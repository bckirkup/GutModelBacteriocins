# Carbon capacity ladder — what actually sets the carbon gradient

Design document for the re-run of the #314 carbon ladder. Written **before**
running, per the work order: identify and measure the variables that set the
carbon field, so the ladder sweeps levers instead of inert knobs.

Status: design + agent-free sensitivity measurements. No production campaign
has been run from this document yet.

## Why the ladder cannot be re-run as originally specified

Two findings from the delivery/ROS campaign (`docs/DELIVERY_ROS_CAMPAIGN.md`)
change what the ladder is allowed to conclude:

1. **Uptake is not the carbon lever at these densities.** Carbon was funded at
   99.5% of demand with zero maintenance shortfall. Carbon-limited growth is set
   by the field concentration relative to `km_carbon = 5e-3 mol/m^3`
   (hard-coded in `src/core/agent.cpp`, no parser key), not by the delivery
   path.
2. **The 99.5% figure was measured in a demand-limited regime**, so the
   pre-fix delivery double-count did not bind on it. Delivery-limited carbon
   funding numbers measured before that fix are overstated and are not reused
   here.

So the ladder must sweep the terms that set the *concentration*.

## Measured sensitivity of the carbon field (agent-free)

Method, so the numbers are not quoted out of context: agent-free run
(`initial_strains` empty, so no uptake feedback and no crowding), domain
20 x 20 x 100 µm at `grid_dx = 4 µm` (5 x 5 x 25 cells), `bio_dt = 60 s`,
3600 s, serial Release + HDF5, `g++-13`. Column samples at x = y = 0,
z = 2 / 26 / 98 µm (cell centres of `iz` = 0 / 6 / 24). `monod` columns are
`C / (km_carbon + C)` with `km_carbon = 5e-3`. Agent-free means these are
*supply-side* field measurements only: they bound what the field can offer, not
what a crowded population would experience.

Imposed profile on (shipped default `carbon_z_gradient = true`,
`carbon_z_lambda = 25 µm`, amplitude and boundary `5e-3`):

| Arm | C(2 µm) | C(26 µm) | C(98 µm) | monod(26 µm) |
|---|---|---|---|---|
| shipped default | 5.000e-3 | 1.750e-3 | 9.04e-5 | 0.2592 |
| mucin liberation x0 | 5.000e-3 | 1.718e-3 | 4.49e-5 | 0.2557 |
| mucin liberation x0.2 | 5.000e-3 | 1.724e-3 | 5.40e-5 | 0.2564 |
| mucin liberation x5 | 5.000e-3 | 1.877e-3 | 2.73e-4 | 0.2730 |
| VBF carbon sink x0 | 5.000e-3 | 1.802e-3 | 1.68e-4 | 0.2649 |
| VBF carbon sink x0.2 | 5.000e-3 | 1.776e-3 | 1.21e-4 | 0.2621 |
| VBF carbon sink x5 | 5.000e-3 | 1.747e-3 | 8.67e-5 | 0.2589 |
| amplitude+boundary x0.25 | 1.250e-3 | 4.61e-4 | 5.63e-5 | 0.0844 |
| amplitude+boundary x0.5 | 2.500e-3 | 8.90e-4 | 6.72e-5 | 0.1511 |
| amplitude+boundary x2 | 1.000e-2 | 3.482e-3 | 1.49e-4 | 0.4105 |
| amplitude+boundary x4 | 2.000e-2 | 6.977e-3 | 3.00e-4 | 0.5825 |
| lambda 12.5 µm | 5.000e-3 | 6.40e-4 | 3.01e-5 | 0.1135 |
| lambda 50 µm | 5.000e-3 | 2.901e-3 | 6.18e-4 | 0.3672 |
| lambda 100 µm | 5.000e-3 | 3.730e-3 | 1.667e-3 | 0.4272 |

Emergent profile (`carbon_z_gradient = false`; the field is then set by
liberation, the VBF sink, and the epithelial boundary alone):

| Arm | C(2 µm) | C(26 µm) | C(98 µm) | C(26)/C(2) |
|---|---|---|---|---|
| default | 5.000e-3 | 4.821e-3 | 4.556e-3 | 0.964 |
| mucin liberation x0 | 5.000e-3 | 4.787e-3 | 4.505e-3 | 0.957 |
| mucin liberation x5 | 5.000e-3 | 4.958e-3 | 4.761e-3 | 0.992 |
| VBF carbon sink x0 | 5.000e-3 | 5.034e-3 | 5.052e-3 | 1.007 |
| VBF carbon sink x5 | 5.000e-3 | 4.726e-3 | 4.345e-3 | 0.945 |

**The carbon gradient in the shipped model is a boundary condition, not an
emergent source/sink balance.** Over a 25-fold span in mucin liberation
(0 to 5x) the Monod factor at 26 µm moves 0.2557 -> 0.2730 (6.8%), and over an
unbounded span in the VBF sink (0 to 5x) it moves 0.2649 -> 0.2589 (2.3%).
The same two terms with the imposed profile disabled produce no stratification
at all (C(26)/C(2) between 0.945 and 1.007 across the same spans).

Why, in timescales at z = 26 µm: diffusion `z^2/D_eff` with
`D_eff = 5e-10 m^2/s` is 1.35 s, while sink depletion `C/vmax` is
1.75e-3 / 5.5e-5 = 32 s and liberation refill `C/rate` is 35 s. Transport is
about 24x faster than either reaction term, the epithelial cell is
Dirichlet-pinned, and the reference profile `G(z)` is restored every step — so
the imposed profile behaves as an inexhaustible reservoir that the shipped
source/sink magnitudes cannot bend. Liberation only reshapes the field at
roughly 100x default, where it becomes comparable to the boundary-supplied
flux (measured earlier: the column flattens to 4.998e-3 at 26 µm and 4.822e-3
at 98 µm).

## Gradient-setting variables, ranked

| Rank | Variable | Key | Default | Effect on monod(26 µm) |
|---|---|---|---|---|
| 1 | carbon reference amplitude + epithelial boundary | `carbon.boundary_conc` (or explicit `carbon.z_amplitude`) | 5e-3 mol/m^3 | 0.084 -> 0.583 over x0.25 to x4 |
| 2 | gradient decay length | `carbon_z_lambda` | 25e-6 m | 0.114 -> 0.427 over 12.5 to 100 µm |
| 3 | imposed vs emergent profile | `carbon_z_gradient` | true | removes stratification entirely when false |
| 4 | mucin liberation | `vbf_mucin_liberation` | 5e-5 mol/m^3/s | 6.8% over x0 to x5; material only at >= ~100x |
| 5 | VBF carbon sink | `vbf_carbon_sink_vmax`, `vbf_carbon_sink_km` | 5.5e-5 mol/m^3/s, 1e-4 mol/m^3 | 2.3% over x0 to x5 |
| — | agent Monod half-saturation | none (`km_carbon` hard-coded 5e-3) | 5e-3 mol/m^3 | sets where on the ladder a given field sits |

## Gaps that block the ladder as designed

1. **Carbon reference amplitude (closed).** `carbon.boundary_conc` now sets both
   the epithelial boundary and the z-gradient reference amplitude by default.
   Capacity arms are reachable from configuration, and
   `carbon.z_amplitude` / `carbon_z_amplitude` provides an explicit override
   when a deliberately decoupled profile is needed.
2. **`km_carbon` is not configurable.** It is assigned in
   `src/core/agent.cpp`; a ladder in `C/km` space can only move `C`.
3. **Dynamic mucin liberation is dimensionally wrong and unusable.**
   `dynamic_mucin_liberation()` returns
   `k_liberation [1/s] * vbf.density [cells/m^3] * Monod(M)`, i.e. 1e-4 * 1e11
   = 1e7 mol/m^3/s at defaults — eleven orders above the static liberation
   term. Measured with `mucin.enabled = true` and no other change: carbon
   reaches 6.3e6 mol/m^3 at 6 µm and 3.4e7 mol/m^3 at 26 µm after one 60 s
   step (physiological scale is 5e-3, and even pure crystalline glucose is
   ~5.5e3), then washes out over roughly ten steps toward baseline. Mucin is
   emptied outside the epithelial layer in one step while carbon receives the
   full unlimited amount, so the two species are not stoichiometrically
   coupled. `mucin.secretion_rate` is additionally inert because liberation is
   Monod-saturated (`M = 1e-2` against `Km_degradation = 1e-3`). The path is
   off by default (`mucin.enabled = false`), so no shipped result is affected —
   but it must be fixed or refused before any mucin-driven carbon arm is run.

## Ladder to run, once the gaps are closed

All arms: serial Release + HDF5, delivery-limited uptake
(`metabolism.uptake_limit = delivery`, after the delivery pricing fix),
`oxygen.k_ROS = 0.0` so carbon is the only limiter under test, population-scale
(order 1e2 agents) rather than per-agent, with the asserted-contrast discipline
that caught #334's inertness.

- **Arm set A — capacity ladder.** Amplitude + boundary jointly at x0.25, x0.5,
  x1, x2, x4 of 5e-3 (monod at the epithelium 0.20 / 0.33 / 0.50 / 0.67 / 0.80).
  Asks: where does carrying capacity become carbon-limited rather than
  space-limited, and does realized growth track `C/(km+C)` at the agent's own
  z?
- **Arm set B — gradient shape.** `carbon_z_lambda` at 12.5 / 25 / 50 / 100 µm
  at default amplitude. Asks: how far off the epithelium can a colony be
  supported, i.e. does the vertical extent of the population track lambda?
- **Arm set C — emergent supply.** `carbon_z_gradient = false` with mucin
  liberation at x0 / x1 / x5 / x100 crossed with the VBF sink at x0 / x1 / x5.
  Asks: can any physiologically defensible liberation/sink pair generate the
  imposed profile's stratification on its own? The agent-free measurements say
  no below ~100x liberation; the arm exists to state that as a result rather
  than an assumption.

Each arm must report, and label: whether it early-terminated (dysbiosis guard
or horizon), whether agents were crowded (contact-limited rather than
carbon-limited), and the grid spacing used, since delivery support and shell
sampling are grid-dependent. The `sqrt(3)` shell artifact is not compensated.

## Measured ladder (28 arms, serial, local)

Run configuration, identical across arms unless named: 4 µm grid,
200 × 200 × 100 µm domain, 100 founders placed in a 0–20 µm z slab and carrying
no plasmids, `bio_dt = 60 s`, 12 h requested horizon,
`metabolism.uptake_limit = delivery` over the 10 µm support,
`oxygen.k_ROS = 0.0`, oxygen species disabled so carbon is the only limiter.
`dysbiosis_threshold` was raised to 1e10: at the shipped 1e8 cells/mL the guard
fires at roughly 400 agents in this volume, i.e. before any arm reaches a
plateau. 27 of 28 arms reached `horizon_reached`; the one exception is named
below. `km_carbon` is still hard-coded, so this is a ladder in `C`, not `C/km`.

### Two silent defects voided the first pass of this ladder

Both were invisible in the numbers, and both are the reason the results below
are the ones to cite:

- The first 20-arm pass was produced by a binary predating the delivery pricing
  fix (#344), so every gradient-on arm priced delivery-limited carbon at
  roughly twice the stored field. It reported 1377 agents at default carbon;
  the corrected binary reports 122. Any carbon number from that pass is void.
- The seeded λ replicates were produced by a binary predating #347, so all
  "replicates" ran at seed 0 and reported a within-λ spread of exactly zero —
  which is what prompted the audit. Exact reproducibility across nominally
  distinct seeds is a defect signature, not a clean result.

The analyzer now refuses to interpret an arm set whose `/run_provenance/git_sha`
is not identical and non-dirty across every arm. All numbers below come from
`a115447`.

### A — capacity: the shipped carbon default sits at the threshold of growth

| amplitude | C_epi | C(26 µm) | Monod at agents | peak N | divisions | boundary export | µ/µmax | rationing |
|---|---|---|---|---|---|---|---|---|
| x0.25 | 1.25e-3 | 4.61e-4 | 0.075 | 100 | 10 | 105 | 0.029 | 1 |
| x0.5 | 2.5e-3 | 8.90e-4 | 0.153 | 100 | 47 | 122 | 0.083 | 1 |
| x1 (default) | 5e-3 | 1.75e-3 | 0.146 | 122 | 194 | 190 | 0.082 | 1 |
| x2 | 1e-2 | 3.46e-3 | 0.318 | 1063 | 1547 | 584 | 0.197 | 1 |
| x4 | 2e-2 | 6.78e-3 | 0.371 | 9046 | 13083 | 4137 | 0.189 | 0.753 |

Capacity is a threshold in amplitude, not a proportionality. Below the shipped
default no arm exceeds its own founder count in 12 h — divisions are real but
are matched by boundary export (x0.25: 10 divisions against 105 exports), so
the patch is not self-sustaining. At the default it barely is: 122 peak against
100 founders, 194 divisions against 190 exports. A 2x amplitude buys 8.7x
capacity and 4x buys 74x. The shipped default is therefore not a mid-range
choice; it sits at the break-even point between division and export, which is
also exactly `km_carbon` (Monod = 0.5 at the epithelium, 0.35 at 26 µm).

Only the top arm is delivery-rationed (0.753, funded fraction 0.910) and it is
also the most contact-limited (1.85 displacement clamps/agent against ≤1.2
elsewhere, 600 bacteriostatic cells), so its µ/µmax of 0.189 is not a clean
carbon read-out — the ladder's growth-fraction claim is therefore made over the
unrationed range (0.029 → 0.197) only. Adjacent unrationed steps are not
separated by more than the seed spread measured in B (x0.5 and x1 are tied at
0.083 within 1%), so no step-by-step monotone claim is made. Maintenance
shortfall and reaction clips were zero in every arm.

### B — gradient shape: λ does not resolve vertical extent

Three seeds per λ (the 25 µm group shares the A x1 arm).

| λ | peak N per seed | p90 agent z (µm) per seed |
|---|---|---|
| 12.5 µm | 100, 100, 100 | 93.7, 91.8, 36.8 |
| 25 µm | 122, 100, 100 | 82.4, 95.2, 72.6 |
| 50 µm | 191, 100, 100 | 80.4, 77.4, 69.6 |
| 100 µm | 472, 128, 343 | 90.5, 69.3, 82.2 |

The single-seed pass of this arm set produced a monotone p90(z) in λ. That
ordering did not survive replication: the worst within-λ seed spread in p90 z is
0.77 of the group mean against a between-λ range of 0.12, so **vertical extent
is not resolved by λ at this scale** and the earlier ordering was noise. With
most arms pinned at the founder floor, the surviving agents' z distribution is
set by transport and export, not by how far carbon reaches.

What does replicate is a floor effect: at λ = 12.5 µm no seed achieves net
growth, at λ = 100 µm every seed does. Between those, one seed of three does.
Since the within-λ spread in peak N reaches 1.09 of the group mean, only that
never/always contrast is asserted; the mean 100 → 314 trend is not separately
claimable.

`R_lam_12p5um_s20260828` terminated at `population_stop` (extinction) at
3.51e4 s, not at the horizon, and its mean µ/µmax of −0.013 is a shrinking
population. It must not be quoted as a horizon value.

### C — emergent supply: the imposed profile is a penalty, not a supply

With `carbon_z_gradient = false`, liberation x0/x1/x5/x100 crossed with the VBF
sink x0/x1/x5. Across every physiologically defensible pair (liberation ≤ x5)
the field stays essentially flat: C(26 µm)/C(epithelium) ≥ 0.934, against 0.350
with the profile imposed. So no defensible source/sink pair reproduces the
imposed stratification — as the agent-free scoping predicted, now stated as a
measured result.

The population consequence is the direction that matters: those flat-field arms
reach peak 628–1291 agents against 122 for the imposed profile at the same
epithelial concentration. Relative to a uniform field at the boundary value,
the prescribed gradient is a **depletion penalty everywhere off the
epithelium**, not a supply mechanism; switching it off raises capacity 5–10x.
Peak N varies 2.06x across the nominally physiological arms, which is within
the seed spread measured in B, so the ordering *within* arm set C is not
interpretable — only the 5–10x contrast against the imposed profile exceeds it.

Liberation at x100 inverts the profile (C(26)/C(epi) = 1.55–1.66) instead of
steepening it, because VBF liberation is spatially uniform and has no
epithelial localization. It is a diagnostic of that fact, not a physiological
arm.

### What this does and does not license

- Carbon amplitude is a real population-scale lever, and the shipped default is
  at the edge of patch self-sustainment. Any campaign that needs a growing
  patch at default carbon needs a larger domain, a shorter mucus turnover, or
  more founders — not a small amplitude nudge.
- The prescribed z-gradient is a boundary condition whose net population effect
  is negative relative to a uniform field. It should not be described as the
  mechanism that supplies carbon to the mucus layer.
- λ, at this domain and horizon, is not a usable knob for colony vertical
  extent; claims about depth of colonization need a design where growth is not
  floored, and need replicates.
- Everything here is single-seed except arm set B. Differences smaller than
  roughly 2x in A or C are not resolved.
- Grid spacing is 4 µm and the delivery support is 10 µm; shell sampling is
  grid-dependent and the `sqrt(3)` shell artifact is not compensated.

## Provenance

The measured ladder above was produced locally and serially at
`a115447` (post-#344, post-#346, plus the #347 seed-parsing fix), with all
28 arms verified to share that single non-dirty sha. Campaign configs, HDF5
outputs, logs, the runner and the analyzer live outside the repository and are
not committed.

Agent-free measurements earlier in this document were produced with an
untracked standalone probe against the serial Release library at `4644120` plus
the delivery pricing fix; they are supply-side scoping numbers, not campaign
outputs, and no HDF5 or campaign artifacts were committed.
