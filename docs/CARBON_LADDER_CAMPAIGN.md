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
| 1 | carbon reference amplitude + epithelial boundary | `carbon.boundary_conc` (amplitude currently has **no** key — see gaps) | 5e-3 mol/m^3 | 0.084 -> 0.583 over x0.25 to x4 |
| 2 | gradient decay length | `carbon_z_lambda` | 25e-6 m | 0.114 -> 0.427 over 12.5 to 100 µm |
| 3 | imposed vs emergent profile | `carbon_z_gradient` | true | removes stratification entirely when false |
| 4 | mucin liberation | `vbf_mucin_liberation` | 5e-5 mol/m^3/s | 6.8% over x0 to x5; material only at >= ~100x |
| 5 | VBF carbon sink | `vbf_carbon_sink_vmax`, `vbf_carbon_sink_km` | 5.5e-5 mol/m^3/s, 1e-4 mol/m^3 | 2.3% over x0 to x5 |
| — | agent Monod half-saturation | none (`km_carbon` hard-coded 5e-3) | 5e-3 mol/m^3 | sets where on the ladder a given field sits |

## Gaps that block the ladder as designed

1. **No parser key for the carbon reference amplitude.** `carbon.boundary_conc`
   re-pins only the epithelial cell; the reference profile keeps
   `initial_conc = 5e-3`. Setting the boundary alone therefore produces a
   non-physical inversion — measured at `carbon.boundary_conc = 1e-3`:
   C(2 µm) = 1.00e-3 *below* C(26 µm) = 1.75e-3, i.e. carbon rising away from
   its source. The amplitude arms in the table above were only reachable by
   setting the field struct directly, which no config can do.
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

## Provenance

Agent-free measurements above were produced with an untracked standalone probe
against the serial Release library at `4644120` plus the delivery pricing fix;
they are supply-side scoping numbers, not campaign outputs, and no HDF5 or
campaign artifacts were committed.
