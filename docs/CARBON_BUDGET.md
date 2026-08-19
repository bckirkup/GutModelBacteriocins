# Carbon Budget — where a run's carbon actually comes from

Why this document exists: the 6-hour AWS calibration (Aug 2026) grew from 654 to
63,004 agents with **no deceleration**, while the grid carbon quantiles were
*identical* at steps 360, 720 and 1080. That is only possible if something
refills carbon as fast as the agents consume it, and for two hours the wrong
mechanism was blamed each time. The accounting described here exists so the
question is answered from an artifact instead of from argument.

Related:

- [`docs/AWS_CALIBRATION_6H.md`](AWS_CALIBRATION_6H.md) — the run that raised this
- [`docs/BRANCHING_FROM_CHECKPOINTS.md`](BRANCHING_FROM_CHECKPOINTS.md) — what to fork and why
- [`docs/MECHANISMS.md`](MECHANISMS.md) — VBF and metabolism mechanisms
- [`docs/WIRING_AUDIT.md`](WIRING_AUDIT.md) — species mass balance and coupling points

---

## 1. The four carbon terms

Carbon enters and leaves the domain through exactly four paths. Three of them
are configured mechanisms; the fourth is a boundary condition, which is why it
was invisible for so long.

| Term | Mechanism | Configured by |
|------|-----------|---------------|
| VBF liberation | mucin-derived carbon released by the background flora | `vbf_mucin_liberation`, `vbf_mucin_z_gradient*`, or the dynamic mucin path |
| VBF consumption | Monod sink representing background-flora competition | `vbf_carbon_sink_vmax`, `vbf_carbon_sink_km` |
| Agent uptake | `FixMetabolism::grow_agent()` writing negative reaction into `chem.reac` | strain `mu_max`, `K_carbon` |
| **Epithelial boundary** | net exchange across the `z=0` face, including the clamp and implicit z-solve boundary exchange | `carbon.boundary_conc` or `carbon_boundary_conc` |

The first three are bounded by parameters. The fourth defaults to a **Dirichlet
reservoir with no inventory behind it**: whatever the interior loses, the clamp
restores, forever. Carbon can instead select a finite-rate `robin` or `flux`
delivery variant at runtime; those variants meter epithelial supply through the
bottom face.

## 2. Measured scale of each term (calibration defaults)

Domain `2 mm × 2 mm × 0.1 mm`, `grid_dx = 2 µm`, so `A = 4e-6 m²`,
`V = 4e-10 m³`.

**VBF liberation.** With the default `mucin_liberation = 5e-5 mol/m³/s` and
`lambda = 25 µm`, depth-averaged over `h = 100 µm`:

```
r̄ = (1/h) ∫₀ʰ 5e-5 · exp(-z/25e-6) dz
  = 5e-5 · (25/100) · (1 - e⁻⁴)
  = 1.227e-5 mol/m³/s
```

As an area flux, and in units a physiologist can check:

```
J_VBF = r̄ · h = 1.227e-9 mol/m²/s
      = 10.6 nmol/cm²/day
      ≈ 1.9 µg glucose-equivalent/cm²/day
```

For scale, human ex vivo colonic mucus grows at ~240 µm/h at 1.3–1.9 wt%
solids, i.e. `~7.5–10.9 mg solids/cm²/day`. The model therefore treats roughly
**0.02%** of secreted mucin as immediately available monosaccharide. Mucin
polymer secretion and free glycan liberation are different quantities, so this
is a scale check and not an equality — but it means the liberation parameter is
conservative rather than generous.

**VBF consumption.** At the observed mean `c ≈ 1.176e-3 mol/m³`:

```
sink = 5.5e-5 · c/(1e-4 + c) ≈ 5.07e-5 mol/m³/s
```

which is **4× larger than the liberation term**. The background flora is a net
carbon *sink* in this configuration, not a subsidy.

**Epithelial boundary.** Fickian scale with `D = 5e-10 m²/s`, `C₀ = 5e-3
mol/m³`, `λ = 25 µm`:

```
J_boundary ≈ D·C₀/λ ≈ 1e-7 mol/m²/s
           ≈ 864 nmol/cm²/day
           ≈ 155 µg glucose-equivalent/cm²/day
```

≈ **81× the entire mucin source**. This is an upper-scale estimate, not the
realized flux — which is precisely why the accounting below records the realized
value.

## 3. Why growth never decelerated

Two consequences follow from the boundary being an unmetered reservoir, and both
were mistaken for other things during the calibration:

1. **The carbon field cannot deplete**, so identical quantiles across a 60×
   population increase are *expected*, not evidence that uptake fails to couple
   back to the field. Uptake does couple: `grow_agent()` writes into
   `chem.reac`, the GPU kernel does the analogous `atomicAdd`, and
   `chem.conc += chem.reac · dt` integrates it.
2. **Realized μ is pinned.** The default config has
   `carbon.boundary_conc == K_carbon == 5e-3 mol/m³`, so the Monod factor is
   *exactly* 0.5 at the richest point in the domain and never improves or
   degrades there. Mean realized μ stayed flat at 0.31 of `mu_max` across the
   whole run. Whether that equality is physiology or coincidence is a modelling
   question, but it has total control over whether anything can grow.

There is no carrying-capacity mechanism anywhere in the live path.
`VBFConfig::carrying_cap` defaults to `1e12` and `VBF::local_capacity()` returns
it, but **no call site applies either** to growth, division, death, placement or
mechanics. Growth is bounded only by nutrients, washout, and soft-sphere
repulsion. Volume exclusion is real but never binding at these densities: at
step 1080 the peak local occupied fraction was 0.28% within 10 µm, with zero
overlapping nearest-neighbour pairs and a minimum surface gap of 6.75 nm.

## 4. Reading the accounting from an artifact

Each summary write carries:

```
summary/step_NNNNNN/nutrient_flux/
  species_names                  # (nspec, 48) char — index → species, self-describing
  boundary_interval              # net mol supplied by the epithelial boundary this interval
  boundary_cumulative            #   … since t=0, survives resumes
  boundary_area_flux_interval    # mol/m²/s — compare directly against literature
  vbf_source_interval            # mol liberated by VBF
  vbf_source_cumulative
  vbf_sink_interval              # mol consumed by VBF (positive = removed)
  vbf_sink_cumulative
  agent_uptake_interval          # mol taken up by agents (positive = removed)
  agent_uptake_cumulative
  reaction_clip_interval         # mass discarded by reaction positivity clipping
  reaction_clip_cumulative
  interval_start_step, interval_end_step
  interval_start_time, interval_end_time
```

All per-species arrays are indexed by chemical species, so `species_names` is
the only lookup needed.

**Interval values are interval values.** They cover
`[interval_start_step, interval_end_step]`, not the whole run — the same
convention as the event counters, adopted after interval kill counts were
misread as a cumulative series. Use the `_cumulative` fields for totals; they
are restored from closed restarts and so are continuous across Spot resumes.

The question this is for:

```
boundary_cumulative[carbon]  vs  vbf_source_cumulative[carbon]
```

If the first dominates, the population is being fed by a boundary condition
rather than by anything with a budget, and growth rates measured from that run
describe the reservoir, not the ecology.

**Sign convention:** `boundary_*` and `vbf_source_*` are signed amounts *added*
to the domain (a negative boundary value means the clamp removed carbon, which
happens when the interior is richer than the boundary). `vbf_sink_*` and
`agent_uptake_*` are positive amounts *removed*.

## 5. How the numbers are obtained (and why it matters that they are cheap)

- **Boundary flux is taken where the boundary is applied.** For the default
  Dirichlet mode, the accounting adds two `O(nx·ny)` contributions per
  diffusing species to the existing boundary channel:
  1. the clamp-discard term
     `(boundary_conc - conc[idx]) · cell_volume` from
     `set_epithelial_boundary()`;
  2. the net z-solve face exchange
     `alpha · (diffusion_boundary - C_first_after_solve) · cell_volume`,
     summed over the x-y face where `diffuse_bounded_z()` applies it.
  This avoids differencing whole-grid inventories. In the
  gradient-preserving path, `diffusion_boundary = 0` for the departure
  field; the resulting exchange is still part of the epithelial boundary
  channel.
- **Finite-rate delivery is post-solve accounting.** Robin solves all `nz`
  cells with `beta = k·dt/dx_z` and records
  `beta·(C_epi - c0_after)·cell_volume` per bottom cell. Flux records
  `J·area·dt` per bottom face. These realized amounts are added directly to
  `boundary_cumulative`; no boundary clamp is applied in either mode.
- **The prescribed z-gradient is deliberately excluded.** Its
  subtract-diffuse-re-add is linear superposition and adds no net mass; the
  departure-field exchange at the z=0 face is included in the z-solve
  face-exchange term above.
- **One expression per mechanism.** VBF totals accumulate inside
  `apply_carbon_source` / `apply_carbon_sink` / `apply_iron_sink` /
  `apply_oxygen_sink` as they compute the reaction they apply, so the audit
  cannot drift from the physics it claims to describe.
- **Uptake is counted where it is produced** — in `grow_agent()` and the GPU
  metabolism kernel — so it is `O(agents)`, not a grid scan.

**MPI treatment is asymmetric, and inverting it is an easy mistake:**

| term | reduction | why |
|------|-----------|-----|
| VBF source/sink | **none** | every rank holds the full grid and computes an identical VBF field; summing would multiply by rank count |
| agent uptake | `Allreduce` sum | ranks hold only their own agents |

The instrumentation is verified to change no trajectory: identical
`FINGERPRINT` and `FINGERPRINT_STOCHASTIC` against `main`. It is safe to leave
on for a whole campaign rather than switching it on for audits.

## 6. Measured boundary flux and resolved carbon overdraw

- The scaled calibration measured a late-run **clamp-discard-only** flux of
  `1.17641933397e-11 mol/m²/s` (`0.101642630 nmol/cm²/day`). This
  clamp-discard-only number is approximately 104.3× below the
  `1.227e-9 mol/m²/s` mucin-source scale and approximately 8500× below the
  `1e-7 mol/m²/s` Fickian upper estimate; those ratios do not describe the
  complete z-solve boundary exchange.
- In the **pre-fix** diagnosis, with the z-gradient disabled, the boundary
  supplied approximately **77%** of reported agent-plus-VBF consumption. With
  the default gradient enabled, the boundary supplied approximately **28%**,
  while the unaccounted reaction-integration positivity clip supplied
  approximately **43%**.
- In the **post-fix** calibration with the default gradient enabled, the
  boundary supplied **29.1%**, VBF/mucin liberation supplied **70.7%**, and
  the residual reaction clip supplied **0.11%** of reported consumption. The
  pre-fix 43% was not carbon redirected to another pathway: it was the
  discarded explicit-sink demand, so removing it leaves a smaller honest
  denominator and raises liberation's share.
- The carbon VBF sink now uses backward-Euler implicit Monod integration and
  reports realized removal, so it cannot demand more carbon than the local
  field contains. The positivity-clip question for that VBF overdraw is
  therefore closed. The `reaction_clip_*` channel remains in every artifact so
  any residual clip is visible.
- Agent-side uptake can still overdraw a cell because realized uptake is not
  yet fed back into growth. That remains an open follow-up and is not fixed
  here.
- **Finite-rate epithelial delivery is implemented.** Set
  `carbon.epithelial_boundary` (or `carbon_epithelial_boundary`) to `robin` or
  `flux`. Robin uses `carbon.epithelial_transfer_coeff` in m/s and
  `carbon.boundary_conc` as `C_epi`; flux uses
  `carbon.epithelial_flux` in mol/m²/s. The default remains bit-compatible
  Dirichlet. These modes cannot be combined with `carbon_z_gradient`, because
  the gradient reference profile is pinned at the epithelial boundary.
- **Whether `boundary_conc == K_carbon` is intentional.** It fixes realized μ at
  half-max near the epithelium by construction.
- **Whether the VBF Monod sink has independent biological calibration** or was
  tuned to hold a target concentration.
- **The available fraction of secreted mucin.** The 0.02% above is the ratio the
  model implies, not a measured quantity.
