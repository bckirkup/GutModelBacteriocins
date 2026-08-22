# Epithelial delivery bracket for a guard-safe 168 h horizon

12 runs, main `f1a7aee` (post-#300/#301), image digest
`sha256:9868a520…fd531`, seed 1001, dysbiosis guard on at 1e8 cells/mL,
168 h requested horizon. Only two things varied: `carbon.epithelial_flux`
at 0.10/0.14/0.18/0.22/0.26/0.30 × J_dir (J_dir = 1.0756e-8 mol/m²/s,
measured), and `kd_b12_btuB` at 1e-6 (colicin silent) and 1e-4 (colicin
biting). Analysis: `bracket_analysis.py`, raw numbers in
`bracket_analysis.json`, run inventory in `bracket_report.json`.

Domain volume 9.6e-5 × 9.6e-5 × 3.0e-4 m³ gives **276.5 agents = 1e8
cells/mL**, so the guard is a ~277-agent ceiling. Each run's
`halt_density_cells_per_mL` agrees with that conversion to the agent, which is
the check that the guard and this analysis are measuring the same quantity.

## Result: the band is empty

| kd | flux | reached | guard | peak density | final n (type1/type2) | verdict |
|---|---|---:|---|---:|---|---|
| 1e-6 | 0.10x | 86.1 h | no | 3.11e7 | 0 / 0 | extinct |
| 1e-6 | 0.14x | 43.9 h | YES | 1.18e8 | 325 / 0 | bloom |
| 1e-6 | 0.18x | 32.9 h | YES | 1.45e8 | 212 / 190 | bloom |
| 1e-6 | 0.22x | 25.1 h | YES | 1.04e8 | 162 / 126 | bloom |
| 1e-6 | 0.26x | 19.2 h | YES | 1.24e8 | 247 / 95 | bloom |
| 1e-6 | 0.30x | 17.2 h | YES | 1.13e8 | 220 / 92 | bloom |
| 1e-4 | 0.10x | 49.3 h | no | 3.11e7 | 0 / 1 | extinct |
| 1e-4 | 0.14x | 44.2 h | YES | 1.42e8 | 392 / 0 | bloom |
| 1e-4 | 0.18x | 33.9 h | YES | 2.02e8 | 473 / 85 | bloom |
| 1e-4 | 0.22x | 36.5 h | YES | 2.17e8 | 585 / 14 | bloom |
| 1e-4 | 0.26x | 20.0 h | YES | 1.18e8 | 262 / 64 | bloom |
| 1e-4 | 0.30x | 18.0 h | YES | 1.22e8 | 326 / 0 | bloom |

**No arm survived the week.** Every arm at or above 0.14x tripped the guard
inside 44 h; both 0.10x arms went extinct instead (the last producer dies
around 60 h, then the susceptibles decline to nothing). The intended answer —
"delivery rate X holds a steady population for 168 h" — does not exist
anywhere in [0.10, 0.30] × J_dir, so the horizon question cannot be settled by
choosing a delivery rate.

## Why: nothing in the model limits growth at density

Per-agent demography, first 17 h (matched across all arms), rates per live
agent per hour:

| flux | divisions | losses (outflow+lysis+colicin) | net |
|---|---:|---:|---:|
| 0.10x | 0.060 | 0.113 | −0.054 |
| 0.14x | 0.093 | 0.113 | −0.020 |
| 0.18x | 0.122 | 0.109 | +0.013 |
| 0.22x | 0.152 | 0.121 | +0.031 |
| 0.26x | 0.191 | 0.104 | +0.088 |
| 0.30x | 0.196 | 0.101 | +0.096 |

Division rate is **linear in delivery** (`div/agent/h ≈ 0.72 × mult − 0.008`,
R² ≈ 0.98) while total loss is **flat at 0.11/agent/h and delivery-independent**
— it is dominated by advective outflow, which does not care how much carbon
arrives. Net growth is therefore a straight line crossing zero once, near
0.16–0.18x, with no term that bends it back down as the population rises.

The reason nothing bends it: **the agents are not the carbon sink.** Cumulative
carbon at the end of each run, as a fraction of epithelial supply:

| flux | VBF sink / supply | agent uptake / supply | funded / demanded |
|---|---:|---:|---:|
| 0.10x | 2.17 | 0.038 | 1.0000 |
| 0.18x | 1.65 | 0.133 | 1.0000 |
| 0.30x | 1.39 | 0.110 | 1.0000 |

Agents book 4–15% of the carbon the epithelium delivers; the VBF eats 139–217%
of it (the excess over 1.0 is mucin liberation, not a budget violation). So
doubling the population barely depresses the field the population feeds on, and
`funded/demanded = 1.0000` in every arm confirms the reason it *cannot*:
`uptake_limit` is still `none`, so per-agent uptake is never capped by what can
physically reach the cell. **There is no density-dependent brake in these runs
at all** — which is exactly the thing the earlier analysis said was missing, now
visible as a population-level consequence rather than a per-step ledger oddity.

Consequences worth being explicit about:

1. A finer bracket between 0.10x and 0.14x would find a **knife edge, not a
   basin**. Successive 12 h windows of the 0.10x arm give net rates of −0.061,
   +0.012, +0.033, −0.044, −0.080, −0.016, +0.028 per agent per hour: at the
   crossover the *sign* of population growth is set by stochastic noise. A flux
   tuned there would give a different answer per seed, and the "coexistence for
   168 h" outcome would be a coin flip rather than a property of the
   configuration.
2. Guard-crossing time is then not a biological readout of the arm. It is
   `log(277/80) / r` with `r` set by delivery, which is why it falls smoothly
   from 44 h to 17 h as flux rises.

## What this says about `kd_b12_btuB`

The two corrinoid regimes are **separable from the delivery calibration and not
separable from the biology**:

- guard-crossing time is nearly identical between them at every flux (44 vs 44,
  34 vs 33, 20 vs 19, 18 vs 17 h) — colicin does not measurably change how fast
  the total population reaches the ceiling;
- but composition differs sharply. At 0.18x the susceptibles end at 85 (kd 1e-4)
  vs 190 (kd 1e-6); at 0.22x, 14 vs 126. Colicin deaths run 0.004–0.006 per
  agent per hour at 1e-4 and exactly 0.000 at 1e-6.

So the choice does not have to be made to calibrate delivery, and it should not
be inherited silently for the campaign: it is the difference between an arm
where the susceptibles persist and one where they are removed.

## Recommendation

The horizon problem is not a delivery-rate problem, so bracketing harder is the
wrong next move. The candidate that would actually create a carrying capacity is
the one already built and unmeasured: **`uptake_limit=sherwood`** (#297), which
caps uptake at the diffusive delivery to each cell and therefore binds harder as
cells crowd. Proposed next run, 6 arms, one seed, 168 h:
`uptake_limit ∈ {none, sherwood}` × `flux ∈ {0.14, 0.18, 0.22}x` at
`kd_b12_btuB = 1e-6`. That simultaneously (a) tests whether Sherwood produces a
sub-guard steady state and at what density, and (b) is the sherwood-vs-none
divergence measurement that was queued and never run.

If Sherwood does not bind at these densities either — the earlier per-cell
estimate said supply beats demand by 127–262x, so it may not — then the model
has **no** density-dependent limitation in the reachable range, and the 7-day
coexistence outcome in the RPS spec is not achievable by parameter choice. The
honest options then are a shorter matched-time horizon (24–36 h, which every
blooming arm supports), or adding a defensible density-dependent mechanism and
saying so explicitly rather than relying on the guard to end runs.
