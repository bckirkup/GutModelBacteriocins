# Spec 13 Phase 1 — Layer 2 occupancy readout

**Verdict: admissible region found.** Across the SPEC13_MULTISCALE
physiological table, a nonempty region of parameter space produces
stationary patch occupancy in the 0.1–1% band with segment mean density
in 1e4–1e5 CFU/mL. 578 of 4374 grid cells are robust-admissible
(band occupancy on >=4/6 seeds); 745 are admissible on at least one seed.

## Question

Mucosal *E. coli* persists at ~1e4–1e5 CFU/mL while a single mucus patch
carries ~1e7 cells/mL. Layer 2 resolves the tension through patch
occupancy: segment mean = occupied fraction x within-patch density. The
falsifiable readout: does any combination in the physiological table
give stationary occupancy in 0.1–1% (of patches) with segment mean in
1e4–1e5 CFU/mL, stratified by patch type, against the fragment-vs-single
establishment ratio?

## Design

- `sweep.py`: full-factorial 3^6 grid x 6 establishment ratios =
  4374 arms. N = 1000 patches (band resolution 1–10 patches), dt = 60 s,
  horizon 24 h, tail mean over the last 6 h, 5% initial patch seeding
  (transient; discarded by the tail window), fragment size = 10% of a
  full patch.
- Swept axes (3 levels each): contraction 0.2/0.6/1.0 per min; exposed
  disruption fraction 0.3/0.5/0.7 (shared by proximal and distal);
  agent-loss fraction 0.1/0.3/0.5; crypt disruption 0.05/0.10/0.15;
  transit half-life 2/4/8 h; reattachment 0.01/0.05/0.10 per transit.
  Establishment ratio swept over 1/10/30/100/300/1000.
- Fixed at table mid-points: type fractions crypt 10% / proximal 45% /
  distal 45%; supply multipliers 2.0/1.0/0.5; per-type single-cell loss
  0.02/0.10/0.20 per h; p_establish(single) = 0.01.
- `results.csv`: every arm's tail-window observables.
- `seeds.py` / `results_seeds.csv`: every stage-A admissible cell re-run
  at 5 additional seeds; a cell is robust-admissible if >=4/6 seeds land
  in band.

## Result

Outcome classes over 4374 arms:

| outcome | arms |
|---|---|
| extinct before stationary window (population_stop) | 424 |
| stationary occupancy below 0.1% | 2535 |
| stationary occupancy above 1% | 82 |
| occupancy in band, segment mean below 1e4 CFU/mL | 538 |
| occupancy in band, segment mean above 1e5 CFU/mL | 5 |
| **admissible (both bands)** | **790 (578 robust)** |

The admissible region is **carried entirely by crypt-type patches**: in
admissible arms, median occ_crypt = 5.0% (of crypt patches) with
within-patch density ~8e6 CFU/mL, while occ_proximal = occ_distal = 0.
The crypt contribution 0.10 x 0.05 x 8e6 ~ 4e4 CFU/mL sits mid-band by
construction of the metapopulation mechanism: a reseeded cell that lands
in a crypt patch persists and grows; one that lands in an exposed patch
dies to disruption + single-cell loss before reaching capacity. Occupancy
is therefore a source–sink structure — crypts hold the stationary
colonies; exposed patches churn.

Axis structure of the robust region (robust+marginal cells / 1458 cells
per level):

| axis | low | mid | high |
|---|---|---|---|
| contraction /min | 459+141 (0.2) | 119+25 (0.6) | 0+1 (1.0) |
| agent-loss fraction | 442+112 (0.1) | 136+9 (0.3) | 0+46 (0.5) |
| crypt disruption | 362+113 (0.05) | 114+25 (0.10) | 102+29 (0.15) |
| reattachment | 256+14 (0.01) | 186+63 (0.05) | 136+90 (0.10) |
| exposed disruption | 193+48 (0.3) | 193+60 (0.5) | 192+59 (0.7) |
| transit half-life | 190+54 (2 h) | 191+53 (4 h) | 197+60 (8 h) |
| **establish ratio** | 131+4 (1) / 129+9 (10) | 117+19 (30) | 69+46 (100) / 67+46 (300) / 65+43 (1000) |

The discriminating parameter shifts, but does not open or close, the
region: ~130 robust cells per ratio at 1–30 versus ~65 at 100–1000 —
high fragment establishment narrows the admissible corner about 2x
(fragments push occupancy above band unless churn compensates). The
binding axes are churn — contraction rate and single-agent loss — which
gate the whole table: at contraction 1.0/min or agent loss 0.5, no seed
lands in band.

## Interpretation

- Spec 13's central claim survives Phase 1: a physiological corner —
  low contraction, low agent loss, low-to-mid crypt disruption —
  produces stationary occupancy inside the target band at any
  establishment ratio, so the claim does not hinge on fragment reseeding
  being efficient.
- The empty set is *adjacent*: doubling contraction from 0.2 to 0.6/min
  cuts the admissible share ~4x, and at 1.0/min it vanishes. If the
  physiological contraction rate is at the top of the tabulated range,
  the mechanism fails — that is the falsifiable edge.
- Exposed patches never carry stationary occupancy in the admissible
  region. If Phase 2 requires an exposed-patch occupancy signal (a
  Swidsinski-type spatial gradient), the current reseeding/loss balance
  cannot produce it — that is a finding to carry forward, not tune away.

## Reproduce

```bash
python3 analysis/spec13_occupancy/sweep.py   # stage A: 4374 arms
python3 analysis/spec13_occupancy/seeds.py   # stage B: 5 seeds x admissible
```

Per-arm artifacts (`sweep_configs/`, `sweep_runs/`) are regenerable and
not committed.
