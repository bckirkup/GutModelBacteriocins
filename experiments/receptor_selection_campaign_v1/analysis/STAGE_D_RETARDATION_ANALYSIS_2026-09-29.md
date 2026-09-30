# Stage D retardation-sensitivity analysis — 33 runs (2026-09-29)

Paired readout of two new arrays plus the existing refine anchors, all on
execution source `29fd563c`, image digest `fafc77ee`, job definition
`gutibm-cuda-campaign:14`:

- endpoint array `f0e1c478-0f49-49af-9a2b-b0b96524f248` (18 children:
  mucin-charge amplitude 0 and 60 x targets {0.25%, 1%, 5%} x seeds
  20261009/1011/1013)
- mid array `11eb5ef8-3fac-4f39-a0a8-2aaf640eb255` (6 children: amplitude
  15 x targets {1%, 5%} x same seeds)
- refine anchors `adfd5af1` reused: indices 0-2 (amp-15 zero-lysis
  carriers), 3-5 (amp-15 @0.25%), 18-20 (plasmid-free nulls)

The mucin-charge amplitude sets ColE1 retardation through
`retardation_from_pI` (pI 9.0): amplitude 0 -> r = 1.2 (literature
D_eff ~3.3e-11 m^2/s), amplitude 15 -> r ~ 13.5 (the deployed Stage D
measurement), amplitude 60 -> r ~ 50.2 (the pI-derived shipped value).
This probe asks whether the Stage D dose-response is a measurement of
the lysis prior or of the transport parameter.

## Integrity

- 33/33 outputs readable; `execution_source_sha` matches `29fd563c` in
  all; `chemistry_placement=device_delivery` in all.
- All 24 new runs succeeded first attempt and terminated on the
  **dysbiosis guard**; paired contrasts use the shared calendar window
  ending at `t_common` = 16800 s per seed.
- n = 3 seeds per cell (screening replicates — medians and ranges, no
  p-values). Same-seed plasmid-free nulls pair every cell; zero-lysis
  carriers give the released-nothing baseline.
- Realized lysis is amplitude-invariant within noise at every target
  (5% arm: 0.0451 / 0.0458 / 0.0444 across amp 0/15/60) — the amplitude
  lever moved toxin transport only, not the lysis rate, so the design
  isolates the transport axis.

## Per-arm medians

| arm | n | divisions | lysis | realized lysis / producer division | kills | kills / lysis | own-end slope (dec/h) |
|---|---|---|---|---|---|---|---|
| amp 0, target 0.0025 | 3 | 18,941 | 20 | 0.00258 | 275 | 13.8 | +0.009 |
| amp 0, target 0.0100 | 3 | 27,020 | 111 | 0.00806 | 1,190 | 11.4 | +0.061 |
| amp 0, target 0.0500 | 3 | 24,603 | 599 | 0.04507 | 4,888 | 8.0 | +0.191 |
| amp 15, target 0 (carrier) | 3 | 19,279 | 0 | 0 | 0 | - | -0.001 |
| amp 15, target 0.0025 | 3 | 26,267 | 32 | 0.00232 | 1,986 | 54.4 | +0.086 |
| amp 15, target 0.0100 | 3 | 23,220 | 117 | 0.00820 | 3,456 | 27.0 | +0.206 |
| amp 15, target 0.0500 | 3 | 28,703 | 1,122 | 0.04583 | 3,721 | 3.3 | +0.597 |
| amp 60, target 0.0025 | 3 | 25,683 | 35 | 0.00233 | 2,681 | 81.2 | +0.137 |
| amp 60, target 0.0100 | 3 | 35,076 | 249 | 0.00878 | 5,967 | 23.7 | +0.448 |
| amp 60, target 0.0500 | 3 | 27,363 | 1,090 | 0.04441 | 2,263 | 2.1 | +0.744 |

## Paired contrast vs same-seed plasmid-free null (shared `t_common` window)

Delta slope (log10 producer:sensitive ratio per hour, treatment minus null):

| target | amp 0 (r = 1.2) | amp 15 (r ~ 13.5) | amp 60 (r ~ 50.2) |
|---|---|---|---|
| 0.0025 | +0.025 [+0.023, +0.029] | +0.073 [+0.069, +0.082] | +0.100 [+0.084, +0.162] |
| 0.0100 | +0.051 [+0.048, +0.053] | +0.166 [+0.160, +0.196] | +0.248 [+0.235, +0.258] |
| 0.0500 | +0.133 [+0.121, +0.133] | +0.372 [+0.363, +0.402] | +0.517 [+0.487, +0.529] |
| 0 (carrier, amp 15) | | +0.014 [+0.014, +0.018] | |

Delta tail-median (decades) follows the same ordering; at 5% it is
+0.046 / +0.333 / +0.568 across amp 0/15/60.

## Readout

- **The gradient moves along the transport axis — strongly.** At every
  target the paired delta slope is monotone in amplitude, with a ~4-5x
  contrast between endpoints (5% arm: +0.133 vs +0.517 dec/h). The Stage
  D dose-response is therefore **not** a measurement of the lysis prior
  alone; it convolves the prior with the retardation value, exactly the
  confound `docs/SPEC13_LYSIS_SELECTION.md` warns about ("a lysis
  coefficient calibrated at current retardation absorbs the transport
  error"). The admissible prior must be reported conditional on
  retardation; the SPEC13 two-retardation-arm design ({50, 1.2}) is
  confirmed necessary, and the deployed amplitude-15 grid is not a
  transport-neutral midpoint.
- **Retardation concentrates the benefit of each lysis.** kills/lysis at
  the sparse end rises with amplitude (13.8 at r = 1.2, 54.4 at r ~ 13.5,
  81.2 at r = 50.2 for the 0.25% target): a strongly retarded toxin is
  released in a tight halo around the producer and repeatedly kills the
  sensitives competing closest to it. At r = 1.2 the same release rate
  spreads so thin that even a 5% lysis rate buys only +0.13 dec/h.
- **The gradient's shape, not only its level, shifts.** At amp 0 the
  response stays near-linear in target (0.025 -> 0.051 -> 0.133; ~2.6
  dec/h per percent at the top); at amp 60 the same target axis gives
  0.100 -> 0.248 -> 0.517, steeper through mid-range. The "knee at
  0.20-0.25%" located in the sub-grid analysis is a property of the
  amplitude-15 transport regime.
- **Direction for calibration:** at fixed realized lysis, more
  retardation means more selection benefit per lysis. A prior selected
  at r ~ 50.2 is not transferable to r = 1.2 — the two endpoints demand
  roughly 4x different lysis rates to express the same producer
  advantage.

## Classification

- **Measured** (29fd563c, seeds 20261009/1011/1013, n = 3 per cell):
  the three-amplitude paired deltas and per-arm rows above.
- **Inferred:** the lysis-prior selection is not transport-invariant;
  it must be parameterized jointly with (or conditional on) the
  retardation value.
- **Hypothesis:** which retardation is physically correct remains open —
  literature D_eff analogues sit near r = 1.2 while the shipped
  pI-derived value gives r ~ 50.2. The 12-arm selection must bracket
  both, as SPEC13 already specifies.

## Caveats

- n = 3 per cell, screening replicates; amp-60 @0.25% shows one high seed
  (+0.162 vs +0.084/+0.100) — wider range, same ordering.
- All arms halt on the dysbiosis guard; contrasts on the shared window
  measure relative advantage over the window, not a fixed-point
  displacement outcome. The patch scale remains a lysis-response
  observable, not the weeks-scale in-vivo endpoint (that needs Layer 2+).
- Carrier baseline (+0.014 dec/h) is retained from the amp-15 refine
  run; carriers emit no toxin so amplitude does not act on them.
