# pI-driven mucin retardation for bacteriocins — design and calibration

Authored by lead. Implementation must follow this spec exactly; the numeric
constants and the functional form are not to be re-derived.

## 1. What exists today

`BICluster.pI` is carried through bursts, `GreensFunctionParams`, HDF5 and MPI
transfer, and is read by **no physics**. Its only functional effect anywhere is
`classify_by_pI()` setting `bclass` when `FixMutation` invents a novel BI, and
`bclass` has no consumer outside I/O and agent transfer. The pI to retardation
coupling that the documented Lethal Core / Lethal Halo concept requires lived
only in `FixBacteriocin::retardation_for_pI`, which was never called and is
deleted in PR #300.

Actual toxin spread is set by the independent hardcoded per-plasmid
`retardation` field:

| Plasmid | pI | `retardation` | `diff_coeff` (m^2/s) | D_eff = D/R (m^2/s) |
|---|---:|---:|---:|---:|
| MccV (microcin V) | 5.0 | 1.2 | 1.0e-10 | 8.3e-11 |
| ColB | 5.4 | 1.5 | 4.0e-11 | 2.7e-11 |
| ColE2 | 6.5 | 3.0 | 3.5e-11 | 1.2e-11 |
| ColIa | 7.2 | 5.0 | 4.0e-11 | 8.0e-12 |
| ColE1 | 9.0 | 50.0 | 4.0e-11 | 8.0e-13 |
| ColM | 9.3 | 60.0 | 5.0e-11 | 8.3e-13 |

## 2. Calibration finding

The six hardcoded retardations are **not** independent numbers: they lie on a
charge-titration sigmoid in pI. Fitting

    R(pI) = R_min + A / (1 + 10^((pI_0 - pI) / w))

to the table above gives R_min = 1.2, A = 60, pI_0 = 8.35, w = 1.0, which
reproduces every library value within about 30%:

| Plasmid | pI | library R | fitted R |
|---|---:|---:|---:|
| MccV | 5.0 | 1.2 | 1.23 |
| ColB | 5.4 | 1.5 | 1.27 |
| ColE2 | 6.5 | 3.0 | 2.04 |
| ColIa | 7.2 | 5.0 | 5.40 |
| ColE1 | 9.0 | 50.0 | 50.2 |
| ColM | 9.3 | 60.0 | 55.1 |

So the library already encodes a pI coupling implicitly; making it explicit is
mostly a refactor of intent, not a change of regime.

Mechanistically the sigmoid should be centred on **net charge**, not on pI in
absolute terms: a protein is cationic when the local pH is below its pI, and
mucin is polyanionic across the whole physiological range (sialic acid and
sulfate pKa about 2.0 to 2.6, so mucin charge does not titrate here). The
parameterisation is therefore in `(pI - pH_local)`:

    R(pI) = R_min + A / (1 + 10^((dz_0 - (pI - pH)) / w))

with defaults `R_min = 1.2`, `A = 60.0`, `dz_0 = 1.35`, `w = 1.0`, and pH
defaulting to 7.0 (colonic mucus). At pH 7.0 this is numerically identical to
the pI_0 = 8.35 fit above, so the defaults reproduce section 1 within 30% while
the form now responds to pH.

## 3. The unresolved constant, and why it becomes a campaign axis

`A` is the maximum charge-driven excess retardation and it has no measured
anchor. The two available anchors disagree by a factor of about 40:

- The library values imply `A = 60`: ColE1 is retarded 50x, D_eff = 8.0e-13.
- The literature parameterisation supplied by the user gives ColE1 retardation
  about 1.2, i.e. essentially charge-inert transport, D_eff about 5.8e-11.

`A = 0` collapses the model to the second anchor (all bacteriocins at R_min,
transport set by molecular weight alone, no core/halo distinction). `A = 60`
is today's behaviour. This is exactly the quantity the RPS spec's retardation
axis should sweep, in place of three arbitrary per-plasmid numbers:

    A in {0, 15, 60}

Report D_eff and the resulting near-field/far-field toxin profile per arm.

## 4. Required implementation

1. Add `retardation_from_pI(Real pI, const MucinChargeConfig&)` as the single
   source of truth, in the same translation unit as `classify_by_pI()`
   (`src/genome/plasmid.{h,cpp}`) so pI-derived quantities stay together.
2. Config struct with the section 2 defaults, parser keys under
   `bacteriocin.mucin_charge.` for `r_min`, `amplitude`, `dz_half`, `width`,
   and `ph`. Follow the repo's new-config-key checklist: flat key plus
   `config_json.cpp`, fixture, `test_input_parser.cpp` assertion,
   `test_config_ingestion.cpp` probe, resolved-config serialisation.
3. `BICluster.retardation` stops being an independent physical parameter. It is
   populated from `retardation_from_pI(pI)` when the cluster is built, and the
   PR #300 per-plasmid `plasmid_overrides.<name>.retardation` continues to win
   over the derived value (that is the documented escape hatch).
4. Mutation: when `FixMutation` drifts `novel.pI`, recompute `retardation` from
   the new pI alongside the existing `bclass` reclassification. This is the
   point of the change — a mutated pI must move transport, not just a label.
5. Checkpoint/restart and MPI transfer already serialise `retardation`, so a
   restored cluster keeps its resolved value; do not recompute on load, or an
   overridden cluster would silently revert.

## 5. Required tests (behaviour, not goldens)

- Graded sensitivity: retardation is monotonically non-decreasing in pI, and
  strictly increasing across the transition (pI 6, 7, 8, 9, 10).
- Bounds: `R_min <= R <= R_min + A` for pI in [3, 12]; finite everywhere.
- `A = 0` makes retardation independent of pI (all clusters at `R_min`).
- Raising `ph` shifts the curve **down**: at fixed pI, higher pH gives *lower*
  retardation, because net positive charge scales with `(pI - pH)` and a protein
  in a more alkaline medium is less cationic. Equivalently, retardation is
  monotonically non-decreasing in `(pI - pH)`, which is the single invariant
  worth asserting.
- Transport consequence, via the Green's function: three pI values produce
  ordered near-field and far-field concentrations, same ordering as an explicit
  retardation sweep (near-field up, far-field down with increasing pI).
- Mutation coupling: a novel BI whose pI drifts upward past the transition
  gets a higher retardation than its parent.
- Override precedence: `plasmid_overrides.ColE1.retardation` wins over the
  derived value, and the resolved config records the override.
- Defaults reproduce the section 2 table within **35%** for all six library
  plasmids (labelled as a change-detector on the calibration, not a physics
  assertion). ColE2 is the one expected outlier at -32%: its library pI of 6.5
  is the pI of the toxin plus Im2 *complex*, not of a single protein, so its
  hardcoded retardation of 3.0 was never set from the same single-protein
  charge argument the curve encodes. Name that in the test as the reason the
  tolerance is 35% rather than tightening the fit around one off-curve point.
