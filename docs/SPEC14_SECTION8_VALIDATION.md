# Spec 14 §8 validation: the in-vitro half of the reversal does not occur

Spec 14 proposes accumulated DNA damage, a lysogen state, and a Hill induction
hazard, motivated by Henrot et al. (bioRxiv 10.64898/2026.05.26.727859): an
E-type endonuclease colicin producer **outcompetes a λ-lysogen within 6 h in
vitro**, yet both persist at comparable levels over 10 days in dixenic mice.
That in-vitro-wins / in-vivo-coexists reversal is what the model is supposed to
explain, and the accompanying brief is explicit that §8 must be run **first,
with no code changes**, because if the reversal falls out of spatial refugia
and density limitation alone then no new mechanism is licensed.

It does not fall out, and the reason is upstream of anything Spec 14 proposes:
**GutIBM does not currently reproduce the in-vitro half.** There is no
exclusion to reverse. This document records what was measured, why, and what it
does and does not license.

Run at `5fb4250`, no code changes. Configs and analysis are in
`experiments/spec14_sec8/`; HDF5 outputs are not committed.

## 1. Result

Well-mixed analogue: 100 µm cube, 2 µm grid, no gradient, no peristalsis, no
crypts, no motility, Dirichlet carbon at 0.05 mol/m³, 60 ColE2 producers + 60
sensitives, 6 h horizon, five paired seeds. Every arm carries a **null**: the
identical founders and seed with the ColE2 plasmid removed from type 1, so that
plasmid carriage cost, corrinoid-limited growth, and crowding all cancel in the
paired difference and only the colicin contributes.

| `b12_initial_conc` | apparent Kd | selection, decades/h | kills per producer lysis | kills/division | final agents | termination |
|---|---|---|---|---|---|---|
| 1e-3 **(shipped)** | 5.01e-4 | **−0.036** | 0.121 | 0.0013 | 28888 | guard |
| 1e-4 | 5.05e-5 | **−0.052** | 1.17 | 0.0112 | 17873 | guard |
| 1e-5 | 5.50e-6 | **+0.090** | 7.38 | 0.092 | 23185 | guard |
| 1e-6 | 1.00e-6 | **+0.081** | 6.29 | 0.091 | 1651 | horizon |
| 1e-9 | 5.00e-7 | +0.000 | 1.5 | — | 120 | horizon |

Medians over five paired seeds; every level's null returns 0.000 ± 0.004
decades/h, which is the control the contrast is read against.

**V1 FAILS at every corrinoid level.** The Spec 14 target is 100× exclusion
within 6 h, i.e. 2 decades in 6 h, i.e. 0.33 decades/h. The best the model
achieves anywhere on the ladder is 0.090 decades/h — at that rate the claimed
endpoint arrives in **22 h, not 6**, and only in a culture that has already
been stopped at the dysbiosis guard four times over. This is a shortfall in
kind, not a near miss to be tightened up.

**V2 is not reported.** The in-vivo arm was started and abandoned: it reached
28,085 agents at 1.2 of its 10 days with the population still climbing, and the
projected serial cost was not defensible against a V1 that had already failed.
An in-vivo arm cannot demonstrate a *reversal* of an exclusion that never
happened, so it was not worth its wall time before V1 is addressed.

## 2. Why: BtuB is corrinoid-saturated at the shipped pool

The shipped result — 299 producer lyses buying 36 kills, a strategy running at
**0.12 kills per certain producer death** — is not a transport failure. ColE2's
retardation is `retardation_from_pI(6.5)` = 2.04, not the 50 that ColE1 (pI 9.0)
carries, so the nuclease colicin is nearly unretarded and a 5e4-molecule burst
over `burst_release_tau` holds ~1.3e-6 mol/m³ at 1 µm and ~2.6e-7 at 5 µm, at or
above `kd_colicinE_btuB = 5e-7` across a plume containing hundreds of cells.
Two orders of magnitude of kills are missing against that estimate.

They are missing into competitive binding. `FixReceptor::toxin_occupancy` is

```text
apparent_kd = kd_colicinE_btuB * (1 + [corrinoid] / kd_b12_btuB)
```

and the shipped field is `b12_initial_conc = 1e-3 mol/m³` (1 µM) against
`kd_b12_btuB = 1e-6` (1 nM). BtuB is corrinoid-occupied 1000-fold, so colicin E
faces an apparent Kd of **5e-4 mol/m³ rather than 5e-7** and the plume that
should be saturating is three orders below half-occupancy. The header comment on
that constant already names it "the key unknown governing colicin-E competition
against the ~1 µM corrinoid pool"; the ladder measures what it governs, and the
answer is a **61× swing in kills per lysis across a 91× swing in apparent Kd**,
which is the competitive model behaving exactly as written.

So the shipped configuration is not one in which colicins matter. That is a
parameter position, not a defect: whether ~1 µM *free* corrinoid is the right
mucus value is a scientific question this campaign does not settle, and nothing
here was retuned to make the producer win.

## 3. Why it still fails where colicin is potent

Fixing potency does not buy exclusion, and the ladder says why in one number:
at `b12 = 1e-4` the producer kills **1.17 sensitives per lysis** — ten times the
shipped efficacy — and is nonetheless the **worst** arm on the ladder at −0.052
decades/h, worse than the near-inert shipped level.

Break-even is not one kill per lysis. A lysis is a certain death that forfeits
the producer's entire future lineage, while a killed sensitive in a
density-capped culture is substantially replaced by its neighbours' growth into
the space it vacated. The ladder brackets the break-even between 1.17 (net
negative) and 7.38 (net positive) kills per lysis, so the strategy needs to be
**several-fold better than one-for-one** merely to stop losing. This is a
structural statement about SOS-release bacteriocins in a saturating batch
culture and it does not depend on the corrinoid value.

Two boundaries close the window from the other side. Corrinoid is also a growth
substrate against `km_b12 = 1e-6`: at 1e-6 the culture reaches 1651 agents
instead of ~29000, and at 1e-9 it does not grow at all — 120 agents, the
founders, unchanged over 6 h. Colicin potency and B12-dependent growth are
driven by the *same pool*, so the lever cannot be pushed to maximum potency
without extinguishing the culture that is supposed to demonstrate exclusion. The
only level that is both potent and growing is 1e-5, and 1e-5 is the +0.090
decades/h that still misses by 3.7×.

## 4. What this licenses

- **Not** Spec 14 Change 1 (accumulated DNA damage) or Change 2 (lysogeny).
  Their motivation is a reversal whose in-vitro half the model does not produce.
  Calibrating a damage/induction mechanism against this target now would tune
  new parameters to compensate for whatever is actually suppressing the
  exclusion — and the parameters in question have correlated provenance (De
  Paepe 2016 and Henrot 2026 share a senior author), so they are the last thing
  that should absorb an unexplained residual.
- **Not** Change 3 (free phage), which was already deferred.
- The corrinoid pool is now a **first-class parameter decision**, ahead of any
  phage work: `b12_initial_conc` moves colicin efficacy by 61× and it is
  currently set where colicins do not matter.
- The break-even result stands on its own and is testable independently of
  phage: an SOS-release bacteriocin must clear several kills per lysis before it
  pays, which is a quantitative constraint on the lysis prior that Spec 13
  selects.

## 5. Method note: do not read a ratio off the final sample

`log10(n1/n2)` at the last output step is not a usable statistic here. Once the
culture is dense the two strains divide in **alternating synchronized waves**,
and the last eight samples at `b12=1e-5` run

```text
1.10  1.12  0.68  1.18  1.27  1.30  0.78  1.44
```

so the final sample reports +0.16 decades where the wave-averaged value is
+0.065. A statistic sampled at one instant of a synchronized population
measures the phase of the wave. Every number in §1 is therefore a tail median
and a fitted slope over multiple waves, and the same correction should be
applied to any other campaign analyzer that quotes a final-step ratio.
