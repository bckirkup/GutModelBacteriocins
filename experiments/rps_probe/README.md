# btuB-null resistant-strain probe

Population-scale check of the claim in `docs/SPEC13_LYSIS_SELECTION.md` §4a
that a three-strain C/S/R system is constructible today through per-strain
`receptor_expression`, and that a `{"BtuB": 0.0}` founder is both colicin-proof
and measurably costly. Both halves are read from source; #334 established that
per-agent reasoning about a hazard does not transfer to a population, so the
claim is not usable until a run asserts it.

Config generator and analysis are committed. Outputs (HDF5, run directories,
logs) are not, per the repository rule on campaign artifacts.

## Arms

12 arms, serial, local, 12 h horizon in a 200×200×100 µm patch at 2× the
shipped carbon amplitude — the corrected ladder
(`docs/CARBON_LADDER_CAMPAIGN.md`) puts capacity at 122 agents from 100
founders at 1×, which leaves no headroom for a fitness difference to appear in.

| Group | Arms | Founders | Question |
|---|---|---|---|
| `A_three_strain_s*` | 3 seeds | 40 C (ColE1) + 40 S + 40 R | does the knockout survive a producer that kills S? |
| `A_null_no_producer_s*` | 1 seed | same three founders, no ColE1 | is the A separation actually colicin? |
| `B_b12_*_s*` | 4 corrinoid levels × 2 seeds | 60 S + 60 R, no producer | what does the knockout cost, and where does the cost become real? |

R is declared `"receptor_expression": {"BtuB": 0.0}`. Type is the strain label
throughout; `PhenoState::RESISTANT` is **not** usable for strain accounting
because the same field is later overwritten by `SOS_INDUCED`/`DEAD`.

## Assertions

- **A1** producer present ⇒ final R/S > 1.
- **A2** the contrast vanishes in the null arm.
- **B** R/S is non-increasing as `b12_initial_conc` falls
  `1e-3 → 1e-4 → 3e-5 → 1e-5`.

A FAIL is a result about the model. Nothing in these configs may be retuned to
make a contrast pass, and every row is reported with its termination cause and
delivery rationing factor so that a crowded or early-terminated arm is not
quoted as an ordinary horizon value.

## Result at `3c166c1`: A1 FAIL, A2 FAIL, B PASS-then-withdrawn

The probe could not answer its question, and why is the finding. Full account
in `docs/SPEC13_LYSIS_SELECTION.md` §4a-result; in brief:

- Outcome is a founder lottery, not a treatment effect. 120–180 founders fall
  to ~30 live within ~2 h, then the patch either escapes (2383–7629 divisions)
  or does not (22–297), and the survivors sweep it. Three identical-treatment
  seeds ended R-only, S-only, and C+S.
- The null arm — no ColE1 anywhere — gave the largest R/S separation in the
  probe (395 vs 0.67 with a producer). That is A2 failing outright.
- The producer treatment is nearly inert: 11 receptor-mediated kills against
  7629 divisions, consistent with the hardcoded `retardation = 50`. A1's FAIL
  therefore says nothing about the `{"BtuB": 0.0}` genotype.
- B's PASS was an artifact of `n3 / max(n2, 1)` with `n2 = 0`; one strain was
  extinct in 8 of 12 arms. `analyze.py` now excludes arms without both strains
  and reports `INSUFFICIENT` rather than averaging a ratio against an extinct
  denominator. Re-run on the same outputs, nothing passes.

A probe with power needs the patch off its self-sustainment threshold (Layer 2
reseeding or Layer 3), ≥5 seeds per arm with a paired within-seed statistic
rather than a ratio of between-arm means, and colicin transport fixed first.

## Running

```bash
python3 experiments/rps_probe/prepare.py    # writes arm directories
# run each arm's config.json with a serial Release + HDF5 build
python3 experiments/rps_probe/analyze.py    # asserted contrasts
```
