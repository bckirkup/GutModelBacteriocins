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

## Running

```bash
python3 experiments/rps_probe/prepare.py    # writes arm directories
# run each arm's config.json with a serial Release + HDF5 build
python3 experiments/rps_probe/analyze.py    # asserted contrasts
```
