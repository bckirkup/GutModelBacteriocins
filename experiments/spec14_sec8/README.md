# Spec 14 §8 validation and the corrinoid ladder

Config generators and analysis for the §8 no-code validation reported in
`docs/SPEC14_SECTION8_VALIDATION.md`. Outputs (HDF5, run directories, logs) are
not committed, per the repository rule on campaign artifacts.

Run at `5fb4250`, serial, local, no code changes.

## Scripts

| File | Role |
|---|---|
| `prepare.py` | the §8 pair itself: in-vitro (6 h, well-mixed 100 µm cube) and in-vivo (10 d, 200×200×100 µm with peristalsis, crypts, motility, gradient), producer vs null, five paired seeds |
| `prepare_corrinoid.py` | the follow-up ladder: `b12_initial_conc` ∈ {1e-3, 1e-4, 1e-5, 1e-6, 1e-9} × {producer, null} × 5 seeds, in-vitro only |
| `analyze.py` | per-arm counts, events, termination cause, provenance SHA; paired producer − null |
| `analyze_corrinoid.py` | per-level paired selection coefficient, kills per lysis, kills per division |

The analysis scripts import `h5py`, which the system `python3` does not have.
Run them with the venv the blueprint's maintenance step builds:

```bash
python/.venv311/bin/python experiments/spec14_sec8/analyze_corrinoid.py --run-root <campaign-root>
```

Run `analyze_corrinoid.py` from this directory; it imports helpers from
`analyze.py` and `prepare_corrinoid.py` imports the shared config blocks from
`prepare.py`, so neither is self-contained by design — the ladder must not be
able to drift from the §8 configuration it extends.

## Arms

`prepare.py` generates paired competition and control arms across `vitro` and `vivo`:
- `producer`: type 1 carries `ColE2` (`BtuB` target, nuclease, `SOS_LYSIS`), type 2 is plasmid-free sensitive.
- `null`: plasmid-free reference null (`plasmids=[]` on both strains).
- `kill_null`: plasmid-carriage-matched kill-null control (type 1 carries `ColE2` with identical metabolic carriage and SOS lysis cost; `kill_rate_colicin = 0`, `kill_rate_microcin = 0`, `sos_cross_induction_rate = 0`). Pairing `producer - kill_null` isolates receptor-mediated killing from plasmid carriage/lysis burden, while `kill_null - null` isolates producer carriage/lysis cost.
- `immunity_null`: carriage- and immunity-matched control (both type 1 and type 2 carry `ColE2`).
- `colIa_producer` and `colIa_kill_null`: non-`BtuB` comparator targeting `CirA` (`ColIa`, paired with its carriage-matched kill-null).
- `colB_producer` and `colB_kill_null`: non-`BtuB` comparator targeting `FepA` (`ColB`, paired with its carriage-matched kill-null). (`ColE1` also targets `BtuB` and is never used as a non-`BtuB` control.)

Optional CLI flags (`--no-b12-growth-essential`, `--btub-competitor-conc`, `--b12-free-fraction`) decouple extracellular `BtuB` competitive corrinoid shielding from obligate `B12` Monod growth limitation.

## Assertions

- **V1** in vitro, the producer excludes the sensitive by 100× (2 decades)
  within 6 h, i.e. ≥ 0.33 decades/h.
- **V2** in vivo, the two persist at comparable levels over 10 days.
- **V3** the V1 contrast is attributable to receptor-mediated colicin kills,
  reported as kills per division and kills per producer lysis.

**V1 FAILED at every corrinoid level** (best 0.090 decades/h, 3.7× short).
**V2 was not obtained** — the in-vivo arm was abandoned at 1.2 of 10 days with
the population still climbing, because an in-vivo arm cannot demonstrate the
reversal of an exclusion that never happened.

A FAIL is a result about the model. Read
`docs/SPEC14_SECTION8_VALIDATION.md` before quoting any row: 30 of 50 ladder
arms ended at the dysbiosis guard rather than the horizon, the 1e-6 and 1e-9
levels are corrinoid-growth-limited rather than comparable cultures, and a
ratio read off the final output step measures the phase of the division wave
rather than selection.
