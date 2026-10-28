# SPEC-13 step 3 — live-patch placement decision

Deliverable: placement decision for the C/S/R selection campaign —
can one `layer2.audit_patch_index` live patch represent the campaign's
selection contrast, and are ≥2 live patches required?

**Decision: ≥2 live patches are required — in practice the Edison
reference of all-ten-live.**  The lineage readout itself works
(per-lineage alive/emigrant trajectories, closed ledgers, lever
response), but a single live patch cannot hold a selection contrast
under any tested condition and requires a spec change to scale past
one (`audit_patch_index` is a single Int).

Build: serial `build-serial/gut_ibm` at 778c2af.  Runs are local/serial
per SPEC13 §7.  Seeds 20261027–20261031 throughout.  All provenance
files report `ledger_closed: true` (36/36 runs).

## Arms

C/S/R founders co-inoculated into the single live patch
(`layer2.audit_patch_index=0`), 10 founders each, `initial_strains`
order → lineages 1–10=C (ColE1 producer), 11–20=S (sensitive),
21–30=R (`receptor_expression {"BtuB": 0.0}` btuB-null phenocopy).
Mutation and conjugation zeroed so lineage identity = strain identity.
Lysis lever `sos_lysis_prob ∈ {0, 0.02, 0.10}`; retardation lever
`lit` (`ColE1 diff_coeff 7e-11, retardation 1.2`) vs `ship`
(shipped 4e-11, 50.2).  Generators: `make_configs.py`,
`make_incub_configs.py`; per-run outputs under `out/`; machine table
in `readout_summary.json` (`analyze.py`).

Three environmental conditions were exercised:

| condition | parameters | purpose |
|---|---|---|
| churn | contraction 0.5/min, crypt disruption 0.10, agent-loss 0.3, shipped advection | physiological baseline (occupancy-study admissible region) |
| persistent | loss channels zeroed, shipped advection | isolate founder-washout from churn |
| incubator | all loss channels damped (`radial_turnover`/`distal_transit` 1e12), b12=1e-5 | give the contrast its best chance; drop the shipped corrinoid shield (SPEC14 §8 landmine) |

## Measured results

### Churn — total extinction, 15/15 runs (paired seeds)

Every strain in every seed is extinct in ≤3 h of patch time
(crypt; last-alive 0.5–3 h) or effectively instantly (prox; last-alive
0–0.5 h — founders gone before the first 30-min census).  Deaths 0,
reseeds 0–1.  Lysis on/off makes no difference — nothing lives long
enough for the lever to act.  Same founder-lottery failure class as
the rps_probe §4a-result, now measured at patch scale with paired
seeds.  There is no selection window to read.

### Persistent — composition by founder-position washout (smoke, 1 seed)

With churn off but shipped advection, the patch persists, but all 10 C
founders are exported <2 h, R is gone ~4 h, S persists and exports 434
cells.  `deaths=0` — strain fate is decided entirely by where founders
happened to sit relative to the outflow, not by biology.

### Incubator — readout discriminates levers; composition still seed-bound (20 runs, 4 arms × 5 seeds)

Per-lineage alive/emigrant trajectories resolve cleanly.  Lever
response measured in the deaths channel:

| arm | deaths by seed (27/28/29/30/31) |
|---|---|
| lit, lysis 0 | 0, 0, 0, 0, 0 |
| lit, lysis 0.10 | 13, 0, 10, 0, 15 |
| ship, lysis 0 | 0, 0, 0, 0, 0 |
| ship, lysis 0.10 | 42, 0, 44, 0, 281 |

- The **lysis null is a clean falsifier**: deaths=0 in every lysis=0
  run.  The **retardation lever is strongly discriminated**: shipped
  ret=50 holds toxin in the patch and produces up to 21× the deaths
  of the literature analogue on the same seed (281 vs 15).  The ledger
  does not decompose lysis self-kills vs receptor-mediated kills —
  attribution is inferential, not measured.
- The **lysis lever does not move composition inside the window**:
  lysis 0 vs 0.10 give near-identical per-strain trajectories within
  seed (e.g. s20261028 identical 0/41/0 at 8 h, s20261030 identical
  0/0/1).  Kill/lysis events (~10–280 over 24 h) are small against
  ~800 births and are compensated by regrowth — the SPEC14 structural
  break-even ("a killed sensitive is largely replaced") shows up at
  patch scale too.
- Outcome is **seed-determined**: seeds 20261028 and 20261030 collapse
  to near-extinction even quiescent; 27/29/31 bloom S to ~500.
- The **bloom guard freezes the signal**: `patch_cells` crosses the
  10×30=300 bound at ~30 min; the audit record shows n_live=488 vs
  n_shadow=57 at 12 h (the live sim carries ~8× the tabulated capacity,
  rel_diff 7.56), the patch is marked halted (`n_halted_crypt=1`), and
  composition is frozen from ~10–12 h into the 24 h run.  The final
  census row roughly doubles the frozen counts (492→984) — a terminal
  emit artifact worth checking before the campaign.
- **reseeds=0 in every run**: `emigrants_by_lineage` resolves fine
  (per-strain emigrant counts above), but no emigrant lineage ever
  returned to a live patch within 24 h.  Patch→patch lineage spread —
  half of the deliverable's question — is unobservable at n=1: the
  `Single`/`Fragment` packets carry lineage out, and there is nowhere
  live for them to land.

### Structural limits (read from `live_config_for`, segment.cpp)

- `audit_patch_index` is one Int — multi-live is not implementable
  without a spec change.
- `advection.crypts_enabled` is forced false in the live patch, so a
  live "crypt" vs "exposed" patch differs only in `supply_mult`
  (carbon boundary 5e-3×mult).  The crypt/exposed *flow* contrast is
  tabulated-only; measured corroboration: prox founders never survive
  the first census under churn.

## Interpretation

Measured: the instrumentation is good — lineage-resolved outputs
discriminate the retardation arm and cleanly separate the null arm.
Also measured: a single live patch either goes extinct (churn), is
decided by founder-position washout (persistent), or saturates and is
frozen by the bloom guard before selection can reorder composition
(incubator) — and the pool→live reseed channel is never exercised at
n=1.

Inferred: the campaign's displacement contrast (Kirkup & Riley-scale,
weeks of turnover) is an ensemble property — sustained residence plus
emigrant→reseed lineage recycling across patches of different type.
No single patch can express that: it is simultaneously the whole
source and the whole sink of every lineage.  The observed freeze at
the bloom bound means even a quiescent single patch yields only ~10 h
of effective dynamics per 24 h run.

Hypothesis (not measured here): with all ten patches live, each patch
still churns (~2–3 h residence) but the ensemble retains lineages via
pool→live reseeds, which the code path already supports
(`dest.live` branch in `pool_phase`).  Whether that retention is fast
enough to keep three strains co-resident for a discrimination window
is the first thing the all-ten-live reference must measure.

## Spec deltas needed for ≥2 live patches

1. `layer2.audit_patch_index` (Int) → a list form
   (`audit_patch_indices` or `n_live_patches`), keeping the single-int
   alias for back-compat.
2. Per-patch-type live advection: `live_config_for` currently forces
   `crypts_enabled=false` for every live patch; the crypt/exposed
   contrast needs the live patch to keep its type's flow profile.
3. Revisit the live-patch bloom bound (baseline 30 founders × 10):
   the live sim's own carrying capacity is ~500 cells, so the guard
   trips well before the live patch is actually saturated.
4. Decompose the ledger `deaths` channel (lysis self-kill vs
   receptor-mediated kill) — currently the lever is measurable but
   not attributable from outputs alone.

## Validation gate status

- Selection signal measured per lineage with paired seeds: **met**
  (5 seeds, crypt+prox × ret × lysis under churn; 4 arms × 5 seeds
  under incubator; signal = extinction/saturation, see above).
- Ledger closed: **met** (36/36 `ledger_closed: true`).
- Config/test additions green on CI: analysis-only files
  (`analysis/spec13_placement/`); no source or schema changes.

## Files

- `make_configs.py` — churn-arm generator (40 configs + manifest).
- `make_incub_configs.py` — incubator-arm generator (20 configs).
- `analyze.py` → `readout_summary.json` — per-arm machine table.
- `configs/`, `out/` — configs and raw outputs (lineage.csv, prov.json,
  ts.csv per run).
