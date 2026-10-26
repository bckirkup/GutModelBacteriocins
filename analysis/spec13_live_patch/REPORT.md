# Spec 13 — Layer-2 mechanism-capable live patch (Edison step 2)

One `layer2.audit_patch_index` slot runs a real `Simulation` inside a
tabulated `MucusSegment`: correct influx/efflux, biological-death
accounting, lineage-resolved outputs, and a minimal segment→chain handoff.
Edison step-1 settled input: founder genotype propagates via
`apply_founder_genotype` (#452); the tabulated layers refuse
mechanism-dynamics configs except this slot.

## Verdicts (measured)

- **V1 PASS — readout matches a Layer-1 reference.** A standalone
  `gut_ibm` run with the live patch's serialized resolved config
  (`configs/l1_ref.json`, built from `ConfigJson::serialize_document` of
  the live cfg, `hdf5` enabled) reproduces the segment's live patch
  **bit-for-bit** through 10 bio-steps: positions, phenotype state,
  `mu_realized`, carbon and iron voxel fields identical at all four module
  boundaries (biology / chemistry / advection / post-physics); RNG engine
  states identical. Lineage-alive counts in
  `layer2.audit_lineage_file` (`out/iso_lineage.csv`) match
  `agents/step_*/lineage_id` in the reference HDF5 with **0 mismatches**.
- **V2 PASS — cell ledger closes.** `out/iso_prov.json`:
  initial=20, births=0, biological deaths=11, distal departures=9,
  washout=9, final stocks=0 → 20 = 11 + 9, closed every step
  (`population_stop` at t=1080 s; deaths verified independent of the
  contraction/washout channels).
- **V3 PASS — ledger closes under coupling arms** (probes in
  `tests/test_layer2_segment.cpp`): lysis-lever arm books
  `deaths=9` under `sos_basal_rate`; influx arm runs to horizon with
  pool reseeds entering the live patch alive.

## Defects found and fixed (this PR)

1. **Shadow-birth leak** — the tabulated twin's `grow_tabulated`
   unconditionally added fictional births to `ledger_.births` beside the
   live sim's real divisions; now called with `ledger=false` for live
   patches (delta from the sim's `divisions` counter is booked instead).
2. **Biological deaths invisible** — new `biological_deaths` ledger
   channel harvested from StepEvents mortality deltas
   (`step_events + cumulative_events`; cumulative alone lags the
   summary-commit boundary) and subtracted in `expected_stocks()`.
3. **Corpse counting** — live `patch_count`/`audit_check`/shadow reseed
   counted DEAD agents retained by the CDI corpse window; all three now
   select `state != DEAD`.
4. **Dead-cell harvest** — contraction hook now selects among alive-only
   agents so a biology-killed cell cannot be booked twice (death +
   fragment departure).
5. **DEAD-before-hook packet export** — `check_washout` marked agents
   DEAD before firing `departure_hook_`; exported packets entered the
   pool as corpses and were reaped unbooked on reseed. Hook now fires
   first at both export sites so packets carry live state.
6. **Carbon spec stale after supply rescale** — `live_config_for`
   rescaled `carbon_boundary_conc` by `tc.supply_mult`, but the parser had
   already baked `chemicals[carbon].boundary_conc` and `initial_conc`
   from the *unscaled* boundary; the live patch's carbon field was half
   the intended crypt-supply profile. The parser's rule is now reapplied
   (`boundary_conc` always; `initial_conc` unless an explicit
   `carbon_z_amplitude` override is set). This was the root cause of the
   readout's growth-rate divergence — `FixMotility` energy taxis scales
   swim speed by `mu_realized`, so a halved carbon field moved motility,
   not just metabolism.

## Outputs and handoff

- `layer2.audit_lineage_file` → per-step, per-lineage
  `alive_count,emigrants_cum` CSV (`kLiveLineageHeader`).
- Provenance gains `"deaths"` and `"audit_live": {patch, live_seed}`
  fields; timeseries gains `deaths_cum`; checkpoint carries
  `/ledger/biological_deaths`.
- Segment→chain handoff: `departure_hook_` singles + harvested fragments
  reach `push_packet` as live `Fragment`/`Single` packets
  (lineage-resolved); `try_deposit` still refuses live destinations
  (spec), and `set_pool_exit_sink → LumenCompartment::add_cells` is the
  chain-facing surface a Layer-3 arm wires without touching
  `ColonicChain`'s deliberate mechanism refusal.

## Per-patch cost (measured)

Same geometry/seed, `total_time=3600 s`, `bio_dt=60 s`, 2-core VM:

| Arm | Wall | Steps | Cost |
|---|---|---|---|
| Tabulated-only segment (audit off) | ~0.00 s | 60 | unmeasurable at this size |
| Segment + 1 live patch | 2.37 s | 60 | **~39.5 ms/step per live patch** |
| Standalone identical sim | 0.75 s | 10 | ~75 ms/step (higher sustained load) |

RSS delta ≈ +19.7 MB per live patch (the patch carries a full
`Simulation` + `ChemicalField`). Cost scales with agent load like the
standalone engine; the segment adds only packet push/harvest overhead.

## Reproduce

```bash
build-serial/gut_ibm analysis/spec13_live_patch/configs/iso_arm.json
build-serial/gut_ibm analysis/spec13_live_patch/configs/l1_ref.json
python3 analysis/spec13_live_patch/run_readout.py \
  out_dir/iso_lineage.csv l1_ref.h5
```
