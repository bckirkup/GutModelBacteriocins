# Receptor-selection staged campaign v1

This GPU-only package separates two provenance identities:

- `model_baseline_sha=ad037ef96f3433a1ad4008994956f89c68a3ed39` is the audited design baseline. It fixes what was reviewed; it is not the commit that must be built after this package is added.
- `execution_source_sha` is the later commit containing this campaign package. That exact commit must be built, and its full SHA must be emitted as HDF5 `/run_provenance/git_sha`.

The package creates inputs and reviewed commands. It never runs the model or calls AWS.

## Scientific design

The campaign is a sequential, explicit-run ladder, not a sweep product:

| Stage | Jobs | Isolated question | Gate |
|---|---:|---|---|
| Q | 2 | identical fresh-image GPU repeats | provenance, `device_delivery`, repeat behavior, runtime |
| A | 24 | corrected corrinoid–BtuB Kd `{1e-6,1e-5,1e-4,1e-3}` mol/m³ × producer/null × 3 paired seeds | receptor engagement |
| B | 12 | toxin-free BtuB-normal versus BtuB-null over corrinoid `{1e-3,1e-4,1e-5,1e-6}` mol/m³ | viable-growth/action intersection |
| C | 12 | ColE1 mucin-charge amplitude `{0,15,60}` (9 producers + 3 shared nulls) | ordered source-centred transport |
| D | 15 | 12 ColE1 carriers: realized lysis targets `{0,1%,2%,5%}` × 3 seeds, plus 3 same-image plasmid-free nulls | graded realized lysis and producer benefit |

Every run is a 100 µm cube, 2 µm grid (50³), 6 h maximum, CUDA enabled, GPU 0, and one Message Passing Interface (MPI) rank. Only Stage C requests grids: every 60 steps, named species `bacteriocin_BtuB`; provenance is every 10 steps. Other stages set grid output to zero.

## Planning generation

In this retrieved tarball (no `.git`), the baseline retrieval manifest verifies the design baseline. Planning artifacts intentionally use:

- `execution_source_sha=EXECUTION_SOURCE_SHA_REQUIRED_AFTER_CAMPAIGN_COMMIT`
- `container_image_digest=REQUIRED_BEFORE_SUBMISSION`

Generate and inspect them without blocking a future commit:

```bash
cd experiments/receptor_selection_campaign_v1
./prepare_and_preflight.sh
# or: python3 prepare.py --clean && python3 preflight.py
```

Planning preflight passes with an explicit image warning. It is not deployment authorization.

## Required deployment workflow

1. Copy this directory into a normal checkout of the canonical repository whose history contains the audited baseline.
2. Commit the campaign package. Do not build the old baseline and do not leave the package uncommitted.
3. Capture `EXECUTION_SOURCE_SHA=$(git rev-parse HEAD)`.
4. Build and push the runtime image from that exact HEAD; resolve its immutable registry digest as `sha256:<64 hex>`.
5. Regenerate operational deployment inputs/manifests. These generated files may remain uncommitted:

```bash
python3 prepare.py --clean --deployment   --execution-source-sha "$EXECUTION_SOURCE_SHA"   --image-digest "sha256:$DIGEST" [recorded promotion options]
python3 preflight.py --deployment --execution-source-sha "$EXECUTION_SOURCE_SHA"
```

`prepare.py --deployment` fails unless it is in a Git checkout, the supplied execution SHA equals `HEAD`, the audited baseline is an ancestor, the provenance-bearing package files are committed and unchanged, and the image digest is immutable. The manifests and every input record both SHAs. Generated deployment artifacts can differ from HEAD because they are operational records created after the image digest is known.

Submit **Q only**. A–D remain gate-locked. A/B selections are `Kd=1e-4 mol/m³` and `B12=1e-3 mol/m³`. The ecological C population gate passed for amplitudes 0/15/60, but the PR416 intrinsic-transport gate is invalid/failed and superseded. Run the revised single-source assay first. Only its pass, combined with ecological C, supports conditional amplitude 15. Regenerate D with that amplitude on the assay's exact execution SHA and digest; any revision change requires a new assay.

## Outputs and analysis

Place returned files at `generated/stage_X/results/<array-index>/output.h5.gz` (uncompressed `.h5` is accepted), then run `python3 analyze.py`. It writes run and paired metrics plus explicit gate and missing-output JSON. The analyzer requires exact equality between `/run_provenance/git_sha` and `execution_source_sha`; a baseline SHA or prefix match does not pass. It also requires Q `chemistry_placement=device_delivery`.

For Stage C transport, also run:

```bash
python3 analyze_transport.py \
  --results-root generated/stage_C/results \
  --output-dir analysis
```

That writes `transport_profiles.csv`, `transport_metrics.json`, `transport_gate.json`, and `transport_profiles.png`. Kill-event coordinates are read from `/provenance` (runtime HDF5 name); the ecological analyzer accepts `/provenance` or legacy `/kill_provenance`.

Missing outputs are blocked, not inferred as zero. Guard halts remain outcomes. Time points are not replicates; screening summaries are seed-level medians/ranges, not p-values. Stage C blocks if named-grid or source-event-coordinate provenance is absent.

## Package files

- `campaign_contract.json` — science/execution contract and two-SHA provenance policy.
- `campaign_decision_record.json` — authoritative, stage-specific revised C-to-D decisions and execution identities.
- `CURSOR_HANDOFF.md` — exact no-shortcuts assay-to-D operator workflow.
- `prepare.py` — explicit planning or deployment generator; submits nothing.
- `preflight.py` — planning and fail-closed deployment validation.
- `analyze.py` — HDF5, provenance, agent, and gate analysis.
- `analyze_transport.py` — Stage C source-centred bacteriocin transport profiles and gate.
- `aws_commands.py` — prints digest-pinned upload/array commands; executes nothing.
- `AWS_HANDOFF.md` — operator procedure and gates.
- `COMPLETION_NOTE.md` — correction scope and validation status.
- `generated/` — planning artifacts in version control; deployment regeneration may be uncommitted.

## Revised C-to-D status

Stage D is blocked until `single_source_transport_assay_v1` passes and a stage-specific record sets `C_transport_gate=true` for the assay execution SHA and digest. Stage D then uses that same identity and 15 jobs. Its P=0 ColE1 carrier remains a mechanistic no-release arm; it is not the plasmid-free null. The analyzer pairs each of the 12 D carriers only to the D plasmid-free null with the same seed and never searches Stage C controls. `approval.json` is historical and not revised-D authorization.
