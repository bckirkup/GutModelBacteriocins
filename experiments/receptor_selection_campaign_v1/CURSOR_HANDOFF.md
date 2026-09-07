# CURSOR handoff: authoritative revised C-to-D workflow

`campaign_decision_record.json` is the machine-readable authority. The existing
`approval.json` is historical and does not authorize the revised assay or Stage D.

## Status (2026-09-07)

- Ecological Stage C population gate: **PASS** (SHA `3f176b26…`).
- Single-source transport assay first Batch attempt (SHA `d7b16c3`, job
  `231c8870…`): **BLOCKED** — science ordered; authentication failed because
  inputs omitted `metabolism.uptake_limit=delivery` (provenance wrote `device`).
  Evidence: `experiments/single_source_transport_assay_v1/analysis/`.
- `prepare_assay.py` now sets delivery. **Do not** record `C_transport_gate` or
  submit Stage D until a new assay pass on a post-fix merge SHA + matching
  image digest.
- Historical `d7b16c3` / `:8` image is **invalidated** for D by this prep fix.

## Exact workflow (current)

1. Merge the delivery-auth fix to canonical `main`. Do not build from an
   uncommitted tree.
2. In a clean checkout of the resulting canonical HEAD, record
   `EXECUTION_SOURCE_SHA=$(git rev-parse HEAD)`.
3. Build and push the runtime image from that exact HEAD. Resolve and record its
   immutable `sha256:<64 hex>` digest. The checkout must remain at the same HEAD.
4. Generate the **12-job** `single_source_transport_assay_v1` deployment package
   with that exact execution SHA and digest; confirm every input has
   `metabolism.uptake_limit=delivery`; review commands; then submit only that
   assay (explicit authorization required).
5. Return all 12 outputs and run `analyze_assay.py`. Missing or
   authentication-invalid outputs are blockers. Do not treat the superseded
   PR416 nearest-source result, or the blocked `d7b16c3` Batch attempt, as a pass.
6. If and only if `assay_gate.json` records an assay pass, record
   `C_transport_gate=true` for this exact execution SHA and digest. In
   combination with the already-passed ecological C population gate, record
   amplitude **15** as the conditional Stage D choice. The assay alone does not
   calibrate amplitude.
7. Regenerate the receptor-selection campaign deployment artifacts with
   `--selected-kd 1e-4 --selected-b12 1e-3 --selected-amplitude 15`, using the
   **same exact** `EXECUTION_SOURCE_SHA` and image digest as the assay. Run
   deployment preflight.
8. Verify Stage D contains exactly 15 explicit jobs: 12 ColE1-carrier jobs (four
   unchanged lysis targets × three seeds, including the distinct `P=0` carrier)
   plus three plasmid-free nulls (one per seed), with no null × lysis factorial.
9. Supply a stage-specific approval that sets `C_transport_gate=true` and matches
   the exact SHA and digest. Only then may Stage D be considered for submission.
   Any code change, new HEAD, or new image digest invalidates this handoff and
   requires rerunning the assay before D.

## Required commands (operator fills registry details)

```bash
# after commit/merge, in a clean canonical checkout
EXECUTION_SOURCE_SHA=$(git rev-parse HEAD)
# build image from exactly $EXECUTION_SOURCE_SHA and resolve IMAGE_DIGEST

cd experiments/single_source_transport_assay_v1
python3 prepare_assay.py --clean --deployment \
  --execution-source-sha "$EXECUTION_SOURCE_SHA" \
  --image-digest "<repository>@sha256:<64 hex>"
# confirm delivery key is present in generated/jobs/*/input.json
# review and submit the 12 assay jobs only; then return outputs
python3 analyze_assay.py

# only after a recorded assay pass
cd ../receptor_selection_campaign_v1
python3 prepare.py --clean --deployment \
  --execution-source-sha "$EXECUTION_SOURCE_SHA" \
  --image-digest "sha256:<same 64 hex>" \
  --selected-kd 1e-4 --selected-b12 1e-3 --selected-amplitude 15
python3 preflight.py --deployment --execution-source-sha "$EXECUTION_SOURCE_SHA"
```
