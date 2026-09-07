# CURSOR handoff: authoritative revised C-to-D workflow

This patch is planning-only. **Do not submit AWS jobs while applying it.**
`campaign_decision_record.json` is the machine-readable authority for the revised handoff. The existing `approval.json` is historical and does not authorize the revised assay or Stage D.

## Exact workflow

1. Apply this patch to canonical `main` at merge SHA `f9a908e0eb5323cfaeea82501a8ca16b1a95910c`, commit it, and merge it. Do not build from the pre-patch SHA or an uncommitted tree.
2. In a clean checkout of the resulting canonical HEAD, record `EXECUTION_SOURCE_SHA=$(git rev-parse HEAD)`.
3. Build and push the runtime image from that exact HEAD. Resolve and record its immutable `sha256:<64 hex>` digest. The checkout must remain at the same HEAD.
4. Generate the **12-job** `single_source_transport_assay_v1` deployment package with that exact execution SHA and digest; run its preflight; review commands; then submit only that assay.
5. Return all 12 outputs and run `analyze_assay.py`. Missing/authentication-invalid outputs are blockers. Do not treat the superseded PR416 intrinsic-transport result as a pass.
6. If and only if `assay_gate.json` records an assay pass, record `C_transport_gate=true` for this exact execution SHA and digest. In combination with the already-passed ecological C population gate, record amplitude **15** as the conditional Stage D choice. The assay alone does not calibrate amplitude.
7. Regenerate the receptor-selection campaign deployment artifacts with `--selected-kd 1e-4 --selected-b12 1e-3 --selected-amplitude 15`, using the **same exact** `EXECUTION_SOURCE_SHA` and image digest as the assay. Run deployment preflight.
8. Verify Stage D contains exactly 15 explicit jobs: 12 ColE1-carrier jobs (four unchanged lysis targets × three seeds, including the distinct `P=0` carrier) plus three plasmid-free nulls (one per seed), with no null × lysis factorial.
9. Supply a stage-specific approval that sets `C_transport_gate=true` and matches the exact SHA and digest. Only then may Stage D be considered for submission. Any code change, new HEAD, or new image digest invalidates this handoff and requires rerunning the assay before D.

## Required commands (operator fills registry details)

```bash
# after commit/merge, in a clean canonical checkout
EXECUTION_SOURCE_SHA=$(git rev-parse HEAD)
# build image from exactly $EXECUTION_SOURCE_SHA and resolve IMAGE_DIGEST

cd experiments/single_source_transport_assay_v1
python3 prepare_assay.py --clean --deployment   --execution-source-sha "$EXECUTION_SOURCE_SHA"   --image-digest "<repository>@sha256:<64 hex>"
python3 preflight_assay.py --deployment --execution-source-sha "$EXECUTION_SOURCE_SHA"
# review and submit the 12 assay jobs only; then return outputs
python3 analyze_assay.py

# only after a recorded assay pass
cd ../receptor_selection_campaign_v1
python3 prepare.py --clean --deployment   --execution-source-sha "$EXECUTION_SOURCE_SHA"   --image-digest "sha256:<same 64 hex>"   --selected-kd 1e-4 --selected-b12 1e-3 --selected-amplitude 15
python3 preflight.py --deployment --execution-source-sha "$EXECUTION_SOURCE_SHA"
```
