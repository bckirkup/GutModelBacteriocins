# CURSOR handoff: corrected adaptive C-to-D workflow

`campaign_decision_record.json` is authoritative. Historical `approval.json` does not authorize the adaptive intrinsic assay, ecological refinement, or Stage D.

## Status (2026-09-07)

- Ecological Stage C population gate at SHA `3f176b26…`: **PASS** for the original 0/15/60 design.
- Historical 12-job intrinsic attempt at `d7b16c30…`: recorded **BLOCKED** only by the inherited `device_delivery` criterion. Its observed `device` provenance is the correct physical GPU placement for a bacteriocin+receptor assay with no metabolism. All then-declared scientific criteria passed. Preserve its artifacts; see the separate adjudication record.
- Adaptive intrinsic assay (issue 420): **PENDING**, 18 jobs (`SS_amp` 0/15/20/30/60 × 3 plus 3 nulls). It requires new runs because d7b16c lacks 20/30.
- Separate ecological refinement (issue 420): **PENDING**, 12 jobs (`C_amp` 15/20/30 × 3 plus 3 same-image nulls). Never submit it through the intrinsic assay package.
- Stage D: **BLOCKED**. Do not alter D selection yet.

## Required intrinsic workflow

1. Use a clean canonical checkout at the exact execution HEAD, descending from PR417 event-time SHA `f9a908e0eb5323cfaeea82501a8ca16b1a95910c`.
2. Build/push from that exact HEAD and resolve `<repository>@sha256:<64 hex>`.
3. Generate the 18-job deployment package with `prepare_assay.py --deployment`. Inputs must contain only `fixes=[bacteriocin,receptor]`, no `metabolism.*`, and expected placement `device`.
4. Run `preflight_assay.py --deployment`; it checks exact HEAD, clean tracked package, PR417 ancestry, immutable image, hashes, IDs, contiguous jobs, configuration, pairing design, and isolated outputs.
5. Review output from `aws_commands_assay.py`. It never executes commands. Submit only with explicit human authorization outside these scripts.
6. Analyze all 18 returned outputs. The original endpoint margins remain on 0>15>60; 20/30 add monotonic interpolation only. Any missing amplitude/time blocks.
7. Keep ecological refinement separate. Do not infer its result from the intrinsic assay.
8. Do not release Stage D or change amplitude selection until both required issue-420 evidence streams and a new explicit decision record support it.

## Planning commands

```bash
cd experiments/single_source_transport_assay_v1
python3 prepare_assay.py --clean
python3 preflight_assay.py
```

Planning placeholder mode may pass with warnings but is never deployable.
