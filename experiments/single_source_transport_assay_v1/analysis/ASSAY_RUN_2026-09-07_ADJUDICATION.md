# Adjudication of the first single-source assay attempt (2026-09-07)

This record is separate from, and does not overwrite, the original blocked artifacts in this directory.

## Scope and evidence

The immutable historical record is `ASSAY_RUN_2026-09-07.md`, `assay_gate.json`, `assay_metrics.json`, and `assay_snapshots.csv` for execution SHA `d7b16c30d4afdd3a317133e4db13ba8686a2027b`, image digest `sha256:786576980ca74d836a4cf006c66a65cfd33fa6d8f7311fdf77cc77384020d9d8`, and Batch job `231c8870-2564-4734-ab8c-1359f9a62afa`.

## Corrected placement interpretation

The original analyzer expected `device_delivery`, inherited from an ecological configuration that exercised metabolism delivery. That criterion was irrelevant to this intrinsic assay, whose declared mechanisms are only bacteriocin and receptor. Physical GPU chemistry placement is `device`; metabolism must not be enabled merely to change the provenance label. The observed `device` placement is therefore the correct authenticated placement under the corrected assay contract. `host` and `host_forced_delivery` would not be acceptable.

## Scientific adjudication

The historical artifacts report that all 12 outputs returned, each producer had one source, source positions and 36 paired times were complete, toxin-free nulls were zero, and the original endpoint criteria passed at every paired time: `r50/r90` ordered 0 > 15 > 60 with the declared margins and near field ordered 0 < 15 < 60. The 12 identical placement mismatch strings were the only blockers. I therefore adjudicate the existing d7b16c HDF5 results as scientifically supportive, and authenticated for physical placement under the corrected `device` contract.

This adjudication does not rewrite the historical gate status and does not authorize `C_transport_gate` or Stage D. The files lack the issue-420 interpolation amplitudes 20 and 30. A new 18-job adaptive intrinsic assay is required for the current formal gate. The separate 12-job ecological 15/20/30 refinement is also pending and must remain separate. Stage D remains blocked and no D amplitude selection changes here.
