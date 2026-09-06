# AWS handoff (commands only; no cloud actions performed)

## Provenance contract

`model_baseline_sha=ad037ef96f3433a1ad4008994956f89c68a3ed39` identifies the audited design baseline. `execution_source_sha` identifies the later commit that contains this package and is actually built. They are deliberately not required to be equal.

Do not build the old baseline merely to preserve the audit label. Do not build from an uncommitted campaign directory.

## Sequential deployment procedure

1. Copy the package into the canonical repository checkout and commit it.
2. Capture the full execution commit:

```bash
EXECUTION_SOURCE_SHA=$(git rev-parse HEAD)
```

3. Build and push the image from exactly that HEAD. Resolve the pushed image to `ECR_REPO@sha256:<64 hex>` and independently verify that the AWS Batch job-definition revision uses that digest URI, not a tag.
4. Generate deployment artifacts only after both immutable identities exist:

```bash
python3 prepare.py --clean --deployment   --execution-source-sha "$EXECUTION_SOURCE_SHA"   --image-digest "sha256:$DIGEST" [recorded promotion options]
python3 preflight.py --deployment --execution-source-sha "$EXECUTION_SOURCE_SHA"
```

The regenerated inputs/manifests may be uncommitted operational artifacts. They record `model_baseline_sha`, `execution_source_sha`, and `container_image_digest`. Deployment generation/preflight fails without Git, if `HEAD` differs, if baseline ancestry cannot be established, if package source files are uncommitted, or if the digest is mutable/placeholder. A tarball plus baseline retrieval manifest is sufficient only for planning.

5. Generate and review commands for **Q only**:

```bash
python3 aws_commands.py --stage Q   --execution-source-sha "$EXECUTION_SOURCE_SHA"   --image-uri "$ECR_REPO@sha256:$DIGEST"   --input-prefix s3://BUCKET/receptor-v1/Q/inputs   --output-prefix s3://BUCKET/receptor-v1/Q/outputs   --job-queue gutibm-gpu-spot   --job-definition gutibm-cuda-digest:REV
```

The command generator also requires a Git checkout at the declared execution HEAD and matching deployment manifests. It only prints `aws s3 cp` and `aws batch submit-job`; never pipe its output directly to a shell.

## Qualification and sequential gates

Q must show GPU on, one rank, `chemistry_placement=device_delivery`, acceptable repeat behavior/runtime, and HDF5 `/run_provenance/git_sha` exactly equal to `execution_source_sha`. A baseline SHA does not satisfy runtime provenance.

Only after review, create approval JSON for the next stage:

```json
{
  "model_baseline_sha": "ad037ef96f3433a1ad4008994956f89c68a3ed39",
  "execution_source_sha": "<full HEAD used for image>",
  "image_digest": "sha256:<64 hex>",
  "Q_pass": true
}
```

For A–D add `--approval-file approval.json`. Each approval must match both SHAs and the digest. Review A, record Kd, regenerate, then approve A for B; repeat for B and C. Do not queue later arrays speculatively. Timeouts, missing arms, Spot reclaim, CPU fallback, digest mismatch, runtime-SHA mismatch, or missing source-centred transport provenance block promotion. Comparable arms use one digest. D reuses C's plasmid-free null analytically; its `P=0` ColE1 carrier is a distinct no-release control.

Before executing any printed command, review bucket paths, queue, job-definition digest, Identity and Access Management (IAM) role, array size, timeout, and approval.
