# GitHub → AWS OIDC setup (dispatch + collector)

One-time setup giving GitHub Actions short-lived AWS credentials with **no
stored keys anywhere**. GitHub mints an OIDC token per workflow run; AWS
trades it for a scoped role. Two roles:

- `github-gutibm-dispatch` — Batch submit/describe on `gutibm-*`
  queues/jobdefs + S3 read/write on `campaign/*` in both
  `gutibm-outputs-994254241749` and `gutibm-inputs-994254241749`
  (receipts, manifests, job inputs).
- `github-gutibm-collector` — read-only S3 on `campaign/*` in both buckets +
  Batch describe/list. Cannot delete, restore, or write anything.

Both require the GitHub environment `aws-campaign`, so every run must pass
that environment's protection rules (e.g. required reviewer approval)
before any AWS credential is issued.

## Step 1 — create the AWS resources (one command)

In the AWS console (account `994254241749`), open **CloudShell** and run:

```bash
curl -sL https://raw.githubusercontent.com/bckirkup/GutModelBacteriocins/main/deploy/aws/github_oidc_stack.yaml -o /tmp/stack.yaml && aws cloudformation deploy --stack-name github-oidc-gutibm --template-file /tmp/stack.yaml --capabilities CAPABILITY_NAMED_IAM --parameter-overrides CreateOidcProvider=no NotificationEmail=YOUR_EMAIL
```

(replace `YOUR_EMAIL`, or drop `NotificationEmail=YOUR_EMAIL` to skip the
monthly Batch+S3 budget alert — capped at $500 with 80%/100% email
notifications).

`CreateOidcProvider=no` is required: the account already has the
`token.actions.githubusercontent.com` OIDC provider (one per URL, created
for the Crusher bridge) and a second create fails on
`AWS::EarlyValidation::ResourceExistenceCheck`. `--template-file` needs a
real file, so the `curl` lands it first — never `<(curl …)`.

The command finishes with two outputs:

- `DispatchRoleArn` → `arn:aws:iam::994254241749:role/github-gutibm-dispatch`
- `CollectorRoleArn` → `arn:aws:iam::994254241749:role/github-gutibm-collector`

To change anything later, rerun the same command — CloudFormation updates
in place.

## Step 2 — create the GitHub environment

Repo → **Settings → Environments → New environment**, name it exactly
`aws-campaign`. Under *Deployment protection* add yourself as a required
reviewer — every dispatch/collect run then waits for your approve click.

Under *Environment variables*, add:

| Name | Value |
|---|---|
| `AWS_DISPATCH_ROLE_ARN` | `arn:aws:iam::994254241749:role/github-gutibm-dispatch` |
| `AWS_COLLECTOR_ROLE_ARN` | `arn:aws:iam::994254241749:role/github-gutibm-collector` |

(ARNs are identifiers, not secrets — variables, not Secrets.)

## Step 3 — run the handshake (no spend)

**Actions → campaign-collect → Run workflow** on `main`. It will pause on
the environment approval; approve it. Expected:

- `aws sts get-caller-identity` shows `assumed-role/github-gutibm-collector`
- the write denial-probe is **denied** (proves the role is correctly
  under-scoped)
- a bounded index of `campaign/` prefixes on both buckets is published as a
  comment on a tracking issue — that issue is the Edison-readable surface

If a step fails with `AccessDenied`, the error names the action — the fix
is a policy edit in `github_oidc_stack.yaml` and a re-deploy, not a
credential change. "Re-run failed jobs" reuses the *old* workflow file —
always start a fresh run after a fix.

## What this does NOT do

- No Edison callback into AWS: Edison consumes the *published index* via
  the repo's GitHub connection (its PAT needs `contents:read` +
  `actions:read+write` to dispatch collector runs itself; every run still
  waits on the `aws-campaign` environment approval) and compact uploads via
  `edison-client`. Raw S3 access for Edison is not possible — third-party
  bearer tokens cannot assume roles.
- No real submits yet: only the collector workflow exists. Wiring an actual
  `submit` behind the dispatch role is a separate, deliberate step.
- Only `campaign/` prefixes are in scope. Older flat campaign prefixes
  (`receptor-selection-v1/`, `stage-d-*`, `practice/`, …) sit outside the
  roles' grant; widen `CampaignPrefix`/add prefixes in the stack if Edison
  should index those too.
- No cleanup rights: deletion/Glacier restoration stays on the
  `PowerUserAccess` admin profile — deliberately.
