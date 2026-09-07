# Campaign package provenance correction — completion note

**Status:** package correction complete; planning and isolated deployment-path validation executed; no AWS calls, image builds, or simulations executed.

## Corrected defect

The original package overloaded one `source_sha` and required the audited baseline commit to be the runtime build. That made it impossible to commit the campaign package before building it. The corrected contract uses:

- `model_baseline_sha=ad037ef96f3433a1ad4008994956f89c68a3ed39` for the audited design baseline;
- `execution_source_sha` for the later commit containing the campaign package and actually used to build the image.

Planning artifacts carry an explicit execution-SHA placeholder. Deployment artifacts cannot be generated in this tarball. In a canonical Git checkout, they require an explicit execution SHA equal to `HEAD`, baseline ancestry, committed/unchanged provenance-bearing package files, and an immutable image digest. Inputs and manifests record both SHAs and the digest. Returned HDF5 must report the exact execution SHA.

## Required next operational steps

Copy package → commit package → capture HEAD → build/push image from HEAD → resolve digest → regenerate deployment artifacts (allowed to remain uncommitted) → deployment preflight → submit Q only. No old-baseline build and no uncommitted-plan workaround is required or permitted.

## Validation scope

The package test suite, planning generation, planning preflight, and analyzer missing-output path passed locally. Deployment correctly fails in this no-Git retrieved snapshot. I also exercised the full two-commit workflow in an isolated synthetic Git history: deployment generation passed, deployment preflight checked 62 runs and passed, and Q command generation printed (but did not execute) a submission command. A positive test against the exact canonical Git history was not possible because the canonical repository fetch required unavailable credentials; this does not weaken the local tarball baseline check, but the real operator must run deployment preflight in the authenticated canonical checkout. No cloud submission, image build, GPU model run, or scientific promotion is part of this correction.

## Revised C-to-D addendum

`campaign_decision_record.json` supersedes the prior implied handoff. A/B passed under `b884cc54b1c0684c079a90195c031e388fe70534`; ecological C passed under `3f176b26c0d18a22a61db218e56106b1355b781e`; the PR416 intrinsic-transport gate is invalid/failed and superseded. The revised assay is pending. D now has 15 jobs, including three formal same-revision plasmid-free nulls, and remains blocked until assay pass. The prior validation statement covering 62 runs is historical; this patch revalidates 65 planning runs. No AWS action is part of this addendum.
