#!/usr/bin/env bash
# Install and preview the gutibm ECR retention policy.
#
# Run this with administrator credentials in CloudShell. The Devin role cannot
# put an ECR lifecycle policy until 06_setup_devin_role.sh is rerun with the
# widened permissions.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
POLICIES="${HERE}/policies"
REPO_NAME="${REPO_NAME:-gutibm}"
POLICY_FILE="${POLICIES}/ecr-lifecycle-gutibm.json"

aws ecr put-lifecycle-policy \
  --repository-name "${REPO_NAME}" \
  --lifecycle-policy-text "file://${POLICY_FILE}" \
  >/dev/null

PREVIEW_STATUS="$(
  aws ecr get-lifecycle-policy-preview \
    --repository-name "${REPO_NAME}" \
    --query 'status' \
    --output text 2>/dev/null || true
)"
if [[ "${PREVIEW_STATUS}" != "IN_PROGRESS" ]]; then
  aws ecr start-lifecycle-policy-preview \
    --repository-name "${REPO_NAME}" \
    >/dev/null
fi

for ((attempt = 0; attempt < 60; attempt++)); do
  PREVIEW_STATUS="$(
    aws ecr get-lifecycle-policy-preview \
      --repository-name "${REPO_NAME}" \
      --query 'status' \
      --output text 2>/dev/null || true
  )"
  case "${PREVIEW_STATUS}" in
    COMPLETE)
      break
      ;;
    FAILED)
      echo "ECR lifecycle policy preview failed." >&2
      exit 1
      ;;
    IN_PROGRESS|"")
      sleep 2
      ;;
    *)
      echo "Unexpected ECR lifecycle policy preview status: ${PREVIEW_STATUS}" >&2
      exit 1
      ;;
  esac
done

if [[ "${PREVIEW_STATUS}" != "COMPLETE" ]]; then
  echo "Timed out waiting for the ECR lifecycle policy preview." >&2
  exit 1
fi

EXPIRE_COUNT="$(
  aws ecr get-lifecycle-policy-preview \
    --repository-name "${REPO_NAME}" \
    --query 'length(previewResults[?action==`EXPIRE`])' \
    --output text
)"
echo "ECR lifecycle policy installed for ${REPO_NAME}."
echo "Lifecycle preview: ${EXPIRE_COUNT} image(s) would expire."
aws ecr get-lifecycle-policy-preview \
  --repository-name "${REPO_NAME}" \
  --query 'previewResults[?action==`EXPIRE`].[imageDigest, imageTags, appliedRulePriority]' \
  --output table
echo "Next: rerun bash deploy/aws/06_setup_devin_role.sh so Devin can prune future leftovers."
