#!/usr/bin/env bash
# One-time administrator setup for the GitHub Actions GPU device-test role.
#
# Run with administrative AWS credentials:
#   bash deploy/aws/08_setup_github_oidc.sh
#
# Store the printed role ARN as repository variable
# GUTIBM_AWS_GPU_DEVICE_ROLE_ARN.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
POLICIES="${HERE}/policies"
AWS_REGION="${AWS_REGION:-us-east-1}"
GITHUB_REPOSITORY="${GITHUB_REPOSITORY:-bckirkup/GutModelBacteriocins}"
ROLE_NAME="${GITHUB_OIDC_ROLE_NAME:-gutibm-github-gpu-device}"
POLICY_NAME="${GITHUB_OIDC_POLICY_NAME:-gutibm-gpu-device}"
THUMBPRINT="${GITHUB_OIDC_THUMBPRINT:-6938fd4d98bab03faadb97b34396831e3780aea1}"

ACCOUNT="$(aws sts get-caller-identity --query Account --output text)"
PROVIDER_ARN="arn:aws:iam::${ACCOUNT}:oidc-provider/token.actions.githubusercontent.com"
if aws iam get-open-id-connect-provider \
  --open-id-connect-provider-arn "${PROVIDER_ARN}" >/dev/null 2>&1; then
  echo "OIDC provider exists: ${PROVIDER_ARN}"
else
  echo "Creating OIDC provider: ${PROVIDER_ARN}"
  aws iam create-open-id-connect-provider \
    --url https://token.actions.githubusercontent.com \
    --client-id-list sts.amazonaws.com \
    --thumbprint-list "${THUMBPRINT}" >/dev/null
fi

TRUST="$(mktemp)"
PERMISSIONS="$(mktemp)"
trap 'rm -f "${TRUST}" "${PERMISSIONS}"' EXIT
sed \
  -e "s|__ACCOUNT_ID__|${ACCOUNT}|g" \
  -e "s|__GITHUB_REPOSITORY__|${GITHUB_REPOSITORY}|g" \
  "${POLICIES}/github-oidc-trust.json" > "${TRUST}"
sed \
  -e "s|__ACCOUNT_ID__|${ACCOUNT}|g" \
  -e "s|__REGION__|${AWS_REGION}|g" \
  "${POLICIES}/github-oidc-permissions.json" > "${PERMISSIONS}"

if aws iam get-role --role-name "${ROLE_NAME}" >/dev/null 2>&1; then
  echo "Updating role trust policy: ${ROLE_NAME}"
  aws iam update-assume-role-policy \
    --role-name "${ROLE_NAME}" \
    --policy-document "file://${TRUST}"
else
  echo "Creating role: ${ROLE_NAME}"
  aws iam create-role \
    --role-name "${ROLE_NAME}" \
    --assume-role-policy-document "file://${TRUST}" \
    --description "GitHub OIDC role for GutIBM GPU device CI" >/dev/null
fi

aws iam put-role-policy \
  --role-name "${ROLE_NAME}" \
  --policy-name "${POLICY_NAME}" \
  --policy-document "file://${PERMISSIONS}"

echo
echo "GitHub repository variable:"
echo "GUTIBM_AWS_GPU_DEVICE_ROLE_ARN=arn:aws:iam::${ACCOUNT}:role/${ROLE_NAME}"
echo "Repository: ${GITHUB_REPOSITORY}"
