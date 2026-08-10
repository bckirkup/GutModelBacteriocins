#!/usr/bin/env bash
# Create the role Devin assumes for GutIBM GPU work, and reduce the
# gutibm-devin user to assume-role-only.
#
# Run this with your own admin credentials, not Devin's.
#
#   GUTIBM_EXTERNAL_ID='<a long random string>' ./06_setup_devin_role.sh
#
# Store the same external id as the GUTIBM_AWS_EXTERNAL_ID secret in Devin,
# mirroring the KINGSANDI_AWS_EXTERNAL_ID setup on TheKingsAndI.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
POLICIES="${HERE}/policies"

USER_NAME="${USER_NAME:-gutibm-devin}"
ROLE_NAME="${ROLE_NAME:-gutibm-devin-role}"
SESSION_HOURS="${SESSION_HOURS:-4}"

if [[ -z "${GUTIBM_EXTERNAL_ID:-}" ]]; then
  echo "GUTIBM_EXTERNAL_ID is required." >&2
  echo "Generate one with: openssl rand -hex 32" >&2
  exit 1
fi

ACCOUNT="$(aws sts get-caller-identity --query Account --output text)"
echo "Account: ${ACCOUNT}"

TRUST="$(mktemp)"
trap 'rm -f "${TRUST}"' EXIT
sed "s|REPLACE_WITH_EXTERNAL_ID|${GUTIBM_EXTERNAL_ID}|" \
  "${POLICIES}/devin-assume-role-trust.json" > "${TRUST}"

if aws iam get-role --role-name "${ROLE_NAME}" >/dev/null 2>&1; then
  echo "Updating trust policy on ${ROLE_NAME}"
  aws iam update-assume-role-policy \
    --role-name "${ROLE_NAME}" \
    --policy-document "file://${TRUST}"
else
  echo "Creating ${ROLE_NAME}"
  aws iam create-role \
    --role-name "${ROLE_NAME}" \
    --assume-role-policy-document "file://${TRUST}" \
    --max-session-duration "$((SESSION_HOURS * 3600))" \
    --description "Scoped GutIBM access assumed by the ${USER_NAME} user"
fi

echo "Attaching permissions to ${ROLE_NAME}"
aws iam put-role-policy \
  --role-name "${ROLE_NAME}" \
  --policy-name gutibm-devin-permissions \
  --policy-document "file://${POLICIES}/devin-role-permissions.json"

echo "Restricting ${USER_NAME} to assume-role-only"
aws iam put-user-policy \
  --user-name "${USER_NAME}" \
  --policy-name gutibm-devin-assume-only \
  --policy-document "file://${POLICIES}/devin-user-assume-only.json"

echo
echo "Inline policies now on ${USER_NAME}:"
aws iam list-user-policies --user-name "${USER_NAME}" --output text
echo
echo "Managed policies still attached to ${USER_NAME} (detach any that grant"
echo "more than assume-role, or the tightening above achieves nothing):"
aws iam list-attached-user-policies --user-name "${USER_NAME}" --output text
echo
echo "Role ARN: arn:aws:iam::${ACCOUNT}:role/${ROLE_NAME}"
echo "Store GUTIBM_AWS_EXTERNAL_ID and GUTIBM_AWS_ROLE_ARN as Devin secrets."
