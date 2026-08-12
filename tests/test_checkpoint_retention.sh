#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
FIXTURE="$ROOT/tests/fixtures/checkpoint_retention.json"
# shellcheck source=../deploy/aws/checkpoint_retention.sh
source "$ROOT/deploy/aws/checkpoint_retention.sh"

fail() {
  echo "test_checkpoint_retention: FAILED: $*" >&2
  exit 1
}

KEYS="$(jq -r '.keys[]' "$FIXTURE")"
LATEST_STEP="$(jq -r '.latest.uri | split("/")[-1]' "$FIXTURE" |
  sed -E 's/^step_([0-9]+)\.h5$/\1/')"
RETAINED="$(
  checkpoint_select_retained_keys "${LATEST_STEP}" 2 360 <<<"${KEYS}"
)"

grep -Fxq 'campaign/calibration/ckpt/step_000810.h5' <<<"${RETAINED}" ||
  fail "latest checkpoint was not retained"
grep -Fxq 'campaign/calibration/ckpt/step_000780.h5' <<<"${RETAINED}" ||
  fail "newest-K checkpoint was not retained"
grep -Fxq 'campaign/calibration/ckpt/step_000360.h5' <<<"${RETAINED}" ||
  fail "360-step archive was not retained"
grep -Fxq 'campaign/calibration/ckpt/step_000720.h5' <<<"${RETAINED}" ||
  fail "720-step archive was not retained"
if grep -Fxq 'campaign/calibration/ckpt/step_000750.h5' <<<"${RETAINED}"; then
  fail "non-retained checkpoint was unexpectedly retained"
fi

LATEST_TARGET='campaign/calibration/ckpt/step_000720.h5'
RETAINED_WITH_OLDER_LATEST="$(
  checkpoint_select_retained_keys 720 2 360 <<<"${KEYS}"
)"
grep -Fxq "${LATEST_TARGET}" <<<"${RETAINED_WITH_OLDER_LATEST}" ||
  fail "latest.json target was not retained"

echo "Checkpoint retention selection self-test passed."
