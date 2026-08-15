#!/usr/bin/env bash
set -euo pipefail

tmp="$(mktemp -d)"
trap 'rm -rf "${tmp}"' EXIT
ready="${tmp}/ready"
observed="${tmp}/observed"
group_script="${tmp}/group.sh"
child_script="${tmp}/child.sh"
cat >"${child_script}" <<'SCRIPT'
#!/usr/bin/env bash
set -eu
trap 'printf "%s\n" "$BASHPID" >> "$OBSERVED"; exit 0' TERM
printf "%s\n" "$BASHPID" >> "$OBSERVED"
while :; do sleep 0.01; done
SCRIPT
chmod +x "${child_script}"
cat >"${group_script}" <<'SCRIPT'
#!/usr/bin/env bash
set -u
observed="$1"
ready="$2"
trap ':' TERM
touch "${ready}.leader"
for _ in 1 2 3; do
  OBSERVED="${observed}" bash "$3" &
done
touch "${ready}.children"
wait || true
SCRIPT
chmod +x "${group_script}"

setsid "${group_script}" "${observed}" "${ready}" "${child_script}" 2>/dev/null &
leader=$!
pgid="${leader}"
for _ in {1..100}; do
  kill -0 "${leader}" 2>/dev/null || exit 1
  [[ "$(ps -o pgid= -p "${leader}" | tr -d ' ')" == "${pgid}" ]] && break
  sleep 0.01
done
for _ in {1..200}; do
  [[ -f "${ready}.children" ]] && break
  sleep 0.01
done
if [[ ! -f "${ready}.children" ]]; then
  echo "process-group smoke failed: children readiness marker missing" >&2
  exit 1
fi
for _ in {1..200}; do
  [[ "$(wc -l < "${observed}")" -ge 3 ]] && break
  sleep 0.01
done
if [[ "$(wc -l < "${observed}")" -lt 3 ]]; then
  echo "process-group smoke failed: three child processes never became ready" >&2
  exit 1
fi
kill -TERM -- "-${pgid}"
touch "${ready}.done"
for _ in {1..200}; do
  [[ "$(wc -l < "${observed}")" -ge 6 ]] && break
  sleep 0.01
done
wait "${leader}" || true
observed_lines="$(wc -l < "${observed}")"
if [[ "${observed_lines}" -ne 6 ]]; then
  echo "process-group smoke failed: expected 3 child signal observations, got ${observed_lines}" >&2
  exit 1
fi

# This proves only the shell-level setsid/negative-PGID signalling contract.
# It does not prove OpenMPI process-group behaviour, AWS Batch Spot reclaim
# timing, end-to-end S3 upload ordering, or production rank termination.
