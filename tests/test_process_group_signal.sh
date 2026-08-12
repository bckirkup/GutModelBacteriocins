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
while [[ ! -f "${ready}.children" ]]; do sleep 0.01; done
while [[ "$(wc -l < "${observed}")" -lt 3 ]]; do sleep 0.01; done
kill -TERM -- "-${pgid}"
touch "${ready}.done"
wait "${leader}" || true
test "$(wc -l < "${observed}")" -eq 6

# This proves only the shell-level setsid/negative-PGID signalling contract.
# It does not prove OpenMPI process-group behaviour, AWS Batch Spot reclaim
# timing, end-to-end S3 upload ordering, or production rank termination.
