#!/usr/bin/env bash

checkpoint_step_from_key() {
  local key="${1##*/}"
  if [[ "${key}" =~ ^step_([0-9]+)\.h5$ ]]; then
    printf '%s\n' "${BASH_REMATCH[1]}"
  fi
}

checkpoint_select_retained_keys() {
  local latest_step="$1"
  local newest_k="$2"
  local archive_interval="$3"
  local key step newest_steps
  local -a keys

  mapfile -t keys
  newest_steps="$(
    for key in "${keys[@]}"; do
      checkpoint_step_from_key "${key}"
    done | sort -nr | head -n "${newest_k}"
  )"

  for key in "${keys[@]}"; do
    step="$(checkpoint_step_from_key "${key}")"
    [[ -n "${step}" ]] || continue
    if [[ "${step}" == "${latest_step}" ]] \
      || ((10#${step} % archive_interval == 0)) \
      || grep -Fxq "${step}" <<<"${newest_steps}"; then
      printf '%s\n' "${key}"
    fi
  done
}
