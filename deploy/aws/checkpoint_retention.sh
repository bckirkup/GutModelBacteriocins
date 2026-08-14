#!/usr/bin/env bash

checkpoint_keys_from_list_response() {
  jq -r '.Contents[]?.Key'
  return 0
}

checkpoint_step_from_key() {
  local key="${1##*/}"
  if [[ "${key}" =~ ^step_([0-9]+)\.h5$ ]]; then
    printf '%d\n' "$((10#${BASH_REMATCH[1]}))"
  elif [[ "${key}" =~ ^step_([0-9]+)$ ]]; then
    printf '%d\n' "$((10#${BASH_REMATCH[1]}))"
  elif [[ "${key}" =~ ^[0-9]+$ ]]; then
    printf '%d\n' "$((10#${key}))"
  fi
  return 0
}

checkpoint_select_retained_keys() {
  local latest_step="$1"
  local newest_k="$2"
  local archive_interval="$3"
  local key step latest_step_number newest_steps
  local -a keys

  mapfile -t keys
  latest_step_number="$(checkpoint_step_from_key "${latest_step}")"
  [[ -n "${latest_step_number}" ]] || return 0
  newest_steps="$(
    for key in "${keys[@]}"; do
      checkpoint_step_from_key "${key}"
    done | sort -nr | head -n "${newest_k}"
  )"

  for key in "${keys[@]}"; do
    step="$(checkpoint_step_from_key "${key}")"
    [[ -n "${step}" ]] || continue
    if [[ "${step}" == "${latest_step_number}" ]] \
      || ((10#${step} % archive_interval == 0)) \
      || grep -Fxq "${step}" <<<"${newest_steps}"; then
      printf '%s\n' "${key}"
    fi
  done
  return 0
}
