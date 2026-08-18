#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cmake_file="${repo_root}/tests/CMakeLists.txt"
workflow_file="${repo_root}/.github/workflows/ci.yml"

mapfile -t gpu_tests < <(
  awk '
    /gutibm_add_test\(/ && /LABELS "[^"]*gpu/ {
      line = $0
      sub(/.*gutibm_add_test\(/, "", line)
      sub(/[^A-Za-z0-9_].*/, "", line)
      print line
    }
    /set_tests_properties\(/ {
      line = $0
      sub(/.*set_tests_properties\(/, "", line)
      sub(/[^A-Za-z0-9_].*/, "", line)
      property_name = line
      in_properties = 1
    }
    in_properties && /LABELS "[^"]*gpu/ {
      print property_name
      in_properties = 0
    }
    in_properties && /^[[:space:]]*\)/ { in_properties = 0 }
  ' "${cmake_file}" | sort -u
)

target_text="$(sed -n '/--target /,/^[[:space:]]*$/p' "${workflow_file}")"
missing=0
for test_name in "${gpu_tests[@]}"; do
  target_name="test_${test_name}"
  if ! grep -Eq "(^|[[:space:]])${target_name}([[:space:]]|$)" <<<"${target_text}"; then
    printf 'GPU target missing from cuda-compile target list: %s\n' "${target_name}" >&2
    missing=1
  fi
done

if (( missing != 0 )); then
  exit 1
fi

printf 'GPU target manifest: %d labelled tests covered by cuda-compile\n' \
  "${#gpu_tests[@]}"
