#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
ENTRYPOINT="${ROOT}/deploy/aws/entry.sh"
DOCKERFILE="${ROOT}/deploy/aws/Dockerfile"

mapfile -t sourced_files < <(
  sed -nE 's/^[[:space:]]*source[[:space:]]+"\$\{SCRIPT_DIR\}\/([^"]+)".*/\1/p' \
    "${ENTRYPOINT}"
)

if [[ "${#sourced_files[@]}" -eq 0 ]]; then
  echo "No SCRIPT_DIR source files found in entrypoint" >&2
  exit 1
fi

for file in "${sourced_files[@]}"; do
  source_path="${ROOT}/deploy/aws/${file}"
  destination="/opt/gutibm/${file}"

  [[ -f "${source_path}" ]] || {
    echo "Missing entrypoint source file: ${source_path}" >&2
    exit 1
  }
  copy_line="$(grep -E \
    "COPY .*deploy/aws/${file}[[:space:]]+${destination}([[:space:]]|$)" \
    "${DOCKERFILE}")"
  grep -Fq -- "--chown=root:root" <<<"${copy_line}" \
    || {
    echo "Dockerfile does not preserve root ownership for sourced file: ${file}" >&2
    exit 1
  }
  grep -Fq -- "--chmod=755" <<<"${copy_line}" || {
    echo "Dockerfile does not preserve executable permissions for sourced file: ${file}" >&2
    exit 1
  }
  [[ -n "${copy_line}" ]] || {
    echo "Dockerfile does not install sourced file: ${file}" >&2
    exit 1
  }
  sed_line="$(grep -E "^[[:space:]]*RUN sed -i 's/\\\\r\\$//'" "${DOCKERFILE}")"
  grep -Fq "${destination}" <<<"${sed_line}" || {
    echo "Dockerfile does not normalize sourced file line endings: ${file}" >&2
    exit 1
  }
done

echo "AWS entrypoint runtime source packaging invariant passed."
