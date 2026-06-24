#!/usr/bin/env bash
# Validate GutIBM example and parser fixture configs are strict JSON (#54).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PYTHON="${PYTHON:-python3}"

mapfile -t FILES < <(
  find "$ROOT/examples" "$ROOT/tests/fixtures" -name '*.json' -type f | sort
)

if [[ ${#FILES[@]} -eq 0 ]]; then
  echo "ERROR: no JSON config files found" >&2
  exit 1
fi

echo "Validating ${#FILES[@]} config JSON file(s)..."
for file in "${FILES[@]}"; do
  "$PYTHON" -m json.tool "$file" > /dev/null
  echo "  OK $(realpath --relative-to="$ROOT" "$file")"
done

echo "All config JSON files are valid."
