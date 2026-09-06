#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mode=planning
execution_sha=""
args=("$@")
for ((i=0; i<${#args[@]}; i++)); do
  case "${args[$i]}" in
    --deployment) mode=deployment ;;
    --execution-source-sha) ((i+=1)); execution_sha="${args[$i]:-}" ;;
    --execution-source-sha=*) execution_sha="${args[$i]#*=}" ;;
  esac
done
python3 prepare.py --clean "$@"
if [[ "$mode" == deployment ]]; then
  [[ -n "$execution_sha" ]] || { echo 'REFUSED: --deployment requires --execution-source-sha' >&2; exit 2; }
  python3 preflight.py --deployment --execution-source-sha "$execution_sha"
else
  python3 preflight.py
fi
python3 analyze.py
