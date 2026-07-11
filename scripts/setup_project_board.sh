#!/usr/bin/env bash
# Create GutIBM GitHub labels, milestones, and a Projects v2 board.
# Requires: gh CLI with `project`, `repo` scopes
# Usage: ./scripts/setup_project_board.sh [--dry-run]

set -euo pipefail

REPO="${REPO:-bckirkup/GutModelBacteriocins}"
OWNER="${REPO%%/*}"
MILESTONE_P6='P6 — HPC & GPU phase 2'
DRY_RUN=false
if [[ "${1:-}" == "--dry-run" ]]; then
  DRY_RUN=true
fi

run() {
  if $DRY_RUN; then
    echo "[dry-run] $*"
  else
    "$@"
  fi
}

echo "=== GutIBM project board setup for ${REPO} ==="

# ── Labels ──────────────────────────────────────────────────────────────
declare -A LABELS=(
  ["track:hdf5"]="1e40af|HDF5 I/O and checkpoint"
  ["track:mpi"]="0e7490|MPI multi-rank correctness"
  ["track:config"]="fb8500|Input parser and configuration"
  ["track:docs"]="6f42c1|Documentation and agent guides"
  ["track:scale"]="8250df|Scaling benchmarks and profiling"
  ["track:science"]="bf3989|Experimental validation and algorithms"
  ["track:ci"]="1d76db|CI coverage and parity gates"
  ["track:gpu"]="5319e7|CUDA GPU acceleration"
  ["board:ready"]="ededed|On project board — ready column"
  ["board:backlog"]="d4d4d4|On project board — backlog column"
)

for name in "${!LABELS[@]}"; do
  IFS='|' read -r color desc <<< "${LABELS[$name]}"
  run gh label create "$name" --repo "$REPO" --color "$color" --description "$desc" 2>/dev/null \
    || echo "  label exists: $name"
done

# ── Milestones ──────────────────────────────────────────────────────────
declare -A MILESTONES=(
  ["${MILESTONE_P6}"]="Post-#152: MPI np4+, GPU CI parity, mechanics/FMM GPU, science regressions.|2026-09-30"
)

for title in "${!MILESTONES[@]}"; do
  IFS='|' read -r desc due <<< "${MILESTONES[$title]}"
  args=(api "repos/${REPO}/milestones" -f title="$title" -f description="$desc")
  if [[ -n "$due" ]]; then
    args+=(-f "due_on=${due}T00:00:00Z")
  fi
  run gh "${args[@]}" 2>/dev/null || echo "  milestone exists: $title"
done

# ── Issue labels + milestones ───────────────────────────────────────────
assign_issue() {
  local num="$1" labels="$2" milestone="$3"
  echo "  issue #${num} → labels: ${labels} milestone: ${milestone}"
  if ! $DRY_RUN; then
  # shellcheck disable=SC2086
    gh issue edit "$num" --repo "$REPO" --add-label "$labels" --milestone "$milestone" 2>/dev/null \
      || gh issue edit "$num" --repo "$REPO" --add-label "$labels" 2>/dev/null \
      || echo "    warning: could not edit #${num}"
  fi
}

echo ""
echo "=== Tagging open issues (Jul 2026 backlog) ==="
assign_issue 154 "track:mpi,board:ready" "${MILESTONE_P6}"
assign_issue 156 "track:mpi,board:ready" "${MILESTONE_P6}"
assign_issue 158 "track:ci,board:ready" "${MILESTONE_P6}"
assign_issue 155 "track:ci,board:ready" "${MILESTONE_P6}"
assign_issue 157 "track:gpu,board:backlog" "${MILESTONE_P6}"
assign_issue 159 "track:science,board:backlog" "${MILESTONE_P6}"
assign_issue 160 "track:science,board:backlog" "${MILESTONE_P6}"
assign_issue 161 "track:ci,board:backlog" "${MILESTONE_P6}"

# ── GitHub Project (v2) ─────────────────────────────────────────────────
echo ""
echo "=== Creating Projects v2 board ==="
if $DRY_RUN; then
  echo "[dry-run] gh project create --owner $OWNER --title 'GutIBM Roadmap' --format board"
else
  if ! gh auth status 2>&1 | grep -q 'project'; then
    echo "WARNING: gh token may lack 'project' scope."
    echo "  Run: gh auth refresh -h github.com -s project"
  fi

  PROJECT_URL=$(gh project create --owner "$OWNER" --title "GutIBM Roadmap" --format board 2>&1) || true
  if [[ -z "$PROJECT_URL" ]]; then
    echo "Could not create project (check 'project' scope). Labels/milestones may still have been created."
    echo "Create the board manually: ${REPO} → Projects → New project → Board"
    exit 0
  fi
  echo "Created: $PROJECT_URL"

  PROJECT_NUM=$(gh project list --owner "$OWNER" --format json --limit 5 \
    | python3 -c "import json,sys; d=json.load(sys.stdin); print(next((p['number'] for p in d['projects'] if p['title']=='GutIBM Roadmap'), ''))")

  if [[ -n "$PROJECT_NUM" ]]; then
  run gh project link "$PROJECT_NUM" --owner "$OWNER" --repo "$REPO"
  for num in 154 155 156 157 158 159 160 161; do
    run gh project item-add "$PROJECT_NUM" --owner "$OWNER" --url "https://github.com/${REPO}/issues/${num}" 2>/dev/null \
      || echo "  could not add issue #${num} to project"
  done
  echo ""
  echo "Open board: gh project view ${PROJECT_NUM} --owner ${OWNER} --web"
  fi
fi

echo ""
echo "Done. See docs/PROJECT_BOARD.md for column layout and PR bundles."
