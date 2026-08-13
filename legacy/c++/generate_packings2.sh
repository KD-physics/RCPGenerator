#!/usr/bin/env bash
set -euo pipefail

# ── CONFIGURATION ────────────────────────────────────────────────────────────
: "${NPROC:=$(nproc)}"          # number of parallel jobs
MASTER=master_list.txt          # path to your master list
SAMPLE_COUNT=${SAMPLE_COUNT:-0} # 0 = run all, >0 = sample that many
# ─────────────────────────────────────────────────────────────────────────────

# ── CLEANUP ON Ctrl+C ───────────────────────────────────────────────────────
trap 'echo "⏹️  Caught SIGINT — killing all jobs…"; kill 0; exit 1' INT
# ─────────────────────────────────────────────────────────────────────────────

echo "▶️  Starting pipeline with up to $NPROC concurrent jobs"
echo "   Master list: $MASTER"
if (( SAMPLE_COUNT > 0 )); then
  TOTAL=$(grep -v '^[[:space:]]*#' "$MASTER" | grep -v '^[[:space:]]*$' | wc -l)
  echo "   Sampling $SAMPLE_COUNT out of $TOTAL lines"
  LAUNCH_LIST=$(mktemp)
  grep -v '^[[:space:]]*#' "$MASTER" | grep -v '^[[:space:]]*$' | shuf -n "$SAMPLE_COUNT" > "$LAUNCH_LIST"
else
  echo "   Shuffling all lines"
  LAUNCH_LIST=$(mktemp)
  grep -v '^[[:space:]]*#' "$MASTER" | grep -v '^[[:space:]]*$' | shuf > "$LAUNCH_LIST"
fi
echo

# ── LAUNCH LOOP ─────────────────────────────────────────────────────────────
while IFS= read -r LINE || [[ -n "$LINE" ]]; do
  [[ -z "$LINE" || "${LINE:0:1}" == "#" ]] && continue

  echo "➤ Launching: $LINE"
  bash "generate_single_case.sh" $LINE &

  while (( $(jobs -r | wc -l) >= NPROC )); do
    sleep 1
  done
done < "$LAUNCH_LIST"

# ── WAIT FOR ALL TO FINISH ───────────────────────────────────────────────────
wait
echo
echo "✅ All jobs completed."
