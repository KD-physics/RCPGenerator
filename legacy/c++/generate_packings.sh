#!/usr/bin/env bash
set -euo pipefail

# ── CONFIGURATION ────────────────────────────────────────────────────────────
NPROC=$(nproc)               # number of parallel jobs = number of logical cores
# NPROC=2   # 
MASTER=master_list.txt  # adjust path if needed
# SCRIPTDIR=''             # where generate_single_case.sh lives
# ─────────────────────────────────────────────────────────────────────────────

# ── CLEANUP ON Ctrl+C ───────────────────────────────────────────────────────
# When user hits Ctrl+C, kill all child processes (jobs) and exit
trap 'echo "⏹️  Caught SIGINT — killing all jobs…"; kill 0; exit 1' INT
# ─────────────────────────────────────────────────────────────────────────────

echo "▶️  Starting pipeline with up to $NPROC concurrent jobs"
echo "   Master list: $MASTER"
echo

# ── LAUNCH LOOP ─────────────────────────────────────────────────────────────
# Read each job line, launch it, and throttle to $NPROC
while IFS= read -r LINE || [[ -n "$LINE" ]]; do
  # skip empty or commented lines
  [[ -z "$LINE" || "${LINE:0:1}" == "#" ]] && continue

  echo "➤ Launching: $LINE"
  # Fire it off
  bash "generate_single_case.sh" $LINE &

  # If we've hit our limit, wait for at least one to finish
  while (( $(jobs -r | wc -l) >= NPROC )); do
    sleep 1
  done
done < "$MASTER"

# ── WAIT FOR ALL TO FINISH ───────────────────────────────────────────────────
wait
echo
echo "✅ All jobs completed."
