#!/usr/bin/env bash
set -euo pipefail
# set -x
# trap 'echo "❌ Error at line $LINENO"; exit 1' ERR

# ─────────────────────────────────────────────────────────────────────────────
# Quick usage‐check:
if [ $# -lt 4 ]; then
  cat <<-USAGE
Usage: $0 ID NDIMS N_PARTICLES DIST [DIST_PARAMS...]
Example: $0 1 3 25000 powerlaw 1.0 20.0 -3.0
USAGE
  exit 1
fi

# ------------------------------------------------------------------
# Usage:
#   ./generate_single_case.sh ID NDIMS N_PARTICLES DIST [DIST_PARAMS...]
#
# Example:
#   ./generate_single_case.sh 1 3 25000 powerlaw 1.0 20.0 -3.0
# ------------------------------------------------------------------

NEIGH=(0    0    1000   1500   3000   3500   8500)

# 1. Read fixed fields
ID="$1";        shift
NDIMS="$1";     shift
NPART="$1";     shift
DIST="$1";      shift

# 2. All remaining args are distribution parameters, in order
DIST_ARGS=("$@")

# 3. Prepare output directory
OUTDIR="/mnt/hda/data/ID_${ID}"
mkdir -p "${OUTDIR}"

# 4. Prepare the provenance log
JOBLOG="${OUTDIR}/job.txt"
echo "=== JOB: ID=${ID} NDIMS=${NDIMS} N=${NPART} DIST=${DIST} ARGS=${DIST_ARGS[*]} ===" \
    > "${JOBLOG}"
echo "" >> "${JOBLOG}"

# 5. Build the common InitializeParticles command
INIT_CMD=( "./InitializeParticles.exe" 
           "--N" "${NPART}"
           "--Ndim" "${NDIMS}"
           "--phi" "0.005"   # fixed phi—change here if you want it driven by master_list
           "--dist" "${DIST}"
         )


# 6. Append distribution-specific flags
case "${DIST}" in
  mono)
    # DIST_ARGS = [ d ]
    INIT_CMD+=( "--d" "${DIST_ARGS[0]}" )
    ;;
  bidisperse)
    # [ d1 d2 p ]
    INIT_CMD+=( "--d1" "${DIST_ARGS[0]}" "--d2" "${DIST_ARGS[1]}" "--p" "${DIST_ARGS[2]}" )
    ;;
  gaussian)
    # [ mu sigma ]
    INIT_CMD+=( "--mu" "${DIST_ARGS[0]}" "--sigma" "${DIST_ARGS[1]}" )
    ;;
  bigaussian)
    # [ mu1 sigma1 mu2 sigma2 p ]
    INIT_CMD+=( "--mu1" "${DIST_ARGS[0]}" "--sigma1" "${DIST_ARGS[1]}" \
                "--mu2" "${DIST_ARGS[2]}" "--sigma2" "${DIST_ARGS[3]}" \
                "--p" "${DIST_ARGS[4]}" )
    ;;
  lognormal)
    # [ mu sigma ]
    INIT_CMD+=( "--mu" "${DIST_ARGS[0]}" "--sigma" "${DIST_ARGS[1]}" )
    ;;
  flat)
    # [ d_min d_max ]
    INIT_CMD+=( "--d_min" "${DIST_ARGS[0]}" "--d_max" "${DIST_ARGS[1]}" )
    ;;
  powerlaw)
    # [ d_min d_max exponent ]
    INIT_CMD+=( "--d_min" "${DIST_ARGS[0]}" "--d_max" "${DIST_ARGS[1]}" \
                "--exponent" "${DIST_ARGS[2]}" )
    ;;
  exponential)
    # [ d_min d_max ]
    INIT_CMD+=( "--d_min" "${DIST_ARGS[0]}" "--d_max" "${DIST_ARGS[1]}" )
    ;;
  *)
    echo "ERROR: Unknown distribution '${DIST}'" >&2
    exit 1
    ;;
esac


# 7. Add box & walls flags (1s and 0s repeated NDIMS times)
# BOX_FLAGS=$(yes 1 | head -n "${NDIMS}" | paste -sd, )
# WALL_FLAGS=$(yes 0 | head -n "${NDIMS}" | paste -sd, )
# INIT_CMD+=( "--box" "${BOX_FLAGS}" "--walls" "${WALL_FLAGS}" )
# 7. Add box & walls flags
BOX_FLAGS=$(printf '1%.0s,' $(seq 1 $NDIMS) | sed 's/,$//')
WALL_FLAGS=$(printf '0%.0s,' $(seq 1 $NDIMS) | sed 's/,$//')
INIT_CMD+=( "--box" "${BOX_FLAGS}" "--walls" "${WALL_FLAGS}" )

# 8. Decide whether to start from init.txt or last save point
if [[ -f "${OUTDIR}/final.txt" ]]; then
  # fully done
  echo "Already have final.txt; skipping." >> "${JOBLOG}"
  exit 0
fi

# find any intermediate 'final_<ITER>.txt' files
INTERMED=( "${OUTDIR}"/final_[0-9]*.txt )
if [[ -e "${INTERMED[0]}" ]]; then
  # pick the one with the highest iteration number
  LAST_SAVE=$(printf "%s\n" "${INTERMED[@]}" \
              | sed -E 's/.*final_([0-9]+)\.txt/\1 \0/' \
              | sort -nr \
              | head -n1 \
              | cut -d' ' -f2-)
  echo "Resuming from intermediate: ${LAST_SAVE}" >> "${JOBLOG}"
  INIT_FILE="${LAST_SAVE}"
else
  # no intermediate saves → fresh initialization
  INIT_FILE="${OUTDIR}/init.txt"
  echo "INIT command:"    >> "${JOBLOG}"
  printf "  %s " "${INIT_CMD[@]}" >> "${JOBLOG}"
  printf "> %s\n" "${INIT_FILE}" >> "${JOBLOG}"

  # run initialization
  "${INIT_CMD[@]}" > "${INIT_FILE}"
fi

# 9. Build & run RCPGenerator command
NEIGHVAL="${NEIGH[$((NDIMS))]}"

# RCP_CMD=( "./RCPGenerator.exe"
#           "--file" "${INIT_FILE}"
#           "--output" "${OUTDIR}/final.txt"
#           "--verbose"
#           "--save-interval" "1000"
#           "--box" "${BOX_FLAGS}"
#           "--walls" "${WALL_FLAGS}"
#         )
RCP_CMD=( "./RCPGenerator.exe"
          "--file" "${INIT_FILE}"
          "--output" "${OUTDIR}/final.txt"
          "--verbose"
          "--box" "${BOX_FLAGS}"
          "--walls" "${WALL_FLAGS}"
          "--NeighborMax" "${NEIGHVAL}"
        )        

echo "" >> "${JOBLOG}"
echo "RCPGenerator command:" >> "${JOBLOG}"
printf "  %s " "${RCP_CMD[@]}" >> "${JOBLOG}"
printf "> %s\n" "${OUTDIR}/summary.txt" >> "${JOBLOG}"

"${RCP_CMD[@]}" > "${OUTDIR}/summary.txt"

# ─── Print “Done” when this case finishes ─────────────────────────────────
# Align columns to match the “Launching” lines
printf "➤ Done:      %s\t%s\t%s\t%s\t%s\n" \
  "$ID" "$NDIMS" "$NPART" "$DIST" "${DIST_ARGS[*]}"

# instead of:
#   ./RCPGenerator.exe … > "${OUTDIR}/summary.txt"
# do:

