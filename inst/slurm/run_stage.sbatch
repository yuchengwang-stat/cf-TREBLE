#!/bin/bash -l
set -euo pipefail
STAGE="$1"; CONFIG="$2"
: "${CFT_HOME:?CFT_HOME not exported -- submit through bin/cftreble}"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

if [[ -z "${R_LIBS_USER:-}" ]]; then
  R_LIBS_USER="$(Rscript -e 'cat(strsplit(Sys.getenv("R_LIBS_USER"), .Platform$path.sep)[[1]][1])')"
fi
export R_LIBS_USER
mkdir -p "$R_LIBS_USER" 2>/dev/null || true

echo "host=$(hostname) stage=$STAGE array=${SLURM_ARRAY_TASK_ID:-none} ntasks=${SLURM_NTASKS_PER_NODE:-?}"
exec Rscript "$CFT_HOME/scripts/stage_${STAGE}.R" --config "$CONFIG"
