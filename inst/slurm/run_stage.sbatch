#!/bin/bash -l
# Generic stage runner.  All resource flags come from `cftreble submit`, which
# reads them out of the config -- nothing is pinned in this file.
#   $1 = stage name    $2 = config path
# CFT_HOME arrives through the exported environment (sbatch --export=ALL default).
set -euo pipefail
STAGE="$1"; CONFIG="$2"
: "${CFT_HOME:?CFT_HOME not exported -- submit through bin/cftreble}"

# SLURM does not set OMP_NUM_THREADS.  Without it OpenMP/BLAS sees all 128 cores
# of the node while the cgroup only allows the few we asked for.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

# On clusters where $HOME is automounted, R only appends the default user library
# if that directory already exists at startup; if the mount has not been
# triggered yet the path is dropped silently and the job dies on the first
# library() call, which looks like a missing package.  Set it explicitly.
if [[ -z "${R_LIBS_USER:-}" ]]; then
  R_LIBS_USER="$(Rscript -e 'cat(strsplit(Sys.getenv("R_LIBS_USER"), .Platform$path.sep)[[1]][1])')"
fi
export R_LIBS_USER
mkdir -p "$R_LIBS_USER" 2>/dev/null || true

echo "host=$(hostname) stage=$STAGE array=${SLURM_ARRAY_TASK_ID:-none} ntasks=${SLURM_NTASKS_PER_NODE:-?}"
# Stage scripts live in inst/scripts; CFT_HOME is the inst directory.
exec Rscript "$CFT_HOME/scripts/stage_${STAGE}.R" --config "$CONFIG"
