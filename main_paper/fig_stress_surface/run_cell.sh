#!/bin/bash
# Run ONE Figure-9 stress cell: Scenario 2, mixture truth, P = HD null covariates.
#
# usage: bash run_cell.sh <N> <HD_EXTRA> [out_root] [n_cores] [n_sim]
#   N        sample size (50 / 100 / 200 / 400)
#   HD_EXTRA number of extra null covariates per submodel (0 / 5 / 20 / 30)
#   out_root where the per-cell folder is created (default: current directory)
#   n_cores  parallel workers (default 64)
#   n_sim    number of replicates (default 500)
#
# Set RSCRIPT if Rscript is not on PATH, e.g.
#   RSCRIPT=~/miniforge3/envs/rsim/bin/Rscript bash run_cell.sh 50 0 out
set -e
N="${1:?usage: run_cell.sh <N> <HD_EXTRA> [out_root] [n_cores] [n_sim]}"
HD="${2:?usage: run_cell.sh <N> <HD_EXTRA> [out_root] [n_cores] [n_sim]}"
ROOT="${3:-$PWD}"
NC="${4:-64}"
NSIM="${5:-500}"
RSCRIPT="${RSCRIPT:-Rscript}"

HERE="$(cd "$(dirname "$0")" && pwd)"
ENGINE="$HERE/../../00_shared/prog1_cdpmm_joint_vs_separate.R"
[ -f "$ENGINE" ] || { echo "engine not found: $ENGINE" >&2; exit 1; }

CELL="n${N}_hd${HD}"
D="$ROOT/$CELL"
mkdir -p "$D"
cd "$D"

export OUT_DIR="$D"
export OMP_NUM_THREADS=1
export SCENARIO=2
export RE_DIST=mixture
export HD_EXTRA="$HD"
export N="$N"
export N_SIM="$NSIM"
export CHAIN=10000
export BURN=5000
export THIN=5
export G_MIX=8
export N_CORES="$NC"

echo "START $(date) N=$N HD_EXTRA=$HD N_SIM=$NSIM cores=$NC" > run.log
"$RSCRIPT" "$ENGINE" "$D" >> run.log 2>&1
echo "DONE $(date)" >> run.log
