#!/bin/bash
#SBATCH -J st2400nor
#SBATCH -p wzhctest
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH -t 3-00:00:00

cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
OUT_DIR="$(pwd)"
export OUT_DIR
export OMP_NUM_THREADS=1
export CMDSTAN="${HOME}/cmdstan"
export TBB_CXX_TYPE=gcc
if [ -f /opt/rh/devtoolset-7/enable ]; then
  source /opt/rh/devtoolset-7/enable
fi
export CXX="${CXX:-/opt/rh/devtoolset-7/root/usr/bin/g++}"
export CC="${CC:-/opt/rh/devtoolset-7/root/usr/bin/gcc}"
export PATH="/opt/rh/devtoolset-7/root/usr/bin:${PATH}"
TBB_DIR=$(find ${HOME}/cmdstan -type d -path '*/stan/lib/stan_math/lib/tbb' 2>/dev/null | head -1)
if [ -n "$TBB_DIR" ]; then
  export LIBRARY_PATH="${TBB_DIR}${LIBRARY_PATH:+:${LIBRARY_PATH}}"
  export LD_LIBRARY_PATH="${TBB_DIR}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
fi

export NSIM=500
export N_SUBJ=400
export N=400
export G_MIX=8
export SCENARIO=2
export RE_DIST=normal
export FIT_SEPARATE=0
export STAN_WARMUP=500
export STAN_SAMPLING=500
export STAN_CHAINS=2
export STAN_PARALLEL_CHAINS=1
export STAN_REFRESH=0

/public/software/VersionHub/R/4.2.0/Rscript code.R "${OUT_DIR}"
echo EXIT=$?
ls -la "${OUT_DIR}"/*.csv 2>/dev/null | head