#!/bin/bash
# Queue the 14 stress cells of Figure 9 at S = 500, 3 cells in parallel.
# (The two P = 0 cells at N = 200/400 are NOT run here: they come from the
#  main-text design, which uses the same specification at S = 500.)
#
# usage: bash queue.sh [out_root] [n_cores_per_cell] [n_sim]
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="${1:-$HERE/out}"
NC="${2:-64}"
NSIM="${3:-500}"
MAXPAR=3

mkdir -p "$ROOT"

CELLS="50:0 50:5 50:20 50:30 \
100:0 100:5 100:20 100:30 \
200:5 200:20 200:30 \
400:5 400:20 400:30"

echo "QUEUE START $(date)  cells=$(echo $CELLS | wc -w)  maxpar=$MAXPAR  n_sim=$NSIM"
for c in $CELLS; do
  N="${c%%:*}"
  HD="${c##*:}"
  while [ "$(jobs -rp | wc -l)" -ge "$MAXPAR" ]; do
    sleep 30
  done
  echo "launch N=$N HD=$HD $(date)"
  bash "$HERE/run_cell.sh" "$N" "$HD" "$ROOT" "$NC" "$NSIM" &
  sleep 10
done
wait
echo "QUEUE DONE $(date)"
