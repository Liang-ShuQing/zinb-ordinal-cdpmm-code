# Figure 9 — `fig:stress_surface`

N × P design-stress surface: three heatmaps of (a) ARMSE, (b) ACP and (c) the largest
absolute bias among null covariates, over sample size `N ∈ {50, 100, 200, 400}` and
`P ∈ {0, 5, 20, 30}` additional null covariates per submodel.

Data-generating design: Scenario 2, mixture random-effects truth (Truth 2), `S = 500`
replicates, `chain = 10000`, `burn = 5000`, `thin = 5`, `G = 8`.

## Files

| File | Purpose |
|------|---------|
| `run_cell.sh` | Run **one** cell at a given `N` and `HD_EXTRA` (default `S = 500`) |
| `queue.sh` | Queue all 14 stress cells, 3 in parallel (3 × 64 = 192 cores) |
| `aggregate_stress_cells.R` | Turn per-cell `joint_param_summary_*.csv` into the 16-row surface table |
| `plot_stress_surface_NP.R` | Draw the 3-panel heatmap from `stress_surface_long.csv` |

The sampler itself is the shared engine
[`00_shared/prog1_cdpmm_joint_vs_separate.R`](../../00_shared/prog1_cdpmm_joint_vs_separate.R),
which supports `HD_EXTRA` to append null (true coefficient 0) standard-normal
covariates to each submodel design. The columns are appended **after** `y1`/`y2` are
generated, so the base data are identical across `P` levels (common random numbers).

## How to reproduce

```bash
# 1. run the cells (each writes its own results in <root>/n<N>_hd<HD>/)
bash run_cell.sh 50 0   "$PWD/out"
bash run_cell.sh 50 5   "$PWD/out"
# ... (see queue.sh for the full list)

# 2. aggregate 14 cells -> 14-row table
Rscript aggregate_stress_cells.R "$PWD/out"

# 3. add the two P = 0 cells at N = 200/400 from the main-text design
#    (same spec, S = 500) and draw the figure
Rscript plot_stress_surface_NP.R "$PWD/out" "$PWD/out"
```

On a cluster replace step 1 by `bash queue.sh`; the scripts honour `RSCRIPT` if
`Rscript` is not on `PATH`.

## Aggregation rule

For every cell the three plotted quantities are computed from
`joint_param_summary_scen2_mixture_joint_vs_sep_n<N>_chain10000[_hd<HD>].csv`:

* **ARMSE** = mean of `RMSE` over the 19 base parameters
  `alpha1–alpha4`, `beta1–beta5`, `gamma1–gamma3`, `r`, `Sigma11`, `Sigma22`,
  `Sigma33`, `rho1–rho3`
* **ACP** = mean of `CP` over the same 19 parameters
* **max |bias| of null coefs** = `max(abs(Bias))` over the appended `P` null
  covariates (`alpha5+`, `beta6+`, `gamma4+`); equal to 0 when `P = 0`

`delta1–delta4` (the category-threshold nuisance functions) and the off-diagonal
entries of `Sigma` are excluded. The same rule reproduces the published values of
all cells.

## Cells

The two `P = 0` cells at `N = 200` and `N = 400` are **not** re-run here: they are
taken from the main-text design (`N_SIM = 500`, identical specification) and only
re-aggregated with the rule above. The other 14 cells are run by `queue.sh`.

## Results (S = 500)

ARMSE / ACP:

| N \ P | +0 | +5 | +20 | +30 |
|-------|----|----|-----|-----|
| 50  | 0.2143 / 0.9572 | 0.2551 / 0.9502 | 0.6089 / 0.8009 | 0.6662 / 0.7720 |
| 100 | 0.1532 / 0.9517 | 0.1619 / 0.9488 | 0.2527 / 0.9164 | 0.4499 / 0.8175 |
| 200 | 0.1125 / 0.9557 | 0.1146 / 0.9566 | 0.1331 / 0.9438 | 0.1543 / 0.9218 |
| 400 | 0.0810 / 0.9529 | 0.0817 / 0.9534 | run by `queue.sh` | run by `queue.sh` |

Largest absolute bias among null coefficients:

| N \ P | +0 | +5 | +20 | +30 |
|-------|----|----|-----|-----|
| 50  | 0 | 0.0108 | 0.0487 | 0.0416 |
| 100 | 0 | 0.0070 | 0.0152 | 0.0205 |
| 200 | 0 | 0.0048 | 0.0082 | 0.0113 |
| 400 | 0 | 0.0048 | run by `queue.sh` | run by `queue.sh` |

## Notes

* `prog1_cdpmm_joint_vs_separate.R` is self-contained (no `source()`): it only needs
  `BayesLogit`, `mvtnorm`, `MCMCpack`, `truncnorm`, `parallel`, `doParallel`,
  `foreach`, `loo`, `coda`.
* On Linux use fork-style `registerDoParallel(cores = n)` (the engine does); avoid
  `makeCluster()` at high core counts.
* 500 replicates × 16 cells is a large job: the 14 re-runs above took ≈ 7,800 core-hours
  on a 192-core node.
