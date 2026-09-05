# `fig:hrs_trace`

HRS joint-analysis trace plots and posterior densities for selected fixed effects and \(r\) (\(G=8\), sparse cluster initialization; three chains).

## Files

- `prog1_cdpmm_joint_vs_separate.R` — HRS joint fit driver
- `plot_diagnostics.R` — diagnostics / trace PDF from saved RDS

## Data (not shipped)

Place your HRS analysis-ready CSV where the driver expects it (in the original project: `Data_model_complete.csv`). Do **not** commit restricted HRS files to a public repo.

## Run (outline)

```bash
# Fit (long; HPC recommended: 3 x 30000, burn 15000, thin 5, G=8)
Rscript prog1_cdpmm_joint_vs_separate.R ./output

# Plots from RDS
Rscript plot_diagnostics.R ./output ./figures
```

Paper figure: `hrs_g8k23_trace.pdf`.
