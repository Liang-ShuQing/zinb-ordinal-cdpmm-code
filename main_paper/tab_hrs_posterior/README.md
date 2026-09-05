# `tab:hrs_posterior`

Posterior summaries for the HRS joint model (\(3\times 30{,}000\) iterations, burn-in \(15{,}000\), thin \(5\), \(G=8\)).

## Files

- `prog1_cdpmm_joint_vs_separate.R` — primary joint fit
- `code.R` — alternate/legacy entry used in some runs
- `plot_diagnostics.R` — optional diagnostics

## Data (not shipped)

Same HRS extract requirement as [`fig_hrs_trace`](../fig_hrs_trace/).

## Run

```bash
Rscript prog1_cdpmm_joint_vs_separate.R ./output
# Then export posterior means / 95% CrI from the saved RDS into the manuscript table.
```

Cross-ref: [`fig_hrs_trace`](../fig_hrs_trace/).
