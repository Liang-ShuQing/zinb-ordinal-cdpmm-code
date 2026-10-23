# `fig:hrs_ppc`

Posterior predictive checks for the HRS joint CDPMM fit:
(a) observed HSPNIT proportions with posterior predictive means and intervals,
(b) posterior predictive distribution of \(\mathrm{cor}(\text{HSPNIT},\text{SHLT})\),
(c) mean HSPNIT by SHLT.

## Files

| File | Panels |
|------|--------|
| `plot_hrs_ppc_ab.R` | (a) + (b) from the saved replicate summary |
| `plot_hrs_ppc_mean_y1_by_y2.R` | (c) from the saved posterior samples |
| `plot_hrs_ppc_abc.R` | combines (a)(b)(c) into the paper figure |

## Inputs (not shipped)

Replicate summaries and thinned posterior samples come from the HRS fit
([`../tab_hrs_posterior/code.R`](../tab_hrs_posterior/); HPC bundle `hrs_ppc_g8k23_xh5_20260915`):
`hrs_ppc_reps_summary.rds`, `combined_samples_slim.rds`.
HRS microdata are not shipped; obtain them under their terms (<https://hrs.isr.umich.edu/>).

## Run

```bash
# (a) + (b)
Rscript plot_hrs_ppc_ab.R ./results /path/to/hrs_ppc_reps_summary.rds

# (c) -> results/hrs_ppc_mean_y1_by_y2.rds
Rscript plot_hrs_ppc_mean_y1_by_y2.R ./results /path/to/project \
  /path/to/combined_samples_slim.rds /path/to/Data_model_complete.csv \
  /path/to/hrs_ppc_reps_summary.rds

# combine into the paper figure
Rscript plot_hrs_ppc_abc.R ./results /path/to/project \
  /path/to/hrs_ppc_reps_summary.rds ./results/hrs_ppc_mean_y1_by_y2.rds
```

Paper figure: `hrs_ppc_abc.pdf`.
