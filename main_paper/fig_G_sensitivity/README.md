# `fig:G_sensitivity`

Sensitivity of ARMSE and ACP to CDPMM truncation level \(G\) (Scenario 2, \(N=200\), \(S=200\)).

## Files

- `plot_G_sensitivity.R` — Monte Carlo / panel generation
- `prog1_cdpmm_joint_vs_separate.R` — helpers
- `combine_G_sensitivity_plot.R` — combine normal/mixture panels into one PDF

## Run

```bash
Rscript plot_G_sensitivity.R ./output
Rscript combine_G_sensitivity_plot.R   # after CSVs/PDFs exist; adjust paths in script if needed
```

Paper figure: `G_sensitivity_armse_acp_scen2_n200_combined.pdf`.
