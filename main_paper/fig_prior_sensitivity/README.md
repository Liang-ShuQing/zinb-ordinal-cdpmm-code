# `fig:prior_sensitivity`

Sensitivity of ARMSE and ACP to fixed-effect prior variance \(\kappa\) (Scenario 2, mixture truth, \(N=200\), \(S=200\)).

## Files

- `plot_prior_sensitivity.R` — entry
- `prog1_cdpmm_joint_vs_separate.R` — helpers
- `combine_prior_sensitivity_plot.R` — combine panels

## Run

```bash
Rscript plot_prior_sensitivity.R ./output
Rscript combine_prior_sensitivity_plot.R
```

Paper figure: `prior_sensitivity_scen2_mixture_n200_combined.pdf`.
