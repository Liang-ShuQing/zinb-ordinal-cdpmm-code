# `fig:ppc_combined`

Posterior predictive checks under Scenario 2 with mixture random-effects truth (\(N=200\); MCMC \(5{,}000\) / burn \(2{,}000\) / thin \(5\)).

## Files

- `ppc_plots.R` — entry
- `prog1_cdpmm_joint_vs_separate.R` — model helpers

## Run

```bash
Rscript ppc_plots.R ./output
```

Typical setting in script/env: Scenario 2, mixture RE, \(N=200\). Output: combined PPC PDF used as `ppc_fig4_5_6_combined.pdf` in the paper.
