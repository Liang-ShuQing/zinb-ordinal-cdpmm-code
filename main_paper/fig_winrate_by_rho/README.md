# `fig:winrate_by_rho`

Model-selection win rates (WAIC/LOOIC) by correlation \(\rho\) (Scenario 1, normal truth, \(N=200\), \(S=128\)).

## Files

- `plot_winrate_by_rho.R` — entry (Monte Carlo + plot)
- `prog1_cdpmm_joint_vs_separate.R` — helpers
- `redraw_winrate_by_rho.R` — optional restyle from saved results (if you keep local RDS)

## Run

```bash
Rscript plot_winrate_by_rho.R ./output
```

Paper figure: `win_rate_waic_looic_by_rho_scen1_normal_n200_nsim128.pdf`.
