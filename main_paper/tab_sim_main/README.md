# `tab:sim_main`

Main Monte Carlo table (joint CDPMM): Scenario 2, mixture random-effects truth, \(N\in\{200,400\}\).

## Files

- `prog1_cdpmm_joint_vs_separate.R` — full simulation driver (HPC recommended, \(S=500\))
- `build_sim_tables.R` — build LaTeX from finished result folders

Sibling folders `tab_sim_scen1_normal`, `tab_sim_scen2_normal`, `tab_sim_scen1` share the same pipeline with different Scenario / RE settings.

## Run (outline)

```bash
# On HPC: set SCENARIO, RE_DIST, N, N_SIM, CHAIN, … then
Rscript prog1_cdpmm_joint_vs_separate.R /path/to/out

# After all cells finish, point build_sim_tables.R at result roots:
Rscript build_sim_tables.R
```

Paper `\input{sim_tables_generated}` → label `tab:sim_main`.
