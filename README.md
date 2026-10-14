# Adaptive potential across temperature and diet in *Drosophila*

Analysis code for Walter et al. A half-sib breeding design (143 sires, 688 dams, 6906 females) was assayed in six treatments: three temperatures (14, 25, 28 °C) crossed with two diets (control, low-calorie). A zero-inflated negative binomial sire model (`k0`, fitted in Stan) estimates the additive genetic variance in relative fitness, `V_A(w)`, in each treatment and the genetic correlations in fitness among treatments.

## Requirements

- R 4.6.0 with the packages pinned in `renv.lock`. Restore them with `renv::restore()` from the project root.
- CmdStan, needed only to refit the model (`cmdstanr::install_cmdstan()`).

## Running

Run the scripts in order, one at a time, from the project root:

```         
Rscript R/01_data.R          # data checks, design tables, descriptive plots
Rscript R/02_fit.R           # fit k0 (skipped if results/data/k0.rds exists)
Rscript R/03_diagnostics.R   # convergence, posterior predictive checks, LOO
Rscript R/04_results.R       # mean fitness, V_A(w), diet contrasts, r_A
Rscript R/05_figures_ms.R    # manuscript Figs 2-4
```

The fitted model can be found at `results/data/k0.rds` (a new fit takes about 20-25 minutes on four cores).

## Layout

| path | contents |
|----|----|
| `R/00_functions.R` | constants and helpers, sourced by every script |
| `R/01_data.R` … `R/05_figures_ms.R` | the pipeline, run in order |
| `stan/k0_zinb.stan` | the model |
| `data/Data_Productivity.csv` | the data, one row per female |
| `results/data/k0.rds` | the fitted model |
| `results/tables/`, `results/plots/` | created by the scripts |
