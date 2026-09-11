# Multiple imputation and debiased machine learning under MAR and MNAR

A simulation study comparing how standard missing-data methods and a cross-fitted debiased machine learning estimator behave when the missingness mechanism and the model specification vary. Everything runs from seed, and the report regenerates from the saved results.

**Status:** code complete, results not yet generated. Run `make quick` first, then `make all`.

## Questions

1. When an outcome is missing at random and the conditional mean is nonlinear, what does a misspecified imputation model cost in bias and in interval coverage?
2. Does a cross-fitted debiased ML estimator recover what a correctly specified parametric analysis recovers, without the analyst knowing the functional form?
3. When missingness depends on the unobserved values, how far from MAR do the data have to be before the conclusion changes?

Coverage of nominal 95 per cent intervals is the headline measure. Bias alone rewards methods that understate uncertainty, and mean imputation shows exactly that pattern.

## Design

### Scenario A: missing outcome, target E[Y]

The outcome is generated as

```
Y = 1 + X1 + 0.5*X2 - 0.75*X1*X3 + 0.5*X1^2 + e
```

with `X1, X2 ~ N(0,1)`, `X3 ~ Bernoulli(0.5)`, `e ~ N(0,1)`, all independent. The truth is therefore `E[Y] = 1.5` in closed form. Covariates are always observed. Three mechanisms govern the response indicator:

| Mechanism | Response model |
|---|---|
| MCAR | `P(R = 1) = 0.7` |
| MAR | `logit P(R = 1 given X) = 0.8 + 0.9*X1 - 0.6*X2 + 0.5*X3` |
| MNAR | `logit P(R = 1 given X, Y) = 0.5*X1 + 0.8*Y` |

Estimators compared:

| Estimator | What it assumes |
|---|---|
| `full_data` | Benchmark analysis before deletion |
| `complete_case` | Mean of observed outcomes |
| `mean_imputation` | Constant fill, analysed as complete data |
| `mi_main` | MI with a linear imputation model, main effects only |
| `mi_correct` | MI with the quadratic and interaction terms included |
| `ipw` | Correct parametric propensity, Hajek weighting |
| `aipw_glm` | Correct propensity, misspecified outcome model, doubly robust |
| `dml_rf` | Random forest nuisances, 5-fold cross-fitting, efficient influence function |

### Scenario B: missing covariate, congeniality

The substantive model contains an interaction:

```
Y = 1 + X1 + 0.5*X2 + 0.8*X1*X2 + e
logit P(R_X1 = 1) = 1.0 - 0.7*Y + 0.4*X2
```

`X1` is missing at random given the observed data, so complete-case analysis is biased because the missingness depends on the outcome. An imputation model for `X1` built from main effects of `Y` and `X2` cannot reproduce the interaction in the analysis model. Estimators compared: `complete_case`, `mean_imputation`, `mice_main`, `mice_interaction` (product term as just another variable), and `smcfcs` (substantive model compatible imputation).

### Sensitivity analysis

Under MNAR none of these estimators is consistent. The study reports a delta-adjustment sweep instead: impute under MAR, shift every imputed outcome by delta, and trace the estimate across a grid of delta values.

## Running it

```bash
Rscript scripts/install-deps.R   # mice, smcfcs, ranger, ggplot2, knitr, rmarkdown
make smoke                       # 3 replications, checks every estimator returns finite output
make quick                       # QUICK=1, about a minute
make all                         # full study, then renders the report
```

Runtime for the full study depends on cores. At `n = 1000` and `n_sim = 500` across three mechanisms, budget roughly 20 to 40 minutes on 4 to 8 cores. `parallel::mclapply` handles the parallelism on macOS and Linux; Windows falls back to sequential execution.

Outputs:

- `results/scenario-a-raw.rds`, `results/scenario-b-raw.rds`, `results/sensitivity-raw.rds`
- `results/scenario-a-summary.csv`, `results/scenario-b-summary.csv`
- `figures/*.png`
- `results/session-info.txt` with the base seed, run settings and package versions
- `report/report.html`

## Reproducibility

`draw_seeds()` generates one seed per replication from a single base seed (`config$base_seed`, default 20271001). Each replication sets its own seed before generating data, so results do not depend on the number of cores or the order in which replications finish. Changing `n`, `n_sim` or `m_imp` changes the results; changing `n_cores` does not.

Every performance measure comes with a Monte Carlo standard error. Coverage of 0.93 against a nominal 0.95 means little without knowing whether the Monte Carlo error is 0.005 or 0.03.

## Structure

```
R/00-config.R            settings, seeds, truth values
R/utils.R                seed handling, Rubin's rules, performance summaries
R/01-dgp.R               data generating processes
R/02-estimators-mean.R   scenario A estimators and the delta adjustment
R/03-estimators-coef.R   scenario B estimators
R/04-run-scenario-a.R    scenario A runner
R/05-run-scenario-b.R    scenario B runner and sensitivity sweep
R/06-summarise.R         tables and figures
scripts/install-deps.R   package installation
scripts/run-all.R        full pipeline
tests/test-smoke.R       fast check before committing to a full run
report/report.Rmd        write-up, rebuilt from saved results
```

Rubin's rules are written out longhand in `R/utils.R` rather than called from a package. The split between within- and between-imputation variance, and the Barnard-Rubin degrees of freedom, are the part worth being able to derive.

## Known limitations

- Standard errors for `ipw` and `aipw_glm` treat the fitted propensity score as known. That is conservative when the propensity model is correct. A sandwich estimator accounting for the first-stage estimation is the obvious extension.
- The debiased ML estimator uses a single learner (random forests). A super learner ensemble over several candidates would be closer to standard practice.
- Only one sample size is studied. Performance of the cross-fitted estimator depends on n through the nuisance convergence rates, so a sweep over n would say more.
- The MNAR sensitivity analysis applies a single constant delta. Pattern mixture models with covariate-dependent departures are the next step.

## References

- Rubin, D. B. (1976). Inference and missing data. *Biometrika* 63(3), 581-592.
- Little, R. J. A. and Rubin, D. B. (2019). *Statistical Analysis with Missing Data*, 3rd edition. Wiley.
- Meng, X.-L. (1994). Multiple-imputation inferences with uncongenial sources of input. *Statistical Science* 9(4), 538-558.
- Bartlett, J. W., Seaman, S. R., White, I. R. and Carpenter, J. R. (2015). Multiple imputation of covariates by fully conditional specification: accommodating the substantive model. *Statistical Methods in Medical Research* 24(4), 462-487.
- van Buuren, S. and Groothuis-Oudshoorn, K. (2011). mice: multivariate imputation by chained equations in R. *Journal of Statistical Software* 45(3).
- Chernozhukov, V. et al. (2018). Double/debiased machine learning for treatment and structural parameters. *The Econometrics Journal* 21(1), C1-C68.
- Kennedy, E. H. (2024). Semiparametric doubly robust targeted double machine learning: a review.
- Morris, T. P., White, I. R. and Crowther, M. J. (2019). Using simulation studies to evaluate statistical methods. *Statistics in Medicine* 38(11), 2074-2102.

## Licence

MIT. See `LICENSE`.
