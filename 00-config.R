## 00-config.R -------------------------------------------------------------
## Global settings for the simulation study. Every script sources this file
## first, so all sample sizes, seeds and truth values live in one place.

config <- list(
  ## Reproducibility --------------------------------------------------------
  ## base_seed generates one seed per replication (see draw_seeds() in
  ## utils.R). Each replication calls set.seed() with its own seed, so results
  ## do not depend on the number of cores or on the order of execution.
  base_seed = 20271001L,

  ## Simulation size --------------------------------------------------------
  n         = 1000L,   # observations per simulated dataset
  n_sim     = 500L,    # replications per scenario
  m_imp     = 20L,     # imputations per multiply imputed analysis
  numit     = 10L,     # smcfcs Gibbs iterations

  ## Debiased ML ------------------------------------------------------------
  n_folds   = 5L,      # cross-fitting folds
  num_trees = 500L,    # ranger trees per nuisance model
  trim      = 0.01,    # propensity scores trimmed to [trim, 1 - trim]

  ## MNAR sensitivity analysis ----------------------------------------------
  delta_grid = seq(-1.5, 1.5, by = 0.25),

  ## Execution --------------------------------------------------------------
  n_cores     = max(1L, parallel::detectCores() - 1L),
  results_dir = "results",
  figures_dir = "figures"
)

## Truth ---------------------------------------------------------------------
## Scenario A. Y = 1 + X1 + 0.5*X2 - 0.75*X1*X3 + 0.5*X1^2 + e, with
## X1, X2 ~ N(0, 1), X3 ~ Bernoulli(0.5), e ~ N(0, 1), all mutually
## independent. Therefore E[Y] = 1 + 0 + 0 - 0 + 0.5 * E[X1^2] = 1.5.
truth_a <- c(theta = 1.5)

## Scenario B. Y = 1 + X1 + 0.5*X2 + 0.8*X1*X2 + e is the substantive model,
## so the coefficients are known by construction.
truth_b <- c(
  "(Intercept)" = 1.0,
  "x1"          = 1.0,
  "x2"          = 0.5,
  "x1:x2"       = 0.8
)

## QUICK=1 shrinks the study to something that finishes in about a minute.
## Use it to check that the pipeline runs before committing to the full study.
if (nzchar(Sys.getenv("QUICK"))) {
  config$n_sim  <- 25L
  config$m_imp  <- 5L
  config$numit  <- 5L
  message("QUICK mode: n_sim = ", config$n_sim, ", m_imp = ", config$m_imp)
}
