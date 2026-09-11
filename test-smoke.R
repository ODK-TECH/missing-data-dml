## test-smoke.R ------------------------------------------------------------
## Runs three replications of everything and checks that each estimator
## returns a finite estimate. Run this before starting the full study:
##
##   Rscript tests/test-smoke.R

suppressPackageStartupMessages({
  library(mice); library(smcfcs); library(ranger)
})

source("R/00-config.R")
source("R/utils.R")
source("R/01-dgp.R")
source("R/02-estimators-mean.R")
source("R/03-estimators-coef.R")
source("R/04-run-scenario-a.R")
source("R/05-run-scenario-b.R")

config$n_sim <- 3L
config$n     <- 400L
config$m_imp <- 3L
config$numit <- 3L
config$n_cores <- 1L

check <- function(df, label) {
  bad <- df[!is.finite(df$estimate) | !is.finite(df$se), ]
  if (nrow(bad) > 0) {
    print(unique(bad[, c("estimator", "target")]))
    stop(label, ": non-finite results above.")
  }
  message(label, ": ", nrow(df), " rows, all finite. Estimators: ",
          paste(sort(unique(df$estimator)), collapse = ", "))
}

res_a <- run_scenario_a(config, mechanisms = c("MCAR", "MAR", "MNAR"))
check(res_a, "Scenario A")

res_b <- run_scenario_b(config)
check(res_b, "Scenario B")

sens <- run_sensitivity(config, n_rep = 2L)
check(sens, "Sensitivity sweep")

## Rubin's rules against a hand-computed example
est <- c(1.0, 1.2, 0.8, 1.1, 0.9)
var <- rep(0.04, 5)
out <- rubin_rules(est, var, df_com = 99)
stopifnot(abs(out$estimate - mean(est)) < 1e-12)
stopifnot(abs(out$se^2 - (0.04 + 1.2 * stats::var(est))) < 1e-12)
message("Rubin's rules: total variance matches Ubar + (1 + 1/m) B.")

message("Smoke test passed.")
