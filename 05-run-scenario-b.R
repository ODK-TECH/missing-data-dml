## 05-run-scenario-b.R -----------------------------------------------------
## Scenario B: missing covariate, target the substantive model coefficients.
## Also holds the MNAR sensitivity sweep, which reuses the scenario A data.

run_one_b <- function(seed, config) {
  set.seed(seed)
  dat <- sim_data_b(config$n)

  res <- rbind(
    safe_call("full_data",        est_b_full_data(dat)),
    safe_call("complete_case",    est_b_complete_case(dat)),
    safe_call("mean_imputation",  est_b_mean_imputation(dat)),
    safe_call("mice_main",        est_b_mice(dat, m = config$m_imp,
                                             include_interaction = FALSE)),
    safe_call("mice_interaction", est_b_mice(dat, m = config$m_imp,
                                             include_interaction = TRUE)),
    safe_call("smcfcs",           est_b_smcfcs(dat, m = config$m_imp,
                                               numit = config$numit))
  )

  res$mechanism     <- "MAR_covariate"
  res$seed          <- seed
  res$response_rate <- response_rate(dat, "x1")
  res
}

run_scenario_b <- function(config) {
  seeds   <- draw_seeds(config$n_sim, config$base_seed + 1L)
  started <- Sys.time()
  message("Scenario B | MAR covariate | ", config$n_sim, " replications")

  reps <- run_replications(
    seeds,
    fun = function(seed, i) try(run_one_b(seed, config), silent = TRUE),
    n_cores = config$n_cores
  )

  out <- bind_replications(reps, label = "Scenario B")
  message("  finished in ",
          round(difftime(Sys.time(), started, units = "mins"), 2), " minutes")
  out
}

## MNAR sensitivity ---------------------------------------------------------
##
## Under MNAR no estimator in this study recovers E[Y], so the question shifts
## from "which method is right" to "how far from MAR would the data have to be
## to overturn the answer". The sweep runs on a smaller number of replications
## because it exists to show the shape of the curve.

run_sensitivity <- function(config, n_rep = 50L) {
  n_rep <- min(n_rep, config$n_sim)
  seeds <- draw_seeds(config$n_sim, config$base_seed)[seq_len(n_rep)]

  message("Sensitivity sweep | ", length(config$delta_grid),
          " delta values | ", n_rep, " replications")

  reps <- run_replications(
    seeds,
    fun = function(seed, i) {
      try({
        set.seed(seed)
        dat <- sim_data_a(config$n, "MNAR")
        rows <- lapply(config$delta_grid, function(d) {
          est_delta_adjusted(dat, m = config$m_imp, delta = d)
        })
        out <- do.call(rbind, rows)
        out$seed <- seed
        out
      }, silent = TRUE)
    },
    n_cores = config$n_cores
  )

  bind_replications(reps, label = "Sensitivity sweep")
}
