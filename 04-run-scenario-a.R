## 04-run-scenario-a.R -----------------------------------------------------
## Scenario A: missing outcome, target E[Y], three mechanisms.

run_one_a <- function(seed, mechanism, config) {
  set.seed(seed)
  dat <- sim_data_a(config$n, mechanism)

  res <- rbind(
    safe_call("full_data",       est_full_data(dat)),
    safe_call("complete_case",   est_complete_case(dat)),
    safe_call("mean_imputation", est_mean_imputation(dat)),
    safe_call("mi_main",         est_mi_mice(dat, m = config$m_imp,
                                             terms = "main")),
    safe_call("mi_correct",      est_mi_mice(dat, m = config$m_imp,
                                             terms = "correct")),
    safe_call("ipw",             est_ipw(dat, trim = config$trim)),
    safe_call("aipw_glm",        est_aipw_glm(dat, trim = config$trim)),
    safe_call("dml_rf",          est_dml(dat,
                                         n_folds   = config$n_folds,
                                         num_trees = config$num_trees,
                                         trim      = config$trim))
  )

  res$mechanism     <- mechanism
  res$seed          <- seed
  res$response_rate <- response_rate(dat, "y")
  res
}

run_scenario_a <- function(config,
                           mechanisms = c("MCAR", "MAR", "MNAR")) {
  seeds <- draw_seeds(config$n_sim, config$base_seed)

  all_res <- lapply(mechanisms, function(mech) {
    started <- Sys.time()
    message("Scenario A | ", mech, " | ", config$n_sim, " replications")

    reps <- run_replications(
      seeds,
      fun = function(seed, i) try(run_one_a(seed, mech, config), silent = TRUE),
      n_cores = config$n_cores
    )

    out <- bind_replications(reps, label = paste("Scenario A", mech))
    message("  finished in ",
            round(difftime(Sys.time(), started, units = "mins"), 2), " minutes")
    out
  })

  do.call(rbind, all_res)
}
