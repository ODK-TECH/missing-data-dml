## run-all.R ---------------------------------------------------------------
## Runs both scenarios and the sensitivity sweep, then writes every result and
## figure to disk. Call it from the repository root:
##
##   Rscript scripts/run-all.R          full study
##   QUICK=1 Rscript scripts/run-all.R  short version for checking the pipeline

suppressPackageStartupMessages({
  library(mice)
  library(smcfcs)
  library(ranger)
  library(ggplot2)
})

source("R/00-config.R")
source("R/utils.R")
source("R/01-dgp.R")
source("R/02-estimators-mean.R")
source("R/03-estimators-coef.R")
source("R/04-run-scenario-a.R")
source("R/05-run-scenario-b.R")
source("R/06-summarise.R")

dir.create(config$results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(config$figures_dir, showWarnings = FALSE, recursive = TRUE)

started <- Sys.time()

## Scenario A ---------------------------------------------------------------
raw_a <- run_scenario_a(config)
saveRDS(raw_a, file.path(config$results_dir, "scenario-a-raw.rds"))

summ_a <- write_summary(
  raw_a, truth_a,
  path   = file.path(config$results_dir, "scenario-a-summary.csv"),
  levels = estimator_levels_a
)

save_plot(plot_bias(summ_a, "Scenario A: bias in the estimated mean of Y",
                    estimator_levels_a),
          file.path(config$figures_dir, "scenario-a-bias.png"))

save_plot(plot_coverage(summ_a, "Scenario A: interval coverage",
                        estimator_levels_a),
          file.path(config$figures_dir, "scenario-a-coverage.png"))

## Scenario B ---------------------------------------------------------------
raw_b <- run_scenario_b(config)
saveRDS(raw_b, file.path(config$results_dir, "scenario-b-raw.rds"))

summ_b <- write_summary(
  raw_b, truth_b,
  path   = file.path(config$results_dir, "scenario-b-summary.csv"),
  levels = estimator_levels_b
)

save_plot(plot_bias(summ_b, "Scenario B: bias in the substantive model coefficients",
                    estimator_levels_b),
          file.path(config$figures_dir, "scenario-b-bias.png"),
          height = 8)

save_plot(plot_coverage(summ_b, "Scenario B: interval coverage",
                        estimator_levels_b),
          file.path(config$figures_dir, "scenario-b-coverage.png"),
          height = 8)

## Sensitivity --------------------------------------------------------------
sens <- run_sensitivity(config)
saveRDS(sens, file.path(config$results_dir, "sensitivity-raw.rds"))

save_plot(plot_sensitivity(sens, truth_a[["theta"]]),
          file.path(config$figures_dir, "sensitivity-delta.png"),
          width = 7, height = 5)

## Provenance ---------------------------------------------------------------
writeLines(
  c(
    paste("Run completed:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("Base seed:", config$base_seed),
    paste("n:", config$n, " n_sim:", config$n_sim, " m_imp:", config$m_imp),
    paste("Elapsed (mins):",
          round(difftime(Sys.time(), started, units = "mins"), 2)),
    "",
    utils::capture.output(utils::sessionInfo())
  ),
  file.path(config$results_dir, "session-info.txt")
)

message("Done. Results in ", config$results_dir,
        ", figures in ", config$figures_dir, ".")
