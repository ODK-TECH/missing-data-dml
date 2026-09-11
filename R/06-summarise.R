## 06-summarise.R ----------------------------------------------------------
## Performance tables and figures. Bias and coverage carry Monte Carlo
## standard errors so that small differences are not over-read.

estimator_levels_a <- c("full_data", "complete_case", "mean_imputation",
                        "mi_main", "mi_correct", "ipw", "aipw_glm", "dml_rf")

estimator_levels_b <- c("full_data", "complete_case", "mean_imputation",
                        "mice_main", "mice_interaction", "smcfcs")

order_estimators <- function(df, levels) {
  df$estimator <- factor(df$estimator, levels = levels)
  df[order(df$estimator), ]
}

## Tables --------------------------------------------------------------------

write_summary <- function(raw, truth, path, levels) {
  summ <- summarise_performance(raw, truth)
  summ <- order_estimators(summ, levels)
  utils::write.csv(summ, path, row.names = FALSE)
  summ
}

## Figures -------------------------------------------------------------------

plot_bias <- function(summ, title, levels) {
  summ <- order_estimators(summ, levels)
  ggplot2::ggplot(summ, ggplot2::aes(x = estimator, y = bias)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    ggplot2::geom_pointrange(
      ggplot2::aes(ymin = bias - 1.96 * bias_mcse,
                   ymax = bias + 1.96 * bias_mcse)
    ) +
    ggplot2::facet_grid(target ~ mechanism, scales = "free_y") +
    ggplot2::coord_flip() +
    ggplot2::labs(title = title, x = NULL,
                  y = "Bias with Monte Carlo interval") +
    ggplot2::theme_minimal(base_size = 11)
}

plot_coverage <- function(summ, title, levels) {
  summ <- order_estimators(summ, levels)
  ggplot2::ggplot(summ, ggplot2::aes(x = estimator, y = coverage)) +
    ggplot2::geom_hline(yintercept = 0.95, linetype = 2, colour = "grey40") +
    ggplot2::geom_pointrange(
      ggplot2::aes(ymin = coverage - 1.96 * cov_mcse,
                   ymax = coverage + 1.96 * cov_mcse)
    ) +
    ggplot2::facet_grid(target ~ mechanism) +
    ggplot2::coord_flip(ylim = c(0, 1)) +
    ggplot2::labs(title = title, x = NULL,
                  y = "Coverage of nominal 95% intervals") +
    ggplot2::theme_minimal(base_size = 11)
}

plot_sensitivity <- function(sens, truth_value) {
  agg <- stats::aggregate(
    cbind(estimate, ci_lo, ci_hi) ~ delta,
    data = sens, FUN = mean
  )
  ggplot2::ggplot(agg, ggplot2::aes(x = delta, y = estimate)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lo, ymax = ci_hi),
                         alpha = 0.15) +
    ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = truth_value, linetype = 2,
                        colour = "grey40") +
    ggplot2::labs(
      x = "Delta applied to imputed outcomes",
      y = "Estimated E[Y]",
      title = "Delta-adjusted multiple imputation under MNAR",
      subtitle = "Dashed line marks the true value"
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

save_plot <- function(plot, path, width = 9, height = 6) {
  ggplot2::ggsave(path, plot, width = width, height = height, dpi = 200)
  invisible(path)
}
