## utils.R 
## Shared helpers: seed handling, Rubin's rules, interval construction and
## the Monte Carlo summaries used to judge every estimator.

## Seeds

#' Draw one seed per replication
#'
#' Drawing the seeds up front and setting them inside each replication keeps
#' the study reproducible whether it runs sequentially or in parallel.
draw_seeds <- function(n_sim, base_seed) {
  set.seed(base_seed)
  sample.int(.Machine$integer.max, n_sim)
}

## Rubin's rules 

#' Combine estimates across multiply imputed datasets
#'
#' @param estimates numeric vector of m point estimates
#' @param variances numeric vector of m squared standard errors
#' @param df_com complete-data degrees of freedom (n minus parameters)
#'
#' Returns the pooled estimate, its standard error, the Barnard-Rubin degrees
#' of freedom, the fraction of missing information and a 95 per cent interval.
#' Written out longhand rather than called from a package, because the
#' between- and within-imputation split is the part worth understanding.
rubin_rules <- function(estimates, variances, df_com = Inf, level = 0.95) {
  m <- length(estimates)
  stopifnot(m >= 2, length(variances) == m)

  q_bar <- mean(estimates)          # pooled point estimate
  u_bar <- mean(variances)          # within-imputation variance
  b     <- stats::var(estimates)    # between-imputation variance
  t_var <- u_bar + (1 + 1 / m) * b  # total variance

  lambda <- ((1 + 1 / m) * b) / t_var          # fraction of missing information
  lambda <- max(lambda, .Machine$double.eps)

  df_old <- (m - 1) / lambda^2
  if (is.finite(df_com)) {
    df_obs <- ((df_com + 1) / (df_com + 3)) * df_com * (1 - lambda)
    df     <- df_old * df_obs / (df_old + df_obs)
  } else {
    df <- df_old
  }

  se   <- sqrt(t_var)
  crit <- stats::qt(1 - (1 - level) / 2, df = df)

  list(
    estimate = q_bar,
    se       = se,
    df       = df,
    fmi      = lambda,
    ci_lo    = q_bar - crit * se,
    ci_hi    = q_bar + crit * se
  )
}

## Intervals 

#' Wald interval on the normal scale
wald_ci <- function(estimate, se, level = 0.95) {
  crit <- stats::qnorm(1 - (1 - level) / 2)
  list(
    estimate = estimate,
    se       = se,
    ci_lo    = estimate - crit * se,
    ci_hi    = estimate + crit * se
  )
}

#' Standard shape for every estimator in this study
as_result <- function(estimator, target, fit) {
  data.frame(
    estimator = estimator,
    target    = target,
    estimate  = fit$estimate,
    se        = fit$se,
    ci_lo     = fit$ci_lo,
    ci_hi     = fit$ci_hi,
    stringsAsFactors = FALSE
  )
}

## Monte Carlo summaries 

#' Performance measures with Monte Carlo standard errors
#'
#' Bias, empirical standard deviation, RMSE, mean model standard error and
#' interval coverage. The Monte Carlo standard errors matter: without them a
#' coverage of 0.93 against 0.95 cannot be told apart from noise.
summarise_performance <- function(results, truth) {
  truth_df <- data.frame(
    target = names(truth),
    truth  = as.numeric(truth),
    stringsAsFactors = FALSE
  )

  merged <- merge(results, truth_df, by = "target")
  merged <- merged[stats::complete.cases(merged$estimate), , drop = FALSE]

  split_key <- interaction(merged$estimator, merged$target,
                           merged$mechanism, drop = TRUE, sep = "|")

  out <- lapply(split(merged, split_key), function(d) {
    n_rep    <- nrow(d)
    err      <- d$estimate - d$truth
    covered  <- d$ci_lo <= d$truth & d$truth <= d$ci_hi
    cov_rate <- mean(covered)

    data.frame(
      mechanism    = d$mechanism[1],
      estimator    = d$estimator[1],
      target       = d$target[1],
      n_rep        = n_rep,
      truth        = d$truth[1],
      mean_est     = mean(d$estimate),
      bias         = mean(err),
      bias_mcse    = stats::sd(d$estimate) / sqrt(n_rep),
      emp_sd       = stats::sd(d$estimate),
      mean_se      = mean(d$se),
      rmse         = sqrt(mean(err^2)),
      coverage     = cov_rate,
      cov_mcse     = sqrt(cov_rate * (1 - cov_rate) / n_rep),
      ci_width     = mean(d$ci_hi - d$ci_lo),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[order(out$mechanism, out$target, out$estimator), ]
}

## Execution 

#' Run replications, in parallel where the platform allows it
run_replications <- function(seeds, fun, n_cores = 1L) {
  if (n_cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_along(seeds), function(i) fun(seeds[i], i),
                       mc.cores = n_cores)
  } else {
    lapply(seq_along(seeds), function(i) fun(seeds[i], i))
  }
}

#' Drop replications that failed and report how many
bind_replications <- function(reps, label = "") {
  ok      <- !vapply(reps, inherits, logical(1), what = "try-error")
  n_fail  <- sum(!ok)
  if (n_fail > 0) {
    message(label, ": ", n_fail, " of ", length(reps),
            " replications failed and were dropped.")
  }
  do.call(rbind, reps[ok])
}

#' Evaluate one estimator, returning an empty row if it fails
#'
#' A single failing estimator should not discard the whole replication, so
#' failures are recorded rather than raised.
safe_call <- function(estimator, expr) {
  out <- try(expr, silent = TRUE)
  if (inherits(out, "try-error")) {
    message("Estimator ", estimator, " failed: ",
            conditionMessage(attr(out, "condition")))
    return(data.frame(
      estimator = estimator,
      target    = NA_character_,
      estimate  = NA_real_,
      se        = NA_real_,
      ci_lo     = NA_real_,
      ci_hi     = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  out
}
