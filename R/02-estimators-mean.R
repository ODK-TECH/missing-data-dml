## 02-estimators-mean.R -----------------------------------------------------
## Estimators of theta = E[Y] when Y is partially observed and X is complete.
## Every function takes the simulated data frame and returns a one-row data
## frame with estimate, standard error and a 95 per cent interval.

## Benchmark -----------------------------------------------------------------

#' Analysis of the data before any values were deleted
est_full_data <- function(dat) {
  y  <- dat$y_full
  n  <- length(y)
  as_result("full_data", "theta",
            wald_ci(mean(y), stats::sd(y) / sqrt(n)))
}

## Ad hoc methods ------------------------------------------------------------

#' Complete-case analysis: the mean of the observed outcomes
est_complete_case <- function(dat) {
  y  <- dat$y[dat$r == 1]
  n  <- length(y)
  as_result("complete_case", "theta",
            wald_ci(mean(y), stats::sd(y) / sqrt(n)))
}

#' Mean imputation, analysed as if the filled values were real
#'
#' The point estimate matches complete-case analysis exactly. The standard
#' error does not, because filling in a constant shrinks the sample variance
#' while the denominator grows to n. Coverage collapses as a result.
est_mean_imputation <- function(dat) {
  y_obs  <- dat$y[dat$r == 1]
  y_fill <- ifelse(dat$r == 1, dat$y, mean(y_obs))
  n      <- length(y_fill)
  as_result("mean_imputation", "theta",
            wald_ci(mean(y_fill), stats::sd(y_fill) / sqrt(n)))
}

## Multiple imputation -------------------------------------------------------

#' Multiple imputation of Y under a normal linear imputation model
#'
#' @param terms "main" uses X1, X2, X3 as main effects, which omits the
#'   quadratic and interaction terms in the true conditional mean. "correct"
#'   adds them. The contrast shows what a misspecified imputation model costs
#'   even when the missingness really is MAR.
est_mi_mice <- function(dat, m, terms = c("main", "correct"), seed = NA) {
  terms <- match.arg(terms)

  work <- dat[, c("y", "x1", "x2", "x3")]
  if (terms == "correct") {
    work$x1sq <- dat$x1^2
    work$x1x3 <- dat$x1 * dat$x3
  }

  meth      <- mice::make.method(work)
  meth[]    <- ""
  meth["y"] <- "norm"

  pred        <- mice::make.predictorMatrix(work)
  pred[, ]    <- 0
  pred["y", setdiff(names(work), "y")] <- 1

  imp <- mice::mice(work, m = m, method = meth, predictorMatrix = pred,
                    printFlag = FALSE, seed = seed)

  imp_list <- mice::complete(imp, action = "all")

  ests <- vapply(imp_list, function(d) mean(d$y), numeric(1))
  vars <- vapply(imp_list, function(d) stats::var(d$y) / nrow(d), numeric(1))

  label <- paste0("mi_", terms)
  as_result(label, "theta",
            rubin_rules(ests, vars, df_com = nrow(work) - 1))
}

## Weighting and doubly robust methods --------------------------------------

#' Inverse probability weighting, Hajek form, with a parametric propensity
#'
#' The standard error treats the fitted propensity as known. That is
#' conservative when the propensity model is correct, and it keeps the code
#' readable. A sandwich estimator accounting for the first stage is the
#' next step if you want exact intervals.
est_ipw <- function(dat, trim) {
  ps_fit <- stats::glm(r ~ x1 + x2 + x3, family = stats::binomial(), data = dat)
  pi_hat <- stats::predict(ps_fit, type = "response")
  pi_hat <- pmin(pmax(pi_hat, trim), 1 - trim)

  w     <- dat$r / pi_hat
  y0    <- ifelse(dat$r == 1, dat$y, 0)
  theta <- sum(w * y0) / sum(w)

  psi <- w * (y0 - theta * dat$r) / mean(w)
  se  <- stats::sd(psi) / sqrt(nrow(dat))

  as_result("ipw", "theta", wald_ci(theta, se))
}

#' Augmented IPW with parametric nuisance models
#'
#' The outcome model omits the quadratic and interaction terms, so it is
#' misspecified while the propensity model is correct. Double robustness
#' should keep the estimator consistent under MAR.
est_aipw_glm <- function(dat, trim) {
  ps_fit <- stats::glm(r ~ x1 + x2 + x3, family = stats::binomial(), data = dat)
  pi_hat <- pmin(pmax(stats::predict(ps_fit, type = "response"), trim), 1 - trim)

  obs    <- dat$r == 1
  om_fit <- stats::lm(y ~ x1 + x2 + x3, data = dat[obs, ])
  m_hat  <- stats::predict(om_fit, newdata = dat)

  resid <- ifelse(obs, dat$y - m_hat, 0)
  psi   <- m_hat + (dat$r / pi_hat) * resid

  theta <- mean(psi)
  se    <- stats::sd(psi) / sqrt(nrow(dat))

  as_result("aipw_glm", "theta", wald_ci(theta, se))
}

#' Cross-fitted debiased machine learning
#'
#' Both nuisance functions come from random forests. Cross-fitting splits the
#' sample into K folds, fits the nuisances on the other K-1 folds and
#' evaluates them on the held-out fold, which removes the own-observation bias
#' that makes naive plug-in estimators fail. The estimating function is the
#' efficient influence function for E[Y] under MAR, so it is Neyman orthogonal
#' and tolerates slow nuisance convergence rates.
est_dml <- function(dat, n_folds, num_trees, trim) {
  n      <- nrow(dat)
  covars <- c("x1", "x2", "x3")
  folds  <- sample(rep_len(seq_len(n_folds), n))

  m_hat  <- numeric(n)
  pi_hat <- numeric(n)

  for (k in seq_len(n_folds)) {
    train <- folds != k
    test  <- folds == k

    ps_fit <- ranger::ranger(
      x = dat[train, covars, drop = FALSE],
      y = factor(dat$r[train], levels = c(0, 1)),
      probability = TRUE, num.trees = num_trees
    )
    pi_hat[test] <- stats::predict(
      ps_fit, data = dat[test, covars, drop = FALSE]
    )$predictions[, "1"]

    train_obs <- train & dat$r == 1
    om_fit <- ranger::ranger(
      x = dat[train_obs, covars, drop = FALSE],
      y = dat$y[train_obs],
      num.trees = num_trees
    )
    m_hat[test] <- stats::predict(
      om_fit, data = dat[test, covars, drop = FALSE]
    )$predictions
  }

  pi_hat <- pmin(pmax(pi_hat, trim), 1 - trim)
  resid  <- ifelse(dat$r == 1, dat$y - m_hat, 0)
  psi    <- m_hat + (dat$r / pi_hat) * resid

  theta <- mean(psi)
  se    <- stats::sd(psi) / sqrt(n)

  as_result("dml_rf", "theta", wald_ci(theta, se))
}

## Sensitivity analysis ------------------------------------------------------

#' Delta-adjusted multiple imputation for departures from MAR
#'
#' Impute under MAR, then shift every imputed value by delta before analysis.
#' delta = 0 returns the MAR answer. The sweep answers a different question
#' from the ones above: how far from MAR would the data have to be before the
#' conclusion changes?
est_delta_adjusted <- function(dat, m, delta, seed = NA) {
  work      <- dat[, c("y", "x1", "x2", "x3")]
  meth      <- mice::make.method(work)
  meth[]    <- ""
  meth["y"] <- "norm"

  pred     <- mice::make.predictorMatrix(work)
  pred[, ] <- 0
  pred["y", c("x1", "x2", "x3")] <- 1

  imp      <- mice::mice(work, m = m, method = meth, predictorMatrix = pred,
                         printFlag = FALSE, seed = seed)
  imp_list <- mice::complete(imp, action = "all")
  missing  <- is.na(work$y)

  ests <- vapply(imp_list, function(d) {
    d$y[missing] <- d$y[missing] + delta
    mean(d$y)
  }, numeric(1))

  vars <- vapply(imp_list, function(d) {
    d$y[missing] <- d$y[missing] + delta
    stats::var(d$y) / nrow(d)
  }, numeric(1))

  out <- rubin_rules(ests, vars, df_com = nrow(work) - 1)
  res <- as_result("delta_adjusted_mi", "theta", out)
  res$delta <- delta
  res
}
