## 01-dgp.R ----------------------------------------------------------------
## Data generating processes. Both scenarios return the fully observed data
## alongside the incomplete data, so bias can be measured against a known
## truth and complete-data benchmarks stay available.

## Scenario A: missing outcome, target E[Y] ---------------------------------
##
## Y depends on X1 through a quadratic term and interacts with X3, so an
## imputation model with main effects only is misspecified. The covariates are
## always observed. Three mechanisms govern the response indicator R:
##
##   MCAR   P(R = 1) = 0.7
##   MAR    logit P(R = 1 | X) = 0.8 + 0.9*X1 - 0.6*X2 + 0.5*X3
##   MNAR   logit P(R = 1 | X, Y) = 0.0 + 0.5*X1 + 0.8*Y
##
## Under MNAR every estimator here assumes the wrong thing, which is the
## point: the MAR-based toolkit breaks and sensitivity analysis takes over.

sim_data_a <- function(n, mechanism = c("MCAR", "MAR", "MNAR")) {
  mechanism <- match.arg(mechanism)

  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  x3 <- stats::rbinom(n, size = 1, prob = 0.5)

  y_full <- 1 + x1 + 0.5 * x2 - 0.75 * x1 * x3 + 0.5 * x1^2 + stats::rnorm(n)

  eta <- switch(
    mechanism,
    MCAR = rep(stats::qlogis(0.7), n),
    MAR  = 0.8 + 0.9 * x1 - 0.6 * x2 + 0.5 * x3,
    MNAR = 0.0 + 0.5 * x1 + 0.8 * y_full
  )

  r <- stats::rbinom(n, size = 1, prob = stats::plogis(eta))

  data.frame(
    y      = ifelse(r == 1, y_full, NA_real_),
    x1     = x1,
    x2     = x2,
    x3     = x3,
    r      = r,
    y_full = y_full
  )
}

## Scenario B: missing covariate, target regression coefficients ------------
##
## The substantive model carries an interaction:
##
##   Y = 1 + X1 + 0.5*X2 + 0.8*X1*X2 + e
##
## X1 goes missing with probability depending on Y and X2, so the missingness
## is MAR given the observed data but complete-case analysis is biased. An
## imputation model for X1 with main effects only is uncongenial with the
## substantive model, which is the gap smcfcs closes.

sim_data_b <- function(n) {
  x1 <- stats::rnorm(n)
  x2 <- stats::rbinom(n, size = 1, prob = 0.5)

  y <- 1 + x1 + 0.5 * x2 + 0.8 * x1 * x2 + stats::rnorm(n)

  eta <- 1.0 - 0.7 * y + 0.4 * x2
  r   <- stats::rbinom(n, size = 1, prob = stats::plogis(eta))

  data.frame(
    y       = y,
    x1      = ifelse(r == 1, x1, NA_real_),
    x2      = x2,
    r       = r,
    x1_full = x1
  )
}

## Reporting helper ----------------------------------------------------------

#' Observed-data proportion, useful for checking that mechanisms are
#' comparable in how much information they remove
response_rate <- function(dat, var) mean(!is.na(dat[[var]]))
