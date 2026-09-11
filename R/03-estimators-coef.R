## 03-estimators-coef.R ----------------------------------------------------
## Estimators of the substantive model Y = b0 + b1*X1 + b2*X2 + b3*X1*X2 + e
## when X1 is partially observed and missingness depends on Y.
##
## This is the congeniality scenario. An imputation model for X1 built from
## main effects of Y and X2 cannot reproduce the interaction in the analysis
## model, so the pooled interaction coefficient is pulled toward zero.
## Substantive model compatible imputation (smcfcs) draws X1 from a
## distribution that respects the analysis model instead.

SM_FORMULA <- y ~ x1 * x2
SM_STRING  <- "y ~ x1 + x2 + x1:x2"

## Helper --------------------------------------------------------------------

#' Turn an lm fit into one row per coefficient
lm_results <- function(fit, estimator) {
  cf <- summary(fit)$coefficients
  ci <- stats::confint(fit)
  data.frame(
    estimator = estimator,
    target    = rownames(cf),
    estimate  = cf[, "Estimate"],
    se        = cf[, "Std. Error"],
    ci_lo     = ci[, 1],
    ci_hi     = ci[, 2],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' Pool a list of lm fits across imputations, coefficient by coefficient
pool_lm_fits <- function(fits, estimator, df_com) {
  coefs <- lapply(fits, stats::coef)
  vars  <- lapply(fits, function(f) diag(stats::vcov(f)))
  names_all <- names(coefs[[1]])

  rows <- lapply(names_all, function(nm) {
    est <- vapply(coefs, function(x) x[[nm]], numeric(1))
    v   <- vapply(vars,  function(x) x[[nm]], numeric(1))
    out <- rubin_rules(est, v, df_com = df_com)
    data.frame(
      estimator = estimator,
      target    = nm,
      estimate  = out$estimate,
      se        = out$se,
      ci_lo     = out$ci_lo,
      ci_hi     = out$ci_hi,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

## Benchmark -----------------------------------------------------------------

est_b_full_data <- function(dat) {
  work <- data.frame(y = dat$y, x1 = dat$x1_full, x2 = dat$x2)
  lm_results(stats::lm(SM_FORMULA, data = work), "full_data")
}

## Ad hoc methods ------------------------------------------------------------

est_b_complete_case <- function(dat) {
  work <- dat[dat$r == 1, c("y", "x1", "x2")]
  lm_results(stats::lm(SM_FORMULA, data = work), "complete_case")
}

est_b_mean_imputation <- function(dat) {
  work    <- dat[, c("y", "x1", "x2")]
  work$x1 <- ifelse(is.na(work$x1), mean(work$x1, na.rm = TRUE), work$x1)
  lm_results(stats::lm(SM_FORMULA, data = work), "mean_imputation")
}

## Multiple imputation -------------------------------------------------------

#' Standard MICE for X1
#'
#' @param include_interaction when TRUE, a Y-by-X2 product term enters the
#'   imputation model as a just-another-variable predictor. This narrows the
#'   gap with the substantive model without closing it.
est_b_mice <- function(dat, m, include_interaction = FALSE, seed = NA) {
  work <- dat[, c("y", "x1", "x2")]
  if (include_interaction) work$yx2 <- dat$y * dat$x2

  meth       <- mice::make.method(work)
  meth[]     <- ""
  meth["x1"] <- "norm"

  pred     <- mice::make.predictorMatrix(work)
  pred[, ] <- 0
  pred["x1", setdiff(names(work), "x1")] <- 1

  imp <- mice::mice(work, m = m, method = meth, predictorMatrix = pred,
                    printFlag = FALSE, seed = seed)

  fits <- lapply(mice::complete(imp, action = "all"), function(d) {
    stats::lm(SM_FORMULA, data = d)
  })

  label <- if (include_interaction) "mice_interaction" else "mice_main"
  pool_lm_fits(fits, label, df_com = nrow(work) - 4)
}

#' Substantive model compatible fully conditional specification
#'
#' smcfcs draws X1 from a distribution proportional to
#' f(Y | X1, X2) * f(X1 | X2), so the interaction in the analysis model enters
#' the imputation step by construction.
est_b_smcfcs <- function(dat, m, numit) {
  work <- dat[, c("y", "x1", "x2")]

  meth <- rep("", ncol(work))
  meth[match("x1", names(work))] <- "norm"

  imps <- smcfcs::smcfcs(
    originaldata = work,
    smtype       = "lm",
    smformula    = SM_STRING,
    method       = meth,
    m            = m,
    numit        = numit,
    rjlimit      = 5000,
    noisy        = FALSE
  )

  fits <- lapply(imps$impDatasets, function(d) stats::lm(SM_FORMULA, data = d))
  pool_lm_fits(fits, "smcfcs", df_com = nrow(work) - 4)
}
