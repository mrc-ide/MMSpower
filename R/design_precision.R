#' Calculate sample size for a target margin of error on prevalence
#'
#' Returns the sample size needed to estimate true prevalence with a given
#' margin of error (half-width of confidence interval).
#'
#' The sample size formula is derived from the delta-method variance of the
#' Rogan-Gladen estimator:
#'
#' \deqn{n = \frac{z^2 \, p_{app}(1 - p_{app})}{MOE^2 \cdot (Se + Sp - 1)^2}}
#'
#' where \eqn{p_{app} = p \cdot Se + (1 - p)(1 - Sp)} is the apparent
#' prevalence implied by the expected true prevalence.  When `sensitivity = 1`
#' and `specificity = 1` this reduces to the standard formula
#' \eqn{n = z^2 p(1-p) / MOE^2}.
#'
#' Note: design uses the normal-approximation MOE; `estimate_prevalence()` uses
#' a Wilson CI, so exact agreement for marginal cases is not guaranteed.
#'
#' @param prevalence Numeric in (0, 1). Expected true prevalence.
#' @param moe Numeric in (0, 0.5). Target margin of error (half-width of
#'   confidence interval) on the true-prevalence scale.
#' @param sensitivity Numeric in (0, 1]. Diagnostic sensitivity; default 1.
#' @param specificity Numeric in (0, 1]. Diagnostic specificity; default 1.
#' @param conf_level Numeric in (0, 1). Confidence level; default 0.95.
#' @param n_per_site Optional integer. Target samples per site (cluster size).
#'   Determines the Kish design effect.  If omitted, no clustering is assumed.
#' @param n_sites Optional integer. Informational: passed through to the
#'   return value.  Not used in DEFF calculation unless `n_per_site` is also
#'   given (in which case sites needed = `ceiling(n / n_per_site)`).
#' @param icc Numeric >= 0. Intraclass correlation coefficient; default 0.
#' @param fpc_N Optional integer. Finite population size.  Reduces required n
#'   via the approximation \eqn{n_{fpc} \approx n_{srs} / (1 + n_{srs}/N)}.
#'
#' @return An object of class `"mms_design"` (a named list) with:
#'   \item{n}{Total sample size required (ceiling, after DEFF and FPC)}
#'   \item{n_base}{Required n ignoring clustering and FPC}
#'   \item{n_sites}{Sites needed when `n_per_site` is supplied; else `NULL`}
#'   \item{n_per_site}{Cluster size used in DEFF calculation}
#'   \item{prevalence}{Expected true prevalence assumed}
#'   \item{apparent_prev}{Apparent prevalence used in the formula}
#'   \item{moe}{Target margin of error}
#'   \item{conf_level}{Confidence level}
#'   \item{sensitivity}{Sensitivity used}
#'   \item{specificity}{Specificity used}
#'   \item{design_effect}{Kish DEFF applied (1 = no clustering)}
#'
#' @references
#' MMS-SD workshop Module 2 (MOE sample sizing) and Module 5 (ICC / design
#' effect).
#'
#' @export
#' @examples
#' # Perfect test: samples needed for ±5 % MOE at 30 % true prevalence
#' design_precision(prevalence = 0.3, moe = 0.05)
#'
#' # Imperfect test (sensitivity 90 %, specificity 95 %)
#' design_precision(0.3, 0.05, sensitivity = 0.9, specificity = 0.95)
#'
#' # With clustering: 10 samples per site, ICC = 0.05
#' design_precision(0.3, 0.05, n_per_site = 10, icc = 0.05)
design_precision <- function(
  prevalence,
  moe,
  sensitivity = 1,
  specificity = 1,
  conf_level  = 0.95,
  n_per_site  = NULL,
  n_sites     = NULL,
  icc         = 0,
  fpc_N       = NULL
) {
  # ---- validation ----
  if (!is.numeric(prevalence) || length(prevalence) != 1 ||
      prevalence <= 0 || prevalence >= 1)
    stop("`prevalence` must be a single number in (0, 1)")
  if (!is.numeric(moe) || length(moe) != 1 || moe <= 0 || moe >= 0.5)
    stop("`moe` must be a single number in (0, 0.5)")
  if (!is.numeric(sensitivity) || sensitivity <= 0 || sensitivity > 1)
    stop("`sensitivity` must be in (0, 1]")
  if (!is.numeric(specificity) || specificity <= 0 || specificity > 1)
    stop("`specificity` must be in (0, 1]")
  correction <- sensitivity + specificity - 1
  if (correction <= 0)
    stop("`sensitivity` + `specificity` must exceed 1")
  if (!is.numeric(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("`conf_level` must be in (0, 1)")

  # ---- base sample size (no clustering, no FPC) ----
  z     <- qnorm((1 + conf_level) / 2)
  p_app <- .apparent_prev(prevalence, sensitivity, specificity)

  # Delta-method: Var(p_true) = p_app*(1-p_app) / (n * correction^2)
  # MOE = z * sqrt(Var(p_true))  =>  n = z^2 * p_app*(1-p_app) / (MOE^2 * correction^2)
  n_continuous <- z^2 * p_app * (1 - p_app) / (moe^2 * correction^2)

  # ---- design effect ----
  design <- .sampling_design(n_per_site = n_per_site, n_sites = n_sites, icc = icc)

  # ---- inflate for clustering ----
  n_inflated <- n_continuous * design$deff

  # ---- deflate for finite-population correction ----
  if (!is.null(fpc_N)) {
    # Approximation: n_fpc = n_srs / (1 + n_srs / N)
    n_inflated <- n_inflated / (1 + n_inflated / fpc_N)
  }

  n_total <- ceiling(n_inflated)
  n_base  <- ceiling(n_continuous)

  n_sites_needed <- if (!is.null(n_per_site)) ceiling(n_total / n_per_site) else NULL

  structure(
    list(
      n             = n_total,
      n_base        = n_base,
      n_sites       = n_sites_needed,
      n_per_site    = n_per_site,
      prevalence    = prevalence,
      apparent_prev = p_app,
      moe           = moe,
      conf_level    = conf_level,
      sensitivity   = sensitivity,
      specificity   = specificity,
      design_effect = design$deff
    ),
    class = "mms_design"
  )
}

#' @export
print.mms_design <- function(x, ...) {
  cat(sprintf(
    "Required n = %d  (%d%% CI, MOE = ±%.1f%%)\n",
    x$n, round(100 * x$conf_level), 100 * x$moe
  ))
  if (!is.null(x$n_per_site) && x$design_effect > 1) {
    cat(sprintf(
      "  Base n (no clustering): %d  |  DEFF: %.2f  |  Sites needed: %s\n",
      x$n_base, x$design_effect,
      if (!is.null(x$n_sites)) as.character(x$n_sites) else "n/a"
    ))
  }
  if (x$sensitivity < 1 || x$specificity < 1)
    cat(sprintf(
      "  Imperfect test (se=%.2f, sp=%.2f): apparent prevalence = %.1f%%\n",
      x$sensitivity, x$specificity, 100 * x$apparent_prev
    ))
  invisible(x)
}
