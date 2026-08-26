#' Estimate true prevalence from observed test positives
#'
#' Returns a bias-corrected prevalence estimate and confidence interval.
#' Apparent prevalence (positives / n) is corrected for imperfect diagnostic
#' sensitivity and specificity using the Rogan-Gladen method.  The Wilson
#' score interval is computed on the apparent prevalence scale and then
#' linearly transformed to the true-prevalence scale.
#'
#' Optional clustering arguments (`n_per_site`, `icc`) widen the CI via the
#' Kish design effect.  Optional `fpc_N` narrows it via a finite-population
#' correction.
#'
#' @param positives Integer. Observed number of positive test results.
#' @param n Integer. Total number of samples tested.
#' @param sensitivity Numeric in (0, 1]. Diagnostic sensitivity; default 1
#'   (perfect test).
#' @param specificity Numeric in (0, 1]. Diagnostic specificity; default 1
#'   (perfect test).
#' @param conf_level Numeric in (0, 1). Confidence level; default 0.95.
#' @param n_per_site Optional integer. Samples per site (cluster size).
#'   Used to compute the Kish design effect together with `icc`.
#' @param n_sites Optional integer. Number of sites / clusters (informational
#'   only; required alongside `n_per_site` for FPC).
#' @param icc Numeric >= 0. Intraclass correlation coefficient; default 0.
#' @param fpc_N Optional integer. Finite population size.  When supplied,
#'   requires `n_per_site` and `n_sites` to compute the sampled fraction.
#'
#' @return An object of class `"mms_estimate"` (a named list) with:
#'   \item{prevalence}{Rogan-Gladen corrected true-prevalence point estimate}
#'   \item{ci}{Named numeric vector `c(lower, upper)` on true-prevalence scale}
#'   \item{conf_level}{Confidence level used}
#'   \item{apparent_prev}{Observed apparent prevalence (positives / n)}
#'   \item{apparent_ci}{Wilson CI on apparent prevalence scale}
#'   \item{n}{Total samples}
#'   \item{positives}{Observed positives}
#'   \item{sensitivity}{Sensitivity used}
#'   \item{specificity}{Specificity used}
#'   \item{design_effect}{Kish DEFF applied (1 = no clustering)}
#'
#' @references
#' Rogan WJ, Gladen B (1978). "Estimating prevalence from the results of a
#' screening test." *Am J Epidemiol* 107(1):71-76.
#'
#' MMS-SD workshop Module 2 (margin-of-error estimation) and Module 5
#' (design effect / ICC).
#'
#' @export
#' @examples
#' # Perfect test: 30 positives from 100 samples
#' estimate_prevalence(positives = 30, n = 100)
#'
#' # Imperfect test (sensitivity 90 %, specificity 95 %)
#' estimate_prevalence(30, 100, sensitivity = 0.9, specificity = 0.95)
#'
#' # With clustering: 10 samples per site, ICC = 0.05
#' estimate_prevalence(30, 100, n_per_site = 10, icc = 0.05)
estimate_prevalence <- function(
  positives,
  n,
  sensitivity = 1,
  specificity = 1,
  conf_level  = 0.95,
  n_per_site  = NULL,
  n_sites     = NULL,
  icc         = 0,
  fpc_N       = NULL
) {
  # ---- validation ----
  if (!is.numeric(positives) || length(positives) != 1 || positives < 0)
    stop("`positives` must be a single non-negative number")
  if (!is.numeric(n) || length(n) != 1 || n <= 0)
    stop("`n` must be a single positive number")
  if (positives > n)
    stop("`positives` cannot exceed `n`")
  if (!is.numeric(sensitivity) || sensitivity <= 0 || sensitivity > 1)
    stop("`sensitivity` must be in (0, 1]")
  if (!is.numeric(specificity) || specificity <= 0 || specificity > 1)
    stop("`specificity` must be in (0, 1]")
  correction <- sensitivity + specificity - 1
  if (correction <= 0)
    stop("`sensitivity` + `specificity` must exceed 1 for Rogan-Gladen correction")
  if (!is.numeric(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("`conf_level` must be in (0, 1)")

  # ---- design effect ----
  design <- .sampling_design(n_per_site = n_per_site, n_sites = n_sites, icc = icc)

  # ---- apparent prevalence and Wilson CI ----
  p_app <- positives / n
  z     <- qnorm((1 + conf_level) / 2)
  n_eff <- n / design$deff   # effective n: smaller → wider CI when DEFF > 1

  w      <- .wilson_ci(p_app, n_eff, z)
  ci_app <- c(lower = w$lower, upper = w$upper)

  # ---- finite-population correction (narrows CI) ----
  if (!is.null(fpc_N)) {
    if (is.null(n_per_site) || is.null(n_sites))
      warning("`fpc_N` supplied without `n_per_site` and `n_sites`; using `n` as sampled fraction")
    n_sampled  <- if (!is.null(n_per_site) && !is.null(n_sites)) n_per_site * n_sites else n
    fpc_factor <- sqrt(1 - n_sampled / fpc_N)
    center     <- (ci_app["lower"] + ci_app["upper"]) / 2
    half_width <- (ci_app["upper"] - ci_app["lower"]) / 2 * fpc_factor
    ci_app     <- c(lower = center - half_width, upper = center + half_width)
  }

  # ---- Rogan-Gladen correction: apparent → true prevalence ----
  p_true  <- .rogan_gladen(p_app,   sensitivity, specificity)
  ci_true <- .rogan_gladen(ci_app,  sensitivity, specificity)

  # clamp to [0, 1]
  p_true  <- max(0, min(1, p_true))
  ci_true <- pmax(0, pmin(1, ci_true))
  names(ci_true) <- c("lower", "upper")

  structure(
    list(
      prevalence    = p_true,
      ci            = ci_true,
      conf_level    = conf_level,
      apparent_prev = p_app,
      apparent_ci   = ci_app,
      n             = n,
      positives     = positives,
      sensitivity   = sensitivity,
      specificity   = specificity,
      design_effect = design$deff
    ),
    class = "mms_estimate"
  )
}

#' @export
print.mms_estimate <- function(x, digits = 1, ...) {
  pct <- function(v) sprintf("%.*f%%", digits, 100 * v)
  cat(sprintf(
    "Prevalence: %s  (%s–%s %d%% CI)\n",
    pct(x$prevalence), pct(x$ci["lower"]), pct(x$ci["upper"]),
    round(100 * x$conf_level)
  ))
  if (x$sensitivity < 1 || x$specificity < 1)
    cat(sprintf(
      "  Apparent prevalence %s corrected via Rogan-Gladen (se=%.2f, sp=%.2f)\n",
      pct(x$apparent_prev), x$sensitivity, x$specificity
    ))
  if (x$design_effect > 1)
    cat(sprintf("  Design effect (DEFF): %.2f\n", x$design_effect))
  invisible(x)
}
