#' Calculate sample size to test prevalence against a threshold
#'
#' Returns the sample size needed to detect (with given power) that true
#' prevalence exceeds a decision threshold, using either a frequentist or
#' Bayesian framework.
#'
#' Pairs with `test_threshold()`.  Corresponds to DRpower / WHO pfhrp2/3
#' protocol sample sizing (M6).
#'
#' @param prevalence Numeric in (0, 1). Expected true prevalence (the
#'   prevalence you want to be able to detect).
#' @param threshold Numeric in (0, 1). Decision threshold.  Typically
#'   `prevalence > threshold` is the hypothesis of interest.
#' @param power Numeric in (0, 1). Target statistical power; default 0.80.
#' @param alpha Numeric. Significance level (type-I error rate); default 0.05.
#' @param sensitivity Numeric in (0, 1]. Diagnostic sensitivity; default 1.
#' @param specificity Numeric in (0, 1]. Diagnostic specificity; default 1.
#' @param method Character. `"frequentist"` (default) or `"bayesian"`.
#' @param alternative Character. `"greater"` (default) or `"less"`.
#' @param prior_alpha Numeric. Shape1 of Beta prior (Bayesian only); default 1.
#' @param prior_beta Numeric. Shape2 of Beta prior (Bayesian only); default 1.
#' @param n_per_site Optional integer. Samples per site (design effect).
#' @param n_sites Optional integer. Number of sites.
#' @param icc Numeric >= 0. Intraclass correlation; default 0.
#' @param fpc_N Optional integer. Finite population size.
#'
#' @return An object of class `"mms_design"` (planned) with:
#'   \item{n}{Total sample size required.}
#'   \item{n_base}{Required n ignoring clustering and FPC.}
#'   \item{power}{Target power used.}
#'   \item{alpha}{Significance level used.}
#'   \item{prevalence}{Expected prevalence assumed.}
#'   \item{threshold}{Decision threshold.}
#'   \item{design_effect}{DEFF applied.}
#'
#' @references
#' MMS-SD workshop Module 6 (DRpower / WHO pfhrp2/3 protocol).
#'
#' @export
design_threshold <- function(
  prevalence,
  threshold,
  power       = 0.80,
  alpha       = 0.05,
  sensitivity = 1,
  specificity = 1,
  method      = c("frequentist", "bayesian"),
  alternative = c("greater", "less"),
  prior_alpha = 1,
  prior_beta  = 1,
  n_per_site  = NULL,
  n_sites     = NULL,
  icc         = 0,
  fpc_N       = NULL
) {
  stop("design_threshold() is not yet implemented")
}
