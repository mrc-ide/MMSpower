#' Test whether prevalence exceeds a decision threshold
#'
#' Given observed positives and sample size, tests whether the true prevalence
#' is above (or below) a specified decision threshold using either a
#' frequentist one-sided test or a Bayesian posterior-probability approach.
#'
#' Corresponds to the WHO pfhrp2/3 deletion surveillance protocol implemented
#' in DRpower (M6).
#'
#' @param positives Integer. Observed positive test results.
#' @param n Integer. Total samples tested.
#' @param threshold Numeric in (0, 1). Decision threshold on true prevalence.
#' @param sensitivity Numeric in (0, 1]. Diagnostic sensitivity; default 1.
#' @param specificity Numeric in (0, 1]. Diagnostic specificity; default 1.
#' @param method Character. `"frequentist"` (default) or `"bayesian"`.
#' @param alternative Character. `"greater"` (default) or `"less"`.
#'   Tests whether prevalence is above or below `threshold`.
#' @param alpha Numeric. Significance level for frequentist test; default 0.05.
#' @param prior_alpha Numeric. Shape1 of Beta prior (Bayesian only); default 1.
#' @param prior_beta Numeric. Shape2 of Beta prior (Bayesian only); default 1.
#' @param conf_level Numeric in (0, 1). CI to report alongside test; default 0.95.
#' @param n_per_site Optional integer. Samples per site (for design effect).
#' @param n_sites Optional integer. Number of sites.
#' @param icc Numeric >= 0. Intraclass correlation; default 0.
#' @param fpc_N Optional integer. Finite population size.
#'
#' @return An object of class `"mms_test_threshold"` (planned) with:
#'   \item{reject}{Logical. `TRUE` if the data support prevalence on the
#'     `alternative` side of `threshold` at level `alpha`.}
#'   \item{p_value}{Frequentist p-value (frequentist method only).}
#'   \item{posterior_prob}{P(prevalence > threshold | data) (Bayesian only).}
#'   \item{prevalence_est}{Rogan-Gladen corrected point estimate.}
#'   \item{ci}{Confidence / credible interval on true prevalence.}
#'   \item{method}{Method used.}
#'
#' @references
#' MMS-SD workshop Module 6 (DRpower / Bayesian WHO pfhrp2/3 protocol).
#'
#' @export
test_threshold <- function(
  positives,
  n,
  threshold,
  sensitivity = 1,
  specificity = 1,
  method      = c("frequentist", "bayesian"),
  alternative = c("greater", "less"),
  alpha       = 0.05,
  prior_alpha = 1,
  prior_beta  = 1,
  conf_level  = 0.95,
  n_per_site  = NULL,
  n_sites     = NULL,
  icc         = 0,
  fpc_N       = NULL
) {
  stop("test_threshold() is not yet implemented")
}
