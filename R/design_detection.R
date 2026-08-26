#' Calculate sample size to detect at least one occurrence of a variant
#'
#' Returns the sample size needed to achieve a target probability of observing
#' at least one positive test result, given an expected (possibly rare) true
#' prevalence.  Accounts for imperfect diagnostic sensitivity and specificity.
#'
#' Corresponds to DRpower presence/detection logic (M6).
#'
#' @param prevalence Numeric in (0, 1). Minimum true prevalence of the variant
#'   to detect (i.e., the prevalence at which detection probability is
#'   evaluated).
#' @param detection_prob Numeric in (0, 1). Target probability of observing
#'   at least one positive; default 0.95.
#' @param sensitivity Numeric in (0, 1]. Diagnostic sensitivity; default 1.
#' @param specificity Numeric in (0, 1]. Diagnostic specificity; default 1.
#' @param n_per_site Optional integer. Samples per site (design effect).
#' @param n_sites Optional integer. Number of sites.
#' @param icc Numeric >= 0. Intraclass correlation; default 0.
#' @param fpc_N Optional integer. Finite population size.
#'
#' @details
#' Without diagnostic error, the probability of detecting ≥1 positive in `n`
#' samples is \eqn{1 - (1 - p)^n}, giving
#' \eqn{n = \log(1 - detection\_prob) / \log(1 - p)}.
#'
#' With imperfect tests, the apparent prevalence \eqn{p_{app}} replaces
#' \eqn{p} in the detection formula.  The true-positive detection probability
#' is adjusted accordingly.
#'
#' @return An object of class `"mms_design"` (planned) with:
#'   \item{n}{Total sample size required.}
#'   \item{n_base}{Required n ignoring clustering and FPC.}
#'   \item{prevalence}{Minimum prevalence to detect.}
#'   \item{detection_prob}{Target detection probability.}
#'   \item{apparent_prev}{Apparent prevalence used in the formula.}
#'   \item{design_effect}{DEFF applied.}
#'
#' @references
#' MMS-SD workshop Module 6 (DRpower presence/detection logic).
#'
#' @export
design_detection <- function(
  prevalence,
  detection_prob = 0.95,
  sensitivity    = 1,
  specificity    = 1,
  n_per_site     = NULL,
  n_sites        = NULL,
  icc            = 0,
  fpc_N          = NULL
) {
  stop("design_detection() is not yet implemented")
}
