#' Test whether prevalence is changing over time
#'
#' Given observed positives and sample sizes at multiple timepoints, tests
#' whether there is a statistically significant trend (increase or decrease)
#' in true prevalence over time.
#'
#' Corresponds to longitudinal power/sample-size formulae (M4) with optional
#' clustering (M5).
#'
#' @param positives Integer vector. Observed positives at each timepoint.
#' @param n Integer vector (same length as `positives`). Samples at each
#'   timepoint.
#' @param timepoints Numeric vector (same length as `positives`). Time values
#'   (e.g. years or rounds).  If omitted, defaults to `seq_along(positives)`.
#' @param sensitivity Numeric in (0, 1]. Diagnostic sensitivity; default 1.
#' @param specificity Numeric in (0, 1]. Diagnostic specificity; default 1.
#' @param alternative Character. `"two.sided"` (default), `"increasing"`, or
#'   `"decreasing"`.
#' @param alpha Numeric. Significance level; default 0.05.
#' @param conf_level Numeric in (0, 1). CI on slope; default 0.95.
#' @param n_per_site Optional integer. Samples per site (design effect, applied
#'   to all timepoints).
#' @param n_sites Optional integer. Number of sites.
#' @param icc Numeric >= 0. Intraclass correlation; default 0.
#' @param fpc_N Optional integer. Finite population size.
#'
#' @return An object of class `"mms_test_trend"` (planned) with:
#'   \item{reject}{Logical. `TRUE` if trend is significant at `alpha`.}
#'   \item{p_value}{P-value for trend.}
#'   \item{slope}{Estimated change in true prevalence per unit time.}
#'   \item{ci_slope}{CI on slope.}
#'   \item{prevalence_est}{Estimated true prevalence at each timepoint.}
#'   \item{method}{Model/method used.}
#'
#' @references
#' MMS-SD workshop Module 4 (power/sample-size formulae for trend detection),
#' Module 5 (ICC / design effect).
#'
#' @export
test_trend <- function(
  positives,
  n,
  timepoints  = NULL,
  sensitivity = 1,
  specificity = 1,
  alternative = c("two.sided", "increasing", "decreasing"),
  alpha       = 0.05,
  conf_level  = 0.95,
  n_per_site  = NULL,
  n_sites     = NULL,
  icc         = 0,
  fpc_N       = NULL
) {
  stop("test_trend() is not yet implemented")
}
