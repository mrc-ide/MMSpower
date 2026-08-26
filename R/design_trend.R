#' Calculate sample size per timepoint to detect a trend over time
#'
#' Returns the per-timepoint sample size needed to detect a specified trend
#' (rate of change) in true prevalence across a given number of rounds or
#' years with given power.
#'
#' Pairs with `test_trend()`.  Corresponds to longitudinal power/sample-size
#' formulae (M4) with ICC adjustment (M5).
#'
#' @param baseline_prevalence Numeric in (0, 1). True prevalence at the first
#'   timepoint.
#' @param trend Numeric. Expected change in true prevalence per unit time
#'   (e.g., -0.05 means 5 percentage-point decline per year).
#' @param n_timepoints Integer >= 2. Number of equally-spaced timepoints
#'   (rounds, years).
#' @param power Numeric in (0, 1). Target statistical power; default 0.80.
#' @param alpha Numeric. Significance level; default 0.05.
#' @param sensitivity Numeric in (0, 1]. Diagnostic sensitivity; default 1.
#' @param specificity Numeric in (0, 1]. Diagnostic specificity; default 1.
#' @param alternative Character. `"two.sided"` (default), `"increasing"`, or
#'   `"decreasing"`.
#' @param n_per_site Optional integer. Samples per site per timepoint (design
#'   effect).
#' @param n_sites Optional integer. Number of sites.
#' @param icc Numeric >= 0. Intraclass correlation; default 0.
#' @param fpc_N Optional integer. Finite population size.
#'
#' @return An object of class `"mms_design"` (planned) with:
#'   \item{n_per_timepoint}{Required samples at each timepoint.}
#'   \item{n_total}{Total samples across all timepoints.}
#'   \item{n_base_per_timepoint}{Per-timepoint n ignoring clustering and FPC.}
#'   \item{power}{Target power.}
#'   \item{baseline_prevalence}{Baseline assumed.}
#'   \item{trend}{Trend per unit time assumed.}
#'   \item{n_timepoints}{Timepoints used.}
#'   \item{design_effect}{DEFF applied.}
#'
#' @references
#' MMS-SD workshop Module 4 (longitudinal power/sample-size formulae),
#' Module 5 (ICC / design effect).
#'
#' @export
design_trend <- function(
  baseline_prevalence,
  trend,
  n_timepoints = 3,
  power        = 0.80,
  alpha        = 0.05,
  sensitivity  = 1,
  specificity  = 1,
  alternative  = c("two.sided", "increasing", "decreasing"),
  n_per_site   = NULL,
  n_sites      = NULL,
  icc          = 0,
  fpc_N        = NULL
) {
  stop("design_trend() is not yet implemented")
}
