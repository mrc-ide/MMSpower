#' Calculate sample size to detect a difference in prevalence between groups
#'
#' Returns the per-group sample size needed to detect a specified difference in
#' true prevalence between two groups (or sites) with given power.
#'
#' Pairs with `test_difference()`.  Corresponds to two-sample power/sample-size
#' formulae (M3, M4) with ICC adjustment (M5).
#'
#' @param prevalence_1 Numeric in (0, 1). Expected true prevalence in group 1.
#' @param prevalence_2 Numeric in (0, 1). Expected true prevalence in group 2.
#' @param power Numeric in (0, 1). Target statistical power; default 0.80.
#' @param alpha Numeric. Significance level; default 0.05.
#' @param sensitivity Numeric in (0, 1]. Diagnostic sensitivity (assumed equal
#'   across groups); default 1.
#' @param specificity Numeric in (0, 1]. Diagnostic specificity; default 1.
#' @param alternative Character. `"two.sided"` (default), `"greater"`, or
#'   `"less"`.
#' @param n_per_site Optional integer. Samples per site (design effect).
#' @param n_sites Optional integer. Number of sites per group.
#' @param icc Numeric >= 0. Intraclass correlation; default 0.
#' @param fpc_N Optional integer. Finite population size (per group).
#'
#' @return An object of class `"mms_design"` (planned) with:
#'   \item{n_per_group}{Per-group sample size required.}
#'   \item{n_total}{Total sample size (both groups combined).}
#'   \item{n_base_per_group}{Per-group n ignoring clustering and FPC.}
#'   \item{power}{Target power.}
#'   \item{alpha}{Significance level.}
#'   \item{prevalence_1}{Group 1 prevalence assumed.}
#'   \item{prevalence_2}{Group 2 prevalence assumed.}
#'   \item{design_effect}{DEFF applied.}
#'
#' @references
#' MMS-SD workshop Modules 3 (hypothesis testing), 4 (power/sample-size
#' formulae), 5 (ICC / design effect).
#'
#' @export
design_difference <- function(
  prevalence_1,
  prevalence_2,
  power       = 0.80,
  alpha       = 0.05,
  sensitivity = 1,
  specificity = 1,
  alternative = c("two.sided", "greater", "less"),
  n_per_site  = NULL,
  n_sites     = NULL,
  icc         = 0,
  fpc_N       = NULL
) {
  stop("design_difference() is not yet implemented")
}
