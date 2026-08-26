#' Test whether prevalence differs between two groups or sites
#'
#' Given observed positives in each of two groups, tests whether the true
#' prevalence differs between them (two-sided) or is higher in one than the
#' other (one-sided).
#'
#' Corresponds to hypothesis-testing and two-sample proportion tests (M3, M4),
#' extended for imperfect diagnostics and optional clustering (M5).
#'
#' @param positives_1 Integer. Observed positives in group 1.
#' @param n_1 Integer. Total samples in group 1.
#' @param positives_2 Integer. Observed positives in group 2.
#' @param n_2 Integer. Total samples in group 2.
#' @param sensitivity Numeric in (0, 1]. Diagnostic sensitivity (assumed equal
#'   across groups); default 1.
#' @param specificity Numeric in (0, 1]. Diagnostic specificity; default 1.
#' @param alternative Character. `"two.sided"` (default), `"greater"`, or
#'   `"less"`.
#' @param alpha Numeric. Significance level; default 0.05.
#' @param conf_level Numeric in (0, 1). CI level on the difference; default 0.95.
#' @param n_per_site Optional integer. Samples per site (design effect).
#'   Applied equally to both groups.
#' @param n_sites Optional integer. Number of sites per group.
#' @param icc Numeric >= 0. Intraclass correlation; default 0.
#' @param fpc_N Optional integer. Finite population size (per group).
#'
#' @return An object of class `"mms_test_difference"` (planned) with:
#'   \item{reject}{Logical. `TRUE` if difference is significant at `alpha`.}
#'   \item{p_value}{P-value.}
#'   \item{prevalence_1}{Estimated true prevalence in group 1.}
#'   \item{prevalence_2}{Estimated true prevalence in group 2.}
#'   \item{difference}{Estimated true prevalence difference (group 1 - group 2).}
#'   \item{ci_difference}{CI on the difference.}
#'   \item{method}{Test method used.}
#'
#' @references
#' MMS-SD workshop Modules 3 (hypothesis testing), 4 (power/sample-size
#' formulae), 5 (ICC-aware design effect).
#'
#' @export
test_difference <- function(
  positives_1,
  n_1,
  positives_2,
  n_2,
  sensitivity = 1,
  specificity = 1,
  alternative = c("two.sided", "greater", "less"),
  alpha       = 0.05,
  conf_level  = 0.95,
  n_per_site  = NULL,
  n_sites     = NULL,
  icc         = 0,
  fpc_N       = NULL
) {
  stop("test_difference() is not yet implemented")
}
