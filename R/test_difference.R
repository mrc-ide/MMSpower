#' Test whether prevalence differs between two groups
#'
#' @description
#' Analysis function: given observed positive counts in two groups, tests
#' whether their true prevalence differs. Returns a z-statistic, p-value,
#' reject / fail-to-reject decision, the estimated prevalence difference, and
#' a confidence interval on that difference.
#'
#' The test is run on the apparent-prevalence scale (observed positives),
#' which -- when sensitivity and specificity are equal across groups -- is
#' equivalent to testing the true-prevalence difference. Point estimates and
#' the confidence interval are transformed back to the true scale via the
#' Rogan-Gladen correction. Clustering is handled per group with a Kish
#' design effect, mirroring \code{test_threshold()}.
#'
#' This is the analysis complement to \code{design_difference()}.
#'
#' \strong{Hypothesis} (\eqn{p_1}, \eqn{p_2} are the true prevalences):
#' \describe{
#'   \item{\code{"two.sided"} (default)}{H0: \eqn{p_1 = p_2} vs
#'     Ha: \eqn{p_1 \ne p_2}.}
#'   \item{\code{"greater"}}{H0: \eqn{p_1 \le p_2} vs Ha: \eqn{p_1 > p_2}.}
#'   \item{\code{"less"}}{H0: \eqn{p_1 \ge p_2} vs Ha: \eqn{p_1 < p_2}.}
#' }
#'
#' @param x1 Integer vector of positive counts per cluster/site in group 1
#'   (a single number if group 1 is one site).
#' @param n1 Integer vector of total samples per cluster/site in group 1
#'   (same length as \code{x1}).
#' @param x2,n2 As \code{x1} / \code{n1}, for group 2.
#' @param alternative Direction of the alternative hypothesis: \code{"two.sided"}
#'   (default), \code{"greater"} (Ha: group 1 > group 2), or \code{"less"}.
#' @param sensitivity Diagnostic sensitivity in (0, 1]; default 1. Assumed
#'   equal across groups.
#' @param specificity Diagnostic specificity in (0, 1]; default 1. Assumed
#'   equal across groups.
#' @param conf_level Confidence level for the reported CI on the difference;
#'   also sets the significance level as alpha = 1 - conf_level. Default 0.95.
#' @param icc Optional intra-cluster correlation. \code{NULL} (default) =
#'   estimate separately for each group from its own clusters. \code{0} =
#'   force SRS (no clustering adjustment). A positive scalar is applied to
#'   both groups.
#' @param fpc_N Optional finite-population size. \code{NULL} (default) = no
#'   FPC. A single value is used as the population size of \emph{each} group;
#'   a length-2 vector gives \code{c(group1, group2)}.
#' @param ci_method Confidence-interval method for the difference:
#'   \code{"wald"} (default) or \code{"newcombe"} (Newcombe 1998 hybrid-score
#'   interval, recommended for small samples or prevalences near 0/1). This
#'   affects only the reported interval, not the hypothesis test, which
#'   always uses the pooled-variance z-statistic.
#'
#' @details
#' \strong{Test statistic.} With \eqn{\hat{p}_{g}} the apparent prevalence in
#' group \eqn{g} and \eqn{m_g} its effective sample size (total n divided by
#' the design effect, then by the squared FPC factor):
#'
#' \deqn{z = \frac{\hat{p}_1 - \hat{p}_2}
#'   {\sqrt{\bar{p}(1-\bar{p})\,(1/m_1 + 1/m_2)}}, \qquad
#'   \bar{p} = \frac{\hat{p}_1 m_1 + \hat{p}_2 m_2}{m_1 + m_2}.}
#'
#' \strong{Confidence interval.} Built on the apparent scale then divided by
#' \eqn{(Se + Sp - 1)} to return to the true scale (the Rogan-Gladen offset
#' cancels in a difference). \code{"wald"} uses the unpooled standard error;
#' \code{"newcombe"} combines per-group Wilson intervals. Bounds are clamped
#' to \eqn{[-1, 1]}.
#'
#' \strong{Design effect.} Estimated per group exactly as in
#' \code{test_threshold()} / \code{estimate_prevalence()}: from the between-
#' vs within-cluster variance when \code{icc = NULL}, or as
#' \eqn{1 + (\bar{n} - 1)\,ICC} for a supplied \code{icc}. A group with fewer
#' than two clusters (or clusters all of size 1) falls back to
#' \eqn{D_{eff} = 1}.
#'
#' @section Equations and sources:
#' Building blocks from the MMS-SD Study Design Workshop
#' (\url{https://mrc-ide.github.io/MMS-SD_workshop/}); the two-group test
#' itself is \strong{not} covered there.
#' \itemize{
#'   \item \emph{One-sample z-test for a proportion}
#'     \eqn{Z = (\hat p - p_0)/\sqrt{\hat p(1-\hat p)/n}} -- Module 3
#'     "Hypothesis testing", slide "Null hypothesis testing" (lecture
#'     slides p. 8). This function uses the \emph{two-sample}
#'     (pooled-variance) form, which the workshop does not cover; standard,
#'     e.g. Fleiss, Levin & Paik (2003), ch. 2.
#'   \item \emph{Design effect} \eqn{D_{eff} = 1 + (\bar n - 1)\,r} and the
#'     \eqn{D_{eff} = \mathrm{Var}_{obs}/\mathrm{Var}_{SRS}} estimator --
#'     Module 5, slides "The Design Effect" / "Why is the ICC useful?"
#'     (pp. 4-5); applied per group.
#'   \item \emph{Wilson score interval} (used by \code{ci_method =
#'     "newcombe"}) and the \emph{Newcombe} hybrid-score difference
#'     interval -- not in the workshop; Newcombe (1998).
#'   \item \emph{Rogan-Gladen correction} of the difference and its CI
#'     (divide by \eqn{Se + Sp - 1}) -- not in the workshop;
#'     Rogan & Gladen (1978).
#'   \item \emph{Finite-population correction} -- not in the workshop;
#'     Cochran (1977).
#' }
#'
#' @return A named list:
#'   \item{difference}{Rogan-Gladen corrected prevalence difference
#'     (group 1 - group 2) on the true scale, computed as
#'     \code{difference_app / (Se + Sp - 1)} and clamped to [-1, 1]. This
#'     is the quantity the confidence interval is centred on; it can differ
#'     slightly from \code{prevalence1 - prevalence2} when a per-group
#'     estimate is clamped.}
#'   \item{difference_app}{Apparent-scale difference
#'     (\eqn{\hat{p}_1 - \hat{p}_2}).}
#'   \item{statistic}{Pooled-variance z-statistic (apparent scale).}
#'   \item{p_value}{p-value for the chosen \code{alternative}.}
#'   \item{reject}{Logical: \code{TRUE} if \code{p_value < 1 - conf_level}.}
#'   \item{alternative}{Alternative hypothesis as supplied.}
#'   \item{prevalence1, prevalence2}{Per-group Rogan-Gladen corrected point
#'     estimates, each clamped to [0, 1].}
#'   \item{apparent_prev1, apparent_prev2}{Per-group apparent prevalences.}
#'   \item{ci_lower, ci_upper}{Confidence interval on the true-scale
#'     difference at \code{conf_level} (two-sided regardless of
#'     \code{alternative}).}
#'   \item{ci_method}{CI method used (as supplied).}
#'   \item{n_total1, n_total2}{Total samples per group.}
#'   \item{n_eff1, n_eff2}{Per-group effective sample size before the FPC
#'     (\code{n_total / deff}).}
#'   \item{n_eff_adj1, n_eff_adj2}{Per-group effective sample size the test
#'     and CI use (\code{n_eff} divided by the squared FPC factor; equals
#'     \code{n_eff} when \code{fpc_N} is \code{NULL}).}
#'   \item{conf_level}{Confidence level as supplied.}
#'   \item{sensitivity, specificity}{As supplied.}
#'   \item{icc_used1, icc_used2}{Per-group ICC applied (estimated or supplied).}
#'   \item{deff1, deff2}{Per-group design effect applied.}
#'   \item{fpc_N}{\code{fpc_N} as supplied, or \code{NULL}.}
#'
#' @references
#' MMS-SD Study Design Workshop, Modules 3 (hypothesis testing) and 5 (ICC /
#' design effect). \url{https://mrc-ide.github.io/MMS-SD_workshop/}
#'
#' Newcombe RG (1998). Interval estimation for the difference between
#' independent proportions: comparison of eleven methods. Statistics in
#' Medicine 17(8):873-890.
#'
#' Fleiss JL, Levin B, Paik MC (2003). Statistical Methods for Rates and
#' Proportions, 3rd ed. Wiley. (Two-sample z-test for proportions.)
#'
#' Rogan WJ, Gladen B (1978). Estimating prevalence from the results of a
#' screening test. American Journal of Epidemiology 107(1):71-76.
#'
#' Cochran WG (1977). Sampling Techniques, 3rd ed. Wiley.
#'
#' @seealso \code{\link{design_difference}}
#'
#' @export
#'
#' @examples
#' # Two single sites: 8/100 vs 20/120
#' test_difference(x1 = 8, n1 = 100, x2 = 20, n2 = 120)
#'
#' # Clustered data, one-sided, imperfect test
#' test_difference(
#'   x1 = c(2, 5, 3), n1 = c(50, 60, 40),
#'   x2 = c(8, 11, 9), n2 = c(55, 70, 45),
#'   alternative = "less",
#'   sensitivity = 0.9, specificity = 0.97
#' )
#'
#' # Newcombe interval on the difference
#' test_difference(x1 = 3, n1 = 80, x2 = 12, n2 = 90, ci_method = "newcombe")
test_difference <- function(x1, n1, x2, n2,
                            alternative = "two.sided",
                            sensitivity = 1,
                            specificity = 1,
                            conf_level  = 0.95,
                            icc         = NULL,
                            fpc_N       = NULL,
                            ci_method   = "wald") {

  # ---- validate the four count vectors ----
  check_counts <- function(x, n, xnm, nnm) {
    if (is.logical(x) || is.logical(n))
      stop("`", xnm, "` and `", nnm, "` must be numeric, not logical (got class `",
           class(x)[1], "` for ", xnm, ", `", class(n)[1], "` for ", nnm, "). ",
           "Note: plain `NA` is logical in R -- filter missing observations first.")
    if (!is.numeric(x) || !is.numeric(n))
      stop("`", xnm, "` and `", nnm, "` must be numeric vectors (got class `",
           class(x)[1], "` for ", xnm, ", `", class(n)[1], "` for ", nnm, ").")
    if (length(x) == 0 || length(n) == 0)
      stop("`", xnm, "` and `", nnm, "` must be non-empty vectors.")
    if (length(x) != length(n))
      stop("`", xnm, "` and `", nnm, "` must have the same length ",
           "(got ", length(x), " and ", length(n), ").")
    if (!all(is.finite(x)) || !all(is.finite(n)))
      stop("`", xnm, "` / `", nnm, "` contain NA, NaN, or infinite values. ",
           "All counts must be finite.")
    bad <- which(x != floor(x))
    if (length(bad))
      stop("`", xnm, "` must contain whole numbers (found ", xnm, "[",
           bad[1], "] = ", x[bad[1]], ").")
    bad <- which(n != floor(n))
    if (length(bad))
      stop("`", nnm, "` must contain whole numbers (found ", nnm, "[",
           bad[1], "] = ", n[bad[1]], ").")
    bad <- which(n <= 0)
    if (length(bad))
      stop("`", nnm, "` must be positive for every cluster (found ", nnm, "[",
           bad[1], "] = ", n[bad[1]], ").")
    bad <- which(x < 0)
    if (length(bad))
      stop("`", xnm, "` must be non-negative (found ", xnm, "[",
           bad[1], "] = ", x[bad[1]], ").")
    bad <- which(x > n)
    if (length(bad))
      stop("`", xnm, "` cannot exceed `", nnm, "` (found ", xnm, "[",
           bad[1], "] = ", x[bad[1]], " > ", nnm, "[",
           bad[1], "] = ", n[bad[1]], ").")
    invisible(TRUE)
  }
  check_counts(x1, n1, "x1", "n1")
  check_counts(x2, n2, "x2", "n2")

  # ---- validate scalar parameters ----
  if (!is.character(alternative) || length(alternative) != 1)
    stop("`alternative` must be a single string: 'greater', 'less', or 'two.sided'.")
  if (!alternative %in% c("greater", "less", "two.sided"))
    stop("`alternative` must be 'greater', 'less', or 'two.sided' (got '",
         alternative, "').")
  if (is.null(sensitivity) || length(sensitivity) != 1 || !is.numeric(sensitivity))
    stop("`sensitivity` must be a single number in (0, 1] (got ",
         if (is.null(sensitivity)) "NULL"
         else if (!is.numeric(sensitivity)) paste0("class `", class(sensitivity)[1], "`")
         else paste0("length = ", length(sensitivity)), "). ",
         "Note: `TRUE`/`FALSE` is logical, not numeric -- pass 1 for a perfect test.")
  if (is.null(specificity) || length(specificity) != 1 || !is.numeric(specificity))
    stop("`specificity` must be a single number in (0, 1] (got ",
         if (is.null(specificity)) "NULL"
         else if (!is.numeric(specificity)) paste0("class `", class(specificity)[1], "`")
         else paste0("length = ", length(specificity)), "). ",
         "Note: `TRUE`/`FALSE` is logical, not numeric -- pass 1 for a perfect test.")
  if (!is.finite(sensitivity) || sensitivity <= 0 || sensitivity > 1)
    stop("`sensitivity` must be in (0, 1] (got ", sensitivity, ").")
  if (!is.finite(specificity) || specificity <= 0 || specificity > 1)
    stop("`specificity` must be in (0, 1] (got ", specificity, ").")
  correction <- sensitivity + specificity - 1
  if (correction <= 0)
    stop("sensitivity + specificity must exceed 1 for the Rogan-Gladen correction ",
         "(got ", sensitivity, " + ", specificity, " = ", sensitivity + specificity, "). ",
         "With se + sp <= 1 the test performs at or below chance.")
  if (correction < 0.1)
    warning("sensitivity + specificity = ", round(sensitivity + specificity, 4),
            " is very close to 1. The Rogan-Gladen adjustment is numerically unstable.")
  if (length(conf_level) != 1 || !is.numeric(conf_level))
    stop("`conf_level` must be a single number (got ",
         if (!is.numeric(conf_level)) paste0("class `", class(conf_level)[1], "`")
         else paste0("length = ", length(conf_level)), ").")
  if (!is.finite(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("`conf_level` must be strictly between 0 and 1 (got ", conf_level, ").")
  if (!is.null(icc) && (length(icc) != 1 || !is.numeric(icc)))
    stop("`icc` must be a single number or NULL (got ",
         if (!is.numeric(icc)) paste0("class `", class(icc)[1], "`")
         else paste0("length = ", length(icc)), "). ",
         "Set `icc = NULL` to estimate ICC from the data.")
  if (!is.null(icc) && (!is.finite(icc) || icc < 0 || icc > 1))
    stop("`icc` must be in [0, 1] (got ", icc, "). ",
         "Set `icc = NULL` to estimate ICC from the data.")
  if (!is.null(fpc_N) && (!is.numeric(fpc_N) || !(length(fpc_N) %in% c(1, 2)) ||
      !all(is.finite(fpc_N)) || any(fpc_N <= 0)))
    stop("`fpc_N` must be NULL, a single finite positive number, or a ",
         "length-2 positive vector c(group1, group2) (got ",
         if (!is.numeric(fpc_N)) paste0("class `", class(fpc_N)[1], "`")
         else paste0("length = ", length(fpc_N)), ").")
  if (!is.character(ci_method) || length(ci_method) != 1)
    stop("`ci_method` must be a single character string: 'wald' or 'newcombe' ",
         "(got class `", class(ci_method)[1], "`, length ", length(ci_method), ").")
  if (!ci_method %in% c("wald", "newcombe"))
    stop("`ci_method` must be 'wald' or 'newcombe' (got '", ci_method, "').")

  # ---- per-group design effect (same logic as test_threshold) ----
  group_deff <- function(x, n, icc_arg) {
    n_clusters <- length(n)
    n_bar      <- mean(n)
    p_hat      <- sum(x) / sum(n)

    if (is.null(icc_arg)) {
      if (n_clusters < 2 || n_bar == 1)
        return(list(icc_used = 0, deff = 1))
      p_i     <- x / n
      var_obs <- stats::var(p_i)
      var_srs <- mean(p_hat * (1 - p_hat) / n)
      deff    <- if (var_srs > 0) max(var_obs / var_srs, 1) else 1
      icc_u   <- min(max((deff - 1) / (n_bar - 1), 0), 1)
      list(icc_used = icc_u, deff = 1 + (n_bar - 1) * icc_u)
    } else if (n_clusters < 2 || n_bar == 1) {
      # A supplied icc has no effect with a single cluster (or clusters all
      # of size 1): the Kish denominator is undefined -- fall back to SRS.
      list(icc_used = 0, deff = 1)
    } else {
      list(icc_used = icc_arg, deff = 1 + (n_bar - 1) * icc_arg)
    }
  }

  g1 <- group_deff(x1, n1, icc)
  g2 <- group_deff(x2, n2, icc)

  n_tot1 <- sum(n1); n_tot2 <- sum(n2)
  p_hat1 <- sum(x1) / n_tot1
  p_hat2 <- sum(x2) / n_tot2
  n_eff1 <- n_tot1 / g1$deff
  n_eff2 <- n_tot2 / g2$deff

  # ---- finite-population correction (per group) ----
  fpc_N1 <- fpc_N2 <- NULL
  if (!is.null(fpc_N)) {
    if (length(fpc_N) == 1) { fpc_N1 <- fpc_N; fpc_N2 <- fpc_N }
    else                    { fpc_N1 <- fpc_N[1]; fpc_N2 <- fpc_N[2] }
  }
  fpc_check <- function(fN, n_tot, n_eff, grp) {
    if (is.null(fN)) return(invisible())
    if (fN < n_tot)
      stop("`fpc_N` for group ", grp, " (", fN, ") is less than that group's ",
           "sample size (", n_tot, "). The population must be at least the sample.")
    if (fN <= n_eff)
      stop("`fpc_N` for group ", grp, " (", fN, ") is <= its effective sample ",
           "size n_eff (", round(n_eff, 2), "); the variance collapses to zero. ",
           "Set `fpc_N = NULL` if no FPC is needed.")
  }
  fpc_check(fpc_N1, n_tot1, n_eff1, 1)
  fpc_check(fpc_N2, n_tot2, n_eff2, 2)

  fpc1 <- if (!is.null(fpc_N1)) sqrt((fpc_N1 - n_eff1) / (fpc_N1 - 1)) else 1
  fpc2 <- if (!is.null(fpc_N2)) sqrt((fpc_N2 - n_eff2) / (fpc_N2 - 1)) else 1
  m1   <- n_eff1 / fpc1^2     # variance-equivalent simple sample size
  m2   <- n_eff2 / fpc2^2

  # ---- pooled-variance z-test on the apparent scale ----
  d_app  <- p_hat1 - p_hat2
  p_pool <- (p_hat1 * m1 + p_hat2 * m2) / (m1 + m2)
  se_null <- sqrt(p_pool * (1 - p_pool) * (1 / m1 + 1 / m2))

  if (se_null == 0)
    stop("Pooled standard error is zero (both groups have apparent prevalence ",
         p_pool, "). There is no variability to test.")

  z_stat <- d_app / se_null

  p_value <- switch(alternative,
    greater   = stats::pnorm(z_stat, lower.tail = FALSE),
    less      = stats::pnorm(z_stat, lower.tail = TRUE),
    two.sided = 2 * stats::pnorm(abs(z_stat), lower.tail = FALSE)
  )

  alpha  <- 1 - conf_level
  reject <- p_value < alpha

  # ---- two-sided CI on the difference (apparent scale, then /correction) ----
  z_ci <- stats::qnorm(1 - alpha / 2)

  if (ci_method == "wald") {
    se_unpooled <- sqrt(p_hat1 * (1 - p_hat1) / m1 + p_hat2 * (1 - p_hat2) / m2)
    lo_app <- d_app - z_ci * se_unpooled
    hi_app <- d_app + z_ci * se_unpooled
  } else {   # newcombe (hybrid score)
    w1 <- .wilson_ci(p_hat1, m1, z_ci)
    w2 <- .wilson_ci(p_hat2, m2, z_ci)
    lo_app <- d_app - sqrt((p_hat1 - w1$lower)^2 + (w2$upper - p_hat2)^2)
    hi_app <- d_app + sqrt((w1$upper - p_hat1)^2 + (p_hat2 - w2$lower)^2)
  }

  ci_lower <- max(-1, min(1, lo_app / correction))
  ci_upper <- max(-1, min(1,  hi_app / correction))

  # ---- point estimates on the true scale ----
  # `difference` is the Rogan-Gladen corrected difference d_app / (Se+Sp-1)
  # (the offset cancels), which is exactly what the CI above is centred on.
  # Deriving it as prevalence1 - prevalence2 instead would diverge from the
  # CI whenever a per-group estimate is clamped to [0, 1] -- common for a
  # low true prevalence with imperfect specificity -- and could place the
  # point estimate outside its own interval.
  rg <- function(p) .rogan_gladen(p, sensitivity, specificity)
  prevalence1 <- max(0, min(1, rg(p_hat1)))
  prevalence2 <- max(0, min(1, rg(p_hat2)))
  difference  <- max(-1, min(1, d_app / correction))

  list(
    difference     = difference,
    difference_app = d_app,
    statistic      = z_stat,
    p_value        = p_value,
    reject         = reject,
    alternative    = alternative,
    prevalence1    = prevalence1,
    prevalence2    = prevalence2,
    apparent_prev1 = p_hat1,
    apparent_prev2 = p_hat2,
    ci_lower       = ci_lower,
    ci_upper       = ci_upper,
    ci_method      = ci_method,
    n_total1       = n_tot1,
    n_total2       = n_tot2,
    n_eff1         = n_eff1,
    n_eff2         = n_eff2,
    n_eff_adj1     = m1,
    n_eff_adj2     = m2,
    conf_level     = conf_level,
    sensitivity    = sensitivity,
    specificity    = specificity,
    icc_used1      = g1$icc_used,
    icc_used2      = g2$icc_used,
    deff1          = g1$deff,
    deff2          = g2$deff,
    fpc_N          = fpc_N
  )
}
