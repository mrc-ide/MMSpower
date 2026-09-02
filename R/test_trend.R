#' Test whether prevalence is changing linearly over time
#'
#' @description
#' Analysis function: given observed positive counts at several timepoints
#' (optionally with multiple clusters per timepoint), fits a weighted
#' least-squares line to the apparent prevalences and tests whether the
#' slope differs from zero. Returns the slope estimate on the true-prevalence
#' scale, a confidence interval, a z-statistic, a p-value, and a
#' reject / fail-to-reject decision.
#'
#' The fit is done on the apparent scale with inverse-variance (binomial)
#' weights; the slope and its CI are divided by \eqn{(Se + Sp - 1)} to return
#' to the true scale. Clustering is folded in as a single design-effect
#' inflation of the weights, mirroring \code{test_threshold()} /
#' \code{test_difference()}.
#'
#' This is the analysis complement to \code{design_trend()}.
#'
#' \strong{Hypothesis} (\eqn{\beta} is the true-prevalence slope):
#' \describe{
#'   \item{\code{"two.sided"} (default)}{H0: \eqn{\beta = 0} vs
#'     Ha: \eqn{\beta \ne 0}.}
#'   \item{\code{"greater"}}{H0: \eqn{\beta \le 0} vs Ha: \eqn{\beta > 0}
#'     (an increasing trend).}
#'   \item{\code{"less"}}{H0: \eqn{\beta \ge 0} vs Ha: \eqn{\beta < 0}.}
#' }
#'
#' @param x Integer vector of positive counts. One entry per timepoint, or
#'   one entry per cluster-by-timepoint (in which case give repeated
#'   \code{time} values).
#' @param n Integer vector of sample sizes, same length as \code{x}.
#' @param time Numeric vector of timepoint coordinates, same length as
#'   \code{x} (e.g. calendar years or survey rounds). If \code{NULL}
#'   (default), \code{seq_along(x)} is used -- only sensible when there is
#'   one entry per timepoint.
#' @param alternative Direction of the alternative hypothesis:
#'   \code{"two.sided"} (default), \code{"greater"} (increasing), or
#'   \code{"less"} (decreasing).
#' @param sensitivity Diagnostic sensitivity in (0, 1]; default 1.
#' @param specificity Diagnostic specificity in (0, 1]; default 1.
#' @param conf_level Confidence level for the CI on the slope; also sets the
#'   significance level as alpha = 1 - conf_level. Default 0.95.
#' @param icc Optional intra-cluster correlation. \code{NULL} (default) =
#'   estimate from the data when any timepoint has two or more clusters
#'   (rows), otherwise treat as 0. \code{0} = force SRS. A positive scalar
#'   is used directly \emph{only} when the data has a cluster structure
#'   (some timepoint has >= 2 rows); with one row per timepoint there is no
#'   cluster size to form a design effect from, so it is ignored and
#'   \code{deff = 1} (as in \code{test_threshold()}).
#' @param fpc_N Optional finite-population size \strong{per timepoint}.
#'   \code{NULL} (default) = no FPC.
#' @param ci_method Reserved for future methods; currently only
#'   \code{"wls"} (weighted least squares) is implemented.
#'
#' @details
#' \strong{Model.} Counts are aggregated to unique timepoints
#' \eqn{(t_j, X_j, N_j)}, giving apparent prevalences
#' \eqn{\hat p_j = X_j / N_j} with variance \eqn{\hat p_j (1 - \hat p_j)/N_j}
#' (a \eqn{(X_j + 0.5)/(N_j + 1)} plug-in is used for the weight when
#' \eqn{X_j} is 0 or \eqn{N_j}). The weighted least-squares slope and its
#' variance are
#'
#' \deqn{\hat\beta_{app} = \frac{\sum_j w_j (t_j - \bar t_w)\,\hat p_j}
#'   {\sum_j w_j (t_j - \bar t_w)^2}, \qquad
#'   \mathrm{Var}(\hat\beta_{app}) = \frac{1}{\sum_j w_j (t_j - \bar t_w)^2},}
#'
#' with weights \eqn{w_j = 1 / \mathrm{Var}(\hat p_j)} divided by the design
#' effect and, if requested, by the squared FPC factor
#' \eqn{(N_{pop} - N_j)/(N_{pop} - 1)}. The z-statistic is
#' \eqn{\hat\beta_{app} / \sqrt{\mathrm{Var}(\hat\beta_{app})}}; the true-scale
#' slope and CI follow by dividing by \eqn{(Se + Sp - 1)}. The CI is
#' matched to \code{alternative}: one-sided \eqn{[L, \infty)} for
#' \code{"greater"}, \eqn{(-\infty, U]} for \code{"less"}, two-sided for
#' \code{"two.sided"}, so \code{reject} agrees with whether 0 is outside it.
#'
#' \strong{Design effect.} With \code{icc = NULL}, a Kish design effect is
#' estimated from the between-cluster spread within each timepoint that has
#' at least two clusters, averaged, and applied as
#' \eqn{1 + (\bar n - 1)\,\widehat{ICC}} where \eqn{\bar n} is the mean
#' cluster size. This is a single uniform inflation -- a full mixed-model
#' treatment of cluster-by-time structure is future work.
#'
#' @section Equations and sources:
#' Building blocks from the MMS-SD Study Design Workshop
#' (\url{https://mrc-ide.github.io/MMS-SD_workshop/}); the trend test
#' itself is \strong{not} covered there.
#' \itemize{
#'   \item The workshop poses "Has prevalence increased over the last 5
#'     years?" as a hypothesis-test question (Module 3 "Hypothesis
#'     testing", slide "Null hypothesis testing", lecture slides p. 8) but
#'     gives no method. This function fits a weighted least-squares line;
#'     the classical count-based alternative is the Cochran-Armitage trend
#'     test (Armitage 1955; Cochran 1954).
#'   \item \emph{Inverse-variance weights} \eqn{w_j = 1/\mathrm{Var}(\hat
#'     p_j)} with \eqn{\mathrm{Var}(\hat p_j) = \hat p_j(1-\hat p_j)/N_j}
#'     -- the binomial variance of a proportion, as in Module 1's Wald
#'     interval (lecture slides p. 11).
#'   \item \emph{Design effect} \eqn{D_{eff} = 1 + (\bar n - 1)\,r} and the
#'     \eqn{\mathrm{Var}_{obs}/\mathrm{Var}_{SRS}} estimator -- Module 5,
#'     slides "The Design Effect" / "Why is the ICC useful?" (pp. 4-5).
#'   \item \emph{Rogan-Gladen correction} of the slope and its CI (divide
#'     by \eqn{Se + Sp - 1}) -- not in the workshop; Rogan & Gladen (1978).
#'   \item \emph{Finite-population correction} -- not in the workshop;
#'     Cochran (1977).
#' }
#'
#' @return A named list:
#'   \item{slope}{Rogan-Gladen corrected slope (true prevalence change per
#'     unit time).}
#'   \item{slope_app}{Apparent-scale slope.}
#'   \item{statistic}{z-statistic for the slope.}
#'   \item{p_value}{p-value for the chosen \code{alternative}.}
#'   \item{reject}{Logical: \code{TRUE} if \code{p_value < 1 - conf_level}.}
#'   \item{alternative}{As supplied.}
#'   \item{ci_lower, ci_upper}{Confidence interval on the true-scale slope
#'     at \code{conf_level}, matched to \code{alternative}: \code{ci_upper}
#'     is \code{Inf} for \code{"greater"}, \code{ci_lower} is \code{-Inf}
#'     for \code{"less"}. So \code{reject} agrees with whether 0 lies
#'     outside the interval.}
#'   \item{prevalence_start_est, prevalence_end_est}{Fitted true prevalence
#'     at the earliest and latest observed timepoint, clamped to [0, 1].}
#'   \item{times}{Sorted unique timepoints.}
#'   \item{apparent_prev}{Apparent prevalence at each unique timepoint.}
#'   \item{x_by_time, n_by_time}{Pooled positives and totals per timepoint.}
#'   \item{n_timepoints}{Number of unique timepoints.}
#'   \item{n_total}{Total samples.}
#'   \item{conf_level}{As supplied.}
#'   \item{sensitivity, specificity}{As supplied.}
#'   \item{icc_used}{ICC applied (estimated or supplied).}
#'   \item{deff}{Design effect applied to the weights.}
#'   \item{fpc_N}{As supplied, or \code{NULL}.}
#'   \item{method}{\code{"wls"}.}
#'
#' @references
#' MMS-SD Study Design Workshop, Modules 1 (Wald interval), 3 (hypothesis
#' testing) and 5 (ICC / design effect).
#' \url{https://mrc-ide.github.io/MMS-SD_workshop/}
#'
#' Armitage P (1955). Tests for linear trends in proportions and
#' frequencies. Biometrics 11(3):375-386.
#'
#' Cochran WG (1954). Some methods for strengthening the common chi-square
#' tests. Biometrics 10(4):417-451.
#'
#' Rogan WJ, Gladen B (1978). Estimating prevalence from the results of a
#' screening test. American Journal of Epidemiology 107(1):71-76.
#'
#' Cochran WG (1977). Sampling Techniques, 3rd ed. Wiley.
#'
#' @seealso \code{\link{design_trend}}
#'
#' @export
#'
#' @examples
#' # Five annual rounds, one survey per round
#' test_trend(
#'   x    = c(5, 8, 9, 14, 17),
#'   n    = c(100, 100, 100, 100, 100),
#'   time = 2019:2023
#' )
#'
#' # Clustered: several sites per round
#' test_trend(
#'   x    = c(2, 3, 4, 5, 6, 8, 7, 10),
#'   n    = c(50, 45, 55, 50, 48, 52, 47, 51),
#'   time = c(1, 1, 2, 2, 3, 3, 4, 4)
#' )
#'
#' # Imperfect test, one-sided (decline)
#' test_trend(x = c(20, 15, 11, 7), n = rep(120, 4), time = 1:4,
#'            alternative = "less", sensitivity = 0.9, specificity = 0.97)
test_trend <- function(x, n, time = NULL,
                       alternative = "two.sided",
                       sensitivity = 1,
                       specificity = 1,
                       conf_level  = 0.95,
                       icc         = NULL,
                       fpc_N       = NULL,
                       ci_method   = "wls") {

  # ---- validate x / n ----
  if (is.logical(x) || is.logical(n))
    stop("`x` and `n` must be numeric, not logical (got class `", class(x)[1],
         "` for x, `", class(n)[1], "` for n). ",
         "Note: plain `NA` is logical in R -- filter missing observations first.")
  if (!is.numeric(x) || !is.numeric(n))
    stop("`x` and `n` must be numeric vectors (got class `", class(x)[1],
         "` for x, `", class(n)[1], "` for n).")
  if (length(x) == 0 || length(x) != length(n))
    stop("`x` and `n` must be non-empty vectors of the same length ",
         "(got ", length(x), " and ", length(n), ").")
  if (!all(is.finite(x)) || !all(is.finite(n)))
    stop("`x` / `n` contain NA, NaN, or infinite values. All counts must be finite.")
  if (any(x != floor(x)))
    stop("`x` must contain whole numbers (found x[",
         which(x != floor(x))[1], "] = ", x[which(x != floor(x))[1]], ").")
  if (any(n != floor(n)))
    stop("`n` must contain whole numbers (found n[",
         which(n != floor(n))[1], "] = ", n[which(n != floor(n))[1]], ").")
  if (any(n <= 0))
    stop("`n` must be positive for every row (found n[",
         which(n <= 0)[1], "] = ", n[which(n <= 0)[1]], ").")
  if (any(x < 0))
    stop("`x` must be non-negative (found x[",
         which(x < 0)[1], "] = ", x[which(x < 0)[1]], ").")
  if (any(x > n))
    stop("`x` cannot exceed `n` (found x[", which(x > n)[1], "] = ",
         x[which(x > n)[1]], " > n[", which(x > n)[1], "] = ",
         n[which(x > n)[1]], ").")

  # ---- validate time ----
  if (is.null(time)) {
    time <- seq_along(x)
  } else {
    if (!is.numeric(time) || length(time) != length(x) || !all(is.finite(time)))
      stop("`time` must be a finite numeric vector the same length as `x` ",
           "(got length ", length(time), " vs ", length(x), ").")
  }
  if (length(unique(time)) < 2)
    stop("`test_trend()` needs at least two distinct `time` values to fit a slope ",
         "(got ", length(unique(time)), ").")

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
         "(got ", sensitivity, " + ", specificity, " = ", sensitivity + specificity, ").")
  if (correction < 0.1)
    warning("sensitivity + specificity = ", round(sensitivity + specificity, 4),
            " is very close to 1. The Rogan-Gladen adjustment is numerically unstable.")
  if (length(conf_level) != 1 || !is.numeric(conf_level) ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("`conf_level` must be a single number strictly between 0 and 1 (got ",
         conf_level, ").")
  if (!is.null(icc) && (length(icc) != 1 || !is.numeric(icc) ||
      !is.finite(icc) || icc < 0 || icc > 1))
    stop("`icc` must be NULL or a single number in [0, 1] (got ",
         if (!is.numeric(icc)) paste0("class `", class(icc)[1], "`") else icc, "). ",
         "Set `icc = NULL` to estimate it from the data.")
  if (!is.null(fpc_N) && (!is.numeric(fpc_N) || length(fpc_N) != 1 ||
      !is.finite(fpc_N) || fpc_N <= 0))
    stop("`fpc_N` must be NULL or a single finite positive number (got ",
         if (length(fpc_N) != 1) paste0("length = ", length(fpc_N)) else fpc_N, ").")
  if (!is.character(ci_method) || length(ci_method) != 1 || ci_method != "wls")
    stop("`ci_method` currently supports only 'wls' (got '", ci_method, "').")

  # ---- aggregate to unique timepoints ----
  ord   <- order(time)
  time  <- time[ord]; x <- x[ord]; n <- n[ord]
  tj    <- sort(unique(time))
  Xj    <- vapply(tj, function(t) sum(x[time == t]), numeric(1))
  Nj    <- vapply(tj, function(t) sum(n[time == t]), numeric(1))
  pj    <- Xj / Nj

  # ---- design effect (per-timepoint Kish, averaged) ----
  # A "cluster" is a row; a timepoint only carries cluster structure when it
  # has >= 2 rows. With one row per timepoint the row n is the timepoint
  # total, not a cluster size, so no design effect can be formed -- fall back
  # to deff = 1 even when an explicit `icc` is supplied (matches the
  # single-cluster convention in test_threshold() / estimate_prevalence()).
  multi_tp    <- tj[vapply(tj, function(t) sum(time == t) >= 2, logical(1))]
  has_clusters <- length(multi_tp) > 0
  n_bar_cl    <- if (has_clusters) mean(n[time %in% multi_tp]) else NA_real_

  if (is.null(icc)) {
    icc_tp <- c()
    for (t in multi_tp) {
      xi <- x[time == t]; ni <- n[time == t]
      if (mean(ni) == 1) next
      p_pool  <- sum(xi) / sum(ni)
      var_obs <- stats::var(xi / ni)
      var_srs <- mean(p_pool * (1 - p_pool) / ni)
      d_tp    <- if (var_srs > 0) max(var_obs / var_srs, 1) else 1
      icc_tp  <- c(icc_tp, min(max((d_tp - 1) / (mean(ni) - 1), 0), 1))
    }
    icc_used <- if (length(icc_tp) > 0) mean(icc_tp) else 0
  } else {
    icc_used <- icc
  }
  deff <- if (icc_used > 0 && has_clusters && n_bar_cl > 1)
    1 + (n_bar_cl - 1) * icc_used else 1

  # ---- FPC per timepoint ----
  if (!is.null(fpc_N)) {
    if (fpc_N < max(Nj))
      stop("`fpc_N` (", fpc_N, ") is less than the largest per-timepoint sample ",
           "size (", max(Nj), "). The population must be at least the sample.")
    if (any(fpc_N <= Nj / deff))
      stop("`fpc_N` (", fpc_N, ") is <= a timepoint's effective sample size; ",
           "the variance collapses to zero. Set `fpc_N = NULL` if no FPC is needed.")
    fpc2 <- (fpc_N - Nj / deff) / (fpc_N - 1)
  } else {
    fpc2 <- rep(1, length(tj))
  }

  # ---- inverse-variance weights on the apparent scale ----
  # plug-in p for the weight when a timepoint is all-negative or all-positive
  p_w   <- ifelse(Xj == 0 | Xj == Nj, (Xj + 0.5) / (Nj + 1), pj)
  var_j <- p_w * (1 - p_w) / Nj * deff * fpc2
  wj    <- 1 / var_j

  # ---- weighted least squares slope ----
  tw   <- sum(wj * tj) / sum(wj)
  pwm  <- sum(wj * pj) / sum(wj)
  sxx_w <- sum(wj * (tj - tw)^2)

  beta_app <- sum(wj * (tj - tw) * (pj - pwm)) / sxx_w
  var_beta <- 1 / sxx_w
  se_beta  <- sqrt(var_beta)
  z_stat   <- beta_app / se_beta

  p_value <- switch(alternative,
    greater   = stats::pnorm(z_stat, lower.tail = FALSE),
    less      = stats::pnorm(z_stat, lower.tail = TRUE),
    two.sided = 2 * stats::pnorm(abs(z_stat), lower.tail = FALSE)
  )
  alpha  <- 1 - conf_level
  reject <- p_value < alpha

  slope_true <- beta_app / correction

  # CI on the slope, matched to `alternative` (like the tests above):
  # "greater" -> [L, Inf), "less" -> (-Inf, U], "two.sided" -> [L, U].
  # The one-sided bound uses z_{1-alpha}; `reject` then agrees with
  # whether 0 lies outside the interval.
  alpha_lo <- switch(alternative,
    two.sided = alpha / 2, greater = alpha, less = 0)
  alpha_hi <- switch(alternative,
    two.sided = alpha / 2, greater = 0,     less = alpha)
  ci_lower <- if (alpha_lo > 0)
    (beta_app - stats::qnorm(1 - alpha_lo) * se_beta) / correction else -Inf
  ci_upper <- if (alpha_hi > 0)
    (beta_app + stats::qnorm(1 - alpha_hi) * se_beta) / correction else Inf

  # ---- fitted endpoints on the true scale ----
  fit_app_start <- pwm + beta_app * (min(tj) - tw)
  fit_app_end   <- pwm + beta_app * (max(tj) - tw)
  rg <- function(p) .rogan_gladen(p, sensitivity, specificity)
  prevalence_start_est <- max(0, min(1, rg(fit_app_start)))
  prevalence_end_est   <- max(0, min(1, rg(fit_app_end)))

  list(
    slope                = slope_true,
    slope_app            = beta_app,
    statistic            = z_stat,
    p_value              = p_value,
    reject               = reject,
    alternative          = alternative,
    ci_lower             = ci_lower,
    ci_upper             = ci_upper,
    prevalence_start_est = prevalence_start_est,
    prevalence_end_est   = prevalence_end_est,
    times                = tj,
    apparent_prev        = pj,
    x_by_time            = Xj,
    n_by_time            = Nj,
    n_timepoints         = length(tj),
    n_total              = sum(Nj),
    conf_level           = conf_level,
    sensitivity          = sensitivity,
    specificity          = specificity,
    icc_used             = icc_used,
    deff                 = deff,
    fpc_N                = fpc_N,
    method               = "wls"
  )
}
