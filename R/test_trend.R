#' Test whether prevalence is changing linearly over time
#'
#' @description
#' Given the number of positives and the number tested at several
#' timepoints (optionally at several sites per timepoint), fits a straight
#' line to the measured prevalences and tests whether its slope differs
#' from zero with the Cochran-Armitage trend test. Returns the estimated
#' change in true prevalence per unit of time with its confidence interval,
#' the p-value and whether the trend is significant. It is the analysis twin
#' of [design_trend()].
#'
#' **Hypotheses** (\eqn{\beta} is the slope of true prevalence):
#' \describe{
#'   \item{`"two.sided"` (default)}{H0: \eqn{\beta = 0} against
#'     H1: \eqn{\beta \ne 0}.}
#'   \item{`"greater"`}{H0: \eqn{\beta \le 0} against H1: \eqn{\beta > 0}
#'     (an increasing trend).}
#'   \item{`"less"`}{H0: \eqn{\beta \ge 0} against H1: \eqn{\beta < 0}
#'     (a decreasing trend).}
#' }
#'
#' @param x Vector of whole-number counts of people who tested positive: one
#'   per timepoint, or one per site and timepoint (give the timepoint of each
#'   in `time`).
#' @param n Vector of whole-number totals tested, the same length as `x`.
#' @param time Numeric vector giving the timepoint of each value of `x` (for
#'   example calendar years or survey rounds). If `NULL` (default),
#'   `1, 2, ...` is used, which only makes sense with one value per
#'   timepoint.
#' @param alternative Direction of the test: `"two.sided"` (default),
#'   `"greater"` (increasing) or `"less"` (decreasing).
#' @param sensitivity Diagnostic sensitivity in (0, 1]; default 1 (the test
#'   never misses a case).
#' @param specificity Diagnostic specificity in (0, 1]; default 1 (the test
#'   never gives a false positive). `sensitivity + specificity` must exceed 1.
#' @param conf_level Confidence level for the interval, in (0, 1); default
#'   0.95. The significance level of the test is 1 - `conf_level`.
#' @param icc Optional. Intra-cluster correlation, in \[0, 1\]. If `NULL`
#'   (default), it is estimated from the timepoints with two or more sites.
#'   Set to `0` for no clustering.
#'
#' @details
#' **Line.** The counts are added up within each timepoint, giving
#' \eqn{X_j} positives of \eqn{N_j} tested at time \eqn{t_j}, and measured
#' prevalence \eqn{\hat p_j = X_j / N_j}. Each timepoint is weighted by its
#' effective number tested, \eqn{w_j = N_j / D_{eff}}. The slope is
#' \deqn{\hat\beta_{app} = \frac{\sum_j w_j (t_j - \bar t)\,\hat p_j}
#'   {\sum_j w_j (t_j - \bar t)^2},}
#' where \eqn{\bar t} is the weighted mean time.
#'
#' **Test.** With no trend every timepoint has the weighted mean prevalence
#' \eqn{\bar p}, and the slope has variance
#' \eqn{\bar p (1 - \bar p) / \sum_j w_j (t_j - \bar t)^2}. The z-statistic
#' is the slope divided by its standard error (the Cochran-Armitage trend
#' test). Without sites it gives the same p-value as R's
#' `prop.trend.test()`.
#'
#' **Interval.** For a possible slope \eqn{\beta}, take the line with that
#' slope through \eqn{\bar p} at \eqn{\bar t}, cut to \[0, 1\]; call its
#' value at each timepoint \eqn{f_j}. If that line were true, the slope
#' would have variance \eqn{V(\beta) = \sum_j w_j (t_j - \bar t)^2 f_j
#' (1 - f_j) / (\sum_j w_j (t_j - \bar t)^2)^2}. The interval holds the
#' slopes with \eqn{(\hat\beta_{app} - \beta)^2 \le z^2 V(\beta)}. At
#' \eqn{\beta = 0} this is the test, so `reject` is `TRUE` exactly when 0
#' is outside the interval. It is matched to `alternative`: one-sided for
#' `"greater"` and `"less"`. The slope and its interval are divided by
#' \eqn{Se + Sp - 1} to give the trend in true prevalence.
#'
#' **Design effect.** It is \eqn{1 + (n_w - 1) ICC}, where
#' \eqn{n_w = \sum n^2 / \sum n} is the mean number of people per site,
#' weighted by site size, and is applied to every timepoint. With
#' `icc = NULL`, an ICC is estimated within each timepoint that has two or
#' more sites (as in [estimate_prevalence()]), and the estimates are
#' averaged; timepoints with one person per site, or with all negative or
#' all positive, are skipped. With one site per timepoint there is nothing
#' to estimate from, and a supplied `icc` is ignored with a warning.
#'
#' @return A list with:
#'   \item{slope}{Estimated change in true prevalence per unit of time.}
#'   \item{ci_lower, ci_upper}{Confidence interval for the slope;
#'     `ci_upper` is `Inf` for `"greater"` and `ci_lower` is `-Inf` for
#'     `"less"`.}
#'   \item{p_value}{p-value for the chosen `alternative`.}
#'   \item{reject}{`TRUE` if `p_value` is below 1 - `conf_level`.}
#'   \item{n_eff}{Effective sample size, `sum(n) / deff`.}
#'   \item{icc}{ICC used: as supplied, or estimated from the data.}
#'   \item{deff}{Design effect used (1 with no clustering).}
#'
#' @references
#' Armitage, P. (1955) Tests for linear trends in proportions and
#' frequencies. *Biometrics* **11**(3), 375-386. \doi{10.2307/3001775}
#'
#' Cochran, W. G. (1954) Some methods for strengthening the common
#' chi-square tests. *Biometrics* **10**(4), 417-451. \doi{10.2307/3001616}
#'
#' Kish, L. (1965) *Survey Sampling*. Wiley. (Design effect.)
#'
#' Rogan, W. J. & Gladen, B. (1978) Estimating prevalence from the results
#' of a screening test. *American Journal of Epidemiology* **107**(1),
#' 71-76. \doi{10.1093/oxfordjournals.aje.a112510}
#'
#' MMS-SD Study Design Workshop. \url{https://mrc-ide.github.io/MMS-SD_workshop/}
#'
#' @seealso [design_trend()]
#'
#' @export
#'
#' @examples
#' # Five annual rounds, one survey per round
#' test_trend(x = c(5, 8, 9, 14, 17), n = rep(100, 5), time = 2019:2023)
#'
#' # Two sites per round, ICC estimated from the data
#' test_trend(
#'   x    = c(2, 3, 4, 5, 6, 8, 7, 10),
#'   n    = c(50, 45, 55, 50, 48, 52, 47, 51),
#'   time = c(1, 1, 2, 2, 3, 3, 4, 4)
#' )
#'
#' # Imperfect diagnostic test, one-sided (a decline)
#' test_trend(x = c(20, 15, 11, 7), n = rep(120, 4), time = 1:4,
#'            alternative = "less", sensitivity = 0.9, specificity = 0.97)
test_trend <- function(x,
                       n,
                       time        = NULL,
                       alternative = c("two.sided", "greater", "less"),
                       sensitivity = 1,
                       specificity = 1,
                       conf_level  = 0.95,
                       icc         = NULL) {

  .check_counts(x, n)
  if (is.null(time)) time <- seq_along(x)
  if (!is.numeric(time) || length(time) != length(x) || !all(is.finite(time)))
    stop("`time` must be a numeric vector the same length as `x`.", call. = FALSE)
  if (length(unique(time)) < 2)
    stop("`time` must have at least two different timepoints.", call. = FALSE)
  .check_test(sensitivity, specificity)
  .check_number(conf_level, "conf_level", 0, 1, lower_open = TRUE, upper_open = TRUE)
  if (!is.null(icc)) .check_number(icc, "icc", 0, 1)
  alternative <- match.arg(alternative)

  # counts added up within each timepoint
  tj <- sort(unique(time))
  Xj <- vapply(tj, function(t) sum(x[time == t]), numeric(1))
  Nj <- vapply(tj, function(t) sum(n[time == t]), numeric(1))
  pj <- Xj / Nj
  if (sum(Xj) == 0 || sum(Xj) == sum(Nj))
    stop("Every person tested was ", if (sum(Xj) == 0) "negative" else "positive",
         ", so there is no trend to test.", call. = FALSE)

  # design effect: a site is a row of the data, so only timepoints with two
  # or more rows have sites to measure clustering from
  multi <- tj[vapply(tj, function(t) sum(time == t) >= 2, logical(1))]
  if (length(multi) == 0) {
    if (!is.null(icc) && icc > 0)
      warning("`icc` was ignored: no timepoint has two or more sites.",
              call. = FALSE)
    icc <- 0
  } else if (is.null(icc)) {
    # average of the ICCs estimated within each timepoint, skipping those
    # with nothing to estimate from
    icc_tp <- unlist(lapply(multi, function(t) {
      xi <- x[time == t]; ni <- n[time == t]
      p  <- sum(xi) / sum(ni)
      if (mean(ni) == 1 || p == 0 || p == 1) return(NULL)
      .estimate_deff(xi, ni)$icc
    }))
    icc <- if (length(icc_tp) > 0) mean(icc_tp) else 0
  }
  rows <- time %in% multi
  deff <- if (icc > 0) .deff(sum(n[rows]^2) / sum(n[rows]), icc) else 1

  # slope, weighting each timepoint by its effective number tested
  wj       <- Nj / deff
  tw       <- sum(wj * tj) / sum(wj)
  pbar     <- sum(wj * pj) / sum(wj)
  dj       <- tj - tw
  sxx_w    <- sum(wj * dj^2)
  beta_app <- sum(wj * dj * pj) / sxx_w

  # Cochran-Armitage test: with no trend the slope has variance
  # pbar (1 - pbar) / sxx_w
  z <- beta_app / sqrt(pbar * (1 - pbar) / sxx_w)
  p_value <- switch(alternative,
    greater   = stats::pnorm(z, lower.tail = FALSE),
    less      = stats::pnorm(z),
    two.sided = 2 * stats::pnorm(abs(z), lower.tail = FALSE))
  alpha <- 1 - conf_level

  # Interval: the slopes b with gap(b) = (beta_app - b)^2 - z^2 var_at(b)
  # <= 0, where var_at(b) is the slope's variance if the line through pbar
  # at tw with slope b (cut to [0, 1]) were true. gap(0) <= 0 exactly when
  # the test does not reject, so each end is searched for beyond 0 when 0
  # is inside and between beta_app and 0 when it is not; `reject` then
  # always agrees with the interval. var_at(b) is at most 0.25 / sxx_w, so
  # the gap is positive beyond beta_app +/- z * sqrt(0.25 / sxx_w).
  var_at <- function(b) {
    f <- pmin(1, pmax(0, pbar + b * dj))
    sum(wj * dj^2 * f * (1 - f)) / sxx_w^2
  }
  gap    <- function(b, zq) (beta_app - b)^2 - zq^2 * var_at(b)
  ci_end <- function(zq, side) {
    far     <- beta_app + side * (zq * sqrt(0.25 / sxx_w) + 1e-12)
    toward0 <- side * (0 - beta_app) > 0
    from    <- if (toward0 && gap(0, zq) <= 0) 0 else beta_app
    to      <- if (toward0 && gap(0, zq) > 0) 0 else far
    if (gap(from, zq) >= 0) return(from)
    stats::uniroot(function(b) gap(b, zq), sort(c(from, to)), tol = 1e-14)$root
  }
  alpha_lo <- switch(alternative, two.sided = alpha / 2, greater = alpha, less = 0)
  alpha_hi <- switch(alternative, two.sided = alpha / 2, greater = 0, less = alpha)
  corr     <- sensitivity + specificity - 1

  list(
    slope    = beta_app / corr,
    ci_lower = if (alpha_lo > 0)
      ci_end(stats::qnorm(alpha_lo, lower.tail = FALSE), -1) / corr else -Inf,
    ci_upper = if (alpha_hi > 0)
      ci_end(stats::qnorm(alpha_hi, lower.tail = FALSE), 1) / corr else Inf,
    p_value  = p_value,
    reject   = p_value < alpha,
    n_eff    = sum(n) / deff,
    icc      = icc,
    deff     = deff
  )
}
