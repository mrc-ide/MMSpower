#' Test whether prevalence is changing linearly over time
#'
#' @description
#' Analysis function: given the number of positives and the number tested at
#' several timepoints (optionally at several sites per timepoint), fits a
#' straight line to the measured prevalences and tests whether its slope
#' differs from zero. Returns the estimated change in true prevalence per
#' unit of time with a confidence interval, the z-statistic, the p-value
#' and whether the trend is significant. It is the analysis twin of
#' \code{design_trend()}.
#'
#' The line is fitted to the measured (apparent) prevalences, weighting each
#' timepoint by the number of people tested there, and the slope is tested
#' with the Cochran-Armitage trend test. The slope and its interval are then
#' divided by \eqn{Se + Sp - 1} to give the trend in true prevalence.
#' Clustering is allowed for with one design effect applied to every
#' timepoint.
#'
#' \strong{Hypotheses} (\eqn{\beta} is the slope of true prevalence):
#' \describe{
#'   \item{\code{"two.sided"} (default)}{H0: \eqn{\beta = 0} against
#'     H1: \eqn{\beta \ne 0}.}
#'   \item{\code{"greater"}}{H0: \eqn{\beta \le 0} against H1:
#'     \eqn{\beta > 0} (an increasing trend).}
#'   \item{\code{"less"}}{H0: \eqn{\beta \ge 0} against H1:
#'     \eqn{\beta < 0} (a decreasing trend).}
#' }
#'
#' @param x Vector of integer counts of people who tested positive: one per
#'   timepoint, or one per site and timepoint (give the timepoint of each in
#'   \code{time}).
#' @param n Vector of integer totals tested, the same length as \code{x}.
#' @param time Numeric vector giving the timepoint of each value of \code{x}
#'   (e.g. calendar years or survey rounds). If \code{NULL} (default),
#'   \code{1, 2, ...} is used, which only makes sense with one value per
#'   timepoint.
#' @param alternative Direction of the test: \code{"two.sided"} (default),
#'   \code{"greater"} (increasing) or \code{"less"} (decreasing).
#' @param sensitivity Diagnostic sensitivity in (0, 1]: the share of people
#'   who truly have the condition that the test detects as positive.
#'   Default 1 (the test never misses a case).
#' @param specificity Diagnostic specificity in (0, 1]: the share of people
#'   who truly do not have the condition that the test correctly reads as
#'   negative. Default 1 (the test never gives a false positive).
#' @param conf_level Confidence level for the interval on the slope, in
#'   \[0.5, 1); default 0.95. The significance level of the test is
#'   1 - \code{conf_level}.
#' @param icc Optional. Intra-cluster correlation, in \[0, 1\]. If \code{NULL}
#'   (default), it is estimated from the timepoints with two or more sites.
#'   Set to \code{0} for no clustering. It is only used when some timepoint
#'   has two or more sites; otherwise it is ignored (with a warning). Values
#'   below about 1.5e-8 are treated as 0.
#' @param fpc_N Optional positive integer. Population size \strong{per
#'   timepoint}, for a finite-population correction; must be larger than the
#'   number tested at every timepoint. \code{NULL} (default) = no
#'   finite-population correction applied.
#' @param method Method for the test and the interval on the slope. Only
#'   \code{"cochran-armitage"} is available.
#'
#' @details
#' \strong{Line.} The counts are added up within each timepoint, giving
#' \eqn{X_j} positives of \eqn{N_j} tested at time \eqn{t_j}, and measured
#' prevalence \eqn{\hat p_j = X_j / N_j}. Each timepoint is weighted by its
#' number tested divided by the design effect, \eqn{w_j = N_j / D_{eff}}
#' (with \code{fpc_N}, also divided by the finite-population factor
#' \eqn{(N - N_j)/(N - 1)}). The slope is
#'
#' \deqn{\hat\beta_{app} = \frac{\sum_j w_j (t_j - \bar t)\,\hat p_j}
#'   {\sum_j w_j (t_j - \bar t)^2},}
#'
#' where \eqn{\bar t} is the weighted mean time.
#'
#' \strong{Test.} If prevalence is not changing, every timepoint has the
#' same prevalence, estimated by the weighted mean \eqn{\bar p} of the
#' \eqn{\hat p_j}. The slope then has variance
#' \eqn{\bar p (1 - \bar p) / \sum_j w_j (t_j - \bar t)^2}, and the
#' z-statistic is \eqn{\hat\beta_{app}} divided by the square root of this
#' variance. This is the Cochran-Armitage trend test. Without sites and
#' \code{fpc_N} it gives the same p-value as R's \code{prop.trend.test()}.
#'
#' \strong{Interval.} For a possible slope \eqn{\beta}, take the line with
#' that slope through \eqn{\bar p} at \eqn{\bar t}; where it goes below 0
#' or above 1, 0 or 1 is used. Call its value at each timepoint \eqn{f_j}.
#' The slope's variance if that line were true is
#' \eqn{V(\beta) = \sum_j w_j (t_j - \bar t)^2 f_j (1 - f_j) /
#' (\sum_j w_j (t_j - \bar t)^2)^2}. The interval holds the slopes with
#' \eqn{(\hat\beta_{app} - \beta)^2 \le z^2 V(\beta)}, and its ends are
#' found numerically. At \eqn{\beta = 0} this is the test above, so
#' \code{reject} agrees with whether 0 is outside the interval. The slope
#' and its interval are divided by \eqn{Se + Sp - 1} to give the trend in
#' true prevalence. The interval matches \code{alternative}:
#' \eqn{[L, \infty)} for \code{"greater"}, \eqn{(-\infty, U]} for
#' \code{"less"} and two-sided for \code{"two.sided"}.
#'
#' \strong{Design effect.} With \code{icc = NULL}, an ICC is estimated within
#' each timepoint that has two or more sites (as in
#' \code{estimate_prevalence()}), and the estimates are averaged. The design
#' effect \eqn{1 + (\bar n - 1) ICC}, with \eqn{\bar n} the mean number of
#' people per site, is then applied to every timepoint.
#'
#' @return A named list with the following fields, in the order returned:
#'   \item{slope}{Estimated change in true prevalence per unit of time.}
#'   \item{slope_app}{Estimated change in measured prevalence per unit of
#'     time.}
#'   \item{statistic}{The z-statistic.}
#'   \item{p_value}{p-value for the chosen \code{alternative}.}
#'   \item{reject}{\code{TRUE} if \code{p_value} is below
#'     1 - \code{conf_level}.}
#'   \item{alternative}{As supplied.}
#'   \item{ci_lower, ci_upper}{Confidence interval for the true slope;
#'     \code{ci_upper} is \code{Inf} for \code{"greater"} and \code{ci_lower}
#'     is \code{-Inf} for \code{"less"}.}
#'   \item{prevalence_start_est, prevalence_end_est}{Fitted true prevalence
#'     at the first and last timepoint, between 0 and 1.}
#'   \item{times}{The timepoints, sorted.}
#'   \item{apparent_prev}{Measured prevalence at each timepoint.}
#'   \item{x_by_time, n_by_time}{Positives and number tested at each
#'     timepoint.}
#'   \item{n_timepoints}{Number of timepoints.}
#'   \item{n_total}{Total number of people tested.}
#'   \item{conf_level}{As supplied.}
#'   \item{sensitivity, specificity}{As supplied.}
#'   \item{icc_used}{ICC used. If \code{icc = NULL}, it is the value
#'     estimated from the data. If \code{icc} was supplied, it is that value,
#'     unless it was below about 1.5e-8 or was ignored (no timepoint with two
#'     or more sites), then it is 0.}
#'   \item{deff}{Design effect applied at every timepoint (1 with no
#'     clustering).}
#'   \item{fpc_N}{As supplied, or \code{NULL}.}
#'   \item{method}{\code{"cochran-armitage"}.}
#'
#' @references
#' MMS-SD Study Design Workshop. \url{https://mrc-ide.github.io/MMS-SD_workshop/}
#'
#' Armitage, P. (1955) Tests for linear trends in proportions and
#' frequencies. \emph{Biometrics} \strong{11}(3), 375-386. \doi{10.2307/3001775}
#'
#' Cochran, W. G. (1954) Some methods for strengthening the common
#' chi-square tests. \emph{Biometrics} \strong{10}(4), 417-451.
#' \doi{10.2307/3001616}
#'
#' Kish, L. (1965) \emph{Survey Sampling}. Wiley. (Design effect.)
#'
#' Cochran, W. G. (1977) \emph{Sampling Techniques}, 3rd ed. Wiley.
#' (Finite-population correction.)
#'
#' Rogan, W. J. & Gladen, B. (1978) Estimating prevalence from the results
#' of a screening test. \emph{American Journal of Epidemiology}
#' \strong{107}(1), 71-76. \doi{10.1093/oxfordjournals.aje.a112510}
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
#' # Clustered: two sites per round
#' test_trend(
#'   x    = c(2, 3, 4, 5, 6, 8, 7, 10),
#'   n    = c(50, 45, 55, 50, 48, 52, 47, 51),
#'   time = c(1, 1, 2, 2, 3, 3, 4, 4)
#' )
#'
#' # Imperfect diagnostic test, one-sided (a decline)
#' test_trend(x = c(20, 15, 11, 7), n = rep(120, 4), time = 1:4,
#'            alternative = "less", sensitivity = 0.9, specificity = 0.97)
test_trend <- function(x, n, time = NULL,
                       alternative = "two.sided",
                       sensitivity = 1,
                       specificity = 1,
                       conf_level  = 0.95,
                       icc         = NULL,
                       fpc_N       = NULL,
                       method      = "cochran-armitage") {

  # ---- validate x and n ----
  # The same checks and messages as estimate_prevalence().
  if (is.logical(x) || is.logical(n))
    stop("`x` and `n` must be numeric, not logical (got class `",
         class(x)[1], "` for x, `", class(n)[1], "` for n). ",
         "Note: a plain `NA` is logical in R. Remove missing observations ",
         "before calling.")
  if (!is.numeric(x) || !is.numeric(n))
    stop("`x` and `n` must be numeric (got class `", class(x)[1], "` for x, ",
         "`", class(n)[1], "` for n).")
  if (length(x) == 0 || length(n) == 0)
    stop("`x` and `n` must each contain at least one value. ",
         "Supply at least one count and total.")
  # For each check below, `i` is the position of the first bad element.
  if (!all(is.finite(x))) {
    i <- which(!is.finite(x))[1]
    stop("`x` contains a missing or infinite value (found x[", i, "] = ",
         x[i], "). The number of positives must be a finite non-negative ",
         "integer.")
  }
  if (!all(is.finite(n))) {
    i <- which(!is.finite(n))[1]
    stop("`n` contains a missing or infinite value (found n[", i, "] = ",
         n[i], "). The number tested must be a finite positive integer.")
  }
  if (length(x) != length(n))
    stop("`x` and `n` must have the same length ",
         "(got length(x) = ", length(x), ", length(n) = ", length(n), "). ",
         "Each value of `x` is a number of positives, and the value of `n` ",
         "in the same position is the number tested.")
  if (any(x < 0)) {
    i <- which(x < 0)[1]
    stop("`x` must be non-negative (found x[", i, "] = ", x[i], "). ",
         "Counts cannot be negative.")
  }
  if (any(n <= 0)) {
    i <- which(n <= 0)[1]
    stop("`n` must be positive everywhere (found n[", i, "] = ", n[i], "). ",
         "Each value must count at least one person tested.")
  }
  if (any(x != floor(x))) {
    i <- which(x != floor(x))[1]
    stop("`x` must contain integers, not decimals (found x[", i, "] = ", x[i], ").")
  }
  if (any(n != floor(n))) {
    i <- which(n != floor(n))[1]
    stop("`n` must contain integers, not decimals (found n[", i, "] = ", n[i], ").")
  }
  if (any(x > n)) {
    i <- which(x > n)[1]
    stop("`x` cannot exceed `n` (found x[", i, "] = ", x[i],
         " > n[", i, "] = ", n[i], "). ",
         "The number of positives cannot be more than the number tested.")
  }

  # ---- validate time ----
  if (is.null(time)) {
    time <- seq_along(x)
  } else {
    if (!is.numeric(time) || length(time) != length(x))
      stop("`time` must be a numeric vector the same length as `x` (got class `",
           class(time)[1], "`, length ", length(time), ", for ", length(x),
           " values of `x`).")
    if (!all(is.finite(time))) {
      i <- which(!is.finite(time))[1]
      stop("`time` contains a missing or infinite value (found time[", i,
           "] = ", time[i], ").")
    }
  }
  if (length(unique(time)) < 2)
    stop("At least two different timepoints are needed to fit a slope (got ",
         length(unique(time)), ").")

  # ---- validate single-value parameters ----
  # Each of the parameters below must be a single value.
  if (length(sensitivity) != 1 || !is.numeric(sensitivity))
    stop("`sensitivity` must be a single number in (0, 1] (got class `",
         class(sensitivity)[1], "`, length ", length(sensitivity), "). ",
         "Note: `TRUE`, `FALSE` and a plain `NA` are logical, not numeric. ",
         "Pass 1 if the test never misses a case.")
  if (length(specificity) != 1 || !is.numeric(specificity))
    stop("`specificity` must be a single number in (0, 1] (got class `",
         class(specificity)[1], "`, length ", length(specificity), "). ",
         "Note: `TRUE`, `FALSE` and a plain `NA` are logical, not numeric. ",
         "Pass 1 if the test never gives a false positive.")
  if (length(conf_level) != 1 || !is.numeric(conf_level))
    stop("`conf_level` must be a single number in (0, 1) (got class `",
         class(conf_level)[1], "`, length ", length(conf_level), ").")
  if (!is.null(icc) && (length(icc) != 1 || !is.numeric(icc)))
    stop("`icc` must be a single number in [0, 1] (got class `",
         class(icc)[1], "`, length ", length(icc), "). ",
         "To estimate ICC from the data, leave `icc = NULL`.")
  if (!is.character(alternative) || length(alternative) != 1)
    stop("`alternative` must be one of 'two.sided', 'greater', or 'less' ",
         "(got class `", class(alternative)[1], "`, length ",
         length(alternative), ").")
  if (!alternative %in% c("two.sided", "greater", "less"))
    stop("`alternative` must be one of 'two.sided', 'greater', or 'less' ",
         "(got '", alternative, "').")
  if (!is.character(method) || length(method) != 1 || method != "cochran-armitage")
    stop("`method` must be 'cochran-armitage', the only method available (got ",
         if (is.character(method) && length(method) == 1) paste0("'", method, "'")
         else paste0("class `", class(method)[1], "`, length ", length(method)),
         ").")

  # Check for NA, NaN and Inf before the range checks below.
  if (!is.finite(sensitivity))
    stop("`sensitivity` must be a single finite number (got ", sensitivity, ").")
  if (!is.finite(specificity))
    stop("`specificity` must be a single finite number (got ", specificity, ").")
  if (!is.finite(conf_level))
    stop("`conf_level` must be a single finite number (got ", conf_level, ").")
  if (!is.null(icc) && !is.finite(icc))
    stop("`icc` must be a single finite number (got ", icc, "). ",
         "To estimate ICC from the data, leave `icc = NULL`.")

  if (sensitivity < 0 || sensitivity > 1)
    stop("`sensitivity` must be in (0, 1] (got ", sensitivity, "). ",
         "`sensitivity` is a diagnostic probability and cannot be ",
         if (sensitivity < 0) "negative." else "greater than 1.")
  if (sensitivity == 0)
    stop("`sensitivity` must be in (0, 1] (got 0). ",
         "A sensitivity of 0 means everyone who truly has the condition ",
         "tests negative.")
  if (specificity < 0 || specificity > 1)
    stop("`specificity` must be in (0, 1] (got ", specificity, "). ",
         "`specificity` is a diagnostic probability and cannot be ",
         if (specificity < 0) "negative." else "greater than 1.")
  if (specificity == 0)
    stop("`specificity` must be in (0, 1] (got 0). ",
         "A specificity of 0 means everyone who truly does not have the ",
         "condition tests positive.")
  # A sum above 1 only by rounding (within about 1.5e-8 of 1) counts as 1.
  correction <- sensitivity + specificity - 1
  if (correction < sqrt(.Machine$double.eps))
    stop("`sensitivity` + `specificity` must exceed 1 for the Rogan-Gladen correction ",
         "(got ", sensitivity, " + ", specificity, " = ", sensitivity + specificity, "). ",
         "The correction divides by `sensitivity` + `specificity` - 1, so a sum ",
         "of 1 or less cannot be corrected.")
  if (correction < 0.1)
    warning("`sensitivity` + `specificity` = ",
            format(sensitivity + specificity, digits = 15),
            ", which is very close to 1. Correcting for this much test error ",
            "makes the estimated slope and its confidence interval very ",
            "wide, and small changes in the assumed sensitivity or specificity ",
            "will change them a lot.")
  if (!is.null(icc) && (icc < 0 || icc > 1))
    stop("`icc` must be in [0, 1] (got ", icc, "). ",
         "`icc` is a correlation and cannot be ",
         if (icc < 0) "negative. " else "greater than 1. ",
         "To estimate ICC from the data, leave `icc = NULL`.")
  if (conf_level <= 0 || conf_level >= 1)
    stop("`conf_level` must be in (0, 1) (got ", conf_level, "). ",
         "For example, use 0.95 for a 95% confidence interval.")
  if (conf_level < 0.5)
    stop("`conf_level` must be at least 0.5 (got ", conf_level, "). ",
         "Below 0.5, when prevalence is not changing, the test would wrongly ",
         "report a trend more than half the time.")
  if (!is.null(fpc_N) &&
      (!is.numeric(fpc_N) || length(fpc_N) != 1 || !is.finite(fpc_N) ||
       fpc_N < 1 || fpc_N != floor(fpc_N)))
    stop("`fpc_N` must be a single finite positive integer (got ",
         if (!is.numeric(fpc_N)) paste0("class `", class(fpc_N)[1], "`")
         else if (length(fpc_N) != 1) paste0("length = ", length(fpc_N))
         else fpc_N,
         "). `fpc_N` is the population size per timepoint. ",
         "Set `fpc_N = NULL` to skip the finite-population correction.")

  # A supplied icc below this tiny threshold (e.g. 1e-12, a rounding
  # leftover) is treated as exactly 0.
  if (!is.null(icc) && icc < sqrt(.Machine$double.eps))
    icc <- 0

  # ---- add up the counts within each timepoint ----
  ord  <- order(time)
  time <- time[ord]; x <- x[ord]; n <- n[ord]
  tj   <- sort(unique(time))
  Xj   <- vapply(tj, function(t) sum(x[time == t]), numeric(1))
  Nj   <- vapply(tj, function(t) sum(n[time == t]), numeric(1))
  pj   <- Xj / Nj

  # ---- design effect ----
  # A site is a row of the data. Only a timepoint with two or more rows has
  # sites to measure clustering from; with one row per timepoint the row is
  # the whole timepoint, and no design effect can be formed.
  multi_tp     <- tj[vapply(tj, function(t) sum(time == t) >= 2, logical(1))]
  has_clusters <- length(multi_tp) > 0
  n_bar_cl     <- if (has_clusters) mean(n[time %in% multi_tp]) else NA_real_
  no_structure <- !has_clusters || n_bar_cl == 1

  if (is.null(icc)) {
    # Estimate an ICC within each timepoint with two or more sites (as in
    # estimate_prevalence), skipping timepoints that carry no information
    # (one person per site, or all negative / all positive), and average.
    icc_tp <- c()
    for (t in multi_tp) {
      xi <- x[time == t]; ni <- n[time == t]
      if (mean(ni) == 1) next
      p_pool  <- sum(xi) / sum(ni)
      var_srs <- mean(p_pool * (1 - p_pool) / ni)
      if (var_srs == 0) next
      d_tp   <- max(stats::var(xi / ni) / var_srs, 1)
      icc_tp <- c(icc_tp, min(max((d_tp - 1) / (mean(ni) - 1), 0), 1))
    }
    icc_used <- if (length(icc_tp) > 0) mean(icc_tp) else 0
  } else if (no_structure) {
    if (icc > 0)
      warning("`icc` = ", icc, " was ignored: no timepoint has two or more ",
              "sites (or every site has one person), so there is no ",
              "clustering to adjust for. `icc_used` is reported as 0.")
    icc_used <- 0
  } else {
    icc_used <- icc
  }
  deff <- if (icc_used > 0 && !no_structure) 1 + (n_bar_cl - 1) * icc_used else 1

  # ---- finite-population correction for each timepoint ----
  if (!is.null(fpc_N)) {
    if (fpc_N < max(Nj))
      stop("`fpc_N` (", fpc_N, ") is smaller than the number tested at a ",
           "timepoint (", max(Nj), "): you cannot test more people than there ",
           "are in the population.")
    if (fpc_N == max(Nj))
      stop("`fpc_N` (", fpc_N, ") equals the number tested at a timepoint: ",
           "the whole population was tested there (a census), so that ",
           "timepoint has no sampling uncertainty. Set `fpc_N = NULL` if no ",
           "finite-population correction is needed.")
    fpc2 <- (fpc_N - Nj) / (fpc_N - 1)
  } else {
    fpc2 <- rep(1, length(tj))
  }

  # ---- all negative or all positive ----
  if (sum(Xj) == 0)
    stop("Every person tested was negative (`x` is 0 everywhere), so the ",
         "measured prevalence is 0 at every timepoint and there is no trend ",
         "to test.")
  if (sum(Xj) == sum(Nj))
    stop("Every person tested was positive (`x` equals `n` everywhere), so ",
         "the measured prevalence is 1 at every timepoint and there is no ",
         "trend to test.")

  # ---- slope, weighting each timepoint by its effective number tested ----
  wj       <- Nj / (deff * fpc2)
  tw       <- sum(wj * tj) / sum(wj)
  pbar     <- sum(wj * pj) / sum(wj)
  dj       <- tj - tw
  sxx_w    <- sum(wj * dj^2)
  beta_app <- sum(wj * dj * pj) / sxx_w

  # ---- Cochran-Armitage test ----
  # With no trend every timepoint has prevalence pbar, so the slope has
  # variance pbar * (1 - pbar) / sxx_w.
  z_stat <- beta_app / sqrt(pbar * (1 - pbar) / sxx_w)

  p_value <- switch(alternative,
    greater   = stats::pnorm(z_stat, lower.tail = FALSE),
    less      = stats::pnorm(z_stat, lower.tail = TRUE),
    two.sided = 2 * stats::pnorm(abs(z_stat), lower.tail = FALSE))
  alpha  <- 1 - conf_level
  reject <- p_value < alpha

  # ---- interval on the slope, matched to `alternative` ----
  # For a slope b, the line through pbar at tw has value pbar + b * dj, cut
  # to [0, 1]; var_at(b) is the slope's variance if that line were true.
  # The interval holds the b with gap(b) = (beta_app - b)^2 - z^2 var_at(b)
  # <= 0. gap(0) <= 0 exactly when the test does not reject, so each end is
  # searched for on the far side of 0 when 0 is inside, and between
  # beta_app and 0 when it is not: `reject` then always agrees with whether
  # 0 is outside the interval. var_at(b) is at most var_max, so the gap is
  # positive beyond beta_app +/- z * sqrt(var_max).
  var_at  <- function(b) {
    f <- pmin(1, pmax(0, pbar + b * dj))
    sum(wj * dj^2 * f * (1 - f)) / sxx_w^2
  }
  var_max <- 0.25 / sxx_w
  gap     <- function(b, z) (beta_app - b)^2 - z^2 * var_at(b)
  ci_end  <- function(z, side) {
    far     <- beta_app + side * (z * sqrt(var_max) + 1e-12)
    toward0 <- side * (0 - beta_app) > 0
    from    <- if (toward0 && gap(0, z) <= 0) 0 else beta_app
    to      <- if (toward0 && gap(0, z) > 0) 0 else far
    if (gap(from, z) >= 0) return(from)
    stats::uniroot(function(b) gap(b, z), sort(c(from, to)), tol = 1e-14)$root
  }
  alpha_lo <- switch(alternative, two.sided = alpha / 2, greater = alpha, less = 0)
  alpha_hi <- switch(alternative, two.sided = alpha / 2, greater = 0,     less = alpha)
  slope_true <- beta_app / correction
  ci_lower <- if (alpha_lo > 0)
    ci_end(stats::qnorm(alpha_lo, lower.tail = FALSE), -1) / correction else -Inf
  ci_upper <- if (alpha_hi > 0)
    ci_end(stats::qnorm(alpha_hi, lower.tail = FALSE),  1) / correction else Inf

  # ---- fitted prevalence at the first and last timepoint ----
  fit_app_start <- pbar + beta_app * (min(tj) - tw)
  fit_app_end   <- pbar + beta_app * (max(tj) - tw)
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
    method               = "cochran-armitage"
  )
}
