#' Test whether prevalence exceeds a decision threshold
#'
#' @description
#' Analysis function: given observed counts, tests whether the true prevalence
#' is above or below a specified action threshold. Returns a z-statistic,
#' p-value, reject/fail-to-reject decision, and a confidence interval.
#'
#' The test works on the apparent-prevalence scale (accounting for imperfect
#' diagnostics via the Rogan-Gladen transform) and adjusts for clustering
#' (design effect) and finite-population corrections in the same way as
#' \code{estimate_prevalence()}.
#'
#' **Hypothesis:**
#' \describe{
#'   \item{`"greater"` (default)}{H0: p <= threshold vs Ha: p > threshold.
#'     Use this when asking "has prevalence exceeded the action threshold?"}
#'   \item{`"less"`}{H0: p >= threshold vs Ha: p < threshold.}
#'   \item{`"two.sided"`}{H0: p = threshold vs Ha: p != threshold.}
#' }
#'
#' @param x Integer vector of positive counts per cluster/site.
#' @param n Integer vector of total samples per cluster/site (same length as x).
#' @param threshold Numeric in (0, 1). Decision threshold on the
#'   true-prevalence scale (e.g. 0.05 for the WHO pfhrp2/3 deletion cutoff).
#' @param alternative Direction of the alternative hypothesis: `"greater"`
#'   (default), `"less"`, or `"two.sided"`.
#' @param sensitivity Diagnostic sensitivity in (0, 1]; default 1.
#' @param specificity Diagnostic specificity in (0, 1]; default 1.
#' @param conf_level Confidence level for the reported CI; also sets the
#'   significance level as alpha = 1 - conf_level. Default 0.95.
#' @param icc Optional. Intra-cluster correlation. `NULL` (default) = estimate
#'   from data. `0` = force SRS (no clustering adjustment).
#' @param fpc_N Optional. Total population size for a finite-population
#'   correction. `NULL` (default) = no FPC.
#' @param ci_method Method for the reported confidence interval on true
#'   prevalence: `"wald"` (default), `"clopper-pearson"`, or
#'   `"agresti-coull"`. Same three methods as `estimate_prevalence()`. This
#'   affects only the reported interval, not the hypothesis test, which
#'   always uses the null-variance z-statistic.
#'
#' @details
#' **Test statistic**: the z-statistic is computed using the \emph{null}
#' variance (evaluated at \code{threshold}, not at \eqn{\hat{p}}) for tighter
#' Type I error control:
#'
#' \deqn{z = \frac{\hat{p}_{app} - \theta_{app}}
#'              {\sqrt{\theta_{app}(1-\theta_{app}) / n_{adj}}}}
#'
#' where \eqn{\theta_{app} = \theta \cdot Se + (1-\theta)(1-Sp)} is the
#' threshold on the apparent-prevalence scale, and
#' \eqn{n_{adj} = n_{eff} / fpc^2} folds in the design effect and FPC.
#'
#' The reported CI is \strong{matched to \code{alternative}}: a lower
#' one-sided interval \eqn{[L, 1]} for \code{"greater"}, an upper one-sided
#' interval \eqn{[0, U]} for \code{"less"}, and the usual two-sided
#' interval for \code{"two.sided"}, so the interval and the test point the
#' same way. \code{ci_method} selects Wald (default), Clopper-Pearson, or
#' Agresti-Coull for the interval shape.
#'
#' Note that the interval and \code{reject} are not \emph{exactly}
#' equivalent near the decision boundary: the test uses the null-variance
#' (score) SE evaluated at \code{threshold}, while the CI uses each
#' method's own SE evaluated at \eqn{\hat p}. They agree in direction and
#' away from the boundary; a score-based (Wilson) interval -- not offered
#' here -- would be the exact test inversion.
#'
#' @section Equations and sources:
#' Direct workshop material (MMS-SD Study Design Workshop,
#' \url{https://mrc-ide.github.io/MMS-SD_workshop/}) plus three extensions:
#' \itemize{
#'   \item \emph{One-sample z-test for a proportion against a fixed value}
#'     \eqn{Z = (\hat p - p_0)/\sqrt{\hat p(1-\hat p)/n}} -- Module 3
#'     "Hypothesis testing", slide "Null hypothesis testing" (lecture
#'     slides p. 8); \eqn{\alpha} and the +/-1.96 critical values, same
#'     slide.
#'   \item \emph{Design effect / effective sample size}
#'     \eqn{D_{eff} = 1 + (\bar n - 1)\,r}, \eqn{N_{eff} = N/D_{eff}} --
#'     Module 5 "Dealing with over-dispersion", slides "The effective
#'     sample size" and "Why is the ICC useful?" (pp. 4-5).
#'   \item \emph{Null-variance z} (SE evaluated at \code{threshold}, not
#'     \eqn{\hat p}) -- a standard refinement of the Module 3 statistic,
#'     not spelled out in the workshop.
#'   \item \emph{Rogan-Gladen correction}, \emph{Clopper-Pearson} and
#'     \emph{Agresti-Coull} intervals, and the \emph{finite-population
#'     correction} -- not in the workshop; Rogan & Gladen (1978),
#'     Clopper & Pearson (1934), Agresti & Coull (1998), Cochran (1977).
#' }
#'
#' @return A named list:
#'   \item{statistic}{z-statistic (on the apparent-prevalence scale)}
#'   \item{p_value}{p-value for the chosen \code{alternative}}
#'   \item{reject}{Logical: \code{TRUE} if \code{p_value < 1 - conf_level}}
#'   \item{threshold}{Decision threshold as supplied}
#'   \item{alternative}{Alternative hypothesis as supplied}
#'   \item{prevalence}{Rogan-Gladen corrected point estimate}
#'   \item{ci_lower}{Lower CI bound on true prevalence. \code{0} when
#'     \code{alternative = "less"} (upper one-sided interval).}
#'   \item{ci_upper}{Upper CI bound on true prevalence. \code{1} when
#'     \code{alternative = "greater"} (lower one-sided interval).}
#'   \item{ci_method}{CI method used (as supplied)}
#'   \item{n_total}{Total samples}
#'   \item{n_eff}{Effective independent sample size before the FPC:
#'     \code{n_total / deff}}
#'   \item{n_eff_adj}{Effective sample size the CI and test statistic actually
#'     use: \code{n_eff} divided by the squared FPC factor (equals
#'     \code{n_eff} when \code{fpc_N} is \code{NULL})}
#'   \item{conf_level}{Confidence level (as supplied)}
#'   \item{sensitivity}{Sensitivity (as supplied)}
#'   \item{specificity}{Specificity (as supplied)}
#'   \item{icc_used}{ICC applied: estimated from data or as supplied}
#'   \item{deff}{Design effect applied}
#'   \item{fpc_N}{\code{fpc_N} as supplied, or \code{NULL}}
#'
#' @references
#' MMS-SD Study Design Workshop, Modules 3 (hypothesis testing) and 5 (ICC /
#' design effect). \url{https://mrc-ide.github.io/MMS-SD_workshop/}
#'
#' Rogan WJ, Gladen B (1978). Estimating prevalence from the results of a
#' screening test. American Journal of Epidemiology 107(1):71-76.
#'
#' Clopper CJ, Pearson ES (1934). The use of confidence or fiducial limits
#' illustrated in the case of the binomial. Biometrika 26(4):404-413.
#'
#' Agresti A, Coull BA (1998). Approximate is better than "exact" for
#' interval estimation of binomial proportions. The American Statistician
#' 52(2):119-126.
#'
#' Cochran WG (1977). Sampling Techniques, 3rd ed. Wiley.
#'
#' @export
#'
#' @examples
#' # Single site: test whether prevalence > 5%
#' test_threshold(x = 8, n = 100, threshold = 0.05)
#'
#' # Multi-site clustered data
#' test_threshold(
#'   x = c(2, 5, 1, 8, 6),
#'   n = c(50, 60, 40, 80, 70),
#'   threshold = 0.05
#' )
#'
#' # Two-sided test with imperfect diagnostic
#' test_threshold(x = 12, n = 100, threshold = 0.10,
#'                alternative = "two.sided",
#'                sensitivity = 0.9, specificity = 0.95)
test_threshold <- function(x,
                           n,
                           threshold,
                           alternative = "greater",
                           sensitivity = 1,
                           specificity = 1,
                           conf_level  = 0.95,
                           icc         = NULL,
                           fpc_N       = NULL,
                           ci_method   = "wald") {

  # ---- validate x and n ----
  if (is.logical(x) || is.logical(n))
    stop("`x` and `n` must be numeric, not logical (got class `",
         class(x)[1], "` for x, `", class(n)[1], "` for n). ",
         "Note: plain `NA` is logical in R -- filter missing observations before calling.")
  if (!is.numeric(x) || !is.numeric(n))
    stop("`x` and `n` must be numeric vectors (got class `", class(x)[1], "` for x, ",
         "`", class(n)[1], "` for n). Supply integer or double counts.")
  if (length(x) == 0 || length(n) == 0)
    stop("`x` and `n` must be non-empty vectors. ",
         "Supply at least one cluster's count and total.")
  if (!all(is.finite(x)))
    stop("`x` contains NA, NaN, or infinite values. ",
         "All counts must be finite non-negative integers.")
  if (!all(is.finite(n)))
    stop("`n` contains NA, NaN, or infinite values. ",
         "All cluster totals must be finite positive integers.")
  if (length(x) != length(n))
    stop("`x` and `n` must have the same length ",
         "(got length(x) = ", length(x), ", length(n) = ", length(n), ").")
  if (any(x != floor(x)))
    stop("`x` must contain whole numbers (found x[",
         which(x != floor(x))[1], "] = ", x[which(x != floor(x))[1]], ").")
  if (any(n != floor(n)))
    stop("`n` must contain whole numbers (found n[",
         which(n != floor(n))[1], "] = ", n[which(n != floor(n))[1]], ").")
  if (any(n <= 0))
    stop("`n` must be positive for every cluster (found n[",
         which(n <= 0)[1], "] = ", n[which(n <= 0)[1]], ").")
  if (any(x < 0))
    stop("`x` must be non-negative (found x[",
         which(x < 0)[1], "] = ", x[which(x < 0)[1]], ").")
  if (any(x > n))
    stop("`x` cannot exceed `n` (found x[", which(x > n)[1], "] = ",
         x[which(x > n)[1]], " > n[", which(x > n)[1], "] = ",
         n[which(x > n)[1]], ").")

  # ---- validate scalar parameters ----
  if (length(threshold) != 1 || !is.numeric(threshold))
    stop("`threshold` must be a single number in (0, 1) (got class `",
         class(threshold)[1], "`, length ", length(threshold), ").")
  if (!is.finite(threshold) || threshold <= 0 || threshold >= 1)
    stop("`threshold` must be strictly between 0 and 1 (got ", threshold, "). ",
         "It represents the action threshold on the true-prevalence scale.")
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
  if (length(conf_level) != 1 || !is.numeric(conf_level))
    stop("`conf_level` must be a single number (got ",
         if (!is.numeric(conf_level)) paste0("class `", class(conf_level)[1], "`")
         else paste0("length = ", length(conf_level)), ").")
  if (!is.null(icc) && (length(icc) != 1 || !is.numeric(icc)))
    stop("`icc` must be a single number or NULL (got ",
         if (!is.numeric(icc)) paste0("class `", class(icc)[1], "`")
         else paste0("length = ", length(icc)), "). ",
         "Set `icc = NULL` to estimate ICC from the data.")
  if (!is.null(fpc_N) && length(fpc_N) != 1)
    stop("`fpc_N` must be a single number or NULL (got length = ", length(fpc_N), ").")

  if (!is.finite(sensitivity))
    stop("`sensitivity` must be a finite number (got ", sensitivity, ").")
  if (sensitivity <= 0 || sensitivity > 1)
    stop("`sensitivity` must be in (0, 1] (got ", sensitivity, ").")
  if (!is.finite(specificity))
    stop("`specificity` must be a finite number (got ", specificity, ").")
  if (specificity <= 0 || specificity > 1)
    stop("`specificity` must be in (0, 1] (got ", specificity, ").")
  correction <- sensitivity + specificity - 1
  if (correction <= 0)
    stop("sensitivity + specificity must exceed 1 for the Rogan-Gladen correction ",
         "(got ", sensitivity, " + ", specificity, " = ", sensitivity + specificity, "). ",
         "With se + sp <= 1 the test performs at or below chance.")
  if (correction < 0.1)
    warning("sensitivity + specificity = ", round(sensitivity + specificity, 4),
            " is very close to 1. The Rogan-Gladen adjustment is numerically unstable.")
  if (!is.finite(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("`conf_level` must be strictly between 0 and 1 (got ", conf_level, ").")
  if (!is.null(icc) && (!is.finite(icc) || icc < 0 || icc > 1))
    stop("`icc` must be in [0, 1] (got ", icc, "). ",
         "Set `icc = NULL` to estimate ICC from the data.")
  if (!is.null(fpc_N) && (!is.numeric(fpc_N) || !is.finite(fpc_N) || fpc_N <= 0))
    stop("`fpc_N` must be a finite positive number (got ", fpc_N, ").")
  if (!is.character(ci_method) || length(ci_method) != 1)
    stop("`ci_method` must be a single character string: ",
         "'wald', 'clopper-pearson', or 'agresti-coull' (got class `",
         class(ci_method)[1], "`, length ", length(ci_method), ").")
  if (!ci_method %in% c("wald", "clopper-pearson", "agresti-coull"))
    stop("`ci_method` must be one of 'wald', 'clopper-pearson', or ",
         "'agresti-coull' (got '", ci_method, "').")

  n_clusters <- length(n)
  n_total    <- sum(n)

  if (!is.null(fpc_N) && fpc_N <= n_total)
    stop("`fpc_N` (", fpc_N, ") must be greater than the total sample size (",
         n_total, "). At `fpc_N = n_total` the whole population was surveyed, ",
         "the sampling variance is zero, and the test statistic is undefined. ",
         "Set `fpc_N = NULL` if no FPC is needed.")

  p_hat <- sum(x) / n_total   # apparent prevalence

  # ---- design effect / ICC ----
  # TODO(review): uses the arithmetic mean cluster size mean(n) (matches
  # Module 5's "average cluster size"); classical Kish for unequal
  # clusters uses the size-weighted sum(n^2)/sum(n). Same open decision as
  # estimate_prevalence() -- keep the two in lockstep. Do not "fix"
  # unilaterally; the team resolves the convention.
  n_bar <- mean(n)

  if (is.null(icc)) {
    if (n_clusters < 2 || n_bar == 1) {
      icc_used <- 0
      deff     <- 1
    } else {
      p_i     <- x / n
      var_obs <- stats::var(p_i)
      var_srs <- mean(p_hat * (1 - p_hat) / n)

      deff <- if (var_srs > 0) var_obs / var_srs else 1
      deff <- max(deff, 1)

      icc_used <- (deff - 1) / (n_bar - 1)
      icc_used <- min(max(icc_used, 0), 1)
      deff     <- 1 + (n_bar - 1) * icc_used
    }
  } else if (n_clusters < 2 || n_bar == 1) {
    # A supplied icc has no effect with a single cluster (or clusters all of
    # size 1): the Kish denominator is undefined, so fall back to SRS -- same
    # as the icc = NULL path above.
    icc_used <- 0
    deff     <- 1
  } else {
    icc_used <- icc
    deff     <- 1 + (n_bar - 1) * icc_used
  }

  n_eff <- n_total / deff

  # FPC on the collected count n_total (the actual sampling fraction), NOT
  # n_eff -- same as estimate_prevalence() / design_precision(). The design
  # effect and the FPC are independent adjustments (Cochran 1977 sec. 2.8;
  # Kish 1965). fpc_N > n_total is guaranteed above.
  fpc   <- if (!is.null(fpc_N)) sqrt((fpc_N - n_total) / (fpc_N - 1)) else 1

  # Variance-equivalent simple sample size (same as estimate_prevalence)
  n_eff_adj <- n_eff / (fpc^2)

  # ---- threshold on apparent scale ----
  theta_app <- .apparent_prev(threshold, sensitivity, specificity)

  # ---- z-statistic using null variance ----
  # Null SE evaluated at theta_app, not p_hat, for correct Type I error control.
  se_null <- sqrt(theta_app * (1 - theta_app) / n_eff_adj)

  if (se_null == 0)
    stop("Null standard error is zero (theta_app = ", theta_app,
         ", n_eff_adj = ", round(n_eff_adj, 2), "). ",
         "`threshold` is validated to (0, 1), so this only happens when the ",
         "FPC drives the effective sample size to infinity, i.e. `fpc_N` is ",
         "just above `n_eff`. Set `fpc_N = NULL`, or use a larger `fpc_N`.")

  z_stat <- (p_hat - theta_app) / se_null

  p_value <- switch(alternative,
    greater   = stats::pnorm(z_stat, lower.tail = FALSE),
    less      = stats::pnorm(z_stat, lower.tail = TRUE),
    two.sided = 2 * stats::pnorm(abs(z_stat), lower.tail = FALSE)
  )

  alpha  <- 1 - conf_level
  reject <- p_value < alpha

  # ---- CI on true prevalence ----
  # `ci_method` controls the interval shape; the interval is matched to
  # `alternative` so it points the same way as the test: for "greater" a
  # lower one-sided interval [L, 1], for "less" an upper one-sided interval
  # [0, U], for "two.sided" the usual two-sided interval. Independent of the
  # hypothesis test above (which uses the null-variance z at `threshold`),
  # so `reject` and "threshold outside the CI" agree in direction but not
  # exactly at the boundary -- the CI uses each method's own SE at p_hat.
  alpha_lo <- switch(alternative,
    two.sided = alpha / 2, greater = alpha, less = 0)
  alpha_hi <- switch(alternative,
    two.sided = alpha / 2, greater = 0,     less = alpha)
  # z multiplier for the closed side(s); NA on an open side.
  z_lo <- if (alpha_lo > 0) stats::qnorm(1 - alpha_lo) else NA_real_
  z_hi <- if (alpha_hi > 0) stats::qnorm(1 - alpha_hi) else NA_real_
  z_ci <- max(z_lo, z_hi, na.rm = TRUE)   # for the AC pseudo-count

  if (ci_method == "wald") {
    se_est    <- sqrt(p_hat * (1 - p_hat) / n_eff_adj)
    ci_lo_app <- if (is.na(z_lo)) 0 else max(p_hat - z_lo * se_est, 0)
    ci_hi_app <- if (is.na(z_hi)) 1 else min(p_hat + z_hi * se_est, 1)

  } else if (ci_method == "clopper-pearson") {
    x_eff     <- p_hat * n_eff_adj   # effective successes (continuous)
    ci_lo_app <- if (is.na(z_lo) || p_hat == 0) 0 else
      max(stats::qbeta(alpha_lo,     x_eff,     n_eff_adj - x_eff + 1), 0)
    ci_hi_app <- if (is.na(z_hi) || p_hat == 1) 1 else
      min(stats::qbeta(1 - alpha_hi, x_eff + 1, n_eff_adj - x_eff), 1)

  } else {   # agresti-coull
    x_eff     <- p_hat * n_eff_adj
    n_tilde   <- n_eff_adj + z_ci^2
    p_tilde   <- (x_eff + z_ci^2 / 2) / n_tilde
    se_tilde  <- sqrt(p_tilde * (1 - p_tilde) / n_tilde)
    ci_lo_app <- if (is.na(z_lo)) 0 else max(p_tilde - z_lo * se_tilde, 0)
    ci_hi_app <- if (is.na(z_hi)) 1 else min(p_tilde + z_hi * se_tilde, 1)
  }

  rg <- function(p) .rogan_gladen(p, sensitivity, specificity)

  prevalence <- max(0, min(1, rg(p_hat)))
  ci_lower   <- max(0, min(1, rg(ci_lo_app)))
  ci_upper   <- max(0, min(1, rg(ci_hi_app)))

  list(
    statistic   = z_stat,
    p_value     = p_value,
    reject      = reject,
    threshold   = threshold,
    alternative = alternative,
    prevalence  = prevalence,
    ci_lower    = ci_lower,
    ci_upper    = ci_upper,
    ci_method   = ci_method,
    n_total     = n_total,
    n_eff       = n_eff,
    n_eff_adj   = n_eff_adj,
    conf_level  = conf_level,
    sensitivity = sensitivity,
    specificity = specificity,
    icc_used    = icc_used,
    deff        = deff,
    fpc_N       = fpc_N
  )
}
