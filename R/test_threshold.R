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
#' The reported CI uses the Wald method at \code{conf_level} (two-sided
#' regardless of \code{alternative}), consistent with
#' \code{estimate_prevalence()}.
#'
#' @return A named list:
#'   \item{statistic}{z-statistic (on the apparent-prevalence scale)}
#'   \item{p_value}{p-value for the chosen \code{alternative}}
#'   \item{reject}{Logical: \code{TRUE} if \code{p_value < 1 - conf_level}}
#'   \item{threshold}{Decision threshold as supplied}
#'   \item{alternative}{Alternative hypothesis as supplied}
#'   \item{prevalence}{Rogan-Gladen corrected point estimate}
#'   \item{ci_lower}{Lower bound of the two-sided Wald CI}
#'   \item{ci_upper}{Upper bound of the two-sided Wald CI}
#'   \item{n_total}{Total samples}
#'   \item{n_eff}{Effective independent sample size: \code{n_total / deff}}
#'   \item{conf_level}{Confidence level (as supplied)}
#'   \item{sensitivity}{Sensitivity (as supplied)}
#'   \item{specificity}{Specificity (as supplied)}
#'   \item{icc_used}{ICC applied: estimated from data or as supplied}
#'   \item{deff}{Design effect applied}
#'   \item{fpc_N}{\code{fpc_N} as supplied, or \code{NULL}}
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
                           fpc_N       = NULL) {

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
  if (is.null(sensitivity) || length(sensitivity) != 1)
    stop("`sensitivity` must be a single number in (0, 1] (got ",
         if (is.null(sensitivity)) "NULL" else paste0("length = ", length(sensitivity)), ").")
  if (is.null(specificity) || length(specificity) != 1)
    stop("`specificity` must be a single number in (0, 1] (got ",
         if (is.null(specificity)) "NULL" else paste0("length = ", length(specificity)), ").")
  if (length(conf_level) != 1)
    stop("`conf_level` must be a single number (got length = ", length(conf_level), ").")
  if (!is.null(icc) && length(icc) != 1)
    stop("`icc` must be a single number or NULL (got length = ", length(icc), ").")
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

  n_clusters <- length(n)
  n_total    <- sum(n)

  if (!is.null(fpc_N) && fpc_N < n_total)
    stop("`fpc_N` = ", fpc_N, " is less than the total sample size = ", n_total, ". ",
         "The population must be at least as large as the sample.")
  if (!is.null(fpc_N) && fpc_N == n_total)
    warning("`fpc_N` equals the total sample size: you have surveyed the entire ",
            "population, so the test statistic is undefined (zero variance). ",
            "Set `fpc_N = NULL` if this is not intended.")

  p_hat <- sum(x) / n_total   # apparent prevalence

  # ---- design effect / ICC ----
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
  } else {
    icc_used <- icc
    deff     <- 1 + (n_bar - 1) * icc_used
  }

  n_eff <- n_total / deff
  fpc   <- if (!is.null(fpc_N)) sqrt((fpc_N - n_eff) / (fpc_N - 1)) else 1

  # Variance-equivalent simple sample size (same as estimate_prevalence)
  n_eff_adj <- n_eff / (fpc^2)

  # ---- threshold on apparent scale ----
  theta_app <- threshold * sensitivity + (1 - threshold) * (1 - specificity)

  # ---- z-statistic using null variance ----
  # Null SE evaluated at theta_app, not p_hat, for correct Type I error control.
  se_null <- sqrt(theta_app * (1 - theta_app) / n_eff_adj)

  if (se_null == 0)
    stop("Null standard error is zero (theta_app = ", theta_app,
         ", n_eff_adj = ", round(n_eff_adj, 2), "). ",
         "This occurs when threshold is 0 or 1, or when n_eff_adj is effectively infinite.")

  z_stat <- (p_hat - theta_app) / se_null

  p_value <- switch(alternative,
    greater   = stats::pnorm(z_stat, lower.tail = FALSE),
    less      = stats::pnorm(z_stat, lower.tail = TRUE),
    two.sided = 2 * stats::pnorm(abs(z_stat), lower.tail = FALSE)
  )

  alpha  <- 1 - conf_level
  reject <- p_value < alpha

  # ---- two-sided Wald CI on true prevalence (same as estimate_prevalence) ----
  z_ci      <- stats::qnorm(1 - alpha / 2)
  se_est    <- sqrt(p_hat * (1 - p_hat) / n_eff_adj)
  ci_lo_app <- max(p_hat - z_ci * se_est, 0)
  ci_hi_app <- min(p_hat + z_ci * se_est, 1)

  rg <- function(p) (p - (1 - specificity)) / correction

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
    n_total     = n_total,
    n_eff       = n_eff,
    conf_level  = conf_level,
    sensitivity = sensitivity,
    specificity = specificity,
    icc_used    = icc_used,
    deff        = deff,
    fpc_N       = fpc_N
  )
}
