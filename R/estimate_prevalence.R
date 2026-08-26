#' Estimate prevalence and its precision
#'
#' Analysis function: given observed data, what is the prevalence and how
#' precise is our estimate? Implements the classical Wald confidence
#' interval (Module 1), generalized to account for clustering via the
#' design effect (Module 5): p_hat +/- z * sqrt(p_hat*(1-p_hat) /
#' (n * Deff)).
#'
#' If `icc` is not supplied, it is estimated directly from the data
#' as the ratio of observed-to-expected variance across clusters (Module
#' 5's design-effect worked example), rather than assumed to be 0 --
#' assuming independence by default would understate uncertainty for any
#' genuinely clustered survey.
#'
#' For a Bayesian treatment of the same problem (posterior distribution
#' over prevalence and ICC jointly, rather than a plug-in point estimate
#' of ICC), see `DRpower::get_prevalence()` -- not wrapped here,
#' since the two approaches answer the question differently and
#' shouldn't be silently mixed behind one interface.
#'
#' @param x Integer vector of positive counts per cluster/site.
#' @param n Integer vector of total samples per cluster/site (same length as x).
#' @param icc Optional. Intra-cluster correlation coefficient. If `NULL`,
#'   ICC is estimated from the data (Module 5 worked-example method). Set
#'   explicitly to `0` to force the simple-random-sampling case with no
#'   clustering adjustment.
#' @param fpc_N Optional. Total population size, for a finite-population
#'   correction. `NULL` = no FPC applied.
#' @param conf_level Confidence level. Defaults to 0.95.
#' @param sensitivity Diagnostic sensitivity in (0, 1]; default 1 (perfect
#'   test). When less than 1, the Rogan-Gladen correction is applied to
#'   convert apparent prevalence to true prevalence. See Details.
#' @param specificity Diagnostic specificity in (0, 1]; default 1 (perfect
#'   test).
#'
#' @details
#' **Rogan-Gladen correction** (applied when `sensitivity < 1` or
#' `specificity < 1`): an imperfect test inflates apparent prevalence via
#' false positives and deflates it via false negatives. The correction
#' unscrambles these two biases:
#'
#' \deqn{p_{\text{true}} = \frac{p_{\text{apparent}} - (1 - \text{specificity})}
#'       {\text{sensitivity} + \text{specificity} - 1}}
#'
#' The correction is a linear transform, so it is applied directly to the
#' Wald CI endpoints as well as the point estimate. For most PCR-based MMS
#' assays, sensitivity and specificity are effectively 1 and no correction
#' is needed.
#'
#' @return A list with:
#'   \item{prevalence}{Point estimate of true prevalence}
#'   \item{ci_lower}{Lower confidence limit}
#'   \item{ci_upper}{Upper confidence limit}
#'   \item{margin_of_error}{Half-width of CI on the true-prevalence scale}
#'   \item{icc_used}{ICC applied (estimated from data or supplied)}
#'   \item{deff}{Design effect applied}
#'   \item{n_total}{Total samples}
#'   \item{n_eff}{Effective sample size (n_total / deff)}
#'
#' @export
#'
#' @examples
#' # Single site / simple random sample (no clustering)
#' estimate_prevalence(x = 8, n = 50)
#'
#' # Multi-site, clustered data -- ICC estimated from the data
#' estimate_prevalence(
#'   x = c(0, 4, 0, 22, 25, 16, 12, 8),
#'   n = c(60, 80, 70, 100, 40, 60, 50, 90)
#' )
#'
#' # Imperfect diagnostic test
#' estimate_prevalence(x = 30, n = 100, sensitivity = 0.9, specificity = 0.95)
estimate_prevalence <- function(x,
                                n,
                                icc         = NULL,
                                fpc_N       = NULL,
                                conf_level  = 0.95,
                                sensitivity = 1,
                                specificity = 1) {

  # Check for NA/NaN/Inf before any comparisons -- otherwise R throws a
  # generic "missing value where TRUE/FALSE needed" with no context.
  if (is.logical(x) || is.logical(n))
    stop("`x` and `n` must be numeric, not logical. ",
         "Did you accidentally pass TRUE/FALSE instead of counts?")
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
         "(got length(x) = ", length(x), ", length(n) = ", length(n), "). ",
         "Each element of `x` is the positive count for one cluster and each ",
         "element of `n` is that cluster's total.")
  if (any(x != floor(x)))
    stop("`x` must contain whole numbers — counts cannot be fractional ",
         "(found x[", which(x != floor(x))[1], "] = ", x[which(x != floor(x))[1]], ").")
  if (any(n != floor(n)))
    stop("`n` must contain whole numbers — sample sizes cannot be fractional ",
         "(found n[", which(n != floor(n))[1], "] = ", n[which(n != floor(n))[1]], ").")
  if (any(n <= 0))
    stop("`n` must be positive for every cluster ",
         "(found n[", which(n <= 0)[1], "] = ", n[which(n <= 0)[1]], "). ",
         "A cluster with zero or negative total is undefined.")
  if (any(x < 0))
    stop("`x` must be non-negative ",
         "(found x[", which(x < 0)[1], "] = ", x[which(x < 0)[1]], "). ",
         "Counts cannot be negative.")
  if (any(x > n))
    stop("`x` cannot exceed `n` ",
         "(found x[", which(x > n)[1], "] = ", x[which(x > n)[1]],
         " > n[", which(x > n)[1], "] = ", n[which(x > n)[1]], "). ",
         "Positive counts cannot exceed the total tested per cluster.")
  # Scalar-length checks: catch vectors and NULL before comparisons fire R's
  # generic "condition has length > 1" or "argument is of length zero" errors.
  if (is.null(sensitivity) || length(sensitivity) != 1)
    stop("`sensitivity` must be a single number in (0, 1] (got ",
         if (is.null(sensitivity)) "NULL" else paste0("length = ", length(sensitivity)), "). ",
         "Leave at the default (1) for a perfect test.")
  if (is.null(specificity) || length(specificity) != 1)
    stop("`specificity` must be a single number in (0, 1] (got ",
         if (is.null(specificity)) "NULL" else paste0("length = ", length(specificity)), "). ",
         "Leave at the default (1) for a perfect test.")
  if (length(conf_level) != 1)
    stop("`conf_level` must be a single number (got length = ", length(conf_level), ").")
  if (!is.null(icc) && length(icc) != 1)
    stop("`icc` must be a single number or NULL (got length = ", length(icc), "). ",
         "Set `icc = NULL` to estimate ICC from the data.")
  if (!is.null(fpc_N) && length(fpc_N) != 1)
    stop("`fpc_N` must be a single number or NULL (got length = ", length(fpc_N), "). ",
         "Set `fpc_N = NULL` to skip the finite-population correction.")

  if (!is.finite(conf_level))
    stop("`conf_level` must be a single finite number (got ", conf_level, ").")
  if (conf_level <= 0 || conf_level >= 1)
    stop("`conf_level` must be strictly between 0 and 1 (got ", conf_level, "). ",
         "Use, e.g., 0.95 for a 95% confidence interval.")
  if (!is.null(icc) && !is.finite(icc))
    stop("`icc` must be a finite number (got ", icc, "). ",
         "To estimate ICC from the data, leave `icc = NULL`.")
  if (!is.null(icc) && (icc < 0 || icc > 1))
    stop("`icc` must be in [0, 1] (got ", icc, "). ",
         "Values outside this range imply negative within-cluster variance, ",
         "which is not possible. To estimate ICC from the data, leave `icc = NULL`.")
  if (!is.null(fpc_N) && (!is.numeric(fpc_N) || !is.finite(fpc_N) || fpc_N <= 0))
    stop("`fpc_N` must be a finite positive number representing the total population size ",
         "(got ", fpc_N, "). Set `fpc_N = NULL` to skip the finite-population correction.")
  if (!is.finite(sensitivity))
    stop("`sensitivity` must be a finite number (got ", sensitivity, ").")
  if (sensitivity <= 0 || sensitivity > 1)
    stop("`sensitivity` must be in (0, 1] (got ", sensitivity, "). ",
         "A sensitivity of 0 means the test never detects true positives.")
  if (!is.finite(specificity))
    stop("`specificity` must be a finite number (got ", specificity, ").")
  if (specificity <= 0 || specificity > 1)
    stop("`specificity` must be in (0, 1] (got ", specificity, "). ",
         "A specificity of 0 means the test always returns a false positive.")
  correction <- sensitivity + specificity - 1
  if (correction <= 0)
    stop("sensitivity + specificity must exceed 1 for the Rogan-Gladen correction ",
         "(got ", sensitivity, " + ", specificity, " = ", sensitivity + specificity, "). ",
         "With se + sp <= 1 the test performs at or below chance and true prevalence ",
         "is not identifiable from apparent prevalence.")
  if (correction < 0.1)
    warning("sensitivity + specificity = ", round(sensitivity + specificity, 4),
            " is very close to 1 (correction = ", round(correction, 4), "). ",
            "The Rogan-Gladen adjustment is numerically unstable -- ",
            "prevalence estimates and CI will be unreliable.")

  n_clusters <- length(n)
  n_total    <- sum(n)

  if (!is.null(fpc_N) && fpc_N < n_total)
    stop("`fpc_N` = ", fpc_N, " is less than the total sample size = ", n_total, ". ",
         "The population must be at least as large as the sample. ",
         "Check your inputs, or set `fpc_N = NULL` to skip the FPC.")
  if (!is.null(fpc_N) && fpc_N == n_total)
    warning("`fpc_N` equals the total sample size (", n_total, "): you have surveyed the ",
            "entire population, so sampling variance is zero and the CI collapses to a point. ",
            "Set `fpc_N = NULL` if this is not intended.")
  p_hat      <- sum(x) / n_total   # apparent prevalence

  # -----------------------------------------------------------------
  # Design effect / ICC (Module 5)
  #
  #   Deff    = Var_obs / Var_SRS
  #   Var_SRS = mean(p_hat*(1-p_hat) / n_i)   expected under independence
  #   Var_obs = sample variance of per-cluster prevalences
  #   ICC     = (Deff - 1) / (n_bar - 1)
  # -----------------------------------------------------------------
  n_bar <- mean(n)

  if (is.null(icc)) {
    if (n_clusters < 2 || n_bar == 1) {
      # Can't estimate ICC: single cluster, or every cluster has exactly 1
      # observation (Kish formula has n_bar-1 in the denominator → div/0).
      icc_used <- 0
      deff     <- 1
    } else {
      p_i     <- x / n
      var_obs <- stats::var(p_i)
      var_srs <- mean(p_hat * (1 - p_hat) / n)

      deff <- if (var_srs > 0) var_obs / var_srs else 1
      deff <- max(deff, 1)   # Deff < 1 not meaningful here

      icc_used <- (deff - 1) / (n_bar - 1)
      icc_used <- min(max(icc_used, 0), 1)
      deff     <- 1 + (n_bar - 1) * icc_used  # keep pair mutually consistent
    }
  } else {
    icc_used <- icc
    deff     <- 1 + (n_bar - 1) * icc_used
  }

  n_eff <- n_total / deff

  # Finite-population correction
  fpc <- if (!is.null(fpc_N)) {
    sqrt((fpc_N - n_eff) / (fpc_N - 1))
  } else {
    1
  }

  # -----------------------------------------------------------------
  # Wald interval on apparent prevalence (Module 1 + Module 5)
  # -----------------------------------------------------------------
  z            <- stats::qnorm(1 - (1 - conf_level) / 2)
  se           <- sqrt(p_hat * (1 - p_hat) / n_eff) * fpc
  moe_apparent <- z * se

  ci_lo_app <- max(p_hat - moe_apparent, 0)
  ci_hi_app <- min(p_hat + moe_apparent, 1)

  # -----------------------------------------------------------------
  # Rogan-Gladen correction: apparent -> true prevalence
  # Identity transform when sensitivity = specificity = 1.
  # -----------------------------------------------------------------
  rg <- function(p) (p - (1 - specificity)) / correction

  prevalence <- max(0, min(1, rg(p_hat)))
  ci_lower   <- max(0, min(1, rg(ci_lo_app)))
  ci_upper   <- max(0, min(1, rg(ci_hi_app)))

  # Derive moe from the actual returned CI, not the theoretical formula.
  # When clamping occurs (prevalence near 0 or 1), the CI endpoints are
  # truncated and moe_apparent/correction would be inconsistent with them.
  moe <- (ci_upper - ci_lower) / 2

  list(
    prevalence      = prevalence,
    ci_lower        = ci_lower,
    ci_upper        = ci_upper,
    margin_of_error = moe,
    icc_used        = icc_used,
    deff            = deff,
    n_total         = n_total,
    n_eff           = n_eff
  )
}
