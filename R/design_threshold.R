#' Calculate sample size to detect a prevalence above (or below) a threshold
#'
#' @description
#' Design function: given a decision threshold and an alternative prevalence
#' to detect, returns the minimum sample size needed to achieve the target
#' statistical power. This is the direct design complement to
#' \code{test_threshold()}: it answers "how many samples do I need to run
#' \code{test_threshold()} with adequate power?"
#'
#' The same three design modes as \code{design_precision()} are supported:
#' \describe{
#'   \item{SRS (no clustering args)}{Returns total \code{n} only.}
#'   \item{Fixed cluster size (\code{n_per_site} supplied)}{Solves for
#'     number of clusters needed.}
#'   \item{Fixed number of clusters (\code{n_sites} supplied)}{Solves
#'     for target samples per cluster via closed-form.}
#' }
#'
#' @param threshold Numeric in (0, 1). The decision threshold on the
#'   true-prevalence scale (the boundary under H0 in \code{test_threshold()}).
#' @param prevalence Numeric in (0, 1). The alternative true prevalence to
#'   detect. Must be > \code{threshold} for \code{alternative = "greater"},
#'   < \code{threshold} for \code{"less"}.
#' @param power Target power (1 - beta); default 0.80.
#' @param alternative Direction of the alternative hypothesis; must match the
#'   \code{alternative} used in the companion \code{test_threshold()} call:
#'   \code{"greater"} (default), \code{"less"}, or \code{"two.sided"}.
#' @param sensitivity Diagnostic sensitivity in (0, 1]; default 1.
#' @param specificity Diagnostic specificity in (0, 1]; default 1.
#' @param conf_level Significance level expressed as a confidence level:
#'   alpha = 1 - conf_level. Default 0.95 (alpha = 0.05).
#' @param n_sites Optional positive integer. Fix the number of clusters;
#'   the function solves for the target samples per cluster (\code{n_per_site}).
#'   Cannot be used with \code{n_per_site}.
#' @param n_per_site Optional positive integer. Fix the samples per cluster;
#'   the function solves for clusters needed (\code{n_sites}).
#'   Cannot be used with \code{n_sites}.
#' @param icc Numeric in \[0, 1\]. Intra-cluster correlation; default 0 (SRS).
#'   Requires one of \code{n_sites} or \code{n_per_site} when \code{icc > 0}.
#' @param fpc_N Optional positive integer. Total population size for a
#'   finite-population correction. \code{NULL} (default) = no FPC.
#'
#' @details
#' **Power formula**: the required SRS sample size for a one-proportion
#' z-test is:
#'
#' \deqn{n = \left(\frac{z_\alpha \sqrt{\theta_{app}(1-\theta_{app})} +
#'   z_\beta \sqrt{p_{app}(1-p_{app})}}
#'   {|p_{app} - \theta_{app}|}\right)^2}
#'
#' where \eqn{\theta_{app}} and \eqn{p_{app}} are the threshold and
#' alternative prevalence on the apparent scale via Rogan-Gladen:
#' \eqn{p_{app} = p \cdot Se + (1-p)(1-Sp)}. For a one-sided test,
#' \eqn{z_\alpha = qnorm(1 - alpha)}; for two-sided,
#' \eqn{z_\alpha = qnorm(1 - alpha/2)}.
#'
#' Clustering and FPC are applied using the same closed-form logic as
#' \code{design_precision()}.
#'
#' @return A named list. Always present:
#'   \item{n}{Total sample size required}
#'   \item{n_eff}{SRS-equivalent independent sample size the design achieves
#'     (the base power-formula n before the design effect and FPC). Clustering
#'     inflates the collected \code{n} above it; the FPC lets it fall below.
#'     Defined the same way as \code{n_eff} in \code{design_precision()} and
#'     \code{estimate_prevalence()}.}
#'   \item{threshold}{Decision threshold as supplied}
#'   \item{prevalence}{Alternative prevalence as supplied}
#'   \item{apparent_prev}{Alternative prevalence on the apparent scale}
#'   \item{threshold_app}{Threshold on the apparent scale}
#'   \item{power}{Target power as supplied}
#'   \item{alternative}{Alternative hypothesis as supplied}
#'   \item{conf_level}{Confidence level as supplied}
#'   \item{sensitivity}{Sensitivity as supplied}
#'   \item{specificity}{Specificity as supplied}
#'   \item{icc}{ICC as supplied}
#'   \item{deff}{Design effect applied}
#'   \item{fpc_N}{\code{fpc_N} as supplied, or \code{NULL}}
#'
#'   Present when clustering is specified:
#'   \item{n_sites}{If \code{n_per_site} supplied: clusters required.
#'     If \code{n_sites} supplied: echoed back. \code{NULL} for SRS.}
#'   \item{n_per_site}{If \code{n_sites} supplied: target samples per cluster.
#'     If \code{n_per_site} supplied: echoed back. \code{NULL} for SRS.}
#'
#' @export
#'
#' @examples
#' # How many samples to detect p=0.10 when threshold=0.05, at 80% power?
#' design_threshold(threshold = 0.05, prevalence = 0.10)
#'
#' # Higher power requires more samples
#' design_threshold(threshold = 0.05, prevalence = 0.10, power = 0.90)
#'
#' # Fixed cluster size: how many sites needed?
#' design_threshold(threshold = 0.05, prevalence = 0.10,
#'                  n_per_site = 20, icc = 0.05)
design_threshold <- function(threshold,
                              prevalence,
                              power       = 0.80,
                              alternative = "greater",
                              sensitivity = 1,
                              specificity = 1,
                              conf_level  = 0.95,
                              n_sites     = NULL,
                              n_per_site  = NULL,
                              icc         = 0,
                              fpc_N       = NULL) {

  # ---- validate scalars ----
  if (length(threshold) != 1 || !is.numeric(threshold))
    stop("`threshold` must be a single number in (0, 1) (got class `",
         class(threshold)[1], "`, length ", length(threshold), ").")
  if (!is.finite(threshold) || threshold <= 0 || threshold >= 1)
    stop("`threshold` must be strictly between 0 and 1 (got ", threshold, ").")
  if (length(prevalence) != 1 || !is.numeric(prevalence))
    stop("`prevalence` must be a single number in (0, 1) (got class `",
         class(prevalence)[1], "`, length ", length(prevalence), ").")
  if (!is.finite(prevalence) || prevalence <= 0 || prevalence >= 1)
    stop("`prevalence` must be strictly between 0 and 1 (got ", prevalence, ").")
  if (length(power) != 1 || !is.numeric(power))
    stop("`power` must be a single number (got ",
         if (!is.numeric(power)) paste0("class `", class(power)[1], "`")
         else paste0("length ", length(power)), ").")
  if (!is.finite(power) || power <= 0 || power >= 1)
    stop("`power` must be strictly between 0 and 1 (got ", power, "). ",
         "Use, e.g., 0.80 for 80% power.")
  if (power < 0.5)
    stop("`power` must be at least 0.5 (got ", power, "). ",
         "Below 0.5 the normal-approximation sample-size formula breaks down: ",
         "z_beta = qnorm(power) turns negative and the required n is no longer ",
         "monotone in power. A design with < 50% power is not meaningful anyway.")
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
  if (length(icc) != 1 || !is.numeric(icc))
    stop("`icc` must be a single number (got ",
         if (!is.numeric(icc)) paste0("class `", class(icc)[1], "`")
         else paste0("length ", length(icc)), "). Use 0 for an unclustered (SRS) design.")

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
         "(got ", sensitivity, " + ", specificity, " = ", sensitivity + specificity, ").")
  if (correction < 0.1)
    warning("sensitivity + specificity = ", round(sensitivity + specificity, 4),
            " is very close to 1. The Rogan-Gladen adjustment is numerically unstable.")
  if (!is.finite(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("`conf_level` must be strictly between 0 and 1 (got ", conf_level, ").")
  if (!is.finite(icc) || icc < 0 || icc > 1)
    stop("`icc` must be in [0, 1] (got ", icc, ").")
  if (!is.null(n_sites) && !is.null(n_per_site))
    stop("Supply at most one of `n_sites` or `n_per_site`, not both.")
  if (!is.null(n_sites) &&
      (!is.numeric(n_sites) || length(n_sites) != 1 || !is.finite(n_sites) ||
       n_sites != floor(n_sites) || n_sites < 1))
    stop("`n_sites` must be a single finite positive integer (got ",
         if (length(n_sites) != 1) paste0("length = ", length(n_sites)) else n_sites, ").")
  if (!is.null(n_per_site) &&
      (!is.numeric(n_per_site) || length(n_per_site) != 1 || !is.finite(n_per_site) ||
       n_per_site != floor(n_per_site) || n_per_site < 1))
    stop("`n_per_site` must be a single finite positive integer (got ",
         if (length(n_per_site) != 1) paste0("length = ", length(n_per_site)) else n_per_site, ").")
  if (!is.null(fpc_N) && (!is.numeric(fpc_N) || length(fpc_N) != 1 ||
      !is.finite(fpc_N) || fpc_N <= 0))
    stop("`fpc_N` must be a single finite positive number (got ",
         if (length(fpc_N) != 1) paste0("length = ", length(fpc_N)) else fpc_N, ").")
  if (icc > 0 && is.null(n_sites) && is.null(n_per_site))
    stop("icc > 0 requires a cluster structure: supply `n_sites` or `n_per_site`.")

  # ---- direction check ----
  if (alternative == "greater" && prevalence <= threshold)
    stop("`prevalence` (", prevalence, ") must be > `threshold` (", threshold,
         ") when alternative = 'greater'.")
  if (alternative == "less" && prevalence >= threshold)
    stop("`prevalence` (", prevalence, ") must be < `threshold` (", threshold,
         ") when alternative = 'less'.")
  if (alternative == "two.sided" &&
      abs(prevalence - threshold) < 1e-8 * max(1, abs(threshold)))
    stop("`prevalence` (", prevalence, ") must differ meaningfully from ",
         "`threshold` (", threshold, ") for alternative = 'two.sided'.")

  # ---- apparent scale (Rogan-Gladen forward transform) ----
  theta_app <- .apparent_prev(threshold,  sensitivity, specificity)
  p1_app    <- .apparent_prev(prevalence, sensitivity, specificity)

  delta <- abs(p1_app - theta_app)
  if (delta < 1e-8)
    stop("`threshold` (", threshold, ") and `prevalence` (", prevalence,
         ") map to effectively identical apparent prevalences ",
         "(|difference| = ", signif(delta, 3), "), often because an imperfect ",
         "test compresses them together. The required sample size explodes ",
         "toward infinity -- move `prevalence` further from `threshold`, or use ",
         "a more accurate test.")

  # ---- SRS power formula ----
  alpha  <- 1 - conf_level
  z_a    <- if (alternative == "two.sided") stats::qnorm(1 - alpha / 2) else
                                            stats::qnorm(1 - alpha)
  z_b    <- stats::qnorm(power)

  se0    <- sqrt(theta_app * (1 - theta_app))
  se1    <- sqrt(p1_app    * (1 - p1_app))

  n_base_cont <- ((z_a * se0 + z_b * se1) / delta)^2

  # ---- design effect (same closed-form logic as design_precision) ----
  # icc == 0 covers every unclustered case: an earlier guard errored if
  # icc > 0 without a cluster structure.
  if (icc == 0) {
    deff   <- 1
    n_cont <- n_base_cont

  } else if (!is.null(n_per_site)) {
    deff   <- 1 + (n_per_site - 1) * icc
    n_cont <- n_base_cont * deff

  } else {
    denom <- n_sites - n_base_cont * icc
    if (denom <= 0) {
      stop(sprintf(paste0(
        "Target power of %.0f%% is unachievable with %d sites and ICC = %.3f.\n",
        "Increase `n_sites`, lower the ICC assumption, or relax the power target."
      ), 100 * power, n_sites, icc))
    }
    n_cont <- n_base_cont * n_sites * (1 - icc) / denom
    deff   <- n_cont / n_base_cont

    if (deff <= 1) {
      stop("n_sites = ", n_sites, " is >= the SRS sample size (~",
           ceiling(n_base_cont), "), so each site would receive < 1 person on ",
           "average -- not a valid cluster design. ",
           "Use n_sites < ", ceiling(n_base_cont), ", or supply `n_per_site` instead.")
    }
  }

  # ---- finite-population correction ----
  if (!is.null(fpc_N)) {
    n_cont <- (n_cont * fpc_N) / (n_cont + fpc_N - 1)
  }

  n_total <- ceiling(n_cont)

  # n_eff: SRS-equivalent independent sample size the design achieves (the
  # base power-formula n before the design effect and FPC). Clustering
  # inflates the collected `n` above it; the FPC lets it fall below.
  # Defined the same way as n_eff in design_precision()/estimate_prevalence().
  n_eff <- ceiling(n_base_cont)

  # ---- distribute across sites ----
  if (!is.null(n_per_site)) {
    n_sites_out    <- ceiling(n_total / n_per_site)
    n_per_site_out <- n_per_site
  } else if (!is.null(n_sites)) {
    n_per_site_out <- ceiling(n_total / n_sites)
    n_sites_out    <- n_sites
  } else {
    n_sites_out    <- NULL
    n_per_site_out <- NULL
  }

  list(
    n             = n_total,
    n_eff         = n_eff,
    n_sites       = n_sites_out,
    n_per_site    = n_per_site_out,
    threshold     = threshold,
    prevalence    = prevalence,
    apparent_prev = p1_app,
    threshold_app = theta_app,
    power         = power,
    alternative   = alternative,
    conf_level    = conf_level,
    sensitivity   = sensitivity,
    specificity   = specificity,
    icc           = icc,
    deff          = deff,
    fpc_N         = fpc_N
  )
}
