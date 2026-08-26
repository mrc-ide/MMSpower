#' Calculate sample size for a target margin of error on prevalence
#'
#' Design function: given a target precision, what sample size do I need?
#' Implements the standard MOE formula (Module 2), generalized to account
#' for diagnostic test accuracy (Rogan-Gladen) and clustering via the
#' design effect (Module 5).
#'
#' @param prevalence Numeric in (0, 1). Expected true prevalence.
#' @param moe Numeric in (0, 0.5). Target margin of error (half-width of
#'   confidence interval) on the true-prevalence scale.
#' @param sensitivity Diagnostic sensitivity in (0, 1]; default 1 (perfect
#'   test). Propagates into the Rogan-Gladen variance adjustment.
#' @param specificity Diagnostic specificity in (0, 1]; default 1.
#' @param conf_level Confidence level; default 0.95.
#' @param n_sites Optional integer. Fixed number of sites/clusters. When
#'   supplied with `icc > 0`, the n_sites/Deff circularity is resolved via
#'   closed-form (see Details). Cannot be used if `n_per_site` is also given.
#' @param n_per_site Optional integer. Fixed samples per site (cluster size).
#'   When supplied with `icc > 0`, Deff is non-circular (cluster size is
#'   known). Cannot be used if `n_sites` is also given.
#' @param icc Numeric in [0, 1]. Intra-cluster correlation; default 0 (SRS).
#'   If `icc > 0`, at least one of `n_sites` or `n_per_site` must be given
#'   -- without a cluster structure, Deff is not computable.
#' @param fpc_N Optional integer. Total population size, for a
#'   finite-population correction. When supplied, the required n is
#'   reduced by the FPC factor: n_adj = n * N / (n + N - 1). NULL (default)
#'   means no FPC applied.
#'
#' @details
#' **Rogan-Gladen variance adjustment**: the MOE formula uses the
#' delta-method variance of the Rogan-Gladen estimator:
#'
#' \deqn{n = \frac{z^2 \, p_{app}(1-p_{app}) \cdot D_{eff}}{MOE^2 \cdot (Se + Sp - 1)^2}}
#'
#' where \eqn{p_{app} = p \cdot Se + (1-p)(1-Sp)} is the apparent prevalence
#' implied by the expected true prevalence and the test characteristics.
#' When \eqn{Se = Sp = 1} this reduces to the standard MOE formula.
#'
#' **n_sites/Deff circularity**: when `n_sites` is fixed, Deff depends on
#' average cluster size \eqn{\bar{n} = n/n_{sites}}, which depends on n.
#' Substituting and solving gives the closed-form:
#'
#' \deqn{n = \frac{n_0 \cdot n_{sites} \cdot (1-ICC)}{n_{sites} - n_0 \cdot ICC}}
#'
#' where \eqn{n_0} is the SRS sample size. This has a solution only when
#' \eqn{n_{sites} > n_0 \cdot ICC}; if not, the target MOE is unachievable
#' regardless of samples per site (more samples increase Deff proportionally).
#' The function stops with the minimum achievable MOE for that site count.
#'
#' When `n_per_site` is fixed instead, Deff is determined directly
#' (non-circular) and the required number of sites follows from ceiling(n /
#' n_per_site). The two routes will give slightly different totals due to
#' integer rounding -- this is expected and documented.
#'
#' @return A list with:
#'   \item{n}{Total sample size required}
#'   \item{n_base}{Required n under SRS with a perfect test}
#'   \item{n_sites}{Sites needed (if `n_per_site` supplied) or as given}
#'   \item{n_per_site}{Samples per site (if `n_sites` supplied) or as given}
#'   \item{prevalence}{True prevalence assumed}
#'   \item{apparent_prev}{Apparent prevalence used in the formula}
#'   \item{moe}{Target MOE}
#'   \item{conf_level}{Confidence level}
#'   \item{sensitivity}{Sensitivity used}
#'   \item{specificity}{Specificity used}
#'   \item{icc}{ICC used}
#'   \item{deff}{Design effect applied}
#'   \item{fpc_N}{Population size used for FPC, or NULL}
#'
#' @export
#'
#' @examples
#' # Simple random sample, perfect test
#' design_precision(prevalence = 0.3, moe = 0.05)
#'
#' # Imperfect test (sensitivity 90%, specificity 95%)
#' design_precision(0.3, 0.05, sensitivity = 0.9, specificity = 0.95)
#'
#' # Fixed cluster size: how many sites do I need?
#' design_precision(0.3, 0.05, n_per_site = 10, icc = 0.05)
#'
#' # Fixed number of sites: how many samples per site?
#' design_precision(0.3, 0.05, n_sites = 50, icc = 0.05)
design_precision <- function(prevalence,
                              moe,
                              sensitivity = 1,
                              specificity = 1,
                              conf_level  = 0.95,
                              n_sites     = NULL,
                              n_per_site  = NULL,
                              icc         = 0,
                              fpc_N       = NULL) {

  # ---- validation ----
  # Reject vectors: all parameters are scalars -- a vector input causes R's
  # generic "the condition has length > 1" error deep in an if(), not ours.
  if (length(prevalence) != 1)
    stop("`prevalence` must be a single number, not a vector (got length ", length(prevalence), ").")
  if (length(moe) != 1)
    stop("`moe` must be a single number, not a vector (got length ", length(moe), ").")
  if (length(icc) != 1)
    stop("`icc` must be a single number, not a vector (got length ", length(icc), ").")
  if (length(conf_level) != 1)
    stop("`conf_level` must be a single number, not a vector (got length ", length(conf_level), ").")

  # Check for NA/NaN/Inf before any comparisons -- otherwise R throws a
  # generic "missing value where TRUE/FALSE needed" with no context.
  if (!is.finite(prevalence))
    stop("`prevalence` must be a single finite number (got ", prevalence, ").")
  if (!is.finite(moe))
    stop("`moe` must be a single finite number (got ", moe, ").")
  if (!is.finite(sensitivity))
    stop("`sensitivity` must be a single finite number (got ", sensitivity, ").")
  if (!is.finite(specificity))
    stop("`specificity` must be a single finite number (got ", specificity, ").")
  if (!is.finite(icc))
    stop("`icc` must be a single finite number (got ", icc, "). ",
         "Use 0 for an unclustered (SRS) design.")
  if (!is.finite(conf_level))
    stop("`conf_level` must be a single finite number (got ", conf_level, ").")

  if (prevalence <= 0 || prevalence >= 1)
    stop("`prevalence` must be strictly between 0 and 1 (got ", prevalence, "). ",
         "Use a value from a pilot study, historical data, or conservative guess.")
  if (moe <= 0 || moe >= 0.5)
    stop("`moe` must be in (0, 0.5) (got ", moe, "). ",
         "`moe` is the target half-width of the confidence interval, e.g. 0.05 for +/-5 pp.")
  if (sensitivity <= 0 || sensitivity > 1)
    stop("`sensitivity` must be in (0, 1] (got ", sensitivity, "). ",
         "A sensitivity of 0 means the test never detects true positives.")
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
            "The Rogan-Gladen adjustment is numerically unstable here -- ",
            "required n will be extremely large and results unreliable.")
  if (icc < 0 || icc > 1)
    stop("`icc` must be in [0, 1] (got ", icc, "). ",
         "Values outside this range imply negative within-cluster variance, which ",
         "is not possible. Use 0 for an unclustered (SRS) design.")
  if (!is.null(n_sites) && !is.null(n_per_site))
    stop("Supply at most one of `n_sites` or `n_per_site`, not both. ",
         "`n_sites` fixes the number of clusters and solves for samples per cluster; ",
         "`n_per_site` fixes the cluster size and solves for the number of clusters.")
  if (!is.null(n_sites) &&
      (!is.numeric(n_sites) || !is.finite(n_sites) || n_sites != floor(n_sites) || n_sites < 1))
    stop("`n_sites` must be a finite positive integer (got ", n_sites, "). ",
         "It represents the number of sampling clusters in your design.")
  if (!is.null(n_per_site) &&
      (!is.numeric(n_per_site) || !is.finite(n_per_site) || n_per_site != floor(n_per_site) || n_per_site < 1))
    stop("`n_per_site` must be a finite positive integer (got ", n_per_site, "). ",
         "It represents the fixed number of individuals sampled per cluster.")
  if (!is.null(fpc_N) && (!is.numeric(fpc_N) || !is.finite(fpc_N) || fpc_N <= 0))
    stop("`fpc_N` must be a finite positive number representing total population size ",
         "(got ", fpc_N, "). Set `fpc_N = NULL` to skip the finite-population correction.")
  if (conf_level <= 0 || conf_level >= 1)
    stop("`conf_level` must be strictly between 0 and 1 (got ", conf_level, "). ",
         "Use, e.g., 0.95 for a 95% confidence interval.")

  if (icc > 0 && is.null(n_sites) && is.null(n_per_site))
    stop("icc > 0 requires a cluster structure to compute the design effect ",
         "(Deff = 1 + (n_bar - 1) * icc, where n_bar = total n / n_sites). ",
         "Supply `n_sites` (fix the number of clusters) or `n_per_site` (fix the ",
         "cluster size), or set `icc = 0` for an unclustered (SRS) design."
    )

  # ---- base sample size (SRS, apparent-prevalence scale) ----
  z           <- stats::qnorm(1 - (1 - conf_level) / 2)
  p_app       <- prevalence * sensitivity + (1 - prevalence) * (1 - specificity)
  n_base_cont <- z^2 * p_app * (1 - p_app) / (moe^2 * correction^2)

  # ---- design effect and total n ----
  if (icc == 0 || (is.null(n_sites) && is.null(n_per_site))) {
    # SRS: no clustering adjustment needed
    deff   <- 1
    n_cont <- n_base_cont

  } else if (!is.null(n_per_site)) {
    # Cluster size fixed -- non-circular: Deff determined directly
    deff   <- 1 + (n_per_site - 1) * icc
    n_cont <- n_base_cont * deff

  } else {
    # Number of sites fixed -- circular: Deff depends on n, n depends on Deff.
    #
    # Substituting Deff = 1 + (n/n_sites - 1)*icc into n = n_base*Deff and
    # solving for n gives the closed-form:
    #   n = n_base * n_sites * (1 - icc) / (n_sites - n_base * icc)
    #
    # A solution exists only when n_sites > n_base * icc. If not, adding
    # more samples per site increases Deff proportionally, so MOE never
    # reaches the target -- it floors at:
    #   MOE_min = z * sqrt(p_app*(1-p_app)*icc / (n_sites * correction^2))
    denom <- n_sites - n_base_cont * icc
    if (denom <= 0) {
      min_moe <- z * sqrt(p_app * (1 - p_app) * icc /
                            (n_sites * correction^2))
      stop(sprintf(paste0(
        "Target MOE of %.1f%% is unachievable with %d sites and ICC = %.3f.\n",
        "Minimum achievable MOE with these settings: %.1f%%.\n",
        "Increase `n_sites`, lower the ICC assumption, or relax the MOE target."
      ), 100 * moe, n_sites, icc, 100 * min_moe))
    }
    n_cont <- n_base_cont * n_sites * (1 - icc) / denom
    deff   <- n_cont / n_base_cont

    # Sanity check: deff must be > 1 when icc > 0 and we have multiple sites.
    # deff <= 1 means n_cont <= n_sites, i.e., average cluster size < 1 --
    # a physically impossible design. This happens when n_sites >= n_base,
    # meaning you have more sites than you'd need people under SRS.
    if (deff <= 1) {
      stop("n_sites = ", n_sites, " is >= the SRS sample size (n_base ≈ ",
           ceiling(n_base_cont), "), so each site would receive < 1 person on ",
           "average -- not a valid cluster design. ",
           "Use n_sites < ", ceiling(n_base_cont), ", or supply `n_per_site` ",
           "to fix the cluster size and solve for the number of sites instead.")
    }
  }

  # ---- finite-population correction ----
  if (!is.null(fpc_N)) {
    n_cont <- (n_cont * fpc_N) / (n_cont + fpc_N - 1)
  }

  n_total <- ceiling(n_cont)

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
    n_base        = ceiling(n_base_cont),
    n_sites       = n_sites_out,
    n_per_site    = n_per_site_out,
    prevalence    = prevalence,
    apparent_prev = p_app,
    moe           = moe,
    conf_level    = conf_level,
    sensitivity   = sensitivity,
    specificity   = specificity,
    icc           = icc,
    deff          = deff,
    fpc_N         = fpc_N
  )
}
