#' Calculate sample size for a target margin of error on prevalence
#'
#' @description
#' Given a target precision (margin of error), returns the minimum total
#' sample size required. Handles imperfect diagnostic tests via the
#' Rogan-Gladen variance adjustment, clustered sampling via the design
#' effect (Kish formula), and finite-population corrections.
#'
#' **Three design modes** -- controlled by `n_sites` and `n_per_site`:
#' \describe{
#'   \item{SRS (`n_sites = NULL`, `n_per_site = NULL`)}{Treats all
#'     observations as independent. Returns total `n` only; `n_sites` and
#'     `n_per_site` are both `NULL` in the output.}
#'   \item{Fixed cluster size (`n_per_site` supplied)}{Deff is determined
#'     directly from the cluster size; solves for the required number of
#'     clusters (`n_sites` in the output).}
#'   \item{Fixed number of clusters (`n_sites` supplied)}{Resolves the Deff/n
#'     circularity via closed-form algebra; returns target samples per
#'     cluster (`n_per_site` in the output).}
#' }
#'
#' @param prevalence Numeric in (0, 1). Expected true prevalence.
#' @param moe Numeric in (0, 0.5). Target margin of error (half-width of
#'   confidence interval) on the true-prevalence scale.
#' @param sensitivity Diagnostic sensitivity in (0, 1]; default 1 (perfect
#'   test). Set below 1 to activate the Rogan-Gladen variance adjustment.
#' @param specificity Diagnostic specificity in (0, 1]; default 1.
#' @param conf_level Confidence level; default 0.95.
#' @param n_sites Optional positive integer. Fix the number of clusters.
#'   The function solves for the required samples per cluster and returns it
#'   as `n_per_site` in the output.
#'   Cannot be used together with `n_per_site`.
#' @param n_per_site Optional positive integer. Fix the samples per cluster.
#'   The function computes Deff directly, then solves for the number of
#'   clusters needed and returns it as `n_sites` in the output.
#'   Cannot be used together with `n_sites`.
#' @param icc Numeric in \[0, 1\]. Intra-cluster correlation; default 0 (SRS).
#'   If `icc > 0`, supply exactly one of `n_sites` or `n_per_site` --
#'   without a cluster structure, Deff is not computable.
#' @param fpc_N Optional positive integer. Total population size, for a
#'   finite-population correction. Reduces the required `n` when the sample
#'   is a non-trivial fraction of the population. `NULL` (default) = no FPC.
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
#' **Fixed number of clusters -- resolving the circularity**: when `n_sites`
#' is fixed, Deff depends on average cluster size \eqn{\bar{n} = n/n_{sites}},
#' which depends on n. Substituting and solving gives the closed-form:
#'
#' \deqn{n = \frac{n_0 \cdot n_{sites} \cdot (1-ICC)}{n_{sites} - n_0 \cdot ICC}}
#'
#' where \eqn{n_0} is the SRS sample size. A solution exists only when
#' \eqn{n_{sites} > n_0 \cdot ICC}; if not, the target MOE is unachievable.
#' The function stops with the minimum achievable MOE for that site count.
#'
#' When `n_per_site` is fixed instead, Deff is non-circular (cluster size is
#' known directly) and the number of sites follows from
#' \code{ceiling(n / n_per_site)}.
#'
#' @return A named list. The following fields are always present:
#'   \item{n}{Total sample size required (ceiling of the continuous solution)}
#'   \item{n_eff}{Effective independent sample size: \code{n / deff}. Accounts
#'     for clustering and the finite-population correction. This is the
#'     equivalent number of independent (SRS) observations your design is worth.}
#'   \item{apparent_prev}{Apparent (observed-test) prevalence implied by
#'     \code{prevalence}, \code{sensitivity}, and \code{specificity}}
#'   \item{moe}{Target MOE (as supplied)}
#'   \item{conf_level}{Confidence level (as supplied)}
#'   \item{sensitivity}{Sensitivity (as supplied)}
#'   \item{specificity}{Specificity (as supplied)}
#'   \item{icc}{ICC (as supplied; 0 for SRS)}
#'   \item{deff}{Design effect applied: 1 for SRS, > 1 for clustered designs}
#'   \item{fpc_N}{\code{fpc_N} as supplied, or \code{NULL}}
#'
#'   The following fields depend on the design mode:
#'   \item{n_sites}{If `n_per_site` was supplied: clusters required
#'     (\code{ceiling(n / n_per_site)}). If `n_sites` was supplied: echoed
#'     back. \code{NULL} for SRS.}
#'   \item{n_per_site}{If `n_sites` was supplied: target samples per cluster
#'     (\code{ceiling(n / n_sites)}). If `n_per_site` was supplied: echoed
#'     back. \code{NULL} for SRS. Note: this is the minimum whole-number
#'     cluster size needed -- actual allocation may differ if your real cluster
#'     sizes vary.}
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
#' # Fixed number of sites: what is the target per-site sample?
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
  # Reject vectors and non-numeric types: all parameters are numeric scalars.
  # A vector causes R's generic "the condition has length > 1" error deep in
  # an if(); a character/logical value slips past the length check and then
  # either crashes in is.finite() or (for logical) coerces silently.
  if (length(prevalence) != 1 || !is.numeric(prevalence))
    stop("`prevalence` must be a single number (got ",
         if (!is.numeric(prevalence)) paste0("class `", class(prevalence)[1], "`")
         else paste0("length ", length(prevalence)), ").")
  if (length(moe) != 1 || !is.numeric(moe))
    stop("`moe` must be a single number (got ",
         if (!is.numeric(moe)) paste0("class `", class(moe)[1], "`")
         else paste0("length ", length(moe)), ").")
  if (length(icc) != 1 || !is.numeric(icc))
    stop("`icc` must be a single number (got ",
         if (!is.numeric(icc)) paste0("class `", class(icc)[1], "`")
         else paste0("length ", length(icc)), "). Use 0 for an unclustered (SRS) design.")
  if (length(conf_level) != 1 || !is.numeric(conf_level))
    stop("`conf_level` must be a single number (got ",
         if (!is.numeric(conf_level)) paste0("class `", class(conf_level)[1], "`")
         else paste0("length ", length(conf_level)), ").")
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
      (!is.numeric(n_sites) || length(n_sites) != 1 || !is.finite(n_sites) || n_sites != floor(n_sites) || n_sites < 1))
    stop("`n_sites` must be a single finite positive integer (got ",
         if (length(n_sites) != 1) paste0("length = ", length(n_sites)) else n_sites, "). ",
         "It represents the number of sampling clusters in your design.")
  if (!is.null(n_per_site) &&
      (!is.numeric(n_per_site) || length(n_per_site) != 1 || !is.finite(n_per_site) || n_per_site != floor(n_per_site) || n_per_site < 1))
    stop("`n_per_site` must be a single finite positive integer (got ",
         if (length(n_per_site) != 1) paste0("length = ", length(n_per_site)) else n_per_site, "). ",
         "It represents the fixed number of individuals sampled per cluster.")
  if (!is.null(fpc_N) && (!is.numeric(fpc_N) || length(fpc_N) != 1 || !is.finite(fpc_N) || fpc_N <= 0))
    stop("`fpc_N` must be a single finite positive number representing total population size (got ",
         if (length(fpc_N) != 1) paste0("length = ", length(fpc_N)) else fpc_N,
         "). Set `fpc_N = NULL` to skip the finite-population correction.")
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
      stop("n_sites = ", n_sites, " is >= the SRS sample size (n_base ~= ",
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

  # n_eff: effective independent sample size after clustering and FPC.
  # Equivalent to the number of independent (SRS) observations this
  # design is worth. Always <= n_total; equals n_total only for SRS with no FPC.
  n_eff <- n_total / deff

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
