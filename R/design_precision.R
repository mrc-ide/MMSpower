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
  if (prevalence <= 0 || prevalence >= 1)
    stop("`prevalence` must be in (0, 1)")
  if (moe <= 0 || moe >= 0.5)
    stop("`moe` must be in (0, 0.5)")
  if (sensitivity <= 0 || sensitivity > 1)
    stop("`sensitivity` must be in (0, 1]")
  if (specificity <= 0 || specificity > 1)
    stop("`specificity` must be in (0, 1]")
  correction <- sensitivity + specificity - 1
  if (correction <= 0)
    stop("`sensitivity` + `specificity` must exceed 1")
  if (icc < 0 || icc > 1)
    stop("`icc` must be in [0, 1]")
  if (!is.null(n_sites) && !is.null(n_per_site))
    stop("Supply at most one of `n_sites` or `n_per_site`, not both")

  # icc > 0 with no cluster structure is underspecified:
  # Deff = 1 + (n_bar - 1)*icc requires knowing n_bar, which requires
  # knowing the cluster structure. Rather than silently assume anything,
  # we stop and ask the caller to be explicit.
  if (icc > 0 && is.null(n_sites) && is.null(n_per_site))
    stop(
      "`icc` > 0 but cluster structure not specified.\n",
      "Supply `n_sites` or `n_per_site` to compute the design effect.\n",
      "To assume simple random sampling (no clustering), set `icc = 0`."
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
