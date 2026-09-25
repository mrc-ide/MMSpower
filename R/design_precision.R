#' Calculate sample size for a target margin of error when estimating prevalence
#'
#' @description
#' Given a target precision (margin of error), returns the minimum total
#' sample size required. Handles imperfect diagnostic tests via the
#' Rogan-Gladen variance adjustment, clustered sampling via the design
#' effect (Kish formula), and finite-population corrections via Cochran's
#' (1977) adjustment.
#'
#' **Three design modes** -- controlled by `n_sites` and `n_per_site`:
#' \describe{
#'   \item{SRS (`n_sites = NULL`, `n_per_site = NULL`)}{Treats all
#'     observations as independent. Returns total `n` only; `n_sites` and
#'     `n_per_site` are both `NULL` in the output.}
#'   \item{Fixed cluster size (`n_per_site` supplied)}{Solves for the
#'     required number of clusters, returned as `n_sites` in the output.}
#'   \item{Fixed number of clusters (`n_sites` supplied)}{Solves for the
#'     target samples per cluster, returned as `n_per_site` in the output.}
#' }
#'
#' @param prevalence Numeric in (0, 1). Expected true prevalence.
#' @param moe Numeric in (0, 0.5). Target margin of error (half-width of
#'   confidence interval) on the true-prevalence scale.
#' @param sensitivity Diagnostic sensitivity in (0, 1]; default 1 (perfect
#'   test). Set below 1 to activate the Rogan-Gladen variance adjustment.
#' @param specificity Diagnostic specificity in (0, 1]; default 1.
#' @param conf_level Confidence level, in (0, 1); default 0.95.
#' @param n_sites Optional positive integer. Fix the number of clusters.
#'   The function solves for the required samples per cluster and returns it
#'   as `n_per_site` in the output.
#'   Cannot be used together with `n_per_site` -- see the "Known limitation"
#'   section below.
#' @param n_per_site Optional positive integer. Fix the samples per cluster.
#'   The function computes Deff directly, then solves for the number of
#'   clusters needed and returns it as `n_sites` in the output.
#'   Cannot be used together with `n_sites` -- see the "Known limitation"
#'   section below.
#' @param icc Numeric in \[0, 1\]. Intra-cluster correlation; default 0 (SRS).
#'   If `icc > 0`, supply exactly one of `n_sites` or `n_per_site` --
#'   without a cluster structure, Deff is not computable. An `icc` below
#'   `sqrt(.Machine$double.eps)` (about 1.5e-8) is treated as 0 (SRS), so a
#'   negligible upstream estimate does not force the clustered code path.
#' @param fpc_N Optional positive integer. Total population size, for a
#'   finite-population correction. Reduces the required `n` when the sample
#'   is a non-trivial fraction of the population. `NULL` (default) = no FPC.
#'
#' @details
#' **Rogan-Gladen variance adjustment**: the Rogan-Gladen transform is
#' linear in the apparent prevalence, so scaling its variance by the
#' (constant) RG slope is exact -- not a delta-method approximation. The
#' MOE formula uses that exact variance:
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
#' \eqn{n_{sites} > n_0 \cdot ICC}; if not, the target MOE is unachievable --
#' adding more samples per site inflates Deff proportionally, so MOE floors
#' at:
#'
#' \deqn{MOE_{min} = z_{1-\alpha/2}\sqrt{\frac{p_{app}(1-p_{app}) \cdot ICC}
#' {n_{sites} \cdot (Se + Sp - 1)^2}}}
#'
#' The function stops and reports this minimum achievable MOE for that site
#' count.
#'
#' When `n_per_site` is fixed instead, Deff is non-circular (cluster size is
#' known directly) and the number of sites follows from
#' \code{ceiling(n / n_per_site)}.
#'
#' **Finite-population correction**: sampling without replacement from a
#' finite population of size \eqn{N} shrinks the sampling variance by the
#' factor \eqn{(N-n)/(N-1)} relative to sampling from an infinite
#' population. Requiring the finite-population design to hit the same
#' target variance as the continuous sample size \eqn{n_{pre}} computed
#' above (after the Rogan-Gladen and design-effect inflation, before the
#' FPC) and solving for \eqn{n} gives Cochran's (1977) adjustment:
#'
#' \deqn{n = \frac{n_{pre}}{1 + (n_{pre} - 1)/N} = \frac{n_{pre} \cdot N}{N + n_{pre} - 1}}
#'
#' applied with \eqn{N} = `fpc_N`. As \eqn{N \to \infty} this converges to
#' \eqn{n_{pre}} (no correction); as \eqn{N} shrinks toward \eqn{n_{pre}},
#' it pulls the required sample down toward \eqn{N} -- you cannot sample
#' more people than exist in the population.
#'
#' @section Equations and sources:
#' Mostly direct workshop material (MMS-SD Study Design Workshop,
#' \url{https://mrc-ide.github.io/MMS-SD_workshop/}):
#' \itemize{
#'   \item \emph{MOE sample size} \eqn{n = z_{1-\alpha/2}^2\,p(1-p)/m^2} --
#'     Module 2 "Sample size calculation based on margin of error", slides
#'     "Deriving the sample size formula" (steps 1-3) and "Worked example"
#'     (lecture slides pp. 9-10).
#'   \item \emph{Margin of error / Wald interval}
#'     \eqn{\hat p \pm z_{1-\alpha/2}\sqrt{\hat p(1-\hat p)/n}} -- Module 1
#'     "Sampling from a population" (p. 11), Module 2 "Precision and sample
#'     size" (p. 9).
#'   \item \emph{Design effect and effective sample size}
#'     \eqn{D_{eff} = 1 + (\bar n - 1)\,r}, \eqn{N_{eff} = N/D_{eff}}, and
#'     the clustered Wald interval
#'     \eqn{\hat p \pm z_{1-\alpha/2}\sqrt{\hat p(1-\hat p)/N \cdot
#'     D_{eff}}} -- Module 5 "Dealing with over-dispersion", slides "The
#'     effective sample size", "Why is the ICC useful?" and "How can we
#'     design multi-cluster studies?" (pp. 4-5).
#'   \item \emph{Rogan-Gladen variance inflation} (the
#'     \eqn{/(Se + Sp - 1)^2} factor) -- \strong{not} in the workshop;
#'     Rogan & Gladen (1978).
#'   \item \emph{Fixed-\code{n_sites} closed-form circularity solve} and the
#'     \emph{finite-population correction} -- not in the workshop
#'     (package's own algebra; Cochran (1977) for the FPC).
#' }
#'
#' @section Known limitation -- n_sites and n_per_site cannot both be fixed:
#' \strong{This is a scope limitation, flagged for a team decision, not a
#' bug.} This function only solves forward for the required sample size
#' \code{n}, so supplying both a fixed number of clusters and a fixed
#' cluster size is rejected -- with both fixed, the total sample size is
#' already determined and there is nothing left to solve for. Someone who
#' already has both numbers fixed (e.g. "I'm running 50 sites of 20 people
#' each -- what MOE do I actually get?") needs the reverse calculation
#' instead: achieved MOE from a given \code{n_sites * n_per_site}. That mode
#' does not exist yet -- it would need to be added, not just unlocked by
#' relaxing this guard.
#'
#' @return A named list with the following fields, in the order returned:
#'   \item{n}{Total sample size required (ceiling of the continuous solution).
#'     With clusters, the total actually collected is
#'     \code{n_sites * n_per_site}, which can be slightly higher than
#'     \code{n} because the per-site size is rounded up; budget for that
#'     figure.}
#'   \item{n_eff}{SRS-equivalent independent sample size the design achieves:
#'     the number of independent observations needed to hit the same \code{moe}
#'     (equal to the base SRS sample size before the design effect and FPC).
#'     Clustering inflates the collected \code{n} above this, before any FPC
#'     adjustment. With no FPC \code{n_eff} \eqn{\le} \code{n}, with equality
#'     only for an unclustered design; \strong{FPC shrinks \code{n} but does
#'     not affect \code{n_eff}, so \code{n_eff} can then exceed \code{n}.}
#'     Rounded up (\code{ceiling()}) to match \code{n}'s rounding, so the
#'     unclustered-equality case holds exactly.}
#'   \item{n_sites}{If `n_per_site` was supplied: clusters required
#'     (\code{ceiling(n / n_per_site)}). If `n_sites` was supplied: echoed
#'     back. \code{NULL} for SRS.}
#'   \item{n_per_site}{If `n_sites` was supplied: target samples per cluster
#'     (\code{ceiling(n / n_sites)}). If `n_per_site` was supplied: echoed
#'     back. \code{NULL} for SRS. Note: this is the minimum whole-number
#'     cluster size needed -- actual allocation may differ if your real cluster
#'     sizes vary.}
#'   \item{prevalence}{Expected true prevalence (as supplied).}
#'   \item{apparent_prev}{Apparent (observed-test) prevalence implied by
#'     \code{prevalence}, \code{sensitivity}, and \code{specificity}}
#'   \item{moe}{Target MOE (as supplied)}
#'   \item{conf_level}{Confidence level (as supplied)}
#'   \item{sensitivity}{Sensitivity (as supplied)}
#'   \item{specificity}{Specificity (as supplied)}
#'   \item{icc}{ICC (as supplied; 0 for SRS)}
#'   \item{deff}{Design effect applied: 1 for SRS, > 1 for clustered
#'     designs. Reported as the effect of the rounded, \emph{fielded}
#'     design, \code{1 + (n_per_site - 1) * icc}. Because \code{n_per_site}
#'     is rounded up (and, with an FPC, the FPC acts as a separate
#'     variance factor), \code{n} is \strong{not} exactly
#'     \code{n_eff * deff} -- treat \code{n} as the headline figure and
#'     \code{deff} / \code{n_eff} as diagnostics. This assumes clusters
#'     come out equal-sized (\code{n_per_site} exactly, not an average)
#'     -- unequal realized cluster sizes in the field will inflate the
#'     true design effect beyond this estimate.}
#'   \item{fpc_N}{\code{fpc_N} as supplied, or \code{NULL}}
#'
#' @references
#' MMS-SD Study Design Workshop, Modules 1 (sampling), 2 (sample size from
#' margin of error) and 5 (ICC / design effect).
#' \url{https://mrc-ide.github.io/MMS-SD_workshop/}
#'
#' Rogan WJ, Gladen B (1978). Estimating prevalence from the results of a
#' screening test. American Journal of Epidemiology 107(1):71-76.
#'
#' Cochran WG (1977). Sampling Techniques, 3rd ed. Wiley.
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
  # Every parameter must be a single number. Without this check, a vector
  # fails later with a cryptic R error, and a logical (TRUE/FALSE) silently
  # coerces to 1/0 instead of failing loudly.
  if (length(prevalence) != 1 || !is.numeric(prevalence))
    stop("`prevalence` must be a single number in (0, 1) (got class `",
         class(prevalence)[1], "`, length ", length(prevalence), ").")
  if (length(moe) != 1 || !is.numeric(moe))
    stop("`moe` must be a single number in (0, 0.5) (got class `",
         class(moe)[1], "`, length ", length(moe), ").")
  if (length(sensitivity) != 1 || !is.numeric(sensitivity))
    stop("`sensitivity` must be a single number in (0, 1] (got class `",
         class(sensitivity)[1], "`, length ", length(sensitivity), "). ",
         "Note: `TRUE`/`FALSE` is logical, not numeric -- pass 1 for a perfect test.")
  if (length(specificity) != 1 || !is.numeric(specificity))
    stop("`specificity` must be a single number in (0, 1] (got class `",
         class(specificity)[1], "`, length ", length(specificity), "). ",
         "Note: `TRUE`/`FALSE` is logical, not numeric -- pass 1 for a perfect test.")
  if (length(conf_level) != 1 || !is.numeric(conf_level))
    stop("`conf_level` must be a single number in (0, 1) (got class `",
         class(conf_level)[1], "`, length ", length(conf_level), ").")
  if (length(icc) != 1 || !is.numeric(icc))
    stop("`icc` must be a single number in [0, 1] (got class `",
         class(icc)[1], "`, length ", length(icc), ").")

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
  if (!is.finite(conf_level))
    stop("`conf_level` must be a single finite number (got ", conf_level, ").")
  if (!is.finite(icc))
    stop("`icc` must be a single finite number (got ", icc, ").")

  if (prevalence < 0 || prevalence > 1)
    stop("`prevalence` must be in (0, 1) (got ", prevalence, "). ",
         "`prevalence` is a fraction of the population and cannot be ",
         if (prevalence < 0) "negative." else "greater than 1.")
  if (prevalence == 0 || prevalence == 1)
    stop("`prevalence` must be in (0, 1) (got ", prevalence, "). ",
         "When the prevalence is 0 or 1, there is no uncertainty to estimate ",
         "and no meaningful sample size. Use a value from a pilot study, ",
         "historical data, or conservative guess.")
  if (moe <= 0)
    stop("`moe` must be a positive target margin of error (got ", moe, "). ",
         "`moe` is the target half-width of the confidence interval; a value ",
         "of 0 or less would demand infinite precision to achieve.")
  if (moe >= 0.5)
    stop("`moe` must be less than 0.5 (got ", moe, "). ",
         "A margin of error of 0.5 or more spans (or exceeds) the entire ",
         "(0, 1) prevalence range, so the target interval carries no ",
         "information. `moe` is a decimal fraction, not a percentage -- ",
         "e.g. use 0.05 to target a precision of +/-5 percentage points.")
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
  correction <- sensitivity + specificity - 1
  if (correction <= 0)
    stop("`sensitivity` + `specificity` must exceed 1 for the Rogan-Gladen correction ",
         "(got ", sensitivity, " + ", specificity, " = ", sensitivity + specificity, "). ",
         "A test with `sensitivity` + `specificity` of 1 or less performs no better ",
         "than random guessing, so its results cannot be corrected into a ",
         "reliable prevalence estimate.")
  if (correction < 0.1)
    warning("`sensitivity` + `specificity` = ", round(sensitivity + specificity, 4),
            " is only ", round(correction, 4), " above the required minimum of 1. ",
            "The Rogan-Gladen adjustment is numerically unstable here -- ",
            "required n will be extremely large and results unreliable.")
  if (icc < 0 || icc > 1)
    stop("`icc` must be in [0, 1] (got ", icc, "). ",
         "`icc` is a correlation and cannot be ",
         if (icc < 0) "negative." else "greater than 1.")
  if (conf_level <= 0 || conf_level >= 1)
    stop("`conf_level` must be in (0, 1) (got ", conf_level, "). ",
         "For example, use 0.95 for a 95% confidence interval.")
  # Both fixed would mean the total n is already determined, leaving nothing
  # to solve for -- see the "Known limitation" roxygen section: this rejects
  # a real use case (achieved MOE from a fixed n_sites * n_per_site) that
  # would need a new reverse mode, not just relaxing this guard.
  if (!is.null(n_sites) && !is.null(n_per_site))
    stop("Supply at most one of `n_sites` or `n_per_site`, not both. ",
         "`n_sites` fixes the number of clusters and solves for samples per cluster; ",
         "`n_per_site` fixes the cluster size and solves for the number of clusters.")
  if (!is.null(n_sites) &&
      (!is.numeric(n_sites) || length(n_sites) != 1 || !is.finite(n_sites) || n_sites != floor(n_sites) || n_sites < 1))
    stop("`n_sites` must be a single finite positive integer (got ",
         if (!is.numeric(n_sites)) paste0("class `", class(n_sites)[1], "`")
         else if (length(n_sites) != 1) paste0("length = ", length(n_sites))
         else n_sites, "). ",
         "`n_sites` is the number of sampling clusters.")
  if (!is.null(n_per_site) &&
      (!is.numeric(n_per_site) || length(n_per_site) != 1 || !is.finite(n_per_site) || n_per_site != floor(n_per_site) || n_per_site < 1))
    stop("`n_per_site` must be a single finite positive integer (got ",
         if (!is.numeric(n_per_site)) paste0("class `", class(n_per_site)[1], "`")
         else if (length(n_per_site) != 1) paste0("length = ", length(n_per_site))
         else n_per_site, "). ",
         "`n_per_site` is the fixed number of individuals sampled per cluster.")
  if (!is.null(fpc_N) && (!is.numeric(fpc_N) || length(fpc_N) != 1 || !is.finite(fpc_N) ||
      fpc_N < 1 || fpc_N != floor(fpc_N)))
    stop("`fpc_N` must be a single finite positive integer (got ",
         if (!is.numeric(fpc_N)) paste0("class `", class(fpc_N)[1], "`")
         else if (length(fpc_N) != 1) paste0("length = ", length(fpc_N))
         else fpc_N,
         "). `fpc_N` is the total population size. ",
         "Set `fpc_N = NULL` to skip the finite-population correction.")

  # A cluster structure cannot be larger than the population it is drawn from.
  if (!is.null(fpc_N) && !is.null(n_sites) && n_sites > fpc_N)
    stop("`n_sites` (", n_sites, ") exceeds `fpc_N` (", fpc_N, "): you cannot ",
         "have more clusters than individuals in the population.")
  if (!is.null(fpc_N) && !is.null(n_per_site) && n_per_site > fpc_N)
    stop("`n_per_site` (", n_per_site, ") exceeds `fpc_N` (", fpc_N, "): a ",
         "cluster cannot be larger than the whole population.")

  # Uses a fuzzy zero-threshold instead of an exact icc == 0 comparison, so
  # an upstream estimate like 1e-12 (floating-point noise, not a real signal)
  # is still treated as SRS rather than forcing the clustered code path.
  icc_is_zero <- icc < sqrt(.Machine$double.eps)
  if (!icc_is_zero && is.null(n_sites) && is.null(n_per_site))
    stop("icc > 0 requires a cluster structure to compute the design effect. ",
         "Supply `n_sites` (fix the number of clusters) or `n_per_site` (fix the ",
         "cluster size), or set `icc = 0` for an unclustered (SRS) design."
    )

  # ---- base sample size (SRS, apparent-prevalence scale) ----
  z           <- stats::qnorm(1 - (1 - conf_level) / 2)
  p_app       <- .apparent_prev(prevalence, sensitivity, specificity)
  n_base_cont <- z^2 * p_app * (1 - p_app) / (moe^2 * correction^2)

  # ---- design effect and total n ----
  # This branch handles n_sites and n_per_site both NULL. If icc > 0 with no
  # cluster structure, the earlier guard would already have rejected it --
  # so reaching here with both NULL guarantees icc_is_zero.
  if (icc_is_zero) {
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
    #   min_moe = z * sqrt(p_app*(1-p_app)*icc / (n_sites * correction^2))
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

    # deff must be > 1 when icc > 0 and we have multiple sites. deff <= 1
    # means n_cont <= n_sites, i.e., average cluster size <= 1 -- a
    # physically impossible design. This happens when n_sites >= n_base,
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

  # n_eff: the SRS-equivalent independent sample size this design achieves,
  # i.e. the number of independent observations needed to hit the same MOE.
  # That is exactly n_base_cont -- clustering inflates the collected `n`
  # above it, before any FPC adjustment. With no FPC, n_eff <= n_total.
  # FPC shrinks n_total but does not affect n_eff, so n_eff can then
  # exceed n_total.
  n_eff <- ceiling(n_base_cont)

  # ---- distribute across sites ----
  if (!is.null(n_per_site)) {
    n_sites_out    <- ceiling(n_total / n_per_site)
    n_per_site_out <- n_per_site
  } else if (!is.null(n_sites)) {
    n_per_site_out <- ceiling(n_total / n_sites)
    n_sites_out    <- n_sites
    # deff was computed from the closed-form solve's continuous (pre-rounding,
    # pre-FPC) n. But as n_per_site_out is rounded up here, and n_total may
    # have shrunk because of FPC -- deff no longer matches the actual
    # per-site design being reported. It is recomputed from n_per_site_out
    # so the output stays self-consistent (deff == 1 + (n_per_site - 1) *
    # icc always holds).
    if (!icc_is_zero) deff <- 1 + (n_per_site_out - 1) * icc
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
