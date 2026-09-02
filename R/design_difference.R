#' Calculate sample size to detect a difference in prevalence between two groups
#'
#' @description
#' Design function: given the expected true prevalence in two groups, returns
#' the minimum \strong{per-group} sample size needed to detect that
#' difference with the target power. This is the design complement to
#' \code{test_difference()}: it answers "how many samples per group do I need
#' to run \code{test_difference()} with adequate power?"
#'
#' Equal allocation between groups is assumed. The same three design modes as
#' \code{design_precision()} / \code{design_threshold()} are supported, and
#' clustering / FPC are applied identically to \emph{each} group:
#' \describe{
#'   \item{SRS (no clustering args)}{Returns per-group \code{n} only.}
#'   \item{Fixed cluster size (\code{n_per_site} supplied)}{Solves for the
#'     number of clusters needed per group.}
#'   \item{Fixed number of clusters (\code{n_sites} supplied)}{Solves for
#'     the target samples per cluster via closed form.}
#' }
#'
#' @param prevalence1 Numeric in (0, 1). Expected true prevalence in group 1.
#' @param prevalence2 Numeric in (0, 1). Expected true prevalence in group 2.
#' @param power Target power (1 - beta); default 0.80.
#' @param alternative Direction of the alternative hypothesis; must match the
#'   \code{alternative} used in the companion \code{test_difference()} call.
#'   \code{"two.sided"} (default) tests Ha: p1 != p2; \code{"greater"} tests
#'   Ha: p1 > p2 (requires \code{prevalence1 > prevalence2}); \code{"less"}
#'   tests Ha: p1 < p2 (requires \code{prevalence1 < prevalence2}).
#' @param sensitivity Diagnostic sensitivity in (0, 1]; default 1. Assumed
#'   equal across groups.
#' @param specificity Diagnostic specificity in (0, 1]; default 1. Assumed
#'   equal across groups.
#' @param conf_level Significance level expressed as a confidence level:
#'   alpha = 1 - conf_level. Default 0.95 (alpha = 0.05).
#' @param n_sites Optional positive integer. Fix the number of clusters per
#'   group; the function solves for the target samples per cluster
#'   (\code{n_per_site}). Cannot be used with \code{n_per_site}.
#' @param n_per_site Optional positive integer. Fix the samples per cluster;
#'   the function solves for the number of clusters needed per group
#'   (\code{n_sites}). Cannot be used with \code{n_sites}.
#' @param icc Numeric in \[0, 1\]. Intra-cluster correlation; default 0
#'   (SRS). Requires one of \code{n_sites} or \code{n_per_site} when
#'   \code{icc > 0}.
#' @param fpc_N Optional positive integer. Population size \strong{per group}
#'   for a finite-population correction. \code{NULL} (default) = no FPC.
#'
#' @details
#' \strong{Power formula.} The per-group SRS sample size for a two-proportion
#' z-test with equal allocation is
#'
#' \deqn{n = \frac{\left(z_\alpha \sqrt{2\,\bar{p}(1-\bar{p})} +
#'   z_\beta \sqrt{p_{1,app}(1-p_{1,app}) + p_{2,app}(1-p_{2,app})}\right)^2}
#'   {(p_{1,app} - p_{2,app})^2}}
#'
#' where \eqn{p_{g,app} = p_g \cdot Se + (1-p_g)(1-Sp)} is each group's
#' prevalence on the apparent scale (Rogan-Gladen forward transform),
#' \eqn{\bar{p} = (p_{1,app} + p_{2,app})/2}, \eqn{z_\beta = qnorm(power)},
#' and \eqn{z_\alpha = qnorm(1 - alpha)} for a one-sided test or
#' \eqn{qnorm(1 - alpha/2)} for two-sided.
#'
#' Because \eqn{p_{1,app} - p_{2,app} = (p_1 - p_2)(Se + Sp - 1)}, an
#' imperfect test shrinks the detectable gap and inflates \code{n}. When
#' \eqn{Se = Sp = 1} the formula reduces to the standard two-proportion
#' result.
#'
#' Clustering and the FPC are then applied to each group with the same
#' closed-form logic as \code{design_threshold()} / \code{design_precision()}.
#'
#' @section Equations and sources:
#' Building blocks from the MMS-SD Study Design Workshop
#' (\url{https://mrc-ide.github.io/MMS-SD_workshop/}); the two-group
#' extension is \strong{not} covered there.
#' \itemize{
#'   \item \emph{One-sample power / sample-size form}
#'     \eqn{n = (z_{1-\beta} + z_{1-\alpha/2})^2\, p(1-p)/(p-p_0)^2} --
#'     Module 4 "Statistical power", slide "Sample size formulae"
#'     (lecture slides p. 7). This function uses the two-proportion
#'     generalisation (pooled + unpooled variance terms, equal
#'     allocation), which the workshop does not derive; standard form,
#'     e.g. Fleiss, Levin & Paik (2003), ch. 2.
#'   \item \emph{Design effect} \eqn{D_{eff} = 1 + (\bar n - 1)\,r} --
#'     Module 5 "Dealing with over-dispersion in multi-cluster studies",
#'     slide "Why is the ICC useful?" (p. 5); applied per group.
#'   \item \emph{Apparent-prevalence (Rogan-Gladen) scale}
#'     \eqn{p_{app} = p\,Se + (1-p)(1-Sp)} -- not in the workshop;
#'     Rogan & Gladen (1978).
#'   \item \emph{Finite-population correction} -- not in the workshop;
#'     Cochran (1977), \emph{Sampling Techniques}.
#' }
#'
#' @return A named list. Always present:
#'   \item{n_per_group}{Per-group sample size required.}
#'   \item{n_total}{Total sample size, both groups (\code{2 * n_per_group}).}
#'   \item{n_eff}{SRS-equivalent independent per-group sample size (the base
#'     power-formula n, before the design effect and FPC). Defined the same
#'     way as \code{n_eff} in \code{design_precision()} /
#'     \code{design_threshold()}.}
#'   \item{prevalence1}{Group 1 true prevalence as supplied.}
#'   \item{prevalence2}{Group 2 true prevalence as supplied.}
#'   \item{apparent_prev1}{Group 1 prevalence on the apparent scale.}
#'   \item{apparent_prev2}{Group 2 prevalence on the apparent scale.}
#'   \item{delta_app}{\code{apparent_prev1 - apparent_prev2}.}
#'   \item{power}{Target power as supplied.}
#'   \item{alternative}{Alternative hypothesis as supplied.}
#'   \item{conf_level}{Confidence level as supplied.}
#'   \item{sensitivity}{Sensitivity as supplied.}
#'   \item{specificity}{Specificity as supplied.}
#'   \item{icc}{ICC as supplied.}
#'   \item{deff}{Design effect applied (same for both groups).}
#'   \item{fpc_N}{\code{fpc_N} as supplied, or \code{NULL}.}
#'
#'   Present when clustering is specified:
#'   \item{n_sites}{If \code{n_per_site} supplied: clusters required per
#'     group. If \code{n_sites} supplied: echoed back. \code{NULL} for SRS.}
#'   \item{n_per_site}{If \code{n_sites} supplied: target samples per
#'     cluster. If \code{n_per_site} supplied: echoed back. \code{NULL} for
#'     SRS.}
#'
#' @references
#' MMS-SD Study Design Workshop, Modules 3 (hypothesis testing), 4 (power /
#' sample-size formulae) and 5 (ICC / design effect).
#' \url{https://mrc-ide.github.io/MMS-SD_workshop/}
#'
#' Fleiss JL, Levin B, Paik MC (2003). Statistical Methods for Rates and
#' Proportions, 3rd ed. Wiley. (Two-sample proportion sample size.)
#'
#' Rogan WJ, Gladen B (1978). Estimating prevalence from the results of a
#' screening test. American Journal of Epidemiology 107(1):71-76.
#'
#' Cochran WG (1977). Sampling Techniques, 3rd ed. Wiley.
#' (Finite-population correction.)
#'
#' @seealso \code{\link{test_difference}}
#'
#' @export
#'
#' @examples
#' # Detect a 10 vs 20 percent prevalence difference at 80% power
#' design_difference(prevalence1 = 0.10, prevalence2 = 0.20)
#'
#' # Higher power needs more samples
#' design_difference(0.10, 0.20, power = 0.90)
#'
#' # Imperfect test widens the requirement
#' design_difference(0.10, 0.20, sensitivity = 0.9, specificity = 0.95)
#'
#' # Fixed cluster size: how many sites per group?
#' design_difference(0.10, 0.20, n_per_site = 30, icc = 0.02)
design_difference <- function(prevalence1,
                              prevalence2,
                              power       = 0.80,
                              alternative = "two.sided",
                              sensitivity = 1,
                              specificity = 1,
                              conf_level  = 0.95,
                              n_sites     = NULL,
                              n_per_site  = NULL,
                              icc         = 0,
                              fpc_N       = NULL) {

  # ---- validate scalars ----
  if (length(prevalence1) != 1 || !is.numeric(prevalence1))
    stop("`prevalence1` must be a single number in (0, 1) (got class `",
         class(prevalence1)[1], "`, length ", length(prevalence1), ").")
  if (!is.finite(prevalence1) || prevalence1 <= 0 || prevalence1 >= 1)
    stop("`prevalence1` must be strictly between 0 and 1 (got ", prevalence1, ").")
  if (length(prevalence2) != 1 || !is.numeric(prevalence2))
    stop("`prevalence2` must be a single number in (0, 1) (got class `",
         class(prevalence2)[1], "`, length ", length(prevalence2), ").")
  if (!is.finite(prevalence2) || prevalence2 <= 0 || prevalence2 >= 1)
    stop("`prevalence2` must be strictly between 0 and 1 (got ", prevalence2, ").")

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
         "monotone in power.")

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
  if (!is.finite(sensitivity) || sensitivity <= 0 || sensitivity > 1)
    stop("`sensitivity` must be in (0, 1] (got ", sensitivity, ").")
  if (!is.finite(specificity) || specificity <= 0 || specificity > 1)
    stop("`specificity` must be in (0, 1] (got ", specificity, ").")
  correction <- sensitivity + specificity - 1
  if (correction <= 0)
    stop("sensitivity + specificity must exceed 1 for the Rogan-Gladen correction ",
         "(got ", sensitivity, " + ", specificity, " = ", sensitivity + specificity, ").")
  if (correction < 0.1)
    warning("sensitivity + specificity = ", round(sensitivity + specificity, 4),
            " is very close to 1. The Rogan-Gladen adjustment is numerically unstable.")

  if (length(conf_level) != 1 || !is.numeric(conf_level))
    stop("`conf_level` must be a single number (got ",
         if (!is.numeric(conf_level)) paste0("class `", class(conf_level)[1], "`")
         else paste0("length = ", length(conf_level)), ").")
  if (!is.finite(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("`conf_level` must be strictly between 0 and 1 (got ", conf_level, ").")

  if (length(icc) != 1 || !is.numeric(icc))
    stop("`icc` must be a single number (got ",
         if (!is.numeric(icc)) paste0("class `", class(icc)[1], "`")
         else paste0("length ", length(icc)), "). Use 0 for an unclustered (SRS) design.")
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
      !is.finite(fpc_N) || fpc_N < 1 || fpc_N != floor(fpc_N)))
    stop("`fpc_N` must be a single finite positive integer (got ",
         if (length(fpc_N) != 1) paste0("length = ", length(fpc_N)) else fpc_N, ").")
  if (icc > 0 && is.null(n_sites) && is.null(n_per_site))
    stop("icc > 0 requires a cluster structure: supply `n_sites` or `n_per_site`.")

  # ---- direction check ----
  if (alternative == "greater" && prevalence1 <= prevalence2)
    stop("`prevalence1` (", prevalence1, ") must be > `prevalence2` (",
         prevalence2, ") when alternative = 'greater'.")
  if (alternative == "less" && prevalence1 >= prevalence2)
    stop("`prevalence1` (", prevalence1, ") must be < `prevalence2` (",
         prevalence2, ") when alternative = 'less'.")
  if (alternative == "two.sided" && abs(prevalence1 - prevalence2) < 1e-8)
    stop("`prevalence1` (", prevalence1, ") and `prevalence2` (", prevalence2,
         ") are effectively equal; there is no difference to power for.")

  # ---- apparent scale (Rogan-Gladen forward transform) ----
  p1_app <- .apparent_prev(prevalence1, sensitivity, specificity)
  p2_app <- .apparent_prev(prevalence2, sensitivity, specificity)
  pbar   <- (p1_app + p2_app) / 2

  delta <- abs(p1_app - p2_app)
  if (delta < 1e-8)
    stop("`prevalence1` (", prevalence1, ") and `prevalence2` (", prevalence2,
         ") map to effectively identical apparent prevalences ",
         "(|difference| = ", signif(delta, 3), "), often because an imperfect ",
         "test compresses them together. The required sample size explodes ",
         "toward infinity -- move the prevalences further apart, or use a more ",
         "accurate test.")

  # ---- SRS per-group power formula ----
  alpha <- 1 - conf_level
  z_a   <- if (alternative == "two.sided") stats::qnorm(1 - alpha / 2) else
                                           stats::qnorm(1 - alpha)
  z_b   <- stats::qnorm(power)

  se0 <- sqrt(2 * pbar * (1 - pbar))
  se1 <- sqrt(p1_app * (1 - p1_app) + p2_app * (1 - p2_app))

  n_base_cont <- ((z_a * se0 + z_b * se1) / delta)^2

  # ---- design effect (same closed-form logic as design_threshold) ----
  # Applied identically to each group. icc == 0 covers every unclustered
  # case: an earlier guard errored if icc > 0 without a cluster structure.
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
        "Target power of %.0f%% is unachievable with %d sites per group and ICC = %.3f.\n",
        "Increase `n_sites`, lower the ICC assumption, or relax the power target."
      ), 100 * power, n_sites, icc))
    }
    n_cont <- n_base_cont * n_sites * (1 - icc) / denom
    deff   <- n_cont / n_base_cont

    if (deff <= 1) {
      stop("n_sites = ", n_sites, " is >= the per-group SRS sample size (~",
           ceiling(n_base_cont), "), so each site would receive < 1 person on ",
           "average -- not a valid cluster design. ",
           "Use n_sites < ", ceiling(n_base_cont), ", or supply `n_per_site` instead.")
    }
  }

  # ---- finite-population correction (per group) ----
  if (!is.null(fpc_N)) {
    n_cont <- (n_cont * fpc_N) / (n_cont + fpc_N - 1)
  }

  n_per_group <- ceiling(n_cont)
  n_total     <- 2L * n_per_group

  # n_eff: SRS-equivalent independent per-group sample size (base
  # power-formula n before the design effect and FPC).
  n_eff <- ceiling(n_base_cont)

  # ---- distribute across sites ----
  if (!is.null(n_per_site)) {
    n_sites_out    <- ceiling(n_per_group / n_per_site)
    n_per_site_out <- n_per_site
  } else if (!is.null(n_sites)) {
    n_per_site_out <- ceiling(n_per_group / n_sites)
    n_sites_out    <- n_sites
  } else {
    n_sites_out    <- NULL
    n_per_site_out <- NULL
  }

  list(
    n_per_group    = n_per_group,
    n_total        = n_total,
    n_eff          = n_eff,
    n_sites        = n_sites_out,
    n_per_site     = n_per_site_out,
    prevalence1    = prevalence1,
    prevalence2    = prevalence2,
    apparent_prev1 = p1_app,
    apparent_prev2 = p2_app,
    delta_app      = p1_app - p2_app,
    power          = power,
    alternative    = alternative,
    conf_level     = conf_level,
    sensitivity    = sensitivity,
    specificity    = specificity,
    icc            = icc,
    deff           = deff,
    fpc_N          = fpc_N
  )
}
