#' Calculate sample size for a target margin of error when estimating prevalence
#'
#' @description
#' Given a target precision (margin of error), returns the minimum total
#' sample size required; or, given a full design, returns the margin of
#' error it gives. Handles imperfect diagnostic tests via the Rogan-Gladen
#' variance adjustment, and clustered sampling via the design effect (Kish
#' formula).
#'
#' **What is solved for** -- leave out the one quantity you want:
#' \describe{
#'   \item{`moe` only}{Returns the number of people, with no sites.}
#'   \item{`moe` and `n_per_site`}{Returns the number of sites.}
#'   \item{`moe` and `n_sites`}{Returns the number of people per site.}
#'   \item{A full design, without `moe`}{Returns the margin of error: give
#'     `n` (no sites), `n_sites` and `n_per_site`, or a vector `n_per_site`
#'     with one size per site (for sites of different sizes).}
#' }
#'
#' @param prevalence Expected true prevalence, in (0, 1).
#' @param moe Target margin of error (half-width of the confidence interval)
#'   on the true-prevalence scale, in (0, 0.5). A decimal fraction, e.g. 0.05
#'   for +/-5 percentage points. Leave out with a full design.
#' @param sensitivity Diagnostic sensitivity in (0, 1]; default 1 (perfect
#'   test).
#' @param specificity Diagnostic specificity in (0, 1]; default 1.
#'   `sensitivity + specificity` must exceed 1.
#' @param conf_level Confidence level, in (0, 1); default 0.95.
#' @param ci_method The interval the study will be analysed with, as in
#'   [estimate_prevalence()]: `"clopper-pearson"` (default, the same default
#'   as [estimate_prevalence()]), `"wald"` or `"agresti-coull"`. See
#'   Details.
#' @param n Optional whole number: the number of people in a full design
#'   with no sites. Cannot be used with `n_sites` or `n_per_site`.
#' @param n_sites Optional whole number. Fixes the number of sites; alone,
#'   the function solves for people per site.
#' @param n_per_site Optional. People per site: one whole number (alone, the
#'   function solves for the number of sites), or a vector with one size per
#'   site, a full design with sites of different sizes (then leave out
#'   `n_sites`).
#' @param icc Intra-cluster correlation, in \[0, 1\]; default 0 (no
#'   clustering). An `icc` above 0 needs `n_sites` or `n_per_site`.
#'
#' @details
#' **Margin of error and the interval method**: the margin of error is the
#' half-width of the interval [estimate_prevalence()] would report if the
#' measured prevalence came out exactly as expected, so a design sized here
#' gives back the target when analysed with the same `ci_method`. For
#' `"wald"` it is the usual closed form below. For `"clopper-pearson"` and
#' `"agresti-coull"` there is no formula: the unclustered size \eqn{n_0} is
#' found by a numerical search, and the clustering steps below then work as
#' for Wald. Clopper-Pearson is the widest of the three, so it needs the
#' most people (341 against Wald's 323 for 30 percent +/- 5); the gap is
#' largest at low prevalence.
#'
#' **Rogan-Gladen variance adjustment**: the Rogan-Gladen transform is
#' linear in the apparent prevalence, so scaling its variance by the
#' (constant) RG slope is exact -- not a delta-method approximation. The
#' unclustered Wald sample size is
#'
#' \deqn{n_0 = \frac{z^2 \, p_{app}(1-p_{app})}{MOE^2 \cdot (Se + Sp - 1)^2}}
#'
#' where \eqn{p_{app} = p \cdot Se + (1-p)(1-Sp)} is the apparent prevalence
#' implied by the expected true prevalence and the test characteristics.
#' When \eqn{Se = Sp = 1} this reduces to the standard MOE formula.
#'
#' **Fixed site size**: the design effect \eqn{D_{eff} = 1 + (n_{per site}
#' - 1) \cdot ICC} is known directly, so \eqn{n = n_0 \cdot D_{eff}} and the
#' number of sites is \code{ceiling(n / n_per_site)}.
#'
#' **Fixed number of sites**: Deff depends on the number of people per site,
#' \eqn{n / n_{sites}}, which depends on n. Substituting and solving gives
#' the closed form:
#'
#' \deqn{n = \frac{n_0 \cdot n_{sites} \cdot (1-ICC)}{n_{sites} - n_0 \cdot ICC}}
#'
#' A solution exists only when \eqn{n_{sites} > n_0 \cdot ICC}. Otherwise
#' the target MOE is unachievable: each site is worth at most \eqn{1/ICC}
#' independent people however many are sampled there, so the MOE floors at
#' the margin of error of \eqn{n_{sites} / ICC} independent people (for
#' Wald, \eqn{z\sqrt{p_{app}(1-p_{app}) ICC / n_{sites}} / (Se + Sp - 1)}),
#' and the function stops and reports this minimum.
#'
#' **Rounding**: the number of sites (or people per site) is rounded up, and
#' the returned `n` is the number actually sampled, `n_sites * n_per_site`,
#' so it can be slightly above the continuous solution. The returned `moe`
#' is that of the rounded design, so it can be slightly below the target.
#'
#' **Full design**: the margin of error is that of \eqn{n / D_{eff}}
#' independent people, for Wald
#' \eqn{z \sqrt{p_{app}(1-p_{app}) D_{eff} / n} / (Se + Sp - 1)}. For sites
#' of different sizes, \eqn{D_{eff} = 1 + (n_w - 1) ICC} with
#' \eqn{n_w = \sum s^2 / \sum s} the size-weighted mean site size, as in
#' [estimate_prevalence()].
#'
#' @return A list with:
#'   \item{n}{Total sample size. With sites, \code{n_sites * n_per_site}.}
#'   \item{n_sites}{Number of sites (\code{NULL} if unclustered).}
#'   \item{n_per_site}{People per site (\code{NULL} if unclustered).}
#'   \item{deff}{Design effect of the returned design,
#'     \code{1 + (n_per_site - 1) * icc}; 1 if unclustered.}
#'   \item{n_eff}{Unclustered sample size \eqn{n_0} for the target
#'     precision, rounded up; for a full design, \code{n / deff}.}
#'   \item{moe}{Margin of error of the returned design: the half-width of
#'     its expected interval.}
#'
#' @references
#' Rogan WJ, Gladen B (1978). Estimating prevalence from the results of a
#' screening test. American Journal of Epidemiology 107(1):71-76.
#'
#' MMS-SD Study Design Workshop, Modules 1 (sampling), 2 (sample size from
#' margin of error) and 5 (ICC / design effect).
#' \url{https://mrc-ide.github.io/MMS-SD_workshop/}
#'
#' @export
#'
#' @examples
#' # Unclustered, perfect test
#' design_precision(prevalence = 0.3, moe = 0.05)
#'
#' # Sized for a Wald interval instead of the default Clopper-Pearson
#' design_precision(prevalence = 0.3, moe = 0.05, ci_method = "wald")
#'
#' # Imperfect test (sensitivity 90%, specificity 95%)
#' design_precision(0.3, 0.05, sensitivity = 0.9, specificity = 0.95)
#'
#' # Fixed site size: how many sites do I need?
#' design_precision(0.3, 0.05, n_per_site = 10, icc = 0.05)
#'
#' # Fixed number of sites: how many people per site?
#' design_precision(0.3, 0.05, n_sites = 50, icc = 0.05)
#'
#' # Full design: what margin of error do 6 sites of 40 give?
#' design_precision(0.2, n_sites = 6, n_per_site = 40, icc = 0.05)
#'
#' # Sites of different sizes
#' design_precision(0.2, n_per_site = c(250, 250, 60, 60, 40, 40), icc = 0.05)
design_precision <- function(prevalence,
                             moe         = NULL,
                             sensitivity = 1,
                             specificity = 1,
                             conf_level  = 0.95,
                             ci_method   = c("clopper-pearson", "wald",
                                             "agresti-coull"),
                             n           = NULL,
                             n_sites     = NULL,
                             n_per_site  = NULL,
                             icc         = 0) {

  .check_number(prevalence, "prevalence", 0, 1, lower_open = TRUE, upper_open = TRUE)
  .check_number(conf_level, "conf_level", 0, 1, lower_open = TRUE, upper_open = TRUE)
  .check_test(sensitivity, specificity)
  .check_number(icc, "icc", 0, 1)
  des <- .site_design(n, n_sites, n_per_site, icc)
  .check_target(des$full, !is.null(moe), "moe", has_default = FALSE)
  ci_method <- match.arg(ci_method)

  alpha <- 1 - conf_level
  z     <- stats::qnorm(1 - alpha / 2)
  p_app <- .apparent_prev(prevalence, sensitivity, specificity)
  corr  <- sensitivity + specificity - 1
  # Margin of error at effective sample size m: the half-width of the
  # interval estimate_prevalence() would give if the measured prevalence
  # came out at p_app. Wald is the usual closed form (not cut at 0 or 1).
  moe_at <- function(m) {
    if (ci_method == "wald") return(z * sqrt(p_app * (1 - p_app) / m) / corr)
    ci <- .correct_ci(.ci_app(p_app, m, alpha, ci_method), sensitivity, specificity)
    (ci[2] - ci[1]) / 2
  }
  moe_of <- function(n, deff) moe_at(n / deff)

  # ---- full design: its margin of error ----
  if (des$full) {
    deff <- .deff(des$n_w, icc)
    return(c(list(n = des$n), .report_sizes(des$sizes),
             list(deff = deff, n_eff = des$n / deff, moe = moe_of(des$n, deff))))
  }

  .check_number(moe, "moe", 0, 0.5, lower_open = TRUE, upper_open = TRUE)

  # unclustered sample size: a closed form for Wald; otherwise the margin
  # of error falls as m grows, so one root search on log(m) finds it
  n0 <- if (ci_method == "wald") z^2 * p_app * (1 - p_app) / (moe^2 * corr^2)
        else if (moe_at(1) <= moe) 1
        else exp(stats::uniroot(function(u) moe_at(exp(u)) - moe, c(0, 1),
                                extendInt = "downX", tol = 1e-10)$root)
  # clustered sample size
  # n0 is the number of independent people needed. Clustering makes people
  # at the same site partly redundant, so each design below scales n0 up by
  # the design effect, deff = 1 + (n_per_site - 1) * icc. n is the
  # (continuous) total sample size; it is rounded to whole sites below.
  if (!is.null(n_per_site)) {
    # Site size fixed: deff is known directly, so n = n0 * deff. The number
    # of sites, ceiling(n / n_per_site), is worked out by .allocate() below.
    n <- n0 * .deff(n_per_site, icc)
  } else if (!is.null(n_sites)) {
    # Number of sites fixed: deff depends on people per site, n / n_sites,
    # which depends on n. Solving n = n0 * deff(n / n_sites) gives a closed
    # form (see Details), which exists only if n_sites > n0 * icc.
    n <- .solve_fixed_sites(n0, n_sites, icc)
    if (is.na(n)) {
      # No solution: each site is worth at most 1 / icc independent people
      # however many are sampled there, so n_sites sites can never supply
      # n0 independent people. Report the best MOE they can reach, i.e.
      # the MOE from an effective sample size of n_sites / icc.
      min_moe <- moe_at(n_sites / icc)
      stop(sprintf(paste0("Target MOE of %.1f%% is unachievable with %d sites ",
                          "and ICC = %.3f; the minimum achievable MOE is %.1f%%."),
                   100 * moe, n_sites, icc, 100 * min_moe), call. = FALSE)
    }
  } else {
    # Unclustered: every person is independent, so n0 is the answer.
    n <- n0
  }

  design <- .allocate(n, n_sites, n_per_site)
  deff   <- if (is.null(design$n_per_site)) 1 else .deff(design$n_per_site, icc)

  list(
    n          = design$n,
    n_sites    = design$n_sites,
    n_per_site = design$n_per_site,
    deff       = deff,
    n_eff      = ceiling(n0),
    moe        = moe_of(design$n, deff)
  )
}
