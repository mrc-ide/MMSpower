#' Estimate prevalence and its precision
#'
#' @description
#' Analysis function: given observed counts, returns a prevalence point
#' estimate and confidence interval. Three CI methods are available:
#'
#' \describe{
#'   \item{`"clopper-pearson"` (default)}{Asymmetric. Exact binomial
#'     interval, via the beta distribution. For a simple random sample, its
#'     coverage is at least the stated confidence level, so it is usually
#'     wider than the Wald interval. It is also the interval that matches
#'     the exact test in \code{test_threshold()}.}
#'   \item{`"wald"`}{Symmetric unless an end is moved to 0 or 1.
#'     \eqn{\hat{p} \pm z \cdot SE}. Can be too narrow for small samples or
#'     prevalence near 0 or 1, and has zero width when no one (or everyone)
#'     tests positive.}
#'   \item{`"agresti-coull"`}{Asymmetric. A Wald interval calculated after
#'     adding \eqn{z^2/2} positive and \eqn{z^2/2} negative cases to the data
#'     (about 1.92 of each at 95 percent). This moves the interval towards
#'     0.5; the reported prevalence is unchanged.}
#' }
#'
#' All methods account for imperfect diagnostic tests (Rogan-Gladen
#' correction) and clustered sampling (Kish design effect).
#'
#' If `icc` is not supplied, it is estimated from the data (the observed
#' variance between sites compared with the variance expected if every
#' person were independent). Setting `icc = 0` treats every person as
#' independent; if people at the same site tend to be similar, this gives an
#' interval that is too narrow.
#'
#' @param x Vector of whole-number counts of people who tested positive, one
#'   per site (a single number for one site).
#' @param n Vector of whole-number totals tested, one per site (same length
#'   as `x`).
#' @param sensitivity Diagnostic sensitivity in (0, 1]; default 1 (the test
#'   never misses a case).
#' @param specificity Diagnostic specificity in (0, 1]; default 1 (the test
#'   never gives a false positive). `sensitivity + specificity` must exceed
#'   1.
#' @param conf_level Confidence level, in (0, 1); default 0.95.
#' @param icc Optional. Intra-cluster correlation, in \[0, 1\]. If `NULL`
#'   (default), ICC is estimated from the data. Set to `0` for no
#'   clustering. A supplied value is also used with a single site.
#' @param ci_method CI method: `"clopper-pearson"` (default), `"wald"`, or
#'   `"agresti-coull"`. See Description. \code{design_precision()} takes the
#'   same `ci_method`, with the same default, so a design sized there gives
#'   back its margin of error here.
#'
#' @details
#' The interval is built in five steps: apparent prevalence, a design
#' effect for clustering, an effective sample size, a confidence interval by
#' the chosen method, and finally the Rogan-Gladen correction for an
#' imperfect test.
#'
#' \strong{1. Apparent prevalence.} The pooled proportion of test
#' positives across all sites,
#'
#' \deqn{\hat{p} = \frac{\sum_i x_i}{\sum_i n_i}}
#'
#' It is corrected for test imperfection only at step 5, so the correction
#' applies identically to the point estimate and to both interval endpoints.
#'
#' \strong{2. Design effect (Kish).} Observations within a site are
#' positively correlated, so a clustered sample carries less information
#' than the same number of independent people. The Kish (1965) design
#' effect scales the variance by
#'
#' \deqn{D_{eff} = 1 + (n_w - 1)\,\rho}
#'
#' where \eqn{n_w = \sum n^2 / \sum n} is the mean number of people per site,
#' weighted by site size (the plain mean when all sites are the same size),
#' and \eqn{\rho} the intra-cluster correlation (ICC).
#'
#' \emph{ICC supplied}: \eqn{\rho} is used directly, also for a single site.
#'
#' \emph{ICC estimated} (\code{icc = NULL}, the default): with two or more
#' sites and an average of more than one person per site, the design effect
#' is first estimated as the ratio of the observed variance of the site
#' proportions to the variance expected under simple random sampling,
#'
#' \deqn{\widehat{D} = \frac{\mathrm{Var}(\hat{p}_i)}{\;\overline{\hat{p}(1 - \hat{p}) / n_i}\;}, \qquad \hat{p}_i = \frac{x_i}{n_i}}
#'
#' and the ICC recovered from \eqn{\widehat{D} = 1 + (n_h - 1)\,\hat{\rho}},
#'
#' \deqn{\hat{\rho} = \frac{\widehat{D} - 1}{n_h - 1}, \qquad
#'   n_h = \frac{k}{\sum_i 1 / n_i}}
#'
#' kept between 0 and 1, where \eqn{k} is the number of sites and
#' \eqn{n_h} the harmonic mean site size. The harmonic mean is the one that
#' links \eqn{\widehat{D}} to the ICC when sites differ in size: the
#' variance of a site's prevalence is \eqn{p(1-p)(1 + (n_i - 1)\rho) / n_i},
#' and averaging over sites gives \eqn{n_h}. \eqn{D_{eff}} is then computed
#' from \eqn{\hat{\rho}} as above, with the size-weighted mean \eqn{n_w},
#' because that is the design effect of the pooled prevalence
#' \eqn{\hat{p}}. With equal sites the two means are the same. With a
#' single site, one person per site, or nobody (or everybody) positive, the
#' ICC is 0. If everyone tested comes from one site, supply \code{icc}.
#'
#' \emph{Few sites.} An ICC estimated from \eqn{k} sites rests on only
#' \eqn{k - 1} degrees of freedom (returned as \code{df}), so the design
#' effect is itself uncertain. The interval allows for this as in Korn and
#' Graubard (1998): see step 4. In simulations with a true ICC of 0.1, 50
#' people per site and a prevalence of 0.3, the 95 percent Clopper-Pearson
#' interval contained the true prevalence 99 percent of the time with 3
#' sites, 96 with 5 and 95 to 96 with 10 or more (Wald: 97, 94 and 94). With
#' 2 sites (1 degree of freedom) the interval is very wide, which is an
#' honest picture of how little 2 sites say. Supplying \code{icc} from
#' other studies gives narrower intervals when it can be trusted.
#'
#' \strong{3. Effective sample size.}
#'
#' \deqn{n_{eff} = \frac{n_{total}}{D_{eff}}}
#'
#' the number of independent observations carrying the same information as
#' the clustered sample.
#'
#' \strong{4. Confidence interval on apparent prevalence.} Let
#' \eqn{z = \Phi^{-1}(1 - \alpha/2)}, \eqn{\alpha = 1 - } \code{conf_level}.
#' When the ICC is estimated, the effective sample size used for the
#' interval is shrunk to
#'
#' \deqn{n_{eff} \left(\frac{z}{t_{k-1}}\right)^2}
#'
#' where \eqn{t_{k-1}} is the matching quantile of a t distribution with
#' \eqn{k - 1} degrees of freedom (Korn and Graubard, 1998). For the Wald
#' interval this is the same as using \eqn{t_{k-1}} in place of \eqn{z}.
#' With \code{icc} supplied, or one site, nothing changes. Below,
#' \eqn{n_{eff}} stands for this adjusted size, and
#' \eqn{x_{eff} = \hat{p}\,n_{eff}} is the number of positives scaled to
#' it. Any end below 0 is set to 0 and any end above 1 is set to 1.
#'
#' \emph{\code{"wald"}} -- the normal-approximation interval,
#'
#' \deqn{\hat{p} \;\pm\; z \sqrt{\frac{\hat{p}(1 - \hat{p})}{n_{eff}}}}
#'
#' For small \eqn{n} or prevalence near 0 or 1 it can be too narrow, and it
#' has zero width when \eqn{\hat{p} = 0} or \eqn{1}. For example, with 30
#' people and a true prevalence of 5 percent, a 95 percent Wald interval
#' contains the true prevalence only about 78 percent of the time.
#'
#' \emph{\code{"clopper-pearson"}} -- the exact binomial interval,
#'
#' \deqn{L = B^{-1}\!\left(\tfrac{\alpha}{2};\; x_{eff},\; n_{eff} - x_{eff} + 1\right)}
#' \deqn{U = B^{-1}\!\left(1 - \tfrac{\alpha}{2};\; x_{eff} + 1,\; n_{eff} - x_{eff}\right)}
#'
#' where \eqn{B^{-1}(q; a, b)} is the \eqn{q}-th quantile of a beta
#' distribution, and \eqn{L = 0} when \eqn{\hat{p} = 0} and \eqn{U = 1} when
#' \eqn{\hat{p} = 1}. \eqn{x_{eff}} need not be a whole number, so the
#' design effect can be applied; coverage of at least the stated level then
#' holds only approximately.
#'
#' \emph{\code{"agresti-coull"}} -- add \eqn{z^2/2} positive and \eqn{z^2/2}
#' negative cases, then take a Wald interval on the adjusted proportion,
#'
#' \deqn{\tilde{n} = n_{eff} + z^2, \qquad \tilde{p} = \frac{x_{eff} + z^2/2}{\tilde{n}}}
#' \deqn{\tilde{p} \;\pm\; z \sqrt{\frac{\tilde{p}(1 - \tilde{p})}{\tilde{n}}}}
#'
#' \strong{5. Rogan-Gladen correction.} An imperfect test inflates apparent
#' prevalence through false positives and deflates it through false
#' negatives. Inverting \eqn{p_{app} = p\,Se + (1 - p)(1 - Sp)} gives true
#' prevalence,
#'
#' \deqn{p = \frac{p_{app} - (1 - Sp)}{Se + Sp - 1}}
#'
#' This same formula is applied to \eqn{\hat{p}} and to both CI endpoints,
#' and the results are kept within \[0, 1\]. If the whole corrected
#' interval lies below 0 (fewer positives than false positives alone would
#' produce) or above 1, the function warns that the data do not fit the
#' assumed sensitivity and specificity.
#'
#' \strong{Margin of error.} \code{moe} is the half-width of the corrected
#' interval, \code{(ci_upper - ci_lower) / 2}. Clopper-Pearson and
#' Agresti-Coull intervals are asymmetric, and so is a Wald interval with an
#' end moved to 0 or 1; report \code{ci_lower} and \code{ci_upper} then.
#'
#' @references
#' MMS-SD Study Design Workshop. \url{https://mrc-ide.github.io/MMS-SD_workshop/}
#'
#' Kish, L. (1965) \emph{Survey Sampling}. Wiley.
#'
#' Clopper, C. J. & Pearson, E. S. (1934) The use of confidence or fiducial
#' limits illustrated in the case of the binomial. \emph{Biometrika}
#' \strong{26}(4), 404-413. \doi{10.1093/biomet/26.4.404}
#'
#' Agresti, A. & Coull, B. A. (1998) Approximate is better than "exact" for
#' interval estimation of binomial proportions. \emph{The American
#' Statistician} \strong{52}(2), 119-126. \doi{10.1080/00031305.1998.10480550}
#'
#' Korn, E. L. & Graubard, B. I. (1998) Confidence intervals for
#' proportions with small expected number of positive counts estimated
#' from survey data. \emph{Survey Methodology} \strong{24}(2), 193-201.
#' (Degrees-of-freedom adjustment.)
#'
#' Rogan, W. J. & Gladen, B. (1978) Estimating prevalence from the results
#' of a screening test. \emph{American Journal of Epidemiology}
#' \strong{107}(1), 71-76. \doi{10.1093/oxfordjournals.aje.a112510}
#'
#' @return A list with:
#'   \item{prevalence}{Point estimate of true prevalence (Rogan-Gladen
#'     corrected).}
#'   \item{ci_lower, ci_upper}{Confidence limits on the true-prevalence
#'     scale.}
#'   \item{moe}{Half-width of the interval, \code{(ci_upper - ci_lower) / 2}.}
#'   \item{n_eff}{Effective sample size, \code{sum(n) / deff}.}
#'   \item{icc}{ICC used: as supplied, or estimated from the data.}
#'   \item{deff}{Design effect applied (1 when there is no clustering).}
#'   \item{df}{Degrees of freedom of the estimated design effect, sites - 1;
#'     \code{Inf} when \code{icc} is supplied or not estimated.}
#'
#' @export
#'
#' @examples
#' # Single site / simple random sample
#' estimate_prevalence(x = 8, n = 50)
#'
#' # Multi-site clustered data -- ICC estimated from the data
#' estimate_prevalence(
#'   x = c(0, 4, 0, 22, 25, 16, 12, 8),
#'   n = c(60, 80, 70, 100, 40, 60, 50, 90)
#' )
#'
#' # Wald interval instead of the default Clopper-Pearson
#' estimate_prevalence(x = 3, n = 30, ci_method = "wald")
#'
#' # Imperfect diagnostic test
#' estimate_prevalence(x = 30, n = 100, sensitivity = 0.9, specificity = 0.95)
#'
#' # Known ICC
#' estimate_prevalence(x = c(5, 8, 3), n = c(40, 40, 40), icc = 0.05)
estimate_prevalence <- function(x,
                                n,
                                sensitivity = 1,
                                specificity = 1,
                                conf_level  = 0.95,
                                icc         = NULL,
                                ci_method   = c("clopper-pearson", "wald",
                                                "agresti-coull")) {

  .check_counts(x, n)
  .check_test(sensitivity, specificity)
  .check_number(conf_level, "conf_level", 0, 1, lower_open = TRUE, upper_open = TRUE)
  if (!is.null(icc)) .check_number(icc, "icc", 0, 1)
  ci_method <- match.arg(ci_method)

  # apparent prevalence, design effect and effective sample size
  n_total <- sum(n)
  p_hat   <- sum(x) / n_total
  clus    <- .estimate_deff(x, n, icc)
  n_eff   <- n_total / clus$deff

  # With an estimated ICC the design effect is itself uncertain: shrink the
  # effective sample size for the interval (Korn-Graubard), which makes the
  # Wald interval use a t value with sites - 1 degrees of freedom.
  alpha <- 1 - conf_level
  n_ci  <- .df_adjust(n_eff, alpha / 2, clus$df)

  # interval on the apparent prevalence
  ci_app <- .ci_app(p_hat, n_ci, alpha, ci_method)

  # Rogan-Gladen correction to the true-prevalence scale, kept in [0, 1]
  ci_true    <- .correct_ci(ci_app, sensitivity, specificity)
  prevalence <- min(max(.rogan_gladen(p_hat, sensitivity, specificity), 0), 1)

  list(
    prevalence = prevalence,
    ci_lower   = ci_true[1],
    ci_upper   = ci_true[2],
    moe        = (ci_true[2] - ci_true[1]) / 2,
    n_eff      = n_eff,
    icc        = clus$icc,
    deff       = clus$deff,
    df         = clus$df
  )
}
