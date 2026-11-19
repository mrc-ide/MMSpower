#' Test whether prevalence is above or below a decision threshold
#'
#' @description
#' Analysis function: given the number of positives and the number tested
#' at one or more sites, tests whether the true prevalence is above (or
#' below) a decision threshold, with an exact binomial test (default) or a
#' Wald z-test. Returns the p-value, whether the result is significant, and
#' the confidence interval for the true prevalence that matches the test:
#' Clopper-Pearson for the exact test, Wald for the Wald test. It is the
#' analysis twin of \code{design_threshold()}.
#'
#' With a perfect diagnostic test, the share of people who test positive is
#' the prevalence, and it is compared with the threshold directly. With an
#' imperfect test it is not: some people with the condition test negative,
#' and some without it test positive. So the threshold is converted to the
#' share who would test positive if the true prevalence were exactly at the
#' threshold, \eqn{\theta \cdot Se + (1 - \theta)(1 - Sp)}, where
#' \eqn{\theta} is the threshold, \eqn{Se} the sensitivity and \eqn{Sp} the
#' specificity. For example, with a threshold of 5 percent, sensitivity 0.9
#' and specificity 0.95: of the 5 percent who have the condition, 0.9 test
#' positive (0.045), and of the 95 percent who do not, 0.05 wrongly test
#' positive (0.0475), so 9.25 percent would test positive. The test then
#' asks whether the share testing positive is far enough above 9.25 percent
#' that chance alone is unlikely to explain it: with 100 people tested, 12
#' positives (12 percent) is not enough (p = 0.21), but 15 (15 percent) is
#' (p = 0.04). This works because, when \eqn{Se + Sp} is above 1, a higher
#' true prevalence always gives a higher share testing positive. The
#' prevalence and its interval are then converted back to the true scale
#' (Rogan-Gladen correction). Clustering is allowed for with a design
#' effect, as in \code{estimate_prevalence()}.
#'
#' \strong{Hypotheses} (\eqn{p} is the true prevalence, \eqn{\theta} the
#' threshold):
#' \describe{
#'   \item{\code{"greater"} (default)}{H0: \eqn{p \le \theta} against
#'     H1: \eqn{p > \theta}: is prevalence above the threshold?}
#'   \item{\code{"less"}}{H0: \eqn{p \ge \theta} against H1:
#'     \eqn{p < \theta}: is prevalence below the threshold?}
#'   \item{\code{"two.sided"}}{H0: \eqn{p = \theta} against H1:
#'     \eqn{p \ne \theta}. A threshold question is usually one-way
#'     (above or below), so \code{"greater"} or \code{"less"} is usually
#'     the right choice.}
#' }
#'
#' @param x Vector of whole-number counts of people who tested positive, one
#'   per site (a single number for one site).
#' @param n Vector of whole-number totals tested, one per site (same length
#'   as \code{x}).
#' @param threshold Decision threshold for the true prevalence, strictly
#'   between 0 and 1 (for example 0.05 for the WHO pfhrp2/3 deletion
#'   threshold).
#' @param alternative Direction of the test: \code{"greater"} (default,
#'   prevalence above the threshold), \code{"less"} (below it) or
#'   \code{"two.sided"} (different from it).
#' @param sensitivity Diagnostic sensitivity in (0, 1]; default 1 (the test
#'   never misses a case).
#' @param specificity Diagnostic specificity in (0, 1]; default 1 (the test
#'   never gives a false positive). \code{sensitivity + specificity} must
#'   exceed 1.
#' @param conf_level Confidence level for the interval, in (0, 1); default
#'   0.95. The significance level of the test is 1 - \code{conf_level}.
#' @param icc Optional. Intra-cluster correlation (ICC), in \[0, 1\]. If
#'   \code{NULL} (default), it is estimated from the sites. Set to \code{0}
#'   for no clustering. A supplied value is also used with only one site
#'   (see Details).
#' @param method \code{"exact"} (default): exact binomial test with the
#'   Clopper-Pearson interval. \code{"wald"}: Wald z-test with the Wald
#'   interval. Each interval is its own test turned around, so
#'   \code{reject} is \code{TRUE} exactly when \code{threshold} is outside
#'   the interval, whichever method is used.
#'
#' @details
#' \strong{Design effect} (returned as \code{deff}). It is
#' \eqn{1 + (n_w - 1) ICC}, where \eqn{n_w = \sum n^2 / \sum n} is the
#' mean number of people per site, weighted by site size. With \code{icc}
#' supplied (for example from earlier surveys), that value is used. Supply
#' \code{icc} with one site only when the site stands for a wider
#' population, such as a district; if the site itself is the population of
#' interest, leave \code{icc = NULL}. With \code{icc = NULL}, the ICC is
#' estimated from the data as in \code{estimate_prevalence()}; this needs
#' at least two sites, and with one site the design effect is 1.
#'
#' \strong{Effective sample size.} \eqn{m} (returned as \code{n_eff}) is the
#' number tested divided by the design effect. The number of positives is
#' scaled the same way, to \eqn{\hat{p} m}, where
#' \eqn{\hat{p} = \sum x / \sum n} is the measured prevalence. Neither need
#' be a whole number: the binomial tail probabilities are worked out
#' through the beta distribution, which equals \code{pbinom()} for whole
#' numbers.
#'
#' \strong{Exact test} (\code{method = "exact"}). Suppose the true
#' prevalence is exactly at the threshold. For \code{"greater"}, the p-value
#' is the chance of at least as many positives as observed; for
#' \code{"less"}, the chance of at most as many. The threshold is first
#' converted to the measured scale, \eqn{\theta_{app} = \theta \cdot Se +
#' (1 - \theta)(1 - Sp)}. With one site and a design effect of 1, the
#' one-sided p-values are exactly those of R's \code{binom.test()}.
#'
#' \strong{Wald test} (\code{method = "wald"}). The normal approximation,
#'
#' \deqn{z = \frac{\hat{p} - \theta_{app}}{\sqrt{\hat{p}(1 - \hat{p}) / m}}}
#'
#' with the standard error worked out at the measured prevalence
#' \eqn{\hat{p}}. The one-sided p-values are the normal tail areas above
#' and below \eqn{z}. It can be too optimistic for small samples or a
#' prevalence near 0; with no positives (or all positive) the standard
#' error is 0, and any threshold on the other side counts as rejected.
#'
#' \strong{Two-sided p-value.} For both methods, twice the smaller of the
#' two one-sided p-values, at most 1. For the exact test this is the
#' "central" two-sided exact test (Fay 2010), the one that matches the
#' two-sided Clopper-Pearson interval. It is cautious: the false alarm rate
#' is at most 1 - \code{conf_level}, and often well below it. It can differ
#' from \code{binom.test()}, which adds up all counts that are less likely
#' than the one observed.
#'
#' \strong{Interval.} The interval is the test turned around: it holds the
#' thresholds the test would not reject, so \code{reject} is \code{TRUE}
#' exactly when \code{threshold} is outside it. It is matched to
#' \code{alternative}: \eqn{[L, 1]} for \code{"greater"}, \eqn{[0, U]} for
#' \code{"less"}, and the two-sided interval (\eqn{\alpha/2} in each tail)
#' for \code{"two.sided"}. For the exact test it is the Clopper-Pearson
#' interval,
#'
#' \deqn{L = B^{-1}(\alpha;\; x_{eff},\; m - x_{eff} + 1), \qquad
#'       U = B^{-1}(1 - \alpha;\; x_{eff} + 1,\; m - x_{eff})}
#'
#' with \eqn{x_{eff} = \hat{p} m} and \eqn{B^{-1}} the beta quantile
#' function; for the Wald test it is \eqn{\hat{p} \pm z_{1-\alpha}
#' \sqrt{\hat{p}(1 - \hat{p}) / m}} (\eqn{\alpha/2} in place of \eqn{\alpha}
#' for two-sided). It is built on the measured scale and then corrected for
#' an imperfect test, with ends kept within \[0, 1\]; if the whole
#' corrected interval lies outside \[0, 1\], the function warns that the
#' data do not fit the assumed sensitivity and specificity.
#'
#' \strong{When to use DRpower instead.} With the ICC supplied, this
#' function keeps the chance of a false alarm at or below the significance
#' level. With \code{icc = NULL} and few sites of very different sizes, the
#' ICC is poorly estimated and the false alarm rate can be well above it
#' (13 to 28 percent at a 5 percent level in simulation); with an ICC
#' estimated from few sites it also has low power. For pfhrp2/3 deletion
#' surveys with several sites and an unknown ICC, where the question is
#' whether prevalence is above the threshold, consider DRpower
#' (\url{https://mrc-ide.github.io/DRpower/}). It estimates prevalence and
#' ICC together and has more power, especially with few sites, but more
#' false alarms when prevalence is close to the threshold.
#'
#' @return A list with:
#'   \item{p_value}{p-value for the chosen \code{alternative}.}
#'   \item{reject}{\code{TRUE} if \code{p_value} is below
#'     1 - \code{conf_level}.}
#'   \item{prevalence}{Estimated true prevalence (Rogan-Gladen corrected,
#'     between 0 and 1).}
#'   \item{ci_lower, ci_upper}{Confidence interval for the true prevalence;
#'     \code{ci_upper} is 1 for \code{"greater"} and \code{ci_lower} is 0
#'     for \code{"less"}.}
#'   \item{n_eff}{Effective sample size, \code{sum(n) / deff}.}
#'   \item{icc}{ICC used: as supplied, or estimated from the data.}
#'   \item{deff}{Design effect used (1 with no clustering).}
#'
#' @references
#' MMS-SD Study Design Workshop. \url{https://mrc-ide.github.io/MMS-SD_workshop/}
#'
#' Clopper, C. J. & Pearson, E. S. (1934) The use of confidence or fiducial
#' limits illustrated in the case of the binomial. \emph{Biometrika}
#' \strong{26}(4), 404-413. \doi{10.1093/biomet/26.4.404}
#'
#' Fay, M. P. (2010) Two-sided exact tests and matching confidence intervals
#' for discrete data. \emph{The R Journal} \strong{2}(1), 53-58.
#' \doi{10.32614/RJ-2010-008}
#'
#' Kish, L. (1965) \emph{Survey Sampling}. Wiley. (Design effect.)
#'
#' Rogan, W. J. & Gladen, B. (1978) Estimating prevalence from the results
#' of a screening test. \emph{American Journal of Epidemiology}
#' \strong{107}(1), 71-76. \doi{10.1093/oxfordjournals.aje.a112510}
#'
#' @seealso \code{\link{design_threshold}}
#'
#' @export
#'
#' @examples
#' # One site: is prevalence above 5 percent?
#' test_threshold(x = 8, n = 100, threshold = 0.05)
#'
#' # Five sites, ICC estimated from the data
#' test_threshold(
#'   x = c(0, 6, 1, 11, 4),
#'   n = c(50, 60, 40, 80, 70),
#'   threshold = 0.05
#' )
#'
#' # Two-sided, imperfect diagnostic test
#' test_threshold(x = 12, n = 100, threshold = 0.10,
#'                alternative = "two.sided",
#'                sensitivity = 0.9, specificity = 0.95)
#'
#' # One site, a clear result: the interval lies above the threshold
#' test_threshold(x = 15, n = 100, threshold = 0.05)
#'
#' # Wald test and interval instead of the exact test
#' test_threshold(x = 15, n = 100, threshold = 0.05, method = "wald")
test_threshold <- function(x,
                           n,
                           threshold,
                           alternative = c("greater", "less", "two.sided"),
                           sensitivity = 1,
                           specificity = 1,
                           conf_level  = 0.95,
                           icc         = NULL,
                           method      = c("exact", "wald")) {

  .check_counts(x, n)
  .check_number(threshold, "threshold", 0, 1, lower_open = TRUE, upper_open = TRUE)
  .check_test(sensitivity, specificity)
  .check_number(conf_level, "conf_level", 0, 1, lower_open = TRUE, upper_open = TRUE)
  if (!is.null(icc)) .check_number(icc, "icc", 0, 1)
  alternative <- match.arg(alternative)
  method      <- match.arg(method)

  # measured prevalence, design effect and effective sample size
  p_hat <- sum(x) / sum(n)
  clus  <- .estimate_deff(x, n, icc)
  m     <- sum(n) / clus$deff
  x_eff <- p_hat * m     # positives scaled to m (not always whole)

  # One-sided p-values against the threshold on the measured scale: the
  # chance of a result at least this high (p_up) or this low (p_dn) if the
  # true prevalence were exactly at the threshold.
  theta_app <- .apparent_prev(threshold, sensitivity, specificity)
  if (method == "exact") {
    p_up <- .binom_upper(x_eff, m, theta_app)   # chance of x_eff or more
    p_dn <- .binom_lower(x_eff, m, theta_app)   # chance of x_eff or fewer
  } else {
    # Wald z-test, standard error at the measured prevalence. If p_hat is
    # 0 or 1 the standard error is 0 and z is -Inf or Inf.
    se   <- sqrt(p_hat * (1 - p_hat) / m)
    z    <- (p_hat - theta_app) / se
    p_up <- stats::pnorm(z, lower.tail = FALSE)
    p_dn <- stats::pnorm(z)
  }
  p_value <- switch(alternative,
    greater   = p_up,
    less      = p_dn,
    two.sided = min(1, 2 * min(p_up, p_dn)))
  alpha <- 1 - conf_level

  # The interval is the matching test turned around: the thresholds it
  # would not reject. alpha goes in the one closed tail for a one-sided
  # test, alpha / 2 in each for two-sided; an open tail runs to 0 or 1.
  alpha_lo <- switch(alternative, greater = alpha, less = 0, two.sided = alpha / 2)
  alpha_hi <- switch(alternative, greater = 0, less = alpha, two.sided = alpha / 2)
  if (method == "exact") {
    # Clopper-Pearson
    ci_app <- c(if (alpha_lo == 0 || p_hat == 0) 0 else
                  stats::qbeta(alpha_lo, x_eff, m - x_eff + 1),
                if (alpha_hi == 0 || p_hat == 1) 1 else
                  stats::qbeta(1 - alpha_hi, x_eff + 1, m - x_eff))
  } else {
    # Wald
    ci_app <- c(if (alpha_lo == 0) 0 else
                  p_hat - stats::qnorm(alpha_lo, lower.tail = FALSE) * se,
                if (alpha_hi == 0) 1 else
                  p_hat + stats::qnorm(alpha_hi, lower.tail = FALSE) * se)
    ci_app <- pmin(pmax(ci_app, 0), 1)
  }
  ci_true <- .correct_ci(ci_app, sensitivity, specificity)

  list(
    p_value    = p_value,
    reject     = p_value < alpha,
    prevalence = min(max(.rogan_gladen(p_hat, sensitivity, specificity), 0), 1),
    ci_lower   = ci_true[1],
    ci_upper   = ci_true[2],
    n_eff      = m,
    icc        = clus$icc,
    deff       = clus$deff
  )
}
