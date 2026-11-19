#' Test whether prevalence differs between two groups
#'
#' @description
#' Analysis function: given the number of positives and the number tested in
#' two groups, tests whether their true prevalence differs. Returns the
#' estimated difference with the confidence interval that matches the test,
#' the p-value and whether the difference is significant. It is the analysis
#' twin of \code{design_difference()}.
#'
#' The test compares the measured (apparent) prevalences, the share who
#' test positive. Because sensitivity and specificity are assumed the same
#' in both groups, the measured prevalences differ exactly when the true
#' ones do, in the same direction, so the p-value does not depend on
#' sensitivity or specificity. Only the reported difference and its
#' interval are corrected for an imperfect test (Rogan-Gladen correction).
#' Clustering is handled with a design effect for each group.
#'
#' \strong{Hypotheses} (\eqn{p_1}, \eqn{p_2} are the true prevalences):
#' \describe{
#'   \item{\code{"two.sided"} (default)}{H0: \eqn{p_1 = p_2} against
#'     H1: \eqn{p_1 \ne p_2}.}
#'   \item{\code{"greater"}}{H0: \eqn{p_1 \le p_2} against H1:
#'     \eqn{p_1 > p_2}.}
#'   \item{\code{"less"}}{H0: \eqn{p_1 \ge p_2} against H1: \eqn{p_1 < p_2}.}
#' }
#'
#' @param x1 Vector of whole-number counts of people who tested positive in
#'   group 1, one per site (a single number if group 1 is one site).
#' @param n1 Vector of whole-number totals tested in group 1, one per site
#'   (same length as \code{x1}).
#' @param x2 Vector of whole-number counts of people who tested positive in
#'   group 2, one per site (a single number if group 2 is one site).
#' @param n2 Vector of whole-number totals tested in group 2, one per site
#'   (same length as \code{x2}).
#' @param alternative Direction of the test: \code{"two.sided"} (default),
#'   \code{"greater"} (group 1 higher) or \code{"less"} (group 1 lower).
#' @param sensitivity Diagnostic sensitivity in (0, 1]; default 1. Assumed
#'   the same in both groups.
#' @param specificity Diagnostic specificity in (0, 1]; default 1. Assumed
#'   the same in both groups. \code{sensitivity + specificity} must exceed
#'   1.
#' @param conf_level Confidence level for the interval on the difference, in
#'   (0, 1); default 0.95. The significance level of the test is
#'   1 - \code{conf_level}.
#' @param icc Optional. Intra-cluster correlation (ICC), in \[0, 1\]. If
#'   \code{NULL} (default), it is estimated separately for each group from
#'   its own sites. Otherwise give one number used for both groups, or two,
#'   \code{c(group 1, group 2)}. A supplied value is also used for a group
#'   with only one site (see Details). Set to \code{0} for no clustering.
#' @param method \code{"score"} (default): pooled z-test with the score
#'   interval. \code{"wald"}: unpooled (Wald) z-test with the Wald interval.
#'   Each interval is its own test turned around, so the interval leaves out
#'   0 exactly when \code{reject} is \code{TRUE}, whichever method is used.
#'
#' @details
#' \strong{Design effect} (returned as \code{deff1} and \code{deff2}). For
#' each group it is \eqn{1 + (n_w - 1) ICC}, where
#' \eqn{n_w = \sum n^2 / \sum n} is the mean number of people per site,
#' weighted by site size. With \code{icc} supplied (for example from
#' earlier surveys), that value is used. Supply \code{icc} for a group with
#' one site only when the site stands for a wider population, such as a
#' district; if the site itself is the population of interest, leave
#' \code{icc = NULL}. With \code{icc = NULL}, the ICC of each group is
#' estimated from its own data as in \code{estimate_prevalence()}; this
#' needs at least two sites, and a group with one site gets a design effect
#' of 1.
#'
#' An ICC estimated from only a few sites is uncertain, and the interval
#' then tends to be too narrow. In simulations with a true ICC of 0.1, 50
#' people per site and prevalences of 0.2 to 0.5, a 95 percent interval
#' contained the true difference about 83 to 85 percent of the time with 2
#' sites per group, 87 to 89 percent with 3 and 90 to 92 percent with 5.
#' Supplying \code{icc} from other studies avoids this.
#'
#' \strong{Effective sample size.} For group \eqn{g}, \eqn{m_g} (returned as
#' \code{n_eff1} and \code{n_eff2}) is the number tested divided by the
#' design effect.
#'
#' \strong{Score test} (\code{method = "score"}). With \eqn{\hat{p}_g} the
#' measured prevalence in group \eqn{g},
#'
#' \deqn{z = \frac{\hat{p}_1 - \hat{p}_2}
#'   {\sqrt{\bar{p}(1-\bar{p})\,(1/m_1 + 1/m_2)}}, \qquad
#'   \bar{p} = \frac{\hat{p}_1 m_1 + \hat{p}_2 m_2}{m_1 + m_2}.}
#'
#' With no clustering this is the same test as R's
#' \code{prop.test(correct = FALSE)}, and the one \code{design_difference()}
#' sizes for.
#'
#' \strong{Wald test} (\code{method = "wald"}). The same, but with the
#' standard error worked out from each group's own prevalence,
#'
#' \deqn{z = \frac{\hat{p}_1 - \hat{p}_2}
#'   {\sqrt{\hat{p}_1(1-\hat{p}_1)/m_1 + \hat{p}_2(1-\hat{p}_2)/m_2}}.}
#'
#' With small groups it can be too optimistic: with 20 people per group and
#' true prevalences of 0.05 and 0.10, a 95 percent Wald interval contains
#' the true difference only about 90 percent of the time (score 97). When
#' each group's measured prevalence is 0 or 1 the standard error is 0.
#'
#' For both, the one-sided p-values are the normal tail areas above and
#' below \eqn{z}, and the two-sided p-value is twice the smaller of them.
#'
#' \strong{Interval.} The interval is the test turned around: it holds the
#' differences \eqn{\delta} the test would not reject, so it leaves out 0
#' exactly when \code{reject} is \code{TRUE}. It matches
#' \code{alternative}: two-sided for \code{"two.sided"}, \eqn{[L, 1]} for
#' \code{"greater"} and \eqn{[-1, U]} for \code{"less"}.
#'
#' For \code{"score"}, the two prevalences are estimated again by maximum
#' likelihood under the condition \eqn{p_1 - p_2 = \delta} (the closed form
#' of Farrington and Manning, 1990), giving \eqn{\tilde{p}_1, \tilde{p}_2},
#' and
#'
#' \deqn{z(\delta) = \frac{\hat{p}_1 - \hat{p}_2 - \delta}
#'   {\sqrt{\tilde{p}_1(1-\tilde{p}_1)/m_1 + \tilde{p}_2(1-\tilde{p}_2)/m_2}}}
#'
#' is set equal to the z value for \code{conf_level} and solved for
#' \eqn{\delta} (Mee, 1984). At \eqn{\delta = 0}, \eqn{\tilde{p}_1 =
#' \tilde{p}_2 = \bar{p}} and this is exactly the score test above. Miettinen
#' and Nurminen (1985) multiply the variance by
#' \eqn{(m_1 + m_2)/(m_1 + m_2 - 1)}; that factor is left out here so that
#' the interval and the test always agree.
#'
#' For \code{"wald"}, the interval is \eqn{\hat{p}_1 - \hat{p}_2 \pm z
#' \sqrt{\hat{p}_1(1 - \hat{p}_1)/m_1 + \hat{p}_2(1 - \hat{p}_2)/m_2}}.
#'
#' The interval is built on the measured scale, then divided by
#' \eqn{Se + Sp - 1} to give an interval for the true difference (the rest
#' of the Rogan-Gladen correction cancels in a difference). Ends outside -1
#' to 1 are set to -1 or 1.
#'
#' @return A list with:
#'   \item{difference}{Estimated true difference, group 1 minus group 2:
#'     the measured difference divided by \eqn{Se + Sp - 1}, between -1 and
#'     1.}
#'   \item{ci_lower, ci_upper}{Confidence interval for the true difference;
#'     \code{ci_upper} is 1 for \code{"greater"} and \code{ci_lower} is -1
#'     for \code{"less"}.}
#'   \item{p_value}{p-value for the chosen \code{alternative}.}
#'   \item{reject}{\code{TRUE} if \code{p_value} is below
#'     1 - \code{conf_level}.}
#'   \item{prevalence1, prevalence2}{Estimated true prevalence in each group
#'     (Rogan-Gladen corrected, between 0 and 1).}
#'   \item{n_eff1, n_eff2}{Effective sample size of each group,
#'     \code{sum(n) / deff}.}
#'   \item{icc1, icc2}{ICC used for each group: as supplied, or estimated
#'     from that group's data.}
#'   \item{deff1, deff2}{Design effect used for each group (1 with no
#'     clustering).}
#'
#' @references
#' MMS-SD Study Design Workshop. \url{https://mrc-ide.github.io/MMS-SD_workshop/}
#'
#' Fleiss, J. L., Levin, B. & Paik, M. C. (2003) \emph{Statistical Methods for
#' Rates and Proportions}, 3rd ed. Wiley. (Two-group z-test.)
#'
#' Mee, R. W. (1984) Confidence bounds for the difference between two
#' probabilities. \emph{Biometrics} \strong{40}(4), 1175-1176.
#'
#' Miettinen, O. & Nurminen, M. (1985) Comparative analysis of two rates.
#' \emph{Statistics in Medicine} \strong{4}(2), 213-226.
#'
#' Farrington, C. P. & Manning, G. (1990) Test statistics and sample size
#' formulae for comparative binomial trials with null hypothesis of
#' non-zero risk difference or non-unity relative risk. \emph{Statistics in
#' Medicine} \strong{9}(12), 1447-1454. (Prevalences estimated under a
#' given difference.)
#'
#' Kish, L. (1965) \emph{Survey Sampling}. Wiley. (Design effect.)
#'
#' Rogan, W. J. & Gladen, B. (1978) Estimating prevalence from the results
#' of a screening test. \emph{American Journal of Epidemiology}
#' \strong{107}(1), 71-76. \doi{10.1093/oxfordjournals.aje.a112510}
#'
#' @seealso \code{\link{design_difference}}
#'
#' @export
#'
#' @examples
#' # One count per group: 8 of 100 versus 20 of 120
#' test_difference(x1 = 8, n1 = 100, x2 = 20, n2 = 120)
#'
#' # Three sites per group, one-sided, imperfect diagnostic test. The sites
#' # differ, so the estimated design effects are above 1 (about 2.6 and 3.3)
#' # and p is 0.07; treating everyone as independent would give 0.007.
#' test_difference(
#'   x1 = c(1, 9, 3), n1 = c(50, 60, 40),
#'   x2 = c(4, 18, 9), n2 = c(55, 70, 45),
#'   alternative = "less",
#'   sensitivity = 0.9, specificity = 0.97
#' )
#'
#' # Wald test and interval instead of the default score versions
#' test_difference(x1 = 3, n1 = 80, x2 = 12, n2 = 90, method = "wald")
#'
#' # A different ICC for each group
#' test_difference(
#'   x1 = c(2, 5, 3), n1 = c(50, 60, 40),
#'   x2 = c(8, 11, 9), n2 = c(55, 70, 45),
#'   icc = c(0.02, 0.08)
#' )
test_difference <- function(x1, n1, x2, n2,
                            alternative = c("two.sided", "greater", "less"),
                            sensitivity = 1,
                            specificity = 1,
                            conf_level  = 0.95,
                            icc         = NULL,
                            method      = c("score", "wald")) {

  .check_counts(x1, n1, "x1", "n1")
  .check_counts(x2, n2, "x2", "n2")
  .check_test(sensitivity, specificity)
  .check_number(conf_level, "conf_level", 0, 1, lower_open = TRUE, upper_open = TRUE)
  if (!is.null(icc)) {
    if (!length(icc) %in% 1:2)
      stop("`icc` must be NULL, one number, or two numbers.", call. = FALSE)
    for (v in icc) .check_number(v, "icc", 0, 1)
  }
  alternative <- match.arg(alternative)
  method      <- match.arg(method)

  # measured prevalence, design effect and effective size for each group
  # (one supplied icc is used for both)
  p1 <- sum(x1) / sum(n1)
  p2 <- sum(x2) / sum(n2)
  g1 <- .estimate_deff(x1, n1, icc[1])
  g2 <- .estimate_deff(x2, n2, icc[length(icc)])
  m1 <- sum(n1) / g1$deff
  m2 <- sum(n2) / g2$deff
  d  <- p1 - p2

  # z-statistic. The score version uses the pooled prevalence (the
  # restricted estimate at delta = 0); the Wald version each group's own.
  if (method == "score") {
    z <- .score_stat(p1, m1, p2, m2, 0)
  } else {
    se <- sqrt(p1 * (1 - p1) / m1 + p2 * (1 - p2) / m2)
    z  <- if (d == 0) 0 else d / se     # +/-Inf when se is 0
  }
  p_value <- switch(alternative,
    greater   = stats::pnorm(z, lower.tail = FALSE),
    less      = stats::pnorm(z),
    two.sided = 2 * stats::pnorm(-abs(z)))
  alpha <- 1 - conf_level

  # The interval is the test turned around: alpha in the one closed tail
  # for a one-sided test, alpha / 2 in each for two-sided; an open tail
  # runs to -1 or 1.
  alpha_lo <- switch(alternative, greater = alpha, less = 0, two.sided = alpha / 2)
  alpha_hi <- switch(alternative, greater = 0, less = alpha, two.sided = alpha / 2)
  z_lo <- if (alpha_lo > 0) stats::qnorm(alpha_lo, lower.tail = FALSE)
  z_hi <- if (alpha_hi > 0) stats::qnorm(alpha_hi, lower.tail = FALSE)
  if (method == "score") {
    ci_app <- c(if (is.null(z_lo)) -1 else .score_bound(p1, m1, p2, m2, z_lo, "lower"),
                if (is.null(z_hi))  1 else .score_bound(p1, m1, p2, m2, z_hi, "upper"))
  } else {
    ci_app <- c(if (is.null(z_lo)) -1 else d - z_lo * se,
                if (is.null(z_hi))  1 else d + z_hi * se)
  }

  # to the true scale: divide by Se + Sp - 1, keep within [-1, 1]
  corr  <- sensitivity + specificity - 1
  ci    <- pmin(pmax(ci_app / corr, -1), 1)
  rg    <- function(p) min(max(.rogan_gladen(p, sensitivity, specificity), 0), 1)

  list(
    difference  = min(max(d / corr, -1), 1),
    ci_lower    = ci[1],
    ci_upper    = ci[2],
    p_value     = p_value,
    reject      = p_value < alpha,
    prevalence1 = rg(p1),
    prevalence2 = rg(p2),
    n_eff1      = m1,
    n_eff2      = m2,
    icc1        = g1$icc,
    icc2        = g2$icc,
    deff1       = g1$deff,
    deff2       = g2$deff
  )
}


# ---- score interval helpers ----

# Prevalences estimated again by maximum likelihood under the condition
# t1 - t2 = delta, for measured prevalences p1, p2 with (possibly
# fractional) effective sizes m1, m2: the closed-form root of the
# Farrington-Manning (1990) cubic. Returns c(t1, t2).
.restricted_prev <- function(p1, m1, p2, m2, delta) {
  th <- m2 / m1
  a  <- 1 + th
  b  <- -(1 + th + p1 + th * p2 + delta * (th + 2))
  cc <- delta^2 + delta * (2 * p1 + th + 1) + p1 + th * p2
  d  <- -p1 * delta * (1 + delta)
  v  <- b^3 / (3 * a)^3 - b * cc / (6 * a^2) + d / (2 * a)
  u  <- sqrt(max(b^2 / (3 * a)^2 - cc / (3 * a), 0))
  # (max/min guard against rounding pushing v / u^3 a hair outside [-1, 1])
  t1 <- if (u == 0) -b / (3 * a) else
    2 * u * cos((pi + acos(max(-1, min(1, v / u^3)))) / 3) - b / (3 * a)
  # keep both prevalences within [0, 1]
  t1 <- min(min(1, 1 + delta), max(max(0, delta), t1))
  c(t1, t1 - delta)
}

# Score statistic for the difference delta (Mee, 1984). At delta = 0 it is
# the pooled z-statistic.
.score_stat <- function(p1, m1, p2, m2, delta) {
  num <- (p1 - p2) - delta
  if (num == 0) return(0)
  t <- .restricted_prev(p1, m1, p2, m2, delta)
  v <- t[1] * (1 - t[1]) / m1 + t[2] * (1 - t[2]) / m2
  if (v <= 0) return(sign(num) * Inf)
  num / sqrt(v)
}

# One end of the score interval: the delta where the statistic equals z
# (lower end) or -z (upper end), or -1 or 1 if it never gets there. The
# statistic falls as delta grows, so each end has a single root.
.score_bound <- function(p1, m1, p2, m2, z, side) {
  d <- p1 - p2
  e <- 1e-10
  f <- function(delta) .score_stat(p1, m1, p2, m2, delta)
  if (side == "lower") {
    if (d <= -1 + e || f(-1 + e) <= z) return(-1)
    stats::uniroot(function(x) f(x) - z, c(-1 + e, d), tol = 1e-12)$root
  } else {
    if (d >= 1 - e || f(1 - e) >= -z) return(1)
    stats::uniroot(function(x) f(x) + z, c(d, 1 - e), tol = 1e-12)$root
  }
}
