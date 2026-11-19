#' Calculate sample size to detect a difference in prevalence between two groups
#'
#' @description
#' Design function: given the expected true prevalence in two groups, returns
#' the sample size needed in each group for \code{test_difference()} (with
#' its default score test) to detect that difference with the target power.
#'
#' Each group can have its own number of sites, people per site and ICC:
#' give one value for both groups, or two, \code{c(group 1, group 2)}. The
#' function solves for one unknown, the same in both groups:
#' \describe{
#'   \item{No clustering (\code{n_sites} and \code{n_per_site} both
#'     \code{NULL}, default)}{Solves for the number of people per group.}
#'   \item{Fixed number of people per site (\code{n_per_site})}{Solves for
#'     the number of sites per group.}
#'   \item{Fixed number of sites (\code{n_sites})}{Solves for the number of
#'     people per site.}
#' }
#' So with one value of \code{n_sites} or \code{n_per_site} the two groups
#' are the same size, and with two values their sizes follow those values.
#' For example, \code{n_sites = c(10, 20)} gives group 2 twice as many people
#' as group 1 (the same number of people per site in both groups), and so
#' does \code{n_per_site = c(10, 20)} (the same number of sites in both
#' groups).
#'
#' @param prevalence1 Expected true prevalence in group 1, in \[0, 1\].
#' @param prevalence2 Expected true prevalence in group 2, in \[0, 1\].
#'   The two prevalences must differ.
#' @param power Target power, in \[0.5, 1); default 0.80.
#' @param alternative Which kind of difference to look for.
#'   \code{"two.sided"} (default) looks for any difference; \code{"greater"}
#'   looks for group 1 being higher (requires
#'   \code{prevalence1 > prevalence2}); \code{"less"} looks for group 1 being
#'   lower (requires \code{prevalence1 < prevalence2}).
#' @param sensitivity Diagnostic sensitivity in (0, 1]; default 1. Assumed
#'   the same in both groups.
#' @param specificity Diagnostic specificity in (0, 1]; default 1. Assumed
#'   the same in both groups. \code{sensitivity + specificity} must exceed
#'   1.
#' @param conf_level Confidence level, in (0, 1); default 0.95. The
#'   significance level of the test is 1 - \code{conf_level}.
#' @param n_sites Optional. Number of sites in each group: one whole number
#'   for both groups, or two, \code{c(group 1, group 2)}. The function solves
#'   for the number of people per site, the same in both groups. Cannot be
#'   used with \code{n_per_site}.
#' @param n_per_site Optional. Number of people per site: one whole number
#'   for both groups, or two, \code{c(group 1, group 2)}. The function solves
#'   for the number of sites, the same in both groups. Cannot be used with
#'   \code{n_sites}.
#' @param icc Intra-cluster correlation, in \[0, 1\]: one number for both
#'   groups, or two, \code{c(group 1, group 2)}. Default 0 (no clustering).
#'   An \code{icc} above 0 needs \code{n_sites} or \code{n_per_site}.
#'
#' @details
#' \strong{Sample size.} Write \eqn{m_1}, \eqn{m_2} for the effective sizes
#' of the two groups (below). The score test in \code{test_difference()}
#' reaches the target power when
#'
#' \deqn{z_\alpha \sqrt{\bar{p}(1-\bar{p})\,(1/m_1 + 1/m_2)} +
#'   z_\beta \sqrt{p_{1,app}(1-p_{1,app})/m_1 + p_{2,app}(1-p_{2,app})/m_2}
#'   = |p_{1,app} - p_{2,app}|}
#'
#' where \eqn{p_{g,app} = p_g \cdot Se + (1-p_g)(1-Sp)} is the prevalence the
#' test would measure in group \eqn{g}, \eqn{\bar{p} = (p_{1,app} m_1 +
#' p_{2,app} m_2)/(m_1 + m_2)} is the pooled prevalence the test uses,
#' \eqn{z_\beta} = \code{qnorm(power)} (0.84 for 80 percent power), and
#' \eqn{z_\alpha} = \code{qnorm(1 - alpha/2)} for a two-sided test (1.96 at
#' 95 percent) or \code{qnorm(1 - alpha)} for a one-sided test (1.64 at 95
#' percent), with \eqn{\alpha = 1 -} \code{conf_level}. With no clustering
#' and two groups of \eqn{n} people this is the formula of R's
#' \code{power.prop.test()}.
#'
#' Because \eqn{p_{1,app} - p_{2,app} = (p_1 - p_2)(Se + Sp - 1)}, an
#' imperfect test shrinks the difference the test can see and increases the
#' sample size.
#'
#' \strong{Effective size.} For group \eqn{g} with \eqn{n_g} people in sites
#' of \eqn{s_g} people, \eqn{m_g = n_g / D_g} with design effect
#' \eqn{D_g = 1 + (s_g - 1)\,ICC_g}. With one unknown \eqn{x}, each design
#' gives \eqn{m_g} directly:
#' \describe{
#'   \item{No clustering}{\eqn{m_g = x}, people per group.}
#'   \item{Fixed \code{n_per_site}}{\eqn{m_g = x \, s_g / D_g}, with \eqn{x}
#'     sites per group.}
#'   \item{Fixed \code{n_sites}}{\eqn{m_g = k_g x / (1 + (x - 1) ICC_g)},
#'     with \eqn{k_g} sites and \eqn{x} people per site.}
#' }
#' Each \eqn{m_g} grows with \eqn{x}, so the left side of the equation
#' above falls as \eqn{x} grows, and the function finds the smallest
#' \eqn{x} that meets it with a single root search (\code{uniroot()}),
#' then rounds up to a whole number.
#'
#' \strong{Highest power.} With \code{n_sites} fixed, adding people to each
#' site cannot raise the power past a limit: each site is worth at most
#' \eqn{1/ICC} independent people, so group \eqn{g} never counts as more than
#' \eqn{k_g / ICC_g} independent people. For example, with ICC 0.05 each
#' site is worth at most 20 independent people, so 8 sites per group give at
#' most 160 independent people; for 10 versus 20 percent that gives 70.9
#' percent power, so 80 percent can never be reached. When the target is out
#' of reach, the function stops and reports the power at that limit.
#'
#' \strong{Rounding.} The number of sites, or people per site, is rounded up
#' to a whole number, the same in both groups. With more sites than people
#' needed, 1 person is sampled at each of just enough of the sites (there is
#' then no clustering); with \code{n_per_site} larger than needed, 1 site
#' per group is used. The people in each group are
#' \code{n_sites * n_per_site}, and the design effects are those of the
#' rounded design.
#'
#' @return A list with:
#'   \item{n_per_group}{People in each group, \code{c(group 1, group 2)}.
#'     With sites, \code{n_sites * n_per_site}.}
#'   \item{n}{Total number of people, \code{sum(n_per_group)}.}
#'   \item{n_sites}{Sites in each group (\code{NULL} with no clustering).}
#'   \item{n_per_site}{People per site in each group (\code{NULL} with no
#'     clustering).}
#'   \item{deff}{Design effect of each group (1 with no clustering).}
#'   \item{n_eff}{People each group would need with no clustering, keeping
#'     the groups in the same ratio as the design.}
#'
#' @references
#' MMS-SD Study Design Workshop. \url{https://mrc-ide.github.io/MMS-SD_workshop/}
#'
#' Fleiss, J. L., Levin, B. & Paik, M. C. (2003) \emph{Statistical Methods for
#' Rates and Proportions}, 3rd ed. Wiley. (Two-group sample size.)
#'
#' Kish, L. (1965) \emph{Survey Sampling}. Wiley. (Design effect.)
#'
#' Rogan, W. J. & Gladen, B. (1978) Estimating prevalence from the results
#' of a screening test. \emph{American Journal of Epidemiology}
#' \strong{107}(1), 71-76. \doi{10.1093/oxfordjournals.aje.a112510}
#'
#' @seealso \code{\link{test_difference}}
#'
#' @export
#'
#' @examples
#' # Detect 10 versus 20 percent prevalence with 80 percent power
#' design_difference(prevalence1 = 0.10, prevalence2 = 0.20)
#'
#' # Higher power needs more people
#' design_difference(0.10, 0.20, power = 0.90)
#'
#' # Imperfect diagnostic test
#' design_difference(0.10, 0.20, sensitivity = 0.9, specificity = 0.95)
#'
#' # Fixed number of people per site: how many sites per group?
#' design_difference(0.10, 0.20, n_per_site = 30, icc = 0.02)
#'
#' # Fixed number of sites: how many people per site?
#' design_difference(0.10, 0.20, n_sites = 30, icc = 0.02)
#'
#' # Twice as many sites in group 2, and more clustering there
#' design_difference(0.10, 0.20, n_sites = c(20, 40), icc = c(0.02, 0.05))
design_difference <- function(prevalence1,
                              prevalence2,
                              power       = 0.80,
                              alternative = c("two.sided", "greater", "less"),
                              sensitivity = 1,
                              specificity = 1,
                              conf_level  = 0.95,
                              n_sites     = NULL,
                              n_per_site  = NULL,
                              icc         = 0) {

  .check_number(prevalence1, "prevalence1", 0, 1)
  .check_number(prevalence2, "prevalence2", 0, 1)
  .check_number(power, "power", 0.5, 1, upper_open = TRUE)
  .check_number(conf_level, "conf_level", 0, 1, lower_open = TRUE, upper_open = TRUE)
  .check_test(sensitivity, specificity)
  # n_sites, n_per_site and icc: one value for both groups, or two
  for (nm in c("n_sites", "n_per_site", "icc"))
    if (!is.null(get(nm)) && !length(get(nm)) %in% 1:2)
      stop("`", nm, "` must have one value, or two (one per group).", call. = FALSE)
  for (v in n_sites)    .check_count(v, "n_sites")
  for (v in n_per_site) .check_count(v, "n_per_site")
  for (v in icc)        .check_number(v, "icc", 0, 1)
  .check_sites(n_sites[1], n_per_site[1], max(icc))
  alternative <- match.arg(alternative)
  if (prevalence1 == prevalence2)
    stop("`prevalence1` and `prevalence2` must differ.", call. = FALSE)
  if (alternative == "greater" && prevalence1 < prevalence2)
    stop("`prevalence1` must be above `prevalence2` for alternative = 'greater'.",
         call. = FALSE)
  if (alternative == "less" && prevalence1 > prevalence2)
    stop("`prevalence1` must be below `prevalence2` for alternative = 'less'.",
         call. = FALSE)

  p     <- .apparent_prev(c(prevalence1, prevalence2), sensitivity, specificity)
  alpha <- 1 - conf_level
  z_a   <- stats::qnorm(1 - if (alternative == "two.sided") alpha / 2 else alpha)
  z_b   <- stats::qnorm(power)
  icc   <- rep_len(icc, 2)

  # Effective sizes c(m1, m2) for a value x of the unknown (see Details)
  if (!is.null(n_sites)) {
    k    <- rep_len(n_sites, 2)
    m_of <- function(x) k * x / .deff(x, icc)
  } else if (!is.null(n_per_site)) {
    s    <- rep_len(n_per_site, 2)
    m_of <- function(x) x * s / .deff(s, icc)
  } else {
    m_of <- function(x) c(x, x)
  }

  # With n_sites fixed, the effective sizes level off at k / icc as people
  # per site grow: stop if even that limit falls short of the target.
  if (!is.null(n_sites) && .diff_gap(k / icc, p, z_a, z_b) > 0)
    stop(sprintf(paste0("Target power of %.0f%% cannot be reached with these ",
                        "sites and ICC; the highest power they can give is ",
                        "%.1f%%. Increase `n_sites` or lower the ICC."),
                 100 * power, 100 * .diff_power(k / icc, p, z_a)), call. = FALSE)

  # smallest x meeting the target, rounded up, and the unclustered size
  x     <- .diff_solve(m_of, p, z_a, z_b)
  ratio <- if (!is.null(n_sites)) k else if (!is.null(n_per_site)) s else c(1, 1)
  x0    <- .diff_solve(function(x) ratio * x, p, z_a, z_b, whole = FALSE)

  if (!is.null(n_sites)) {
    sites <- k
    per   <- rep(x, 2)
    # Fewer people needed than there are sites (k * x0 people unclustered,
    # x0 < 1): with 1 person per site there is no clustering, so use 1
    # person at each of just enough of the sites.
    if (x0 < 1) sites <- ceiling(k * x0)
  } else if (!is.null(n_per_site)) {
    sites <- rep(x, 2)
    per   <- s
  } else {
    sites <- per <- NULL
  }
  n_pg <- if (is.null(sites)) rep(x, 2) else sites * per

  list(
    n_per_group = n_pg,
    n           = sum(n_pg),
    n_sites     = sites,
    n_per_site  = per,
    deff        = if (is.null(per)) c(1, 1) else .deff(per, icc),
    n_eff       = ceiling(ratio * x0)
  )
}

# How far effective sizes m = c(m1, m2) fall short of the target power: the
# left side of the sample-size equation (see Details) minus the right side.
# Negative when the target is met. An infinite m (no limit) adds no
# variance, and the pooled prevalence then is that group's.
.diff_gap <- function(m, p, z_a, z_b) {
  w   <- if (any(is.infinite(m))) as.numeric(is.infinite(m)) else m
  pb  <- sum(p * w) / sum(w)
  inv <- ifelse(is.infinite(m), 0, 1 / m)
  z_a * sqrt(pb * (1 - pb) * sum(inv)) + z_b * sqrt(sum(p * (1 - p) * inv)) -
    abs(p[1] - p[2])
}

# Power of the score test at effective sizes m (inverting the equation).
.diff_power <- function(m, p, z_a) {
  w   <- if (any(is.infinite(m))) as.numeric(is.infinite(m)) else m
  pb  <- sum(p * w) / sum(w)
  inv <- ifelse(is.infinite(m), 0, 1 / m)
  stats::pnorm((abs(p[1] - p[2]) - z_a * sqrt(pb * (1 - pb) * sum(inv))) /
                 sqrt(sum(p * (1 - p) * inv)))
}

# Smallest x with .diff_gap(m_of(x)) <= 0. The gap falls as x grows, so a
# single root search finds it; with whole = TRUE it is the smallest whole
# number, at least 1.
.diff_solve <- function(m_of, p, z_a, z_b, whole = TRUE) {
  gap <- function(x) .diff_gap(m_of(x), p, z_a, z_b)
  x   <- stats::uniroot(gap, c(1e-8, 1), extendInt = "downX", tol = 1e-10)$root
  if (!whole) return(x)
  x <- max(1, ceiling(x))
  # uniroot stops within its tolerance of the root, so check one below
  if (x > 1 && gap(x - 1) <= 0) x - 1 else x
}
