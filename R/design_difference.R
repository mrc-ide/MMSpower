#' Calculate sample size to detect a difference in prevalence between two groups
#'
#' @description
#' Design function: given the expected true prevalence in two groups, returns
#' the sample size needed in each group for \code{test_difference()} (with
#' its default score test) to detect that difference with the target power;
#' or, given a full design for both groups, returns the power it gives.
#'
#' Each group has its own arguments, \code{n1}, \code{n_sites1},
#' \code{n_per_site1} for group 1 and \code{n2}, \code{n_sites2},
#' \code{n_per_site2} for group 2, as in \code{test_difference()}. Each
#' group-2 argument defaults to its group-1 value, so for two groups with the
#' same design only the group-1 arguments are needed.
#'
#' \strong{What is solved for} -- leave out the one quantity you want, the
#' same in both groups:
#' \describe{
#'   \item{No sites}{Returns the number of people per group.}
#'   \item{\code{n_per_site1} (and \code{n_per_site2})}{Returns the number
#'     of sites per group, the same in both.}
#'   \item{\code{n_sites1} (and \code{n_sites2})}{Returns the number of
#'     people per site, the same in both.}
#'   \item{A full design for both groups}{Returns the power. For each group
#'     give \code{n} (no sites), \code{n_sites} and \code{n_per_site}, or a
#'     vector \code{n_per_site} with one size per site (for sites of
#'     different sizes). Leave \code{power} at its default.}
#' }
#' When solving, the groups' sizes follow their fixed values: for example,
#' \code{n_sites1 = 10, n_sites2 = 20} gives group 2 twice as many people as
#' group 1.
#'
#' @param prevalence1 Expected true prevalence in group 1, in \[0, 1\].
#' @param prevalence2 Expected true prevalence in group 2, in \[0, 1\].
#'   The two prevalences must differ.
#' @param power Target power, in \[0.5, 1); default 0.80. Not used, and not
#'   to be given, with a full design.
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
#' @param conf_level Confidence level, in \[0.5, 1); default 0.95. The
#'   significance level of the test is 1 - \code{conf_level}. Below 0.5 the
#'   one-sided z value turns negative and the sample size shrinks towards
#'   nonsense, so it is not allowed.
#' @param n1,n2 Optional whole numbers: the number of people in each group,
#'   for a full design with no sites. \code{n2} defaults to \code{n1}.
#' @param n_sites1,n_sites2 Optional whole numbers: the number of sites in
#'   each group. Alone, the function solves for the number of people per
#'   site. \code{n_sites2} defaults to \code{n_sites1}.
#' @param n_per_site1,n_per_site2 Optional: people per site in each group.
#'   One whole number (alone, the function solves for the number of sites),
#'   or a vector with one size per site, a full design with sites of
#'   different sizes (then leave out that group's \code{n_sites}).
#'   \code{n_per_site2} defaults to \code{n_per_site1}.
#' @param icc Intra-cluster correlation, in \[0, 1\]: one number for both
#'   groups, or two, \code{c(group 1, group 2)}. Default 0 (no clustering).
#'   An \code{icc} above 0 needs sites in that group.
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
#' \eqn{D_g = 1 + (s_g - 1)\,ICC_g}. For sites of different sizes, \eqn{s_g}
#' is the size-weighted mean site size \eqn{\sum s^2 / \sum s}, as in
#' \code{test_difference()}. With one unknown \eqn{x}, each design gives
#' \eqn{m_g} directly:
#' \describe{
#'   \item{No clustering}{\eqn{m_g = x}, people per group.}
#'   \item{Fixed \code{n_per_site1}, \code{n_per_site2}}{\eqn{m_g = x \,
#'     s_g / D_g}, with \eqn{x} sites per group.}
#'   \item{Fixed \code{n_sites1}, \code{n_sites2}}{\eqn{m_g = k_g x / (1 +
#'     (x - 1) ICC_g)}, with \eqn{k_g} sites and \eqn{x} people per site.}
#' }
#' Each \eqn{m_g} grows with \eqn{x}, so the left side of the equation
#' above falls as \eqn{x} grows, and the function finds the smallest
#' \eqn{x} that meets it with a single root search (\code{uniroot()}),
#' then rounds up to a whole number.
#'
#' \strong{Power.} For a full design, or the rounded design that was solved
#' for, the power is the equation above solved for \eqn{z_\beta}:
#' \eqn{\Phi((|p_{1,app} - p_{2,app}| - z_\alpha \sqrt{\bar{p}(1-\bar{p})
#' (1/m_1 + 1/m_2)}) / \sqrt{p_{1,app}(1-p_{1,app})/m_1 +
#' p_{2,app}(1-p_{2,app})/m_2})}.
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
#' to a whole number, the same in both groups. With \code{n_sites} fixed,
#' every site is used, with at least 1 person each, even if fewer people
#' would do (as in the other design functions); with \code{n_per_site}
#' larger than needed, 1 site per group is used. The people in each group are
#' \code{n_sites * n_per_site}, and the design effects are those of the
#' rounded design.
#'
#' @return A list with:
#'   \item{n_per_group}{People in each group, \code{c(group 1, group 2)}.
#'     With sites, \code{n_sites * n_per_site}.}
#'   \item{n}{Total number of people, \code{sum(n_per_group)}.}
#'   \item{n_sites}{Sites in each group (\code{NULL} with no clustering).}
#'   \item{n_per_site}{People per site in each group (\code{NULL} with no
#'     clustering); a list of each group's site sizes if they differ within
#'     a group.}
#'   \item{deff}{Design effect of each group (1 with no clustering).}
#'   \item{n_eff}{People each group would need with no clustering, keeping
#'     the groups in the same ratio as the design; for a full design,
#'     \code{n_per_group / deff}.}
#'   \item{power}{Power of the returned design.}
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
#' design_difference(0.10, 0.20, n_per_site1 = 30, icc = 0.02)
#'
#' # Fixed number of sites: how many people per site?
#' design_difference(0.10, 0.20, n_sites1 = 30, icc = 0.02)
#'
#' # Twice as many sites in group 2, and more clustering there
#' design_difference(0.10, 0.20, n_sites1 = 20, n_sites2 = 40,
#'                   icc = c(0.02, 0.05))
#'
#' # Full design: the power of 200 people per group
#' design_difference(0.10, 0.20, n1 = 200)
#'
#' # Full design with sites of different sizes in each group
#' design_difference(0.10, 0.20,
#'                   n_per_site1 = c(80, 60, 40, 40), n_per_site2 = c(100, 50, 50),
#'                   icc = 0.02)
design_difference <- function(prevalence1,
                              prevalence2,
                              power       = 0.80,
                              alternative = c("two.sided", "greater", "less"),
                              sensitivity = 1,
                              specificity = 1,
                              conf_level  = 0.95,
                              n1          = NULL,
                              n2          = n1,
                              n_sites1    = NULL,
                              n_sites2    = n_sites1,
                              n_per_site1 = NULL,
                              n_per_site2 = n_per_site1,
                              icc         = 0) {

  .check_number(prevalence1, "prevalence1", 0, 1)
  .check_number(prevalence2, "prevalence2", 0, 1)
  .check_number(conf_level, "conf_level", 0.5, 1, upper_open = TRUE)
  .check_test(sensitivity, specificity)
  if (!length(icc) %in% 1:2)
    stop("`icc` must have one value, or two (one per group).", call. = FALSE)
  for (v in icc) .check_number(v, "icc", 0, 1)
  icc <- rep_len(icc, 2)
  g1  <- .site_design(n1, n_sites1, n_per_site1, icc[1], "1")
  g2  <- .site_design(n2, n_sites2, n_per_site2, icc[2], "2")
  if (g1$full != g2$full)
    stop("Give a full design for both groups, or for neither.", call. = FALSE)
  .check_target(g1$full, !missing(power), "power")
  if (!g1$full) .check_number(power, "power", 0.5, 1, upper_open = TRUE)
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

  # ---- full design: its power ----
  if (g1$full) {
    if (is.null(g1$sizes) != is.null(g2$sizes))
      stop("Give sites for both groups, or for neither.", call. = FALSE)
    n_pg <- c(g1$n, g2$n)
    deff <- .deff(c(g1$n_w, g2$n_w), icc)
    r1   <- .report_sizes(g1$sizes)
    r2   <- .report_sizes(g2$sizes)
    per  <- if (is.null(g1$sizes)) NULL
            else if (length(r1$n_per_site) == 1 && length(r2$n_per_site) == 1)
              c(r1$n_per_site, r2$n_per_site)
            else list(r1$n_per_site, r2$n_per_site)
    return(list(
      n_per_group = n_pg,
      n           = sum(n_pg),
      n_sites     = c(r1$n_sites, r2$n_sites),
      n_per_site  = per,
      deff        = deff,
      n_eff       = n_pg / deff,
      power       = .diff_power(n_pg / deff, p, z_a)
    ))
  }

  # ---- solve for the one unknown, the same in both groups ----
  mode <- function(g_sites, g_per)
    if (!is.null(g_sites)) "sites" else if (!is.null(g_per)) "per" else "none"
  if (mode(n_sites1, n_per_site1) != mode(n_sites2, n_per_site2))
    stop("Both groups must leave out the same thing: fix `n_sites1` and ",
         "`n_sites2`, or `n_per_site1` and `n_per_site2`, or neither.",
         call. = FALSE)
  z_b <- stats::qnorm(power)

  # Effective sizes c(m1, m2) for a value x of the unknown (see Details)
  if (!is.null(n_sites1)) {
    k    <- c(n_sites1, n_sites2)
    m_of <- function(x) k * x / .deff(x, icc)
  } else if (!is.null(n_per_site1)) {
    s    <- c(n_per_site1, n_per_site2)
    m_of <- function(x) x * s / .deff(s, icc)
  } else {
    m_of <- function(x) c(x, x)
  }

  # With n_sites fixed, the effective sizes level off at k / icc as people
  # per site grow: stop if even that limit falls short of the target.
  if (!is.null(n_sites1) && .diff_gap(k / icc, p, z_a, z_b) > 0)
    stop(sprintf(paste0("Target power of %.0f%% cannot be reached with these ",
                        "sites and ICC; the highest power they can give is ",
                        "%.1f%%. Increase the number of sites or lower the ICC."),
                 100 * power, 100 * .diff_power(k / icc, p, z_a)), call. = FALSE)

  # smallest x meeting the target, rounded up, and the unclustered size
  x     <- .diff_solve(m_of, p, z_a, z_b)
  ratio <- if (!is.null(n_sites1)) k else if (!is.null(n_per_site1)) s else c(1, 1)
  x0    <- .diff_solve(function(x) ratio * x, p, z_a, z_b, whole = FALSE)

  if (!is.null(n_sites1)) {
    sites <- k
    per   <- rep(x, 2)
  } else if (!is.null(n_per_site1)) {
    sites <- rep(x, 2)
    per   <- s
  } else {
    sites <- per <- NULL
  }
  n_pg <- if (is.null(sites)) rep(x, 2) else sites * per
  deff <- if (is.null(per)) c(1, 1) else .deff(per, icc)

  list(
    n_per_group = n_pg,
    n           = sum(n_pg),
    n_sites     = sites,
    n_per_site  = per,
    deff        = deff,
    n_eff       = ceiling(ratio * x0),
    power       = .diff_power(n_pg / deff, p, z_a)
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
