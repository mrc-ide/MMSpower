#' Calculate sample size to detect at least one case
#'
#' @description
#' Given the expected true prevalence of a (possibly rare) variant or
#' condition, returns the number of people needed so that the chance of
#' finding at least one case reaches a target (`detection_prob`). Supply `n`
#' instead to run it in reverse: the function then reports the chance of
#' finding at least one case with that many people.
#'
#' **Three design modes** -- controlled by `n_sites` and `n_per_site`:
#' \describe{
#'   \item{Unclustered (`n_sites = NULL`, `n_per_site = NULL`)}{Returns the
#'     number of people.}
#'   \item{Fixed site size (`n_per_site` supplied)}{Solves for the number of
#'     sites.}
#'   \item{Fixed number of sites (`n_sites` supplied)}{Solves for the number
#'     of people per site.}
#' }
#'
#' @param prevalence Expected true prevalence, in (0, 1): the lowest
#'   prevalence the study should be able to detect.
#' @param detection_prob Target chance of finding at least one case, in
#'   (0, 1); default 0.80. Not used when `n` is supplied.
#' @param n Optional whole number. If supplied, the function runs in reverse
#'   and reports the chance of finding at least one case among `n` people.
#' @param sensitivity Diagnostic sensitivity in (0, 1]; default 1 (the test
#'   never misses a case).
#' @param n_sites Optional whole number. Fixes the number of sites; the
#'   function solves for people per site. Cannot be used with `n_per_site`.
#' @param n_per_site Optional whole number. Fixes the number of people per
#'   site; the function solves for the number of sites. Cannot be used with
#'   `n_sites`.
#' @param icc Intra-cluster correlation, in \[0, 1); default 0 (no
#'   clustering). An `icc` above 0 needs `n_sites` or `n_per_site`.
#'
#' @details
#' **Sample size.** Each person tested has a chance
#' \eqn{q = prevalence \times sensitivity} of being a true case whom the
#' test detects. The chance that no case is found among \eqn{n} people is
#' \eqn{(1 - q)^n}, so the number of people needed is
#'
#' \deqn{n_0 = \frac{\log(1 - detection\_prob)}{\log(1 - q)}}
#'
#' rounded up (Cannon 2001, equation 4.6). Positives are assumed to be
#' confirmed, for example by a second test (Cannon 2001), so a false
#' positive is not counted as finding a case, and specificity is not used.
#'
#' **Clustering.** Cases bunch together at sites, so the chance \eqn{q} is
#' taken to vary between sites as a beta distribution with mean \eqn{q} and
#' intra-cluster correlation ICC, as in DRpower. The chance that no case is
#' found among \eqn{m} people at one site is
#'
#' \deqn{P_0(m) = \frac{B(\alpha, \beta + m)}{B(\alpha, \beta)}, \qquad
#'   \alpha = q \left(\frac{1}{ICC} - 1\right), \quad
#'   \beta = (1 - q) \left(\frac{1}{ICC} - 1\right),}
#'
#' and at \eqn{k} sites it is \eqn{P_0(m)^k}. With `n_per_site` fixed, the
#' number of sites is
#' \eqn{\lceil \log(1 - detection\_prob) / \log P_0(n_{per site}) \rceil}.
#' With `n_sites` fixed there is no formula for the people per site, so
#' the smallest whole number that reaches `detection_prob` is found by a
#' numerical search. \eqn{P_0(m)} always falls towards 0 as \eqn{m} grows,
#' but very slowly when the ICC is high and cases are rare; if even
#' \eqn{10^{15}} people per site fall short, the function stops.
#'
#' **Sites and rounding.** With `n_sites` fixed, every site is used, with at
#' least 1 person each, even if fewer people would do (as in the other
#' design functions). A site larger than needed gives 1 site of that size.
#' With sites, the returned `n` is the number actually sampled,
#' `n_sites * n_per_site`.
#'
#' @return A list with:
#'   \item{n}{Total sample size (or `n` as supplied, in reverse). With
#'     sites, `n_sites * n_per_site`.}
#'   \item{n_sites}{Number of sites (`NULL` if unclustered).}
#'   \item{n_per_site}{People per site (`NULL` if unclustered).}
#'   \item{deff}{Design effect: `n` divided by the number of people with no
#'     clustering who give the same chance of finding a case; 1 if
#'     unclustered.}
#'   \item{n_eff}{Number of people needed with no clustering, \eqn{n_0}
#'     rounded up (in reverse, the number with no clustering who give the
#'     same chance of finding a case).}
#'   \item{detection_prob}{The target (forward), or the chance achieved
#'     with `n` people (reverse).}
#'
#' @references
#' Cannon, R. M. (2001) Sense and sensitivity: designing surveys based on an
#' imperfect test. *Preventive Veterinary Medicine* **49**(3-4), 141-163.
#' \doi{10.1016/S0167-5877(01)00184-2} (Section 4.)
#'
#' Verity, R. & Ruybal, S. DRpower: Study design and analysis for pfhrp2/3
#' deletion prevalence studies. R package version 1.0.4 (functions
#' `get_power_presence()` and `get_sample_size_presence()`).
#' \url{https://github.com/mrc-ide/DRpower}
#'
#' MMS-SD Study Design Workshop, Module 7.
#' \url{https://mrc-ide.github.io/MMS-SD_workshop/}
#'
#' @export
#'
#' @examples
#' # 80 percent chance of finding a case at 5 percent prevalence
#' design_detection(prevalence = 0.05)
#'
#' # Imperfect test
#' design_detection(prevalence = 0.05, sensitivity = 0.8)
#'
#' # Reverse: the chance of finding a 1 percent variant with 139 people
#' design_detection(prevalence = 0.01, n = 139)
#'
#' # Fixed site size: how many sites?
#' design_detection(prevalence = 0.05, n_per_site = 10, icc = 0.05)
#'
#' # Fixed number of sites: how many people per site?
#' design_detection(prevalence = 0.01, n_sites = 5, icc = 0.1)
design_detection <- function(prevalence,
                             detection_prob = 0.80,
                             n              = NULL,
                             sensitivity    = 1,
                             n_sites        = NULL,
                             n_per_site     = NULL,
                             icc            = 0) {

  .check_number(prevalence, "prevalence", 0, 1, lower_open = TRUE, upper_open = TRUE)
  if (is.null(n))
    .check_number(detection_prob, "detection_prob", 0, 1,
                  lower_open = TRUE, upper_open = TRUE)
  .check_number(sensitivity, "sensitivity", 0, 1, lower_open = TRUE)
  .check_number(icc, "icc", 0, 1, upper_open = TRUE)
  .check_sites(n_sites, n_per_site, icc)
  .check_count(n, "n")

  q <- prevalence * sensitivity
  # log chance of no case at k sites of m people, and the number of people
  # with no clustering who have the same chance
  log_p0  <- function(k, m) k * .log_p0_site(m, q, icc)
  n_equiv <- function(lp0) lp0 / log1p(-q)

  # ---- reverse: chance of detection with n people ----
  if (!is.null(n)) {
    if (!is.null(n_sites) && n < n_sites)
      stop("`n` must be at least `n_sites`.", call. = FALSE)
    # people spread as evenly as the sites allow
    lp0 <- if (!is.null(n_sites))
      log_p0(n %% n_sites, n %/% n_sites + 1) +
        log_p0(n_sites - n %% n_sites, n %/% n_sites)
    else if (!is.null(n_per_site))
      log_p0(n %/% n_per_site, n_per_site) + log_p0(1, n %% n_per_site)
    else
      log_p0(1, n)
    design <- .allocate(n, n_sites, n_per_site)
    return(list(
      n              = n,
      n_sites        = design$n_sites,
      n_per_site     = design$n_per_site,
      deff           = n / n_equiv(lp0),
      n_eff          = ceiling(n_equiv(lp0)),
      detection_prob = -expm1(lp0)
    ))
  }

  # ---- forward ----
  log_miss <- log1p(-detection_prob)   # log chance of missing every case allowed
  n0       <- ceiling(log_miss / log1p(-q))
  if (!is.null(n_per_site)) {
    # every site misses every case with the same chance
    design <- list(n_sites    = ceiling(log_miss / log_p0(1, n_per_site)),
                   n_per_site = n_per_site)
  } else if (!is.null(n_sites)) {
    design <- list(n_sites    = n_sites,
                   n_per_site = .detection_per_site(n_sites, log_p0, log_miss))
  } else {
    design <- list(n_sites = NULL, n_per_site = NULL)
  }
  n_total <- if (is.null(n_sites) && is.null(n_per_site)) n0
             else design$n_sites * design$n_per_site

  list(
    n              = n_total,
    n_sites        = design$n_sites,
    n_per_site     = design$n_per_site,
    deff           = if (is.null(design$n_sites)) 1 else
      n_total / n_equiv(log_p0(design$n_sites, design$n_per_site)),
    n_eff          = n0,
    detection_prob = detection_prob
  )
}

# Log chance that no case is found among m people at one site, when the
# chance q varies between sites as a beta distribution with mean q and
# intra-cluster correlation icc (DRpower's get_power_presence()); with icc
# 0 every person is independent. m need not be a whole number.
.log_p0_site <- function(m, q, icc) {
  if (icc == 0) return(m * log1p(-q))
  a <- q * (1 / icc - 1)
  b <- (1 - q) * (1 / icc - 1)
  lbeta(a, b + m) - lbeta(a, b)
}

# Smallest whole number of people per site at n_sites sites with log chance
# of no case at most log_miss. That chance falls as people are added, so one
# root search on log(m) finds it. It falls very slowly when the ICC is high
# and cases are rare, so the search stops at 1e15 people per site (beyond
# which whole numbers are not exact).
.detection_per_site <- function(n_sites, log_p0, log_miss) {
  f <- function(u) log_p0(n_sites, exp(u)) - log_miss
  if (f(0) <= 0) return(1)
  u_max <- log(1e15)
  if (f(u_max) > 0)
    stop("`detection_prob` cannot be reached with ", n_sites, " sites at this ",
         "ICC: even 10^15 people per site give ",
         format(100 * -expm1(log_p0(n_sites, 1e15)), digits = 3),
         " percent. Increase `n_sites` or use a lower `icc`.", call. = FALSE)
  m <- ceiling(exp(stats::uniroot(f, c(0, u_max), tol = 1e-10)$root))
  # uniroot stops within its tolerance of the root, so check one below
  if (m > 1 && f(log(m - 1)) <= 0) m - 1 else m
}
