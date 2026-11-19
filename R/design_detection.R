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
#' @param icc Intra-cluster correlation, in \[0, 1\]; default 0 (no
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
#' and at \eqn{k} sites it is \eqn{P_0(m)^k}. The function finds the
#' smallest number of sites (or people per site) that reaches
#' `detection_prob`. With an ICC of 1 everyone at a site has the same
#' result, so \eqn{k} sites give at most \eqn{1 - (1 - q)^k}. An ICC below
#' about 1.5e-8 is treated as 0, where the beta formula loses accuracy.
#'
#' **Too many sites.** One person per site gives the chance of `n_sites`
#' people with no clustering. If fewer people are needed, fewer sites are
#' used, with 1 person at each, and a warning says so.
#'
#' **Limit and rounding.** A design that needs more than 2 million people
#' stops with an error. With sites, the returned `n` is the number actually
#' sampled, `n_sites * n_per_site`.
#'
#' @return A list with:
#'   \item{n}{Total sample size (or `n` as supplied, in reverse). With
#'     sites, `n_sites * n_per_site`.}
#'   \item{n_sites}{Number of sites (`NULL` if unclustered).}
#'   \item{n_per_site}{People per site (`NULL` if unclustered).}
#'   \item{deff}{Design effect: `n` divided by the number of people with no
#'     clustering who give the same chance of finding a case; 1 if
#'     unclustered.}
#'   \item{n_eff}{Unclustered sample size \eqn{n_0} in a very large
#'     population, rounded up.}
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
  .check_sites(n_sites, n_per_site, icc)
  .check_count(n, "n")

  q <- prevalence * sensitivity
  # Below this the beta formula loses accuracy, and the answer is the same
  # as with no clustering.
  clustered <- icc >= sqrt(.Machine$double.eps)

  # Tiny allowance for rounding noise, so a target exactly equal to the
  # chance of a design is reached by that design (not one person more).
  tol <- 1e-14

  # log chance of no case among people at sites; log1p keeps a tiny q from
  # being rounded away (1 - 1e-17 is stored as exactly 1)
  log_p0_sites <- function(k, m) k * .log_p0_site(m, q, icc)
  n_equiv      <- function(log_p0) log_p0 / log1p(-q)   # unclustered people

  # ---- reverse: chance of detection with n people ----
  if (!is.null(n)) {
    if (!clustered) {
      detection <- -expm1(n * log1p(-q))
      n_eff     <- n
    } else {
      if (!is.null(n_sites) && n < n_sites)
        stop("`n` must be at least `n_sites`.", call. = FALSE)
      # people spread as evenly as the sites allow
      log_p0 <- if (!is.null(n_sites))
        log_p0_sites(n %% n_sites, n %/% n_sites + 1) +
          log_p0_sites(n_sites - n %% n_sites, n %/% n_sites)
      else
        log_p0_sites(n %/% n_per_site, n_per_site) +
          (if (n %% n_per_site > 0) log_p0_sites(1, n %% n_per_site) else 0)
      detection <- -expm1(log_p0)
      n_eff     <- n_equiv(log_p0)
    }
    design <- .allocate(n, n_sites, n_per_site)
    return(list(
      n              = n,
      n_sites        = design$n_sites,
      n_per_site     = design$n_per_site,
      deff           = n / n_eff,
      # the tiny allowance stops rounding noise adding a person
      n_eff          = ceiling(n_eff * (1 - tol)),
      detection_prob = detection
    ))
  }

  # ---- forward: unclustered sample size ----
  n0 <- log1p(-detection_prob) / log1p(-q)
  too_many <- function(needed)
    stop("These settings need more than 2 million people (about ",
         format(signif(ceiling(needed), 7), big.mark = ",", scientific = FALSE),
         ").", call. = FALSE)
  if (n0 > 2e6) too_many(n0)
  needed <- ceiling(n0)
  # One fewer is enough when rounding noise pushed n0 just above a whole number.
  if (needed > 1 && (needed - 1) * log1p(-q) <= log1p(-detection_prob) * (1 - tol))
    needed <- needed - 1

  # more sites than people needed: 1 person at just enough sites
  if (!is.null(n_sites) && n_sites > needed) {
    warning("The study needs ", needed, " people, fewer than the ", n_sites,
            " sites given; the design uses ", needed, " sites with 1 person ",
            "at each.", call. = FALSE)
    n_sites <- needed
  }

  # ---- clustered sample size ----
  deff <- 1
  if (clustered) {
    # log chance of missing every case that the target allows (the tiny
    # allowance stops rounding noise rejecting an exact hit)
    log_miss <- log1p(-detection_prob) * (1 - tol)
    if (!is.null(n_per_site)) {
      # each site misses every case with the same chance; a site larger than
      # needed gives 1 site (as in design_precision())
      k <- ceiling(log_miss / log_p0_sites(1, n_per_site))
      design <- list(n_sites = k, n_per_site = n_per_site)
    } else {
      # smallest number per site, by bisection, up to 2 million people
      reaches <- function(m) log_p0_sites(n_sites, m) <= log_miss
      m_max   <- floor(2e6 / n_sites)
      if (!reaches(m_max)) {
        if (icc < 1)
          stop("These settings need more than 2 million people.", call. = FALSE)
        stop("`detection_prob` cannot be reached with ", n_sites, " sites and ",
             "ICC = 1: they give at most ",
             format(100 * -expm1(log_p0_sites(n_sites, 1)), digits = 4),
             " percent, however many people each site has.", call. = FALSE)
      }
      lo <- 0; hi <- m_max
      while (hi - lo > 1) {
        mid <- floor((lo + hi) / 2)
        if (reaches(mid)) hi <- mid else lo <- mid
      }
      design <- list(n_sites = n_sites, n_per_site = hi)
    }
    design$n <- design$n_sites * design$n_per_site
    if (design$n > 2e6) too_many(design$n)
    deff <- design$n /
      n_equiv(log_p0_sites(design$n_sites, design$n_per_site))
  } else {
    design <- .allocate(needed, n_sites, n_per_site)
  }

  list(
    n              = design$n,
    n_sites        = design$n_sites,
    n_per_site     = design$n_per_site,
    deff           = deff,
    n_eff          = needed,
    detection_prob = detection_prob
  )
}

# Log chance that no case is found among m people at one site, when the
# chance q varies between sites as a beta distribution with mean q and
# intra-cluster correlation icc (DRpower's get_power_presence()). With icc 1
# everyone at a site has the same result, so it is log(1 - q) for any m.
.log_p0_site <- function(m, q, icc) {
  if (icc == 1) return(log1p(-q))
  a <- q * (1 / icc - 1)
  b <- (1 - q) * (1 / icc - 1)
  lbeta(a, b + m) - lbeta(a, b)
}
