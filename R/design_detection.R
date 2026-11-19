#' Calculate sample size to detect at least one case
#'
#' @description
#' Design function: given the expected true prevalence of a (possibly rare)
#' variant or condition, returns the number of people needed so that the
#' chance of finding at least one case reaches a target
#' (\code{detection_prob}). Supply \code{n} instead to run it in reverse: the
#' function then reports the chance of finding at least one case
#' with that many people.
#'
#' Three designs are supported:
#' \describe{
#'   \item{No clustering (neither \code{n_sites} nor \code{n_per_site})}{
#'     Returns the number of people.}
#'   \item{Fixed number of people per site (\code{n_per_site})}{Solves for
#'     the number of sites needed.}
#'   \item{Fixed number of sites (\code{n_sites})}{Solves for the number of
#'     people needed per site.}
#' }
#'
#' @param prevalence Expected true prevalence, in (0, 1): the lowest
#'   prevalence the study should be able to detect.
#' @param detection_prob Target chance of finding at least one case, in
#'   (0, 1); default 0.80. Not used when \code{n} is supplied.
#' @param n Optional positive integer. If supplied, the function runs in
#'   reverse and reports the chance of finding at least one case
#'   among \code{n} people; \code{detection_prob} is then not used.
#' @param sensitivity Diagnostic sensitivity in (0, 1]: the share of people
#'   who truly have the condition that the test detects as positive.
#'   Default 1 (the test never misses a case).
#' @param n_sites Optional positive integer. Fixes the number of sites; the
#'   function solves for the number of people per site. Cannot be used
#'   together with \code{n_per_site}.
#' @param n_per_site Optional positive integer. Fixes the number of people
#'   per site; the function solves for the number of sites. Cannot be used
#'   together with \code{n_sites}. If one site already has more people than
#'   needed, the function stops with an error; use \code{n_sites = 1} to find
#'   how many people that one site needs.
#' @param icc Intra-cluster correlation (ICC): how alike people in the same
#'   site are, from 0 (not at all) to 1 (identical). Default 0 (no
#'   clustering). An \code{icc} above 0 needs \code{n_sites} or
#'   \code{n_per_site}, and cannot be used together with \code{fpc_N}. A
#'   value below about 1.5e-8 is treated as 0: the calculation with sites
#'   loses accuracy there, and the answer is the same as with no clustering.
#' @param fpc_N Optional. Population size for a finite-population
#'   calculation: a positive integer. \code{NULL} (default) means the
#'   population is treated as very large. Only used with \code{icc} 0.
#'
#' @details
#' \strong{Sample size.} Each person tested has a chance
#' \eqn{q = prevalence \times sensitivity} of being a true case whom the
#' test detects. The chance that no case is found among \eqn{n} people is
#' \eqn{(1 - q)^n}, so
#' the chance of finding at least one is \eqn{1 - (1 - q)^n}, and the number
#' of people needed is
#'
#' \deqn{n_0 = \frac{\log(1 - detection\_prob)}{\log(1 - q)}}
#'
#' rounded up (Cannon 2001, equation 4.6). Positives are assumed to be
#' confirmed, for example by a second test (Cannon 2001), so a false
#' positive is not counted as finding a case, and specificity is not used.
#'
#' A design that needs more than 2 million people stops with an error (as
#' in \code{design_threshold()}).
#'
#' \strong{Clustering.} People at the same site tend to resemble each other,
#' so cases bunch together: some sites have many and others none. To allow
#' for this, the chance \eqn{q} is taken to vary from site to site, with
#' average \eqn{q} and intra-cluster correlation ICC (a beta distribution, as
#' in the DRpower package). The chance that no case is found among
#' \eqn{m} people at one site is then
#'
#' \deqn{P_0(m) = \frac{B(\alpha, \beta + m)}{B(\alpha, \beta)}, \qquad
#'   \alpha = q \left(\frac{1}{ICC} - 1\right), \quad
#'   \beta = (1 - q) \left(\frac{1}{ICC} - 1\right),}
#'
#' where \eqn{B} is the beta function, and the chance of finding at least one
#' case at \eqn{k} sites of \eqn{m} people is \eqn{1 - P_0(m)^k}. With
#' \code{n_per_site} fixed, the function finds the smallest number of sites
#' that reaches \code{detection_prob}; with \code{n_sites} fixed, the
#' smallest number of people per site. With an ICC of 1 everyone at a site
#' has the same result, so \eqn{P_0(m) = 1 - q} however many people each
#' site has, and \eqn{k} sites can give at most \eqn{1 - (1 - q)^k}; the
#' function stops if the target is above that. With an ICC below 1, more
#' people per site always raise the chance, but slowly when the ICC is high.
#'
#' With \code{n_sites} fixed, one person per site already gives the chance
#' of \code{n_sites} people with no clustering. If the study needs fewer
#' people than there are sites, fewer sites are used, with 1 person at
#' each, and a warning says so (as in \code{design_threshold()}).
#'
#' \strong{Finite population.} With \code{fpc_N} (the population size
#' \eqn{N}), the population is taken to contain \eqn{D} cases, \eqn{N \times
#' prevalence} rounded down, and people are drawn without replacement. The
#' chance that no case is found among \eqn{n} people is then
#'
#' \deqn{\sum_k \frac{\binom{D}{k}\binom{N-D}{n-k}}{\binom{N}{n}} (1 - Se)^k}
#'
#' (the hypergeometric distribution: \eqn{k} cases in the sample, each missed
#' by the test with probability \eqn{1 - Se}; Cannon 2001, equation 4.1), and \eqn{n_0} is the smallest
#' \eqn{n} that reaches \code{detection_prob}. Testing everyone finds a case
#' with probability \eqn{1 - (1 - Se)^D}, so a target above that cannot be
#' reached. The calculation with clustering assumes a very large population,
#' so \code{fpc_N} cannot be used together with an \code{icc} above 0.
#'
#' @return A named list with the following fields, in the order returned:
#'   \item{n}{Number of people needed (forward), or as supplied (reverse).
#'     With sites and an \code{icc} above 0, this is \code{n_sites *
#'     n_per_site}.}
#'   \item{n_eff}{\eqn{n_0} with no finite-population calculation, rounded
#'     up: the number of people that would be needed with no clustering and
#'     a very large population. In reverse mode, the number of people with
#'     no clustering who would give the same chance of finding a case as the
#'     \code{n} people given, rounded up.}
#'   \item{n_sites}{Number of sites: solved for when \code{n_per_site} is
#'     given, echoed when \code{n_sites} is given. If the study needs fewer
#'     people than there are sites, fewer sites are used, with 1 person at
#'     each, and a warning says so. \code{NULL} with no clustering.}
#'   \item{n_per_site}{People per site: solved for when \code{n_sites} is
#'     given, echoed when \code{n_per_site} is given. \code{NULL} with no
#'     clustering. With \code{icc} 0 the design can collect more people than
#'     \code{n}, because the site numbers are rounded up.}
#'   \item{detection_prob}{Target chance of finding at least one case
#'     (forward), or the chance achieved with \code{n} people
#'     (reverse).}
#'   \item{prevalence}{The true prevalence, as supplied.}
#'   \item{q}{Chance that one person is a true case and tests positive,
#'     \code{prevalence * sensitivity}.}
#'   \item{n_cases}{Number of cases \eqn{D} assumed in the population when
#'     \code{fpc_N} is given; \code{NULL} otherwise.}
#'   \item{sensitivity}{Sensitivity, as supplied.}
#'   \item{icc}{ICC, as supplied.}
#'   \item{deff}{Design effect: 1 with no clustering. With sites and an
#'     \code{icc} above 0, \code{n} divided by the number of people with no
#'     clustering who would give the same chance of finding a case.}
#'   \item{fpc_N}{\code{fpc_N}, as supplied, or \code{NULL}.}
#'   \item{mode}{\code{"solve_n"} (forward) or \code{"solve_detection_prob"}
#'     (reverse).}
#'
#' @references
#' MMS-SD Study Design Workshop. \url{https://mrc-ide.github.io/MMS-SD_workshop/}
#'
#' Verity, R. & Ruybal, S. DRpower: Study design and analysis for pfhrp2/3
#' deletion prevalence studies. R package, version 1.0.4 (functions
#' \code{get_power_presence()} and \code{get_sample_size_presence()}).
#' \url{https://github.com/mrc-ide/DRpower}
#'
#' Cannon, R. M. (2001) Sense and sensitivity: designing surveys based on an
#' imperfect test. \emph{Preventive Veterinary Medicine} \strong{49}(3-4),
#' 141-163. \doi{10.1016/S0167-5877(01)00184-2} (Section 4: sample size
#' with an imperfect sensitivity and confirmed positives, equations 4.1 and
#' 4.6.)
#'
#' @export
#'
#' @examples
#' # How many people to have an 80 percent chance of finding at least one
#' # case when the prevalence is 5 percent?
#' design_detection(prevalence = 0.05)
#'
#' # A higher target needs more people
#' design_detection(prevalence = 0.05, detection_prob = 0.95)
#'
#' # Imperfect diagnostic test
#' design_detection(prevalence = 0.05, sensitivity = 0.8)
#'
#' # Reverse: the chance of finding a 1 percent variant with 139 people
#' design_detection(prevalence = 0.01, n = 139)
#'
#' # Population of 200 people
#' design_detection(prevalence = 0.05, fpc_N = 200)
#'
#' # Fixed number of people per site: how many sites?
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
                             icc            = 0,
                             fpc_N          = NULL) {

  solve_n <- is.null(n)

  # Numbers in messages are printed in full (100000, not 1e+05), whatever
  # the `digits` and `scipen` options.
  num <- function(v) format(v, scientific = FALSE, digits = 15)
  # ...except numbers far beyond any real study, shown as 1e+200.
  huge <- function(v) if (v < 1e16) num(v) else format(v, digits = 15)
  # A value within rounding of a whole number (such as 100 * 1.1, stored as
  # 110.00000000000001) prints as that whole number, so say why it fails.
  near_int <- function(v, nm)
    if (is.finite(v) && v != round(v) &&
        abs(v - round(v)) <= 8 * .Machine$double.eps * abs(v))
      paste0(", which is not exactly an integer; use round(", nm, ")") else ""

  # ---- validate single-value parameters ----
  # Each of the parameters below must be a single value. `detection_prob`
  # is only used, and so only checked, when solving for n.
  if (length(prevalence) != 1 || !is.numeric(prevalence))
    stop("`prevalence` must be a single number in (0, 1) (got class `",
         class(prevalence)[1], "`, length ", length(prevalence), ").")
  if (solve_n && (length(detection_prob) != 1 || !is.numeric(detection_prob)))
    stop("`detection_prob` must be a single number in (0, 1) (got class `",
         class(detection_prob)[1], "`, length ", length(detection_prob), ").")
  if (length(sensitivity) != 1 || !is.numeric(sensitivity))
    stop("`sensitivity` must be a single number in (0, 1] (got class `",
         class(sensitivity)[1], "`, length ", length(sensitivity), "). ",
         "Note: `TRUE`, `FALSE` and a plain `NA` are logical, not numeric. ",
         "Pass 1 if the test never misses a case.")
  if (length(icc) != 1 || !is.numeric(icc))
    stop("`icc` must be a single number in [0, 1] (got class `",
         class(icc)[1], "`, length ", length(icc), ").")

  # Check for NA, NaN and Inf before the range checks below.
  if (!is.finite(prevalence))
    stop("`prevalence` must be a single finite number (got ", prevalence, ").")
  if (solve_n && !is.finite(detection_prob))
    stop("`detection_prob` must be a single finite number (got ",
         detection_prob, ").")
  if (!is.finite(sensitivity))
    stop("`sensitivity` must be a single finite number (got ", sensitivity, ").")
  if (!is.finite(icc))
    stop("`icc` must be a single finite number (got ", icc, ").")

  if (prevalence < 0 || prevalence > 1)
    stop("`prevalence` must be in (0, 1) (got ", num(prevalence), "). ",
         "`prevalence` is a fraction of the population and cannot be ",
         if (prevalence < 0) "negative." else "greater than 1.")
  if (prevalence == 0 || prevalence == 1)
    stop("`prevalence` must be in (0, 1) (got ", num(prevalence), "). ",
         if (prevalence == 0) "With a prevalence of 0 there is no case to find."
         else "With a prevalence of 1 everyone is a case, so there is no sample size to plan.")
  if (solve_n && (detection_prob <= 0 || detection_prob >= 1))
    stop("`detection_prob` must be in (0, 1) (got ", num(detection_prob), "). ",
         "For example, use 0.80 for an 80 percent chance of finding at least ",
         "one case.")
  if (sensitivity < 0 || sensitivity > 1)
    stop("`sensitivity` must be in (0, 1] (got ", num(sensitivity), "). ",
         "`sensitivity` is a diagnostic probability and cannot be ",
         if (sensitivity < 0) "negative." else "greater than 1.")
  if (sensitivity == 0)
    stop("`sensitivity` must be in (0, 1] (got 0). ",
         "A sensitivity of 0 means everyone who truly has the condition ",
         "tests negative.")
  if (icc < 0 || icc > 1)
    stop("`icc` must be in [0, 1] (got ", num(icc), "). ",
         "`icc` is a correlation and cannot be ",
         if (icc < 0) "negative." else "greater than 1.")

  # ---- site and population inputs ----
  if (!solve_n &&
      (!is.numeric(n) || length(n) != 1 || !is.finite(n) ||
       n != floor(n) || n < 1))
    stop("`n` must be a single finite positive integer (got ",
         if (!is.numeric(n)) paste0("class `", class(n)[1], "`")
         else if (length(n) != 1) paste0("length = ", length(n))
         else paste0(num(n), near_int(n, "n")), "). `n` is the number of ",
         "people tested.")
  if (!solve_n && n > 2^53)
    stop("`n` (", huge(n), ") is larger than 2^53 = 9007199254740992, the ",
         "largest whole number R can store exactly. No real study tests this ",
         "many people.")
  if (!is.null(n_sites) && !is.null(n_per_site))
    stop("Supply at most one of `n_sites` or `n_per_site`, not both. ",
         "`n_sites` fixes the number of sites and solves for people per site; ",
         "`n_per_site` fixes the number of people per site and solves for the ",
         "number of sites.")
  if (!is.null(n_sites) &&
      (!is.numeric(n_sites) || length(n_sites) != 1 || !is.finite(n_sites) ||
       n_sites != floor(n_sites) || n_sites < 1))
    stop("`n_sites` must be a single finite positive integer (got ",
         if (!is.numeric(n_sites)) paste0("class `", class(n_sites)[1], "`")
         else if (length(n_sites) != 1) paste0("length = ", length(n_sites))
         else paste0(num(n_sites), near_int(n_sites, "n_sites")), "). ",
         "`n_sites` is the number of sites.")
  if (!is.null(n_per_site) &&
      (!is.numeric(n_per_site) || length(n_per_site) != 1 || !is.finite(n_per_site) ||
       n_per_site != floor(n_per_site) || n_per_site < 1))
    stop("`n_per_site` must be a single finite positive integer (got ",
         if (!is.numeric(n_per_site)) paste0("class `", class(n_per_site)[1], "`")
         else if (length(n_per_site) != 1) paste0("length = ", length(n_per_site))
         else paste0(num(n_per_site), near_int(n_per_site, "n_per_site")), "). ",
         "`n_per_site` is the fixed number of people sampled per site.")
  if (!is.null(fpc_N) &&
      (!is.numeric(fpc_N) || length(fpc_N) != 1 || !is.finite(fpc_N) ||
       fpc_N < 1 || fpc_N != floor(fpc_N)))
    stop("`fpc_N` must be a single finite positive integer (got ",
         if (!is.numeric(fpc_N)) paste0("class `", class(fpc_N)[1], "`")
         else if (length(fpc_N) != 1) paste0("length = ", length(fpc_N))
         else paste0(num(fpc_N), near_int(fpc_N, "fpc_N")),
         "). `fpc_N` is the population size. ",
         "Set `fpc_N = NULL` to treat the population as very large.")
  if (!is.null(fpc_N) && fpc_N > 2^53)
    stop("`fpc_N` (", huge(fpc_N), ") is larger than 2^53 = 9007199254740992, ",
         "the largest whole number R can store exactly. No real population is ",
         "this large: set `fpc_N = NULL`.")
  if (!is.null(n_sites) && n_sites > 2^53)
    stop("`n_sites` (", huge(n_sites), ") is larger than 2^53 = 9007199254740992, ",
         "the largest whole number R can store exactly. No real study has ",
         "this many sites.")
  if (!is.null(n_per_site) && n_per_site > 2^53)
    stop("`n_per_site` (", huge(n_per_site), ") is larger than 2^53 = ",
         "9007199254740992, the largest whole number R can store exactly. No ",
         "real site has this many people.")

  # The number of sites, or people per site, cannot exceed the population.
  if (!is.null(fpc_N) && !is.null(n_sites) && n_sites > fpc_N)
    stop("`n_sites` (", num(n_sites), ") exceeds `fpc_N` (", num(fpc_N),
         "): you cannot have more sites than people in the population.")
  if (!is.null(fpc_N) && !is.null(n_per_site) && n_per_site > fpc_N)
    stop("`n_per_site` (", num(n_per_site), ") exceeds `fpc_N` (", num(fpc_N),
         "): a site cannot have more people than there are in the whole ",
         "population.")

  # A supplied icc below this tiny threshold (e.g. 1e-12, a rounding
  # leftover) is treated as exactly 0.
  icc_is_zero <- icc < sqrt(.Machine$double.eps)
  if (!icc_is_zero && is.null(n_sites) && is.null(n_per_site))
    stop("`icc` is greater than 0 but no sites were given (`n_sites` and ",
         "`n_per_site` are both `NULL`), so the design effect cannot be ",
         "computed. Supply `n_sites` (fix the number of sites) or `n_per_site` ",
         "(fix the number of people per site), or set `icc = 0` for a design ",
         "without clustering.")
  # The calculation with clustering assumes a very large population.
  if (!icc_is_zero && !is.null(fpc_N))
    stop("`fpc_N` cannot be used together with an `icc` above 0 (got icc = ",
         num(icc), "): the calculation with sites assumes a very large population. ",
         "Set `fpc_N = NULL`, or set `icc = 0`.")

  # Reverse mode: with n_sites fixed and icc above 0, every site needs at
  # least one person. (With n_per_site fixed, a number below one full site
  # is fine: forward mode can return one.)
  if (!solve_n && !is.null(n_sites) && !icc_is_zero && n < n_sites)
    stop("`n` (", num(n), ") is smaller than `n_sites` (", num(n_sites), "): there ",
         "would be fewer than one person per site. Increase `n`, or lower ",
         "`n_sites`.")
  # Testing the whole population (n = fpc_N) is allowed.
  if (!solve_n && !is.null(fpc_N) && n > fpc_N)
    stop("`n` (", num(n), ") is larger than `fpc_N` (", num(fpc_N), "): you cannot ",
         "test more people than there are in the population.")

  # From here on these are plain numbers: a name or a 1 x 1 matrix would
  # otherwise leak into the result (as in design_threshold()).
  prevalence  <- as.numeric(prevalence)
  sensitivity <- as.numeric(sensitivity)
  icc         <- as.numeric(icc)
  if (solve_n)  detection_prob <- as.numeric(detection_prob)
  if (!solve_n) n <- as.numeric(n)
  if (!is.null(n_sites))    n_sites    <- as.numeric(n_sites)
  if (!is.null(n_per_site)) n_per_site <- as.numeric(n_per_site)
  if (!is.null(fpc_N))      fpc_N      <- as.numeric(fpc_N)

  # ---- chance that one person is a case and tests positive ----
  q <- prevalence * sensitivity

  # ---- cases in a finite population ----
  # Rounded down, so the calculation never assumes more cases than the
  # prevalence gives. The small tolerance stops a product such as
  # 200 * 0.05 = 9.999999999999998 being rounded down to 9.
  if (!is.null(fpc_N)) {
    n_cases <- floor(fpc_N * prevalence + sqrt(.Machine$double.eps))
    if (n_cases < 1)
      stop("`fpc_N` (", num(fpc_N), ") x `prevalence` (", num(prevalence), ") is less ",
           "than 1: a population of ", num(fpc_N), " at this prevalence contains ",
           "no case to find. Use a larger population or a higher prevalence.")
  } else {
    n_cases <- NULL
  }

  # Number of people with no clustering (and a very large population) whose
  # log chance of missing every case is log_p0.
  n_equiv <- function(log_p0) log_p0 / log1p(-q)

  # ======================================================================
  # Reverse mode: n given -> chance of detection achieved.
  # ======================================================================
  if (!solve_n) {
    if (icc_is_zero) {
      deff  <- 1
      n_eff <- n
      # 1 - (1 - q)^n, worked out so that a tiny q is not lost to rounding
      detection <- if (is.null(fpc_N)) -expm1(n * log1p(-q))
                   else .detect_prob_finite(n, fpc_N, n_cases, sensitivity)
    } else {
      # The n people are spread as evenly as the sites allow: with n_sites
      # fixed, n %% n_sites sites get one extra person; with n_per_site
      # fixed, the last site holds the people left over.
      if (!is.null(n_sites)) {
        base   <- n %/% n_sites
        extra  <- n %% n_sites
        log_p0 <- extra * .log_p0_site(base + 1, q, icc) +
                  (n_sites - extra) * .log_p0_site(base, q, icc)
      } else {
        full   <- n %/% n_per_site
        rest   <- n %% n_per_site
        log_p0 <- full * .log_p0_site(n_per_site, q, icc) +
                  (if (rest > 0) .log_p0_site(rest, q, icc) else 0)
      }
      detection <- -expm1(log_p0)
      n_eff     <- n_equiv(log_p0)
      deff      <- n / n_eff
    }

    if (!is.null(n_per_site)) {
      n_sites_out    <- ceiling(n / n_per_site)
      n_per_site_out <- n_per_site
    } else if (!is.null(n_sites)) {
      n_sites_out    <- n_sites
      n_per_site_out <- ceiling(n / n_sites)
    } else {
      n_sites_out    <- NULL
      n_per_site_out <- NULL
    }

    return(list(
      n              = n,
      # The tiny allowance stops rounding noise adding a person.
      n_eff          = ceiling(n_eff * (1 - 1e-12)),
      n_sites        = n_sites_out,
      n_per_site     = n_per_site_out,
      detection_prob = detection,
      prevalence     = prevalence,
      q              = q,
      n_cases        = n_cases,
      sensitivity    = sensitivity,
      icc            = icc,
      deff           = deff,
      fpc_N          = fpc_N,
      mode           = "solve_detection_prob"
    ))
  }

  # ======================================================================
  # Forward mode: solve for n.
  # ======================================================================

  # One site larger than the number of people needed: it would be only
  # partly sampled (as in design_threshold()).
  stop_one_site <- function(needed)
    stop("`n_per_site` (", num(n_per_site), ") is larger than the number of ",
         "people needed (", num(needed), "), so one site would be more than ",
         "enough. Use a smaller `n_per_site`, or set `n_sites = 1` to find ",
         "how many people one site needs.")

  # More sites than people needed: one person at each site (no clustering
  # then) is already enough, so the design uses fewer sites, with 1 person
  # at each (as in design_threshold()).
  warn_fewer_sites <- function(needed)
    warning("With this prevalence and detection probability, the study ",
            "needs ", num(needed), if (needed == 1) " person" else " people",
            ", fewer than the ", num(n_sites), " sites given (`n_sites`). ",
            "The design uses ", num(needed),
            if (needed == 1) " site" else " sites", " with 1 person at each.",
            call. = FALSE)

  # No real study tests more than 2 million people, so larger designs stop
  # (the same limit as design_threshold()). `about` is the number needed,
  # when it is known (7 digits, so 2000002 does not print as 2,000,000).
  # Site advice is given only when the sites are what push the design past
  # the limit: with no clustering (`sites = FALSE`) it cannot help.
  stop_limit <- function(about = NULL, sites = TRUE) {
    fixes <- c("a higher prevalence", "a lower `detection_prob`")
    if (sites && !icc_is_zero && !is.null(n_per_site)) fixes <- c(fixes, "a smaller `n_per_site`")
    if (sites && !icc_is_zero && !is.null(n_sites))    fixes <- c(fixes, "more sites or a lower ICC")
    if (sensitivity < 1)                               fixes <- c(fixes, "a more sensitive test")
    stop("These settings need more than 2 million people",
         if (!is.null(about))
           paste0(" (about ", format(signif(ceiling(about), 7), big.mark = ",",
                                     scientific = FALSE), ")"),
         ". Use ", paste(fixes[-length(fixes)], collapse = ", "), ", or ",
         fixes[length(fixes)], ".")
  }

  # ---- number of people with no clustering ----
  # log1p(-x) is log(1 - x) without rounding 1 - x to 1 when x is tiny
  # (1 - 1e-17 is stored as exactly 1).
  n_inf <- log1p(-detection_prob) / log1p(-q)
  # With fpc_N the number needed can be smaller, so it is checked after the
  # finite-population search instead.
  if (is.null(fpc_N) && n_inf > 2e6) stop_limit(n_inf, sites = FALSE)
  # n_eff: people needed with no clustering and a very large population.
  n_eff <- ceiling(n_inf)

  # ---- with clustering ----
  if (!icc_is_zero) {
    # Log of the chance of missing every case that the target allows. The
    # tiny allowance stops rounding noise rejecting a design that reaches
    # the target exactly.
    log_miss <- log1p(-detection_prob) * (1 - 1e-12)

    if (!is.null(n_per_site)) {
      # Number of people per site fixed: each site misses every case with
      # the same chance, so k sites miss with that chance to the power k.
      log_p0_site <- .log_p0_site(n_per_site, q, icc)
      n_sites_out    <- ceiling(log_miss / log_p0_site)
      n_per_site_out <- n_per_site
      # One site is enough: find the fewest people that one site needs.
      if (n_sites_out == 1) {
        lo <- 0; hi <- n_per_site
        while (hi - lo > 1) {
          mid <- floor((lo + hi) / 2)
          if (.log_p0_site(mid, q, icc) <= log_miss) hi <- mid else lo <- mid
        }
        if (hi < n_per_site) stop_one_site(hi)
      }

    } else if (n_sites > n_eff) {
      # One person per site gives the chance of n_sites people with no
      # clustering, more than needed: use n_eff sites with 1 person at each.
      # Checked before the search, which would stop at 1 per site anyway.
      warn_fewer_sites(n_eff)
      n_sites_out    <- n_eff
      n_per_site_out <- 1

    } else {
      # Number of sites fixed: find the smallest number of people per site.
      # The chance of detection rises with every person added, so a
      # bisection search finds it. The search stops at 2 million people in
      # all (n_sites is at most n_eff here, so at least 1 per site).
      reaches <- function(m) n_sites * .log_p0_site(m, q, icc) <= log_miss
      m_max <- floor(2e6 / n_sites)
      if (!reaches(m_max)) {
        # With an ICC of 1 no number per site is enough; below 1 the target
        # is reached, but only with more than 2 million people.
        if (icc < 1) stop_limit()
        best <- -expm1(n_sites * .log_p0_site(m_max, q, icc))
        stop("A detection probability of ", .pct(detection_prob),
             " percent cannot be reached with ", num(n_sites),
             if (n_sites == 1) " site" else " sites", " and ICC = 1: with an ",
             "ICC of 1 everyone at a site has the same result, so these sites ",
             "can give at most ", .pct(best), " percent, however many people ",
             "each site has. Increase `n_sites`, lower the ICC, or lower ",
             "`detection_prob`.")
      }
      lo <- 0; hi <- m_max
      while (hi - lo > 1) {
        mid <- floor((lo + hi) / 2)
        if (reaches(mid)) hi <- mid else lo <- mid
      }
      n_sites_out    <- n_sites
      n_per_site_out <- hi
    }

    n_total <- n_sites_out * n_per_site_out
    if (n_total > 2e6) stop_limit(n_total)
    deff    <- n_total / n_equiv(n_sites_out *
                                 .log_p0_site(n_per_site_out, q, icc))

    return(list(
      n              = n_total,
      n_eff          = n_eff,
      n_sites        = n_sites_out,
      n_per_site     = n_per_site_out,
      detection_prob = detection_prob,
      prevalence     = prevalence,
      q              = q,
      n_cases        = n_cases,
      sensitivity    = sensitivity,
      icc            = icc,
      deff           = deff,
      fpc_N          = fpc_N,
      mode           = "solve_n"
    ))
  }

  # ---- no clustering ----
  if (is.null(fpc_N)) {
    n_total <- ceiling(n_inf)
  } else {
    # Even testing everyone misses all D cases with probability (1 - Se)^D.
    max_census <- .detect_prob_finite(fpc_N, fpc_N, n_cases, sensitivity)
    if (detection_prob > max_census)
      stop("A detection probability of ", round(100 * detection_prob, 1),
           " percent is not achievable in a population of ", num(fpc_N), " with ",
           num(n_cases), " cases and `sensitivity` = ", num(sensitivity), ": even ",
           "testing everyone finds a case with probability ",
           round(100 * max_census, 1), " percent. Lower `detection_prob`, or ",
           "use a more sensitive test.")
    # Smallest n that reaches the target. The chance of detection rises
    # with n, so a bisection search finds it.
    lo <- 0; hi <- fpc_N
    while (hi - lo > 1) {
      mid <- (lo + hi) %/% 2
      if (.detect_prob_finite(mid, fpc_N, n_cases, sensitivity) >= detection_prob)
        hi <- mid else lo <- mid
    }
    n_total <- hi
    if (n_total > 2e6) stop_limit(n_total)
  }

  # ---- distribute across sites (icc 0) ----
  if (!is.null(n_per_site) && n_per_site > n_total) stop_one_site(n_total)
  if (!is.null(n_per_site)) {
    n_sites_out    <- ceiling(n_total / n_per_site)
    n_per_site_out <- n_per_site
  } else if (!is.null(n_sites) && n_sites > n_total) {
    warn_fewer_sites(n_total)
    n_sites_out    <- n_total
    n_per_site_out <- 1
  } else if (!is.null(n_sites)) {
    n_per_site_out <- ceiling(n_total / n_sites)
    n_sites_out    <- n_sites
  } else {
    n_sites_out    <- NULL
    n_per_site_out <- NULL
  }

  # Rounding the number of sites or people per site up can push the design
  # past the population. That only happens close to a census.
  if (!is.null(fpc_N) && !is.null(n_sites_out) &&
      n_sites_out * n_per_site_out > fpc_N)
    warning("The rounded design (", num(n_sites_out), " sites x ", num(n_per_site_out),
            " = ", num(n_sites_out * n_per_site_out), ") exceeds `fpc_N` = ",
            num(fpc_N), ": the required n of ", num(n_total), " is close to a census of ",
            "the population. Consider surveying the whole population instead.")

  list(
    n              = n_total,
    n_eff          = n_eff,
    n_sites        = n_sites_out,
    n_per_site     = n_per_site_out,
    detection_prob = detection_prob,
    prevalence     = prevalence,
    q              = q,
    n_cases        = n_cases,
    sensitivity    = sensitivity,
    icc            = icc,
    deff           = 1,
    fpc_N          = fpc_N,
    mode           = "solve_n"
  )
}

# Log of the chance that no case is found among m people at one site, when
# the chance q varies between sites as a beta distribution with mean q and
# intra-cluster correlation icc (the beta-binomial model of DRpower's
# get_power_presence()). Needs icc above 0. With icc 1 everyone at a site
# has the same result, so the chance is 1 - q whatever m is.
.log_p0_site <- function(m, q, icc) {
  if (icc == 1) return(log1p(-q))
  a <- q * (1 / icc - 1)
  b <- (1 - q) * (1 / icc - 1)
  lbeta(a, b + m) - lbeta(a, b)
}

# A probability as a percentage for messages, with enough digits that a
# value just below a target does not print the same as the target.
.pct <- function(p) format(100 * p, digits = 7)

# Chance of finding at least one case when n people are drawn without
# replacement from a population of N that contains D cases, and each case
# in the sample tests positive with probability se. The number of cases k
# in the sample is hypergeometric; all k are missed with probability
# (1 - se)^k. Terms where (1 - se)^k is below 1e-20 are left out: they
# cannot change the result.
.detect_prob_finite <- function(n, N, D, se) {
  k_max <- min(n, D)
  if (se < 1) k_max <- min(k_max, ceiling(log(1e-20) / log1p(-se)))
  k <- 0:(if (se == 1) 0 else k_max)
  1 - sum(stats::dhyper(k, D, N - D, n) * (1 - se)^k)
}
