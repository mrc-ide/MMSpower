#' Calculate sample size to detect a prevalence above or below a threshold
#'
#' @description
#' Design function: given a decision threshold and the true prevalence
#' expected in the population, returns the number of people needed for
#' \code{test_threshold()} (with its default exact binomial test) to show
#' that the prevalence is above (or below) the threshold with the target
#' power.
#'
#' Three designs are supported:
#' \describe{
#'   \item{No clustering (\code{n_sites} and \code{n_per_site} both
#'     \code{NULL}, default)}{Returns the number of people.}
#'   \item{Fixed number of people per site (\code{n_per_site})}{Solves for
#'     the number of sites needed.}
#'   \item{Fixed number of sites (\code{n_sites})}{Solves for the number of
#'     people needed per site.}
#' }
#'
#' @param threshold Decision threshold for the true prevalence, strictly
#'   between 0 and 1 (for example 0.05 for the WHO pfhrp2/3 deletion
#'   threshold). This is the value \code{test_threshold()} tests against.
#' @param prevalence The true prevalence expected in the population, in
#'   (0, 1). It must be above \code{threshold} for
#'   \code{alternative = "greater"} and below it for \code{"less"}.
#' @param power Target power, in (0, 1); default 0.80.
#' @param alternative Which side of the threshold to look for, the same as
#'   in the \code{test_threshold()} call. \code{"greater"} (default) looks
#'   for a prevalence above the threshold, \code{"less"} for one below it,
#'   and \code{"two.sided"} for either.
#' @param sensitivity Diagnostic sensitivity in (0, 1]; default 1 (the test
#'   never misses a case).
#' @param specificity Diagnostic specificity in (0, 1]; default 1 (the test
#'   never gives a false positive). \code{sensitivity + specificity} must
#'   exceed 1.
#' @param conf_level Confidence level, in (0, 1); default 0.95. The
#'   significance level of the test is 1 - \code{conf_level}.
#' @param n_sites Optional whole number. Fixes the number of sites; the
#'   function solves for the number of people per site. Cannot be used with
#'   \code{n_per_site}.
#' @param n_per_site Optional whole number. Fixes the number of people per
#'   site; the function solves for the number of sites. Cannot be used with
#'   \code{n_sites}.
#' @param icc Intra-cluster correlation, in \[0, 1\]; default 0 (no
#'   clustering). An \code{icc} above 0 needs \code{n_sites} or
#'   \code{n_per_site}.
#'
#' @details
#' \strong{Exact power.} The power of a design is the chance that the exact
#' binomial test in \code{test_threshold()} rejects when the true prevalence
#' is \code{prevalence}. For \eqn{n} people and \code{alternative =
#' "greater"}, the test rejects when the number of positives \eqn{X} is at
#' least the critical count \eqn{k}, the smallest count with
#' \eqn{P(X \ge k \mid \theta_{app}) < \alpha}; the power is then
#' \eqn{P(X \ge k \mid p_{app})}. \code{"less"} is the mirror image, and
#' \code{"two.sided"} adds the two tails, each at \eqn{\alpha / 2}. This is
#' worked out from the binomial distribution itself, not from a normal
#' approximation, which can undersize the study badly at low prevalence.
#'
#' \strong{Saw-tooth power and the choice of n.} With an exact test, power
#' does not always go up when people are added: counts are whole numbers,
#' so the critical count jumps, and power goes up and down in a saw-tooth.
#' For example, with \code{threshold = 0.10}, \code{prevalence = 0.05},
#' \code{alternative = "less"} and 95 percent confidence, 179 people give
#' 81.3 percent power, 182 to 190 people all give less than 80 percent, and
#' every number from 203 on gives at least 80 percent. The function returns
#' the larger number (203): the smallest \eqn{n} from which every larger
#' \eqn{n} also reaches the target power. This is the safe choice, as a few
#' extra (or a few lost) people never drop the study below the target.
#' Sample sizes up to 100,000 are checked. A design needing more than that
#' without clustering (a prevalence very close to the threshold) gives a
#' warning that reports the power 100,000 people would reach, and returns
#' \code{NA} for every value.
#'
#' \strong{Imperfect test.} The test is worked out on the measured scale: a
#' true prevalence \eqn{p} is measured as
#' \eqn{p_{app} = p \cdot Se + (1-p)(1-Sp)}, and the same for the threshold.
#' The measured gap is the true gap times \eqn{Se + Sp - 1}, so an imperfect
#' test needs more people.
#'
#' \strong{Clustering.} The exact search above gives \eqn{n_0}, the number
#' of independent people needed (returned as \code{n_eff}). Clustering is
#' then allowed for with the design effect, as in
#' \code{design_precision()}: with \code{n_per_site} fixed,
#' \eqn{n = n_0 \cdot D_{eff}} with \eqn{D_{eff} = 1 + (n_{per site} - 1)
#' ICC}; with \code{n_sites} fixed,
#'
#' \deqn{n = \frac{n_0 \cdot n_{sites} \cdot (1-ICC)}{n_{sites} - n_0 \cdot ICC}}
#'
#' which exists only when \eqn{n_{sites} > n_0 \cdot ICC} (each site is
#' worth at most \eqn{1/ICC} independent people, however many are sampled
#' there); otherwise the function stops. The number of sites (or people per
#' site) is then rounded up, and \code{n} is \code{n_sites * n_per_site}.
#' Scaling the exact unclustered answer by the design effect is an
#' approximation: the effective sample size \eqn{n / D_{eff}} of the
#' returned design is at least \eqn{n_0}, but need not sit on the
#' saw-tooth exactly where \eqn{n_0} does. The power assumes the analysis
#' uses the same ICC.
#'
#' @return A list with (every value \code{NA}, with a warning, if more than
#'   100,000 people are needed without clustering):
#'   \item{n}{Total number of people needed. With sites,
#'     \code{n_sites * n_per_site}.}
#'   \item{n_sites}{Number of sites (\code{NULL} with no clustering).}
#'   \item{n_per_site}{People per site (\code{NULL} with no clustering).}
#'   \item{deff}{Design effect of the returned design,
#'     \code{1 + (n_per_site - 1) * icc} (1 with no clustering).}
#'   \item{n_eff}{Number of independent people needed: the exact answer
#'     with no clustering.}
#'
#' @references
#' MMS-SD Study Design Workshop. \url{https://mrc-ide.github.io/MMS-SD_workshop/}
#'
#' Kish, L. (1965) \emph{Survey Sampling}. Wiley. (Design effect.)
#'
#' Rogan, W. J. & Gladen, B. (1978) Estimating prevalence from the results
#' of a screening test. \emph{American Journal of Epidemiology}
#' \strong{107}(1), 71-76. \doi{10.1093/oxfordjournals.aje.a112510}
#'
#' @seealso \code{\link{test_threshold}}
#'
#' @export
#'
#' @examples
#' # Show that prevalence is above a 5 percent threshold when it is really
#' # 10 percent, with 80 percent power
#' design_threshold(threshold = 0.05, prevalence = 0.10)
#'
#' # Higher power needs more people
#' design_threshold(threshold = 0.05, prevalence = 0.10, power = 0.90)
#'
#' # Below the threshold instead
#' design_threshold(threshold = 0.10, prevalence = 0.05, alternative = "less")
#'
#' # Fixed number of people per site: how many sites?
#' design_threshold(threshold = 0.05, prevalence = 0.10,
#'                  n_per_site = 20, icc = 0.05)
design_threshold <- function(threshold,
                             prevalence,
                             power       = 0.80,
                             alternative = c("greater", "less", "two.sided"),
                             sensitivity = 1,
                             specificity = 1,
                             conf_level  = 0.95,
                             n_sites     = NULL,
                             n_per_site  = NULL,
                             icc         = 0) {

  .check_number(threshold, "threshold", 0, 1, lower_open = TRUE, upper_open = TRUE)
  .check_number(prevalence, "prevalence", 0, 1, lower_open = TRUE, upper_open = TRUE)
  .check_number(power, "power", 0, 1, lower_open = TRUE, upper_open = TRUE)
  .check_number(conf_level, "conf_level", 0, 1, lower_open = TRUE, upper_open = TRUE)
  .check_test(sensitivity, specificity)
  .check_sites(n_sites, n_per_site, icc)
  alternative <- match.arg(alternative)
  if (prevalence == threshold)
    stop("`prevalence` and `threshold` must differ.", call. = FALSE)
  if (alternative == "greater" && prevalence < threshold)
    stop("`prevalence` must be above `threshold` for alternative = 'greater'.",
         call. = FALSE)
  if (alternative == "less" && prevalence > threshold)
    stop("`prevalence` must be below `threshold` for alternative = 'less'.",
         call. = FALSE)

  # exact power for every unclustered sample size up to n_max
  n_max     <- 100000
  theta_app <- .apparent_prev(threshold, sensitivity, specificity)
  p_app     <- .apparent_prev(prevalence, sensitivity, specificity)
  pw <- .threshold_power(seq_len(n_max), theta_app, p_app, alternative,
                         1 - conf_level)

  # n0: the smallest n from which every larger n reaches the target power
  below <- which(pw < power)
  if (length(below) > 0 && max(below) == n_max) {
    warning(sprintf(paste0("These settings need more than %s people without ",
                           "clustering (%s people give %.1f%% power), so no ",
                           "design is returned. Use a prevalence further from ",
                           "the threshold, a lower power, or a more accurate test."),
                    formatC(n_max, format = "d", big.mark = ","),
                    formatC(n_max, format = "d", big.mark = ","),
                    100 * pw[n_max]), call. = FALSE)
    return(list(n = NA_real_, n_sites = NA_real_, n_per_site = NA_real_,
                deff = NA_real_, n_eff = NA_real_))
  }
  n0 <- if (length(below) == 0) 1 else max(below) + 1

  # clustered sample size, as in design_precision()
  if (!is.null(n_per_site)) {
    n <- n0 * .deff(n_per_site, icc)
  } else if (!is.null(n_sites)) {
    n <- .solve_fixed_sites(n0, n_sites, icc)
    if (is.na(n))
      stop("Target power is unreachable with ", n_sites, " sites and ICC = ",
           icc, ": the sites are worth at most ", floor(n_sites / icc),
           " independent people, but ", n0, " are needed. Increase `n_sites` ",
           "or lower the target power.", call. = FALSE)
  } else {
    n <- n0
  }

  design <- .allocate(n, n_sites, n_per_site)

  list(
    n          = design$n,
    n_sites    = design$n_sites,
    n_per_site = design$n_per_site,
    deff       = if (is.null(design$n_per_site)) 1 else .deff(design$n_per_site, icc),
    n_eff      = n0
  )
}

# Exact power of the binomial test in test_threshold() for each sample size
# in `n` (a vector), testing against measured threshold `theta` when the
# measured prevalence is really `p`. The test rejects for counts whose
# one-sided p-value is below `a` (alpha, or alpha / 2 for each tail of a
# two-sided test); power is the chance of such a count.
.threshold_power <- function(n, theta, p, alternative, alpha) {
  # P(X >= k) and P(X <= k) for X ~ Binomial(n, q)
  up <- function(k, q) stats::pbinom(k - 1, n, q, lower.tail = FALSE)
  dn <- function(k, q) stats::pbinom(k, n, q)

  # Smallest k with P(X >= k | theta) < a. qbinom() gives it to within one,
  # and the two steps after it make it exact.
  upper_crit <- function(a) {
    k <- stats::qbinom(a, n, theta, lower.tail = FALSE) + 1
    k <- k + (up(k, theta) >= a)
    k - (up(k - 1, theta) < a)
  }
  # Largest k with P(X <= k | theta) < a (-1 if there is none).
  lower_crit <- function(a) {
    k <- stats::qbinom(a, n, theta) - 1
    k <- k - (dn(k, theta) >= a)
    k + (dn(k + 1, theta) < a)
  }

  switch(alternative,
    greater   = up(upper_crit(alpha), p),
    less      = dn(lower_crit(alpha), p),
    two.sided = up(upper_crit(alpha / 2), p) + dn(lower_crit(alpha / 2), p))
}
