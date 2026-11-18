#' Calculate sample size to detect a prevalence above or below a threshold
#'
#' @description
#' Design function: given a decision threshold and the true prevalence
#' expected in the population, returns the number of people needed for
#' \code{test_threshold()} to show that the prevalence is above (or below)
#' the threshold with the target power.
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
#' @param power Target power, between 1 - \code{conf_level} and 1; default
#'   0.80. Values below 0.5 give a warning.
#' @param alternative Which side of the threshold to look for, the same as
#'   in the \code{test_threshold()} call. \code{"greater"} (default) looks
#'   for a prevalence above the threshold, \code{"less"} for one below it,
#'   and \code{"two.sided"} for either.
#' @param sensitivity Diagnostic sensitivity in (0, 1]: the share of people
#'   who truly have the condition that the test detects as positive.
#'   Default 1 (the test never misses a case).
#' @param specificity Diagnostic specificity in (0, 1]: the share of people
#'   who truly do not have the condition that the test correctly reads as
#'   negative. Default 1 (the test never gives a false positive).
#' @param conf_level Confidence level, in \[0.5, 1); default 0.95. The
#'   significance level of the test is 1 - \code{conf_level}.
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
#'   \code{n_per_site}, and the value given is then used exactly, however
#'   small. Without sites, a tiny value below about 1.5e-8, such as one left
#'   over from rounding in an earlier calculation, is treated as 0.
#' @param fpc_N Optional. Population size for a finite-population
#'   correction: a positive integer. \code{NULL} (default) means no
#'   finite-population correction.
#'
#' @details
#' \strong{Exact power.} The power of a design is the chance that
#' \code{test_threshold()}, which uses an exact binomial test, finds the
#' difference when the true prevalence is \code{prevalence}. It is worked out
#' from the binomial distribution itself, for every possible number of
#' positives, not from a normal approximation.
#'
#' \strong{Why two numbers.} With an exact test, power does not always go up
#' when people are added: the counts are whole numbers, so the cut-off for a
#' significant result jumps, and power goes up and down in a saw-tooth.
#' For example, with \code{threshold = 0.10}, \code{prevalence = 0.05},
#' \code{alternative = "less"} and 95 percent confidence, 179 people give
#' 81.3 percent power, 182 to 190 people all give less than 80 percent, and
#' every number from 203 on gives at least 80 percent. \code{n} is the safe
#' number (203): the smallest number from which every larger number also
#' reaches the target power, so a few extra people never drop the study
#' below it. \code{n_first} is the smallest number that reaches it (179).
#'
#' \strong{Sites and population size.} These enter as in
#' \code{test_threshold()}: the effective size is \eqn{n / D_{eff}}, with the
#' design effect \eqn{D_{eff} = 1 + (n_{per site} - 1) ICC}, and with
#' \code{fpc_N} (the population size \eqn{N}) it is divided further by
#' \eqn{(N - n)/(N - 1)}. With \code{n_per_site} fixed the search is over
#' the number of sites; with \code{n_sites} fixed it is over the number of
#' people per site. The power assumes the analysis uses the same ICC.
#' With sites the power is approximate: in simulations of 12 designs with a
#' target of 0.80, \code{test_threshold()} reached 0.78 to 0.82 with the ICC
#' supplied, and as low as 0.74 when the ICC was estimated from the data
#' (\code{icc = NULL}).
#'
#' \strong{Imperfect test.} The test is worked out on the measured scale: a
#' true prevalence \eqn{p} is measured as
#' \eqn{p_{app} = p \cdot Se + (1-p)(1-Sp)}, and the same for the threshold.
#' The measured gap is the true gap times \eqn{Se + Sp - 1}, so an imperfect
#' test needs more people.
#'
#' @return A named list with the following fields, in the order returned:
#'   \item{n}{Total number of people needed: the smallest number from which
#'     every larger number also reaches the target power (see Details). With
#'     sites, \code{n_sites * n_per_site}.}
#'   \item{n_first}{The smallest total number of people that reaches the
#'     target power. Some numbers between \code{n_first} and \code{n} fall
#'     below it.}
#'   \item{n_eff}{The number of people needed with no sites and no
#'     \code{fpc_N}, worked out the same way as \code{n}. Clustering raises
#'     \code{n} above it and \code{fpc_N} can lower \code{n} below it. With
#'     sites, only whole sites are checked, so with a small ICC \code{n} can
#'     also come out a little below it.}
#'   \item{n_sites}{Number of sites: solved for when \code{n_per_site} is
#'     given, echoed when \code{n_sites} is given. If the study needs fewer
#'     people than there are sites, fewer sites are used, with 1 person at
#'     each, and a warning says so. \code{NULL} with no clustering.}
#'   \item{n_per_site}{People per site: solved for when \code{n_sites} is
#'     given, echoed when \code{n_per_site} is given. \code{NULL} with no
#'     clustering.}
#'   \item{threshold}{The threshold, as supplied.}
#'   \item{prevalence}{The true prevalence, as supplied.}
#'   \item{apparent_prev}{The prevalence the test would measure.}
#'   \item{threshold_app}{The threshold on the measured scale.}
#'   \item{power}{Target power, as supplied.}
#'   \item{alternative}{Which side of the threshold was looked for, as
#'     supplied.}
#'   \item{conf_level}{Confidence level, as supplied.}
#'   \item{sensitivity}{Sensitivity, as supplied.}
#'   \item{specificity}{Specificity, as supplied.}
#'   \item{icc}{ICC, as supplied.}
#'   \item{deff}{Design effect of the returned design,
#'     \code{1 + (n_per_site - 1) * icc} (1 with no clustering).}
#'   \item{fpc_N}{\code{fpc_N}, as supplied, or \code{NULL}.}
#'
#' @references
#' MMS-SD Study Design Workshop. \url{https://mrc-ide.github.io/MMS-SD_workshop/}
#'
#' Clopper, C. J. & Pearson, E. S. (1934) The use of confidence or fiducial
#' limits illustrated in the case of the binomial. \emph{Biometrika}
#' \strong{26}(4), 404-413. \doi{10.1093/biomet/26.4.404}
#'
#' Kish, L. (1965) \emph{Survey Sampling}. Wiley. (Design effect.)
#'
#' Cochran, W. G. (1977) \emph{Sampling Techniques}, 3rd ed. Wiley.
#' (Finite-population correction.)
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
                              alternative = "greater",
                              sensitivity = 1,
                              specificity = 1,
                              conf_level  = 0.95,
                              n_sites     = NULL,
                              n_per_site  = NULL,
                              icc         = 0,
                              fpc_N       = NULL) {

  # ---- validate single-value parameters ----
  # Each of the parameters below must be a single value.
  if (length(threshold) != 1 || !is.numeric(threshold))
    stop("`threshold` must be a single number in (0, 1) (got class `",
         class(threshold)[1], "`, length ", length(threshold), ").")
  if (length(prevalence) != 1 || !is.numeric(prevalence))
    stop("`prevalence` must be a single number in (0, 1) (got class `",
         class(prevalence)[1], "`, length ", length(prevalence), ").")
  if (length(power) != 1 || !is.numeric(power))
    stop("`power` must be a single number in (0, 1) (got class `",
         class(power)[1], "`, length ", length(power), ").")
  if (length(sensitivity) != 1 || !is.numeric(sensitivity))
    stop("`sensitivity` must be a single number in (0, 1] (got class `",
         class(sensitivity)[1], "`, length ", length(sensitivity), "). ",
         "Note: `TRUE`, `FALSE` and a plain `NA` are logical, not numeric. ",
         "Pass 1 if the test never misses a case.")
  if (length(specificity) != 1 || !is.numeric(specificity))
    stop("`specificity` must be a single number in (0, 1] (got class `",
         class(specificity)[1], "`, length ", length(specificity), "). ",
         "Note: `TRUE`, `FALSE` and a plain `NA` are logical, not numeric. ",
         "Pass 1 if the test never gives a false positive.")
  if (length(conf_level) != 1 || !is.numeric(conf_level))
    stop("`conf_level` must be a single number in [0.5, 1) (got class `",
         class(conf_level)[1], "`, length ", length(conf_level), ").")
  if (length(icc) != 1 || !is.numeric(icc))
    stop("`icc` must be a single number in [0, 1] (got class `",
         class(icc)[1], "`, length ", length(icc), ").")
  if (!is.character(alternative) || length(alternative) != 1)
    stop("`alternative` must be one of 'greater', 'less', or 'two.sided' ",
         "(got class `", class(alternative)[1], "`, length ",
         length(alternative), ").")
  if (!alternative %in% c("greater", "less", "two.sided"))
    stop("`alternative` must be one of 'greater', 'less', or 'two.sided' ",
         "(got '", alternative, "').")

  # Check for NA, NaN and Inf before the range checks below.
  if (!is.finite(threshold))
    stop("`threshold` must be a single finite number (got ", threshold, ").")
  if (!is.finite(prevalence))
    stop("`prevalence` must be a single finite number (got ", prevalence, ").")
  if (!is.finite(power))
    stop("`power` must be a single finite number (got ", power, ").")
  if (!is.finite(sensitivity))
    stop("`sensitivity` must be a single finite number (got ", sensitivity, ").")
  if (!is.finite(specificity))
    stop("`specificity` must be a single finite number (got ", specificity, ").")
  if (!is.finite(conf_level))
    stop("`conf_level` must be a single finite number (got ", conf_level, ").")
  if (!is.finite(icc))
    stop("`icc` must be a single finite number (got ", icc, ").")

  # Numbers in messages are printed in full (100000, not 1e+05).
  num <- function(v) format(v, scientific = FALSE, digits = 15)
  # ...except numbers far beyond any real study, shown as 1e+200.
  huge <- function(v) if (v < 1e16) num(v) else format(v, digits = 15)
  # A value within rounding of a whole number (such as 100 * 1.1, stored as
  # 110.00000000000001) prints as that whole number, so say why it fails.
  near_int <- function(v, nm)
    if (is.finite(v) && v != round(v) &&
        abs(v - round(v)) <= 8 * .Machine$double.eps * abs(v))
      paste0(", which is not exactly an integer; use round(", nm, ")") else ""

  # `threshold` as in test_threshold().
  if (threshold < 0 || threshold > 1)
    stop("`threshold` must be in (0, 1) (got ", num(threshold), "). ",
         "`threshold` is a prevalence and cannot be ",
         if (threshold < 0) "negative." else "greater than 1.")
  if (threshold == 0 || threshold == 1)
    stop("`threshold` must be in (0, 1) (got ", num(threshold), "). ",
         "No prevalence can be below 0 or above 1, so there is nothing to ",
         "test. Use the prevalence that would lead to action, for example 0.05.")
  if (prevalence < 0 || prevalence > 1)
    stop("`prevalence` must be in (0, 1) (got ", num(prevalence), "). ",
         "`prevalence` is a fraction of the population and cannot be ",
         if (prevalence < 0) "negative." else "greater than 1.")
  if (prevalence == 0 || prevalence == 1)
    stop("`prevalence` must be in (0, 1) (got ", num(prevalence), "). ",
         "Use a value strictly between 0 and 1.")
  if (power >= 1)
    stop("`power` must be in (0, 1) (got ", num(power), "). ",
         "For example, use 0.80 for 80 percent power.")
  if (sensitivity < 0 || sensitivity > 1)
    stop("`sensitivity` must be in (0, 1] (got ", num(sensitivity), "). ",
         "`sensitivity` is a diagnostic probability and cannot be ",
         if (sensitivity < 0) "negative." else "greater than 1.")
  if (sensitivity == 0)
    stop("`sensitivity` must be in (0, 1] (got 0). ",
         "A sensitivity of 0 means everyone who truly has the condition ",
         "tests negative.")
  if (specificity < 0 || specificity > 1)
    stop("`specificity` must be in (0, 1] (got ", num(specificity), "). ",
         "`specificity` is a diagnostic probability and cannot be ",
         if (specificity < 0) "negative." else "greater than 1.")
  if (specificity == 0)
    stop("`specificity` must be in (0, 1] (got 0). ",
         "A specificity of 0 means everyone who truly does not have the ",
         "condition tests positive.")
  # A sum above 1 only by rounding (within 1e-12 of 1, e.g. 0.3 +
  # 0.7000000000000004) counts as 1.
  correction <- sensitivity + specificity - 1
  if (correction < 1e-12)
    stop("`sensitivity` + `specificity` must exceed 1 ",
         "(got ", num(sensitivity), " + ", num(specificity), " = ",
         # 12 digits, so a sum that is 1 up to rounding (within 1e-12) reads 1
         format(sensitivity + specificity, digits = 12, scientific = FALSE), "). ",
         "The difference the test can see is the true difference times ",
         "`sensitivity` + `specificity` - 1, so with a sum of 1 or less the ",
         "test cannot see the difference, or sees it reversed.")
  if (icc < 0 || icc > 1)
    stop("`icc` must be in [0, 1] (got ", num(icc), "). ",
         "`icc` is a correlation and cannot be ",
         if (icc < 0) "negative." else "greater than 1.")
  if (conf_level >= 1)
    stop("`conf_level` must be in [0.5, 1) (got ", num(conf_level), "). ",
         "For example, use 0.95 for a 95% confidence interval.")
  if (conf_level < 0.5)
    stop("`conf_level` must be at least 0.5 (got ", num(conf_level), "). ",
         "Below 0.5, when the prevalence is exactly at the threshold, the ",
         "test can wrongly find it above or below the threshold more than ",
         "half the time.")
  # Power can never be below the false-alarm rate: even with the prevalence
  # exactly at the threshold the test reports a difference that often. A
  # power equal to it only up to rounding (within 1e-12, e.g. 1 - 0.8 =
  # 0.19999999999999996) counts as equal.
  false_alarm <- 1 - conf_level
  if (power <= false_alarm + 1e-12)
    stop("`power` must be greater than 1 - `conf_level` = ",
         signif(false_alarm, 6), " (got ", num(power), "). Even when the ",
         "prevalence is exactly at the threshold, the test reports a difference ",
         signif(100 * false_alarm, 6), " percent of the time, so a lower power ",
         "is not meaningful.")
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
         "). `fpc_N` is the total population size. ",
         "Set `fpc_N = NULL` to skip the finite-population correction.")
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

  # From here on every number is a plain number: a 1 x 1 matrix, a name or a
  # class would otherwise leak into the arithmetic and the result.
  threshold   <- as.numeric(threshold);   prevalence  <- as.numeric(prevalence)
  power       <- as.numeric(power);       conf_level  <- as.numeric(conf_level)
  sensitivity <- as.numeric(sensitivity); specificity <- as.numeric(specificity)
  icc         <- as.numeric(icc)
  if (!is.null(n_sites))    n_sites    <- as.numeric(n_sites)
  if (!is.null(n_per_site)) n_per_site <- as.numeric(n_per_site)
  if (!is.null(fpc_N))      fpc_N      <- as.numeric(fpc_N)

  # The number of sites, or people per site, cannot exceed the population.
  if (!is.null(fpc_N) && !is.null(n_sites) && n_sites > fpc_N)
    stop("`n_sites` (", num(n_sites), ") exceeds `fpc_N` (", num(fpc_N),
         "): you cannot have more sites than people in the population.")
  if (!is.null(fpc_N) && !is.null(n_per_site) && n_per_site > fpc_N)
    stop("`n_per_site` (", num(n_per_site), ") exceeds `fpc_N` (", num(fpc_N),
         "): a site cannot have more people than there are in the whole ",
         "population.")

  # Without sites there is no clustering, so a tiny icc there (e.g. 1e-12, a
  # rounding leftover) counts as 0. With sites, any icc above 0 is used: with
  # very large sites even a tiny icc changes the design effect.
  no_sites <- is.null(n_sites) && is.null(n_per_site)
  if (no_sites && icc >= sqrt(.Machine$double.eps))
    stop("`icc` is greater than 0 but no sites were given (`n_sites` and ",
         "`n_per_site` are both `NULL`), so the design effect cannot be ",
         "computed. Supply `n_sites` (fix the number of sites) or `n_per_site` ",
         "(fix the number of people per site), or set `icc = 0` for a design ",
         "without clustering.")
  icc_is_zero <- no_sites || icc == 0

  # ---- direction of the test ----
  # Checked after "equal", so that the advice to switch the alternative is
  # only given when switching would work.
  if (abs(prevalence - threshold) < 1e-8)
    stop("`prevalence` (", num(prevalence), ") and `threshold` (", num(threshold),
         ") are equal, or nearly so, so there is no difference to detect. ",
         "Use a prevalence further from the threshold.")
  if (alternative == "greater" && prevalence < threshold)
    stop("`prevalence` (", num(prevalence), ") must be greater than `threshold` (",
         num(threshold), ") when alternative = 'greater'. Use alternative = 'less'.")
  if (alternative == "less" && prevalence > threshold)
    stop("`prevalence` (", num(prevalence), ") must be less than `threshold` (",
         num(threshold), ") when alternative = 'less'. Use alternative = 'greater'.")

  # ---- measured (apparent) prevalences ----
  theta_app <- .apparent_prev(threshold,  sensitivity, specificity)
  p1_app    <- .apparent_prev(prevalence, sensitivity, specificity)

  # An imperfect test can bring two different prevalences very close
  # together on the measured scale.
  delta <- abs(p1_app - theta_app)
  if (delta < 1e-8) {
    # Enough digits that the two numbers never print the same.
    digits <- 6
    repeat {
      shown <- c(format(p1_app, digits = digits, scientific = FALSE),
                 format(theta_app, digits = digits, scientific = FALSE))
      if (shown[1] != shown[2] || digits >= 15) break
      digits <- digits + 1
    }
    stop("The measured prevalence and threshold differ by only ",
         format(signif(delta, 3), scientific = FALSE),
         " (", shown[1], " and ", shown[2], "), so the required sample size ",
         "would be extremely large. Use a prevalence further from the ",
         "threshold, or a more accurate test.")
  }

  # ---- warnings, once every check above has passed ----
  # A sum of 1.1 that comes out a hair below it by rounding (0.15 + 0.95)
  # does not warn.
  if (correction < 0.1 - 1e-12)
    warning("`sensitivity` + `specificity` = ", num(sensitivity + specificity),
            ", which is below 1.1. The test then barely tells apart ",
            "people with and without the condition, so the required n is very ",
            "large, and small changes in the assumed sensitivity or ",
            "specificity will change n a lot.")
  if (power < 0.5)
    warning("`power` = ", num(power), " is below 0.5, so the study would be more ",
            "likely to miss a real difference than to find it.")

  # Advice when the target cannot be reached.
  advice <- if (!is.null(n_sites))
    "Increase `n_sites`, lower the ICC, or lower the target power." else
    "Use a prevalence further from the threshold, or lower the target power."

  # ---- starting point: the normal-approximation formula ----
  # The exact power is searched for below; this formula only says roughly
  # how far to look.
  # 1 - 0.95 is stored as 0.050000000000000044; rounded to 15 digits it is
  # exactly 0.05, so an exact p-value of 0.05 is not counted as below it.
  alpha  <- signif(1 - conf_level, 15)
  z_a    <- if (alternative == "two.sided") stats::qnorm(1 - alpha / 2) else
                                            stats::qnorm(1 - alpha)
  z_b    <- stats::qnorm(power)

  se0    <- sqrt(theta_app * (1 - theta_app))
  se1    <- sqrt(p1_app    * (1 - p1_app))

  n_base_cont <- ((z_a * se0 + z_b * se1) / delta)^2

  # ---- exact search ----
  # Without sites: the number of people, n = u.
  # With `n_per_site` fixed: the number of sites, n = u * n_per_site.
  # With `n_sites` fixed: the number per site, n = n_sites * u.
  # (An icc above 0 with no sites was stopped above.)
  if (is.null(n_sites) && is.null(n_per_site)) {
    unit_n <- function(u) u
    unit_d <- function(u) rep(1, length(u))
    u0     <- n_base_cont
  } else if (!is.null(n_per_site)) {
    d      <- if (icc_is_zero) 1 else 1 + (n_per_site - 1) * icc
    unit_n <- function(u) u * n_per_site
    unit_d <- function(u) rep(d, length(u))
    u0     <- n_base_cont * d / n_per_site
  } else {
    unit_n <- function(u) u * n_sites
    unit_d <- if (icc_is_zero) function(u) rep(1, length(u)) else
      function(u) 1 + (u - 1) * icc
    # With an ICC, the effective size n / D_eff can never pass
    # n_sites / icc, however many people each site has. When the formula says
    # that is not enough, it gives no answer (u0 = Inf), but the exact test
    # can still reach the target, so the search below still looks.
    u0     <- if (icc_is_zero) n_base_cont / n_sites else
      if (n_sites - n_base_cont * icc <= 0) Inf else
      n_base_cont * (1 - icc) / (n_sites - n_base_cont * icc)
  }
  sites_limit <- !icc_is_zero && !is.null(n_sites)
  # The exact search checks every possible design one by one, so it is
  # limited to designs of up to about 2 million people.
  if (is.finite(u0) && unit_n(u0) > 2e6)
    stop("These settings need more than 2 million people (about ",
         format(signif(unit_n(u0), 3), big.mark = ",", scientific = FALSE),
         "), too many to work out exactly. Use a prevalence further from the ",
         "threshold, or a more accurate test.")
  # More sites than people needed: one person at each site (no clustering
  # then) is already enough, so the design uses fewer sites, as in
  # design_difference(). Checked before the search, which would otherwise
  # start with one person at each of the n_sites sites.
  fewer_sites <- FALSE
  if (!is.null(n_sites)) {
    one_each <- .thr_search(function(u) u, function(u) rep(1, length(u)),
                            n_base_cont, if (is.null(fpc_N)) Inf else fpc_N - 1,
                            fpc_N, theta_app, p1_app, alternative, alpha, power)
    if (!is.null(one_each) && one_each$safe < n_sites) {
      fewer_sites <- TRUE
      warning("With this prevalence, threshold and power, the study needs ",
              num(one_each$safe), if (one_each$safe == 1) " person" else " people",
              ", fewer than the ", num(n_sites), " sites given (`n_sites`). ",
              "The design uses ", num(one_each$safe),
              if (one_each$safe == 1) " site" else " sites", " with 1 person at each.")
    }
  }
  # With fpc_N the test needs at least one person left untested.
  u_max <- if (is.null(fpc_N)) Inf else floor((fpc_N - 1) / unit_n(1))
  if (!fewer_sites && u_max < 1)
    stop("`fpc_N` (", num(fpc_N), ") leaves no room for a design: the test ",
         "needs at least one person in the population who is not tested.")
  # When the formula gives no answer, look up to the number per site at which
  # the effective size is 99.9 percent of its limit n_sites / icc (more people
  # per site add almost nothing), and no further than 2 million people.
  if (!is.finite(u0))
    u_max <- min(u_max, floor(2e6 / n_sites),
                 max(1, ceiling(999 * (1 - icc) / icc)))

  if (fewer_sites) {
    found  <- one_each
    unit_n <- function(u) u
  } else {
    found <- .thr_search(unit_n, unit_d, u0, u_max, fpc_N, theta_app, p1_app,
                         alternative, alpha, power)
  }
  if (is.null(found)) {
    if (is.null(fpc_N) && sites_limit)
      stop("Target power of ", num(100 * power), " percent cannot be reached ",
           "with ", num(n_sites), if (n_sites == 1) " site" else " sites",
           " and ICC = ", num(icc), ". ", advice)
    stop("Target power of ", num(100 * power), " percent cannot be reached",
         if (is.null(fpc_N)) "" else " within the population (`fpc_N`)",
         " with these settings. ", advice)
  }
  u_safe  <- found$safe
  n_total <- unit_n(u_safe)
  n_first <- unit_n(found$first)

  # One site larger than the number of people needed: it would be only partly
  # sampled, so its size would not be the real site size in the design
  # effect (as in design_difference()).
  if (!is.null(n_per_site) && u_safe == 1) {
    needed <- .thr_search(function(u) u, function(u) rep(d, length(u)),
                          n_base_cont * d, if (is.null(fpc_N)) Inf else fpc_N - 1,
                          fpc_N, theta_app, p1_app, alternative, alpha, power)$safe
    if (n_per_site > needed)
      stop("`n_per_site` (", num(n_per_site), ") is larger than the number of ",
           "people needed (", num(needed), "), so one site would be more than ",
           "enough. Use a smaller `n_per_site`, or set `n_sites = 1` to find ",
           "how many people one site needs.")
  }

  # n_eff: the number of people needed with no sites and no fpc_N (same
  # exact search). Clustering raises `n` above it; an FPC, or checking only
  # whole sites, can bring `n` below it.
  n_eff <- .thr_search(function(u) u, function(u) rep(1, length(u)),
                       n_base_cont, Inf, NULL, theta_app, p1_app,
                       alternative, alpha, power)$safe

  # ---- sites ----
  if (!is.null(n_per_site)) {
    n_sites_out    <- u_safe
    n_per_site_out <- n_per_site
    deff           <- unit_d(u_safe)
  } else if (fewer_sites) {
    n_sites_out    <- u_safe
    n_per_site_out <- 1
    deff           <- 1
  } else if (!is.null(n_sites)) {
    n_sites_out    <- n_sites
    n_per_site_out <- u_safe
    deff           <- unit_d(u_safe)
  } else {
    n_sites_out    <- NULL
    n_per_site_out <- NULL
    deff           <- 1
  }

  list(
    n             = n_total,
    n_first       = n_first,
    n_eff         = n_eff,
    n_sites       = n_sites_out,
    n_per_site    = n_per_site_out,
    threshold     = threshold,
    prevalence    = prevalence,
    apparent_prev = p1_app,
    threshold_app = theta_app,
    power         = power,
    alternative   = alternative,
    conf_level    = conf_level,
    sensitivity   = sensitivity,
    specificity   = specificity,
    icc           = icc,
    deff          = deff,
    fpc_N         = fpc_N
  )
}

# Exact power of test_threshold() for designs of total size n (vector) with
# design effect d (vector), as test_threshold() would run it with the same
# icc and fpc_N: effective size m = n / d (times the FPC factor), each
# possible count X = 0..n of positives counts as X * m / n, and the test
# rejects for the counts whose exact p-value is below alpha. Power is the
# chance of such a count when the measured prevalence is p, worked out on
# the effective size in the same way.
.thr_power <- function(n, d, fpc_N, theta, p, alternative, alpha) {
  m <- n / d
  if (!is.null(fpc_N)) m <- m * (fpc_N - 1) / (fpc_N - n)
  s <- m / n
  # Smallest whole X whose upper-tail p-value is below a (or n + 1 if none):
  # the p-value goes down as X goes up, so halve the range each step.
  first_upper <- function(a) {
    lo <- rep(0, length(n)); hi <- n + 1
    while (any(hi - lo > 1)) {
      mid <- floor((lo + hi) / 2)
      ok  <- .binom_upper(mid * s, m, theta) < a
      hi  <- ifelse(ok, mid, hi); lo <- ifelse(ok, lo, mid)
    }
    ifelse(.binom_upper(0, m, theta) < a, 0, hi)
  }
  # Largest whole X whose lower-tail p-value is below a (or -1 if none).
  last_lower <- function(a) {
    lo <- rep(-1, length(n)); hi <- n
    while (any(hi - lo > 1)) {
      mid <- ceiling((lo + hi) / 2)
      ok  <- .binom_lower(mid * s, m, theta) < a
      lo  <- ifelse(ok, mid, lo); hi <- ifelse(ok, hi, mid)
    }
    ifelse(.binom_lower(n * s, m, theta) < a, n, lo)
  }
  up <- function(a) { k <- first_upper(a)
    ifelse(k > n, 0, .binom_upper(k * s, m, p)) }
  dn <- function(a) { k <- last_lower(a)
    ifelse(k < 0, 0, .binom_lower(k * s, m, p)) }
  switch(alternative,
    greater   = up(alpha),
    less      = dn(alpha),
    two.sided = up(alpha / 2) + dn(alpha / 2))
}

# Search the design unit u (people, sites, or people per site) for the
# first u that reaches the target power, and the smallest u from which every
# larger u also reaches it. Exact power goes up and down as u grows, so both
# are reported. "Every larger u" is checked over a long stretch: the search
# stops only once the target has held for a quarter of the way again past
# the last failure (and at least 200 more steps). The search starts at 1, so
# the first u is never missed. NULL if the target is never reached by u_max
# (or by 4 times the normal-approximation answer u0).
.thr_search <- function(unit_n, unit_d, u0, u_max, fpc_N, theta, p,
                        alternative, alpha, power) {
  pw_at <- function(u)
    .thr_power(unit_n(u), unit_d(u), fpc_N, theta, p, alternative, alpha)
  # Look no further than 4 times the normal-approximation answer.
  u_max <- min(u_max, max(ceiling(4 * u0), ceiling(u0) + 1000))
  lo <- 1
  first <- NA; last_bad <- 0
  chunk <- 5000
  repeat {
    hi <- min(u_max, lo + chunk - 1)
    u  <- lo:hi
    ok <- pw_at(u) >= power
    if (is.na(first) && any(ok)) first <- u[which(ok)[1]]
    if (any(!ok)) last_bad <- u[max(which(!ok))]
    if (!is.na(first) && hi - last_bad >= max(200, 0.25 * last_bad))
      return(list(first = first, safe = last_bad + 1))
    if (hi == u_max)
      return(if (is.na(first) || last_bad == u_max) NULL else
        list(first = first, safe = last_bad + 1))
    lo <- hi + 1
  }
}
