#' Calculate sample size to detect a difference in prevalence between two groups
#'
#' @description
#' Design function: given the expected true prevalence in two groups, returns
#' the sample size needed \strong{per group} to detect that difference with
#' the target power.
#'
#' The two groups get the same number of people. Clustering and the
#' finite-population correction are applied to each group in the same way.
#' Three designs are supported:
#' \describe{
#'   \item{No clustering (\code{n_sites} and \code{n_per_site} both
#'     \code{NULL}, the default)}{Returns the number of people per group.}
#'   \item{Fixed number of people per site (\code{n_per_site})}{Solves for
#'     the number of sites needed per group.}
#'   \item{Fixed number of sites (\code{n_sites})}{Solves for the number of
#'     people needed per site.}
#' }
#'
#' @param prevalence1 Expected true prevalence in group 1, in \[0, 1\].
#' @param prevalence2 Expected true prevalence in group 2, in \[0, 1\].
#'   The two prevalences must differ; either can be 0 or 1.
#' @param power Target power, between 1 - \code{conf_level} and 1; default
#'   0.80. Values below 0.5 give a warning.
#' @param alternative Which kind of difference to look for.
#'   \code{"two.sided"} (default) looks for any difference; \code{"greater"}
#'   looks for group 1 being higher (requires
#'   \code{prevalence1 > prevalence2}); \code{"less"} looks for group 1 being
#'   lower (requires \code{prevalence1 < prevalence2}).
#' @param sensitivity Diagnostic sensitivity in (0, 1]: the share of people
#'   who truly have the condition that the test detects as positive.
#'   Default 1 (the test never misses a case). Assumed the same in both
#'   groups.
#' @param specificity Diagnostic specificity in (0, 1]: the share of people
#'   who truly do not have the condition that the test correctly reads as
#'   negative. Default 1 (the test never gives a false positive). Assumed the
#'   same in both groups.
#' @param conf_level Confidence level, in \[0.5, 1); default 0.95. The
#'   significance level of the test is 1 - \code{conf_level}.
#' @param n_sites Optional positive integer. Fixes the number of sites per
#'   group; the function solves for the number of people per site. Cannot be
#'   used together with \code{n_per_site}.
#' @param n_per_site Optional positive integer. Fixes the number of people
#'   per site; the function solves for the number of sites per group. Cannot
#'   be used together with \code{n_sites}. If \code{n_per_site} is more than
#'   the number of people a group needs, one site is already enough and the
#'   function stops with an error; use \code{n_sites = 1} to find how many
#'   people that one site needs.
#' @param icc Intra-cluster correlation (ICC): how alike people in the same
#'   site are, from 0 (not at all) to 1 (identical). Default 0 (no
#'   clustering). An \code{icc} above 0 needs \code{n_sites} or
#'   \code{n_per_site}, and the value given is then used exactly, however
#'   small. Without sites, a tiny value below about 1.5e-8, such as one left
#'   over from rounding in an earlier calculation, is treated as 0.
#' @param fpc_N Optional positive integer. Population size \strong{per
#'   group}, for a finite-population correction. \code{NULL} (default) = no
#'   finite-population correction applied.
#'
#' @details
#' \strong{Sample size per group.} With no clustering, the number of people
#' needed per group for a two-proportion z-test is
#'
#' \deqn{n_0 = \frac{\left(z_\alpha \sqrt{2\,\bar{p}(1-\bar{p})} +
#'   z_\beta \sqrt{p_{1,app}(1-p_{1,app}) + p_{2,app}(1-p_{2,app})}\right)^2}
#'   {(p_{1,app} - p_{2,app})^2}}
#'
#' where \eqn{p_{g,app} = p_g \cdot Se + (1-p_g)(1-Sp)} is the prevalence the
#' test would measure in group \eqn{g}, \eqn{\bar{p} = (p_{1,app} +
#' p_{2,app})/2}, \eqn{z_\beta} = \code{qnorm(power)} (0.84 for 80 percent
#' power), and \eqn{z_\alpha} = \code{qnorm(1 - alpha/2)} for a two-sided
#' test (1.96 at 95 percent) or \code{qnorm(1 - alpha)} for a one-sided test
#' (1.64 at 95 percent), with \eqn{\alpha = 1 -} \code{conf_level}. This is
#' the same formula as R's \code{power.prop.test()}.
#'
#' Because \eqn{p_{1,app} - p_{2,app} = (p_1 - p_2)(Se + Sp - 1)}, an
#' imperfect test shrinks the difference the test can see and increases
#' \eqn{n_0}. When \eqn{Se = Sp = 1} the measured and true prevalences are
#' the same.
#'
#' \strong{Clustering.} The design effect \eqn{D_{eff} = 1 + (\bar{n} -
#' 1)\,ICC}, with \eqn{\bar{n}} the number of people per site, is applied to
#' each group. With \code{n_per_site} fixed, \eqn{n = n_0 \cdot D_{eff}}. With
#' \code{n_sites} fixed, \eqn{D_{eff}} itself depends on \eqn{n}, and solving
#' \eqn{n = n_0 (1 + (n/n_{sites} - 1) ICC)} gives
#'
#' \deqn{n = \frac{n_0 \cdot n_{sites} \cdot (1-ICC)}{n_{sites} - n_0 \cdot ICC}}
#'
#' This has a solution only when \eqn{n_{sites} > n_0 \cdot ICC}. Otherwise
#' the target power cannot be reached, and the function stops and reports
#' the maximum power those sites can give. People in the same site
#' tend to resemble each other, so after a point each extra person adds
#' almost nothing new. Each site is worth at most \eqn{1 / ICC} independent
#' people, so a group can never count as more than \eqn{n_{sites} / ICC}
#' independent people. For example, with ICC 0.05 each site is worth at
#' most 20 independent people, so 8 sites per group give at most 160
#' independent people. For 10 versus 20 percent, 160 people per group give
#' 70.9 percent power, so 80 percent can never be reached.
#'
#' The number of people per site is rounded up, and the design effect
#' (returned as \code{deff}) is then worked out again from that rounded
#' number, so that
#' \eqn{D_{eff} = 1 + (n_{per site} - 1) ICC} describes the design reported.
#'
#' \strong{Finite-population correction.} With \code{fpc_N} (the population
#' size \eqn{N} of each group), sampling without replacement reduces the
#' variance by the factor \eqn{(N - n)/(N - 1)}, which gives
#' \eqn{n = n_1 N / (N + n_1 - 1)}, where \eqn{n_1} is the size before the
#' correction (Cochran, 1977). With \code{n_sites} fixed, the correction and
#' the design effect both depend on \eqn{n}, so they are solved together: the
#' equation \eqn{n_0 (1 + (n/n_{sites} - 1) ICC)(N - n)/(N - 1) = n} is a
#' quadratic, solved for its root between 0 and \eqn{N}. With a
#' finite population the target is therefore always reachable. If there are
#' too few sites for this ICC to reach the target power on their own, the
#' correction can still make it reachable, but often only by testing a large
#' share of the population; a warning says so.
#'
#' \code{fpc_N} counts people, not sites. If a study samples nearly all of
#' the sites in a group, the returned \eqn{n} is larger than necessary,
#' because the correction does not remove the variation between sites.
#'
#' @return A named list with the following fields, in the order returned:
#'   \item{n_per_group}{Number of people needed in each group.}
#'   \item{n}{Total number of people across both groups,
#'     \code{2 * n_per_group}. With sites, every site gets the same whole
#'     number of people, rounded up, so the study collects
#'     \code{2 * n_sites * n_per_site} people, which can be more than
#'     \code{n}. Budget for that number. The extra is usually small, but can
#'     almost double the study when sites have only 1 or 2 people each: for
#'     example, a group needing 199 people across 198 sites gets 2 per site,
#'     so 396 are collected in that group instead of 199.}
#'   \item{n_eff}{\eqn{n_0} rounded up: the number of people per group that
#'     would be needed with no clustering and no finite-population
#'     correction.}
#'   \item{n_sites}{Sites per group: solved for when \code{n_per_site} is
#'     given, echoed when \code{n_sites} is given. If the study needs fewer
#'     people per group than \code{n_sites}, only that many sites are used,
#'     with 1 person at each, and a warning says so. \code{NULL} with no
#'     clustering.}
#'   \item{n_per_site}{People per site: solved for when \code{n_sites} is
#'     given, echoed when \code{n_per_site} is given. \code{NULL} with no
#'     clustering.}
#'   \item{prevalence1, prevalence2}{The true prevalences, as supplied.}
#'   \item{apparent_prev1, apparent_prev2}{The prevalences the test would
#'     measure in each group.}
#'   \item{difference}{The true difference, \code{prevalence1 - prevalence2}.}
#'   \item{difference_app}{The difference the test would measure,
#'     \code{apparent_prev1 - apparent_prev2}.}
#'   \item{power}{Target power, as supplied.}
#'   \item{alternative}{Which kind of difference was looked for, as
#'     supplied.}
#'   \item{conf_level}{Confidence level, as supplied.}
#'   \item{sensitivity}{Sensitivity, as supplied.}
#'   \item{specificity}{Specificity, as supplied.}
#'   \item{icc}{ICC, as supplied.}
#'   \item{deff}{Design effect: how much clustering inflates the number of
#'     people needed (1 with no clustering). With sites, it is worked out for
#'     the people actually collected, \code{n_sites * n_per_site} per group.
#'     Dividing that number by \code{deff} gives at least \code{n_eff}
#'     (without \code{fpc_N}); dividing \code{n_per_group} by \code{deff} can
#'     give slightly less, because \code{n_per_group} comes before the
#'     rounding up to whole sites.}
#'   \item{fpc_N}{\code{fpc_N}, as supplied, or \code{NULL}.}
#'
#' @references
#' MMS-SD Study Design Workshop. \url{https://mrc-ide.github.io/MMS-SD_workshop/}
#'
#' Fleiss, J. L., Levin, B. & Paik, M. C. (2003) \emph{Statistical Methods for
#' Rates and Proportions}, 3rd ed. Wiley. (Two-group sample size.)
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
design_difference <- function(prevalence1,
                              prevalence2,
                              power       = 0.80,
                              alternative = "two.sided",
                              sensitivity = 1,
                              specificity = 1,
                              conf_level  = 0.95,
                              n_sites     = NULL,
                              n_per_site  = NULL,
                              icc         = 0,
                              fpc_N       = NULL) {

  # ---- validate single-value parameters ----
  # Each of the parameters below must be a single value.
  if (length(prevalence1) != 1 || !is.numeric(prevalence1))
    stop("`prevalence1` must be a single number in [0, 1] (got class `",
         class(prevalence1)[1], "`, length ", length(prevalence1), ").")
  if (length(prevalence2) != 1 || !is.numeric(prevalence2))
    stop("`prevalence2` must be a single number in [0, 1] (got class `",
         class(prevalence2)[1], "`, length ", length(prevalence2), ").")
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
    stop("`alternative` must be one of 'two.sided', 'greater', or 'less' ",
         "(got class `", class(alternative)[1], "`, length ",
         length(alternative), ").")
  if (!alternative %in% c("two.sided", "greater", "less"))
    stop("`alternative` must be one of 'two.sided', 'greater', or 'less' ",
         "(got '", alternative, "').")

  # Check for NA, NaN and Inf before the range checks below.
  if (!is.finite(prevalence1))
    stop("`prevalence1` must be a single finite number (got ", prevalence1, ").")
  if (!is.finite(prevalence2))
    stop("`prevalence2` must be a single finite number (got ", prevalence2, ").")
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

  for (nm in c("prevalence1", "prevalence2")) {
    p <- get(nm)
    if (p < 0 || p > 1)
      stop("`", nm, "` must be in [0, 1] (got ", num(p), "). `", nm, "` is a ",
           "fraction of the population and cannot be ",
           if (p < 0) "negative." else "greater than 1.")
  }
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
         num(sensitivity + specificity), "). ",
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
         "Below 0.5, when there is no real difference, the test would ",
         "wrongly report one more than half the time.")
  # Power can never be below the false-alarm rate: even with no real
  # difference the test reports one that often. A power equal to it only up
  # to rounding (within 1e-12, e.g. 1 - 0.8 = 0.19999999999999996) counts
  # as equal.
  false_alarm <- 1 - conf_level
  if (power <= false_alarm + 1e-12)
    stop("`power` must be greater than 1 - `conf_level` = ",
         signif(false_alarm, 6), " (got ", num(power), "). Even when there is no ",
         "real difference, the test reports one ", signif(100 * false_alarm, 6),
         " percent of the time, so a lower power is not meaningful.")
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
         else num(n_sites), "). ",
         "`n_sites` is the number of sites per group.")
  if (!is.null(n_per_site) &&
      (!is.numeric(n_per_site) || length(n_per_site) != 1 || !is.finite(n_per_site) ||
       n_per_site != floor(n_per_site) || n_per_site < 1))
    stop("`n_per_site` must be a single finite positive integer (got ",
         if (!is.numeric(n_per_site)) paste0("class `", class(n_per_site)[1], "`")
         else if (length(n_per_site) != 1) paste0("length = ", length(n_per_site))
         else num(n_per_site), "). ",
         "`n_per_site` is the fixed number of people sampled per site.")
  if (!is.null(fpc_N) &&
      (!is.numeric(fpc_N) || length(fpc_N) != 1 || !is.finite(fpc_N) ||
       fpc_N < 1 || fpc_N != floor(fpc_N)))
    stop("`fpc_N` must be a single finite positive integer (got ",
         if (!is.numeric(fpc_N)) paste0("class `", class(fpc_N)[1], "`")
         else if (length(fpc_N) != 1) paste0("length = ", length(fpc_N))
         else num(fpc_N),
         "). `fpc_N` is the population size per group. ",
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
  prevalence1 <- as.numeric(prevalence1); prevalence2 <- as.numeric(prevalence2)
  power       <- as.numeric(power);       conf_level  <- as.numeric(conf_level)
  sensitivity <- as.numeric(sensitivity); specificity <- as.numeric(specificity)
  icc         <- as.numeric(icc)
  if (!is.null(n_sites))    n_sites    <- as.numeric(n_sites)
  if (!is.null(n_per_site)) n_per_site <- as.numeric(n_per_site)
  if (!is.null(fpc_N))      fpc_N      <- as.numeric(fpc_N)

  # The number of sites, or people per site, cannot exceed the population.
  if (!is.null(fpc_N) && !is.null(n_sites) && n_sites > fpc_N)
    stop("`n_sites` (", format(n_sites, scientific = FALSE), ") exceeds `fpc_N` (",
         format(fpc_N, scientific = FALSE), "): you cannot have more sites than ",
         "people in the population.")
  if (!is.null(fpc_N) && !is.null(n_per_site) && n_per_site > fpc_N)
    stop("`n_per_site` (", format(n_per_site, scientific = FALSE), ") exceeds `fpc_N` (",
         format(fpc_N, scientific = FALSE), "): a site cannot have more people ",
         "than there are in the whole population.")

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
  if (abs(prevalence1 - prevalence2) < 1e-8)
    stop("`prevalence1` (", num(prevalence1), ") and `prevalence2` (", num(prevalence2),
         ") are equal, or nearly so, so there is no difference to detect. ",
         "Use two different prevalences.")
  if (alternative == "greater" && prevalence1 < prevalence2)
    stop("`prevalence1` (", num(prevalence1), ") must be greater than `prevalence2` (",
         num(prevalence2), ") when alternative = 'greater'. Use alternative = 'less'.")
  if (alternative == "less" && prevalence1 > prevalence2)
    stop("`prevalence1` (", num(prevalence1), ") must be less than `prevalence2` (",
         num(prevalence2), ") when alternative = 'less'. Use alternative = 'greater'.")

  # Print two different numbers with enough digits that they never look
  # the same (79.9999993 next to a target of 80 would show as 80).
  format_apart <- function(a, b, digits) {
    repeat {
      fa <- format(a, digits = digits, scientific = FALSE)
      fb <- format(b, digits = digits, scientific = FALSE)
      if (fa != fb || digits >= 15) return(c(fa, fb))
      digits <- digits + 1
    }
  }

  # ---- apparent (measured) prevalences ----
  p1_app <- .apparent_prev(prevalence1, sensitivity, specificity)
  p2_app <- .apparent_prev(prevalence2, sensitivity, specificity)
  pbar   <- (p1_app + p2_app) / 2
  delta  <- abs(p1_app - p2_app)
  # An imperfect test can bring two different prevalences very close
  # together on the measured scale.
  if (delta < 1e-8) {
    shown <- format_apart(p1_app, p2_app, digits = 6)
    stop("The measured prevalences differ by only ",
         format(signif(delta, 3), scientific = FALSE),
         " (", shown[1], " and ", shown[2], "), so the ",
         "required sample size would be extremely large. Use ",
         "prevalences further apart, or a more accurate test.")
  }

  # ---- warnings, once every check above has passed ----
  # A sum of 1.1 that comes out a hair below it by rounding (0.15 + 0.95)
  # does not warn.
  if (correction < 0.1 - 1e-12)
    warning("`sensitivity` + `specificity` = ", num(sensitivity + specificity),
            ", which is very close to 1. The test then barely tells apart ",
            "people with and without the condition, so the required n is very ",
            "large, and small changes in the assumed sensitivity or ",
            "specificity will change n a lot.")
  if (power < 0.5)
    warning("`power` = ", num(power), " is below 0.5, so the study would be more ",
            "likely to miss a real difference than to find it.")

  # ---- sample size per group with no clustering ----
  alpha <- 1 - conf_level
  z_a   <- if (alternative == "two.sided") stats::qnorm(1 - alpha / 2) else
                                           stats::qnorm(1 - alpha)
  z_b   <- stats::qnorm(power)
  se0   <- sqrt(2 * pbar * (1 - pbar))
  se1   <- sqrt(p1_app * (1 - p1_app) + p2_app * (1 - p2_app))
  n_base_cont <- ((z_a * se0 + z_b * se1) / delta)^2
  # Only a `conf_level` so close to 1 that R treats z as infinite gets here.
  if (!is.finite(n_base_cont))
    stop("The required sample size is too large for R to compute. Use a lower ",
         "`conf_level`.")

  # ---- design effect and n per group ----
  # If icc is above 0, `n_sites` or `n_per_site` was given (checked above). With
  # n_sites fixed, the finite-population correction is part of the solve;
  # otherwise it is applied afterwards.
  fpc_applied <- FALSE
  # TRUE when the sites alone could not reach the target and only the
  # finite-population correction makes it reachable (warned about below).
  fpc_rescue  <- FALSE
  # TRUE when more sites were given than people needed, so some are dropped.
  sites_cut   <- FALSE

  if (icc_is_zero) {
    deff   <- 1
    n_cont <- n_base_cont

  } else if (!is.null(n_per_site)) {
    deff   <- 1 + (n_per_site - 1) * icc
    n_cont <- n_base_cont * deff

  } else {
    # Number of sites fixed: deff depends on n, and n depends on deff.
    if (is.null(fpc_N)) {
      # Solving n = n_base * (1 + (n / n_sites - 1) * icc) for n gives
      #   n = n_base * n_sites * (1 - icc) / (n_sites - n_base * icc),
      # which needs n_sites > n_base * icc. The highest power these sites
      # can give is where n_base * icc = n_sites exactly.
      denom <- n_sites - n_base_cont * icc
      if (denom <= 0) {
        z_b_max   <- (delta * sqrt(n_sites / icc) - z_a * se0) / se1
        max_power <- stats::pnorm(z_b_max)
        shown <- format_apart(100 * power, 100 * max_power, digits = 7)
        # A maximum like 2e-16 percent is shown as "less than 0.0001".
        if (100 * max_power < 1e-4) shown[2] <- "less than 0.0001"
        stop("Target power of ", shown[1], " percent cannot be ",
             "reached with ", format(n_sites, scientific = FALSE),
             if (n_sites == 1) " site" else " sites", " per group and ICC = ",
             num(icc), ". The highest power ",
             if (n_sites == 1) "this site" else "these sites", " can give is ",
             shown[2], " percent. ",
             "Increase `n_sites`, lower the ICC, or lower the target power.")
      }
      n_cont <- n_base_cont * n_sites * (1 - icc) / denom

    } else {
      # With an FPC, solve n_base * (1 + (n / n_sites - 1) * icc) *
      # (N - n) / (N - 1) = n for n. This is the quadratic A*n^2 + B*n + C = 0
      # with A > 0 and C < 0, so it has exactly one positive root, below N.
      # The root is taken in a form that stays accurate for very large N.
      # With icc = 1, C = 0: the roots are 0 and -B / A. (0 people needed,
      # possible with 0 versus 1 prevalence, stays 0.)
      A <- n_base_cont * icc / n_sites
      B <- n_base_cont * (1 - icc) - A * fpc_N + fpc_N - 1
      C <- -n_base_cont * (1 - icc) * fpc_N
      q <- -(B + (if (B < 0) -1 else 1) * sqrt(B^2 - 4 * A * C)) / 2
      n_cont <- if (n_base_cont == 0) 0
                else if (C == 0) max(0, -B / A)
                else max(q / A, C / q)
      fpc_applied <- TRUE
      fpc_rescue  <- n_sites <= n_base_cont * icc
    }
    deff <- 1 + (n_cont / n_sites - 1) * icc

    # With icc at least 0, deff is below 1 only when the solve gives fewer
    # than one person per site: more sites than people needed with no
    # clustering.
    if (deff < 1) {
      n_srs <- if (is.null(fpc_N) || n_base_cont == 0) n_base_cont
               else (n_base_cont * fpc_N) / (n_base_cont + fpc_N - 1)
      # Use 1 person per site, so no clustering: the whole number of people
      # needed with no clustering (198.96 -> 199), or n_sites if that is
      # smaller. Extra sites are dropped below, with a warning.
      n_cont <- min(n_sites, ceiling(n_srs * (1 - 1e-12)))
      deff   <- 1
    }
  }

  # ---- finite-population correction (per group) ----
  # (0 people needed, possible with 0 versus 1 prevalence, stays 0.)
  if (!is.null(fpc_N) && !fpc_applied && n_cont > 0)
    n_cont <- (n_cont * fpc_N) / (n_cont + fpc_N - 1)

  # Above a trillion people, rounding error in the calculation itself is
  # more than one person, so the answer could not be given exactly.
  if (n_cont > 1e12) {
    fixes <- c("prevalences further apart", "a lower `power` or `conf_level`")
    if (!icc_is_zero && !is.null(n_per_site)) fixes <- c(fixes, "a smaller `n_per_site`")
    if (!icc_is_zero && !is.null(n_sites))    fixes <- c(fixes, "more sites or a lower ICC")
    if (correction < 1)                       fixes <- c(fixes, "a more accurate test")
    stop("The required sample size is more than a trillion (1000000000000) ",
         "people per group, too large to work out to the nearest person. Use ",
         paste(fixes[-length(fixes)], collapse = ", "), ", or ",
         fixes[length(fixes)], ".")
  }

  # Round up to whole people. Shrinking by one part in a trillion first stops
  # rounding error just above a whole number (300.00000000000006) adding a
  # whole person; it could only remove a real person above a trillion
  # people, which is stopped above.
  # At least 1 person per group: with 0 versus 1 prevalence the formula can
  # give 0.
  # Never more people than the population: rounding error in the correction
  # could otherwise give 2 people from a population of 1.
  shrink      <- 1 - 1e-12
  n_per_group <- max(1, ceiling(n_cont * shrink))
  if (!is.null(fpc_N)) n_per_group <- min(n_per_group, fpc_N)
  n_total     <- 2 * n_per_group

  # n_eff: people per group needed with no clustering and no FPC.
  n_eff <- max(1, ceiling(n_base_cont * shrink))

  # ---- distribute across sites ----
  if (!is.null(n_per_site)) {
    # A site larger than a whole group would be only partly sampled, so its
    # size would not be the real site size used in the design effect.
    if (n_per_site > n_per_group)
      stop("`n_per_site` (", format(n_per_site, scientific = FALSE), ") is ",
           "larger than the number of people needed per group (",
           format(n_per_group, scientific = FALSE), "), so one site would be ",
           "more than enough. Use a smaller `n_per_site`, or set ",
           "`n_sites = 1` to find how many people one site needs.")
    n_sites_out    <- ceiling(n_per_group / n_per_site)
    n_per_site_out <- n_per_site
  } else if (!is.null(n_sites)) {
    n_sites_out <- n_sites
    # More sites than people needed: 1 person per site is enough, so only
    # as many sites as people needed are used.
    if (n_sites > n_per_group) {
      n_sites_out <- n_per_group
      sites_cut   <- TRUE
      one <- n_per_group == 1
      warning("With these prevalences and this power, the study needs ",
              format(n_per_group, scientific = FALSE),
              if (one) " person" else " people", " per group, fewer than the ",
              format(n_sites, scientific = FALSE), " sites given (`n_sites`). ",
              "The design uses ", format(n_per_group, scientific = FALSE),
              if (one) " site" else " sites", " with 1 person at each.")
    }
    n_per_site_out <- ceiling(n_per_group / n_sites_out)
    # Work deff out again from the rounded people per site, so that
    # deff == 1 + (n_per_site - 1) * icc describes the reported design.
    if (!icc_is_zero) deff <- 1 + (n_per_site_out - 1) * icc
  } else {
    n_sites_out    <- NULL
    n_per_site_out <- NULL
  }

  # Rounding the number of sites or people per site up can push the design
  # past the population. That only happens close to a census.
  if (!is.null(fpc_N) && !is.null(n_sites_out) &&
      n_sites_out * n_per_site_out > fpc_N)
    warning("The rounded design (",
            format(n_sites_out, scientific = FALSE), " sites x ",
            format(n_per_site_out, scientific = FALSE), " = ",
            format(n_sites_out * n_per_site_out, scientific = FALSE),
            " per group) exceeds `fpc_N` = ", format(fpc_N, scientific = FALSE),
            ": the required n of ", format(n_per_group, scientific = FALSE),
            " per group is close to a census of the population. Consider ",
            "surveying the whole population instead.")

  # Without the population, these sites could not reach the target at all,
  # so the design can be very large.
  if (fpc_rescue)
    warning("With ", format(n_sites, scientific = FALSE),
            if (n_sites == 1) " site" else " sites", " per group and ICC = ",
            num(icc), ", the target power can only be reached because of ",
            "the finite-population correction: the design tests ",
            format(n_per_group, scientific = FALSE), " of the ",
            format(fpc_N, scientific = FALSE), " people in each group (",
            format(n_per_site_out, scientific = FALSE), " per site). Check that ",
            "`fpc_N` is right", if (sites_cut) "." else ", or use more sites.")

  list(
    n_per_group    = n_per_group,
    n              = n_total,
    n_eff          = n_eff,
    n_sites        = n_sites_out,
    n_per_site     = n_per_site_out,
    prevalence1    = prevalence1,
    prevalence2    = prevalence2,
    apparent_prev1 = p1_app,
    apparent_prev2 = p2_app,
    difference     = prevalence1 - prevalence2,
    difference_app = p1_app - p2_app,
    power          = power,
    alternative    = alternative,
    conf_level     = conf_level,
    sensitivity    = sensitivity,
    specificity    = specificity,
    icc            = icc,
    deff           = deff,
    fpc_N          = fpc_N
  )
}
