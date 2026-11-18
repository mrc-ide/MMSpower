#' Test whether prevalence is above or below a decision threshold
#'
#' @description
#' Analysis function: given the number of positives and the number tested
#' at one or more sites, tests whether the true prevalence is above (or
#' below) a decision threshold, with an exact binomial test. Returns the
#' p-value, whether the result is significant, and a confidence interval
#' for the true prevalence that agrees with the test. It is the analysis
#' twin of \code{design_threshold()}.
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
#' @param x Vector of integer counts of people who tested positive, one per
#'   site (a single number for one site).
#' @param n Vector of integer totals tested, one per site (same length as
#'   \code{x}).
#' @param threshold Decision threshold for the true prevalence, strictly
#'   between 0 and 1 (for example 0.05 for the WHO pfhrp2/3 deletion
#'   threshold).
#' @param alternative Direction of the test: \code{"greater"} (default,
#'   prevalence above the threshold), \code{"less"} (below it) or
#'   \code{"two.sided"} (different from it).
#' @param sensitivity Diagnostic sensitivity in (0, 1]: the share of people
#'   who truly have the condition that the test detects as positive.
#'   Default 1 (the test never misses a case).
#' @param specificity Diagnostic specificity in (0, 1]: the share of people
#'   who truly do not have the condition that the test correctly reads as
#'   negative. Default 1 (the test never gives a false positive).
#' @param conf_level Confidence level for the interval on the true
#'   prevalence, in \[0.5, 1); default 0.95. The significance level of the
#'   test is 1 - \code{conf_level}.
#' @param icc Optional. Intra-cluster correlation (ICC): how alike people in
#'   the same site are, from 0 (not at all) to 1 (identical). If
#'   \code{NULL} (default), it is estimated from the sites. A number is
#'   used exactly as given, however small: with very large sites even a
#'   tiny ICC changes the design effect. It is also used with only one site
#'   (see Details). Set to \code{0} for no clustering.
#' @param fpc_N Optional. Population size for a finite-population
#'   correction: a positive integer, larger than the total number tested.
#'   \code{NULL} (default) means no finite-population correction.
#' @param ci_method Method for the interval on the true prevalence:
#'   \code{"clopper-pearson"} (default), the interval that matches the
#'   exact test, so that \code{reject} is \code{TRUE} exactly when
#'   \code{threshold} is outside it; \code{"agresti-coull"} (a Wald interval
#'   after adding the same number of extra positives and negatives); or
#'   \code{"wald"} (the measured prevalence plus or minus z standard
#'   errors). \code{"agresti-coull"} and \code{"wald"} can occasionally leave
#'   out the threshold when the test does not reject, or include it when it
#'   does, and \code{"wald"} can be too narrow with small samples or a
#'   prevalence near 0 (see Details). It changes only the interval, not the
#'   test.
#'
#' @details
#' \strong{Design effect} (returned as \code{deff}). It is
#' \eqn{1 + (n_w - 1) ICC}, where \eqn{n_w = \sum n^2 / \sum n} is the
#' mean number of people per site, weighted by site size (the plain mean
#' when all sites are the same size). With \code{icc} supplied (for example from earlier
#' surveys), that value is used. Supply \code{icc} with one site only when
#' the site stands for a wider population, such as a district; if the site
#' itself is the population of interest, leave \code{icc = NULL}. With
#' \code{icc = NULL}, the design effect is worked out from the data as in
#' \code{estimate_prevalence()}: the variance of the site prevalences
#' divided by the variance expected if everyone were independent. This
#' needs at least two sites; with one site the design effect is set to 1.
#' The ICC that goes with it (returned as \code{icc_used}) is kept between
#' 0 and 1, so the design effect is between 1 and \eqn{n_w}. With one
#' person per site the design effect is 1.
#'
#' \strong{Effective sample size.} \eqn{m} (returned as \code{n_eff_adj})
#' is the number tested divided by the design effect, then divided by the
#' finite-population factor \eqn{(N - n)/(N - 1)} when \code{fpc_N} is
#' given (\eqn{n} the number tested, \eqn{N} the population size). The
#' number of positives is scaled the same way, to \eqn{\hat{p} m}
#' (returned as \code{statistic}), where \eqn{\hat{p} = \sum x / \sum n}
#' is the measured prevalence.
#'
#' \strong{Test.} Suppose the true prevalence is exactly at the threshold. For
#' \code{"greater"}, the p-value is the chance of at least as many positives
#' as observed; for \code{"less"}, the chance of at most as many. The
#' threshold is first converted to the measured scale:
#' \eqn{threshold \times Se + (1 - threshold)(1 - Sp)}. With one site and no
#' \code{fpc_N}, the one-sided p-values are exactly those of R's
#' \code{binom.test()}.
#'
#' \strong{Two-sided p-value.} Twice the smaller of the two one-sided
#' p-values, at most 1. This matches the two-sided Clopper-Pearson
#' interval (the "central" two-sided exact test; Fay 2010, R Journal 2(1)).
#' It is cautious: the false alarm rate is at most
#' 1 - \code{conf_level}, and often well below it. It can differ from
#' \code{binom.test()}, which adds up all counts that are less likely than
#' the one observed.
#'
#' \strong{Interval.} The interval is matched to \code{alternative}:
#' \eqn{[L, 1]} for \code{"greater"}, \eqn{[0, U]} for \code{"less"}, and
#' the two-sided interval for \code{"two.sided"}. The Clopper-Pearson
#' interval (default) holds the thresholds the test would not reject,
#' so \code{reject} is \code{TRUE} exactly when \code{threshold} is outside
#' it. It is built on the measured scale and then corrected for an
#' imperfect test, with ends below 0 or above 1 set to 0 or 1.
#'
#' \code{"agresti-coull"} adds \eqn{z^2/2} positives and \eqn{z^2/2}
#' negatives (about 1.92 of each for a two-sided 95 percent interval) and
#' then works out a Wald interval. It can disagree with \code{reject}
#' near the threshold.
#'
#' \code{"wald"} is \eqn{\hat{p} \pm z \sqrt{\hat{p}(1 - \hat{p})/m}}. It
#' can disagree with \code{reject}, and with a prevalence near 0 it can be
#' too narrow: with 100 people and a true prevalence of 5 percent, a
#' two-sided 95 percent interval contains the true prevalence only about 88
#' percent of the time (Clopper-Pearson 98, Agresti-Coull 97). When the
#' interval has zero width, a warning is given. That happens at a measured
#' prevalence of 0 unless \code{alternative} is \code{"greater"}, and at 1
#' unless it is \code{"less"}.
#'
#' When \eqn{Se + Sp} is below 1.1 a warning is given: the correction divides by
#' \eqn{Se + Sp - 1}, so the prevalence and its interval become very wide.
#'
#' \strong{When to use DRpower instead.} This function allows for imperfect
#' diagnostic tests (sensitivity and specificity below 1), a finite
#' population (\code{fpc_N}), and questions about prevalence below the
#' threshold or different from it. With the ICC supplied, it keeps the
#' chance of a false alarm at or below the significance level. With
#' \code{icc = NULL} and few sites of very different sizes, the ICC is
#' poorly estimated and the false alarm rate can be well above it (13 to 28
#' percent at a 5 percent level in simulation); with an ICC estimated from
#' few sites it also has low power. For pfhrp2/3 deletion surveys with several sites
#' and an unknown ICC, where the question is whether prevalence is above
#' the threshold, consider DRpower
#' (\url{https://mrc-ide.github.io/DRpower/}). It estimates prevalence and
#' ICC together and has more power, especially with few sites, but more
#' false alarms when prevalence is close to the threshold.
#'
#' @return A named list with the following fields, in the order returned:
#'   \item{statistic}{Effective number of positives the exact test uses: the
#'     measured prevalence times \code{n_eff_adj}. It equals \code{sum(x)}
#'     when the design effect is 1 and \code{fpc_N} is \code{NULL}.}
#'   \item{p_value}{p-value for the chosen \code{alternative}.}
#'   \item{reject}{\code{TRUE} if \code{p_value} is below
#'     1 - \code{conf_level}.}
#'   \item{threshold}{Decision threshold, as supplied.}
#'   \item{alternative}{Direction of the test, as supplied.}
#'   \item{prevalence}{Estimated true prevalence (Rogan-Gladen corrected,
#'     between 0 and 1).}
#'   \item{ci_lower, ci_upper}{Confidence interval for the true prevalence;
#'     \code{ci_upper} is 1 for \code{"greater"} and \code{ci_lower} is 0
#'     for \code{"less"}.}
#'   \item{ci_method}{Interval method, as supplied.}
#'   \item{n_total}{Total number of people tested, \code{sum(n)}.}
#'   \item{n_eff}{Effective sample size before the finite-population
#'     correction: \code{n_total / deff}.}
#'   \item{n_eff_adj}{Effective sample size the test and interval use, after
#'     the finite-population correction (equals \code{n_eff} when
#'     \code{fpc_N} is \code{NULL}).}
#'   \item{conf_level}{Confidence level, as supplied.}
#'   \item{sensitivity, specificity}{As supplied.}
#'   \item{icc_used}{ICC used: the value supplied or, if \code{icc = NULL},
#'     the value estimated from the data (0 with one site, or one person per
#'     site).}
#'   \item{deff}{Design effect used (1 with no clustering).}
#'   \item{fpc_N}{\code{fpc_N}, as supplied, or \code{NULL}.}
#'
#' @references
#' MMS-SD Study Design Workshop. \url{https://mrc-ide.github.io/MMS-SD_workshop/}
#'
#' Clopper, C. J. & Pearson, E. S. (1934) The use of confidence or fiducial
#' limits illustrated in the case of the binomial. \emph{Biometrika}
#' \strong{26}(4), 404-413. \doi{10.1093/biomet/26.4.404}
#'
#' Agresti, A. & Coull, B. A. (1998) Approximate is better than "exact" for
#' interval estimation of binomial proportions. \emph{The American
#' Statistician} \strong{52}(2), 119-126. \doi{10.2307/2685469}
#'
#' Fay, M. P. (2010) Two-sided exact tests and matching confidence intervals
#' for discrete data. \emph{The R Journal} \strong{2}(1), 53-58.
#' \doi{10.32614/RJ-2010-008}
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
test_threshold <- function(x,
                           n,
                           threshold,
                           alternative = "greater",
                           sensitivity = 1,
                           specificity = 1,
                           conf_level  = 0.95,
                           icc         = NULL,
                           fpc_N       = NULL,
                           ci_method   = "clopper-pearson") {

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

  # ---- validate x and n ----
  # The same checks and messages as test_difference() and
  # estimate_prevalence().
  if (is.logical(x) || is.logical(n))
    stop("`x` and `n` must be numeric, not logical (got class `",
         class(x)[1], "` for x, `", class(n)[1], "` for n). ",
         "Note: a plain `NA` is logical in R. Remove missing observations ",
         "before calling.")
  if (!is.numeric(x) || !is.numeric(n))
    stop("`x` and `n` must be numeric (got class `", class(x)[1], "` for x, ",
         "`", class(n)[1], "` for n).")
  if (length(x) == 0 || length(n) == 0)
    stop("`x` and `n` must each contain at least one value. ",
         "Supply at least one site's count and total.")
  # For each check below, `i` is the position of the first bad element.
  if (!all(is.finite(x))) {
    i <- which(!is.finite(x))[1]
    stop("`x` contains a missing or infinite value (found x[", i, "] = ",
         num(x[i]), "). The number of positives at each site must be a ",
         "finite non-negative integer.")
  }
  if (!all(is.finite(n))) {
    i <- which(!is.finite(n))[1]
    stop("`n` contains a missing or infinite value (found n[", i, "] = ",
         num(n[i]), "). The number tested at each site must be a finite ",
         "positive integer.")
  }
  if (length(x) != length(n))
    stop("`x` and `n` must have the same length ",
         "(got length(x) = ", length(x), ", length(n) = ", length(n), "). ",
         "Each value of `x` is the number of positives at one site, and ",
         "each value of `n` is the number tested at that site.")
  if (any(x < 0)) {
    i <- which(x < 0)[1]
    stop("`x` must be non-negative (found x[", i, "] = ", num(x[i]), "). ",
         "Counts cannot be negative.")
  }
  if (any(n <= 0)) {
    i <- which(n <= 0)[1]
    stop("`n` must be positive for every site (found n[", i, "] = ",
         num(n[i]), "). Each site must have tested at least one person.")
  }
  if (any(x != floor(x))) {
    i <- which(x != floor(x))[1]
    stop("`x` must contain integers, not decimals (found x[", i, "] = ",
         num(x[i]), near_int(x[i], "x"), ").")
  }
  if (any(n != floor(n))) {
    i <- which(n != floor(n))[1]
    stop("`n` must contain integers, not decimals (found n[", i, "] = ",
         num(n[i]), near_int(n[i], "n"), ").")
  }
  if (any(x > n)) {
    i <- which(x > n)[1]
    stop("`x` cannot exceed `n` (found x[", i, "] = ", num(x[i]),
         " > n[", i, "] = ", num(n[i]), "). The number of positives cannot ",
         "be more than the number tested at that site.")
  }
  # Above 2^53 R cannot store every integer exactly, and no real site tests
  # this many people.
  if (any(n > 2^53)) {
    i <- which(n > 2^53)[1]
    stop("`n` contains a value larger than 2^53 = 9007199254740992, ",
         "the largest whole number R can store exactly (found n[", i,
         "] = ", huge(n[i]), "). No real site tests this many people.")
  }

  # ---- validate single-value parameters ----
  # Each of the parameters below must be a single value.
  if (length(threshold) != 1 || !is.numeric(threshold))
    stop("`threshold` must be a single number in (0, 1) (got class `",
         class(threshold)[1], "`, length ", length(threshold), ").")
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
  if (!is.null(icc) && (length(icc) != 1 || !is.numeric(icc)))
    stop("`icc` must be a single number in [0, 1] (got class `",
         class(icc)[1], "`, length ", length(icc), "). ",
         "To estimate ICC from the data, leave `icc = NULL`.")
  if (!is.character(alternative) || length(alternative) != 1)
    stop("`alternative` must be one of 'greater', 'less', or 'two.sided' ",
         "(got class `", class(alternative)[1], "`, length ",
         length(alternative), ").")
  if (!alternative %in% c("greater", "less", "two.sided"))
    stop("`alternative` must be one of 'greater', 'less', or 'two.sided' ",
         "(got '", alternative, "').")
  if (!is.character(ci_method) || length(ci_method) != 1)
    stop("`ci_method` must be one of 'clopper-pearson', 'agresti-coull', or ",
         "'wald' (got class `", class(ci_method)[1], "`, length ", length(ci_method), ").")
  if (!ci_method %in% c("clopper-pearson", "agresti-coull", "wald"))
    stop("`ci_method` must be one of 'clopper-pearson', 'agresti-coull', or ",
         "'wald' (got '", ci_method, "').")

  # Check for NA, NaN and Inf before the range checks below.
  if (!is.finite(threshold))
    stop("`threshold` must be a single finite number (got ", num(threshold), ").")
  if (!is.finite(sensitivity))
    stop("`sensitivity` must be a single finite number (got ", num(sensitivity), ").")
  if (!is.finite(specificity))
    stop("`specificity` must be a single finite number (got ", num(specificity), ").")
  if (!is.finite(conf_level))
    stop("`conf_level` must be a single finite number (got ", num(conf_level), ").")
  if (!is.null(icc) && !is.finite(icc))
    stop("`icc` must be a single finite number (got ", num(icc), "). ",
         "To estimate ICC from the data, leave `icc = NULL`.")

  if (threshold < 0 || threshold > 1)
    stop("`threshold` must be in (0, 1) (got ", num(threshold), "). ",
         "`threshold` is a prevalence and cannot be ",
         if (threshold < 0) "negative." else "greater than 1.")
  if (threshold == 0 || threshold == 1)
    stop("`threshold` must be in (0, 1) (got ", num(threshold), "). ",
         "No prevalence can be below 0 or above 1, so there is nothing to ",
         "test. Use the prevalence that would lead to action, for example 0.05.")
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
    stop("`sensitivity` + `specificity` must exceed 1 for the Rogan-Gladen correction ",
         "(got ", num(sensitivity), " + ", num(specificity), " = ",
         # 12 digits, so a sum that is 1 up to rounding (within 1e-12) reads 1
         format(sensitivity + specificity, digits = 12, scientific = FALSE), "). ",
         "The correction divides by `sensitivity` + `specificity` - 1, so a sum ",
         "of 1 or less cannot be corrected.")
  if (!is.null(icc) && (icc < 0 || icc > 1))
    stop("`icc` must be in [0, 1] (got ", num(icc), "). ",
         "`icc` is a correlation and cannot be ",
         if (icc < 0) "negative. " else "greater than 1. ",
         "To estimate ICC from the data, leave `icc = NULL`.")
  if (conf_level >= 1)
    stop("`conf_level` must be in [0.5, 1) (got ", num(conf_level), "). ",
         "For example, use 0.95 for a 95% confidence interval.")
  if (conf_level < 0.5)
    stop("`conf_level` must be at least 0.5 (got ", num(conf_level), "). ",
         "Below 0.5, when the prevalence is exactly at the threshold, the ",
         "test can wrongly find it above or below the threshold more than ",
         "half the time.")
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

  # From here on these are plain numbers: a name or a 1 x 1 matrix would
  # otherwise leak into the result.
  threshold   <- as.numeric(threshold)
  sensitivity <- as.numeric(sensitivity)
  specificity <- as.numeric(specificity)
  conf_level  <- as.numeric(conf_level)
  if (!is.null(icc)) icc <- as.numeric(icc)
  if (!is.null(fpc_N)) fpc_N <- as.numeric(fpc_N)

  n_sites <- length(n)
  n_total <- sum(n)

  if (!is.null(fpc_N) && fpc_N < n_total)
    stop("`fpc_N` (", fpc_N, ") is smaller than the total sample size (",
         n_total, "): you cannot test more people than there are in the ",
         "population.")
  if (!is.null(fpc_N) && fpc_N == n_total)
    stop("`fpc_N` (", fpc_N, ") equals the total sample size: the whole ",
         "population was tested (a census), so there is no sampling ",
         "uncertainty and no confidence interval to compute. ",
         "Set `fpc_N = NULL` if no finite-population correction is needed.")

  # Warnings are collected and given only once every check above has
  # passed, so a bad input gives just its error.
  warn_later <- character(0)

  # A sum of 1.1 that comes out a hair below it by rounding (0.15 + 0.95)
  # does not warn.
  if (correction < 0.1 - 1e-12)
    warn_later <- c(warn_later, paste0(
      "`sensitivity` + `specificity` = ",
      format(sensitivity + specificity, digits = 15),
      ", which is below 1.1. Correcting for this much test error ",
      "makes the estimated prevalence and its confidence interval very ",
      "wide, and small changes in the assumed sensitivity or specificity ",
      "will change them a lot."))

  p_hat <- sum(x) / n_total   # measured prevalence

  # ---- design effect (as in estimate_prevalence) ----
  n_bar <- mean(n)   # mean number of people per site
  # Size-weighted mean site size: the pooled prevalence has variance
  # proportional to 1 + (n_w - 1) * ICC, and n_w = n_bar when sites are equal.
  n_w   <- sum(n^2) / sum(n)

  if (is.null(icc)) {
    if (n_sites < 2 || n_bar == 1) {
      # A single site (no between-site variation to measure) or every site
      # has one person (no clustering).
      icc_used <- 0
      deff     <- 1
    } else {
      var_obs <- stats::var(x / n)
      # Variance expected if everyone were independent, using the overall
      # prevalence.
      var_srs <- mean(p_hat * (1 - p_hat) / n)
      deff    <- if (var_srs > 0) max(var_obs / var_srs, 1) else 1
      icc_used <- min(max((deff - 1) / (n_bar - 1), 0), 1)
      deff     <- 1 + (n_w - 1) * icc_used
    }
  } else {
    # A supplied ICC is used as given, also with one site: people at one
    # site still resemble each other. With one person per site this gives
    # deff = 1.
    icc_used <- icc
    deff     <- 1 + (n_w - 1) * icc_used
  }

  # The Wald interval has zero width at a measured prevalence of 0 unless
  # `alternative` is "greater" (then it is [0, 1]), and at 1 unless it is
  # "less" (then it is [0, 1]).
  wald_zero_width <- (p_hat == 0 && alternative != "greater") ||
                     (p_hat == 1 && alternative != "less")
  if (ci_method == "wald" && wald_zero_width)
    warn_later <- c(warn_later, paste0(
      "The measured prevalence is ", num(p_hat), ", so the Wald interval ",
      "has zero width. Use `ci_method = \"clopper-pearson\"` for an interval ",
      "that allows for the uncertainty."))

  # Every check has passed: give the warnings collected above.
  for (w in warn_later) warning(w, call. = FALSE)

  n_eff <- n_total / deff

  # Finite-population correction, using the number of people tested (as in
  # estimate_prevalence): the effective size is divided by
  # (N - n_total) / (N - 1). fpc_N > n_total is guaranteed above.
  f         <- if (!is.null(fpc_N)) (fpc_N - n_total) / (fpc_N - 1) else 1
  n_eff_adj <- n_eff / f

  # ---- threshold on the measured scale ----
  theta_app <- .apparent_prev(threshold, sensitivity, specificity)

  # ---- exact binomial test ----
  # Effective size m = n_eff_adj and effective number of positives
  # p_hat * m (both may be fractions; with one site and no fpc_N they are
  # n_total and sum(x)). The exact binomial tail probabilities are worked
  # out through the beta distribution, which equals pbinom() for whole
  # numbers (as in Clopper-Pearson).
  m     <- n_eff_adj
  x_eff <- p_hat * m
  p_up  <- .binom_upper(x_eff, m, theta_app)   # chance of x_eff or more
  p_dn  <- .binom_lower(x_eff, m, theta_app)   # chance of x_eff or fewer

  p_value <- switch(alternative,
    greater   = p_up,
    less      = p_dn,
    two.sided = min(1, 2 * min(p_up, p_dn))
  )

  # 1 - 0.95 is stored as 0.050000000000000044; rounded to 15 digits it is
  # exactly 0.05, so an exact p-value of 0.05 is not counted as below it.
  alpha  <- signif(1 - conf_level, 15)
  reject <- p_value < alpha

  # ---- interval on the true prevalence, matched to `alternative` ----
  # "greater" -> [L, 1], "less" -> [0, U], "two.sided" -> [L, U] with
  # alpha / 2 on each side. Clopper-Pearson (default) is the inversion
  # of the test above, so `reject` is TRUE exactly when `threshold` is
  # outside it.
  alpha_lo <- switch(alternative,
    two.sided = alpha / 2, greater = alpha, less = 0)
  alpha_hi <- switch(alternative,
    two.sided = alpha / 2, greater = 0,     less = alpha)
  # z for each closed side; NA on an open side. Written this way, not as
  # qnorm(1 - alpha), so it stays exact when `conf_level` is very close to 1.
  z_lo <- if (alpha_lo > 0) stats::qnorm(alpha_lo, lower.tail = FALSE) else NA_real_
  z_hi <- if (alpha_hi > 0) stats::qnorm(alpha_hi, lower.tail = FALSE) else NA_real_
  z_ci <- max(z_lo, z_hi, na.rm = TRUE)   # for the Agresti-Coull added cases

  if (ci_method == "clopper-pearson") {
    ci_lo_app <- if (is.na(z_lo) || p_hat == 0) 0 else
      stats::qbeta(alpha_lo,     x_eff,     m - x_eff + 1)
    ci_hi_app <- if (is.na(z_hi) || p_hat == 1) 1 else
      stats::qbeta(1 - alpha_hi, x_eff + 1, m - x_eff)

  } else if (ci_method == "wald") {
    se_wald   <- sqrt(p_hat * (1 - p_hat) / m)
    ci_lo_app <- if (is.na(z_lo)) 0 else max(p_hat - z_lo * se_wald, 0)
    ci_hi_app <- if (is.na(z_hi)) 1 else min(p_hat + z_hi * se_wald, 1)

  } else {   # agresti-coull
    n_tilde   <- m + z_ci^2
    p_tilde   <- (x_eff + z_ci^2 / 2) / n_tilde
    se_tilde  <- sqrt(p_tilde * (1 - p_tilde) / n_tilde)
    ci_lo_app <- if (is.na(z_lo)) 0 else max(p_tilde - z_lo * se_tilde, 0)
    ci_hi_app <- if (is.na(z_hi)) 1 else min(p_tilde + z_hi * se_tilde, 1)
  }

  rg <- function(p) .rogan_gladen(p, sensitivity, specificity)

  prevalence <- max(0, min(1, rg(p_hat)))
  ci_lower   <- max(0, min(1, rg(ci_lo_app)))
  ci_upper   <- max(0, min(1, rg(ci_hi_app)))

  list(
    statistic   = x_eff,
    p_value     = p_value,
    reject      = reject,
    threshold   = threshold,
    alternative = alternative,
    prevalence  = prevalence,
    ci_lower    = ci_lower,
    ci_upper    = ci_upper,
    ci_method   = ci_method,
    n_total     = n_total,
    n_eff       = n_eff,
    n_eff_adj   = n_eff_adj,
    conf_level  = conf_level,
    sensitivity = sensitivity,
    specificity = specificity,
    icc_used    = icc_used,
    deff        = deff,
    fpc_N       = fpc_N
  )
}
