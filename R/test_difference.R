#' Test whether prevalence differs between two groups
#'
#' @description
#' Analysis function: given the number of positives and the number tested in
#' two groups, tests whether their true prevalence differs. Returns the
#' estimated difference with a confidence interval, the z-statistic, the
#' p-value and whether the difference is significant. It is the analysis
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
#' @param x1 Vector of integer counts of people who tested positive in
#'   group 1, one per site (a single number if group 1 is one site).
#' @param n1 Vector of integer totals tested in group 1, one per site
#'   (same length as \code{x1}).
#' @param x2 Vector of integer counts of people who tested positive in
#'   group 2, one per site (a single number if group 2 is one site).
#' @param n2 Vector of integer totals tested in group 2, one per site
#'   (same length as \code{x2}).
#' @param alternative Direction of the test: \code{"two.sided"} (default),
#'   \code{"greater"} (group 1 higher) or \code{"less"} (group 1 lower).
#' @param sensitivity Diagnostic sensitivity in (0, 1]: the share of people
#'   who truly have the condition that the test detects as positive.
#'   Default 1 (the test never misses a case). Assumed the same in both
#'   groups.
#' @param specificity Diagnostic specificity in (0, 1]: the share of people
#'   who truly do not have the condition that the test correctly reads as
#'   negative. Default 1 (the test never gives a false positive). Assumed the
#'   same in both groups.
#' @param conf_level Confidence level for the interval on the difference, in
#'   \[0.5, 1); default 0.95. The significance level of the test is
#'   1 - \code{conf_level}.
#' @param icc Optional. Intra-cluster correlation (ICC): how alike people in
#'   the same site are, from 0 (not at all) to 1 (identical). If
#'   \code{NULL} (default), it is estimated separately for each group
#'   from its own sites. Otherwise give one number used for both groups, or
#'   two, \code{c(group 1, group 2)}. Each is used exactly as given, however
#'   small: with very large sites even a tiny ICC changes the design effect.
#'   It is also used for a group with only one site (see Details). Set to
#'   \code{0} for no clustering.
#' @param fpc_N Optional. Population size for a finite-population
#'   correction: one positive integer used for each group, or two,
#'   \code{c(group 1, group 2)}. Each must be larger than the number tested
#'   in its group. \code{NULL} (default) means no finite-population
#'   correction.
#' @param ci_method Method for the interval on the difference:
#'   \code{"score"} (default), whose interval leaves out 0 exactly when
#'   the test finds a difference; \code{"newcombe"} (Newcombe's hybrid
#'   score interval, built from a Wilson interval for each group); or
#'   \code{"wald"} (the measured difference plus or minus z standard
#'   errors). \code{"newcombe"} and \code{"wald"} can occasionally leave
#'   out 0 when the test finds no difference, or include 0 when it does,
#'   and \code{"wald"} can be too narrow with small groups (see Details).
#'   It changes only the interval, not the test.
#'
#' @details
#' \strong{Design effect} (returned as \code{deff1} and \code{deff2}). With
#' \code{icc} supplied, it is \eqn{1 + (\bar{n} - 1) ICC} for each group,
#' where \eqn{\bar{n}} is the mean number of people per site. This holds
#' for a group with a single site too: people at the same site tend to be
#' alike, so one site tells you less about the whole group than the same
#' number of people sampled separately. With \code{icc = NULL}, it is
#' estimated for each group as in \code{estimate_prevalence()}: the observed
#' variance of the site prevalences divided by the variance expected if
#' everyone were independent, floored at 1, with the implied ICC kept
#' between 0 and 1. An ICC cannot be estimated from a single site, so a group
#' with one site then gets a design effect of 1; if its people all come from
#' one site, supply \code{icc}. With one person per site the design effect
#' is 1.
#'
#' An ICC estimated from only a few sites is uncertain, and the interval
#' then tends to be too narrow. In simulations with a true ICC of 0.1, 50
#' people per site and prevalences of 0.2 to 0.5, a 95 percent interval
#' contained the true difference about 83 to 85 percent of the time with 2
#' sites per group, 87 to 89 percent with 3 and 90 to 92 percent with 5.
#' Supplying \code{icc} from other studies avoids this.
#'
#' \strong{Effective sample size.} For group \eqn{g}, \eqn{m_g} (returned as
#' \code{n_eff_adj1} and \code{n_eff_adj2}) is the number
#' tested divided by the design effect, then divided by the
#' finite-population factor \eqn{(N - n)/(N - 1)} when \code{fpc_N} is
#' given (\eqn{n} the number tested, \eqn{N} the population size).
#'
#' \strong{Test.} With \eqn{\hat{p}_g} the measured prevalence in group
#' \eqn{g},
#'
#' \deqn{z = \frac{\hat{p}_1 - \hat{p}_2}
#'   {\sqrt{\bar{p}(1-\bar{p})\,(1/m_1 + 1/m_2)}}, \qquad
#'   \bar{p} = \frac{\hat{p}_1 m_1 + \hat{p}_2 m_2}{m_1 + m_2}.}
#'
#' With no clustering and no \code{fpc_N}, this is the same test as R's
#' \code{prop.test(correct = FALSE)}.
#'
#' \strong{Interval.} Built on the measured scale, then divided by
#' \eqn{Se + Sp - 1} to give an interval for the true difference (the rest
#' of the Rogan-Gladen correction cancels in a difference). When
#' \eqn{Se + Sp} is below 1.1 a warning is given: dividing by so small a
#' number makes the difference and its interval very wide.
#'
#' \code{"score"} contains every difference \eqn{\delta} that the same kind
#' of test would not reject. For each \eqn{\delta}, the two prevalences are
#' estimated again under the condition \eqn{p_1 - p_2 = \delta}, giving
#' \eqn{\tilde{p}_1, \tilde{p}_2}, and
#'
#' \deqn{z(\delta) = \frac{\hat{p}_1 - \hat{p}_2 - \delta}
#'   {\sqrt{\tilde{p}_1(1-\tilde{p}_1)/m_1 + \tilde{p}_2(1-\tilde{p}_2)/m_2}}}
#'
#' is compared with the z value for \code{conf_level} (Mee, 1984). At
#' \eqn{\delta = 0} this is exactly the test above, so the interval
#' excludes 0 when, and only when, \code{reject} is \code{TRUE}. Miettinen
#' and Nurminen (1985) multiply the variance by
#' \eqn{(m_1 + m_2)/(m_1 + m_2 - 1)}; that factor is left out here so that
#' the interval and the test always agree. \code{"newcombe"} combines a
#' Wilson interval for each group and can disagree with \code{reject} near
#' the significance cut-off.
#'
#' \code{"wald"} is \eqn{\hat{p}_1 - \hat{p}_2 \pm z \sqrt{\hat{p}_1(1 -
#' \hat{p}_1)/m_1 + \hat{p}_2(1 - \hat{p}_2)/m_2}}. It can disagree with
#' \code{reject}, as the test uses the pooled prevalence \eqn{\bar{p}} and
#' this interval does not. With small groups it can be too narrow: with 20
#' people per group and true prevalences of 0.05 and 0.10, a 95 percent
#' interval contains the true difference only about 90 percent of the time
#' (score 97, Newcombe 99).
#' When each group's measured prevalence is 0 or 1 it has zero width, and a
#' warning is given.
#'
#' The interval matches \code{alternative}: two-sided for
#' \code{"two.sided"}, \eqn{[L, 1]} for \code{"greater"} and \eqn{[-1, U]}
#' for \code{"less"}. Ends outside -1 to 1 are set to -1 or 1.
#'
#' \code{difference} is the measured difference divided by
#' \eqn{Se + Sp - 1}, the value the interval is built around. It can differ
#' from \code{prevalence1 - prevalence2} when a corrected group prevalence is
#' set to 0 or 1.
#'
#' @return A named list with the following fields, in the order returned:
#'   \item{difference}{Estimated true difference, group 1 minus group 2,
#'     between -1 and 1.}
#'   \item{difference_app}{Measured difference,
#'     \code{apparent_prev1 - apparent_prev2}: the difference the test is
#'     done on. Equals \code{difference} when \code{sensitivity} and
#'     \code{specificity} are both 1.}
#'   \item{statistic}{The z-statistic.}
#'   \item{p_value}{p-value for the chosen \code{alternative}.}
#'   \item{reject}{\code{TRUE} if \code{p_value} is below
#'     1 - \code{conf_level}.}
#'   \item{alternative}{Direction of the test, as supplied.}
#'   \item{prevalence1, prevalence2}{Estimated true prevalence in each group
#'     (Rogan-Gladen corrected, between 0 and 1).}
#'   \item{apparent_prev1, apparent_prev2}{Measured prevalence in each
#'     group.}
#'   \item{ci_lower, ci_upper}{Confidence interval for the true difference;
#'     \code{ci_upper} is 1 for \code{"greater"} and \code{ci_lower} is -1
#'     for \code{"less"}.}
#'   \item{ci_method}{Interval method, as supplied.}
#'   \item{n_total1, n_total2}{Number of people tested in each group.}
#'   \item{n_eff1, n_eff2}{Effective sample size of each group before the
#'     finite-population correction: \code{n_total / deff}.}
#'   \item{n_eff_adj1, n_eff_adj2}{Effective sample size each group
#'     contributes to the test and interval, after the finite-population
#'     correction (equals \code{n_eff} when \code{fpc_N} is \code{NULL}).}
#'   \item{conf_level}{Confidence level, as supplied (as a plain number).}
#'   \item{sensitivity, specificity}{As supplied (as plain numbers).}
#'   \item{icc_used1, icc_used2}{ICC used for each group. If
#'     \code{icc = NULL}, it is the value estimated from that group's data
#'     (0 for a group with one site, or one person per site). If \code{icc}
#'     was supplied, it is that group's value.}
#'   \item{deff1, deff2}{Design effect used for each group (1 with no
#'     clustering).}
#'   \item{fpc_N}{\code{fpc_N}, as supplied (as plain numbers), or
#'     \code{NULL}.}
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
#' Newcombe, R. G. (1998) Interval estimation for the difference between
#' independent proportions: comparison of eleven methods. \emph{Statistics
#' in Medicine} \strong{17}(8), 873-890.
#' \doi{10.1002/(SICI)1097-0258(19980430)17:8<873::AID-SIM779>3.0.CO;2-I}
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
#' # Newcombe or Wald interval instead of the default score interval
#' test_difference(x1 = 3, n1 = 80, x2 = 12, n2 = 90, ci_method = "newcombe")
#' test_difference(x1 = 3, n1 = 80, x2 = 12, n2 = 90, ci_method = "wald")
#'
#' # A different ICC for each group
#' test_difference(
#'   x1 = c(2, 5, 3), n1 = c(50, 60, 40),
#'   x2 = c(8, 11, 9), n2 = c(55, 70, 45),
#'   icc = c(0.02, 0.08)
#' )
test_difference <- function(x1, n1, x2, n2,
                            alternative = "two.sided",
                            sensitivity = 1,
                            specificity = 1,
                            conf_level  = 0.95,
                            icc         = NULL,
                            fpc_N       = NULL,
                            ci_method   = "score") {

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

  # ---- validate the counts in each group ----
  # The same checks as estimate_prevalence(), plus a limit of 2^53 people
  # per site, with the argument names of each group.
  check_counts <- function(x, n, xn, nn) {
    if (is.logical(x) || is.logical(n))
      stop("`", xn, "` and `", nn, "` must be numeric, not logical (got class `",
           class(x)[1], "` for ", xn, ", `", class(n)[1], "` for ", nn, "). ",
           "Note: a plain `NA` is logical in R. Remove missing observations ",
           "before calling.")
    if (!is.numeric(x) || !is.numeric(n))
      stop("`", xn, "` and `", nn, "` must be numeric (got class `", class(x)[1],
           "` for ", xn, ", `", class(n)[1], "` for ", nn, ").")
    if (length(x) == 0 || length(n) == 0)
      stop("`", xn, "` and `", nn, "` must each contain at least one value. ",
           "Supply at least one site's count and total.")
    # For each check below, `i` is the position of the first bad element.
    if (!all(is.finite(x))) {
      i <- which(!is.finite(x))[1]
      stop("`", xn, "` contains a missing or infinite value (found ", xn, "[", i,
           "] = ", num(x[i]), "). The number of positives at each site must be a ",
           "finite non-negative integer.")
    }
    if (!all(is.finite(n))) {
      i <- which(!is.finite(n))[1]
      stop("`", nn, "` contains a missing or infinite value (found ", nn, "[", i,
           "] = ", num(n[i]), "). The number tested at each site must be a finite ",
           "positive integer.")
    }
    if (length(x) != length(n))
      stop("`", xn, "` and `", nn, "` must have the same length (got length(", xn,
           ") = ", length(x), ", length(", nn, ") = ", length(n), "). Each value ",
           "of `", xn, "` is the number of positives at one site, and each ",
           "value of `", nn, "` is the number tested at that site.")
    if (any(x < 0)) {
      i <- which(x < 0)[1]
      stop("`", xn, "` must be non-negative (found ", xn, "[", i, "] = ", num(x[i]),
           "). Counts cannot be negative.")
    }
    if (any(n <= 0)) {
      i <- which(n <= 0)[1]
      stop("`", nn, "` must be positive for every site (found ", nn, "[", i,
           "] = ", num(n[i]), "). Each site must have tested at least one person.")
    }
    if (any(x != floor(x))) {
      i <- which(x != floor(x))[1]
      stop("`", xn, "` must contain integers, not decimals (found ", xn, "[", i,
           "] = ", num(x[i]), near_int(x[i], xn), ").")
    }
    if (any(n != floor(n))) {
      i <- which(n != floor(n))[1]
      stop("`", nn, "` must contain integers, not decimals (found ", nn, "[", i,
           "] = ", num(n[i]), near_int(n[i], nn), ").")
    }
    if (any(x > n)) {
      i <- which(x > n)[1]
      stop("`", xn, "` cannot exceed `", nn, "` (found ", xn, "[", i, "] = ", num(x[i]),
           " > ", nn, "[", i, "] = ", num(n[i]), "). The number of positives cannot ",
           "be more than the number tested at that site.")
    }
    # Above 2^53 R cannot store every integer exactly, and no real site tests
    # this many people. (Much larger counts would also overflow the score
    # interval's arithmetic.)
    if (any(n > 2^53)) {
      i <- which(n > 2^53)[1]
      stop("`", nn, "` contains a value larger than 2^53 = 9007199254740992, ",
           "the largest whole number R can store exactly (found ", nn, "[", i,
           "] = ", huge(n[i]), "). No real site tests this many people.")
    }
    invisible(TRUE)
  }
  check_counts(x1, n1, "x1", "n1")
  check_counts(x2, n2, "x2", "n2")

  # ---- validate single-value parameters ----
  # Each must be a single value, except `icc` and `fpc_N`, which can have
  # one value for each group.
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
  if (!is.null(icc) && (!(length(icc) %in% c(1, 2)) || !is.numeric(icc)))
    stop("`icc` must be one or two numbers in [0, 1]: one for both groups, ",
         "or c(group 1, group 2) (got class `", class(icc)[1], "`, length ",
         length(icc), "). To estimate ICC from the data, leave `icc = NULL`.")
  if (!is.character(alternative) || length(alternative) != 1)
    stop("`alternative` must be one of 'two.sided', 'greater', or 'less' ",
         "(got class `", class(alternative)[1], "`, length ",
         length(alternative), ").")
  if (!alternative %in% c("two.sided", "greater", "less"))
    stop("`alternative` must be one of 'two.sided', 'greater', or 'less' ",
         "(got '", alternative, "').")
  if (!is.character(ci_method) || length(ci_method) != 1)
    stop("`ci_method` must be one of 'score', 'newcombe' or 'wald' (got class `",
         class(ci_method)[1], "`, length ", length(ci_method), ").")
  if (!ci_method %in% c("score", "newcombe", "wald"))
    stop("`ci_method` must be one of 'score', 'newcombe' or 'wald' (got '",
         ci_method, "').")

  # Check for NA, NaN and Inf before the range checks below.
  if (!is.finite(sensitivity))
    stop("`sensitivity` must be a single finite number (got ", num(sensitivity), ").")
  if (!is.finite(specificity))
    stop("`specificity` must be a single finite number (got ", num(specificity), ").")
  if (!is.finite(conf_level))
    stop("`conf_level` must be a single finite number (got ", num(conf_level), ").")
  # How a bad `icc` value is shown: "(got 1.5)" for one value, "(found
  # icc[2] = 1.5)" for the first bad one of two.
  icc_got <- function(bad) {
    i <- which(bad)[1]
    if (length(icc) == 1) paste0("(got ", num(icc), ")")
    else paste0("(found icc[", i, "] = ", num(icc[i]), ")")
  }
  if (!is.null(icc) && !all(is.finite(icc)))
    stop("`icc` must be finite ", icc_got(!is.finite(icc)), ". ",
         "To estimate ICC from the data, leave `icc = NULL`.")

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
  if (!is.null(icc) && any(icc < 0 | icc > 1)) {
    bad <- icc < 0 | icc > 1
    stop("`icc` must be in [0, 1] ", icc_got(bad), ". ",
         "`icc` is a correlation and cannot be ",
         if (icc[which(bad)[1]] < 0) "negative. " else "greater than 1. ",
         "To estimate ICC from the data, leave `icc = NULL`.")
  }
  if (conf_level >= 1)
    stop("`conf_level` must be in [0.5, 1) (got ", num(conf_level), "). ",
         "For example, use 0.95 for a 95% confidence interval.")
  if (conf_level < 0.5)
    stop("`conf_level` must be at least 0.5 (got ", num(conf_level), "). ",
         "Below 0.5, when there is no real difference, the test would ",
         "wrongly report one more than half the time.")
  fpc_what <- paste0("`fpc_N` is the population size of each group: one value ",
                     "for both groups, or c(group 1, group 2). Set `fpc_N = NULL` ",
                     "to skip the finite-population correction.")
  if (!is.null(fpc_N)) {
    if (!is.numeric(fpc_N) || !(length(fpc_N) %in% c(1, 2)))
      stop("`fpc_N` must be one or two finite positive integers (got ",
           if (!is.numeric(fpc_N)) paste0("class `", class(fpc_N)[1], "`")
           else paste0("length = ", length(fpc_N)), "). ", fpc_what)
    bad <- !is.finite(fpc_N)
    bad[!bad] <- fpc_N[!bad] < 1 | fpc_N[!bad] != floor(fpc_N[!bad])
    # "(got 100.5)" for one value, "(found fpc_N[2] = 100.5)" for the first
    # bad one of two, as for `icc` above.
    if (any(bad)) {
      i <- which(bad)[1]
      stop("`fpc_N` must be one or two finite positive integers ",
           if (length(fpc_N) == 1) "(got " else paste0("(found fpc_N[", i, "] = "),
           num(fpc_N[i]), near_int(fpc_N[i], "fpc_N"), "). ", fpc_what)
    }
  }

  # From here on these are plain numbers: a name or a 1 x 1 matrix would
  # otherwise leak into the result.
  sensitivity <- as.numeric(sensitivity)
  specificity <- as.numeric(specificity)
  conf_level  <- as.numeric(conf_level)
  if (!is.null(icc)) icc <- as.numeric(icc)
  if (!is.null(fpc_N)) fpc_N <- as.numeric(fpc_N)

  # Warnings are collected and given only once every check below has
  # passed, so a bad input gives just its error.
  warn_later <- character(0)
  # A sum of 1.1 that comes out a hair below it by rounding (0.15 + 0.95)
  # does not warn.
  if (correction < 0.1 - 1e-12)
    warn_later <- c(warn_later, paste0(
      "`sensitivity` + `specificity` = ", num(sensitivity + specificity),
      ", which is very close to 1. Correcting for this much test error ",
      "makes the estimated difference and its confidence interval very ",
      "wide, and small changes in the assumed sensitivity or specificity ",
      "will change them a lot."))

  # ---- design effect for each group ----
  # icc_g is this group's supplied ICC, or NULL to estimate it (as in
  # estimate_prevalence).
  group_deff <- function(x, n, icc_g) {
    n_bar <- mean(n)   # mean number of people per site
    p_hat <- sum(x) / sum(n)
    if (is.null(icc_g)) {
      # A single site (no between-site variation to measure) or every
      # site has one person (no clustering).
      if (length(n) < 2 || n_bar == 1)
        return(list(icc_used = 0, deff = 1))
      var_obs <- stats::var(x / n)
      # Variance expected if everyone were independent, using the overall
      # prevalence of the group.
      var_srs <- mean(p_hat * (1 - p_hat) / n)
      deff    <- if (var_srs > 0) max(var_obs / var_srs, 1) else 1
      icc_u   <- min(max((deff - 1) / (n_bar - 1), 0), 1)
      list(icc_used = icc_u, deff = 1 + (n_bar - 1) * icc_u)
    } else {
      # A supplied ICC is used as given, also with one site: people at
      # one site still resemble each other. With one person per site this
      # gives deff = 1.
      list(icc_used = icc_g, deff = 1 + (n_bar - 1) * icc_g)
    }
  }
  # One supplied value is used for both groups.
  g1 <- group_deff(x1, n1, if (is.null(icc)) NULL else icc[1])
  g2 <- group_deff(x2, n2, if (is.null(icc)) NULL else icc[length(icc)])

  n_tot1 <- sum(n1); n_tot2 <- sum(n2)
  p_hat1 <- sum(x1) / n_tot1
  p_hat2 <- sum(x2) / n_tot2
  n_eff1 <- n_tot1 / g1$deff
  n_eff2 <- n_tot2 / g2$deff

  # ---- finite-population correction for each group ----
  fpc_N1 <- if (is.null(fpc_N)) NULL else fpc_N[1]
  fpc_N2 <- if (is.null(fpc_N)) NULL else fpc_N[length(fpc_N)]
  fpc_check <- function(fN, n_tot, grp) {
    if (is.null(fN)) return(invisible())
    if (fN < n_tot)
      stop("`fpc_N` for group ", grp, " (", num(fN),
           ") is smaller than the number tested in that group (",
           num(n_tot), "): you cannot test more people than ",
           "there are in the population.")
    if (fN == n_tot)
      stop("`fpc_N` for group ", grp, " (", num(fN),
           ") equals the number tested in that group: the whole group was ",
           "tested (a census), so its prevalence is known exactly, and this ",
           "test needs a sample from each group. Set `fpc_N = NULL` if no ",
           "finite-population correction is needed.")
  }
  fpc_check(fpc_N1, n_tot1, "1")
  fpc_check(fpc_N2, n_tot2, "2")
  # Finite-population correction factor, f = (N - n_total) / (N - 1), using
  # the number of people tested. 1 when fpc_N is NULL.
  f1 <- if (is.null(fpc_N1)) 1 else (fpc_N1 - n_tot1) / (fpc_N1 - 1)
  f2 <- if (is.null(fpc_N2)) 1 else (fpc_N2 - n_tot2) / (fpc_N2 - 1)
  m1 <- n_eff1 / f1
  m2 <- n_eff2 / f2

  # ---- pooled z-test on the measured prevalences ----
  d_app   <- p_hat1 - p_hat2
  p_pool  <- (p_hat1 * m1 + p_hat2 * m2) / (m1 + m2)
  # 1 - p_pool worked out directly: when one group's effective size is
  # billions of times the other's, p_pool can round to exactly 1 even
  # though the other group's prevalence is below 1.
  q_pool  <- ((1 - p_hat1) * m1 + (1 - p_hat2) * m2) / (m1 + m2)
  se_null <- sqrt(p_pool * q_pool * (1 / m1 + 1 / m2))
  if (se_null == 0)
    stop("Both groups have a measured prevalence of ", num(p_pool), ", so there is ",
         "no variation and the test statistic cannot be calculated.")
  # The Wald standard error is 0 only when each group's measured
  # prevalence is 0 or 1 (and they differ, or the stop above would have
  # been reached).
  se_wald <- sqrt(p_hat1 * (1 - p_hat1) / m1 + p_hat2 * (1 - p_hat2) / m2)
  if (ci_method == "wald" && se_wald == 0)
    warn_later <- c(warn_later, paste0(
      "Each group's measured prevalence is 0 or 1 (", num(p_hat1), " and ",
      num(p_hat2), "), so the Wald interval has zero width. Use ",
      "`ci_method = \"score\"` for an interval that allows for the ",
      "uncertainty."))
  # Every check has passed: give the warnings collected above.
  for (w in warn_later) warning(w, call. = FALSE)
  z_stat <- d_app / se_null

  p_value <- switch(alternative,
    greater   = stats::pnorm(z_stat, lower.tail = FALSE),
    less      = stats::pnorm(z_stat, lower.tail = TRUE),
    two.sided = 2 * stats::pnorm(abs(z_stat), lower.tail = FALSE))
  alpha  <- 1 - conf_level
  reject <- p_value < alpha

  # ---- interval on the difference, matched to `alternative` ----
  # "greater" -> [L, 1], "less" -> [-1, U], "two.sided" -> [L, U]. A
  # one-sided bound uses z at 1 - alpha.
  alpha_lo <- switch(alternative, two.sided = alpha / 2, greater = alpha, less = 0)
  alpha_hi <- switch(alternative, two.sided = alpha / 2, greater = 0,     less = alpha)
  # Written this way, not as qnorm(1 - alpha), so it stays exact when
  # `conf_level` is very close to 1.
  z_lo <- if (alpha_lo > 0) stats::qnorm(alpha_lo, lower.tail = FALSE) else NA_real_
  z_hi <- if (alpha_hi > 0) stats::qnorm(alpha_hi, lower.tail = FALSE) else NA_real_
  z_ci <- max(z_lo, z_hi, na.rm = TRUE)

  if (ci_method == "score") {
    lo_app <- if (is.na(z_lo)) -Inf else
      .score_bound(p_hat1, m1, p_hat2, m2, z_lo, "lower")
    hi_app <- if (is.na(z_hi))  Inf else
      .score_bound(p_hat1, m1, p_hat2, m2, z_hi, "upper")
  } else if (ci_method == "wald") {
    lo_app <- if (is.na(z_lo)) -Inf else d_app - z_lo * se_wald
    hi_app <- if (is.na(z_hi))  Inf else d_app + z_hi * se_wald
  } else {
    w1 <- .wilson_ci(p_hat1, m1, z_ci)
    w2 <- .wilson_ci(p_hat2, m2, z_ci)
    lo_app <- if (is.na(z_lo)) -Inf else
      d_app - sqrt((p_hat1 - w1$lower)^2 + (w2$upper - p_hat2)^2)
    hi_app <- if (is.na(z_hi))  Inf else
      d_app + sqrt((w1$upper - p_hat1)^2 + (p_hat2 - w2$lower)^2)
  }

  # Back to the true scale: divide by the correction, then move any end
  # below -1 up to -1, and any end above 1 down to 1.
  ci_lower <- if (is.infinite(lo_app)) -1 else max(-1, min(1, lo_app / correction))
  ci_upper <- if (is.infinite(hi_app))  1 else max(-1, min(1, hi_app / correction))

  # ---- estimates on the true scale ----
  # `difference` is the measured difference divided by the correction, the
  # value the interval is built around (prevalence1 - prevalence2 can differ
  # from it when a group prevalence is set to 0 or 1).
  rg <- function(p) .rogan_gladen(p, sensitivity, specificity)
  prevalence1 <- max(0, min(1, rg(p_hat1)))
  prevalence2 <- max(0, min(1, rg(p_hat2)))
  difference  <- max(-1, min(1, d_app / correction))

  list(
    difference     = difference,
    difference_app = d_app,
    statistic      = z_stat,
    p_value        = p_value,
    reject         = reject,
    alternative    = alternative,
    prevalence1    = prevalence1,
    prevalence2    = prevalence2,
    apparent_prev1 = p_hat1,
    apparent_prev2 = p_hat2,
    ci_lower       = ci_lower,
    ci_upper       = ci_upper,
    ci_method      = ci_method,
    n_total1       = n_tot1,
    n_total2       = n_tot2,
    n_eff1         = n_eff1,
    n_eff2         = n_eff2,
    n_eff_adj1     = m1,
    n_eff_adj2     = m2,
    conf_level     = conf_level,
    sensitivity    = sensitivity,
    specificity    = specificity,
    icc_used1      = g1$icc_used,
    icc_used2      = g2$icc_used,
    deff1          = g1$deff,
    deff2          = g2$deff,
    fpc_N          = fpc_N
  )
}


# ---- score interval helpers ----

# Prevalences estimated again under the condition p1 - p2 = delta, by
# maximum likelihood (the closed form of Farrington and Manning, 1990). m1
# and m2 can be fractional (effective sizes).
#
# The same answer can be found four ways: as given, with the groups
# swapped, with positives and negatives swapped, or both. In extreme data
# (a prevalence within about 1e-6 of 0 or 1, or one group's effective size
# billions of times the other's) one of them can lose precision or land on
# the wrong root, so all four are tried. Any clearly less likely than the
# best (a wrong root) is dropped; of the rest, the one that solved for the
# value closest to 0 is kept, as arithmetic is most precise there.
.restricted_prev <- function(p1, m1, p2, m2, delta) {
  tries <- list(
    .restricted_prev_one(p1, m1, p2, m2, delta),
    rev(.restricted_prev_one(p2, m2, p1, m1, -delta)),
    1 - .restricted_prev_one(1 - p1, m1, 1 - p2, m2, -delta),
    1 - rev(.restricted_prev_one(1 - p2, m2, 1 - p1, m1, delta)))
  ll     <- vapply(tries, function(t) .restricted_loglik(p1, m1, p2, m2, t), 0)
  solved <- c(tries[[1]][1], tries[[2]][2], 1 - tries[[3]][1], 1 - tries[[4]][2])
  ok     <- which(ll >= max(ll) - 1e-10 * (1 + abs(max(ll))))
  tries[[ok[which.min(solved[ok])]]]
}

# Log-likelihood of prevalences t = c(t1, t2) for the measured prevalences
# p1, p2 with effective sizes m1, m2 (a 0 count contributes 0).
.restricted_loglik <- function(p1, m1, p2, m2, t) {
  term <- function(w, x) if (w == 0) 0 else w * log(x)
  m1 * (term(p1, t[1]) + term(1 - p1, 1 - t[1])) +
    m2 * (term(p2, t[2]) + term(1 - p2, 1 - t[2]))
}

# Solves the Farrington-Manning cubic for t1, the prevalence of group 1,
# once; t2 = t1 - delta.
.restricted_prev_one <- function(p1, m1, p2, m2, delta) {
  th <- m2 / m1
  a  <- 1 + th
  b  <- -(1 + th + p1 + th * p2 + delta * (th + 2))
  cc <- delta^2 + delta * (2 * p1 + th + 1) + p1 + th * p2
  d  <- -p1 * delta * (1 + delta)
  v  <- b^3 / (3 * a)^3 - b * cc / (6 * a^2) + d / (2 * a)
  # (Published versions give u the sign of v; that makes no difference
  # to the result.)
  u  <- sqrt(max(b^2 / (3 * a)^2 - cc / (3 * a), 0))
  # Guard against rounding pushing v / u^3 a hair outside [-1, 1].
  t1 <- if (u == 0) -b / (3 * a) else
    2 * u * cos((pi + acos(max(-1, min(1, v / u^3)))) / 3) - b / (3 * a)
  # The closed form loses precision when a prevalence is tiny (around 1e-6
  # or less), so polish it with Newton steps on the same cubic, kept inside
  # the range where both prevalences are between 0 and 1.
  lo <- max(0, delta)
  hi <- min(1, 1 + delta)
  t1 <- min(hi, max(lo, t1))
  for (i in 1:50) {
    slope <- (3 * a * t1 + 2 * b) * t1 + cc
    if (slope == 0) break
    t_new <- t1 - (((a * t1 + b) * t1 + cc) * t1 + d) / slope
    t_new <- min(hi, max(lo, t_new))
    done  <- abs(t_new - t1) <= 1e-15 * abs(t1)
    t1    <- t_new
    if (done) break
  }
  c(t1, min(1, max(0, t1 - delta)))
}

# Score statistic for the difference delta (Mee, 1984). At delta = 0 it is
# the pooled z-statistic of the test.
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
# search is precise to 1e-12, or to a billionth of the standard error when
# that is smaller, so it stays accurate when huge samples make the
# interval very narrow.
.score_bound <- function(p1, m1, p2, m2, z, side) {
  d <- p1 - p2
  e <- 1e-10
  f <- function(delta) .score_stat(p1, m1, p2, m2, delta)
  se  <- sqrt(p1 * (1 - p1) / m1 + p2 * (1 - p2) / m2)
  tol <- if (se > 0) min(1e-12, se * 1e-9) else 1e-12
  if (side == "lower") {
    if (d <= -1 + e || f(-1 + e) <= z) return(-1)
    stats::uniroot(function(x) f(x) - z, c(-1 + e, d), tol = tol)$root
  } else {
    if (d >= 1 - e || f(1 - e) >= -z) return(1)
    stats::uniroot(function(x) f(x) + z, c(d, 1 - e), tol = tol)$root
  }
}
