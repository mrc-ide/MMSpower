#' Calculate sample size to detect a difference in prevalence between two groups
#'
#' @description
#' Design function: given the expected true prevalence in two groups, returns
#' the sample size needed in each group to detect that difference with the
#' target power.
#'
#' Each group can have its own number of sites, people per site, ICC and
#' population size: give one value for both groups, or two,
#' \code{c(group 1, group 2)}. The function solves for one unknown, the
#' same in both groups:
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
#' @param n_sites Optional. Number of sites in each group: one positive
#'   integer for both groups, or two, \code{c(group 1, group 2)}. The
#'   function solves for the number of people per site, the same in both
#'   groups. Cannot be used together with \code{n_per_site}.
#' @param n_per_site Optional. Number of people per site: one positive
#'   integer for both groups, or two, \code{c(group 1, group 2)}. The
#'   function solves for the number of sites, the same in both groups.
#'   Cannot be used together with \code{n_sites}. If one site in each group
#'   already has more people than needed, the function stops with an error;
#'   use \code{n_sites = 1} to find how many people that one site needs.
#' @param icc Intra-cluster correlation (ICC): how alike people in the same
#'   site are, from 0 (not at all) to 1 (identical). One number for both
#'   groups, or two, \code{c(group 1, group 2)}. Default 0 (no clustering).
#'   An \code{icc} above 0 needs \code{n_sites} or \code{n_per_site}, and the
#'   value given is then used exactly, however small. Without sites, a tiny
#'   value below about 1.5e-8, such as one left over from rounding in an
#'   earlier calculation, is treated as 0.
#' @param fpc_N Optional. Population size of each group, for a
#'   finite-population correction: one positive integer for both groups, or
#'   two, \code{c(group 1, group 2)}. \code{NULL} (default) means no
#'   finite-population correction.
#'
#' @details
#' \strong{Sample size.} Write \eqn{m_1}, \eqn{m_2} for the effective sizes
#' of the two groups (below). The test reaches the target power when
#'
#' \deqn{z_\alpha \sqrt{\bar{p}(1-\bar{p})\,(1/m_1 + 1/m_2)} +
#'   z_\beta \sqrt{p_{1,app}(1-p_{1,app})/m_1 + p_{2,app}(1-p_{2,app})/m_2}
#'   = |p_{1,app} - p_{2,app}|}
#'
#' where \eqn{p_{g,app} = p_g \cdot Se + (1-p_g)(1-Sp)} is the prevalence the
#' test would measure in group \eqn{g}, \eqn{\bar{p} = (p_{1,app} m_1 +
#' p_{2,app} m_2)/(m_1 + m_2)}, \eqn{z_\beta} = \code{qnorm(power)} (0.84 for
#' 80 percent power), and \eqn{z_\alpha} = \code{qnorm(1 - alpha/2)} for a
#' two-sided test (1.96 at 95 percent) or \code{qnorm(1 - alpha)} for a
#' one-sided test (1.64 at 95 percent), with \eqn{\alpha = 1 -}
#' \code{conf_level}. The function finds the smallest value of the unknown
#' (people per group, sites, or people per site) for which the left side is
#' no larger than the right, by a numerical search, and rounds it up to a
#' whole number.
#'
#' With no clustering and two groups of \eqn{n} people this gives
#'
#' \deqn{n_0 = \frac{\left(z_\alpha \sqrt{2\,\bar{p}(1-\bar{p})} +
#'   z_\beta \sqrt{p_{1,app}(1-p_{1,app}) + p_{2,app}(1-p_{2,app})}\right)^2}
#'   {(p_{1,app} - p_{2,app})^2}}
#'
#' the same formula as R's \code{power.prop.test()}.
#'
#' Because \eqn{p_{1,app} - p_{2,app} = (p_1 - p_2)(Se + Sp - 1)}, an
#' imperfect test shrinks the difference the test can see and increases the
#' sample size. When \eqn{Se = Sp = 1} the measured and true prevalences are
#' the same.
#'
#' \strong{Effective size.} For group \eqn{g} with \eqn{n_g} people,
#' \eqn{m_g = n_g / D_g}, where \eqn{D_g = 1 + (\bar{n}_g - 1)\,ICC_g} is the
#' design effect and \eqn{\bar{n}_g} the number of people per site. With
#' \code{fpc_N} (the population size \eqn{N_g}), sampling without
#' replacement reduces the variance by the factor \eqn{(N_g - n_g)/(N_g - 1)},
#' so \eqn{m_g} is divided by it as well (Cochran, 1977).
#'
#' \strong{Highest power.} With \code{n_sites} fixed and no \code{fpc_N},
#' adding people to each site cannot raise the power past a limit. People in
#' the same site tend to resemble each other, so after a point each extra
#' person adds almost nothing new. Each site is worth at most \eqn{1 / ICC}
#' independent people, so a group can never count as more than
#' \eqn{n_{sites} / ICC} independent people. For example, with ICC 0.05 each
#' site is worth at most 20 independent people, so 8 sites per group give at
#' most 160 independent people. For 10 versus 20 percent, 160 people per
#' group give 70.9 percent power, so 80 percent can never be reached. When
#' the target cannot be reached, the function stops and reports the highest
#' power the sites can give.
#'
#' With \code{fpc_N}, a group whose whole population is tested has no
#' sampling uncertainty left, so with the same design in both groups the
#' target is always reachable. When the groups differ, one group can reach
#' its whole population before the target is reached; the function then
#' stops and reports the highest power at that point. If there are too few
#' sites for the ICC to reach the target power on their own, the
#' correction can still make it reachable, but often only by testing a large
#' share of the population; a warning says so.
#'
#' \code{fpc_N} counts people, not sites. If a study samples nearly all of
#' the sites in a group, the returned sample size is larger than necessary,
#' because the correction does not remove the variation between sites.
#'
#' \strong{Rounding.} The people needed in each group, and the number of
#' sites or people per site, are rounded up. Each design effect (returned as
#' \code{deff1}, \code{deff2}) is then worked out again from the rounded
#' number of people per site, so that
#' \eqn{D_g = 1 + (n_{per site, g} - 1) ICC_g} describes the design reported.
#'
#' @return A named list with the following fields, in the order returned:
#'   \item{n_per_group1, n_per_group2}{Number of people needed in each
#'     group. They are equal when \code{n_sites}, \code{n_per_site},
#'     \code{icc} and \code{fpc_N} each have one value.}
#'   \item{n}{Total number of people across both groups,
#'     \code{n_per_group1 + n_per_group2}. With sites, every site in a group
#'     gets the same whole number of people, rounded up, so the study
#'     collects \code{n_sites1 * n_per_site1 + n_sites2 * n_per_site2}
#'     people, which can be more than \code{n}. Budget for that number. The
#'     extra is usually small, but can almost double the study when sites
#'     have only 1 or 2 people each: for example, a group needing 199 people
#'     across 198 sites gets 2 per site, so 396 are collected in that group
#'     instead of 199.}
#'   \item{n_eff1, n_eff2}{The number of people each group would need with
#'     no clustering and no finite-population correction, keeping the
#'     groups in the same ratio as the design (for example, twice as many
#'     people in group 2 with \code{n_sites = c(10, 20)}).}
#'   \item{n_sites1, n_sites2}{Sites in each group: solved for when
#'     \code{n_per_site} is given, as given when \code{n_sites} is given. If
#'     the study needs fewer people in a group than its sites, fewer sites are
#'     used, with 1 person at each, and a warning says so. \code{NULL} with
#'     no clustering.}
#'   \item{n_per_site1, n_per_site2}{People per site in each group: solved
#'     for when \code{n_sites} is given, as given when \code{n_per_site} is
#'     given. \code{NULL} with no clustering.}
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
#'   \item{icc}{ICC, as supplied (one or two values).}
#'   \item{deff1, deff2}{Design effect of each group: how much clustering
#'     inflates the number of people needed (1 with no clustering). With
#'     sites, it is worked out for the people actually collected,
#'     \code{n_sites * n_per_site} in that group.}
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
#'
#' # Twice as many sites in group 2, and more clustering there
#' design_difference(0.10, 0.20, n_sites = c(20, 40), icc = c(0.02, 0.05))
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
  if (!(length(icc) %in% c(1, 2)) || !is.numeric(icc))
    stop("`icc` must be one or two numbers in [0, 1]: one for both groups, ",
         "or c(group 1, group 2) (got class `", class(icc)[1], "`, length ",
         length(icc), ").")
  if (!is.character(alternative) || length(alternative) != 1)
    stop("`alternative` must be one of 'two.sided', 'greater', or 'less' ",
         "(got class `", class(alternative)[1], "`, length ",
         length(alternative), ").")
  if (!alternative %in% c("two.sided", "greater", "less"))
    stop("`alternative` must be one of 'two.sided', 'greater', or 'less' ",
         "(got '", alternative, "').")

  # Numbers in messages are printed in full (100000, not 1e+05).
  num <- function(v) format(v, scientific = FALSE, digits = 15)
  # ...except numbers far beyond any real study, shown as 1e+200.
  huge <- function(v) if (v < 1e16) num(v) else format(v, digits = 15)
  # How a bad value of a one-or-two-value argument is shown: "(got 1.5)"
  # for one value, "(found icc[2] = 1.5)" for the first bad one of two.
  got <- function(v, nm, bad) {
    i <- which(bad)[1]
    if (length(v) == 1) paste0("(got ", num(v), ")")
    else paste0("(found ", nm, "[", i, "] = ", num(v[i]), ")")
  }

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
  if (!all(is.finite(icc)))
    stop("`icc` must be finite ", got(icc, "icc", !is.finite(icc)), ".")

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
         # 12 digits, so a sum that is 1 up to rounding (within 1e-12) reads 1
         format(sensitivity + specificity, digits = 12, scientific = FALSE), "). ",
         "The difference the test can see is the true difference times ",
         "`sensitivity` + `specificity` - 1, so with a sum of 1 or less the ",
         "test cannot see the difference, or sees it reversed.")
  if (any(icc < 0 | icc > 1)) {
    bad <- icc < 0 | icc > 1
    stop("`icc` must be in [0, 1] ", got(icc, "icc", bad), ". ",
         "`icc` is a correlation and cannot be ",
         if (icc[which(bad)[1]] < 0) "negative." else "greater than 1.")
  }
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

  # n_sites, n_per_site and fpc_N: one or two positive integers each.
  check_counts <- function(v, nm, what, tail = "") {
    if (is.null(v)) return(invisible())
    if (!is.numeric(v) || !(length(v) %in% c(1, 2)))
      stop("`", nm, "` must be one or two finite positive integers (got ",
           if (!is.numeric(v)) paste0("class `", class(v)[1], "`")
           else paste0("length = ", length(v)), "). ", what, tail)
    bad <- !is.finite(v)
    bad[!bad] <- v[!bad] < 1 | v[!bad] != floor(v[!bad])
    if (any(bad))
      stop("`", nm, "` must be one or two finite positive integers ",
           got(v, nm, bad), ". ", what, tail)
  }
  check_counts(n_sites, "n_sites",
               "`n_sites` is the number of sites in each group: one value for both groups, or c(group 1, group 2).")
  check_counts(n_per_site, "n_per_site",
               "`n_per_site` is the number of people sampled per site: one value for both groups, or c(group 1, group 2).")
  check_counts(fpc_N, "fpc_N",
               "`fpc_N` is the population size of each group: one value for both groups, or c(group 1, group 2).",
               " Set `fpc_N = NULL` to skip the finite-population correction.")
  # How the i-th of one or two values is named in a message.
  lab <- function(v, nm, i) if (length(v) == 1) paste0("`", nm, "`") else paste0("`", nm, "[", i, "]`")
  for (i in seq_along(fpc_N))
    if (fpc_N[i] > 2^53)
      stop(lab(fpc_N, "fpc_N", i), " (", huge(fpc_N[i]), ") is larger than 2^53 = ",
           "9007199254740992, the largest whole number R can store exactly. No ",
           "real population is this large: set `fpc_N = NULL`.")
  for (i in seq_along(n_sites))
    if (n_sites[i] > 2^53)
      stop(lab(n_sites, "n_sites", i), " (", huge(n_sites[i]), ") is larger than 2^53 = ",
           "9007199254740992, the largest whole number R can store exactly. No ",
           "real study has this many sites.")
  for (i in seq_along(n_per_site))
    if (n_per_site[i] > 2^53)
      stop(lab(n_per_site, "n_per_site", i), " (", huge(n_per_site[i]), ") is larger ",
           "than 2^53 = 9007199254740992, the largest whole number R can store ",
           "exactly. No real site has this many people.")

  # From here on every number is a plain number: a 1 x 1 matrix, a name or a
  # class would otherwise leak into the arithmetic and the result.
  prevalence1 <- as.numeric(prevalence1); prevalence2 <- as.numeric(prevalence2)
  power       <- as.numeric(power);       conf_level  <- as.numeric(conf_level)
  sensitivity <- as.numeric(sensitivity); specificity <- as.numeric(specificity)
  icc         <- as.numeric(icc)
  if (!is.null(n_sites))    n_sites    <- as.numeric(n_sites)
  if (!is.null(n_per_site)) n_per_site <- as.numeric(n_per_site)
  if (!is.null(fpc_N))      fpc_N      <- as.numeric(fpc_N)

  # One value is used for both groups: from here on each is a pair.
  icc_g <- rep_len(icc, 2)
  k_g   <- if (is.null(n_sites))    NULL else rep_len(n_sites, 2)
  s_g   <- if (is.null(n_per_site)) NULL else rep_len(n_per_site, 2)
  N_g   <- if (is.null(fpc_N))      NULL else rep_len(fpc_N, 2)

  # The number of sites, or people per site, cannot exceed the population.
  # (With one value of each, the message names no group.)
  grp_of <- function(i, a, b) if (length(a) == 1 && length(b) == 1) "" else paste0(" for group ", i)
  for (i in 1:2) {
    if (!is.null(N_g) && !is.null(k_g) && k_g[i] > N_g[i])
      stop("`n_sites`", grp_of(i, n_sites, fpc_N), " (", num(k_g[i]),
           ") exceeds `fpc_N`", grp_of(i, n_sites, fpc_N), " (", num(N_g[i]),
           "): you cannot have more sites than people in the population.")
    if (!is.null(N_g) && !is.null(s_g) && s_g[i] > N_g[i])
      stop("`n_per_site`", grp_of(i, n_per_site, fpc_N), " (", num(s_g[i]),
           ") exceeds `fpc_N`", grp_of(i, n_per_site, fpc_N), " (", num(N_g[i]),
           "): a site cannot have more people than there are in the whole population.")
  }

  # Without sites there is no clustering, so a tiny icc there (e.g. 1e-12, a
  # rounding leftover) counts as 0. With sites, any icc above 0 is used: with
  # very large sites even a tiny icc changes the design effect.
  no_sites <- is.null(n_sites) && is.null(n_per_site)
  if (no_sites && any(icc >= sqrt(.Machine$double.eps)))
    stop("`icc` is greater than 0 but no sites were given (`n_sites` and ",
         "`n_per_site` are both `NULL`), so the design effect cannot be ",
         "computed. Supply `n_sites` (fix the number of sites) or `n_per_site` ",
         "(fix the number of people per site), or set `icc = 0` for a design ",
         "without clustering.")
  if (no_sites) icc_g <- c(0, 0)

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

  # ---- the requirement ----
  # For effective sizes m = c(m1, m2) (Inf for a group whose whole
  # population is tested), the ratio of the left side of the equation in
  # Details to the right side. The target power is reached when it is at
  # most 1. With m1 = m2 = m it is sqrt(n_base_cont / m).
  p_app <- c(p1_app, p2_app)
  ratio <- function(m) {
    if (any(m == 0)) return(Inf)
    inv <- ifelse(is.infinite(m), 0, 1 / m)
    if (all(inv == 0)) return(0)
    if (any(inv == 0)) {
      # A fully tested group carries all the weight in the pooled prevalence.
      pb <- p_app[inv == 0][1]; qb <- 1 - pb
    } else {
      pb <- sum(p_app * m) / sum(m)
      # 1 - pb worked out directly, so it does not round to 0 when pb is
      # within rounding of 1.
      qb <- sum((1 - p_app) * m) / sum(m)
    }
    v0 <- pb * qb * sum(inv)
    # (A prevalence of 0 or 1 adds nothing, even when 1 / m is huge.)
    pq <- p_app * (1 - p_app)
    v1 <- sum(pq[pq > 0] * inv[pq > 0])
    # (z is 0 for a one-sided test at conf_level 0.5, or power 0.5: then
    # its term is 0, even when the variance is huge.)
    t0 <- if (z_a == 0) 0 else z_a * sqrt(v0)
    t1 <- if (z_b == 0 || v1 == 0) 0 else z_b * sqrt(v1)
    (t0 + t1) / delta
  }

  # People in each group, and their effective sizes, for a value x of the
  # unknown. `how` is "none" (x people per group), "per_site" (x sites of
  # n_per_site people), "sites" (n_sites sites of x people) or "one_each"
  # (a share x of the sites, 1 person at each). `cluster` and `fpc` FALSE
  # leave out the design effect or the finite-population correction.
  people <- function(x, how) {
    switch(how,
      none     = c(x, x),
      per_site = x * s_g,
      sites    = k_g * x,
      one_each = x * k_g)
  }
  eff <- function(x, how, cluster = TRUE, fpc = TRUE) {
    n <- people(x, how)
    d <- if (!cluster || how %in% c("none", "one_each")) c(1, 1)
         else if (how == "per_site") 1 + (s_g - 1) * icc_g
         else 1 + (x - 1) * icc_g
    m <- n / d
    if (fpc && !is.null(N_g)) {
      full <- n >= N_g
      m[full]  <- Inf
      m[!full] <- m[!full] * (N_g[!full] - 1) / (N_g[!full] - n[!full])
    }
    m
  }
  # The largest x before the first group reaches its whole population (Inf
  # without fpc_N).
  x_full <- function(how) {
    if (is.null(N_g)) return(Inf)
    min(N_g / switch(how, none = c(1, 1), per_site = s_g, sites = k_g, one_each = k_g))
  }

  # The smallest x in (lo, hi] that meets the requirement, found by halving
  # the interval until it is as small as R can store. ratio() falls as x
  # grows. Returns NA when even hi does not meet it.
  solve_x <- function(how, lo, hi, cluster = TRUE, fpc = TRUE) {
    r <- function(x) ratio(eff(x, how, cluster, fpc))
    if (is.infinite(hi)) {
      # Double until the requirement is met. (Past the largest number R can
      # store this gives Inf, which the halving below keeps, and the
      # trillion check then stops.)
      hi <- max(1, lo * 2)
      while (is.finite(hi) && r(hi) > 1) { lo <- hi; hi <- hi * 2 }
    } else if (r(hi) > 1) {
      return(NA_real_)
    }
    repeat {
      mid <- lo + (hi - lo) / 2
      if (mid <= lo || mid >= hi) break
      if (r(mid) <= 1) hi <- mid else lo <- mid
    }
    hi
  }

  # Rounds up to whole people. Shrinking by one part in a trillion first
  # stops rounding error just above a whole number (300.00000000000006)
  # adding a whole person; it could only remove a real person above a
  # trillion people, which is stopped below.
  shrink <- 1 - 1e-12
  up     <- function(v) max(1, ceiling(v * shrink))

  # Describing the design in messages: "9 sites per group" for one value,
  # "10 sites in group 1 and 20 in group 2" for two.
  same <- function(v) is.null(v) || v[1] == v[2]
  sites_text <- function(k)
    if (same(k)) paste0(num(k[1]), if (k[1] == 1) " site" else " sites", " per group")
    else paste0(num(k[1]), if (k[1] == 1) " site" else " sites", " in group 1 and ",
                num(k[2]), " in group 2")
  # Sites and ICC together: "8 sites per group and ICC = 0.05" when both
  # have one value; otherwise "6 sites in group 1 and 8 in group 2, with
  # ICC = 0.05 in group 1 and 0.1 in group 2" (`join` is "with" or "and").
  design_text <- function(join) {
    icc_part <- if (same(icc_g)) num(icc_g[1])
                else paste0(num(icc_g[1]), " in group 1 and ", num(icc_g[2]), " in group 2")
    if (same(k_g) && same(icc_g)) paste0(sites_text(k_g), " and ICC = ", icc_part)
    else paste0(sites_text(k_g), ", ", join, " ICC = ", icc_part)
  }
  # The highest power with effective sizes m (Inf for a fully tested group).
  power_at <- function(m) {
    inv <- ifelse(is.infinite(m), 0, 1 / m)
    if (any(inv == 0)) { pb <- p_app[inv == 0][1]; qb <- 1 - pb }
    else { pb <- sum(p_app * m) / sum(m); qb <- sum((1 - p_app) * m) / sum(m) }
    pq <- p_app * (1 - p_app)
    v1 <- sum(pq[pq > 0] * inv[pq > 0])
    # With prevalences of 0 and 1 the measured difference never varies, so
    # the test always or never finds it.
    if (v1 == 0) return(as.numeric(delta > z_a * sqrt(pb * qb * sum(inv))))
    stats::pnorm((delta - z_a * sqrt(pb * qb * sum(inv))) / sqrt(v1))
  }
  max_power_text <- function(max_power) {
    shown <- format_apart(100 * power, 100 * max_power, digits = 7)
    # A maximum like 2e-16 percent is shown as "less than 0.0001".
    if (100 * max_power < 1e-4) shown[2] <- "less than 0.0001"
    shown
  }
  # The group that reaches its whole population first.
  first_full <- function(how) {
    per <- N_g / switch(how, none = c(1, 1), per_site = s_g, sites = k_g)
    which.min(per)
  }
  stop_full <- function(how, x_cap, fix) {
    g     <- first_full(how)
    shown <- max_power_text(power_at(eff(x_cap, how)))
    stop("Target power of ", shown[1], " percent cannot be reached: group ", g,
         " has only ", num(N_g[g]),
         if (N_g[g] == 1) " person (`fpc_N`), and testing that person"
         else " people (`fpc_N`), and testing all of them",
         " gives at most ", shown[2],
         " percent power. ", fix)
  }

  # ---- solve for the unknown ----
  # TRUE when the sites alone could not reach the target and only the
  # finite-population correction makes it reachable (warned about below).
  fpc_rescue <- FALSE
  # TRUE when more sites were given than people needed, so some are dropped.
  sites_cut  <- FALSE

  if (no_sites) {
    how <- "none"
    x   <- solve_x(how, 0, x_full(how))
    if (is.na(x))
      stop_full(how, x_full(how), "Check `fpc_N`, or lower the target power.")
    n_pg   <- people(x, how)
    x0     <- solve_x(how, 0, Inf, cluster = FALSE, fpc = FALSE)
    n_eff  <- vapply(people(x0, how), up, 0)
    k_out  <- NULL
    s_out  <- NULL
    deff   <- c(1, 1)

  } else if (!is.null(s_g)) {
    how <- "per_site"
    x   <- solve_x(how, 0, x_full(how))
    if (is.na(x))
      stop_full(how, x_full(how), "Use a different `n_per_site`, check `fpc_N`, or lower the target power.")
    n_pg  <- people(x, how)
    x0    <- solve_x(how, 0, Inf, cluster = FALSE, fpc = FALSE)
    n_eff <- vapply(people(x0, how), up, 0)
    deff  <- 1 + (s_g - 1) * icc_g
    s_out <- s_g

  } else {
    how <- "sites"
    # Without fpc_N, is the target reachable at all with these sites? Each
    # site counts as at most 1 / ICC independent people.
    m_inf <- ifelse(icc_g == 0, Inf, k_g / icc_g)
    unreachable <- ratio(m_inf) >= 1
    if (is.null(N_g) && unreachable) {
      shown <- max_power_text(power_at(m_inf))
      stop("Target power of ", shown[1], " percent cannot be ",
           "reached with ", design_text("with"),
           ". The highest power ",
           if (same(k_g) && k_g[1] == 1) "this site" else "these sites",
           " can give is ", shown[2], " percent. ",
           "Increase `n_sites`, lower the ICC, or lower the target power.")
    }
    fpc_rescue <- !is.null(N_g) && unreachable
    # One person per site has no clustering. If that already meets the
    # target, there are more sites than needed: use 1 person at each of
    # just enough of them (the same share of each group's sites).
    if (ratio(eff(1, how)) <= 1) {
      t      <- solve_x("one_each", 0, 1)
      n_pg   <- people(t, "one_each")
      k_used <- pmin(k_g, vapply(n_pg, up, 0))
      sites_cut <- any(k_used < k_g)
      x <- 1   # every site used has 1 person
    } else {
      x <- solve_x(how, 1, x_full(how))
      if (is.na(x))
        stop_full(how, x_full(how), "Use different `n_sites`, check `fpc_N`, or lower the target power.")
      n_pg <- people(x, how)
    }
    x0    <- solve_x(if (same(k_g)) "none" else "one_each", 0, Inf,
                     cluster = FALSE, fpc = FALSE)
    n_eff <- vapply(people(x0, if (same(k_g)) "none" else "one_each"), up, 0)
  }

  # Above a trillion people, rounding error in the calculation itself is
  # more than one person, so the answer could not be given exactly.
  if (max(n_pg) > 1e12) {
    clustered <- !no_sites && any(icc_g > 0)
    fixes <- c("prevalences further apart", "a lower `power` or `conf_level`")
    if (clustered && !is.null(n_per_site)) fixes <- c(fixes, "a smaller `n_per_site`")
    if (clustered && !is.null(n_sites))    fixes <- c(fixes, "more sites or a lower ICC")
    if (correction < 1)                    fixes <- c(fixes, "a more accurate test")
    stop("The required sample size is more than a trillion (1000000000000) ",
         "people per group, too large to work out to the nearest person. Use ",
         paste(fixes[-length(fixes)], collapse = ", "), ", or ",
         fixes[length(fixes)], ".")
  }

  # At least 1 person per group: with 0 versus 1 prevalence the formula can
  # give 0. Never more people than the population: rounding error in the
  # correction could otherwise give 2 people from a population of 1.
  n_pg <- vapply(n_pg, up, 0)
  if (!is.null(N_g)) n_pg <- pmin(n_pg, N_g)

  # ---- whole sites ----
  if (how == "per_site") {
    # A site larger than a whole group would be only partly sampled, so its
    # size would not be the real site size used in the design effect.
    if (all(s_g > n_pg)) {
      if (same(s_g) && same(n_pg))
        stop("`n_per_site` (", num(s_g[1]), ") is larger than the number of ",
             "people needed per group (", num(n_pg[1]), "), so one site would be ",
             "more than enough. Use a smaller `n_per_site`, or set ",
             "`n_sites = 1` to find how many people one site needs.")
      stop("`n_per_site` (", num(s_g[1]), " and ", num(s_g[2]), ") is larger than ",
           "the number of people needed in each group (", num(n_pg[1]), " and ",
           num(n_pg[2]), "), so one site per group would be more than enough. Use a ",
           "smaller `n_per_site`, or set `n_sites = 1` to find how many people ",
           "one site needs.")
    }
    k_out <- rep(up(x), 2)
  } else if (how == "sites") {
    if (sites_cut) {
      k_out <- k_used
      s_out <- c(1, 1)
      deff  <- c(1, 1)
      if (same(k_g) && same(n_pg)) {
        one <- n_pg[1] == 1
        warning("With these prevalences and this power, the study needs ",
                num(n_pg[1]), if (one) " person" else " people", " per group, ",
                "fewer than the ", num(k_g[1]), " sites given (`n_sites`). ",
                "The design uses ", num(n_pg[1]), if (one) " site" else " sites",
                " with 1 person at each.")
      } else {
        warning("With these prevalences and this power, the study needs ",
                num(n_pg[1]), " people in group 1 and ", num(n_pg[2]), " in group 2, ",
                "fewer than the ", num(k_g[1]), " and ", num(k_g[2]), " sites given ",
                "(`n_sites`). The design uses ", num(k_out[1]), " and ",
                num(k_out[2]), " sites, with 1 person at each.")
      }
    } else {
      k_out <- k_g
      s_out <- rep(up(x), 2)
      # Work deff out again from the rounded people per site, so that
      # deff == 1 + (n_per_site - 1) * icc describes the reported design.
      deff  <- 1 + (s_out - 1) * icc_g
    }
  }

  # Rounding the number of sites or people per site up can push the design
  # past the population. That only happens close to a census.
  if (!is.null(N_g) && !is.null(k_out)) {
    over <- k_out * s_out > N_g
    if (any(over)) {
      if (same(k_out) && same(s_out) && same(N_g))
        warning("The rounded design (", num(k_out[1]), " sites x ", num(s_out[1]),
                " = ", num(k_out[1] * s_out[1]), " per group) exceeds `fpc_N` = ",
                num(N_g[1]), ": the required n of ", num(n_pg[1]),
                " per group is close to a census of the population. Consider ",
                "surveying the whole population instead.")
      else for (g in which(over))
        warning("The rounded design for group ", g, " (", num(k_out[g]), " sites x ",
                num(s_out[g]), " = ", num(k_out[g] * s_out[g]), ") exceeds its `fpc_N` = ",
                num(N_g[g]), ": the required n of ", num(n_pg[g]), " is close to a ",
                "census of that group. Consider surveying the whole group instead.")
    }
  }

  # Without the population, these sites could not reach the target at all,
  # so the design can be very large.
  if (fpc_rescue) {
    tested <- if (same(n_pg) && same(N_g))
      paste0(num(n_pg[1]), " of the ", num(N_g[1]), " people in each group (",
             num(s_out[1]), " per site)")
    else
      paste0(num(n_pg[1]), " of the ", num(N_g[1]), " people in group 1 and ",
             num(n_pg[2]), " of the ", num(N_g[2]), " in group 2 (",
             num(s_out[1]), " per site)")
    warning("With ", design_text("and"), ", the target power ",
            "can only be reached because of the finite-population correction: ",
            "the design tests ", tested, ". Check that `fpc_N` is right",
            if (sites_cut) "." else ", or use more sites.")
  }

  list(
    n_per_group1   = n_pg[1],
    n_per_group2   = n_pg[2],
    n              = sum(n_pg),
    n_eff1         = n_eff[1],
    n_eff2         = n_eff[2],
    n_sites1       = k_out[1],
    n_sites2       = k_out[2],
    n_per_site1    = s_out[1],
    n_per_site2    = s_out[2],
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
    deff1          = deff[1],
    deff2          = deff[2],
    fpc_N          = fpc_N
  )
}
