#' Calculate sample size per timepoint to detect a linear trend in prevalence
#'
#' @description
#' Design function: given the prevalence at the first timepoint and either
#' the prevalence at the last timepoint or a slope, returns the number of
#' people needed \strong{per timepoint} to detect a linear trend with the
#' target power. It is the design twin of \code{test_trend()}.
#'
#' Supply \code{n_per_timepoint} instead to run it in reverse: the function
#' then reports the power that number of people per timepoint gives.
#'
#' Every timepoint gets the same number of people. Clustering and the
#' finite-population correction are applied to each timepoint in the same
#' way. Three designs are supported:
#' \describe{
#'   \item{No clustering (neither \code{n_sites} nor \code{n_per_site})}{
#'     Returns the number of people per timepoint.}
#'   \item{Fixed number of people per site (\code{n_per_site})}{Solves for
#'     the number of sites needed per timepoint.}
#'   \item{Fixed number of sites (\code{n_sites})}{Solves for the number of
#'     people needed per site.}
#' }
#'
#' @param prevalence_start Expected true prevalence at the first timepoint,
#'   in (0, 1).
#' @param prevalence_end Expected true prevalence at the last timepoint, in
#'   (0, 1). Supply exactly one of \code{prevalence_end} or \code{slope}.
#' @param slope Change in true prevalence per unit of time. Supply exactly
#'   one of \code{prevalence_end} or \code{slope}.
#' @param n_timepoints Integer, 2 or more. Number of equally spaced
#'   timepoints, used when \code{times} is not given (the timepoints are then
#'   0, 1, ..., \code{n_timepoints} - 1).
#' @param times Optional numeric vector of at least two timepoints, for
#'   unequal spacing or real units (e.g. calendar years). Used instead of
#'   \code{n_timepoints}.
#' @param n_per_timepoint Optional positive integer. If supplied, the
#'   function runs in reverse and reports the power this number of people per
#'   timepoint gives; \code{power} is then not used.
#' @param power Target power, in \[0.5, 1); default 0.80.
#' @param alternative Direction of the test; must match the
#'   \code{alternative} used later in \code{test_trend()}.
#'   \code{"two.sided"} (default) tests for any trend; \code{"greater"} tests
#'   for an increase; \code{"less"} tests for a decrease.
#' @param sensitivity Diagnostic sensitivity in (0, 1]: the share of people
#'   who truly have the condition that the test detects as positive.
#'   Default 1 (the test never misses a case).
#' @param specificity Diagnostic specificity in (0, 1]: the share of people
#'   who truly do not have the condition that the test correctly reads as
#'   negative. Default 1 (the test never gives a false positive).
#' @param conf_level Confidence level, in \[0.5, 1); default 0.95. The
#'   significance level of the test is 1 - \code{conf_level}.
#' @param n_sites Optional positive integer. Fixes the number of sites per
#'   timepoint; the function solves for the number of people per site. Cannot
#'   be used together with \code{n_per_site}.
#' @param n_per_site Optional positive integer. Fixes the number of people
#'   per site; the function solves for the number of sites per timepoint.
#'   Cannot be used together with \code{n_sites}.
#' @param icc Intra-cluster correlation, in \[0, 1\]; default 0 (no
#'   clustering). An \code{icc} above 0 needs \code{n_sites} or
#'   \code{n_per_site}. Values below about 1.5e-8 are treated as 0.
#' @param fpc_N Optional positive integer. Population size \strong{per
#'   timepoint}, for a finite-population correction. \code{NULL} (default) =
#'   no finite-population correction applied.
#'
#' @details
#' \strong{Sample size per timepoint.} The test fits a straight line through
#' the measured (apparent) prevalences at times \eqn{t_1, \dots, t_T}. With
#' \eqn{n} people per timepoint, the variance of its slope is approximately
#'
#' \deqn{\mathrm{Var}(\hat\beta_{app}) = \frac{\bar\sigma^2_{app}}{n \, S_{xx}},
#'   \qquad S_{xx} = \sum_i (t_i - \bar t)^2,}
#'
#' where \eqn{\bar\sigma^2_{app} = \bar p_{app}(1 - \bar p_{app})} uses the
#' average of the first and last measured prevalences, \eqn{\bar p_{app}}.
#' Detecting a measured slope \eqn{\beta_{app}} with the target power needs
#'
#' \deqn{n_0 = \frac{(z_\alpha + z_\beta)^2 \, \bar\sigma^2_{app}}
#'   {S_{xx} \, \beta_{app}^2},}
#'
#' where \eqn{z_\beta} is the standard normal quantile at \code{power} and
#' \eqn{z_\alpha} the standard normal quantile at \eqn{1 - \alpha/2}
#' (two-sided) or \eqn{1 - \alpha} (one-sided), with
#' \eqn{\alpha = 1 -} \code{conf_level}. The measured slope is the true slope
#' times \eqn{Se + Sp - 1}, so an imperfect test shrinks the trend the test
#' can see and increases \eqn{n_0}.
#'
#' \strong{Clustering.} The design effect \eqn{D_{eff} = 1 + (\bar{n} -
#' 1)\,ICC}, with \eqn{\bar{n}} the number of people per site, is applied to
#' each timepoint. With \code{n_per_site} fixed, \eqn{n = n_0 D_{eff}}. With
#' \code{n_sites} fixed, \eqn{D_{eff}} depends on \eqn{n}, and solving gives
#' \eqn{n = n_0 \cdot n_{sites}(1 - ICC) / (n_{sites} - n_0 \cdot ICC)}. This
#' has a solution only when \eqn{n_{sites} > n_0 \cdot ICC}; otherwise the
#' function stops and reports the maximum power those sites can give. The
#' number of people per site is rounded up, and \code{deff} is then worked
#' out again from that rounded number.
#'
#' \strong{Finite-population correction.} With \code{fpc_N} (the population
#' size \eqn{N} at each timepoint), \eqn{n = n_1 N / (N + n_1 - 1)}, where
#' \eqn{n_1} is the size before the correction. With \code{n_sites} fixed,
#' the correction and the design effect are solved together (a quadratic with
#' one positive root, below \eqn{N}), so the target is always reachable.
#'
#' \strong{Reverse mode.} With \code{n_per_timepoint} given, the same
#' equations are run backwards: the effective size is
#' \eqn{n / (D_{eff} \cdot (N - n)/(N - 1))} (the last factor only with
#' \code{fpc_N}), and the power is that of the slope test with that size, so
#' a forward result put back in gives the target power.
#'
#' Solving for the number of timepoints is not supported: fix
#' \code{n_timepoints} (or \code{times}) and solve for the number of people.
#'
#' @return A named list with the following fields, in the order returned:
#'   \item{n_per_timepoint}{People per timepoint: solved for, or as supplied
#'     in reverse mode.}
#'   \item{n}{Total number of people over all timepoints,
#'     \code{n_per_timepoint * n_timepoints}.}
#'   \item{n_eff}{Effective number of people per timepoint with no
#'     clustering and no finite-population correction: \eqn{n_0} rounded up
#'     in forward mode; in reverse mode the effective size of the supplied
#'     design (\code{Inf} when the whole population is tested, which gives a
#'     power of 1).}
#'   \item{n_timepoints}{Number of timepoints.}
#'   \item{times}{The timepoints used.}
#'   \item{power}{Target power (forward), or the power achieved (reverse).}
#'   \item{prevalence_start}{As supplied.}
#'   \item{prevalence_end}{As supplied, or implied by \code{slope}.}
#'   \item{slope}{True change in prevalence per unit of time.}
#'   \item{slope_app}{Measured change in prevalence per unit of time.}
#'   \item{sxx}{\eqn{S_{xx} = \sum (t_i - \bar t)^2}.}
#'   \item{sigma2_app}{\eqn{\bar\sigma^2_{app}}, the binomial variance at the
#'     average measured prevalence.}
#'   \item{alternative}{As supplied.}
#'   \item{conf_level}{As supplied.}
#'   \item{sensitivity}{As supplied.}
#'   \item{specificity}{As supplied.}
#'   \item{icc}{As supplied.}
#'   \item{deff}{Design effect applied at each timepoint (1 with no
#'     clustering).}
#'   \item{fpc_N}{As supplied, or \code{NULL}.}
#'   \item{mode}{\code{"solve_n"} (forward) or \code{"solve_power"}
#'     (reverse).}
#'   \item{n_sites}{Sites per timepoint: solved for when \code{n_per_site} is
#'     given, echoed when \code{n_sites} is given. \code{NULL} with no
#'     clustering.}
#'   \item{n_per_site}{People per site: solved for when \code{n_sites} is
#'     given, echoed when \code{n_per_site} is given. \code{NULL} with no
#'     clustering.}
#'
#' @references
#' MMS-SD Study Design Workshop. \url{https://mrc-ide.github.io/MMS-SD_workshop/}
#'
#' Dupont, W. D. & Plummer, W. D. (1998) Power and sample size calculations
#' for studies involving linear regression. \emph{Controlled Clinical
#' Trials} \strong{19}(6), 589-601.
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
#' # Detect prevalence rising from 10 to 20 percent over 5 annual rounds
#' design_trend(prevalence_start = 0.10, prevalence_end = 0.20, n_timepoints = 5)
#'
#' # The same trend given as a slope
#' design_trend(prevalence_start = 0.10, slope = 0.025, n_timepoints = 5)
#'
#' # Reverse mode: the power of 150 people per round
#' design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5,
#'              n_per_timepoint = 150)
#'
#' # Clustered: how many sites per round?
#' design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5,
#'              n_per_site = 25, icc = 0.03)
design_trend <- function(prevalence_start,
                         prevalence_end  = NULL,
                         slope           = NULL,
                         n_timepoints    = NULL,
                         times           = NULL,
                         n_per_timepoint = NULL,
                         power           = 0.80,
                         alternative     = "two.sided",
                         sensitivity     = 1,
                         specificity     = 1,
                         conf_level      = 0.95,
                         n_sites         = NULL,
                         n_per_site      = NULL,
                         icc             = 0,
                         fpc_N           = NULL) {

  solve_n <- is.null(n_per_timepoint)

  # ---- validate single-value parameters ----
  # Each of the parameters below must be a single value.
  if (length(prevalence_start) != 1 || !is.numeric(prevalence_start))
    stop("`prevalence_start` must be a single number in (0, 1) (got class `",
         class(prevalence_start)[1], "`, length ", length(prevalence_start), ").")
  if (is.null(prevalence_end) && is.null(slope))
    stop("Supply exactly one of `prevalence_end` or `slope` to define the trend.")
  if (!is.null(prevalence_end) && !is.null(slope))
    stop("Supply exactly one of `prevalence_end` or `slope`, not both.")
  if (!is.null(prevalence_end) && (length(prevalence_end) != 1 || !is.numeric(prevalence_end)))
    stop("`prevalence_end` must be a single number in (0, 1) (got class `",
         class(prevalence_end)[1], "`, length ", length(prevalence_end), ").")
  if (!is.null(slope) && (length(slope) != 1 || !is.numeric(slope)))
    stop("`slope` must be a single number (got class `", class(slope)[1],
         "`, length ", length(slope), ").")
  if (solve_n && (length(power) != 1 || !is.numeric(power)))
    stop("`power` must be a single number in [0.5, 1) (got class `",
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
  if (!is.finite(prevalence_start))
    stop("`prevalence_start` must be a single finite number (got ", prevalence_start, ").")
  if (!is.null(prevalence_end) && !is.finite(prevalence_end))
    stop("`prevalence_end` must be a single finite number (got ", prevalence_end, ").")
  if (!is.null(slope) && !is.finite(slope))
    stop("`slope` must be a single finite number (got ", slope, ").")
  if (solve_n && !is.finite(power))
    stop("`power` must be a single finite number (got ", power, ").")
  if (!is.finite(sensitivity))
    stop("`sensitivity` must be a single finite number (got ", sensitivity, ").")
  if (!is.finite(specificity))
    stop("`specificity` must be a single finite number (got ", specificity, ").")
  if (!is.finite(conf_level))
    stop("`conf_level` must be a single finite number (got ", conf_level, ").")
  if (!is.finite(icc))
    stop("`icc` must be a single finite number (got ", icc, ").")

  for (nm in c("prevalence_start", "prevalence_end")) {
    p <- get(nm)
    if (is.null(p)) next
    if (p < 0 || p > 1)
      stop("`", nm, "` must be in (0, 1) (got ", p, "). `", nm, "` is a ",
           "fraction of the population and cannot be ",
           if (p < 0) "negative." else "greater than 1.")
    if (p == 0 || p == 1)
      stop("`", nm, "` must be in (0, 1) (got ", p, "). Use a value from a ",
           "pilot study, historical data, or a conservative guess.")
  }
  if (solve_n && (power <= 0 || power >= 1))
    stop("`power` must be in [0.5, 1) (got ", power, "). ",
         "For example, use 0.80 for 80 percent power.")
  if (solve_n && power < 0.5)
    stop("`power` must be at least 0.5 (got ", power, "). Below 0.5 the ",
         "sample-size formula stops increasing with power, so it no longer ",
         "gives a sensible answer.")
  if (sensitivity < 0 || sensitivity > 1)
    stop("`sensitivity` must be in (0, 1] (got ", sensitivity, "). ",
         "`sensitivity` is a diagnostic probability and cannot be ",
         if (sensitivity < 0) "negative." else "greater than 1.")
  if (sensitivity == 0)
    stop("`sensitivity` must be in (0, 1] (got 0). ",
         "A sensitivity of 0 means everyone who truly has the condition ",
         "tests negative.")
  if (specificity < 0 || specificity > 1)
    stop("`specificity` must be in (0, 1] (got ", specificity, "). ",
         "`specificity` is a diagnostic probability and cannot be ",
         if (specificity < 0) "negative." else "greater than 1.")
  if (specificity == 0)
    stop("`specificity` must be in (0, 1] (got 0). ",
         "A specificity of 0 means everyone who truly does not have the ",
         "condition tests positive.")
  # A sum above 1 only by rounding (within about 1.5e-8 of 1) counts as 1.
  correction <- sensitivity + specificity - 1
  if (correction < sqrt(.Machine$double.eps))
    stop("`sensitivity` + `specificity` must exceed 1 for the Rogan-Gladen correction ",
         "(got ", sensitivity, " + ", specificity, " = ", sensitivity + specificity, "). ",
         "The correction divides by `sensitivity` + `specificity` - 1, so a sum ",
         "of 1 or less cannot be corrected.")
  if (correction < 0.1)
    warning("`sensitivity` + `specificity` = ",
            format(sensitivity + specificity, digits = 15),
            ", which is very close to 1. Correcting for this much test error ",
            "makes the required n very large, and small changes in the assumed ",
            "sensitivity or specificity will change n a lot.")
  if (icc < 0 || icc > 1)
    stop("`icc` must be in [0, 1] (got ", icc, "). ",
         "`icc` is a correlation and cannot be ",
         if (icc < 0) "negative." else "greater than 1.")
  if (conf_level <= 0 || conf_level >= 1)
    stop("`conf_level` must be in [0.5, 1) (got ", conf_level, "). ",
         "For example, use 0.95 for a 95% confidence interval.")
  if (conf_level < 0.5)
    stop("`conf_level` must be at least 0.5 (got ", conf_level, "). Below 0.5 ",
         "the sample-size formula stops increasing with the confidence level, ",
         "so it no longer gives a sensible answer. Pass a confidence level such ",
         "as 0.95, not a significance level such as 0.05.")

  # ---- timepoints ----
  if (!is.null(times)) {
    if (!is.numeric(times) || length(times) < 2)
      stop("`times` must be a numeric vector of at least two timepoints (got ",
           "class `", class(times)[1], "`, length ", length(times), ").")
    if (!all(is.finite(times))) {
      i <- which(!is.finite(times))[1]
      stop("`times` contains a missing or infinite value (found times[", i,
           "] = ", times[i], ").")
    }
    if (diff(range(times)) <= 0)
      stop("`times` must contain at least two different values (got all ",
           times[1], ").")
    if (!is.null(n_timepoints))
      warning("Both `times` and `n_timepoints` were supplied; `n_timepoints` ",
              "is ignored (the timepoints come from `times`).")
    tvec <- as.numeric(times)
  } else {
    if (is.null(n_timepoints))
      stop("Supply `n_timepoints` (equally spaced timepoints) or `times`.")
    if (length(n_timepoints) != 1 || !is.numeric(n_timepoints) ||
        !is.finite(n_timepoints) || n_timepoints != floor(n_timepoints) ||
        n_timepoints < 2)
      stop("`n_timepoints` must be a single integer of 2 or more (got ",
           if (!is.numeric(n_timepoints)) paste0("class `", class(n_timepoints)[1], "`")
           else if (length(n_timepoints) != 1) paste0("length = ", length(n_timepoints))
           else n_timepoints, "), or supply `times`.")
    tvec <- 0:(n_timepoints - 1)
  }
  T_pts  <- length(tvec)
  t_span <- diff(range(tvec))
  sxx    <- sum((tvec - mean(tvec))^2)

  # ---- site and population inputs ----
  if (!solve_n &&
      (!is.numeric(n_per_timepoint) || length(n_per_timepoint) != 1 ||
       !is.finite(n_per_timepoint) || n_per_timepoint != floor(n_per_timepoint) ||
       n_per_timepoint < 1))
    stop("`n_per_timepoint` must be a single finite positive integer (got ",
         if (!is.numeric(n_per_timepoint)) paste0("class `", class(n_per_timepoint)[1], "`")
         else if (length(n_per_timepoint) != 1) paste0("length = ", length(n_per_timepoint))
         else n_per_timepoint, ").")
  if (!is.null(n_sites) && !is.null(n_per_site))
    stop("Supply at most one of `n_sites` or `n_per_site`, not both. ",
         "`n_sites` fixes the number of clusters and solves for samples per cluster; ",
         "`n_per_site` fixes the cluster size and solves for the number of clusters.")
  if (!is.null(n_sites) &&
      (!is.numeric(n_sites) || length(n_sites) != 1 || !is.finite(n_sites) ||
       n_sites != floor(n_sites) || n_sites < 1))
    stop("`n_sites` must be a single finite positive integer (got ",
         if (!is.numeric(n_sites)) paste0("class `", class(n_sites)[1], "`")
         else if (length(n_sites) != 1) paste0("length = ", length(n_sites))
         else n_sites, "). ",
         "`n_sites` is the number of sampling clusters per timepoint.")
  if (!is.null(n_per_site) &&
      (!is.numeric(n_per_site) || length(n_per_site) != 1 || !is.finite(n_per_site) ||
       n_per_site != floor(n_per_site) || n_per_site < 1))
    stop("`n_per_site` must be a single finite positive integer (got ",
         if (!is.numeric(n_per_site)) paste0("class `", class(n_per_site)[1], "`")
         else if (length(n_per_site) != 1) paste0("length = ", length(n_per_site))
         else n_per_site, "). ",
         "`n_per_site` is the fixed number of individuals sampled per cluster.")
  if (!is.null(fpc_N) &&
      (!is.numeric(fpc_N) || length(fpc_N) != 1 || !is.finite(fpc_N) ||
       fpc_N < 1 || fpc_N != floor(fpc_N)))
    stop("`fpc_N` must be a single finite positive integer (got ",
         if (!is.numeric(fpc_N)) paste0("class `", class(fpc_N)[1], "`")
         else if (length(fpc_N) != 1) paste0("length = ", length(fpc_N))
         else fpc_N,
         "). `fpc_N` is the population size per timepoint. ",
         "Set `fpc_N = NULL` to skip the finite-population correction.")

  # A cluster structure cannot be larger than the population it is drawn from.
  if (!is.null(fpc_N) && !is.null(n_sites) && n_sites > fpc_N)
    stop("`n_sites` (", n_sites, ") exceeds `fpc_N` (", fpc_N, "): you cannot ",
         "have more clusters than individuals in the population.")
  if (!is.null(fpc_N) && !is.null(n_per_site) && n_per_site > fpc_N)
    stop("`n_per_site` (", n_per_site, ") exceeds `fpc_N` (", fpc_N, "): a ",
         "cluster cannot be larger than the whole population.")

  # A supplied icc below this tiny threshold (e.g. 1e-12, a rounding
  # leftover) is treated as exactly 0.
  icc_is_zero <- icc < sqrt(.Machine$double.eps)
  if (!icc_is_zero && is.null(n_sites) && is.null(n_per_site))
    stop("`icc` is greater than 0 but there is no cluster structure (neither ",
         "`n_sites` nor `n_per_site` was supplied), so the design effect cannot ",
         "be computed. ",
         "Supply `n_sites` (fix the number of clusters) or `n_per_site` (fix the ",
         "cluster size), or set `icc = 0` for an unclustered design.")

  # Reverse mode: with n_sites fixed and icc above 0, fewer people than sites
  # would give a design effect below 1. (With n_per_site fixed the design
  # effect does not depend on n, so a number below one full site is fine:
  # forward mode can return one.)
  n <- n_per_timepoint
  if (!solve_n && !is.null(n_sites) && !icc_is_zero && n < n_sites)
    stop("`n_per_timepoint` (", n, ") is smaller than `n_sites` (", n_sites,
         "): there would be fewer than one person per site. Increase ",
         "`n_per_timepoint`, or lower `n_sites`.")
  # Testing the whole population (n = fpc_N) is allowed: there is then no
  # sampling uncertainty and the power is 1. Forward mode can return it when
  # the required n rounds up to the population size.
  if (!solve_n && !is.null(fpc_N) && n > fpc_N)
    stop("`n_per_timepoint` (", n, ") is larger than `fpc_N` (", fpc_N,
         "): you cannot test more people than there are in the population.")

  # ---- true slope and last prevalence ----
  if (!is.null(prevalence_end)) {
    slope_true <- (prevalence_end - prevalence_start) / t_span
    p_end      <- prevalence_end
  } else {
    slope_true <- slope
    p_end      <- prevalence_start + slope_true * t_span
    if (p_end <= 0 || p_end >= 1)
      stop("`slope` = ", slope, " over a time span of ", t_span,
           " takes the prevalence to ", signif(p_end, 4), ", outside (0, 1). ",
           "Reduce the slope or the time span.")
  }

  # ---- direction of the test ----
  if (alternative == "greater" && slope_true <= 0)
    stop("alternative = 'greater' needs an increasing trend, but the slope is ",
         signif(slope_true, 4), ".")
  if (alternative == "less" && slope_true >= 0)
    stop("alternative = 'less' needs a decreasing trend, but the slope is ",
         signif(slope_true, 4), ".")
  if (alternative == "two.sided" && abs(slope_true) < 1e-12)
    stop("The slope is zero, or nearly so, so there is no trend to detect.")

  # ---- measured (apparent) scale ----
  p_start_app <- .apparent_prev(prevalence_start, sensitivity, specificity)
  p_end_app   <- .apparent_prev(p_end,            sensitivity, specificity)
  slope_app   <- (p_end_app - p_start_app) / t_span
  pbar_app    <- (p_start_app + p_end_app) / 2
  sigma2_app  <- pbar_app * (1 - pbar_app)
  if (abs(slope_app) < 1e-10)
    stop("The measured slope is only ", signif(slope_app, 3), ", so the ",
         "required sample size would be astronomically large. Use a larger ",
         "trend, or a more accurate test.")

  alpha <- 1 - conf_level
  z_a   <- if (alternative == "two.sided") stats::qnorm(1 - alpha / 2) else
                                           stats::qnorm(1 - alpha)

  # ======================================================================
  # Reverse mode: n per timepoint given -> power achieved. Runs the forward
  # equations backwards: effective size = n / (deff * FPC factor), with
  # deff from n itself when n_sites is fixed.
  # ======================================================================
  if (!solve_n) {
    deff <- if (icc_is_zero) 1
            else if (!is.null(n_per_site)) 1 + (n_per_site - 1) * icc
            else 1 + (n / n_sites - 1) * icc
    f_n      <- if (is.null(fpc_N)) 1 else (fpc_N - n) / (fpc_N - 1)
    n_eff_pt <- n / (deff * f_n)

    ncp <- abs(slope_app) * sqrt(n_eff_pt * sxx / sigma2_app)
    power_out <- if (alternative == "two.sided")
      stats::pnorm(ncp - z_a) + stats::pnorm(-ncp - z_a)
    else
      stats::pnorm(ncp - z_a)

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
      n_per_timepoint  = as.numeric(n),
      n                = as.numeric(n) * T_pts,
      n_eff            = ceiling(n_eff_pt),
      n_timepoints     = T_pts,
      times            = tvec,
      power            = power_out,
      prevalence_start = prevalence_start,
      prevalence_end   = p_end,
      slope            = slope_true,
      slope_app        = slope_app,
      sxx              = sxx,
      sigma2_app       = sigma2_app,
      alternative      = alternative,
      conf_level       = conf_level,
      sensitivity      = sensitivity,
      specificity      = specificity,
      icc              = icc,
      deff             = deff,
      fpc_N            = fpc_N,
      mode             = "solve_power",
      n_sites          = n_sites_out,
      n_per_site       = n_per_site_out
    ))
  }

  # ======================================================================
  # Forward mode: solve for n per timepoint.
  # ======================================================================
  z_b <- stats::qnorm(power)
  n_base_cont <- (z_a + z_b)^2 * sigma2_app / (sxx * slope_app^2)

  # Reaching the first branch with both n_sites and n_per_site NULL means
  # icc is 0: an icc above 0 without them was rejected above. The
  # fixed-n_sites branch folds the FPC into its own solve; every other
  # branch has it applied afterwards.
  fpc_applied <- FALSE

  if (icc_is_zero) {
    deff   <- 1
    n_cont <- n_base_cont

  } else if (!is.null(n_per_site)) {
    deff   <- 1 + (n_per_site - 1) * icc
    n_cont <- n_base_cont * deff

  } else {
    # Number of sites fixed: deff depends on n, and n depends on deff.
    if (is.null(fpc_N)) {
      # Solving n = n_base * (1 + (n / n_sites - 1) * icc) for n needs
      # n_sites > n_base * icc. The highest power these sites can give is
      # where n_base * icc = n_sites exactly.
      denom <- n_sites - n_base_cont * icc
      if (denom <= 0) {
        z_b_max   <- abs(slope_app) * sqrt(n_sites * sxx / (icc * sigma2_app)) - z_a
        max_power <- stats::pnorm(z_b_max)
        stop(sprintf(paste0(
          "Target power of %.1f percent is unachievable with %d sites per timepoint ",
          "and ICC = %.3f.\n",
          "Maximum achievable power with these settings: %.1f percent.\n",
          "Increase `n_sites`, lower the ICC assumption, or relax the power target."
        ), 100 * power, n_sites, icc, 100 * max_power))
      }
      n_cont <- n_base_cont * n_sites * (1 - icc) / denom

    } else {
      # With an FPC, solve n_base * (1 + (n / n_sites - 1) * icc) *
      # (N - n) / (N - 1) = n for n. This is the quadratic A*n^2 + B*n + C = 0
      # with A > 0 and C < 0, so it has exactly one positive root, below N.
      # The root is taken in a form that stays accurate for very large N.
      A <- n_base_cont * icc / n_sites
      B <- n_base_cont * (1 - icc) - A * fpc_N + fpc_N - 1
      C <- -n_base_cont * (1 - icc) * fpc_N
      q <- -(B + (if (B < 0) -1 else 1) * sqrt(B^2 - 4 * A * C)) / 2
      n_cont <- max(q / A, C / q)
      fpc_applied <- TRUE
    }
    deff <- 1 + (n_cont / n_sites - 1) * icc

    # deff <= 1 means fewer than one person per site on average: there are
    # more sites than people needed with no clustering.
    if (deff <= 1) {
      n_srs <- if (is.null(fpc_N)) n_base_cont
               else (n_base_cont * fpc_N) / (n_base_cont + fpc_N - 1)
      stop("n_sites = ", n_sites, " is >= the number of people needed per ",
           "timepoint with no clustering (about ", ceiling(n_srs), "), so each ",
           "site would receive fewer than 1 person on average, which is not a ",
           "valid cluster design. Use n_sites < ", ceiling(n_srs), ", or supply ",
           "`n_per_site` to fix the cluster size and solve for the number of ",
           "sites instead.")
    }
  }

  # ---- finite-population correction (per timepoint) ----
  if (!is.null(fpc_N) && !fpc_applied)
    n_cont <- (n_cont * fpc_N) / (n_cont + fpc_N - 1)

  n_per_tp <- ceiling(n_cont)
  n_total  <- n_per_tp * T_pts
  n_eff    <- ceiling(n_base_cont)

  # ---- distribute across sites ----
  if (!is.null(n_per_site)) {
    n_sites_out    <- ceiling(n_per_tp / n_per_site)
    n_per_site_out <- n_per_site
  } else if (!is.null(n_sites)) {
    n_per_site_out <- ceiling(n_per_tp / n_sites)
    n_sites_out    <- n_sites
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
    warning("The rounded design (", n_sites_out, " sites x ", n_per_site_out,
            " = ", n_sites_out * n_per_site_out, " per timepoint) exceeds ",
            "`fpc_N` = ", fpc_N, ": the required n of ", n_per_tp, " per ",
            "timepoint is close to a census of the population. Consider ",
            "surveying the whole population instead.")

  list(
    n_per_timepoint  = n_per_tp,
    n                = n_total,
    n_eff            = n_eff,
    n_timepoints     = T_pts,
    times            = tvec,
    power            = power,
    prevalence_start = prevalence_start,
    prevalence_end   = p_end,
    slope            = slope_true,
    slope_app        = slope_app,
    sxx              = sxx,
    sigma2_app       = sigma2_app,
    alternative      = alternative,
    conf_level       = conf_level,
    sensitivity      = sensitivity,
    specificity      = specificity,
    icc              = icc,
    deff             = deff,
    fpc_N            = fpc_N,
    mode             = "solve_n",
    n_sites          = n_sites_out,
    n_per_site       = n_per_site_out
  )
}
