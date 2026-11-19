#' Calculate sample size per timepoint to detect a linear trend in prevalence
#'
#' @description
#' Given the prevalence at the first timepoint and either the prevalence at
#' the last timepoint or a slope, returns the number of people needed
#' **per timepoint** to detect a linear trend with the target power, using
#' the test in [test_trend()]. Supply `n_per_timepoint` instead to run it in
#' reverse: the function then reports the power that number of people per
#' timepoint gives. Every timepoint gets the same design.
#'
#' **Three design modes** -- controlled by `n_sites` and `n_per_site`, at
#' each timepoint:
#' \describe{
#'   \item{Unclustered (`n_sites = NULL`, `n_per_site = NULL`)}{Returns the
#'     number of people per timepoint.}
#'   \item{Fixed site size (`n_per_site` supplied)}{Solves for the number of
#'     sites per timepoint.}
#'   \item{Fixed number of sites (`n_sites` supplied)}{Solves for the number
#'     of people per site.}
#' }
#'
#' @param prevalence_start Expected true prevalence at the first timepoint,
#'   in (0, 1).
#' @param prevalence_end Expected true prevalence at the last timepoint, in
#'   (0, 1). Supply exactly one of `prevalence_end` and `slope`.
#' @param slope Change in true prevalence per unit of time. Supply exactly
#'   one of `prevalence_end` and `slope`.
#' @param n_timepoints Whole number, 2 or more: equally spaced timepoints
#'   0, 1, ..., `n_timepoints` - 1. Supply exactly one of `n_timepoints` and
#'   `times`.
#' @param times Numeric vector of at least two timepoints, for unequal
#'   spacing or real units (for example calendar years).
#' @param n_per_timepoint Optional whole number. If supplied, the function
#'   runs in reverse and reports the power this many people per timepoint
#'   give; `power` is then not used.
#' @param power Target power, in \[0.5, 1); default 0.80.
#' @param alternative Direction of the test, as used later in
#'   [test_trend()]: `"two.sided"` (default), `"greater"` (an increase) or
#'   `"less"` (a decrease).
#' @param sensitivity Diagnostic sensitivity in (0, 1]; default 1.
#' @param specificity Diagnostic specificity in (0, 1]; default 1.
#'   `sensitivity + specificity` must exceed 1.
#' @param conf_level Confidence level, in (0, 1); default 0.95. The
#'   significance level of the test is 1 - `conf_level`.
#' @param n_sites Optional whole number. Fixes the number of sites per
#'   timepoint; the function solves for people per site. Cannot be used with
#'   `n_per_site`.
#' @param n_per_site Optional whole number. Fixes the number of people per
#'   site; the function solves for the number of sites per timepoint. Cannot
#'   be used with `n_sites`.
#' @param icc Intra-cluster correlation, in \[0, 1\]; default 0 (no
#'   clustering). An `icc` above 0 needs `n_sites` or `n_per_site`.
#'
#' @details
#' **Sample size per timepoint.** The test fits a straight line through the
#' measured prevalences at times \eqn{t_1, \dots, t_T}. With \eqn{n} people
#' per timepoint, its slope has variance about
#' \eqn{\bar\sigma^2_{app} / (n S_{xx})}, where
#' \eqn{S_{xx} = \sum_i (t_i - \bar t)^2} and
#' \eqn{\bar\sigma^2_{app} = \bar p_{app}(1 - \bar p_{app})} uses the
#' average of the first and last measured prevalences. Detecting a measured
#' slope \eqn{\beta_{app}} with the target power needs
#'
#' \deqn{n_0 = \frac{(z_\alpha + z_\beta)^2 \, \bar\sigma^2_{app}}
#'   {S_{xx} \, \beta_{app}^2},}
#'
#' with \eqn{z_\beta} the normal quantile at `power` and \eqn{z_\alpha} at
#' \eqn{1 - \alpha/2} (two-sided) or \eqn{1 - \alpha} (one-sided). The
#' measured slope is the true slope times \eqn{Se + Sp - 1}, so an imperfect
#' test needs more people (Dupont & Plummer 1998).
#'
#' **Clustering**, as in [design_precision()]: with `n_per_site` fixed,
#' \eqn{n = n_0 D_{eff}} with \eqn{D_{eff} = 1 + (n_{per site} - 1) ICC};
#' with `n_sites` fixed,
#' \eqn{n = n_0 n_{sites}(1 - ICC) / (n_{sites} - n_0 ICC)}, which exists
#' only when \eqn{n_{sites} > n_0 ICC}; otherwise the function stops and
#' reports the highest power those sites can give. If fewer people are
#' needed than there are sites, fewer sites are used, with 1 person at
#' each, and a warning says so. With sites, `n_per_timepoint` is the number
#' actually sampled, `n_sites * n_per_site`.
#'
#' **Reverse mode.** The same equation run backwards: the power of the
#' slope test with \eqn{n / D_{eff}} people per timepoint.
#'
#' @return A list with:
#'   \item{n_per_timepoint}{People per timepoint (or as supplied, in
#'     reverse). With sites, `n_sites * n_per_site`.}
#'   \item{n}{Total over all timepoints.}
#'   \item{n_sites}{Sites per timepoint (`NULL` if unclustered).}
#'   \item{n_per_site}{People per site (`NULL` if unclustered).}
#'   \item{deff}{Design effect at each timepoint; 1 if unclustered.}
#'   \item{n_eff}{Unclustered people per timepoint \eqn{n_0}, rounded up
#'     (in reverse, `n_per_timepoint / deff`).}
#'   \item{power}{The target (forward), or the power achieved (reverse).}
#'
#' @references
#' Dupont, W. D. & Plummer, W. D. (1998) Power and sample size calculations
#' for studies involving linear regression. *Controlled Clinical Trials*
#' **19**(6), 589-601.
#'
#' Kish, L. (1965) *Survey Sampling*. Wiley. (Design effect.)
#'
#' Rogan, W. J. & Gladen, B. (1978) Estimating prevalence from the results
#' of a screening test. *American Journal of Epidemiology* **107**(1),
#' 71-76. \doi{10.1093/oxfordjournals.aje.a112510}
#'
#' MMS-SD Study Design Workshop. \url{https://mrc-ide.github.io/MMS-SD_workshop/}
#'
#' @seealso [test_trend()]
#'
#' @export
#'
#' @examples
#' # Prevalence rising from 10 to 20 percent over 5 annual rounds
#' design_trend(prevalence_start = 0.10, prevalence_end = 0.20, n_timepoints = 5)
#'
#' # The same trend given as a slope
#' design_trend(prevalence_start = 0.10, slope = 0.025, n_timepoints = 5)
#'
#' # Reverse: the power of 150 people per round
#' design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5,
#'              n_per_timepoint = 150)
#'
#' # Fixed site size: how many sites per round?
#' design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5,
#'              n_per_site = 25, icc = 0.03)
design_trend <- function(prevalence_start,
                         prevalence_end  = NULL,
                         slope           = NULL,
                         n_timepoints    = NULL,
                         times           = NULL,
                         n_per_timepoint = NULL,
                         power           = 0.80,
                         alternative     = c("two.sided", "greater", "less"),
                         sensitivity     = 1,
                         specificity     = 1,
                         conf_level      = 0.95,
                         n_sites         = NULL,
                         n_per_site      = NULL,
                         icc             = 0) {

  .check_number(prevalence_start, "prevalence_start", 0, 1,
                lower_open = TRUE, upper_open = TRUE)
  if (is.null(prevalence_end) == is.null(slope))
    stop("Supply exactly one of `prevalence_end` and `slope`.", call. = FALSE)
  if (!is.null(prevalence_end))
    .check_number(prevalence_end, "prevalence_end", 0, 1,
                  lower_open = TRUE, upper_open = TRUE)
  else .check_number(slope, "slope")
  if (is.null(n_timepoints) == is.null(times))
    stop("Supply exactly one of `n_timepoints` and `times`.", call. = FALSE)
  .check_count(n_timepoints, "n_timepoints", min = 2)
  if (!is.null(times) && (!is.numeric(times) || length(times) < 2 ||
                          !all(is.finite(times)) || diff(range(times)) == 0))
    stop("`times` must be numbers with at least two different timepoints.",
         call. = FALSE)
  .check_count(n_per_timepoint, "n_per_timepoint")
  if (is.null(n_per_timepoint))
    .check_number(power, "power", 0.5, 1, upper_open = TRUE)
  .check_test(sensitivity, specificity)
  .check_number(conf_level, "conf_level", 0, 1, lower_open = TRUE, upper_open = TRUE)
  .check_sites(n_sites, n_per_site, icc)
  alternative <- match.arg(alternative)

  # timepoints and the trend
  tvec   <- if (!is.null(times)) as.numeric(times) else 0:(n_timepoints - 1)
  t_span <- diff(range(tvec))
  sxx    <- sum((tvec - mean(tvec))^2)
  p_end  <- if (!is.null(prevalence_end)) prevalence_end
            else prevalence_start + slope * t_span
  if (p_end <= 0 || p_end >= 1)
    stop("`slope` takes the prevalence at the last timepoint to ",
         signif(p_end, 4), ", outside (0, 1).", call. = FALSE)
  if (p_end == prevalence_start)
    stop("The prevalence must change between the first and last timepoint.",
         call. = FALSE)
  if (alternative == "greater" && p_end < prevalence_start)
    stop("The prevalence must rise for alternative = 'greater'.", call. = FALSE)
  if (alternative == "less" && p_end > prevalence_start)
    stop("The prevalence must fall for alternative = 'less'.", call. = FALSE)

  # measured slope and binomial variance (average of the two ends)
  p_start_app <- .apparent_prev(prevalence_start, sensitivity, specificity)
  p_end_app   <- .apparent_prev(p_end, sensitivity, specificity)
  slope_app   <- (p_end_app - p_start_app) / t_span
  pbar_app    <- (p_start_app + p_end_app) / 2
  sigma2_app  <- pbar_app * (1 - pbar_app)
  alpha <- 1 - conf_level
  z_a   <- stats::qnorm(if (alternative == "two.sided") 1 - alpha / 2 else 1 - alpha)
  power_at <- function(n_eff) {
    ncp <- abs(slope_app) * sqrt(n_eff * sxx / sigma2_app)
    stats::pnorm(ncp - z_a) +
      (if (alternative == "two.sided") stats::pnorm(-ncp - z_a) else 0)
  }
  n_tp <- length(tvec)

  # ---- reverse: power with n_per_timepoint people ----
  if (!is.null(n_per_timepoint)) {
    n <- n_per_timepoint
    if (!is.null(n_sites) && n < n_sites)
      stop("`n_per_timepoint` must be at least `n_sites`.", call. = FALSE)
    deff <- if (!is.null(n_per_site)) .deff(n_per_site, icc)
            else if (!is.null(n_sites)) .deff(n / n_sites, icc)
            else 1
    design <- .allocate(n, n_sites, n_per_site)
    return(list(
      n_per_timepoint = n,
      n               = n * n_tp,
      n_sites         = design$n_sites,
      n_per_site      = design$n_per_site,
      deff            = deff,
      n_eff           = n / deff,
      power           = power_at(n / deff)
    ))
  }

  # ---- forward: unclustered sample size per timepoint ----
  n0 <- (z_a + stats::qnorm(power))^2 * sigma2_app / (sxx * slope_app^2)

  # more sites than people needed: 1 person at just enough sites
  if (!is.null(n_sites) && n_sites > ceiling(n0)) {
    warning("Each timepoint needs ", ceiling(n0), " people, fewer than the ",
            n_sites, " sites given; the design uses ", ceiling(n0), " sites ",
            "with 1 person at each.", call. = FALSE)
    n_sites <- ceiling(n0)
  }

  # clustered sample size, as in design_precision()
  if (!is.null(n_per_site)) {
    n <- n0 * .deff(n_per_site, icc)
  } else if (!is.null(n_sites)) {
    n <- .solve_fixed_sites(n0, n_sites, icc)
    if (is.na(n))
      stop("The target power cannot be reached with ", n_sites, " sites per ",
           "timepoint and ICC = ", icc, ": the highest power these sites can ",
           "give is ", format(100 * power_at(n_sites / icc), digits = 3),
           " percent.", call. = FALSE)
  } else {
    n <- n0
  }

  design <- .allocate(n, n_sites, n_per_site)
  list(
    n_per_timepoint = design$n,
    n               = design$n * n_tp,
    n_sites         = design$n_sites,
    n_per_site      = design$n_per_site,
    deff            = if (is.null(design$n_per_site)) 1 else .deff(design$n_per_site, icc),
    n_eff           = ceiling(n0),
    power           = power
  )
}
