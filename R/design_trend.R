#' Calculate sample size per timepoint to detect a linear trend in prevalence
#'
#' @description
#' Design function: given a starting prevalence and either an ending
#' prevalence or a slope, returns the minimum \strong{per-timepoint} sample
#' size needed to detect a linear time trend with the target power. Supply
#' \code{n} instead of \code{power} to run the calculation in reverse and
#' report the power a given per-timepoint sample size achieves.
#'
#' This is the design complement to \code{test_trend()}. Equal sample size at
#' every timepoint is assumed. The same three design modes as
#' \code{design_precision()} / \code{design_difference()} are supported, and
#' clustering / FPC are applied to \emph{each} timepoint's sample:
#' \describe{
#'   \item{SRS (no clustering args)}{Returns per-timepoint \code{n} only.}
#'   \item{Fixed cluster size (\code{n_per_site} supplied)}{Solves for the
#'     number of clusters needed per timepoint.}
#'   \item{Fixed number of clusters (\code{n_sites} supplied)}{Solves for
#'     the target samples per cluster via closed form.}
#' }
#'
#' @param prevalence_start Numeric in (0, 1). True prevalence at the first
#'   timepoint.
#' @param prevalence_end Numeric in (0, 1). True prevalence at the last
#'   timepoint. Supply exactly one of \code{prevalence_end} or \code{slope}.
#' @param slope Numeric. Change in true prevalence per unit time. Supply
#'   exactly one of \code{prevalence_end} or \code{slope}.
#' @param n_timepoints Integer >= 2. Number of equally spaced timepoints,
#'   used when \code{times} is not given (timepoints are then
#'   \code{0, 1, ..., n_timepoints - 1}).
#' @param times Optional numeric vector (length >= 2) of timepoint
#'   coordinates, if they are not equally spaced or you want explicit units
#'   (e.g. calendar years). Overrides \code{n_timepoints}.
#' @param n Optional positive integer. If supplied, the function reports the
#'   power achieved by this per-timepoint sample size instead of solving for
#'   \code{n}. When \code{n} is supplied, \code{power} is ignored.
#' @param power Target power (1 - beta); default 0.80.
#' @param alternative Direction of the alternative hypothesis; must match the
#'   \code{alternative} used in the companion \code{test_trend()} call.
#'   \code{"two.sided"} (default) tests slope != 0; \code{"greater"} tests
#'   slope > 0 (an increasing trend); \code{"less"} tests slope < 0.
#' @param sensitivity Diagnostic sensitivity in (0, 1]; default 1.
#' @param specificity Diagnostic specificity in (0, 1]; default 1.
#' @param conf_level Significance level expressed as a confidence level:
#'   alpha = 1 - conf_level. Default 0.95 (alpha = 0.05).
#' @param n_sites Optional positive integer. Fix the number of clusters per
#'   timepoint; the function solves for the target samples per cluster
#'   (\code{n_per_site}). Cannot be used with \code{n_per_site}.
#' @param n_per_site Optional positive integer. Fix the samples per cluster;
#'   the function solves for the number of clusters needed per timepoint
#'   (\code{n_sites}). Cannot be used with \code{n_sites}.
#' @param icc Numeric in \[0, 1\]. Intra-cluster correlation; default 0
#'   (SRS). Requires one of \code{n_sites} or \code{n_per_site} when
#'   \code{icc > 0}.
#' @param fpc_N Optional positive integer. Population size \strong{per
#'   timepoint} for a finite-population correction. \code{NULL} (default) =
#'   no FPC.
#'
#' @details
#' \strong{Slope z-test.} A weighted least-squares line through the apparent
#' prevalences at times \eqn{t_1, \dots, t_T} has slope estimator variance
#' approximately
#'
#' \deqn{\mathrm{Var}(\hat\beta_{app}) = \frac{\bar\sigma^2_{app}}{n \, S_{xx}},
#'   \qquad S_{xx} = \sum_i (t_i - \bar t)^2,}
#'
#' where \eqn{n} is the per-timepoint sample size and
#' \eqn{\bar\sigma^2_{app} = \bar p_{app}(1 - \bar p_{app})} is the binomial
#' variance at the mean apparent prevalence
#' \eqn{\bar p_{app} = (p_{start,app} + p_{end,app})/2}. Detecting a slope
#' \eqn{\beta_{app}} with power \eqn{1-\beta} needs
#'
#' \deqn{n = \frac{(z_\alpha + z_\beta)^2 \, \bar\sigma^2_{app}}
#'   {S_{xx} \, \beta_{app}^2},}
#'
#' with \eqn{z_\alpha = qnorm(1 - alpha)} (one-sided) or
#' \eqn{qnorm(1 - alpha/2)} (two-sided), \eqn{z_\beta = qnorm(power)}.
#'
#' The apparent slope relates to the true slope by
#' \eqn{\beta_{app} = \beta_{true} (Se + Sp - 1)}, so an imperfect test
#' shrinks the detectable trend and inflates \code{n}.
#'
#' Clustering and the FPC are then applied to each timepoint's sample with
#' the same closed-form logic as \code{design_threshold()} /
#' \code{design_difference()}.
#'
#' Solving for the \emph{number} of timepoints is not implemented in this
#' version; fix \code{n_timepoints} (or \code{times}) and solve for \code{n}.
#'
#' @section Equations and sources:
#' Building blocks from the MMS-SD Study Design Workshop
#' (\url{https://mrc-ide.github.io/MMS-SD_workshop/}); the trend / slope
#' model is \strong{not} covered there.
#' \itemize{
#'   \item \emph{Power form}
#'     \eqn{n = (z_{1-\beta} + z_{1-\alpha})^2\, \sigma^2/(S_{xx}\,\beta^2)}
#'     is the regression-slope analogue of Module 4's one-proportion
#'     result \eqn{n = (z_{1-\beta} + z_{1-\alpha/2})^2\,p(1-p)/(p-p_0)^2}
#'     (slide "Sample size formulae", lecture slides p. 7). Linear-trend
#'     sample size is not in the workshop; standard form, e.g. Dupont &
#'     Plummer (1998).
#'   \item \emph{Slope variance} \eqn{\mathrm{Var}(\hat\beta) =
#'     \sigma^2/(n\,S_{xx})}, \eqn{S_{xx} = \sum (t_i-\bar t)^2} --
#'     ordinary least squares; not workshop material.
#'   \item \emph{Design effect} \eqn{D_{eff} = 1 + (\bar n - 1)\,r} --
#'     Module 5 "Dealing with over-dispersion", slide "Why is the ICC
#'     useful?" (p. 5); applied per timepoint.
#'   \item \emph{Apparent-prevalence (Rogan-Gladen) scale}
#'     \eqn{\beta_{app} = \beta_{true}(Se + Sp - 1)} -- not in the
#'     workshop; Rogan & Gladen (1978).
#'   \item \emph{Finite-population correction} -- not in the workshop;
#'     Cochran (1977).
#' }
#'
#' @return A named list. Always present:
#'   \item{n_per_timepoint}{Per-timepoint sample size required (forward mode)
#'     or as supplied (reverse mode).}
#'   \item{n_total}{Total sample size across all timepoints
#'     (\code{n_per_timepoint * n_timepoints}).}
#'   \item{n_eff}{SRS-equivalent independent per-timepoint sample size. In
#'     forward mode: the base formula n, before the design effect and FPC.
#'     In reverse mode: the SRS-equivalent of the \emph{supplied} \code{n}
#'     after removing the FPC and design effect (the independent-sample size
#'     with the same slope-detection power) -- a different quantity, not
#'     comparable across modes.}
#'   \item{n_timepoints}{Number of timepoints.}
#'   \item{times}{Timepoint coordinates used.}
#'   \item{power}{Target power (forward mode) or power achieved by \code{n}
#'     (reverse mode).}
#'   \item{prevalence_start}{As supplied.}
#'   \item{prevalence_end}{As supplied, or implied by \code{slope}.}
#'   \item{slope}{True-scale slope (prevalence change per unit time).}
#'   \item{slope_app}{Apparent-scale slope.}
#'   \item{sxx}{\eqn{\sum (t_i - \bar t)^2}.}
#'   \item{sigma2_app}{Binomial variance at the mean apparent prevalence.}
#'   \item{alternative}{As supplied.}
#'   \item{conf_level}{As supplied.}
#'   \item{sensitivity}{As supplied.}
#'   \item{specificity}{As supplied.}
#'   \item{icc}{As supplied.}
#'   \item{deff}{Design effect applied.}
#'   \item{fpc_N}{As supplied, or \code{NULL}.}
#'   \item{mode}{\code{"solve_n"} or \code{"solve_power"}.}
#'
#'   Present when clustering is specified:
#'   \item{n_sites}{If \code{n_per_site} supplied: clusters required per
#'     timepoint. If \code{n_sites} supplied: echoed back. \code{NULL} for
#'     SRS.}
#'   \item{n_per_site}{If \code{n_sites} supplied: target samples per
#'     cluster. If \code{n_per_site} supplied: echoed back. \code{NULL} for
#'     SRS.}
#'
#' @references
#' MMS-SD Study Design Workshop, Modules 4 (statistical power) and 5 (ICC /
#' design effect). \url{https://mrc-ide.github.io/MMS-SD_workshop/}
#'
#' Dupont WD, Plummer WD (1998). Power and sample size calculations for
#' studies involving linear regression. Controlled Clinical Trials
#' 19(6):589-601.
#'
#' Rogan WJ, Gladen B (1978). Estimating prevalence from the results of a
#' screening test. American Journal of Epidemiology 107(1):71-76.
#'
#' Cochran WG (1977). Sampling Techniques, 3rd ed. Wiley.
#'
#' @seealso \code{\link{test_trend}}
#'
#' @export
#'
#' @examples
#' # Detect prevalence rising from 10% to 20% over 5 annual rounds
#' design_trend(prevalence_start = 0.10, prevalence_end = 0.20, n_timepoints = 5)
#'
#' # Same trend expressed as a slope
#' design_trend(prevalence_start = 0.10, slope = 0.025, n_timepoints = 5)
#'
#' # Reverse mode: power from 150 samples per round
#' design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5, n = 150)
#'
#' # Clustered: how many sites per round?
#' design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5,
#'              n_per_site = 25, icc = 0.03)
design_trend <- function(prevalence_start,
                         prevalence_end = NULL,
                         slope          = NULL,
                         n_timepoints   = NULL,
                         times          = NULL,
                         n              = NULL,
                         power          = 0.80,
                         alternative    = "two.sided",
                         sensitivity    = 1,
                         specificity    = 1,
                         conf_level     = 0.95,
                         n_sites        = NULL,
                         n_per_site     = NULL,
                         icc            = 0,
                         fpc_N          = NULL) {

  # ---- validate prevalence_start ----
  if (length(prevalence_start) != 1 || !is.numeric(prevalence_start))
    stop("`prevalence_start` must be a single number in (0, 1) (got class `",
         class(prevalence_start)[1], "`, length ", length(prevalence_start), ").")
  if (!is.finite(prevalence_start) || prevalence_start <= 0 || prevalence_start >= 1)
    stop("`prevalence_start` must be strictly between 0 and 1 (got ",
         prevalence_start, ").")

  # ---- resolve the effect: exactly one of prevalence_end / slope ----
  if (is.null(prevalence_end) && is.null(slope))
    stop("Supply exactly one of `prevalence_end` or `slope` to define the trend.")
  if (!is.null(prevalence_end) && !is.null(slope))
    stop("Supply exactly one of `prevalence_end` or `slope`, not both.")
  if (!is.null(prevalence_end)) {
    if (length(prevalence_end) != 1 || !is.numeric(prevalence_end))
      stop("`prevalence_end` must be a single number in (0, 1) (got class `",
           class(prevalence_end)[1], "`, length ", length(prevalence_end), ").")
    if (!is.finite(prevalence_end) || prevalence_end <= 0 || prevalence_end >= 1)
      stop("`prevalence_end` must be strictly between 0 and 1 (got ",
           prevalence_end, ").")
  }
  if (!is.null(slope) && (length(slope) != 1 || !is.numeric(slope) || !is.finite(slope)))
    stop("`slope` must be a single finite number (got ",
         if (!is.numeric(slope)) paste0("class `", class(slope)[1], "`")
         else paste0("length ", length(slope)), ").")

  # ---- resolve timepoints ----
  if (!is.null(times)) {
    if (!is.numeric(times) || length(times) < 2 || !all(is.finite(times)))
      stop("`times` must be a numeric vector of length >= 2 with finite values.")
    if (diff(range(times)) <= 0)
      stop("`times` must span a positive range (all values are equal).")
    tvec <- as.numeric(times)
  } else {
    if (is.null(n_timepoints) ||
        length(n_timepoints) != 1 || !is.numeric(n_timepoints) ||
        !is.finite(n_timepoints) || n_timepoints != floor(n_timepoints) ||
        n_timepoints < 2)
      stop("`n_timepoints` must be a single integer >= 2 (or supply `times`).")
    tvec <- 0:(n_timepoints - 1)
  }
  T_pts  <- length(tvec)
  t_span <- diff(range(tvec))
  sxx    <- sum((tvec - mean(tvec))^2)

  # ---- reverse vs forward ----
  solve_n <- is.null(n)
  if (!solve_n &&
      (!is.numeric(n) || length(n) != 1 || !is.finite(n) ||
       n != floor(n) || n < 1))
    stop("`n` must be a single finite positive integer (got ",
         if (length(n) != 1) paste0("length = ", length(n)) else n, ").")

  # ---- validate power ----
  if (solve_n) {
    if (length(power) != 1 || !is.numeric(power))
      stop("`power` must be a single number (got ",
           if (!is.numeric(power)) paste0("class `", class(power)[1], "`")
           else paste0("length ", length(power)), ").")
    if (!is.finite(power) || power <= 0 || power >= 1)
      stop("`power` must be strictly between 0 and 1 (got ", power, ").")
    if (power < 0.5)
      stop("`power` must be at least 0.5 (got ", power, "). ",
           "Below 0.5 the normal-approximation formula breaks down.")
  }

  # ---- validate remaining scalars ----
  if (!is.character(alternative) || length(alternative) != 1)
    stop("`alternative` must be a single string: 'greater', 'less', or 'two.sided'.")
  if (!alternative %in% c("greater", "less", "two.sided"))
    stop("`alternative` must be 'greater', 'less', or 'two.sided' (got '",
         alternative, "').")
  if (is.null(sensitivity) || length(sensitivity) != 1 || !is.numeric(sensitivity))
    stop("`sensitivity` must be a single number in (0, 1] (got ",
         if (is.null(sensitivity)) "NULL"
         else if (!is.numeric(sensitivity)) paste0("class `", class(sensitivity)[1], "`")
         else paste0("length = ", length(sensitivity)), "). ",
         "Note: `TRUE`/`FALSE` is logical, not numeric -- pass 1 for a perfect test.")
  if (is.null(specificity) || length(specificity) != 1 || !is.numeric(specificity))
    stop("`specificity` must be a single number in (0, 1] (got ",
         if (is.null(specificity)) "NULL"
         else if (!is.numeric(specificity)) paste0("class `", class(specificity)[1], "`")
         else paste0("length = ", length(specificity)), "). ",
         "Note: `TRUE`/`FALSE` is logical, not numeric -- pass 1 for a perfect test.")
  if (!is.finite(sensitivity) || sensitivity <= 0 || sensitivity > 1)
    stop("`sensitivity` must be in (0, 1] (got ", sensitivity, ").")
  if (!is.finite(specificity) || specificity <= 0 || specificity > 1)
    stop("`specificity` must be in (0, 1] (got ", specificity, ").")
  correction <- sensitivity + specificity - 1
  if (correction <= 0)
    stop("sensitivity + specificity must exceed 1 for the Rogan-Gladen correction ",
         "(got ", sensitivity, " + ", specificity, " = ", sensitivity + specificity, ").")
  if (correction < 0.1)
    warning("sensitivity + specificity = ", round(sensitivity + specificity, 4),
            " is very close to 1. The Rogan-Gladen adjustment is numerically unstable.")
  if (length(conf_level) != 1 || !is.numeric(conf_level))
    stop("`conf_level` must be a single number (got ",
         if (!is.numeric(conf_level)) paste0("class `", class(conf_level)[1], "`")
         else paste0("length = ", length(conf_level)), ").")
  if (!is.finite(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("`conf_level` must be strictly between 0 and 1 (got ", conf_level, ").")
  if (length(icc) != 1 || !is.numeric(icc))
    stop("`icc` must be a single number (got ",
         if (!is.numeric(icc)) paste0("class `", class(icc)[1], "`")
         else paste0("length ", length(icc)), "). Use 0 for an unclustered (SRS) design.")
  if (!is.finite(icc) || icc < 0 || icc > 1)
    stop("`icc` must be in [0, 1] (got ", icc, ").")
  if (!is.null(n_sites) && !is.null(n_per_site))
    stop("Supply at most one of `n_sites` or `n_per_site`, not both.")
  if (!is.null(n_sites) &&
      (!is.numeric(n_sites) || length(n_sites) != 1 || !is.finite(n_sites) ||
       n_sites != floor(n_sites) || n_sites < 1))
    stop("`n_sites` must be a single finite positive integer (got ",
         if (length(n_sites) != 1) paste0("length = ", length(n_sites)) else n_sites, ").")
  if (!is.null(n_per_site) &&
      (!is.numeric(n_per_site) || length(n_per_site) != 1 || !is.finite(n_per_site) ||
       n_per_site != floor(n_per_site) || n_per_site < 1))
    stop("`n_per_site` must be a single finite positive integer (got ",
         if (length(n_per_site) != 1) paste0("length = ", length(n_per_site)) else n_per_site, ").")
  if (!is.null(fpc_N) && (!is.numeric(fpc_N) || length(fpc_N) != 1 ||
      !is.finite(fpc_N) || fpc_N <= 0))
    stop("`fpc_N` must be a single finite positive number (got ",
         if (length(fpc_N) != 1) paste0("length = ", length(fpc_N)) else fpc_N, ").")
  if (icc > 0 && is.null(n_sites) && is.null(n_per_site))
    stop("icc > 0 requires a cluster structure: supply `n_sites` or `n_per_site`.")

  # ---- reverse mode: supplied n must fit the cluster layout ----
  # Without this, n < n_sites gives an average cluster size below 1, so the
  # Kish design effect drops below 1 and the reported power is overstated.
  if (!solve_n && !is.null(n_sites) && n < n_sites)
    stop("`n` (", n, ") is smaller than `n_sites` (", n_sites, "): the supplied ",
         "per-timepoint sample cannot provide even one person per site. ",
         "Increase `n`, or lower `n_sites`.")
  if (!solve_n && !is.null(n_per_site) && n < n_per_site)
    stop("`n` (", n, ") is smaller than `n_per_site` (", n_per_site, "): the ",
         "supplied per-timepoint sample is less than a single cluster. ",
         "Increase `n`, or lower `n_per_site`.")

  # ---- true-scale slope and endpoints ----
  if (!is.null(prevalence_end)) {
    slope_true <- (prevalence_end - prevalence_start) / t_span
    p_end      <- prevalence_end
  } else {
    slope_true <- slope
    p_end      <- prevalence_start + slope_true * t_span
    if (p_end <= 0 || p_end >= 1)
      stop("`slope` = ", slope, " over a time span of ", t_span,
           " drives prevalence to ", signif(p_end, 4),
           ", outside (0, 1). Reduce the slope or the time span.")
  }

  # ---- direction check ----
  if (alternative == "greater" && slope_true <= 0)
    stop("alternative = 'greater' expects an increasing trend, but the implied ",
         "slope is ", signif(slope_true, 4), " (<= 0).")
  if (alternative == "less" && slope_true >= 0)
    stop("alternative = 'less' expects a decreasing trend, but the implied ",
         "slope is ", signif(slope_true, 4), " (>= 0).")
  if (alternative == "two.sided" && abs(slope_true) < 1e-12)
    stop("The implied slope is effectively zero; there is no trend to power for.")

  # ---- apparent scale ----
  p_start_app <- .apparent_prev(prevalence_start, sensitivity, specificity)
  p_end_app   <- .apparent_prev(p_end,            sensitivity, specificity)
  slope_app   <- (p_end_app - p_start_app) / t_span
  pbar_app    <- (p_start_app + p_end_app) / 2
  sigma2_app  <- pbar_app * (1 - pbar_app)

  if (abs(slope_app) < 1e-10)
    stop("The apparent-scale slope is effectively zero (", signif(slope_app, 3),
         "), usually because an imperfect test compresses the trend. The ",
         "required sample size explodes -- use a larger trend or a better test.")

  alpha <- 1 - conf_level
  z_a   <- if (alternative == "two.sided") stats::qnorm(1 - alpha / 2) else
                                           stats::qnorm(1 - alpha)

  # ---- base per-timepoint sample size ----
  n_base_cont <- if (solve_n) {
    z_b <- stats::qnorm(power)
    (z_a + z_b)^2 * sigma2_app / (sxx * slope_app^2)
  } else {
    NA_real_   # reverse mode: no base-n solve
  }

  # ======================================================================
  # REVERSE MODE: n supplied -> report achieved power.
  # ======================================================================
  if (!solve_n) {
    if (icc == 0) {
      deff <- 1
    } else if (!is.null(n_per_site)) {
      deff <- 1 + (n_per_site - 1) * icc
    } else {
      deff <- 1 + (n / n_sites - 1) * icc
    }
    deff <- max(deff, 1)   # Kish deff is >= 1 by construction

    # Undo forward mode's operations in reverse order: forward applied the
    # FPC to the already-deff-inflated n, so invert the FPC on the actual
    # collected n first, then remove the design effect.
    n_cont_full <- n
    if (!is.null(fpc_N)) {
      if (n >= fpc_N)
        stop("`n` (", n, ") must be smaller than `fpc_N` (", fpc_N, ").")
      n_cont_full <- n * (fpc_N - 1) / (fpc_N - n)
    }
    n_eff_pt <- n_cont_full / deff

    ncp   <- abs(slope_app) * sqrt(n_eff_pt * sxx / sigma2_app)
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
      n_per_timepoint  = as.integer(n),
      n_total          = as.integer(n) * T_pts,
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
  # FORWARD MODE: solve for per-timepoint n. Clustering via the same
  # closed-form block as design_threshold() / design_difference().
  # ======================================================================
  if (icc == 0) {
    deff   <- 1
    n_cont <- n_base_cont

  } else if (!is.null(n_per_site)) {
    deff   <- 1 + (n_per_site - 1) * icc
    n_cont <- n_base_cont * deff

  } else {
    denom <- n_sites - n_base_cont * icc
    if (denom <= 0) {
      stop(sprintf(paste0(
        "Target power of %.0f%% is unachievable with %d sites per timepoint and ICC = %.3f.\n",
        "Increase `n_sites`, lower the ICC assumption, or relax the power target."
      ), 100 * power, n_sites, icc))
    }
    n_cont <- n_base_cont * n_sites * (1 - icc) / denom
    deff   <- n_cont / n_base_cont

    if (deff <= 1) {
      stop("n_sites = ", n_sites, " is >= the per-timepoint SRS sample size (~",
           ceiling(n_base_cont), "), so each site would receive < 1 person on ",
           "average -- not a valid cluster design. ",
           "Use n_sites < ", ceiling(n_base_cont), ", or supply `n_per_site` instead.")
    }
  }

  if (!is.null(fpc_N)) {
    n_cont <- (n_cont * fpc_N) / (n_cont + fpc_N - 1)
  }

  n_per_timepoint <- ceiling(n_cont)
  n_total         <- n_per_timepoint * T_pts
  n_eff           <- ceiling(n_base_cont)

  if (!is.null(n_per_site)) {
    n_sites_out    <- ceiling(n_per_timepoint / n_per_site)
    n_per_site_out <- n_per_site
  } else if (!is.null(n_sites)) {
    n_per_site_out <- ceiling(n_per_timepoint / n_sites)
    n_sites_out    <- n_sites
  } else {
    n_sites_out    <- NULL
    n_per_site_out <- NULL
  }

  list(
    n_per_timepoint  = n_per_timepoint,
    n_total          = n_total,
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
