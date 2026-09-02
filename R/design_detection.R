#' Calculate sample size to detect at least one occurrence of a variant
#'
#' @description
#' Design function: given an expected (possibly rare) true prevalence,
#' returns the minimum sample size needed so that the probability of
#' observing at least one positive test result reaches a target
#' (\code{detection_prob}). Supply \code{n} instead of \code{detection_prob}
#' to run the calculation in reverse and report the detection probability a
#' given sample size achieves.
#'
#' The same three design modes as \code{design_precision()} /
#' \code{design_threshold()} are supported:
#' \describe{
#'   \item{SRS (no clustering args)}{Returns total \code{n} only.}
#'   \item{Fixed cluster size (\code{n_per_site} supplied)}{Solves for the
#'     number of clusters needed.}
#'   \item{Fixed number of clusters (\code{n_sites} supplied)}{Solves for
#'     target samples per cluster via closed form.}
#' }
#'
#' @param prevalence Numeric in (0, 1). Minimum true prevalence of the
#'   variant to detect (the prevalence at which detection probability is
#'   evaluated).
#' @param detection_prob Numeric in (0, 1). Target probability of observing
#'   at least one positive; default 0.95. Ignored when \code{n} is supplied.
#' @param n Optional positive integer. If supplied, the function solves the
#'   inverse problem: it reports the detection probability achieved by this
#'   total sample size and does not solve for \code{n}. When \code{n} is
#'   supplied, \code{detection_prob} is ignored.
#' @param sensitivity Numeric in (0, 1]. Diagnostic sensitivity; default 1.
#' @param specificity Numeric in (0, 1]. Diagnostic specificity; default 1.
#'   \strong{Currently validated but NOT used in the calculation} -- see the
#'   "Known limitation" section below.
#' @param n_sites Optional positive integer. Fix the number of clusters; the
#'   function solves for the target samples per cluster
#'   (\code{n_per_site}). Cannot be used with \code{n_per_site}.
#' @param n_per_site Optional positive integer. Fix the samples per cluster;
#'   the function solves for the number of clusters needed (\code{n_sites}).
#'   Cannot be used with \code{n_sites}.
#' @param icc Numeric in \[0, 1\]. Intra-cluster correlation; default 0
#'   (SRS). Requires one of \code{n_sites} or \code{n_per_site} when
#'   \code{icc > 0}.
#' @param fpc_N Optional positive integer. Total population size for a
#'   finite-population correction. \code{NULL} (default) = no FPC.
#'
#' @details
#' \strong{Detection formula.} With a perfect test the number of positives
#' in \code{n} independent samples is binomial, so the probability of at
#' least one positive is \eqn{1 - (1 - q)^n} and the required SRS sample
#' size is
#'
#' \deqn{n = \frac{\log(1 - detection\_prob)}{\log(1 - q)}.}
#'
#' \code{q} is the probability that a single randomly drawn sample is a
#' \emph{genuine} positive that the assay detects:
#'
#' \deqn{q = prevalence \times sensitivity.}
#'
#' This is the standard basis used in surveillance / "freedom from disease"
#' sample-size calculations (Cannon & Roe 1982; Cannon 2001): the
#' calculation is about the chance of catching a \emph{true} occurrence, and
#' specificity is handled separately (e.g. by confirmatory testing) rather
#' than being allowed to shrink the required \code{n}.
#'
#' Clustering and the finite-population correction are then applied with the
#' same closed-form logic as \code{design_threshold()} /
#' \code{design_precision()}.
#'
#' @section Known limitation -- specificity is not yet used:
#' \strong{This is a deliberate v1 simplification that must be revisited.}
#' \code{design_detection()} currently sets \eqn{q = prevalence \times
#' sensitivity} and ignores \code{specificity} entirely. Two questions are
#' unresolved and need a team decision:
#' \enumerate{
#'   \item Whether "detection" should mean "detect a true occurrence"
#'     (current behaviour, \eqn{q = p \cdot Se}) or "observe any positive
#'     result" (the apparent-prevalence basis \eqn{q = p \cdot Se +
#'     (1 - p)(1 - Sp)} used by every other function in this package via
#'     \code{.apparent_prev()}).
#'   \item If false positives are in scope, how a "detection" driven by a
#'     false positive should be reported / discounted.
#' }
#' Until this is settled, results with \code{specificity < 1} should be
#' treated as provisional. See the README for the tracking note.
#'
#' @return A named list. Always present:
#'   \item{n}{Total sample size required (forward mode) or as supplied
#'     (reverse mode).}
#'   \item{n_eff}{SRS-equivalent independent sample size. In forward mode:
#'     the base detection-formula n before the design effect and FPC,
#'     defined the same way as \code{n_eff} in \code{design_precision()} /
#'     \code{design_threshold()}. In reverse mode: the SRS-equivalent of the
#'     \emph{supplied} \code{n} after removing the design effect and FPC
#'     (i.e. the independent-sample size whose detection probability equals
#'     the reported one) -- a different quantity, not comparable across
#'     modes.}
#'   \item{detection_prob}{Target detection probability (forward mode) or the
#'     probability achieved by \code{n} (reverse mode).}
#'   \item{prevalence}{Minimum prevalence to detect, as supplied.}
#'   \item{q}{Per-sample true-positive probability used in the formula
#'     (\code{prevalence * sensitivity}).}
#'   \item{sensitivity}{Sensitivity as supplied.}
#'   \item{specificity}{Specificity as supplied (not used; see above).}
#'   \item{icc}{ICC as supplied.}
#'   \item{deff}{Design effect applied.}
#'   \item{fpc_N}{\code{fpc_N} as supplied, or \code{NULL}.}
#'   \item{mode}{\code{"solve_n"} or \code{"solve_detection_prob"}.}
#'
#'   Present when clustering is specified:
#'   \item{n_sites}{If \code{n_per_site} supplied: clusters required. If
#'     \code{n_sites} supplied: echoed back. \code{NULL} for SRS.}
#'   \item{n_per_site}{If \code{n_sites} supplied: target samples per
#'     cluster. If \code{n_per_site} supplied: echoed back. \code{NULL} for
#'     SRS.}
#'
#' @references
#' Cannon RM, Roe RT (1982). Livestock Disease Surveys: A Field Manual for
#' Veterinarians. Australian Bureau of Animal Health.
#'
#' Cannon RM (2001). Sense and sensitivity -- designing surveys based on an
#' imperfect test. Preventive Veterinary Medicine 49(3):141-163.
#'
#' MMS-SD workshop Module 6 (DRpower presence/detection logic).
#'
#' @export
#'
#' @examples
#' # How many samples to be 95% sure of seeing >=1 positive when p = 0.02?
#' design_detection(prevalence = 0.02)
#'
#' # Imperfect sensitivity inflates the requirement
#' design_detection(prevalence = 0.02, sensitivity = 0.8)
#'
#' # Reverse mode: what detection probability does n = 150 buy me?
#' design_detection(prevalence = 0.02, n = 150)
#'
#' # Fixed cluster size: how many sites needed?
#' design_detection(prevalence = 0.02, n_per_site = 20, icc = 0.05)
design_detection <- function(prevalence,
                             detection_prob = 0.95,
                             n              = NULL,
                             sensitivity    = 1,
                             specificity    = 1,
                             n_sites        = NULL,
                             n_per_site     = NULL,
                             icc            = 0,
                             fpc_N          = NULL) {

  # ---- validate scalars ----
  if (length(prevalence) != 1 || !is.numeric(prevalence))
    stop("`prevalence` must be a single number in (0, 1) (got class `",
         class(prevalence)[1], "`, length ", length(prevalence), ").")
  if (!is.finite(prevalence) || prevalence <= 0 || prevalence >= 1)
    stop("`prevalence` must be strictly between 0 and 1 (got ", prevalence, ").")

  solve_n <- is.null(n)

  if (length(detection_prob) != 1 || !is.numeric(detection_prob))
    stop("`detection_prob` must be a single number (got ",
         if (!is.numeric(detection_prob))
           paste0("class `", class(detection_prob)[1], "`")
         else paste0("length ", length(detection_prob)), ").")
  if (!is.finite(detection_prob) || detection_prob <= 0 || detection_prob >= 1)
    stop("`detection_prob` must be strictly between 0 and 1 (got ",
         detection_prob, "). Use, e.g., 0.95 for 95% detection probability.")

  if (!solve_n &&
      (!is.numeric(n) || length(n) != 1 || !is.finite(n) ||
       n != floor(n) || n < 1))
    stop("`n` must be a single finite positive integer (got ",
         if (length(n) != 1) paste0("length = ", length(n)) else n, ").")

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

  # ---- reverse-mode: supplied n must be consistent with the cluster layout ----
  # Without this, n < n_sites gives an average cluster size below 1, so the
  # Kish design effect drops below 1, inflating n_eff above the n actually
  # collected and (with fpc_N) driving the FPC inversion denominator negative.
  if (!solve_n && !is.null(n_sites) && n < n_sites)
    stop("`n` (", n, ") is smaller than `n_sites` (", n_sites, "): the supplied ",
         "sample size cannot provide even one person per site. Increase `n`, ",
         "or lower `n_sites`.")
  if (!solve_n && !is.null(n_per_site) && n < n_per_site)
    stop("`n` (", n, ") is smaller than `n_per_site` (", n_per_site, "): the ",
         "supplied sample size is less than a single cluster. Increase `n`, or ",
         "lower `n_per_site`.")

  # ---- per-sample true-positive probability ----
  # KNOWN LIMITATION: specificity is deliberately omitted here. `q` is the
  # probability a sample is a genuine positive AND the assay flags it. See
  # the "Known limitation" roxygen section and the README tracking note --
  # this must be revisited with the team before the function is treated as
  # validated.
  q <- prevalence * sensitivity
  if (q <= 0 || q >= 1)
    stop("`prevalence * sensitivity` = ", signif(q, 4),
         " is outside (0, 1); the detection formula is undefined.")

  # ---- base SRS size from the geometric detection formula ----
  # n = log(1 - detection_prob) / log(1 - q)
  n_base_cont <- log(1 - detection_prob) / log(1 - q)

  # ======================================================================
  # REVERSE MODE: n supplied -> report the detection probability it buys.
  # deff is directly computable here (no circularity).
  # ======================================================================
  if (!solve_n) {
    if (icc == 0) {
      deff <- 1
    } else if (!is.null(n_per_site)) {
      deff <- 1 + (n_per_site - 1) * icc
    } else {
      # n_sites fixed: implied samples per cluster is n / n_sites
      deff <- 1 + (n / n_sites - 1) * icc
    }
    deff   <- max(deff, 1)   # Kish deff is >= 1 by construction

    n_cont <- n / deff

    # undo the finite-population correction to recover the SRS-equivalent n
    if (!is.null(fpc_N)) {
      if (n >= fpc_N)
        stop("`n` (", n, ") must be smaller than `fpc_N` (", fpc_N, ").")
      n_cont <- n_cont * (fpc_N - 1) / (fpc_N - n_cont)
    }

    achieved <- 1 - (1 - q)^n_cont

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
      n              = as.numeric(n),
      n_eff          = ceiling(n_cont),
      detection_prob = achieved,
      prevalence     = prevalence,
      q              = q,
      sensitivity    = sensitivity,
      specificity    = specificity,
      icc            = icc,
      deff           = deff,
      fpc_N          = fpc_N,
      mode           = "solve_detection_prob",
      n_sites        = n_sites_out,
      n_per_site     = n_per_site_out
    ))
  }

  # ======================================================================
  # FORWARD MODE: solve for n. Design effect via the same closed-form
  # logic as design_threshold() / design_precision() (copied inline;
  # slated for extraction into .solve_cluster_design()).
  # ======================================================================
  # icc == 0 covers every unclustered case: an earlier guard errored if
  # icc > 0 without a cluster structure.
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
        "Detection probability of %.0f%% is unachievable with %d sites and ICC = %.3f.\n",
        "Increase `n_sites`, lower the ICC assumption, or relax `detection_prob`."
      ), 100 * detection_prob, n_sites, icc))
    }
    n_cont <- n_base_cont * n_sites * (1 - icc) / denom
    deff   <- n_cont / n_base_cont

    if (deff <= 1) {
      stop("n_sites = ", n_sites, " is >= the SRS sample size (~",
           ceiling(n_base_cont), "), so each site would receive < 1 person on ",
           "average -- not a valid cluster design. ",
           "Use n_sites < ", ceiling(n_base_cont), ", or supply `n_per_site` instead.")
    }
  }

  # ---- finite-population correction ----
  if (!is.null(fpc_N)) {
    # A finite population caps the detection probability: even a full census
    # of fpc_N samples gives only 1 - (1 - q)^fpc_N. If the target exceeds
    # that, the FPC would otherwise silently return a truncated n that never
    # reaches detection_prob.
    max_achievable <- 1 - (1 - q)^fpc_N
    if (detection_prob > max_achievable)
      stop("A detection probability of ", round(100 * detection_prob, 1),
           "% is not achievable in a population of ", fpc_N,
           ": even a full census reaches at most ",
           round(100 * max_achievable, 1), "% (prevalence ", prevalence,
           ", sensitivity ", sensitivity, "). Lower `detection_prob`, or drop `fpc_N`.")
    n_cont <- (n_cont * fpc_N) / (n_cont + fpc_N - 1)
  }

  n_total <- ceiling(n_cont)

  # n_eff: SRS-equivalent independent sample size (base detection-formula n
  # before the design effect and FPC). Same definition as n_eff in
  # design_precision() / design_threshold().
  n_eff <- ceiling(n_base_cont)

  # ---- distribute across sites ----
  if (!is.null(n_per_site)) {
    n_sites_out    <- ceiling(n_total / n_per_site)
    n_per_site_out <- n_per_site
  } else if (!is.null(n_sites)) {
    n_per_site_out <- ceiling(n_total / n_sites)
    n_sites_out    <- n_sites
  } else {
    n_sites_out    <- NULL
    n_per_site_out <- NULL
  }

  list(
    n              = n_total,
    n_eff          = n_eff,
    detection_prob = detection_prob,
    prevalence     = prevalence,
    q              = q,
    sensitivity    = sensitivity,
    specificity    = specificity,
    icc            = icc,
    deff           = deff,
    fpc_N          = fpc_N,
    mode           = "solve_n",
    n_sites        = n_sites_out,
    n_per_site     = n_per_site_out
  )
}
