#' Estimate total study cost
#'
#' @description
#' Calculates the total cost of a surveillance or prevalence study from a
#' per-sample cost and an optional fixed cost per site. Designed to pair
#' directly with the output of \code{design_precision()} or
#' \code{design_detection()}: pass their \code{n} and \code{n_sites} return
#' values straight into \code{estimate_cost()}.
#'
#' All monetary inputs and outputs are in the same currency unit -- the
#' function does no currency conversion.
#'
#' @param n Positive integer. Total number of samples (e.g. the \code{n}
#'   element returned by \code{design_precision()} or
#'   \code{design_detection()}).
#' @param cost_per_sample Non-negative number. Cost per individual sample
#'   (e.g. cost of a rapid diagnostic test, lab processing fee, or combined
#'   per-person field cost).
#' @param n_sites Optional positive integer. Number of study sites. When
#'   supplied together with \code{cost_per_site}, a fixed site-level cost is
#'   added to the total. Set to \code{NULL} (default) to omit site costs.
#' @param cost_per_site Non-negative number. Fixed cost per site (e.g.
#'   travel, equipment setup, staff training). Default 0. Ignored when
#'   \code{n_sites} is \code{NULL}.
#'
#' @return A named list of class \code{"mms_cost"}:
#'   \item{total_cost}{Total study cost (\code{sampling_cost + site_cost})}
#'   \item{sampling_cost}{\code{n * cost_per_sample}}
#'   \item{site_cost}{\code{n_sites * cost_per_site} (0 when \code{n_sites}
#'     is \code{NULL})}
#'   \item{n}{Total samples as supplied}
#'   \item{n_sites}{\code{n_sites} as supplied, or \code{NULL}}
#'   \item{cost_per_sample}{\code{cost_per_sample} as supplied}
#'   \item{cost_per_site}{\code{cost_per_site} as supplied}
#'
#' @section Equations and sources:
#' This is a bookkeeping identity, not a statistical formula:
#' \deqn{total = n \times cost\_per\_sample + n_{sites} \times cost\_per\_site.}
#' The MMS-SD Study Design Workshop
#' (\url{https://mrc-ide.github.io/MMS-SD_workshop/}) does not give a cost
#' equation; it frames cost as the feasibility side of the design
#' trade-off -- Module 7 "Designing a study for multiple end-points"
#' ("balance power vs feasibility of target sample sizes") and Module 2
#' ("Feedback from sample size calculation: what can you afford? What is
#' logistically feasible?", lecture slides p. 10). This function just
#' operationalises that trade-off so a design's \code{n} / \code{n_sites}
#' can be costed directly.
#'
#' @references
#' MMS-SD Study Design Workshop, Modules 2 (sample size / feasibility) and
#' 7 (designing for multiple end-points).
#' \url{https://mrc-ide.github.io/MMS-SD_workshop/}
#'
#' @export
#'
#' @examples
#' # 323 samples at $15 each
#' estimate_cost(n = 323, cost_per_sample = 15)
#'
#' # With 10 sites and $200 fixed cost per site
#' estimate_cost(n = 323, cost_per_sample = 15, n_sites = 10, cost_per_site = 200)
#'
#' # Pipe directly from design_precision (once that function is available)
#' \dontrun{
#' des <- design_precision(prevalence = 0.15, moe = 0.05)
#' estimate_cost(n = des$n, cost_per_sample = 12)
#' }
estimate_cost <- function(
  n,
  cost_per_sample,
  n_sites       = NULL,
  cost_per_site = 0
) {

  # ---- validate n ----
  if (length(n) != 1 || !is.numeric(n))
    stop("`n` must be a single positive integer (got class `", class(n)[1],
         "`, length ", length(n), ").")
  if (!is.finite(n) || n <= 0 || n != floor(n))
    stop("`n` must be a finite positive integer (got ", n, ").")

  # ---- validate cost_per_sample ----
  if (length(cost_per_sample) != 1 || !is.numeric(cost_per_sample))
    stop("`cost_per_sample` must be a single non-negative number (got class `",
         class(cost_per_sample)[1], "`, length ", length(cost_per_sample), ").")
  if (!is.finite(cost_per_sample) || cost_per_sample < 0)
    stop("`cost_per_sample` must be a finite non-negative number (got ",
         cost_per_sample, ").")

  # ---- validate n_sites ----
  if (!is.null(n_sites)) {
    if (length(n_sites) != 1 || !is.numeric(n_sites))
      stop("`n_sites` must be a single positive integer or NULL (got class `",
           class(n_sites)[1], "`, length ", length(n_sites), ").")
    if (!is.finite(n_sites) || n_sites <= 0 || n_sites != floor(n_sites))
      stop("`n_sites` must be a finite positive integer (got ", n_sites, ").")
    if (n_sites > n)
      stop("`n_sites` (", n_sites, ") cannot exceed `n` (", n,
           "): that would mean fewer than 1 sample per site on average.")
  }

  # ---- validate cost_per_site ----
  if (length(cost_per_site) != 1 || !is.numeric(cost_per_site))
    stop("`cost_per_site` must be a single non-negative number (got class `",
         class(cost_per_site)[1], "`, length ", length(cost_per_site), ").")
  if (!is.finite(cost_per_site) || cost_per_site < 0)
    stop("`cost_per_site` must be a finite non-negative number (got ",
         cost_per_site, ").")

  # ---- compute ----
  sampling_cost <- n * cost_per_sample
  site_cost     <- if (!is.null(n_sites)) n_sites * cost_per_site else 0

  structure(
    list(
      total_cost      = sampling_cost + site_cost,
      sampling_cost   = sampling_cost,
      site_cost       = site_cost,
      n               = n,
      n_sites         = n_sites,
      cost_per_sample = cost_per_sample,
      cost_per_site   = cost_per_site
    ),
    class = "mms_cost"
  )
}

#' @export
print.mms_cost <- function(x, ...) {
  fmt <- function(v) format(round(v, 2), big.mark = ",", nsmall = 2)
  cat(sprintf("Total cost: %s\n", fmt(x$total_cost)))
  cat(sprintf("  Sampling (%s samples x %s each): %s\n",
              format(x$n, big.mark = ","),
              fmt(x$cost_per_sample),
              fmt(x$sampling_cost)))
  if (!is.null(x$n_sites) && x$site_cost > 0)
    cat(sprintf("  Site fixed costs (%s sites x %s each): %s\n",
                format(x$n_sites, big.mark = ","),
                fmt(x$cost_per_site),
                fmt(x$site_cost)))
  invisible(x)
}
