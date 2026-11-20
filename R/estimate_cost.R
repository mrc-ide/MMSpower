#' Estimate total study cost
#'
#' @description
#' Costs an MMS study from its design choices, as a cost for each stage a
#' sample passes through plus a per-site cost:
#'
#' \itemize{
#'   \item \strong{Stage costs}: for each of up to three stages, the number
#'     of people or samples reaching that stage times that stage's cost per
#'     person or sample. The stages are:
#'     \describe{
#'       \item{screened}{people screened for infection at enrolment: consent,
#'         interview, a rapid test or microscopy.}
#'       \item{collected}{people a sample is collected from: a dried blood
#'         spot or blood draw, storage and shipping.}
#'       \item{sequenced}{samples taken through the molecular laboratory:
#'         extraction, PCR or sequencing, data management.}
#'     }
#'   \item \strong{Site cost} = the sum over sites of (fixed cost per site +
#'     transport cost per site). The fixed cost (training, equipment,
#'     administration) and the transport cost may each be given once for all
#'     sites or as one value per site.
#' }
#'
#' The stages are nested: everyone collected was first screened, and every
#' sample sequenced was first collected. So \code{n_collected} is a subset of
#' \code{n_screened}, not in addition to it (for example, only those who test
#' positive), and \code{n_sequenced} is a subset of \code{n_collected} (for
#' example, after lost or failed samples). A person is counted at every stage
#' they reach, and each stage's cost is the cost of that stage alone, not of
#' everything up to it. Give only the stages a study has; at least one is
#' needed. A study costed with a single all-in cost per sample can give that
#' as \code{n_sequenced} and \code{cost_per_sequenced}.
#'
#' Pairs with the design functions: their \code{$n} is the number of samples
#' needed at analysis, so it is \code{n_sequenced}; scale it up for losses
#' and for suspected cases that turn out negative to get \code{n_collected}
#' and \code{n_screened}. Take \code{n_sites} from the design function's
#' result too. To cost several regions, call \code{estimate_cost()} once per
#' region and add the totals. All monetary inputs and outputs are in one
#' currency unit.
#'
#' @param n_screened,n_collected,n_sequenced Optional whole numbers. The
#'   number of people or samples reaching each stage (see Description). Each
#'   is no larger than the one before it. For \code{n_sequenced}, \code{$n}
#'   from any design function (for \code{design_difference()}, the total
#'   across both groups).
#' @param cost_per_screened,cost_per_collected,cost_per_sequenced
#'   Non-negative numbers. The cost of each stage per person or sample
#'   reaching it. Give each one with its \code{n_} argument, or neither.
#' @param n_sites Optional whole number. Total number of sites to cost:
#'   \code{$n_sites} from a design function (for \code{design_difference()},
#'   \code{sum(n_sites)}). \code{NULL} (default) costs the stages only.
#' @param fixed_cost_per_site Non-negative fixed cost per site: one number
#'   for every site, or one value per site (length \code{n_sites}). Default
#'   0.
#' @param transport_cost_per_site Non-negative transport cost per site: one
#'   number for every site, or one value per site. Default 0.
#' @param budget Optional positive number. Total money available.
#'
#' @return A list with:
#'   \item{total_cost}{\code{variable_cost + site_cost}.}
#'   \item{variable_cost}{The sum of \code{stage_cost}.}
#'   \item{stage_cost}{Named vector of the cost of each stage
#'     (\code{screened}, \code{collected}, \code{sequenced}): its number
#'     times its cost, or 0 for a stage not given.}
#'   \item{site_cost}{Sum over sites of fixed plus transport cost; 0 when
#'     \code{n_sites} is \code{NULL}.}
#'   \item{budget_remaining}{\code{budget - total_cost} (negative if over
#'     budget), or \code{NULL} when no \code{budget} is given.}
#'
#' @references
#' MMS-SD Study Design Workshop, budget-officer activity; Modules 2 and 7.
#' \url{https://mrc-ide.github.io/MMS-SD_workshop/}
#'
#' @export
#'
#' @examples
#' # One all-in cost per sample
#' estimate_cost(n_sequenced = 323, cost_per_sequenced = 50)
#'
#' # Three stages: 1,950 screened, the 780 positives sampled, 700 sequenced
#' # after losses; 12 sites at 5,000 fixed + 1,000 transport each
#' estimate_cost(n_screened  = 1950, cost_per_screened  = 2,
#'               n_collected = 780,  cost_per_collected = 15,
#'               n_sequenced = 700,  cost_per_sequenced = 40,
#'               n_sites = 12, fixed_cost_per_site = 5000,
#'               transport_cost_per_site = 1000)
#'
#' # A different fixed cost at each of 3 sites, checked against a budget
#' estimate_cost(n_sequenced = 50, cost_per_sequenced = 10, n_sites = 3,
#'               fixed_cost_per_site = c(100, 200, 300), budget = 1500)
#'
#' # Cost a design, allowing 10% of samples to be lost and 40% of suspected
#' # cases to be positive
#' des <- design_precision(prevalence = 0.15, moe = 0.05, n_per_site = 30, icc = 0.05)
#' n_collected <- ceiling(des$n / 0.9)
#' estimate_cost(n_screened  = ceiling(n_collected / 0.4), cost_per_screened  = 2,
#'               n_collected = n_collected,                cost_per_collected = 15,
#'               n_sequenced = des$n,                      cost_per_sequenced = 40,
#'               n_sites = des$n_sites, fixed_cost_per_site = 5000,
#'               transport_cost_per_site = 1200)
#'
#' # Cost a two-group design: design_difference() gives one n_sites per
#' # group, so add them up
#' des2 <- design_difference(0.10, 0.20, n_per_site1 = 30, icc = 0.05)
#' estimate_cost(n_sequenced = des2$n, cost_per_sequenced = 50,
#'               n_sites = sum(des2$n_sites), fixed_cost_per_site = 5000)
estimate_cost <- function(n_screened              = NULL,
                          cost_per_screened       = NULL,
                          n_collected             = NULL,
                          cost_per_collected      = NULL,
                          n_sequenced             = NULL,
                          cost_per_sequenced      = NULL,
                          n_sites                 = NULL,
                          fixed_cost_per_site     = 0,
                          transport_cost_per_site = 0,
                          budget                  = NULL) {

  # stages: each needs its number and its cost together, or neither
  stages <- c("screened", "collected", "sequenced")
  n_stage <- list(n_screened, n_collected, n_sequenced)
  cost_stage <- list(cost_per_screened, cost_per_collected, cost_per_sequenced)
  for (i in seq_along(stages)) {
    n_nm <- paste0("n_", stages[i])
    cost_nm <- paste0("cost_per_", stages[i])
    if (is.null(n_stage[[i]]) != is.null(cost_stage[[i]]))
      stop("Give `", n_nm, "` and `", cost_nm, "` together, or neither.",
           call. = FALSE)
    .check_count(n_stage[[i]], n_nm)
    if (!is.null(cost_stage[[i]])) .check_number(cost_stage[[i]], cost_nm, 0)
  }
  given <- !vapply(n_stage, is.null, logical(1))
  if (!any(given))
    stop("Give at least one stage: `n_screened`, `n_collected` or `n_sequenced`, ",
         "each with its cost.", call. = FALSE)
  # nested stages: no stage can have more than the one before it
  n_given <- unlist(n_stage[given])
  if (any(diff(n_given) > 0))
    stop("Each stage is a subset of the one before it, so `n_screened` >= ",
         "`n_collected` >= `n_sequenced`.", call. = FALSE)

  .check_count(n_sites, "n_sites")
  if (!is.null(budget)) .check_number(budget, "budget", 0, lower_open = TRUE)
  # per-site costs: one value for every site, or one per site
  for (nm in c("fixed_cost_per_site", "transport_cost_per_site")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) < 1 || !all(is.finite(v)) || any(v < 0))
      stop("`", nm, "` must be non-negative numbers.", call. = FALSE)
    if (is.null(n_sites) && any(v != 0))
      stop("`", nm, "` needs `n_sites`.", call. = FALSE)
    if (!is.null(n_sites) && !length(v) %in% c(1, n_sites))
      stop("`", nm, "` must have one value, or one per site.", call. = FALSE)
  }

  stage_cost <- stats::setNames(
    mapply(function(n, cost) if (is.null(n)) 0 else n * cost, n_stage, cost_stage),
    stages)
  variable_cost <- sum(stage_cost)
  site_cost <- if (is.null(n_sites)) 0 else
    sum(rep_len(fixed_cost_per_site, n_sites) + rep_len(transport_cost_per_site, n_sites))
  total_cost <- variable_cost + site_cost

  list(
    total_cost       = total_cost,
    variable_cost    = variable_cost,
    stage_cost       = stage_cost,
    site_cost        = site_cost,
    budget_remaining = if (!is.null(budget)) budget - total_cost
  )
}
