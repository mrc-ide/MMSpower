#' Estimate total study cost
#'
#' @description
#' Costs an MMS study from its design choices, as a variable cost plus a
#' per-site cost:
#'
#' \itemize{
#'   \item \strong{Variable cost} = total samples enrolled times the cost
#'     per sample (collection, laboratory testing, consumables, data
#'     management).
#'   \item \strong{Site cost} = the sum over sites of (fixed cost per site +
#'     transport cost per site). The fixed cost (training, equipment,
#'     administration) and the transport cost may each be given once for all
#'     sites or as one value per site.
#' }
#'
#' Pairs with the design functions: take \code{n} and \code{n_sites} from a
#' design function's result. To explore trade-offs, change \code{n} or
#' \code{n_sites} and call it again. To cost several regions, call
#' \code{estimate_cost()} once per region and add the totals. All monetary
#' inputs and outputs are in one currency unit.
#'
#' @param n Whole number. Total number of samples to be enrolled: \code{$n}
#'   from any design function (for \code{design_difference()}, the total
#'   across both groups).
#' @param cost_per_sample Non-negative number. Variable cost per enrolled
#'   sample.
#' @param n_sites Optional whole number. Total number of sites to cost:
#'   \code{$n_sites} from a design function (for \code{design_difference()},
#'   \code{sum(n_sites)}). \code{NULL} (default) costs the samples only.
#' @param fixed_cost_per_site Non-negative fixed cost per site: one number
#'   for every site, or one value per site (length \code{n_sites}). Default
#'   0.
#' @param transport_cost_per_site Non-negative transport cost per site: one
#'   number for every site, or one value per site. Default 0.
#' @param budget Optional positive number. Total money available.
#'
#' @return A list with:
#'   \item{total_cost}{\code{variable_cost + site_cost}.}
#'   \item{variable_cost}{\code{n * cost_per_sample}.}
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
#' # Variable cost only
#' estimate_cost(n = 323, cost_per_sample = 50)
#'
#' # Same costs at every site (5,000 fixed + 1,000 transport per site)
#' estimate_cost(n = 600, cost_per_sample = 50, n_sites = 12,
#'               fixed_cost_per_site = 5000, transport_cost_per_site = 1000)
#'
#' # A different fixed cost at each of 3 sites, checked against a budget
#' estimate_cost(n = 50, cost_per_sample = 10, n_sites = 3,
#'               fixed_cost_per_site = c(100, 200, 300), budget = 1500)
#'
#' # Cost a design
#' des <- design_precision(prevalence = 0.15, moe = 0.05, n_per_site = 30, icc = 0.05)
#' estimate_cost(n = des$n, cost_per_sample = 50, n_sites = des$n_sites,
#'               fixed_cost_per_site = 5000, transport_cost_per_site = 1200)
#'
#' # Cost a two-group design: design_difference() gives one n_sites per
#' # group, so add them up
#' des2 <- design_difference(0.10, 0.20, n_per_site = 30, icc = 0.05)
#' estimate_cost(n = des2$n, cost_per_sample = 50, n_sites = sum(des2$n_sites),
#'               fixed_cost_per_site = 5000)
estimate_cost <- function(n,
                          cost_per_sample,
                          n_sites                 = NULL,
                          fixed_cost_per_site     = 0,
                          transport_cost_per_site = 0,
                          budget                  = NULL) {

  .check_count(n, "n")
  .check_number(cost_per_sample, "cost_per_sample", 0)
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

  variable_cost <- n * cost_per_sample
  site_cost <- if (is.null(n_sites)) 0 else
    sum(rep_len(fixed_cost_per_site, n_sites) + rep_len(transport_cost_per_site, n_sites))
  total_cost <- variable_cost + site_cost

  list(
    total_cost       = total_cost,
    variable_cost    = variable_cost,
    site_cost        = site_cost,
    budget_remaining = if (!is.null(budget)) budget - total_cost
  )
}
