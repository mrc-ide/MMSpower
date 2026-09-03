#' Estimate total study cost
#'
#' @description
#' Costs an MMS study from its design choices, using the budget structure
#' taught in the MMS-SD workshop budget-officer activity:
#'
#' \itemize{
#'   \item \strong{Total variable cost} = total samples enrolled x cost per
#'     sample (collection, laboratory testing, consumables, data
#'     management).
#'   \item \strong{Total fixed cost} = for each region, number of health
#'     facilities x (fixed cost per facility + transport cost per facility).
#'     Both the fixed cost per facility (training, equipment,
#'     administration) and the transport cost per facility may be given once
#'     for all regions or as a per-region named vector.
#' }
#'
#' Provide the design: total samples, health facilities per region, and the
#' cost rates, and \code{estimate_cost()} returns the fixed / variable /
#' total breakdown, a per-region table, and (if you pass a \code{budget})
#' whether you are within it. To explore trade-offs, change \code{n} or
#' \code{n_sites} and call it again.
#'
#' Pairs with the design functions: for a single-region design pass
#' \code{design_precision()} / \code{design_detection()}'s \code{n} and
#' \code{n_sites} straight in; for a multi-region budget, split their
#' \code{n_sites} across regions yourself (a named \code{n_sites} vector).
#'
#' Please note that all monetary inputs and outputs are in one currency
#' unit; the function does no conversion.
#'
#' @param n Positive integer. Total number of samples to be enrolled. All
#'   five design functions report this as \code{$n}:
#'   \code{design_precision()}, \code{design_threshold()},
#'   \code{design_detection()}, \code{design_difference()} and
#'   \code{design_trend()} (for the last two, \code{$n} is the study total --
#'   summed across both groups / all timepoints).
#' @param cost_per_sample Non-negative number. Variable cost per enrolled
#'   sample.
#' @param n_sites Number of health facilities. One of:
#'   \itemize{
#'     \item a single positive integer -- one region;
#'     \item a named positive-integer vector -- the facility count in each
#'       region, e.g. \code{c(North = 10, South = 8, West = 12)};
#'     \item \code{NULL} (the default) -- no health facilities are costed,
#'       so the result is the sample (variable) cost only.
#'   }
#'   A single-region count can be taken straight from a design function's
#'   \code{$n_sites} (reported per group by \code{design_difference()} and
#'   per timepoint by \code{design_trend()}).
#' @param fixed_cost_per_site Non-negative. Fixed cost per health facility
#'   (training, equipment, administration). Either a single number applied to
#'   every region, or a named vector with the same region names as
#'   \code{n_sites}, e.g. \code{c(North = 5000, South = 4000)}. Default 0.
#' @param transport_cost_per_site Non-negative. Transport cost per health
#'   facility. Either a single number applied to every region, or a named
#'   vector with the same region names as \code{n_sites}, e.g.
#'   \code{c(North = 1000, South = 1500, West = 800)}. Default 0.
#' @param budget Optional positive number. Total money available. When
#'   supplied, the result reports \code{budget_remaining} and
#'   \code{over_budget}.
#'
#' @return A named list of class \code{"mms_cost"}:
#'   \item{total_cost}{\code{total_fixed_cost + total_variable_cost}.}
#'   \item{total_variable_cost}{\code{n * cost_per_sample}.}
#'   \item{total_fixed_cost}{\eqn{\sum_{r} n_{sites,r} \times
#'     (fixed\_cost\_per\_site_r + transport\_cost\_per\_site_r)}; 0 when
#'     \code{n_sites} is \code{NULL}.}
#'   \item{by_region}{Data frame with one row per region -- \code{region},
#'     \code{n_sites}, \code{fixed_cost_per_site},
#'     \code{transport_cost_per_site}, \code{fixed_subtotal}
#'     (\code{n_sites * fixed_cost_per_site}), \code{transport_subtotal}
#'     (\code{n_sites * transport_cost_per_site}), \code{region_total}
#'     (their sum). \code{NULL} when \code{n_sites} is \code{NULL}.}
#'   \item{n}{Total samples as supplied.}
#'   \item{n_sites}{Total health facilities across all regions (0 when
#'     \code{n_sites} is \code{NULL}).}
#'   \item{cost_per_sample, fixed_cost_per_site}{As supplied.}
#'   \item{budget}{As supplied, or \code{NULL}.}
#'   \item{budget_remaining}{\code{budget - total_cost} (negative if over);
#'     \code{NULL} when no \code{budget}.}
#'   \item{over_budget}{Logical: \code{total_cost > budget}; \code{NULL}
#'     when no \code{budget}.}
#'
#' @section Equations and sources:
#' Straight from the MMS-SD Study Design Workshop budget-officer group
#' activity (\url{https://mrc-ide.github.io/MMS-SD_workshop/}):
#' \deqn{\mathrm{Total\ Fixed} = \sum_{r} n_{HF,r} \times
#'   (\mathrm{fixed\ cost\ per\ HF}_r + \mathrm{transport\ cost\ per\ HF}_r)}
#' \deqn{\mathrm{Total\ Variable} = \mathrm{samples\ enrolled} \times
#'   \mathrm{cost\ per\ sample}}
#' The function frames cost as the feasibility side of the design
#' trade-off; it operationalises those formulae so a design's \code{n} /
#' \code{n_sites} can be costed and checked against a budget.
#'
#' @references
#' MMS-SD Study Design Workshop, budget-officer activity; Modules 2 and 7.
#' \url{https://mrc-ide.github.io/MMS-SD_workshop/}
#'
#' @seealso The design functions that produce the \code{n} / \code{n_sites}
#'   this costs: \code{design_precision()}, \code{design_threshold()},
#'   \code{design_detection()}, \code{design_difference()},
#'   \code{design_trend()}.
#'
#' @concept cost
#' @concept budget
#' @concept study design
#'
#' @export
#'
#' @examples
#' # Variable cost only
#' estimate_cost(n = 323, cost_per_sample = 50)
#'
#' # Single region, workshop rates ($5,000 fixed + $1,000 transport per HF)
#' estimate_cost(
#'   n = 600, cost_per_sample = 50,
#'   n_sites = 12,
#'   fixed_cost_per_site = 5000,
#'   transport_cost_per_site = 1000
#' )
#'
#' # Multi-region, per-region transport costs, checked against a budget
#' estimate_cost(
#'   n = 800, cost_per_sample = 50,
#'   n_sites = c(North = 10, South = 8, West = 12),
#'   fixed_cost_per_site = 5000,
#'   transport_cost_per_site = c(North = 1000, South = 1500, West = 800),
#'   budget = 300000
#' )
#'
#' # Fixed cost can also vary by region (e.g. higher setup in the North)
#' estimate_cost(
#'   n = 400, cost_per_sample = 50,
#'   n_sites = c(North = 6, South = 6),
#'   fixed_cost_per_site = c(North = 6000, South = 4500),
#'   transport_cost_per_site = c(North = 1000, South = 1500)
#' )
#'
#' # Pipe from a design function
#' \dontrun{
#' des <- design_precision(prevalence = 0.15, moe = 0.05, n_per_site = 30, icc = 0.05)
#' estimate_cost(n = des$n, cost_per_sample = 50, n_sites = des$n_sites,
#'               fixed_cost_per_site = 5000, transport_cost_per_site = 1200)
#' }
estimate_cost <- function(n,
                          cost_per_sample,
                          n_sites                 = NULL,
                          fixed_cost_per_site     = 0,
                          transport_cost_per_site = 0,
                          budget                  = NULL) {

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

  # ---- validate the per-facility costs (scalar or per-region vector) ----
  # Full region-name matching happens later, once we know the regions; here
  # just reject non-numeric / negative / non-finite values.
  for (nm in c("fixed_cost_per_site", "transport_cost_per_site")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) < 1 || !all(is.finite(v)) || any(v < 0))
      stop("`", nm, "` must be a finite non-negative number, or a named ",
           "non-negative vector of per-region costs.")
  }

  # ---- validate budget ----
  if (!is.null(budget) && (length(budget) != 1 || !is.numeric(budget) ||
      !is.finite(budget) || budget <= 0))
    stop("`budget` must be a single finite positive number or NULL (got ",
         if (length(budget) != 1) paste0("length = ", length(budget)) else budget, ").")

  # ---- variable cost ----
  # Coerce to double first: an integer n times an integer unit cost can
  # overflow to NA.
  total_variable_cost <- as.double(n) * cost_per_sample

  # ================================================================
  # No facilities: variable cost only.
  # ================================================================
  if (is.null(n_sites)) {
    if (any(transport_cost_per_site != 0) || any(fixed_cost_per_site != 0))
      warning("`fixed_cost_per_site` / `transport_cost_per_site` are ignored ",
              "because `n_sites` is NULL. Supply `n_sites` to include ",
              "per-facility fixed costs.")

    total_fixed_cost <- 0
    by_region        <- NULL
    n_sites_total    <- 0

  } else {
    # ---- validate n_sites ----
    if (!is.numeric(n_sites) || length(n_sites) < 1 || !all(is.finite(n_sites)) ||
        any(n_sites <= 0) || any(n_sites != floor(n_sites)))
      stop("`n_sites` must be a positive integer, or a named positive-integer ",
           "vector of per-region facility counts.")
    if (length(n_sites) > 1 && is.null(names(n_sites)))
      stop("A multi-region `n_sites` must be named, e.g. ",
           "c(North = 10, South = 8). Got an unnamed length-", length(n_sites),
           " vector.")

    regions <- if (is.null(names(n_sites))) "(all)" else names(n_sites)

    # ---- resolve per-region costs ----
    # For both fixed and transport: an unnamed single value is "same for
    # every region"; anything named (even length 1) is a per-region lookup
    # that must supply a value for each region in `n_sites`.
    fixed_vec     <- .resolve_region_cost(fixed_cost_per_site,
                                          "fixed_cost_per_site", regions)
    transport_vec <- .resolve_region_cost(transport_cost_per_site,
                                          "transport_cost_per_site", regions)

    n_sites_vec <- as.double(n_sites)

    fixed_subtotal     <- n_sites_vec * fixed_vec
    transport_subtotal <- n_sites_vec * transport_vec
    region_total       <- fixed_subtotal + transport_subtotal

    by_region <- data.frame(
      region                  = regions,
      n_sites                 = n_sites_vec,
      fixed_cost_per_site      = unname(fixed_vec),
      transport_cost_per_site  = unname(transport_vec),
      fixed_subtotal           = fixed_subtotal,
      transport_subtotal       = transport_subtotal,
      region_total             = region_total,
      row.names = NULL,
      stringsAsFactors = FALSE
    )

    total_fixed_cost <- sum(region_total)
    n_sites_total    <- sum(n_sites_vec)
  }

  total_cost <- total_fixed_cost + total_variable_cost

  budget_remaining <- if (!is.null(budget)) budget - total_cost else NULL
  over_budget      <- if (!is.null(budget)) total_cost > budget    else NULL

  structure(
    list(
      total_cost          = total_cost,
      total_variable_cost = total_variable_cost,
      total_fixed_cost    = total_fixed_cost,
      by_region           = by_region,
      n                   = n,
      n_sites             = n_sites_total,
      cost_per_sample     = cost_per_sample,
      fixed_cost_per_site  = fixed_cost_per_site,
      budget              = budget,
      budget_remaining    = budget_remaining,
      over_budget         = over_budget
    ),
    class = "mms_cost"
  )
}

#' @export
print.mms_cost <- function(x, ...) {
  # scientific = FALSE: otherwise format() renders round large numbers
  # (e.g. 1e7) in scientific notation, defeating the big.mark grouping.
  money <- function(v) format(round(v, 2), big.mark = ",", nsmall = 2,
                              scientific = FALSE)
  int   <- function(v) format(v, big.mark = ",", scientific = FALSE)

  cat("Study cost estimate\n")
  cat(sprintf("  Variable:  %s samples x %s        = %s\n",
              int(x$n), money(x$cost_per_sample), money(x$total_variable_cost)))

  if (!is.null(x$by_region)) {
    cat(sprintf("  Fixed:     %s health facilities         = %s\n",
                int(x$n_sites), money(x$total_fixed_cost)))
    br <- x$by_region
    for (i in seq_len(nrow(br))) {
      cat(sprintf("    %-10s %s HF x (%s + %s) = %s\n",
                  br$region[i], int(br$n_sites[i]),
                  money(br$fixed_cost_per_site[i]),
                  money(br$transport_cost_per_site[i]),
                  money(br$region_total[i])))
    }
  }

  cat(sprintf("  %s\n", strrep("-", 48)))
  cat(sprintf("  Total                                = %s\n", money(x$total_cost)))

  if (!is.null(x$budget)) {
    if (isTRUE(x$over_budget)) {
      cat(sprintf("  Budget %s -> OVER by %s\n",
                  money(x$budget), money(-x$budget_remaining)))
    } else {
      cat(sprintf("  Budget %s -> within by %s\n",
                  money(x$budget), money(x$budget_remaining)))
      if (x$cost_per_sample > 0)
        cat(sprintf("    slack ~= %s more samples\n",
                    int(floor(x$budget_remaining / x$cost_per_sample))))
      if (!is.null(x$by_region)) {
        per_hf <- mean(x$by_region$fixed_cost_per_site) +
                  mean(x$by_region$transport_cost_per_site)
        if (per_hf > 0)
          cat(sprintf("    slack ~= %s more health facilities (at average per-HF cost)\n",
                      int(floor(x$budget_remaining / per_hf))))
      }
    }
  }

  invisible(x)
}

# Expand a per-facility cost input to one value per region.
#
#   `x` unnamed, length 1 -- recycled to every region.
#   `x` named             -- looked up by region name; every region in
#                            `regions` must be present (extra names ignored).
#
# Numeric / finite / non-negative checks happen in `estimate_cost()` before
# this is called, so only the shape / name matching is enforced here.
.resolve_region_cost <- function(x, name, regions) {
  if (length(x) == 1 && is.null(names(x)))
    return(stats::setNames(rep(as.double(x), length(regions)), regions))

  if (is.null(names(x)))
    stop("`", name, "` has ", length(x), " values but must be named to ",
         "match `n_sites` (e.g. c(",
         paste0(regions, " = ...", collapse = ", "),
         ")); or give a single number to apply to every region.")

  missing_r <- setdiff(regions, names(x))
  if (length(missing_r) > 0)
    stop("`", name, "` is missing a value for region(s): ",
         paste(missing_r, collapse = ", "), ".")

  stats::setNames(as.double(x[regions]), regions)
}
