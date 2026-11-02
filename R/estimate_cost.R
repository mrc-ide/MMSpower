#' Estimate total study cost
#'
#' @description
#' Costs an MMS study from its design choices, as a variable cost plus a
#' fixed cost:
#'
#' \itemize{
#'   \item \strong{Total variable cost} = total samples enrolled times the
#'     cost per sample (collection, laboratory testing, consumables, data
#'     management).
#'   \item \strong{Total fixed cost} = the sum over sites of
#'     (fixed cost per site + transport cost per site). The fixed
#'     cost (training, equipment, administration) and the transport cost may
#'     each be given once for all sites or as one value per site.
#'     When the costs are the same at every site, this is the number of
#'     sites times (fixed cost + transport cost).
#' }
#'
#' Provide the design: total samples, number of sites, and the
#' costs, and \code{estimate_cost()} returns the fixed, variable and
#' total costs, a per-site table, and (if you pass a \code{budget})
#' whether you are within it. To explore trade-offs, change \code{n} or
#' \code{n_sites} and call it again.
#'
#' Pairs with the design functions: take \code{n} and \code{n_sites} from a
#' design function's result, as described under each argument below. To
#' cost several regions, call \code{estimate_cost()} once per region and add
#' the totals.
#'
#' All monetary inputs and outputs are in one currency unit; the function
#' does no conversion.
#'
#' @param n Positive integer. Total number of samples to be enrolled. All
#'   five design functions report this as \code{$n}:
#'   \code{design_precision()}, \code{design_threshold()},
#'   \code{design_detection()}, \code{design_difference()} and
#'   \code{design_trend()} (for the last two, \code{$n} is the study total,
#'   summed across both groups or all timepoints). With sites, every site
#'   gets the same whole number of people, so the number actually collected
#'   can be more than \code{$n}. To budget for what will be collected, pass
#'   \code{n_sites * n_per_site}, times the number of groups or timepoints
#'   for \code{design_difference()} and \code{design_trend()}.
#' @param cost_per_sample Non-negative number. Variable cost per enrolled
#'   sample.
#' @param n_sites Positive integer. Total number of sites to cost. Most
#'   design functions report it as \code{$n_sites}, but
#'   \code{design_difference()} reports it per group, so pass
#'   \code{2 * n_sites}, and \code{design_trend()} reports it per timepoint,
#'   so pass the total for the study, which depends on whether the same
#'   sites are visited at every timepoint. \code{NULL} (the default) means
#'   no sites are costed, so the result is the sample (variable) cost
#'   only.
#' @param fixed_cost_per_site Non-negative number or vector. Fixed cost per
#'   site (training, equipment, administration). Either a single
#'   number applied to every site, or one value per site (length
#'   \code{n_sites}), for example \code{c(100, 200, 300)} for 3 sites.
#'   Values are matched to sites by position; names are ignored.
#'   Default 0.
#' @param transport_cost_per_site Non-negative number or vector. Transport
#'   cost per site. Either a single number applied to every
#'   site, or one value per site (length \code{n_sites}). Values are
#'   matched to sites by position; names are ignored. Default 0.
#' @param budget Optional positive number. Total money available. When
#'   supplied, the result reports \code{budget_remaining} and
#'   \code{over_budget}.
#'
#' @return A named list of class \code{"mms_cost"}:
#'   \item{total_cost}{\code{total_fixed_cost + total_variable_cost}.}
#'   \item{total_variable_cost}{\code{n * cost_per_sample}.}
#'   \item{total_fixed_cost}{\eqn{\sum_{i} (fixed\_cost\_per\_site_i +
#'     transport\_cost\_per\_site_i)} over the \code{n_sites} sites; 0
#'     when \code{n_sites} is \code{NULL}.}
#'   \item{by_site}{Data frame with one row per site and columns
#'     \code{site} (1 to \code{n_sites}), \code{fixed_cost_per_site},
#'     \code{transport_cost_per_site} and \code{site_total} (their sum).
#'     \code{NULL} when \code{n_sites} is \code{NULL}.}
#'   \item{n}{Total samples, as supplied (as a plain number).}
#'   \item{n_sites}{Number of sites as supplied (as a plain number; 0 when
#'     \code{n_sites} is \code{NULL}).}
#'   \item{cost_per_sample}{As supplied (as a plain number).}
#'   \item{fixed_cost_per_site, transport_cost_per_site}{As supplied (single
#'     number or per-site vector).}
#'   \item{budget}{As supplied (as a plain number), or \code{NULL}.}
#'   \item{budget_remaining}{\code{budget - total_cost}, worked out in whole
#'     cents (negative if over), so rounding error in decimal costs never
#'     shows as a difference. \code{NULL} when no \code{budget}.}
#'   \item{over_budget}{Logical: \code{TRUE} when \code{budget_remaining}
#'     is below 0; \code{NULL} when no \code{budget}.}
#'
#' @references
#' MMS-SD Study Design Workshop, budget-officer activity; Modules 2 and 7.
#' \url{https://mrc-ide.github.io/MMS-SD_workshop/}
#'
#' @seealso The design functions that produce the \code{n} and
#'   \code{n_sites} this costs: \code{design_precision()},
#'   \code{design_threshold()}, \code{design_detection()},
#'   \code{design_difference()}, \code{design_trend()}.
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
#' # Same costs at every site (5,000 fixed + 1,000 transport per site)
#' estimate_cost(
#'   n = 600, cost_per_sample = 50,
#'   n_sites = 12,
#'   fixed_cost_per_site = 5000,
#'   transport_cost_per_site = 1000
#' )
#'
#' # A different fixed cost at each of 3 sites, checked against a budget
#' estimate_cost(
#'   n = 50, cost_per_sample = 10,
#'   n_sites = 3,
#'   fixed_cost_per_site = c(100, 200, 300),
#'   budget = 1500
#' )
#'
#' # Feed a design function's output in, budgeting for the people collected
#' \dontrun{
#' des <- design_precision(prevalence = 0.15, moe = 0.05, n_per_site = 30, icc = 0.05)
#' estimate_cost(n = des$n_sites * des$n_per_site, cost_per_sample = 50,
#'               n_sites = des$n_sites,
#'               fixed_cost_per_site = 5000, transport_cost_per_site = 1200)
#' }
estimate_cost <- function(n,
                          cost_per_sample,
                          n_sites                 = NULL,
                          fixed_cost_per_site     = 0,
                          transport_cost_per_site = 0,
                          budget                  = NULL) {

  # Numbers in messages are printed in full (100000, not 1e+05).
  num <- function(v) format(v, scientific = FALSE, digits = 15)
  # A value a hair off an integer (100 * 1.1 = 110.00000000000001) prints
  # as one, so the message says so.
  # Only for a whole number off by rounding error (100 * 1.1 =
  # 110.00000000000001), not for values like 1e10 + 3.7 or 1e-12
  near_int <- function(v, nm)
    if (is.finite(v) && v != round(v) &&
        abs(v - round(v)) <= 8 * .Machine$double.eps * abs(v))
      paste0(", which is not exactly an integer; use round(", nm, ")") else ""

  # ---- validate n ----
  if (length(n) != 1 || !is.numeric(n))
    stop("`n` must be a single positive integer (got class `", class(n)[1],
         "`, length ", length(n), ").")
  if (!is.finite(n) || n <= 0 || n != floor(n))
    stop("`n` must be a finite positive integer (got ", num(n), near_int(n, "n"), ").")

  # ---- validate cost_per_sample ----
  if (length(cost_per_sample) != 1 || !is.numeric(cost_per_sample))
    stop("`cost_per_sample` must be a single non-negative number (got class `",
         class(cost_per_sample)[1], "`, length ", length(cost_per_sample), ").")
  if (!is.finite(cost_per_sample) || cost_per_sample < 0)
    stop("`cost_per_sample` must be a finite non-negative number (got ",
         num(cost_per_sample), ").")

  # ---- validate the per-site costs (one value, or one per site) ----
  # The length is matched against n_sites later, once n_sites is validated;
  # here just reject non-numeric, negative or non-finite values.
  per_site_costs <- list(fixed_cost_per_site     = fixed_cost_per_site,
                         transport_cost_per_site = transport_cost_per_site)
  for (nm in names(per_site_costs)) {
    v <- per_site_costs[[nm]]
    if (!is.numeric(v) || length(v) < 1)
      stop("`", nm, "` must be a finite non-negative number, or a vector of ",
           "non-negative costs with one value per site (got class `",
           class(v)[1], "`, length ", length(v), ").")
    # Name the first bad value, so it can be found in a long vector.
    if (!all(is.finite(v)) || any(v < 0)) {
      i <- which(!is.finite(v) | v < 0)[1]
      stop("`", nm, "` must be a finite non-negative number, or a vector of ",
           "non-negative costs with one value per site (",
           if (length(v) == 1) paste0("got ", num(v)) else
             paste0("found ", nm, "[", i, "] = ", num(v[i])), ").")
    }
  }

  # ---- validate budget ----
  # The type is checked first, so text such as "100" is reported by its
  # class rather than printed as if it were the number 100.
  if (!is.null(budget)) {
    if (length(budget) != 1 || !is.numeric(budget))
      stop("`budget` must be a single positive number or NULL (got class `",
           class(budget)[1], "`, length ", length(budget), ").")
    if (!is.finite(budget) || budget <= 0)
      stop("`budget` must be a finite positive number (got ", num(budget), ").")
  }

  # From here on these are plain numbers: a name or a 1 x 1 matrix would
  # otherwise be carried into the totals and the result.
  n               <- as.numeric(n)
  cost_per_sample <- as.numeric(cost_per_sample)
  if (!is.null(budget)) budget <- as.numeric(budget)

  # ---- variable cost ----
  # Coerce to double first: an integer n times an integer unit cost can
  # overflow to NA.
  total_variable_cost <- as.double(n) * cost_per_sample

  # ================================================================
  # No sites: variable cost only.
  # ================================================================
  if (is.null(n_sites)) {
    if (any(transport_cost_per_site != 0) || any(fixed_cost_per_site != 0))
      warning("`fixed_cost_per_site` and `transport_cost_per_site` are ignored ",
              "because `n_sites` is NULL. Supply `n_sites` to include ",
              "per-site costs.")

    total_fixed_cost <- 0
    by_site          <- NULL
    n_sites_out      <- 0

  } else {
    # ---- validate n_sites ----
    # The region hint only when several values were given, the likely sign
    # of one count per region.
    if (length(n_sites) != 1 || !is.numeric(n_sites))
      stop("`n_sites` must be a single positive integer or NULL (got class `",
           class(n_sites)[1], "`, length ", length(n_sites), ").",
           if (length(n_sites) > 1)
             paste0(" To cost several regions, call estimate_cost() once per ",
                    "region and add the totals."))
    if (!is.finite(n_sites) || n_sites <= 0 || n_sites != floor(n_sites))
      stop("`n_sites` must be a finite positive integer (got ", num(n_sites),
           near_int(n_sites, "n_sites"), ").")
    n_sites <- as.numeric(n_sites)

    # ---- match each per-site cost to n_sites ----
    for (nm in names(per_site_costs)) {
      len <- length(per_site_costs[[nm]])
      if (len != 1 && len != n_sites)
        stop("`", nm, "` has ", len, " values but `n_sites` is ", num(n_sites),
             ". Give a single number for every site, or one value per ",
             "site (", num(n_sites), if (n_sites == 1) " value)." else " values).")
    }

    # Mirror of the NULL-n_sites warning above: sites given but no costs.
    if (all(fixed_cost_per_site == 0) && all(transport_cost_per_site == 0))
      warning("`n_sites` supplied but both per-site costs are 0. ",
              "Supply `fixed_cost_per_site` or `transport_cost_per_site` to ",
              "cost the sites.")

    # A single number is recycled to every site.
    fixed_vec     <- rep_len(as.double(fixed_cost_per_site), n_sites)
    transport_vec <- rep_len(as.double(transport_cost_per_site), n_sites)

    by_site <- data.frame(
      site                    = seq_len(n_sites),
      fixed_cost_per_site     = fixed_vec,
      transport_cost_per_site = transport_vec,
      site_total              = fixed_vec + transport_vec
    )

    total_fixed_cost <- sum(by_site$site_total)
    n_sites_out      <- n_sites
  }

  total_cost <- total_fixed_cost + total_variable_cost

  # Money is compared in whole cents, as it is printed. Whole numbers are
  # exact in R (up to 2^53 cents, about 9e13), so the tiny rounding errors
  # in decimal costs (7 * 1.1 is not exactly 7.7) never change the verdict
  # or show as a difference. (round(x, 2) itself is not exact above about
  # 1e12.)
  budget_remaining <- NULL
  over_budget      <- NULL
  if (!is.null(budget)) {
    cents <- function(v) round(v * 100)
    budget_remaining <- (cents(budget) - cents(total_cost)) / 100
    over_budget      <- budget_remaining < 0
  }

  structure(
    list(
      total_cost = total_cost,
      total_variable_cost = total_variable_cost,
      total_fixed_cost = total_fixed_cost,
      by_site = by_site,
      n = n,
      n_sites = n_sites_out,
      cost_per_sample = cost_per_sample,
      fixed_cost_per_site = fixed_cost_per_site,
      transport_cost_per_site = transport_cost_per_site,
      budget = budget,
      budget_remaining = budget_remaining,
      over_budget = over_budget
    ),
    class = "mms_cost"
  )
}

#' @export
print.mms_cost <- function(x, ...) {
  # Always 2 decimals with comma grouping, never scientific notation, and
  # not affected by options(digits).
  # (+ 0 turns a -0, e.g. from a cost of -0, into 0.)
  money <- function(v) formatC(round(v, 2) + 0, format = "f", digits = 2, big.mark = ",")
  int   <- function(v) format(v, big.mark = ",", scientific = FALSE)
  # "1 site", "2 sites"
  count <- function(k, word) paste(int(k), if (k == 1) word else paste0(word, "s"))

  # Costs that differ by site are listed one line per site, up to
  # this many; the rest are counted in a note (all are in $by_site).
  max_sites_shown <- 10

  # Lines that end in an amount are kept as a label and an amount, then
  # printed with the labels padded to the widest label and the amounts
  # right-aligned, so the = signs and the decimal points each line up.
  # A line with no amount (NA) is printed as it is.
  lines <- list()
  add <- function(label, amount = NA_character_)
    lines[[length(lines) + 1]] <<- c(label, amount)

  add(sprintf("Variable:  %s x %s", count(x$n, "sample"), money(x$cost_per_sample)),
      money(x$total_variable_cost))

  if (!is.null(x$by_site)) {
    add(sprintf("Fixed:     %s", count(x$n_sites, "site")), money(x$total_fixed_cost))
    bs <- x$by_site
    same_costs <- length(unique(bs$fixed_cost_per_site)) == 1 &&
                  length(unique(bs$transport_cost_per_site)) == 1
    add("  (fixed + transport cost per site)")
    if (same_costs) {
      add(sprintf("  %s x (%s + %s)", count(x$n_sites, "site"),
                  money(bs$fixed_cost_per_site[1]),
                  money(bs$transport_cost_per_site[1])))
    } else {
      shown <- min(nrow(bs), max_sites_shown)
      for (i in seq_len(shown))
        add(sprintf("  site %-4s (%s + %s)", bs$site[i],
                    money(bs$fixed_cost_per_site[i]),
                    money(bs$transport_cost_per_site[i])),
            money(bs$site_total[i]))
      if (nrow(bs) > shown)
        add(sprintf("  ... %s more %s not shown (all %s are in $by_site)",
                    int(nrow(bs) - shown),
                    if (nrow(bs) - shown == 1) "site" else "sites",
                    int(nrow(bs))))
    }
  }

  add("Total", money(x$total_cost))

  has_amount <- vapply(lines, function(l) !is.na(l[2]), TRUE)
  label_w    <- max(nchar(vapply(lines[has_amount], `[`, "", 1)))
  amount_w   <- max(nchar(vapply(lines[has_amount], `[`, "", 2)))
  # The dashed line above the total is as wide as the widest line.
  rule_w     <- max(label_w + 3 + amount_w,
                    nchar(vapply(lines[!has_amount], `[`, "", 1)))

  cat("Study cost estimate\n")
  for (k in seq_along(lines)) {
    l <- lines[[k]]
    # The dashed line above the total spans the aligned lines.
    if (k == length(lines))
      cat("  ", strrep("-", rule_w), "\n", sep = "")
    if (is.na(l[2]))
      cat("  ", l[1], "\n", sep = "")
    else
      cat("  ", formatC(l[1], width = -label_w), " = ",
          formatC(l[2], width = amount_w), "\n", sep = "")
  }

  if (!is.null(x$budget)) {
    if (isTRUE(x$over_budget)) {
      cat(sprintf("  Budget %s -> over budget by %s\n",
                  money(x$budget), money(-x$budget_remaining)))
    } else {
      cat(sprintf("  Budget %s -> under budget by %s\n",
                  money(x$budget), money(x$budget_remaining)))

      # Re-express the leftover money as "how much more of one thing you
      # could buy with it" -- each line spends the whole remainder on that
      # item alone, so they are alternatives, not a combined plan. Skipped
      # when nothing is left (exactly on budget). The leftover is already
      # in whole cents. Rounding error in the division (0.29 / 0.01 =
      # 28.999999999999996) is undone by checking whether one more item
      # fits, up to rounding. An item the leftover cannot buy even one of is
      # left out.
      left <- x$budget_remaining
      how_many <- function(cost) {
        tiny <- left * 4 * .Machine$double.eps
        k <- floor(left / cost)
        if ((k + 1) * cost <= left + tiny) k <- k + 1
        k
      }
      if (left > 0) {
        extra <- character(0)
        if (x$cost_per_sample > 0) {
          k <- how_many(x$cost_per_sample)
          if (k > 0)
            extra <- c(extra, sprintf("%s more %s", int(k),
                                      if (k == 1) "sample" else "samples"))
        }
        if (!is.null(x$by_site)) {
          per_site <- x$total_fixed_cost / x$n_sites
          if (per_site > 0) {
            k <- how_many(per_site)
            if (k > 0)
              extra <- c(extra, sprintf(
                "%s more %s (at the average per-site cost)", int(k),
                if (k == 1) "site" else "sites"))
          }
        }
        if (length(extra) > 0) {
          cat(if (length(extra) > 1)
                "    leftover budget would cover roughly one of:\n"
              else "    leftover budget would cover roughly:\n")
          for (i in seq_along(extra))
            cat(sprintf("      %s%s\n", extra[i],
                        if (i < length(extra)) ", or" else ""))
        }
      }
    }
  }

  invisible(x)
}
