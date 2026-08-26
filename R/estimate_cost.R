#' Estimate total study cost
#'
#' Simple arithmetic cost model: per-sample cost multiplied by sample size,
#' plus optional fixed costs per site.  No course-module dependency.
#'
#' @param n Integer. Total number of samples.
#' @param cost_per_sample Numeric. Cost per sample (any currency unit).
#' @param n_sites Optional integer. Number of sites.
#' @param cost_per_site Numeric. Fixed cost per site (e.g. travel, setup);
#'   default 0.  Ignored when `n_sites` is `NULL`.
#'
#' @return An object of class `"mms_cost"` (a named list) with:
#'   \item{total_cost}{Total cost (sampling + site fixed costs)}
#'   \item{sampling_cost}{n * cost_per_sample}
#'   \item{site_cost}{n_sites * cost_per_site (0 if `n_sites` is `NULL`)}
#'   \item{n}{Samples}
#'   \item{n_sites}{Sites (may be `NULL`)}
#'   \item{cost_per_sample}{Cost per sample used}
#'   \item{cost_per_site}{Cost per site used}
#'
#' @export
#' @examples
#' # 323 samples at $15 each
#' estimate_cost(n = 323, cost_per_sample = 15)
#'
#' # With 10 sites and $200 fixed cost per site
#' estimate_cost(n = 323, cost_per_sample = 15, n_sites = 10, cost_per_site = 200)
estimate_cost <- function(
  n,
  cost_per_sample,
  n_sites       = NULL,
  cost_per_site = 0
) {
  if (!is.numeric(n) || n <= 0)          stop("`n` must be a positive number")
  if (!is.numeric(cost_per_sample) || cost_per_sample < 0)
    stop("`cost_per_sample` must be non-negative")

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
  fmt <- function(v) format(round(v), big.mark = ",")
  cat(sprintf("Total cost: %s\n", fmt(x$total_cost)))
  cat(sprintf("  Sampling (%s samples x %s): %s\n",
    fmt(x$n), fmt(x$cost_per_sample), fmt(x$sampling_cost)))
  if (!is.null(x$n_sites) && x$site_cost > 0)
    cat(sprintf("  Site fixed costs (%s sites x %s): %s\n",
      fmt(x$n_sites), fmt(x$cost_per_site), fmt(x$site_cost)))
  invisible(x)
}
