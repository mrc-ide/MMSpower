# Internal helpers shared across MMSpower functions.
# None of these are exported.

# Wilson score confidence interval for a proportion.
# Returns list(lower, upper).
.wilson_ci <- function(p, n, z) {
  denom  <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  margin <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  list(lower = center - margin, upper = center + margin)
}

# True prevalence -> apparent prevalence (forward Rogan-Gladen).
.apparent_prev <- function(prevalence, sensitivity, specificity) {
  prevalence * sensitivity + (1 - prevalence) * (1 - specificity)
}

# Apparent prevalence -> true prevalence (inverse Rogan-Gladen).
.rogan_gladen <- function(p_apparent, sensitivity, specificity) {
  (p_apparent - (1 - specificity)) / (sensitivity + specificity - 1)
}

# Compute design effect (DEFF) from ICC and cluster size.
#
# Returns a list:
#   deff       -- Kish design effect (>= 1); 1 means no clustering
#   n_per_site -- cluster size used (may be NULL)
#   n_sites    -- number of sites supplied (may be NULL)
#
# Planned but not yet implemented:
#   allocation = "pps"  -- proportional-to-size allocation
#   weighting  = "site" -- site-level weighting
.sampling_design <- function(
  n_per_site = NULL,
  n_sites    = NULL,
  icc        = 0,
  allocation = c("equal", "pps"),
  weighting  = c("none", "site")
) {
  allocation <- match.arg(allocation)
  weighting  <- match.arg(weighting)

  if (allocation == "pps")  stop("PPS allocation is not yet implemented")
  if (weighting  == "site") stop("Site weighting is not yet implemented")

  deff <- 1
  if (!is.null(n_per_site) && n_per_site > 1 && icc > 0) {
    deff <- 1 + (n_per_site - 1) * icc   # Kish (1965) formula
  }

  list(deff = deff, n_per_site = n_per_site, n_sites = n_sites)
}
