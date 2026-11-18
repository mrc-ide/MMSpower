# Internal helpers shared across MMSpower functions.
# None of these are exported.

# ---- input checks ----

# Stop unless `x` is a single finite number between `lower` and `upper`.
# `lower_open` / `upper_open` exclude that end of the range.
.check_number <- function(x, name, lower = -Inf, upper = Inf,
                          lower_open = FALSE, upper_open = FALSE) {
  ok <- is.numeric(x) && length(x) == 1 && is.finite(x) &&
    (if (lower_open) x > lower else x >= lower) &&
    (if (upper_open) x < upper else x <= upper)
  if (!ok)
    stop("`", name, "` must be a single number in ",
         if (lower_open) "(" else "[", lower, ", ", upper,
         if (upper_open) ")" else "]", ".", call. = FALSE)
}

# Stop unless `x` is NULL or a single whole number of at least `min`.
.check_count <- function(x, name, min = 1) {
  if (is.null(x)) return(invisible())
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x != round(x) || x < min)
    stop("`", name, "` must be a single whole number of at least ", min, ".",
         call. = FALSE)
}

# Checks shared by every function that takes a diagnostic test.
.check_test <- function(sensitivity, specificity) {
  .check_number(sensitivity, "sensitivity", 0, 1, lower_open = TRUE)
  .check_number(specificity, "specificity", 0, 1, lower_open = TRUE)
  if (sensitivity + specificity <= 1)
    stop("`sensitivity` + `specificity` must be greater than 1.", call. = FALSE)
}

# Checks shared by the design functions: at most one of n_sites and
# n_per_site, and an icc above 0 only with one of them.
.check_sites <- function(n_sites, n_per_site, icc) {
  .check_count(n_sites, "n_sites")
  .check_count(n_per_site, "n_per_site")
  .check_number(icc, "icc", 0, 1)
  if (!is.null(n_sites) && !is.null(n_per_site))
    stop("Supply at most one of `n_sites` and `n_per_site`.", call. = FALSE)
  if (icc > 0 && is.null(n_sites) && is.null(n_per_site))
    stop("`icc` > 0 needs `n_sites` or `n_per_site`.", call. = FALSE)
}

# ---- diagnostic test ----

# True prevalence -> apparent prevalence (forward Rogan-Gladen).
.apparent_prev <- function(prevalence, sensitivity, specificity) {
  prevalence * sensitivity + (1 - prevalence) * (1 - specificity)
}

# Apparent prevalence -> true prevalence (inverse Rogan-Gladen).
.rogan_gladen <- function(p_apparent, sensitivity, specificity) {
  (p_apparent - (1 - specificity)) / (sensitivity + specificity - 1)
}

# ---- clustering ----

# Design effect (Kish) for sites of `n_per_site` people.
.deff <- function(n_per_site, icc) {
  1 + (n_per_site - 1) * icc
}

# Total sample size with `n_sites` fixed, given the unclustered size `n0`.
# Solves n = n0 * deff(n / n_sites, icc) for n. NA when no number of people
# per site is enough (n_sites <= n0 * icc).
.solve_fixed_sites <- function(n0, n_sites, icc) {
  denom <- n_sites - n0 * icc
  if (denom <= 0) return(NA_real_)
  n0 * n_sites * (1 - icc) / denom
}

# Sites and people per site from a total sample size `n`, rounded up to
# whole sites and whole people. Returns list(n_sites, n_per_site, n), where
# n = n_sites * n_per_site is the number actually sampled.
.allocate <- function(n, n_sites = NULL, n_per_site = NULL) {
  if (!is.null(n_per_site)) n_sites <- ceiling(n / n_per_site)
  else if (!is.null(n_sites)) n_per_site <- ceiling(n / n_sites)
  else return(list(n_sites = NULL, n_per_site = NULL, n = ceiling(n)))
  list(n_sites = n_sites, n_per_site = n_per_site, n = n_sites * n_per_site)
}

# ---- intervals ----

# Wilson score confidence interval for a proportion.
# Returns list(lower, upper).
.wilson_ci <- function(p, n, z) {
  denom  <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  margin <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  list(lower = center - margin, upper = center + margin)
}
