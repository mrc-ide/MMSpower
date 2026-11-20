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

# Stop unless `x` (positives) and `n` (people tested) are equal-length
# vectors of whole numbers, one per site, with 0 <= x <= n and n >= 1.
.check_counts <- function(x, n, x_name = "x", n_name = "n") {
  whole <- function(v) is.numeric(v) && length(v) > 0 && all(is.finite(v)) &&
    all(v == round(v))
  if (!whole(x) || !whole(n))
    stop("`", x_name, "` and `", n_name, "` must be vectors of whole numbers.",
         call. = FALSE)
  if (length(x) != length(n))
    stop("`", x_name, "` and `", n_name, "` must have the same length.",
         call. = FALSE)
  if (any(n < 1))
    stop("Every value of `", n_name, "` must be at least 1.", call. = FALSE)
  if (any(x < 0 | x > n))
    stop("Every value of `", x_name, "` must be between 0 and `", n_name, "`.",
         call. = FALSE)
}

# Checks shared by every function that takes a diagnostic test.
.check_test <- function(sensitivity, specificity) {
  .check_number(sensitivity, "sensitivity", 0, 1, lower_open = TRUE)
  .check_number(specificity, "specificity", 0, 1, lower_open = TRUE)
  if (sensitivity + specificity <= 1)
    stop("`sensitivity` + `specificity` must be greater than 1.", call. = FALSE)
}

# Checks and reads the design given to a design function (`suffix` "1" or
# "2" names a group's arguments in design_difference()). A design is either
# solved for, or given in full:
#   - nothing: solve for the number of people (no sites)
#   - n: a full design with no sites
#   - n_per_site alone: solve for the number of sites
#   - n_sites alone: solve for the number of people per site
#   - n_sites and n_per_site, or a vector n_per_site of one size per site:
#     a full design with sites
# Returns list(full, n, n_w, sizes): for a full design the number of people,
# the size-weighted mean site size sum(s^2) / sum(s) (1 with no sites) and
# the site sizes (NULL with no sites); NULL entries otherwise.
.site_design <- function(n, n_sites, n_per_site, icc, suffix = "") {
  nm <- function(x) paste0("`", x, suffix, "`")
  .check_count(n, paste0("n", suffix))
  .check_count(n_sites, paste0("n_sites", suffix))
  if (!is.null(n_per_site) &&
      (!is.numeric(n_per_site) || length(n_per_site) < 1 ||
       !all(is.finite(n_per_site)) || any(n_per_site != round(n_per_site)) ||
       any(n_per_site < 1)))
    stop(nm("n_per_site"), " must be whole numbers of at least 1.", call. = FALSE)
  if (!is.null(n) && (!is.null(n_sites) || !is.null(n_per_site)))
    stop(nm("n"), " is for a design without sites; with sites, give ",
         nm("n_sites"), " and ", nm("n_per_site"), ".", call. = FALSE)
  if (length(n_per_site) > 1 && !is.null(n_sites))
    stop("Give ", nm("n_sites"), " with one ", nm("n_per_site"), ", or a vector ",
         nm("n_per_site"), " with one size per site, not both.", call. = FALSE)
  if (icc > 0 && is.null(n_sites) && is.null(n_per_site))
    stop("`icc` > 0 needs ", nm("n_sites"), " or ", nm("n_per_site"), ".",
         call. = FALSE)

  sizes <- if (length(n_per_site) > 1) n_per_site
           else if (!is.null(n_sites) && !is.null(n_per_site)) rep(n_per_site, n_sites)
  if (!is.null(n)) return(list(full = TRUE, n = n, n_w = 1, sizes = NULL))
  if (!is.null(sizes))
    return(list(full = TRUE, n = sum(sizes), n_w = sum(sizes^2) / sum(sizes),
                sizes = sizes))
  list(full = FALSE, n = NULL, n_w = NULL, sizes = NULL)
}

# Stop if both the target and a full design are given, or (for a target
# with no default) neither. `given` says whether the target was typed in;
# a default does not count.
.check_target <- function(full, given, target, has_default = TRUE) {
  if (full && given)
    stop("Give either `", target, "` or a full design, not both: with a full ",
         "design the function returns `", target, "`.", call. = FALSE)
  if (!full && !given && !has_default)
    stop("Give `", target, "`, or a full design (`n`, or `n_sites` and ",
         "`n_per_site`) to get the `", target, "` it gives.", call. = FALSE)
}

# Sites and people per site to report for a full design: one size if all
# sites are the same, otherwise the vector of sizes.
.report_sizes <- function(sizes) {
  if (is.null(sizes)) return(list(n_sites = NULL, n_per_site = NULL))
  list(n_sites    = as.numeric(length(sizes)),
       n_per_site = if (all(sizes == sizes[1])) sizes[1] else sizes)
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

# ICC and design effect for observed counts `x` out of `n`, one per site.
# A supplied `icc` is used as given. With `icc = NULL` it is estimated from
# the ratio R of the between-site variance of the site prevalences to the
# variance expected if everyone were independent. Under the usual model
# R = 1 + (n_h - 1) * icc, where n_h is the harmonic mean site size, so that
# is the mean used to recover the ICC. With one site, or one person per
# site, there is nothing to estimate from and it is 0. The design effect of
# the pooled prevalence uses the size-weighted mean site size
# sum(n^2) / sum(n) (both means are the plain mean when all sites are the
# same size). `df` is the degrees of freedom of the estimated variance,
# sites - 1, or Inf when the ICC was supplied or not estimated.
# Returns list(icc, deff, df).
.estimate_deff <- function(x, n, icc = NULL) {
  n_w <- sum(n^2) / sum(n)
  df  <- Inf
  if (is.null(icc)) {
    icc <- 0
    n_h <- 1 / mean(1 / n)
    if (length(n) > 1 && n_h > 1) {
      df      <- length(n) - 1
      p_hat   <- sum(x) / sum(n)
      var_obs <- stats::var(x / n)
      var_srs <- mean(p_hat * (1 - p_hat) / n)
      if (var_srs > 0)
        icc <- min(max((var_obs / var_srs - 1) / (n_h - 1), 0), 1)
    }
  }
  list(icc = icc, deff = .deff(n_w, icc), df = df)
}

# Effective sample size for an interval or test with tail area `a`, when
# the design effect was estimated with `df` degrees of freedom: n_eff
# shrunk by (z / t)^2 (Korn & Graubard, 1998). For a Wald interval this is
# the same as using a t value in place of z. Unchanged when df is Inf.
.df_adjust <- function(n_eff, a, df) {
  n_eff * (stats::qnorm(a, lower.tail = FALSE) / stats::qt(a, df, lower.tail = FALSE))^2
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

# ---- exact binomial tails ----
# Counts may be fractional (people divided by the design effect). For whole
# numbers these are exactly pbinom(); for fractions they are the usual
# smooth extension through the beta distribution (as in Clopper-Pearson).

# P(X >= x) for X ~ Binomial(m, p)
.binom_upper <- function(x, m, p) {
  out <- stats::pbeta(p, pmax(x, 1e-300), pmax(m - x + 1, 1e-300))
  out[x <= 0] <- 1
  out
}

# P(X <= x) for X ~ Binomial(m, p)
.binom_lower <- function(x, m, p) {
  out <- stats::pbeta(p, pmax(x + 1, 1e-300), pmax(m - x, 1e-300),
                      lower.tail = FALSE)
  out[x >= m] <- 1
  out
}

# ---- intervals ----

# Interval c(lower, upper) on the apparent-prevalence scale for measured
# prevalence p_hat and effective sample size n (not always whole), by
# `method`, kept within [0, 1]. Shared by estimate_prevalence() and
# design_precision().
.ci_app <- function(p_hat, n, alpha, method) {
  z <- stats::qnorm(1 - alpha / 2)
  x <- p_hat * n
  ci <- switch(method,
    wald = p_hat + c(-1, 1) * z * sqrt(p_hat * (1 - p_hat) / n),
    "clopper-pearson" = c(
      if (p_hat == 0) 0 else stats::qbeta(alpha / 2, x, n - x + 1),
      if (p_hat == 1) 1 else stats::qbeta(1 - alpha / 2, x + 1, n - x)),
    "agresti-coull" = {
      n_t <- n + z^2
      p_t <- (x + z^2 / 2) / n_t
      p_t + c(-1, 1) * z * sqrt(p_t * (1 - p_t) / n_t)
    })
  pmin(pmax(ci, 0), 1)
}

# Interval c(lower, upper) on the apparent-prevalence scale -> true scale
# (Rogan-Gladen), kept within [0, 1]. Warns when the whole corrected
# interval lies outside [0, 1]: the data then do not fit the assumed test
# accuracy (e.g. fewer positives than false positives alone would give).
.correct_ci <- function(ci_app, sensitivity, specificity) {
  ci <- .rogan_gladen(ci_app, sensitivity, specificity)
  if (ci[2] < 0 || ci[1] > 1)
    warning("The corrected interval lies entirely outside [0, 1], so the data ",
            "do not fit the assumed sensitivity and specificity.", call. = FALSE)
  pmin(pmax(ci, 0), 1)
}
