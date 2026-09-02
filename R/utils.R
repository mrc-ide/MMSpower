# Internal helpers shared across MMSpower functions.
# None of these are exported.

# Wilson score confidence interval for a proportion.
# Returns list(lower, upper). Used by test_difference() (Newcombe interval).
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
