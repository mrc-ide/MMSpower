# Hand-checked reference cases for design_detection()
#
# Detection formula (perfect / SRS):
#   q      = prevalence * sensitivity        [true-positive prob per sample]
#   n_base = log(1 - detection_prob) / log(1 - q)
#   n      = ceiling(n_base)
#
# Clustering / FPC use the same closed-form logic as design_threshold().
#
# NOTE: specificity is deliberately unused in v1 (see the "Known limitation"
# section in the roxygen and the README tracking note). Tests below assert
# that specificity has no effect; update them when that decision is made.


# ---------------------------------------------------------------------------
# Round 1 -- core forward-mode math (SRS)
# ---------------------------------------------------------------------------

test_that("DD-R1-1: SRS, perfect test -- hand-checked n", {
  # p = 0.02, detection_prob = 0.95, Se = 1
  # q = 0.02
  # n_base = log(0.05) / log(0.98) = -2.995732 / -0.0202027 = 148.283
  res <- design_detection(prevalence = 0.02)

  n_expected <- ceiling(log(1 - 0.95) / log(1 - 0.02))
  expect_equal(n_expected, 149)
  expect_equal(res$n, 149)
  expect_equal(res$n_eff, 149)
  expect_equal(res$q, 0.02)
  expect_equal(res$deff, 1)
  expect_equal(res$mode, "solve_n")
  expect_null(res$n_sites)
  expect_null(res$n_per_site)
  expect_null(res$fpc_N)
})

test_that("DD-R1-2: rarer variant needs more samples", {
  r_common <- design_detection(prevalence = 0.05)
  r_rare   <- design_detection(prevalence = 0.005)
  expect_gt(r_rare$n, r_common$n)

  expect_equal(r_rare$n, ceiling(log(0.05) / log(1 - 0.005)))
})

test_that("DD-R1-3: higher detection_prob needs more samples", {
  r95 <- design_detection(prevalence = 0.02, detection_prob = 0.95)
  r99 <- design_detection(prevalence = 0.02, detection_prob = 0.99)
  expect_gt(r99$n, r95$n)
  expect_equal(r99$n, ceiling(log(1 - 0.99) / log(0.98)))  # 228
  expect_equal(r99$n, 228)
})

test_that("DD-R1-4: imperfect sensitivity inflates the requirement", {
  r_perfect   <- design_detection(prevalence = 0.02, sensitivity = 1)
  r_imperfect <- design_detection(prevalence = 0.02, sensitivity = 0.8)
  expect_gt(r_imperfect$n, r_perfect$n)
  expect_equal(r_imperfect$q, 0.016)
  expect_equal(r_imperfect$n, ceiling(log(0.05) / log(1 - 0.016)))  # 186
  expect_equal(r_imperfect$n, 186)
})

test_that("DD-R1-5: specificity is currently ignored (tracks known limitation)", {
  r_hi <- design_detection(prevalence = 0.02, specificity = 1)
  r_lo <- design_detection(prevalence = 0.02, specificity = 0.90)
  expect_equal(r_hi$n, r_lo$n)
  expect_equal(r_hi$q, r_lo$q)
  # specificity is still echoed back unchanged
  expect_equal(r_lo$specificity, 0.90)
})


# ---------------------------------------------------------------------------
# Round 2 -- reverse mode (n supplied -> detection probability)
# ---------------------------------------------------------------------------

test_that("DD-R2-1: reverse mode reports achieved detection probability", {
  # p = 0.02, n = 150, Se = 1  ->  1 - 0.98^150
  res <- design_detection(prevalence = 0.02, n = 150)
  expect_equal(res$mode, "solve_detection_prob")
  expect_equal(res$n, 150L)
  expect_equal(res$detection_prob, 1 - 0.98^150, tolerance = 1e-12)
  expect_gt(res$detection_prob, 0.95)   # n=150 exceeds the n=149 needed for 0.95
})

test_that("DD-R2-2: forward / reverse round-trip at the boundary", {
  fwd <- design_detection(prevalence = 0.02, detection_prob = 0.95)   # n = 149
  rev <- design_detection(prevalence = 0.02, n = fwd$n)
  expect_gte(rev$detection_prob, 0.95)

  # one fewer sample drops below the target
  rev_minus <- design_detection(prevalence = 0.02, n = fwd$n - 1)
  expect_lt(rev_minus$detection_prob, 0.95)
})

test_that("DD-R2-3: reverse mode, detection_prob argument is ignored", {
  a <- design_detection(prevalence = 0.02, n = 150)
  b <- design_detection(prevalence = 0.02, n = 150, detection_prob = 0.5)
  expect_equal(a$detection_prob, b$detection_prob)
})

test_that("DD-R2-4: reverse mode with fixed cluster size deflates effective n", {
  # deff = 1 + (20 - 1) * 0.05 = 1.95
  res <- design_detection(prevalence = 0.02, n = 300,
                          n_per_site = 20, icc = 0.05)
  expect_equal(res$deff, 1.95)
  expect_equal(res$detection_prob, 1 - 0.98^(300 / 1.95), tolerance = 1e-10)
  expect_equal(res$n_sites, ceiling(300 / 20))
  expect_equal(res$n_per_site, 20)
})


# ---------------------------------------------------------------------------
# Round 3 -- clustering (forward mode)
# ---------------------------------------------------------------------------

test_that("DD-R3-1: fixed n_per_site inflates n by the Kish deff", {
  r_srs     <- design_detection(prevalence = 0.02)
  r_cluster <- design_detection(prevalence = 0.02, n_per_site = 20, icc = 0.05)
  n_base <- log(0.05) / log(0.98)                      # continuous base, ~148.28
  expect_equal(r_cluster$deff, 1 + 19 * 0.05)          # 1.95
  expect_equal(r_cluster$n, ceiling(n_base * 1.95))    # deff applied to the
                                                       # continuous base, not n_eff
  expect_gt(r_cluster$n, r_srs$n)
  expect_equal(r_cluster$n_eff, r_srs$n_eff)           # SRS-equivalent unchanged
  expect_equal(r_cluster$n_sites, ceiling(r_cluster$n / 20))
  expect_equal(r_cluster$n_per_site, 20)
})

test_that("DD-R3-2: fixed n_sites solved via closed form, deff > 1", {
  res <- design_detection(prevalence = 0.02, n_sites = 50, icc = 0.05)
  n_base <- log(0.05) / log(0.98)
  denom  <- 50 - n_base * 0.05
  n_cont <- n_base * 50 * (1 - 0.05) / denom
  expect_equal(res$n, ceiling(n_cont))
  expect_gt(res$deff, 1)
  expect_equal(res$n_sites, 50)
  expect_equal(res$n_per_site, ceiling(res$n / 50))
})

test_that("DD-R3-3: icc = 0 with a cluster structure is plain SRS", {
  r_srs  <- design_detection(prevalence = 0.02)
  r_zero <- design_detection(prevalence = 0.02, n_per_site = 20, icc = 0)
  expect_equal(r_zero$deff, 1)
  expect_equal(r_zero$n, r_srs$n)
})

test_that("DD-R3-4: n_sites >= SRS size is rejected (each site < 1 person)", {
  expect_error(
    design_detection(prevalence = 0.02, n_sites = 500, icc = 0.05),
    "not a valid cluster design"
  )
})

test_that("DD-R3-5: n_sites below n_base*icc is flagged as unachievable", {
  expect_error(
    design_detection(prevalence = 0.02, n_sites = 100, icc = 0.9),
    "unachievable"
  )
})


# ---------------------------------------------------------------------------
# Round 4 -- finite-population correction
# ---------------------------------------------------------------------------

test_that("DD-R4-1: FPC shrinks the required n, leaves n_eff alone", {
  r_inf <- design_detection(prevalence = 0.02)
  r_fpc <- design_detection(prevalence = 0.02, fpc_N = 500)
  n_base <- log(0.05) / log(0.98)
  n_cont <- n_base * 500 / (n_base + 500 - 1)
  expect_equal(r_fpc$n, ceiling(n_cont))
  expect_lt(r_fpc$n, r_inf$n)
  expect_equal(r_fpc$n_eff, r_inf$n_eff)
  expect_equal(r_fpc$fpc_N, 500)
})

test_that("DD-R4-2: reverse mode rejects n >= fpc_N", {
  expect_error(
    design_detection(prevalence = 0.02, n = 600, fpc_N = 500),
    "smaller than"
  )
})


# ---------------------------------------------------------------------------
# Round 5 -- input validation (what/why/fix errors)
# ---------------------------------------------------------------------------

test_that("DD-R5-1: prevalence must be a single number in (0, 1)", {
  expect_error(design_detection(prevalence = 0),      "strictly between 0 and 1")
  expect_error(design_detection(prevalence = 1),      "strictly between 0 and 1")
  expect_error(design_detection(prevalence = -0.1),   "strictly between 0 and 1")
  expect_error(design_detection(prevalence = "0.02"), "must be a single number")
  expect_error(design_detection(prevalence = c(0.02, 0.03)), "length 2")
})

test_that("DD-R5-2: detection_prob must be a single number in (0, 1)", {
  expect_error(design_detection(prevalence = 0.02, detection_prob = 1.5),
               "strictly between 0 and 1")
  expect_error(design_detection(prevalence = 0.02, detection_prob = 0),
               "strictly between 0 and 1")
  expect_error(design_detection(prevalence = 0.02, detection_prob = TRUE),
               "must be a single number")
})

test_that("DD-R5-3: n must be a finite positive integer when supplied", {
  expect_error(design_detection(prevalence = 0.02, n = 2.5),  "positive integer")
  expect_error(design_detection(prevalence = 0.02, n = 0),    "positive integer")
  expect_error(design_detection(prevalence = 0.02, n = -5),   "positive integer")
  expect_error(design_detection(prevalence = 0.02, n = Inf),  "positive integer")
})

test_that("DD-R5-4: sensitivity / specificity guards", {
  expect_error(design_detection(prevalence = 0.02, sensitivity = TRUE),
               "must be a single number")
  expect_error(design_detection(prevalence = 0.02, sensitivity = 0),
               "must be in \\(0, 1\\]")
  expect_error(design_detection(prevalence = 0.02, sensitivity = 1.2),
               "must be in \\(0, 1\\]")
  expect_error(design_detection(prevalence = 0.02, specificity = 0),
               "must be in \\(0, 1\\]")
  expect_error(design_detection(prevalence = 0.02, specificity = 1.5),
               "must be in \\(0, 1\\]")
})

test_that("DD-R5-5: icc guards and cluster-structure requirement", {
  expect_error(design_detection(prevalence = 0.02, icc = "high"),
               "must be a single number")
  expect_error(design_detection(prevalence = 0.02, icc = -0.1),
               "must be in \\[0, 1\\]")
  expect_error(design_detection(prevalence = 0.02, icc = 0.05),
               "requires a cluster structure")
})

test_that("DD-R5-6: n_sites / n_per_site guards", {
  expect_error(
    design_detection(prevalence = 0.02, n_sites = 10, n_per_site = 20, icc = 0.05),
    "at most one of"
  )
  expect_error(design_detection(prevalence = 0.02, n_sites = 3.5, icc = 0.05),
               "positive integer")
  expect_error(design_detection(prevalence = 0.02, n_per_site = -4, icc = 0.05),
               "positive integer")
})

test_that("DD-R5-7: fpc_N guard (finite positive integer)", {
  expect_error(design_detection(prevalence = 0.02, fpc_N = -10),
               "positive integer")
  expect_error(design_detection(prevalence = 0.02, fpc_N = 500.5),
               "positive integer")
  expect_error(design_detection(prevalence = 0.02, fpc_N = c(100, 200)),
               "length = 2")
})

test_that("DD-R5-8: a cluster structure larger than fpc_N is rejected", {
  expect_error(
    design_detection(prevalence = 0.02, n_sites = 100, fpc_N = 50, icc = 0.05),
    "more clusters than individuals"
  )
  expect_error(
    design_detection(prevalence = 0.02, n_per_site = 100, fpc_N = 50, icc = 0.05),
    "cannot be larger than the whole population"
  )
})


# ---------------------------------------------------------------------------
# Round 6 -- monotonicity / sanity sweeps
# ---------------------------------------------------------------------------

test_that("DD-R6-1: n decreases monotonically as prevalence rises", {
  ps <- c(0.005, 0.01, 0.02, 0.05, 0.10)
  ns <- vapply(ps, function(p) design_detection(prevalence = p)$n, numeric(1))
  expect_true(all(diff(ns) < 0))
})

test_that("DD-R6-2: achieved detection probability rises monotonically with n", {
  ns <- c(50, 100, 150, 200, 300)
  ds <- vapply(ns, function(k)
    design_detection(prevalence = 0.02, n = k)$detection_prob, numeric(1))
  expect_true(all(diff(ds) > 0))
})

test_that("DD-R6-3: larger icc inflates n for a fixed cluster size", {
  r_lo <- design_detection(prevalence = 0.02, n_per_site = 20, icc = 0.02)
  r_hi <- design_detection(prevalence = 0.02, n_per_site = 20, icc = 0.10)
  expect_gt(r_hi$n, r_lo$n)
  expect_gt(r_hi$deff, r_lo$deff)
})


# ---------------------------------------------------------------------------
# Round 7 -- code-review fixes (2026-09-02): reverse-mode consistency,
# finite-population feasibility, return-type consistency
# ---------------------------------------------------------------------------

test_that("DD-R7-1: reverse mode rejects n below the cluster layout", {
  expect_error(
    design_detection(prevalence = 0.02, n = 50, n_sites = 100, icc = 0.5),
    "smaller than `n_sites`"
  )
  expect_error(
    design_detection(prevalence = 0.02, n = 10, n_per_site = 20, icc = 0.05),
    "smaller than `n_per_site`"
  )
})

test_that("DD-R7-2: reverse-mode deff never drops below 1", {
  # n_sites just under n: average cluster size ~1, deff must still be >= 1
  res <- design_detection(prevalence = 0.02, n = 120, n_sites = 100, icc = 0.5)
  expect_gte(res$deff, 1)
  expect_lte(res$n_eff, res$n)            # effective n cannot exceed collected n
})

test_that("DD-R7-3: forward mode errors when fpc_N cannot reach the target", {
  # census of 50 people at prevalence 0.02 gives at most 1 - 0.98^50 ~ 0.636
  expect_error(
    design_detection(prevalence = 0.02, fpc_N = 50),
    "not achievable in a population"
  )
  # a population large enough to reach 95% is fine
  expect_silent(design_detection(prevalence = 0.02, fpc_N = 5000))
})

test_that("DD-R7-4: `n` has a consistent (numeric) type in both modes", {
  fwd <- design_detection(prevalence = 0.02)
  rev <- design_detection(prevalence = 0.02, n = 150)
  expect_type(fwd$n, "double")
  expect_type(rev$n, "double")
})
