# Hand-checked reference cases for design_difference()
#
# Per-group two-proportion z-test sample size (equal allocation):
#   p_g_app = p_g * Se + (1 - p_g) * (1 - Sp)          [apparent scale]
#   pbar    = (p1_app + p2_app) / 2
#   delta   = |p1_app - p2_app|
#   se0     = sqrt(2 * pbar * (1 - pbar))
#   se1     = sqrt(p1_app*(1-p1_app) + p2_app*(1-p2_app))
#   z_a     = qnorm(1 - alpha)   one-sided;   qnorm(1 - alpha/2)   two-sided
#   z_b     = qnorm(power)
#   n       = ceiling( ((z_a*se0 + z_b*se1) / delta)^2 )
#
# Clustering / FPC: same closed-form block as design_threshold(), applied
# to each group.

nbase <- function(p1, p2, power = 0.80, alt = "two.sided",
                  se = 1, sp = 1, conf = 0.95) {
  p1a <- p1 * se + (1 - p1) * (1 - sp)
  p2a <- p2 * se + (1 - p2) * (1 - sp)
  pb  <- (p1a + p2a) / 2
  d   <- abs(p1a - p2a)
  za  <- if (alt == "two.sided") qnorm(1 - (1 - conf) / 2) else qnorm(conf)
  zb  <- qnorm(power)
  ((za * sqrt(2 * pb * (1 - pb)) +
    zb * sqrt(p1a * (1 - p1a) + p2a * (1 - p2a))) / d)^2
}


# ---------------------------------------------------------------------------
# Round 1 -- core forward-mode math (SRS)
# ---------------------------------------------------------------------------

test_that("DF-R1-1: SRS, perfect test, two-sided 80% power -- hand-checked", {
  res <- design_difference(prevalence1 = 0.10, prevalence2 = 0.20)

  expect_equal(res$n_per_group, ceiling(nbase(0.10, 0.20)))   # 199
  expect_equal(res$n_per_group, 199)
  expect_equal(res$n, 398)
  expect_equal(res$n_eff, 199)
  expect_equal(res$deff, 1)
  expect_equal(res$alternative, "two.sided")
  expect_equal(res$apparent_prev1, 0.10, tolerance = 1e-12)
  expect_equal(res$apparent_prev2, 0.20, tolerance = 1e-12)
  expect_equal(res$delta_app, -0.10, tolerance = 1e-12)
  expect_null(res$n_sites)
  expect_null(res$n_per_site)
})

test_that("DF-R1-2: higher power needs more samples", {
  r80 <- design_difference(0.10, 0.20, power = 0.80)
  r90 <- design_difference(0.10, 0.20, power = 0.90)
  expect_gt(r90$n_per_group, r80$n_per_group)
  expect_equal(r90$n_per_group, ceiling(nbase(0.10, 0.20, power = 0.90)))
})

test_that("DF-R1-3: one-sided test needs fewer samples than two-sided", {
  r_2s <- design_difference(0.20, 0.10, alternative = "two.sided")
  r_1s <- design_difference(0.20, 0.10, alternative = "greater")
  expect_lt(r_1s$n_per_group, r_2s$n_per_group)
  expect_equal(r_1s$n_per_group,
               ceiling(nbase(0.20, 0.10, alt = "greater")))
})

test_that("DF-R1-4: wider prevalence gap needs fewer samples", {
  r_narrow <- design_difference(0.10, 0.13)
  r_wide   <- design_difference(0.10, 0.30)
  expect_gt(r_narrow$n_per_group, r_wide$n_per_group)
})

test_that("DF-R1-5: larger conf_level (smaller alpha) needs more samples", {
  r95 <- design_difference(0.10, 0.20)
  r99 <- design_difference(0.10, 0.20, conf_level = 0.99)
  expect_gt(r99$n_per_group, r95$n_per_group)
})

test_that("DF-R1-6: result is symmetric in the two prevalences (two-sided)", {
  a <- design_difference(0.10, 0.20)
  b <- design_difference(0.20, 0.10)
  expect_equal(a$n_per_group, b$n_per_group)
})


# ---------------------------------------------------------------------------
# Round 2 -- imperfect diagnostics (Rogan-Gladen)
# ---------------------------------------------------------------------------

test_that("DF-R2-1: imperfect test shifts apparent prevalences and inflates n", {
  r_perfect   <- design_difference(0.10, 0.20)
  r_imperfect <- design_difference(0.10, 0.20, sensitivity = 0.9, specificity = 0.95)

  expect_equal(r_imperfect$apparent_prev1, 0.135, tolerance = 1e-12)
  expect_equal(r_imperfect$apparent_prev2, 0.22,  tolerance = 1e-12)
  # apparent gap = true gap * (Se + Sp - 1)
  expect_equal(r_imperfect$delta_app, -0.10 * 0.85, tolerance = 1e-12)
  expect_gt(r_imperfect$n_per_group, r_perfect$n_per_group)
  expect_equal(r_imperfect$n_per_group,
               ceiling(nbase(0.10, 0.20, se = 0.9, sp = 0.95)))
})

test_that("DF-R2-2: se + sp <= 1 is rejected", {
  expect_error(
    design_difference(0.10, 0.20, sensitivity = 0.5, specificity = 0.5),
    "must exceed 1"
  )
})

test_that("DF-R2-3: near-identical apparent prevalences explode -> stop", {
  # tiny true gap * an imperfect test collapses the apparent gap below 1e-8
  expect_error(
    design_difference(0.10, 0.10000005, sensitivity = 0.55, specificity = 0.55),
    "explode|identical apparent"
  )
})


# ---------------------------------------------------------------------------
# Round 3 -- clustering
# ---------------------------------------------------------------------------

test_that("DF-R3-1: fixed n_per_site inflates n by the Kish deff", {
  r_srs     <- design_difference(0.10, 0.20)
  r_cluster <- design_difference(0.10, 0.20, n_per_site = 30, icc = 0.02)

  expect_equal(r_cluster$deff, 1 + 29 * 0.02)          # 1.58
  expect_equal(r_cluster$n_per_group,
               ceiling(nbase(0.10, 0.20) * 1.58))
  expect_gt(r_cluster$n_per_group, r_srs$n_per_group)
  expect_equal(r_cluster$n_eff, r_srs$n_eff)           # SRS-equivalent unchanged
  expect_equal(r_cluster$n_sites, ceiling(r_cluster$n_per_group / 30))
  expect_equal(r_cluster$n_per_site, 30)
  expect_equal(r_cluster$n, 2 * r_cluster$n_per_group)
})

test_that("DF-R3-2: fixed n_sites solved via closed form, deff > 1", {
  res <- design_difference(0.10, 0.20, n_sites = 40, icc = 0.02)
  nb    <- nbase(0.10, 0.20)
  denom <- 40 - nb * 0.02
  n_cont <- nb * 40 * (1 - 0.02) / denom
  expect_equal(res$n_per_group, ceiling(n_cont))
  expect_gt(res$deff, 1)
  expect_equal(res$n_sites, 40)
  expect_equal(res$n_per_site, ceiling(res$n_per_group / 40))
})

test_that("DF-R3-3: icc = 0 with a cluster structure is plain SRS", {
  r_srs  <- design_difference(0.10, 0.20)
  r_zero <- design_difference(0.10, 0.20, n_per_site = 30, icc = 0)
  expect_equal(r_zero$deff, 1)
  expect_equal(r_zero$n_per_group, r_srs$n_per_group)
})

test_that("DF-R3-4: larger icc inflates n for fixed cluster size", {
  r_lo <- design_difference(0.10, 0.20, n_per_site = 30, icc = 0.01)
  r_hi <- design_difference(0.10, 0.20, n_per_site = 30, icc = 0.05)
  expect_gt(r_hi$n_per_group, r_lo$n_per_group)
  expect_gt(r_hi$deff, r_lo$deff)
})

test_that("DF-R3-5: n_sites >= per-group SRS size is rejected", {
  expect_error(
    design_difference(0.10, 0.20, n_sites = 5000, icc = 0.02),
    "not a valid cluster design"
  )
})

test_that("DF-R3-6: n_sites below n_base*icc is flagged unachievable", {
  expect_error(
    design_difference(0.10, 0.20, n_sites = 3, icc = 0.9),
    "unachievable"
  )
})


# ---------------------------------------------------------------------------
# Round 4 -- finite-population correction
# ---------------------------------------------------------------------------

test_that("DF-R4-1: FPC shrinks per-group n, leaves n_eff alone", {
  r_inf <- design_difference(0.10, 0.20)
  r_fpc <- design_difference(0.10, 0.20, fpc_N = 400)
  nb     <- nbase(0.10, 0.20)
  n_cont <- nb * 400 / (nb + 400 - 1)
  expect_equal(r_fpc$n_per_group, ceiling(n_cont))
  expect_lt(r_fpc$n_per_group, r_inf$n_per_group)
  expect_equal(r_fpc$n_eff, r_inf$n_eff)
  expect_equal(r_fpc$fpc_N, 400)
})


# ---------------------------------------------------------------------------
# Round 5 -- input validation
# ---------------------------------------------------------------------------

test_that("DF-R5-1: prevalence guards", {
  expect_error(design_difference(0, 0.2),        "strictly between 0 and 1")
  expect_error(design_difference(0.2, 1),        "strictly between 0 and 1")
  expect_error(design_difference("0.1", 0.2),    "must be a single number")
  expect_error(design_difference(c(0.1, 0.15), 0.2), "length 2")
})

test_that("DF-R5-2: power guards", {
  expect_error(design_difference(0.1, 0.2, power = 1.2),  "strictly between 0 and 1")
  expect_error(design_difference(0.1, 0.2, power = 0.3),  "at least 0.5")
  expect_error(design_difference(0.1, 0.2, power = TRUE), "must be a single number")
})

test_that("DF-R5-3: alternative guards and direction checks", {
  expect_error(design_difference(0.1, 0.2, alternative = "bigger"),
               "'greater', 'less', or 'two.sided'")
  expect_error(design_difference(0.1, 0.2, alternative = "greater"),
               "must be > `prevalence2`")
  expect_error(design_difference(0.2, 0.1, alternative = "less"),
               "must be < `prevalence2`")
})

test_that("DF-R5-4: sensitivity / specificity guards", {
  expect_error(design_difference(0.1, 0.2, sensitivity = TRUE),
               "must be a single number")
  expect_error(design_difference(0.1, 0.2, sensitivity = 0),
               "must be in \\(0, 1\\]")
  expect_error(design_difference(0.1, 0.2, specificity = 1.5),
               "must be in \\(0, 1\\]")
})

test_that("DF-R5-5: icc / cluster-structure guards", {
  expect_error(design_difference(0.1, 0.2, icc = -0.1),
               "must be in \\[0, 1\\]")
  expect_error(design_difference(0.1, 0.2, icc = 0.05),
               "requires a cluster structure")
  expect_error(
    design_difference(0.1, 0.2, n_sites = 10, n_per_site = 20, icc = 0.02),
    "at most one of"
  )
  expect_error(design_difference(0.1, 0.2, n_per_site = 2.5, icc = 0.02),
               "positive integer")
})

test_that("DF-R5-6: fpc_N guard (finite positive integer)", {
  expect_error(design_difference(0.1, 0.2, fpc_N = -5),      "positive integer")
  expect_error(design_difference(0.1, 0.2, fpc_N = c(1, 2)), "length = 2")
  expect_error(design_difference(0.1, 0.2, fpc_N = 402.7),   "positive integer")
  expect_silent(design_difference(0.1, 0.2, fpc_N = 400))
})


# ---------------------------------------------------------------------------
# Round 6 -- monotonicity sweeps
# ---------------------------------------------------------------------------

test_that("DF-R6-1: n_per_group decreases as the gap widens", {
  p2s <- c(0.12, 0.15, 0.20, 0.30, 0.50)
  ns  <- vapply(p2s, function(p2)
    design_difference(0.10, p2)$n_per_group, numeric(1))
  expect_true(all(diff(ns) < 0))
})

test_that("DF-R6-2: n_per_group rises monotonically with power", {
  pw <- c(0.6, 0.7, 0.8, 0.9, 0.95)
  ns <- vapply(pw, function(p)
    design_difference(0.10, 0.20, power = p)$n_per_group, numeric(1))
  expect_true(all(diff(ns) > 0))
})
