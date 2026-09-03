# Hand-checked reference cases for design_trend()
#
# Slope z-test, per-timepoint sample size:
#   p_g_app    = p_g * Se + (1 - p_g) * (1 - Sp)
#   slope_app  = (p_end_app - p_start_app) / t_span      [= slope_true * (Se+Sp-1)]
#   pbar_app   = (p_start_app + p_end_app) / 2
#   sigma2_app = pbar_app * (1 - pbar_app)
#   Sxx        = sum((t - mean(t))^2)
#   z_a        = qnorm(1 - alpha)  one-sided ; qnorm(1 - alpha/2)  two-sided
#   z_b        = qnorm(power)
#   n          = ceiling( (z_a + z_b)^2 * sigma2_app / (Sxx * slope_app^2) )
#
# Clustering / FPC: same closed-form block as design_threshold().

nbase_trend <- function(p0, p1, tvec, power = 0.80, alt = "two.sided",
                        se = 1, sp = 1, conf = 0.95) {
  tspan <- diff(range(tvec))
  sxx   <- sum((tvec - mean(tvec))^2)
  p0a   <- p0 * se + (1 - p0) * (1 - sp)
  p1a   <- p1 * se + (1 - p1) * (1 - sp)
  slope_app  <- (p1a - p0a) / tspan
  pbar_app   <- (p0a + p1a) / 2
  sigma2_app <- pbar_app * (1 - pbar_app)
  za <- if (alt == "two.sided") qnorm(1 - (1 - conf) / 2) else qnorm(conf)
  zb <- qnorm(power)
  (za + zb)^2 * sigma2_app / (sxx * slope_app^2)
}


# ---------------------------------------------------------------------------
# Round 1 -- core forward-mode math
# ---------------------------------------------------------------------------

test_that("DT-R1-1: 10% -> 20% over 5 rounds, perfect test -- hand-checked", {
  res <- design_trend(prevalence_start = 0.10, prevalence_end = 0.20,
                      n_timepoints = 5)

  expect_equal(res$times, 0:4)
  expect_equal(res$sxx, 10)
  expect_equal(res$slope, 0.025, tolerance = 1e-12)
  expect_equal(res$slope_app, 0.025, tolerance = 1e-12)
  expect_equal(res$sigma2_app, 0.15 * 0.85, tolerance = 1e-12)
  expect_equal(res$n_per_timepoint, ceiling(nbase_trend(0.10, 0.20, 0:4)))  # 161
  expect_equal(res$n_per_timepoint, 161)
  expect_equal(res$n, 161 * 5)
  expect_equal(res$n_eff, 161)
  expect_equal(res$deff, 1)
  expect_equal(res$mode, "solve_n")
})

test_that("DT-R1-2: slope input matches the equivalent prevalence_end input", {
  a <- design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5)
  b <- design_trend(0.10, slope = 0.025, n_timepoints = 5)
  expect_equal(a$n_per_timepoint, b$n_per_timepoint)
  expect_equal(b$prevalence_end, 0.20, tolerance = 1e-12)
})

test_that("DT-R1-3: higher power needs more samples", {
  r80 <- design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5, power = 0.80)
  r90 <- design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5, power = 0.90)
  expect_gt(r90$n_per_timepoint, r80$n_per_timepoint)
  expect_equal(r90$n_per_timepoint,
               ceiling(nbase_trend(0.10, 0.20, 0:4, power = 0.90)))
})

test_that("DT-R1-4: a steeper trend needs fewer samples", {
  r_shallow <- design_trend(0.10, prevalence_end = 0.13, n_timepoints = 5)
  r_steep   <- design_trend(0.10, prevalence_end = 0.30, n_timepoints = 5)
  expect_gt(r_shallow$n_per_timepoint, r_steep$n_per_timepoint)
})

test_that("DT-R1-5: more timepoints (larger Sxx) needs fewer samples per round", {
  r5 <- design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5)
  r9 <- design_trend(0.10, prevalence_end = 0.20, n_timepoints = 9)
  expect_gt(r9$sxx, r5$sxx)
  expect_lt(r9$n_per_timepoint, r5$n_per_timepoint)
})

test_that("DT-R1-6: one-sided test needs fewer samples than two-sided", {
  r_2s <- design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5,
                       alternative = "two.sided")
  r_1s <- design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5,
                       alternative = "greater")
  expect_lt(r_1s$n_per_timepoint, r_2s$n_per_timepoint)
})

test_that("DT-R1-7: custom unequally-spaced times are honoured", {
  res <- design_trend(0.10, prevalence_end = 0.20, times = c(0, 1, 4))
  expect_equal(res$times, c(0, 1, 4))
  expect_equal(res$sxx, sum((c(0, 1, 4) - mean(c(0, 1, 4)))^2))
  expect_equal(res$slope, (0.20 - 0.10) / 4, tolerance = 1e-12)
})


# ---------------------------------------------------------------------------
# Round 2 -- imperfect diagnostics
# ---------------------------------------------------------------------------

test_that("DT-R2-1: imperfect test shrinks the apparent slope and inflates n", {
  r_perfect   <- design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5)
  r_imperfect <- design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5,
                              sensitivity = 0.9, specificity = 0.95)
  expect_equal(r_imperfect$slope_app, 0.025 * 0.85, tolerance = 1e-12)
  expect_gt(r_imperfect$n_per_timepoint, r_perfect$n_per_timepoint)
  expect_equal(r_imperfect$n_per_timepoint,
               ceiling(nbase_trend(0.10, 0.20, 0:4, se = 0.9, sp = 0.95)))
})

test_that("DT-R2-2: se + sp <= 1 is rejected", {
  expect_error(
    design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5,
                 sensitivity = 0.4, specificity = 0.6),
    "must exceed 1"
  )
})


# ---------------------------------------------------------------------------
# Round 3 -- reverse mode (n supplied -> power)
# ---------------------------------------------------------------------------

test_that("DT-R3-1: reverse mode reports achieved power", {
  fwd <- design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5, power = 0.80)
  rev <- design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5,
                      n = fwd$n_per_timepoint)
  expect_equal(rev$mode, "solve_power")
  expect_gt(rev$power, 0.80)                     # ceiling(n) exceeds the target
  expect_lt(rev$power, 0.82)

  rev_lo <- design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5,
                         n = fwd$n_per_timepoint - 1)
  expect_lt(rev_lo$power, 0.80)
})

test_that("DT-R3-2: power rises monotonically with n", {
  ns <- c(50, 100, 161, 250, 400)
  pw <- vapply(ns, function(k)
    design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5, n = k)$power,
    numeric(1))
  expect_true(all(diff(pw) > 0))
})


# ---------------------------------------------------------------------------
# Round 4 -- clustering
# ---------------------------------------------------------------------------

test_that("DT-R4-1: fixed n_per_site inflates n by the Kish deff", {
  r_srs     <- design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5)
  r_cluster <- design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5,
                            n_per_site = 25, icc = 0.03)
  expect_equal(r_cluster$deff, 1 + 24 * 0.03)         # 1.72
  expect_equal(r_cluster$n_per_timepoint,
               ceiling(nbase_trend(0.10, 0.20, 0:4) * 1.72))
  expect_gt(r_cluster$n_per_timepoint, r_srs$n_per_timepoint)
  expect_equal(r_cluster$n_eff, r_srs$n_eff)
  expect_equal(r_cluster$n_sites, ceiling(r_cluster$n_per_timepoint / 25))
})

test_that("DT-R4-2: fixed n_sites solved via closed form, deff > 1", {
  res <- design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5,
                      n_sites = 30, icc = 0.03)
  nb    <- nbase_trend(0.10, 0.20, 0:4)
  denom <- 30 - nb * 0.03
  n_cont <- nb * 30 * (1 - 0.03) / denom
  expect_equal(res$n_per_timepoint, ceiling(n_cont))
  expect_gt(res$deff, 1)
  expect_equal(res$n_sites, 30)
})

test_that("DT-R4-3: icc = 0 with cluster structure is plain SRS", {
  r_srs  <- design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5)
  r_zero <- design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5,
                         n_per_site = 25, icc = 0)
  expect_equal(r_zero$deff, 1)
  expect_equal(r_zero$n_per_timepoint, r_srs$n_per_timepoint)
})

test_that("DT-R4-4: n_sites >= per-timepoint SRS size is rejected", {
  expect_error(
    design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5,
                 n_sites = 5000, icc = 0.03),
    "not a valid cluster design"
  )
})


# ---------------------------------------------------------------------------
# Round 5 -- finite-population correction
# ---------------------------------------------------------------------------

test_that("DT-R5-1: FPC shrinks per-timepoint n, leaves n_eff alone", {
  r_inf <- design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5)
  r_fpc <- design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5, fpc_N = 400)
  nb     <- nbase_trend(0.10, 0.20, 0:4)
  n_cont <- nb * 400 / (nb + 400 - 1)
  expect_equal(r_fpc$n_per_timepoint, ceiling(n_cont))
  expect_lt(r_fpc$n_per_timepoint, r_inf$n_per_timepoint)
  expect_equal(r_fpc$n_eff, r_inf$n_eff)
})


# ---------------------------------------------------------------------------
# Round 6 -- input validation
# ---------------------------------------------------------------------------

test_that("DT-R6-1: prevalence_start guards", {
  expect_error(design_trend(0, prevalence_end = 0.2, n_timepoints = 5),
               "strictly between 0 and 1")
  expect_error(design_trend("0.1", prevalence_end = 0.2, n_timepoints = 5),
               "must be a single number")
})

test_that("DT-R6-2: exactly one of prevalence_end / slope", {
  expect_error(design_trend(0.1, n_timepoints = 5),
               "exactly one of")
  expect_error(design_trend(0.1, prevalence_end = 0.2, slope = 0.02, n_timepoints = 5),
               "not both")
})

test_that("DT-R6-3: slope that leaves (0, 1) is rejected", {
  expect_error(design_trend(0.1, slope = 0.5, n_timepoints = 5),
               "outside \\(0, 1\\)")
})

test_that("DT-R6-4: timepoints guards", {
  expect_error(design_trend(0.1, prevalence_end = 0.2),
               "n_timepoints` must be a single integer >= 2")
  expect_error(design_trend(0.1, prevalence_end = 0.2, n_timepoints = 1),
               "integer >= 2")
  expect_error(design_trend(0.1, prevalence_end = 0.2, times = c(3, 3, 3)),
               "positive range")
  expect_error(design_trend(0.1, prevalence_end = 0.2, times = 2),
               "length >= 2")
})

test_that("DT-R6-5: direction checks", {
  expect_error(
    design_trend(0.20, prevalence_end = 0.10, n_timepoints = 5,
                 alternative = "greater"),
    "increasing trend"
  )
  expect_error(
    design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5,
                 alternative = "less"),
    "decreasing trend"
  )
})

test_that("DT-R6-6: power / n / icc guards", {
  expect_error(design_trend(0.1, prevalence_end = 0.2, n_timepoints = 5, power = 0.3),
               "at least 0.5")
  expect_error(design_trend(0.1, prevalence_end = 0.2, n_timepoints = 5, n = 2.5),
               "positive integer")
  expect_error(design_trend(0.1, prevalence_end = 0.2, n_timepoints = 5, icc = 0.05),
               "requires a cluster structure")
})


# ---------------------------------------------------------------------------
# Round 7 -- decreasing trend
# ---------------------------------------------------------------------------

test_that("DT-R7-1: decreasing trend is symmetric to the increasing one", {
  up   <- design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5)
  down <- design_trend(0.20, prevalence_end = 0.10, n_timepoints = 5)
  expect_equal(up$n_per_timepoint, down$n_per_timepoint)
  expect_equal(down$slope, -0.025, tolerance = 1e-12)
  expect_lt(down$slope_app, 0)
})

test_that("DT-R7-2: 'less' accepts a genuine decline", {
  res <- design_trend(0.20, prevalence_end = 0.10, n_timepoints = 5,
                      alternative = "less")
  expect_equal(res$alternative, "less")
  expect_true(is.finite(res$n_per_timepoint))
})


# ---------------------------------------------------------------------------
# Round 8 -- code-review fixes (2026-09-02): reverse-mode consistency
# ---------------------------------------------------------------------------

test_that("DT-R8-1: reverse mode rejects n below the cluster layout", {
  expect_error(
    design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5,
                 n = 40, n_sites = 80, icc = 0.03),
    "smaller than `n_sites`"
  )
  expect_error(
    design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5,
                 n = 10, n_per_site = 25, icc = 0.03),
    "smaller than `n_per_site`"
  )
})

test_that("DT-R8-2: reverse-mode one-sided power has no spurious opposite tail", {
  r <- design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5,
                    n = 120, alternative = "greater")
  # recompute the non-centrality by hand and compare to the one-sided form
  p0a <- 0.10; p1a <- 0.20
  sxx <- 10; sigma2 <- 0.15 * 0.85
  slope_app <- (p1a - p0a) / 4
  ncp <- abs(slope_app) * sqrt(120 * sxx / sigma2)
  z_a <- qnorm(0.95)
  expect_equal(r$power, pnorm(ncp - z_a), tolerance = 1e-9)
  # the two-sided call at the same n is strictly less powered
  r2 <- design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5,
                     n = 120, alternative = "two.sided")
  expect_lt(r2$power, r$power)
})

test_that("DT-R8-3: forward/reverse round-trip agrees with FPC + n_per_site", {
  fwd <- design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5,
                      power = 0.80, n_per_site = 25, icc = 0.03, fpc_N = 600)
  rev <- design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5,
                      n = fwd$n_per_timepoint,
                      n_per_site = 25, icc = 0.03, fpc_N = 600)
  # power recovered from the design's own n should sit at ~the target
  expect_equal(rev$power, 0.80, tolerance = 0.03)
})

test_that("DT-R8-4: forward/reverse round-trip agrees with FPC + n_sites", {
  # the case the fresh review flagged: n_sites (circular solve) + fpc_N
  fwd <- design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5,
                      power = 0.80, n_sites = 20, icc = 0.05, fpc_N = 250)
  rev <- design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5,
                      n = fwd$n_per_timepoint,
                      n_sites = 20, icc = 0.05, fpc_N = 250)
  expect_equal(rev$power, 0.80, tolerance = 0.03)   # was ~0.88 before the fix
  expect_equal(rev$deff, fwd$deff, tolerance = 0.05)  # modes describe one design alike
})

test_that("DT-R8-5: fpc_N must be a whole number; `times` + `n_timepoints` warns", {
  expect_error(design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5,
                            fpc_N = 500.5), "positive integer")
  expect_silent(design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5,
                             fpc_N = 500))
  expect_warning(
    design_trend(0.10, prevalence_end = 0.20, times = c(0, 1, 2), n_timepoints = 9),
    "n_timepoints` is ignored"
  )
})

test_that("DT-R8-6: `n_per_timepoint` is numeric in both modes", {
  fwd <- design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5)
  rev <- design_trend(0.10, prevalence_end = 0.20, n_timepoints = 5, n = 150)
  expect_type(fwd$n_per_timepoint, "double")
  expect_type(rev$n_per_timepoint, "double")
})
