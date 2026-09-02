# Hand-checked reference cases for test_trend()
#
# Weighted least squares of apparent prevalence on time:
#   (t_j, X_j, N_j) aggregated per unique time
#   p_j  = X_j / N_j          (weight plug-in (X_j+.5)/(N_j+1) when X_j in {0, N_j})
#   w_j  = 1 / (p_w*(1-p_w)/N_j * deff * fpc2_j)
#   beta_app = sum(w_j (t_j - t_w)(p_j - p_w)) / sum(w_j (t_j - t_w)^2)
#   Var(beta_app) = 1 / sum(w_j (t_j - t_w)^2)
#   slope (true) = beta_app / (Se + Sp - 1)

# independent reimplementation for cross-checking
wls_ref <- function(x, n, tt, se = 1, sp = 1, deff = 1) {
  tj <- sort(unique(tt))
  Xj <- vapply(tj, function(t) sum(x[tt == t]), numeric(1))
  Nj <- vapply(tj, function(t) sum(n[tt == t]), numeric(1))
  pj <- Xj / Nj
  pw <- ifelse(Xj == 0 | Xj == Nj, (Xj + 0.5) / (Nj + 1), pj)
  wj <- 1 / (pw * (1 - pw) / Nj * deff)
  tw <- sum(wj * tj) / sum(wj)
  pwm <- sum(wj * pj) / sum(wj)
  sxx <- sum(wj * (tj - tw)^2)
  b   <- sum(wj * (tj - tw) * (pj - pwm)) / sxx
  corr <- se + sp - 1
  list(slope = b / corr, slope_app = b, se = sqrt(1 / sxx) / corr,
       z = b / sqrt(1 / sxx))
}


# ---------------------------------------------------------------------------
# Round 1 -- core WLS slope, perfect test, one survey per timepoint
# ---------------------------------------------------------------------------

test_that("TT-R1-1: matches an independent WLS computation", {
  x <- c(5, 8, 9, 14, 17); n <- rep(100, 5); tt <- 2019:2023
  res <- test_trend(x, n, tt)
  ref <- wls_ref(x, n, tt)

  expect_equal(res$slope_app, ref$slope_app, tolerance = 1e-9)
  expect_equal(res$slope,     ref$slope,     tolerance = 1e-9)
  expect_equal(res$statistic, ref$z,         tolerance = 1e-9)
  expect_equal(res$ci_lower, ref$slope - qnorm(0.975) * ref$se, tolerance = 1e-9)
  expect_equal(res$ci_upper, ref$slope + qnorm(0.975) * ref$se, tolerance = 1e-9)
  expect_equal(res$method, "wls")
  expect_equal(res$times, 2019:2023)
  expect_equal(res$n_timepoints, 5)
  expect_equal(res$n_total, 500)
})

test_that("TT-R1-2: a clear rising trend is detected", {
  res <- test_trend(x = c(3, 6, 10, 15, 22), n = rep(150, 5), time = 1:5)
  expect_gt(res$slope, 0)
  expect_true(res$reject)
  expect_lt(res$p_value, 0.01)
  expect_gt(res$ci_lower, 0)              # CI entirely positive
})

test_that("TT-R1-3: flat data gives ~zero slope, large p, CI spanning 0", {
  res <- test_trend(x = c(10, 11, 9, 10, 10), n = rep(100, 5), time = 1:5)
  expect_lt(abs(res$slope), 0.01)
  expect_gt(res$p_value, 0.3)
  expect_false(res$reject)
  expect_lt(res$ci_lower, 0)
  expect_gt(res$ci_upper, 0)
})

test_that("TT-R1-4: fitted endpoints bracket the observed trend direction", {
  res <- test_trend(x = c(4, 7, 12, 16), n = rep(120, 4), time = 1:4)
  expect_lt(res$prevalence_start_est, res$prevalence_end_est)
  expect_gt(res$prevalence_start_est, 0)
  expect_lt(res$prevalence_end_est, 1)
})

test_that("TT-R1-5: default time is seq_along(x)", {
  a <- test_trend(x = c(5, 8, 12), n = rep(100, 3))
  b <- test_trend(x = c(5, 8, 12), n = rep(100, 3), time = 1:3)
  expect_equal(a$slope, b$slope)
  expect_equal(a$times, 1:3)
})

test_that("TT-R1-6: reversing time flips the slope sign, keeps |z|", {
  x <- c(3, 6, 10, 15, 22); n <- rep(150, 5)
  up   <- test_trend(x, n, time = 1:5)
  down <- test_trend(x, n, time = 5:1)
  expect_equal(up$slope, -down$slope, tolerance = 1e-9)
  expect_equal(abs(up$statistic), abs(down$statistic), tolerance = 1e-9)
})


# ---------------------------------------------------------------------------
# Round 2 -- one-sided alternatives
# ---------------------------------------------------------------------------

test_that("TT-R2-1: one-sided p-values relate to the two-sided one", {
  x <- c(3, 6, 10, 15, 22); n <- rep(150, 5); tt <- 1:5
  two <- test_trend(x, n, tt, alternative = "two.sided")
  gt  <- test_trend(x, n, tt, alternative = "greater")
  lt  <- test_trend(x, n, tt, alternative = "less")
  expect_equal(gt$p_value, two$p_value / 2, tolerance = 1e-9)   # slope > 0 observed
  expect_equal(lt$p_value, 1 - gt$p_value, tolerance = 1e-9)
  expect_true(gt$reject)
  expect_false(lt$reject)
})

test_that("TT-R2-2: the slope CI is matched to `alternative`", {
  x <- c(3, 6, 10, 15, 22); n <- rep(150, 5); tt <- 1:5
  two <- test_trend(x, n, tt, alternative = "two.sided")
  gt  <- test_trend(x, n, tt, alternative = "greater")
  lt  <- test_trend(x, n, tt, alternative = "less")

  expect_equal(gt$ci_upper, Inf)     # [L, Inf)
  expect_equal(lt$ci_lower, -Inf)    # (-Inf, U]
  # one-sided bound uses z_{1-alpha}, tighter than two-sided z_{1-alpha/2}
  expect_gt(gt$ci_lower, two$ci_lower)
  expect_lt(lt$ci_upper, two$ci_upper)
  # rising trend rejected "greater" -> the one-sided lower bound clears 0
  expect_true(gt$reject)
  expect_gt(gt$ci_lower, 0)
})


# ---------------------------------------------------------------------------
# Round 3 -- Rogan-Gladen correction
# ---------------------------------------------------------------------------

test_that("TT-R3-1: imperfect test rescales the slope by 1 / (Se + Sp - 1)", {
  x <- c(5, 8, 9, 14, 17); n <- rep(100, 5); tt <- 1:5
  a <- test_trend(x, n, tt)
  b <- test_trend(x, n, tt, sensitivity = 0.9, specificity = 0.95)
  corr <- 0.9 + 0.95 - 1
  expect_equal(b$slope_app, a$slope_app, tolerance = 1e-12)   # apparent fit unchanged
  expect_equal(b$slope, a$slope_app / corr, tolerance = 1e-9)
  expect_equal(b$statistic, a$statistic, tolerance = 1e-12)   # z invariant
  # CI on the true scale is the apparent-scale CI divided by the correction
  expect_equal(b$ci_lower * corr, a$ci_lower, tolerance = 1e-9)
  expect_equal(b$ci_upper * corr, a$ci_upper, tolerance = 1e-9)
})

test_that("TT-R3-2: se + sp <= 1 is rejected", {
  expect_error(
    test_trend(x = c(5, 8, 12), n = rep(100, 3), time = 1:3,
               sensitivity = 0.5, specificity = 0.5),
    "must exceed 1"
  )
})


# ---------------------------------------------------------------------------
# Round 4 -- clustering (multiple sites per timepoint)
# ---------------------------------------------------------------------------

test_that("TT-R4-1: between-site spread inflates deff and widens the CI", {
  # 3 sites per round, 4 rounds; sites disagree within a round
  x  <- c(1, 8, 3,  2, 10, 4,  4, 12, 6,  6, 15, 8)
  n  <- rep(50, 12)
  tt <- rep(1:4, each = 3)

  r_srs  <- test_trend(x, n, tt, icc = 0)
  r_clus <- test_trend(x, n, tt)            # icc estimated from data

  expect_equal(r_srs$deff, 1)
  expect_gt(r_clus$deff, 1)
  expect_gt(r_clus$icc_used, 0)
  # wider interval once clustering is acknowledged
  expect_gt(r_clus$ci_upper - r_clus$ci_lower, r_srs$ci_upper - r_srs$ci_lower)
  # slope point estimate is unchanged by a uniform weight rescale
  expect_equal(r_clus$slope, r_srs$slope, tolerance = 1e-9)
})

test_that("TT-R4-2: a supplied scalar icc is used directly", {
  x  <- c(2, 6, 3, 7, 4, 9, 6, 11)
  n  <- rep(40, 8)
  tt <- rep(1:4, each = 2)
  res <- test_trend(x, n, tt, icc = 0.05)
  expect_equal(res$icc_used, 0.05)
  expect_equal(res$deff, 1 + (mean(n) - 1) * 0.05)
})

test_that("TT-R4-3: one row per timepoint => deff = 1 even with icc = NULL", {
  res <- test_trend(x = c(5, 8, 9, 14), n = rep(100, 4), time = 1:4)
  expect_equal(res$deff, 1)
  expect_equal(res$icc_used, 0)
})


# ---------------------------------------------------------------------------
# Round 5 -- finite-population correction
# ---------------------------------------------------------------------------

test_that("TT-R5-1: fpc_N tightens the slope CI", {
  x <- c(5, 8, 9, 14, 17); n <- rep(100, 5); tt <- 1:5
  r_inf <- test_trend(x, n, tt)
  r_fpc <- test_trend(x, n, tt, fpc_N = 300)
  expect_lt(r_fpc$ci_upper - r_fpc$ci_lower, r_inf$ci_upper - r_inf$ci_lower)
})

test_that("TT-R5-2: fpc_N <= the largest timepoint sample errors", {
  expect_error(
    test_trend(x = c(5, 8, 9), n = c(100, 100, 100), time = 1:3, fpc_N = 80),
    "greater than the largest"
  )
  expect_error(
    test_trend(x = c(5, 8, 9), n = c(100, 100, 100), time = 1:3, fpc_N = 100),
    "greater than the largest"
  )
})


# ---------------------------------------------------------------------------
# Round 6 -- input validation
# ---------------------------------------------------------------------------

test_that("TT-R6-1: x / n vector guards", {
  expect_error(test_trend(x = c(5, 8), n = c(100, 100, 100), time = 1:2),
               "same length")
  expect_error(test_trend(x = c(5.5, 8), n = c(100, 100), time = 1:2),
               "whole numbers")
  expect_error(test_trend(x = c(120, 8), n = c(100, 100), time = 1:2),
               "cannot exceed")
  expect_error(test_trend(x = c(5, 8), n = c(0, 100), time = 1:2),
               "must be positive")
  expect_error(test_trend(x = NA, n = 100, time = 1),
               "logical")
})

test_that("TT-R6-2: time guards", {
  expect_error(test_trend(x = c(5, 8, 9), n = rep(100, 3), time = c(1, 2)),
               "same length as `x`")
  expect_error(test_trend(x = c(5, 8, 9), n = rep(100, 3), time = c(2, 2, 2)),
               "at least two distinct")
})

test_that("TT-R6-3: scalar-parameter guards", {
  x <- c(5, 8, 9); n <- rep(100, 3); tt <- 1:3
  expect_error(test_trend(x, n, tt, alternative = "up"),
               "'greater', 'less', or 'two.sided'")
  expect_error(test_trend(x, n, tt, sensitivity = TRUE),
               "must be a single number")
  expect_error(test_trend(x, n, tt, specificity = 1.3),
               "must be in \\(0, 1\\]")
  expect_error(test_trend(x, n, tt, conf_level = 1),
               "strictly between 0 and 1")
  expect_error(test_trend(x, n, tt, icc = 1.4),
               "in \\[0, 1\\]")
  expect_error(test_trend(x, n, tt, icc = "x"),
               "NULL or a single number")
  expect_error(test_trend(x, n, tt, ci_method = "ols"),
               "only 'wls'")
  expect_error(test_trend(x, n, tt, fpc_N = -3),
               "finite positive number")
})


# ---------------------------------------------------------------------------
# Round 7 -- round-trip against design_trend()
# ---------------------------------------------------------------------------

test_that("TT-R7-1: at design_trend()'s n, data on the assumed line rejects", {
  d <- design_trend(prevalence_start = 0.10, prevalence_end = 0.20,
                    n_timepoints = 5, power = 0.80)
  n_pt <- d$n_per_timepoint
  tt   <- 0:4
  # prevalence exactly on the assumed straight line at each round
  p_line <- 0.10 + (0.20 - 0.10) / 4 * tt
  x <- round(p_line * n_pt)
  res <- test_trend(x = x, n = rep(n_pt, 5), time = tt)
  expect_true(res$reject)
  expect_lt(res$p_value, 0.05)
  expect_equal(res$slope, 0.025, tolerance = 0.02)   # recovers ~the design slope
})

test_that("TT-R7-2: far below that n, the trend is not detected", {
  d <- design_trend(prevalence_start = 0.10, prevalence_end = 0.20,
                    n_timepoints = 5, power = 0.80)
  n_pt <- max(5, round(d$n_per_timepoint / 6))
  tt   <- 0:4
  p_line <- 0.10 + 0.025 * tt
  x <- round(p_line * n_pt)
  res <- test_trend(x = x, n = rep(n_pt, 5), time = tt)
  expect_false(res$reject)
})


# ---------------------------------------------------------------------------
# Round 8 -- code-review fixes (2026-09-02): design effect only with a
# real cluster structure
# ---------------------------------------------------------------------------

test_that("TT-R8-1: explicit icc is ignored with one row per timepoint", {
  # no cluster structure -> the row n is the timepoint total, not a
  # cluster size, so deff must stay 1 (matches test_threshold convention)
  res <- test_trend(x = c(5, 8, 9, 14, 17), n = rep(100, 5),
                    time = 1:5, icc = 0.02)
  expect_equal(res$deff, 1)
  expect_equal(res$icc_used, 0.02)     # echoed back, but not applied

  # identical result to not passing icc at all
  plain <- test_trend(x = c(5, 8, 9, 14, 17), n = rep(100, 5), time = 1:5)
  expect_equal(res$statistic, plain$statistic)
  expect_equal(res$ci_lower, plain$ci_lower)
})

test_that("TT-R8-2: deff uses the clustered-rows cluster size, not the global mean", {
  # timepoints 1-2: three sites of 50; timepoint 3: one pooled row of 150
  x  <- c(2, 5, 3,  4, 7, 5,  18)
  n  <- c(50, 50, 50, 50, 50, 50, 150)
  tt <- c(1, 1, 1, 2, 2, 2, 3)
  res <- test_trend(x, n, tt, icc = 0.04)
  # cluster size is 50 (the multi-row timepoints), not mean(n) = 64.3
  expect_equal(res$deff, 1 + (50 - 1) * 0.04)
})

test_that("TT-R8-3: multi-site data still forms deff > 1 from between-site spread", {
  x  <- c(1, 8, 3,  2, 10, 4,  4, 12, 6,  6, 15, 8)
  n  <- rep(50, 12)
  tt <- rep(1:4, each = 3)
  res <- test_trend(x, n, tt)             # icc estimated
  expect_gt(res$deff, 1)
  expect_gt(res$icc_used, 0)
})
