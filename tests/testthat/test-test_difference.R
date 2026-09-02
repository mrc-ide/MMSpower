# Hand-checked reference cases for test_difference()
#
# Pooled-variance two-proportion z-test on the apparent scale:
#   p_hat_g = sum(x_g) / sum(n_g)
#   m_g     = n_total_g / deff_g / fpc_g^2       [effective sample size]
#   p_pool  = (p1*m1 + p2*m2) / (m1 + m2)
#   se_null = sqrt(p_pool*(1-p_pool)*(1/m1 + 1/m2))
#   z       = (p1 - p2) / se_null
#
# CI on the difference: built on the apparent scale, then / (Se + Sp - 1).


# ---------------------------------------------------------------------------
# Round 1 -- core two-proportion z-test, perfect test, SRS
# ---------------------------------------------------------------------------

test_that("TD-R1-1: two single sites, perfect test -- hand-checked z and p", {
  res <- test_difference(x1 = 8, n1 = 100, x2 = 20, n2 = 120)

  p1 <- 0.08; p2 <- 20 / 120
  m1 <- 100;  m2 <- 120
  p_pool  <- (p1 * m1 + p2 * m2) / (m1 + m2)
  se_null <- sqrt(p_pool * (1 - p_pool) * (1 / m1 + 1 / m2))
  z_exp   <- (p1 - p2) / se_null

  expect_equal(res$statistic, z_exp, tolerance = 1e-9)
  expect_equal(res$p_value, 2 * pnorm(abs(z_exp), lower.tail = FALSE),
               tolerance = 1e-9)
  expect_equal(res$difference_app, p1 - p2, tolerance = 1e-12)
  expect_equal(res$difference, p1 - p2, tolerance = 1e-12)   # perfect test
  expect_equal(res$deff1, 1)
  expect_equal(res$deff2, 1)
  expect_equal(res$n_eff_adj1, 100)
  expect_equal(res$n_eff_adj2, 120)
  expect_false(res$reject)          # p ~ 0.0547
})

test_that("TD-R1-2: agrees with prop.test (no continuity correction)", {
  res <- test_difference(x1 = 30, n1 = 200, x2 = 15, n2 = 180)
  pt  <- prop.test(c(30, 15), c(200, 180), correct = FALSE)

  # prop.test reports chi-square = z^2
  expect_equal(res$statistic^2, unname(pt$statistic), tolerance = 1e-8)
  expect_equal(res$p_value, pt$p.value, tolerance = 1e-8)
})

test_that("TD-R1-3: one-sided p-values are half the two-sided (sign permitting)", {
  two  <- test_difference(x1 = 30, n1 = 200, x2 = 15, n2 = 180,
                          alternative = "two.sided")
  gt   <- test_difference(x1 = 30, n1 = 200, x2 = 15, n2 = 180,
                          alternative = "greater")
  lt   <- test_difference(x1 = 30, n1 = 200, x2 = 15, n2 = 180,
                          alternative = "less")
  expect_equal(gt$p_value, two$p_value / 2, tolerance = 1e-9)  # p1 > p2 observed
  expect_equal(lt$p_value, 1 - gt$p_value, tolerance = 1e-9)
})

test_that("TD-R1-4: identical groups give z = 0, p = 1, CI spanning 0", {
  res <- test_difference(x1 = 10, n1 = 100, x2 = 10, n2 = 100)
  expect_equal(res$statistic, 0)
  expect_equal(res$p_value, 1)
  expect_false(res$reject)
  expect_lt(res$ci_lower, 0)
  expect_gt(res$ci_upper, 0)
})

test_that("TD-R1-5: a large clear difference is rejected", {
  res <- test_difference(x1 = 5, n1 = 200, x2 = 60, n2 = 200)
  expect_true(res$reject)
  expect_lt(res$p_value, 1e-6)
  expect_lt(res$ci_upper, 0)        # group1 < group2, CI entirely negative
})


# ---------------------------------------------------------------------------
# Round 2 -- confidence interval on the difference
# ---------------------------------------------------------------------------

test_that("TD-R2-1: Wald CI matches the closed form", {
  res <- test_difference(x1 = 8, n1 = 100, x2 = 20, n2 = 120)
  p1 <- 0.08; p2 <- 20 / 120
  se <- sqrt(p1 * (1 - p1) / 100 + p2 * (1 - p2) / 120)
  z  <- qnorm(0.975)
  expect_equal(res$ci_lower, (p1 - p2) - z * se, tolerance = 1e-9)
  expect_equal(res$ci_upper, (p1 - p2) + z * se, tolerance = 1e-9)
})

test_that("TD-R2-2: Newcombe CI differs from Wald and stays in [-1, 1]", {
  w <- test_difference(x1 = 3, n1 = 80, x2 = 12, n2 = 90, ci_method = "wald")
  n <- test_difference(x1 = 3, n1 = 80, x2 = 12, n2 = 90, ci_method = "newcombe")
  expect_equal(n$ci_method, "newcombe")
  expect_false(isTRUE(all.equal(w$ci_lower, n$ci_lower)))
  expect_gte(n$ci_lower, -1)
  expect_lte(n$ci_upper, 1)
})

test_that("TD-R2-3: wider conf_level gives a wider interval", {
  r95 <- test_difference(x1 = 8, n1 = 100, x2 = 20, n2 = 120, conf_level = 0.95)
  r99 <- test_difference(x1 = 8, n1 = 100, x2 = 20, n2 = 120, conf_level = 0.99)
  expect_lt(r99$ci_lower, r95$ci_lower)
  expect_gt(r99$ci_upper, r95$ci_upper)
})


# ---------------------------------------------------------------------------
# Round 3 -- Rogan-Gladen correction
# ---------------------------------------------------------------------------

test_that("TD-R3-1: imperfect test rescales the difference by (Se + Sp - 1)", {
  se <- 0.9; sp <- 0.95; corr <- se + sp - 1
  res <- test_difference(x1 = 8, n1 = 100, x2 = 20, n2 = 120,
                         sensitivity = se, specificity = sp)

  d_app <- 0.08 - 20 / 120
  expect_equal(res$difference_app, d_app, tolerance = 1e-12)
  expect_equal(res$difference, d_app / corr, tolerance = 1e-9)
  expect_equal(res$ci_lower * corr,
               test_difference(x1 = 8, n1 = 100, x2 = 20, n2 = 120)$ci_lower,
               tolerance = 1e-9)
})

test_that("TD-R3-2: RG rescaling does not change the z-statistic or p-value", {
  a <- test_difference(x1 = 8, n1 = 100, x2 = 20, n2 = 120)
  b <- test_difference(x1 = 8, n1 = 100, x2 = 20, n2 = 120,
                       sensitivity = 0.9, specificity = 0.95)
  expect_equal(a$statistic, b$statistic, tolerance = 1e-12)
  expect_equal(a$p_value,   b$p_value,   tolerance = 1e-12)
})

test_that("TD-R3-3: se + sp <= 1 is rejected", {
  expect_error(
    test_difference(x1 = 8, n1 = 100, x2 = 20, n2 = 120,
                    sensitivity = 0.5, specificity = 0.5),
    "must exceed 1"
  )
})

test_that("TD-R3-4: with a clamped per-group estimate, difference stays in its CI", {
  # group 1: 1/100 apparent -> rg(0.01) < 0 with Sp = 0.97, clamps to 0
  res <- test_difference(x1 = 1, n1 = 100, x2 = 8, n2 = 100,
                         sensitivity = 0.9, specificity = 0.97)
  corr <- 0.9 + 0.97 - 1
  expect_equal(res$prevalence1, 0)                       # clamped
  # difference is d_app / correction, NOT prevalence1 - prevalence2
  expect_equal(res$difference, res$difference_app / corr, tolerance = 1e-9)
  expect_false(isTRUE(all.equal(res$difference,
                                res$prevalence1 - res$prevalence2)))
  # and it lies inside the reported interval
  expect_gte(res$difference, res$ci_lower)
  expect_lte(res$difference, res$ci_upper)
})


# ---------------------------------------------------------------------------
# Round 4 -- clustering (per-group design effect)
# ---------------------------------------------------------------------------

test_that("TD-R4-1: clustered data inflates per-group deff, widens CI", {
  x1 <- c(2, 5, 3, 9); n1 <- c(50, 60, 40, 70)
  x2 <- c(8, 11, 9, 4); n2 <- c(55, 70, 45, 60)

  r_srs  <- test_difference(x1, n1, x2, n2, icc = 0)
  r_clus <- test_difference(x1, n1, x2, n2)            # icc estimated

  expect_equal(r_srs$deff1, 1)
  expect_gte(r_clus$deff1, 1)
  expect_gte(r_clus$deff2, 1)
  # at least one group shows real between-cluster spread -> deff > 1 somewhere
  expect_true(r_clus$deff1 > 1 || r_clus$deff2 > 1)
  expect_lt(r_clus$n_eff_adj1, r_srs$n_eff_adj1 + 1e-9)
})

test_that("TD-R4-2: a supplied scalar icc is applied to both groups", {
  x1 <- c(2, 5, 3); n1 <- c(50, 50, 50)
  x2 <- c(8, 6, 7); n2 <- c(40, 40, 40)
  res <- test_difference(x1, n1, x2, n2, icc = 0.05)
  expect_equal(res$icc_used1, 0.05)
  expect_equal(res$icc_used2, 0.05)
  expect_equal(res$deff1, 1 + (mean(n1) - 1) * 0.05)
  expect_equal(res$deff2, 1 + (mean(n2) - 1) * 0.05)
})

test_that("TD-R4-3: single-site groups ignore a supplied icc (deff = 1)", {
  res <- test_difference(x1 = 8, n1 = 100, x2 = 20, n2 = 120, icc = 0.1)
  expect_equal(res$deff1, 1)
  expect_equal(res$deff2, 1)
  expect_equal(res$icc_used1, 0)
})


# ---------------------------------------------------------------------------
# Round 5 -- finite-population correction
# ---------------------------------------------------------------------------

test_that("TD-R5-1: fpc_N raises the variance-equivalent n and tightens the CI", {
  r_inf <- test_difference(x1 = 30, n1 = 200, x2 = 40, n2 = 200)
  r_fpc <- test_difference(x1 = 30, n1 = 200, x2 = 40, n2 = 200, fpc_N = 500)
  # FPC lowers the sampling variance -> larger effective (variance-equivalent) n
  expect_gt(r_fpc$n_eff_adj1, r_inf$n_eff_adj1)
  expect_gt(r_fpc$n_eff_adj2, r_inf$n_eff_adj2)
  # ... and a tighter interval (difference is negative here, so lower bound rises)
  expect_gt(r_fpc$ci_lower, r_inf$ci_lower)
  expect_lt(r_fpc$ci_upper, r_inf$ci_upper)
})

test_that("TD-R5-2: length-2 fpc_N applies per group", {
  res <- test_difference(x1 = 30, n1 = 200, x2 = 40, n2 = 200,
                         fpc_N = c(400, 5000))
  # group 1 is a bigger fraction of its population -> more FPC inflation
  expect_gt(res$n_eff_adj1, res$n_eff_adj2)
  expect_gt(res$n_eff_adj1, 200)
  expect_gt(res$n_eff_adj2, 200)
})

test_that("TD-R5-3: fpc_N below a group's sample size errors", {
  expect_error(
    test_difference(x1 = 30, n1 = 200, x2 = 40, n2 = 200, fpc_N = 150),
    "less than that group"
  )
})


# ---------------------------------------------------------------------------
# Round 6 -- input validation
# ---------------------------------------------------------------------------

test_that("TD-R6-1: count-vector guards", {
  expect_error(test_difference(x1 = 8, n1 = c(100, 100), x2 = 20, n2 = 120),
               "same length")
  expect_error(test_difference(x1 = 8.5, n1 = 100, x2 = 20, n2 = 120),
               "whole numbers")
  expect_error(test_difference(x1 = 120, n1 = 100, x2 = 20, n2 = 120),
               "cannot exceed")
  expect_error(test_difference(x1 = -1, n1 = 100, x2 = 20, n2 = 120),
               "non-negative")
  expect_error(test_difference(x1 = 8, n1 = 0, x2 = 20, n2 = 120),
               "must be positive")
  expect_error(test_difference(x1 = NA, n1 = 100, x2 = 20, n2 = 120),
               "logical")
})

test_that("TD-R6-2: scalar-parameter guards", {
  expect_error(test_difference(8, 100, 20, 120, alternative = "bigger"),
               "'greater', 'less', or 'two.sided'")
  expect_error(test_difference(8, 100, 20, 120, sensitivity = TRUE),
               "must be a single number")
  expect_error(test_difference(8, 100, 20, 120, specificity = 1.4),
               "must be in \\(0, 1\\]")
  expect_error(test_difference(8, 100, 20, 120, conf_level = 1),
               "strictly between 0 and 1")
  expect_error(test_difference(8, 100, 20, 120, icc = 1.5),
               "must be in \\[0, 1\\]")
  expect_error(test_difference(8, 100, 20, 120, icc = "x"),
               "single number or NULL")
  expect_error(test_difference(8, 100, 20, 120, ci_method = "exact"),
               "'wald' or 'newcombe'")
  expect_error(test_difference(8, 100, 20, 120, fpc_N = c(1, 2, 3)),
               "length-2")
})


# ---------------------------------------------------------------------------
# Round 7 -- round-trip against design_difference()
# ---------------------------------------------------------------------------

test_that("TD-R7-1: at design_difference()'s n, data at the assumed rates rejects", {
  d <- design_difference(prevalence1 = 0.10, prevalence2 = 0.20, power = 0.80)
  n <- d$n_per_group
  # observed counts sitting exactly on the assumed prevalences
  res <- test_difference(x1 = round(0.10 * n), n1 = n,
                         x2 = round(0.20 * n), n2 = n)
  expect_true(res$reject)
  expect_lt(res$p_value, 0.05)
})

test_that("TD-R7-2: just below that n, the effect is borderline / not rejected", {
  d <- design_difference(prevalence1 = 0.10, prevalence2 = 0.20, power = 0.80)
  n <- round(d$n_per_group / 3)      # far under-powered
  res <- test_difference(x1 = round(0.10 * n), n1 = n,
                         x2 = round(0.20 * n), n2 = n)
  expect_false(res$reject)
})
