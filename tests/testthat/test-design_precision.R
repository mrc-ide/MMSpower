# Hand-checked reference cases for design_precision()
#
# Formula: n = z^2 * p_app * (1 - p_app) / (moe^2 * correction^2)
#   where p_app = p*se + (1-p)*(1-sp)   [apparent prevalence]
#         correction = se + sp - 1       [Rogan-Gladen denominator]
#         z = qnorm(0.975) ≈ 1.959964
#
# Case 1  (perfect test, p=0.30, moe=0.05):
#   p_app = 0.30, correction = 1
#   n_cont = 1.959964^2 * 0.30 * 0.70 / 0.05^2 = 322.68 → ceiling → 323
#
# Case 2  (se=0.90, sp=0.95, p=0.30, moe=0.05):
#   p_app = 0.30*0.90 + 0.70*0.05 = 0.305, correction = 0.85
#   n_cont = 1.959964^2 * 0.305 * 0.695 / (0.0025 * 0.7225) = 450.82 → ceiling → 451
#
# Case 3  (n_per_site=10, icc=0.05, p=0.30, moe=0.05):
#   deff = 1 + (10-1)*0.05 = 1.45
#   n_cont = 322.68 * 1.45 = 467.89 → ceiling → 468
#   n_sites = ceiling(468 / 10) = 47
#
# Case 4  (n_sites=50, icc=0.05, p=0.30, moe=0.05):
#   denom = 50 - 322.68*0.05 = 33.87
#   n_cont = 322.68 * 50 * 0.95 / 33.87 = 452.6 → ceiling → 453
#   n_per_site = ceiling(453 / 50) = 10
#
# Case 5  (n_sites=10, icc=0.05, p=0.30, moe=0.05):
#   denom = 10 - 322.68*0.05 = -6.13 < 0 → infeasible
#   min_moe = z * sqrt(0.3*0.7*0.05 / (10*1)) = 1.959964*sqrt(0.00105) ≈ 6.35%

test_that("Case 1: perfect test, SRS returns correct n and deff", {
  res <- design_precision(prevalence = 0.3, moe = 0.05)

  expect_equal(res$n,            323)
  expect_equal(res$deff,         1)
  expect_equal(res$apparent_prev, 0.3, tolerance = 1e-6)
  expect_null(res$n_sites)
  expect_null(res$n_per_site)
})

test_that("Case 2: imperfect test inflates n via Rogan-Gladen variance", {
  res <- design_precision(0.3, 0.05, sensitivity = 0.9, specificity = 0.95)

  # apparent prevalence = 0.30*0.90 + 0.70*0.05 = 0.305
  expect_equal(res$apparent_prev, 0.305, tolerance = 1e-6)
  expect_equal(res$n, 451)
  expect_gt(res$n, 323)   # imperfect test always needs more samples
})

test_that("Case 3: fixed n_per_site inflates n by deff = 1.45", {
  res <- design_precision(0.3, 0.05, n_per_site = 10, icc = 0.05)

  expect_equal(res$deff, 1.45, tolerance = 1e-6)
  expect_equal(res$n, 468)
  expect_gt(res$n, 323)

  # sites needed = ceiling(468 / 10) = 47
  expect_equal(res$n_sites,    47)
  expect_equal(res$n_per_site, 10)
})

test_that("Case 4: fixed n_sites resolves circularity via closed-form", {
  res <- design_precision(0.3, 0.05, n_sites = 50, icc = 0.05)

  expect_equal(res$n,         453)
  expect_equal(res$n_sites,    50)
  expect_equal(res$n_per_site, 10)   # ceiling(453/50) = 10
})

test_that("Case 5: infeasible n_sites produces informative error", {
  # denom = 10 - 322.68*0.05 < 0 → unachievable
  expect_error(
    design_precision(0.3, 0.05, n_sites = 10, icc = 0.05),
    "unachievable"
  )
})

test_that("return list contains all expected fields", {
  res <- design_precision(0.2, 0.05)
  expect_named(res, c("n", "n_base", "n_sites", "n_per_site", "prevalence",
                       "apparent_prev", "moe", "conf_level", "sensitivity",
                       "specificity", "icc", "deff"),
               ignore.order = FALSE)
})

test_that("input validation catches bad arguments", {
  expect_error(design_precision(0,   0.05),  "`prevalence`")
  expect_error(design_precision(0.3, 0),     "`moe`")
  expect_error(design_precision(0.3, 0.05, sensitivity = 0.2, specificity = 0.2),
               "must exceed 1")
  expect_error(design_precision(0.3, 0.05, icc = 0.05),
               "cluster structure")
  expect_error(design_precision(0.3, 0.05, n_sites = 50, n_per_site = 10, icc = 0.05),
               "at most one")
})
