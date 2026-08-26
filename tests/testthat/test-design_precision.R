# Hand-checked reference cases for design_precision()
#
# Formula: n = z^2 * p_app * (1 - p_app) / (moe^2 * correction^2)
#   where p_app = p*se + (1-p)*(1-sp)   [apparent prevalence]
#         correction = se + sp - 1       [Rogan-Gladen denominator]
#         z = qnorm(0.975) ≈ 1.959964
#
# Case 1  (perfect test, p=0.30, moe=0.05):
#   p_app = 0.30, correction = 1
#   n = 1.959964^2 * 0.30 * 0.70 / 0.05^2 = 322.68 → ceiling → 323
#
# Case 2  (se=0.90, sp=0.95, p=0.30, moe=0.05):
#   p_app = 0.30*0.90 + 0.70*0.05 = 0.305, correction = 0.85
#   n = 1.959964^2 * 0.305 * 0.695 / (0.0025 * 0.7225) = 450.82 → ceiling → 451
#
# Case 3  (n_per_site=10, icc=0.05, p=0.30, moe=0.05):
#   DEFF = 1 + (10-1)*0.05 = 1.45
#   n_continuous = 322.68 * 1.45 = 467.89 → ceiling → 468

test_that("Case 1: perfect test returns correct n", {
  res <- design_precision(prevalence = 0.3, moe = 0.05)

  expect_equal(res$n,             323)
  expect_equal(res$design_effect, 1)
  expect_equal(res$apparent_prev, 0.3, tolerance = 1e-6)
})

test_that("Case 2: imperfect test inflates n via Rogan-Gladen variance", {
  res <- design_precision(0.3, 0.05, sensitivity = 0.9, specificity = 0.95)

  # apparent prevalence = 0.30*0.90 + 0.70*0.05 = 0.305
  expect_equal(res$apparent_prev, 0.305, tolerance = 1e-6)

  expect_equal(res$n, 451)
  expect_gt(res$n, 323)   # imperfect test always needs more samples
})

test_that("Case 3: clustering inflates n by DEFF = 1.45", {
  res <- design_precision(0.3, 0.05, n_per_site = 10, icc = 0.05)

  expect_equal(res$design_effect, 1.45, tolerance = 1e-6)
  expect_equal(res$n, 468)
  expect_gt(res$n, 323)

  # sites needed = ceiling(468 / 10) = 47
  expect_equal(res$n_sites, 47)
})

test_that("returns mms_design object with expected fields", {
  res <- design_precision(0.2, 0.05)
  expect_s3_class(res, "mms_design")
  expect_named(res, c("n", "n_base", "n_sites", "n_per_site", "prevalence",
                       "apparent_prev", "moe", "conf_level", "sensitivity",
                       "specificity", "design_effect"))
})

test_that("input validation catches bad arguments", {
  expect_error(design_precision(0,    0.05), "`prevalence`")
  expect_error(design_precision(0.3,  0),    "`moe`")
  expect_error(design_precision(0.3,  0.05, sensitivity = 0.2, specificity = 0.2),
               "must exceed 1")
})
