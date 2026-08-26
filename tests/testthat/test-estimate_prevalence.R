# Hand-checked reference cases for estimate_prevalence()
#
# Key formulas:
#   Apparent prevalence p_app = positives / n
#   Wilson CI on p_app with effective n = n / DEFF
#   True prevalence = (p_app - (1 - sp)) / (se + sp - 1)      [Rogan-Gladen]
#   True CI = Rogan-Gladen applied to apparent CI endpoints
#
# Case 1  (perfect test, x=30, n=100):
#   p_app = 0.300, correction = 1
#   Wilson 95% CI on p_app: [0.2189, 0.3958]
#   p_true = 0.300, true CI = same
#
# Case 2  (se=0.9, sp=0.95, x=30, n=100):
#   p_app = 0.300, correction = 0.85
#   p_true = (0.30 - 0.05) / 0.85 = 0.2941
#   true CI lower = (0.2189 - 0.05) / 0.85 = 0.1988
#   true CI upper = (0.3958 - 0.05) / 0.85 = 0.4069
#
# Case 3  (clustering, n_per_site=10, icc=0.05, x=30, n=100):
#   DEFF = 1 + (10-1)*0.05 = 1.45
#   n_eff = 100/1.45 ≈ 68.97  →  CI wider than Case 1

test_that("Case 1: perfect test returns correct prevalence and Wilson CI", {
  res <- estimate_prevalence(positives = 30, n = 100)

  expect_equal(res$prevalence,    0.3,      tolerance = 1e-6)
  expect_equal(res$apparent_prev, 0.3,      tolerance = 1e-6)
  expect_equal(res$design_effect, 1,        tolerance = 1e-6)

  # Wilson 95% CI (hand-checked: z = qnorm(0.975) ≈ 1.959964)
  expect_equal(res$ci[["lower"]], 0.2189, tolerance = 5e-4)
  expect_equal(res$ci[["upper"]], 0.3958, tolerance = 5e-4)
})

test_that("Case 2: imperfect test applies Rogan-Gladen correction", {
  res <- estimate_prevalence(30, 100, sensitivity = 0.9, specificity = 0.95)

  expect_equal(res$apparent_prev, 0.3, tolerance = 1e-6)

  # p_true = (0.30 - 0.05) / 0.85
  expect_equal(res$prevalence, 0.25 / 0.85, tolerance = 1e-6)

  # CI endpoints: Rogan-Gladen applied to apparent Wilson CI bounds
  expect_equal(res$ci["lower"], (res$apparent_ci["lower"] - 0.05) / 0.85,
    tolerance = 1e-6)
  expect_equal(res$ci["upper"], (res$apparent_ci["upper"] - 0.05) / 0.85,
    tolerance = 1e-6)

  # Approximate values for readability
  expect_equal(res$ci[["lower"]], 0.1988, tolerance = 5e-4)
  expect_equal(res$ci[["upper"]], 0.4069, tolerance = 5e-4)
})

test_that("Case 3: clustering with ICC widens CI and reports correct DEFF", {
  res_plain   <- estimate_prevalence(30, 100)
  res_cluster <- estimate_prevalence(30, 100, n_per_site = 10, icc = 0.05)

  # DEFF = 1 + (10-1)*0.05 = 1.45
  expect_equal(res_cluster$design_effect, 1.45, tolerance = 1e-6)

  # Point estimate unchanged — clustering doesn't bias the estimate
  expect_equal(res_cluster$prevalence, res_plain$prevalence, tolerance = 1e-6)

  # CI width should be wider with clustering
  width_plain   <- res_plain$ci["upper"]   - res_plain$ci["lower"]
  width_cluster <- res_cluster$ci["upper"] - res_cluster$ci["lower"]
  expect_gt(width_cluster, width_plain)
})

test_that("returns mms_estimate object with expected fields", {
  res <- estimate_prevalence(10, 50)
  expect_s3_class(res, "mms_estimate")
  expect_named(res, c("prevalence", "ci", "conf_level", "apparent_prev",
                       "apparent_ci", "n", "positives", "sensitivity",
                       "specificity", "design_effect"))
})

test_that("input validation catches bad arguments", {
  expect_error(estimate_prevalence(-1,  100), "`positives`")
  expect_error(estimate_prevalence(50,  30),  "cannot exceed")
  expect_error(estimate_prevalence(30,  100, sensitivity = 0.3, specificity = 0.3),
               "must exceed 1")
})
