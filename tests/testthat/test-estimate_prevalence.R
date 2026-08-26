# Hand-checked reference cases for estimate_prevalence()
#
# Formula: Wald CI on apparent prevalence, then Rogan-Gladen to true scale.
#
#   p_hat   = sum(x) / sum(n)
#   n_eff   = n_total / deff
#   se      = sqrt(p_hat*(1-p_hat) / n_eff)
#   moe_app = z * se          where z = qnorm(0.975) = 1.959964
#   CI_app  = [p_hat - moe_app, p_hat + moe_app]
#   p_true  = (p_hat - (1-sp)) / (se_test + sp - 1)    [Rogan-Gladen]
#   CI_true = Rogan-Gladen applied to CI_app endpoints
#
# Case 1  (single site, perfect test, x=30, n=100):
#   p_hat = 0.3, deff = 1, n_eff = 100
#   se = sqrt(0.3*0.7/100) = 0.04583
#   moe = 1.959964 * 0.04583 = 0.08982
#   CI = [0.2102, 0.3898]
#
# Case 2  (supplied icc=0.05, 10 uniform sites of 10 each):
#   deff = 1 + 9*0.05 = 1.45, n_eff = 100/1.45 = 68.97
#   se = sqrt(0.3*0.7/68.97) = 0.05519
#   moe = 1.959964 * 0.05519 = 0.10817
#   CI = [0.1918, 0.4082]
#
# Case 3  (imperfect test, se=0.9, sp=0.95, x=30, n=100):
#   p_apparent = 0.3, correction = 0.85
#   p_true = (0.3 - 0.05) / 0.85 = 0.2941
#   moe_apparent = 0.08982  (same as Case 1)
#   moe_true = 0.08982 / 0.85 = 0.10567
#   CI_true = [(0.2102-0.05)/0.85, (0.3898-0.05)/0.85]
#           = [0.1885, 0.3998]

test_that("Case 1: single site, perfect test -- Wald CI", {
  res <- estimate_prevalence(x = 30, n = 100)

  expect_equal(res$prevalence,     0.3,    tolerance = 1e-6)
  expect_equal(res$deff,           1,      tolerance = 1e-6)
  expect_equal(res$icc_used,       0,      tolerance = 1e-6)
  expect_equal(res$n_total,        100)
  expect_equal(res$n_eff,          100,    tolerance = 1e-6)

  # Wald 95% CI (hand-checked)
  expect_equal(res$ci_lower,       0.2102, tolerance = 5e-4)
  expect_equal(res$ci_upper,       0.3898, tolerance = 5e-4)
  expect_equal(res$margin_of_error, 0.08982, tolerance = 5e-4)
})

test_that("Case 2: supplied ICC widens CI via design effect", {
  res <- estimate_prevalence(
    x   = rep(3, 10),
    n   = rep(10, 10),
    icc = 0.05
  )

  expect_equal(res$prevalence, 0.3,  tolerance = 1e-6)
  expect_equal(res$deff,       1.45, tolerance = 1e-6)
  expect_equal(res$n_eff,      100 / 1.45, tolerance = 1e-4)

  # CI should be wider than the no-clustering case
  res_plain <- estimate_prevalence(x = 30, n = 100, icc = 0)
  expect_gt(res$ci_upper - res$ci_lower,
            res_plain$ci_upper - res_plain$ci_lower)

  expect_equal(res$ci_lower, 0.1918, tolerance = 5e-4)
  expect_equal(res$ci_upper, 0.4082, tolerance = 5e-4)
})

test_that("Case 3: imperfect test applies Rogan-Gladen to estimate and CI", {
  res <- estimate_prevalence(x = 30, n = 100,
                             sensitivity = 0.9, specificity = 0.95)

  # p_true = (0.3 - 0.05) / 0.85
  expect_equal(res$prevalence, 0.25 / 0.85, tolerance = 1e-6)

  # MOE inflated by 1/correction
  res_perfect <- estimate_prevalence(x = 30, n = 100)
  expect_equal(res$margin_of_error,
               res_perfect$margin_of_error / 0.85,
               tolerance = 1e-6)

  expect_equal(res$ci_lower, 0.1885, tolerance = 5e-4)
  expect_equal(res$ci_upper, 0.3998, tolerance = 5e-4)
})

test_that("ICC estimated from data is non-negative and finite", {
  res <- estimate_prevalence(
    x = c(0, 4, 0, 22, 25, 16, 12, 8),
    n = c(60, 80, 70, 100, 40, 60, 50, 90)
  )
  expect_true(is.finite(res$icc_used))
  expect_gte(res$icc_used, 0)
  expect_lte(res$icc_used, 1)
  expect_gte(res$deff, 1)
  expect_equal(res$n_total, 550)
})

test_that("icc = 0 forces no clustering adjustment", {
  res_null <- estimate_prevalence(x = 30, n = 100, icc = 0)
  res_no   <- estimate_prevalence(x = 30, n = 100, icc = NULL)  # single site

  expect_equal(res_null$deff, 1, tolerance = 1e-6)
  # Both single-site paths should agree
  expect_equal(res_null$ci_lower, res_no$ci_lower, tolerance = 1e-10)
})

test_that("perfect test is identity: se=sp=1 changes nothing", {
  res_default <- estimate_prevalence(x = 30, n = 100)
  res_explicit <- estimate_prevalence(x = 30, n = 100,
                                      sensitivity = 1, specificity = 1)
  expect_equal(res_default$prevalence, res_explicit$prevalence)
  expect_equal(res_default$ci_lower,   res_explicit$ci_lower)
  expect_equal(res_default$ci_upper,   res_explicit$ci_upper)
})

test_that("input validation catches bad arguments", {
  expect_error(estimate_prevalence(x = c(3, 5), n = c(10)),
               "length")
  expect_error(estimate_prevalence(x = 60, n = 50),
               "cannot exceed")
  expect_error(estimate_prevalence(x = 10, n = 50,
                                   sensitivity = 0.3, specificity = 0.3),
               "must exceed 1")
  # negative counts
  expect_error(estimate_prevalence(x = -1, n = 100),
               "non-negative")
  # conf_level out of range
  expect_error(estimate_prevalence(x = 30, n = 100, conf_level = 0), "`conf_level`")
  expect_error(estimate_prevalence(x = 30, n = 100, conf_level = 1), "`conf_level`")
  # icc out of range when supplied
  expect_error(estimate_prevalence(x = 30, n = 100, icc =  2),  "`icc`")
  expect_error(estimate_prevalence(x = 30, n = 100, icc = -0.1), "`icc`")
  # fpc_N must be positive
  expect_error(estimate_prevalence(x = 30, n = 100, fpc_N =   0), "`fpc_N`")
  expect_error(estimate_prevalence(x = 30, n = 100, fpc_N = -10), "`fpc_N`")
  # fpc_N smaller than sample
  expect_error(estimate_prevalence(x = c(3,3), n = c(10,10), fpc_N = 5),
               "at least as large")
})

test_that("round-trip: estimate_prevalence recovers design_precision MOE (SRS)", {
  # design_precision(0.3, 0.05) -> n = 323
  # observe exactly 30% -> x = round(0.3 * 323) = 97
  res <- estimate_prevalence(x = 97, n = 323)
  expect_equal(res$margin_of_error, 0.05, tolerance = 5e-3)
})

test_that("round-trip: estimate_prevalence recovers design_precision MOE (clustered)", {
  # design_precision(0.3, 0.05, n_per_site=10, icc=0.05) -> 47 sites of 10
  # observe 30% per site -> x = 3 per site
  res <- estimate_prevalence(x = rep(3, 47), n = rep(10, 47), icc = 0.05)
  expect_equal(res$margin_of_error, 0.05, tolerance = 5e-3)
})
