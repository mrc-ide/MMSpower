# Hand-checked reference cases for estimate_prevalence()
#
# p_hat = sum(x) / sum(n), n_eff = sum(n) / deff, z = qnorm(0.975)
# Wald: p_hat +/- z * sqrt(p_hat * (1 - p_hat) / n_eff), then Rogan-Gladen
#   (p - (1 - sp)) / (se + sp - 1) on the estimate and both ends.
#
# x = 30, n = 100: CI = 0.3 +/- 0.089817 = [0.210183, 0.389817]
# se = 0.9, sp = 0.95: prevalence = 0.25 / 0.85 = 0.294118,
#   CI = [0.160183, 0.339817] / 0.85 = [0.188451, 0.399785]
# 10 sites of 10, 3 positive each, icc = 0.05: deff = 1.45,
#   CI = 0.3 +/- 0.108154 = [0.191846, 0.408154]

test_that("one site, perfect test gives the hand-checked Wald interval", {
  res <- estimate_prevalence(x = 30, n = 100)
  expect_equal(res$prevalence, 0.3)
  expect_equal(res$ci_lower, 0.210183, tolerance = 1e-4)
  expect_equal(res$ci_upper, 0.389817, tolerance = 1e-4)
  expect_equal(res$moe, 0.089817, tolerance = 1e-4)
  expect_equal(res$deff, 1)
})

test_that("an imperfect test corrects the estimate and both ends", {
  res <- estimate_prevalence(x = 30, n = 100, sensitivity = 0.9, specificity = 0.95)
  expect_equal(res$prevalence, 0.294118, tolerance = 1e-4)
  expect_equal(res$ci_lower, 0.188451, tolerance = 1e-4)
  expect_equal(res$ci_upper, 0.399785, tolerance = 1e-4)
})

test_that("a supplied icc widens the interval through the design effect", {
  res <- estimate_prevalence(x = rep(3, 10), n = rep(10, 10), icc = 0.05)
  expect_equal(res$deff, 1.45)
  expect_equal(res$n_eff, 100 / 1.45)
  expect_equal(res$ci_lower, 0.191846, tolerance = 1e-4)
  expect_equal(res$ci_upper, 0.408154, tolerance = 1e-4)
})

test_that("icc estimated from the data (8-site example)", {
  # var of site prevalences 0.042917; expected under SRS 0.002101;
  # ratio 20.43 -> icc = 19.43 / 67.75 = 0.2867; size-weighted site size 74,
  # deff = 1 + 73 * 0.2867 = 21.93
  res <- estimate_prevalence(
    x = c(0, 4, 0, 22, 25, 16, 12, 8),
    n = c(60, 80, 70, 100, 40, 60, 50, 90)
  )
  expect_equal(res$prevalence, 0.158182, tolerance = 1e-4)
  expect_equal(res$icc, 0.286728, tolerance = 1e-4)
  expect_equal(res$deff, 21.931152, tolerance = 1e-4)
  expect_equal(res$ci_lower, 0.015363, tolerance = 1e-4)
  expect_equal(res$ci_upper, 0.301000, tolerance = 1e-4)
})

test_that("identical sites, or a single site, give an estimated icc of 0", {
  expect_equal(estimate_prevalence(rep(3, 10), rep(10, 10))$icc, 0)
  expect_equal(estimate_prevalence(30, 100)$deff, 1)
})

test_that("clopper-pearson matches binom.test", {
  res <- estimate_prevalence(x = 30, n = 100, ci_method = "clopper-pearson")
  expect_equal(c(res$ci_lower, res$ci_upper), binom.test(30, 100)$conf.int[1:2])
})

test_that("agresti-coull matches the hand-check", {
  # n~ = 103.84, p~ = 0.30740, CI = 0.30740 +/- 0.088747
  res <- estimate_prevalence(x = 30, n = 100, ci_method = "agresti-coull")
  expect_equal(res$ci_lower, 0.218651, tolerance = 1e-4)
  expect_equal(res$ci_upper, 0.396146, tolerance = 1e-4)
})

test_that("no positives: the interval is kept within [0, 1]", {
  wald <- estimate_prevalence(x = 0, n = 50)
  expect_equal(c(wald$ci_lower, wald$ci_upper), c(0, 0))
  cp <- estimate_prevalence(x = 0, n = 50, ci_method = "clopper-pearson")
  expect_equal(cp$ci_lower, 0)
  expect_gt(cp$ci_upper, 0)
})

test_that("data that do not fit the test accuracy warn", {
  # 1 positive in 200 is far fewer than the 5% false positives expected
  expect_warning(estimate_prevalence(x = 1, n = 200, specificity = 0.95),
                 "do not fit")
})

test_that("the n from design_precision gives back the target moe", {
  # design_precision(0.3, 0.05) -> 323
  expect_equal(estimate_prevalence(x = 97, n = 323)$moe, 0.05, tolerance = 0.01)
  # design_precision(0.3, 0.05, n_per_site = 10, icc = 0.05) -> 47 sites of 10
  res <- estimate_prevalence(x = rep(3, 47), n = rep(10, 47), icc = 0.05)
  expect_equal(res$moe, 0.05, tolerance = 0.01)
})

test_that("invalid inputs are rejected", {
  expect_error(estimate_prevalence(x = c(1, 2), n = 10), "same length")
  expect_error(estimate_prevalence(x = 11, n = 10), "between 0")
  expect_error(estimate_prevalence(x = 1.5, n = 10), "whole numbers")
  expect_error(estimate_prevalence(x = 1, n = 10, icc = 2), "icc")
  expect_error(estimate_prevalence(x = 1, n = 10, ci_method = "exact"))
})
