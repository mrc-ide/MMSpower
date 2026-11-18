# Reference cases for test_threshold()
#
# With one site and deff = 1 the one-sided p-values and the
# Clopper-Pearson interval are those of binom.test().
# Imperfect test, threshold 0.05, se 0.9, sp 0.95: measured threshold
# 0.0925; 12/100 gives p = 0.21, 15/100 gives p = 0.04.

test_that("one site, perfect test matches binom.test (greater and less)", {
  res <- test_threshold(x = 8, n = 100, threshold = 0.05)
  ref <- binom.test(8, 100, p = 0.05, alternative = "greater")
  expect_equal(res$p_value, ref$p.value)
  expect_equal(c(res$ci_lower, res$ci_upper), ref$conf.int[1:2])
  expect_equal(res$prevalence, 0.08)

  res <- test_threshold(x = 2, n = 100, threshold = 0.05, alternative = "less")
  ref <- binom.test(2, 100, p = 0.05, alternative = "less")
  expect_equal(res$p_value, ref$p.value)
  expect_equal(c(res$ci_lower, res$ci_upper), ref$conf.int[1:2])
})

test_that("two-sided p-value is twice the smaller one-sided p-value", {
  g <- test_threshold(x = 12, n = 100, threshold = 0.10)$p_value
  l <- test_threshold(x = 12, n = 100, threshold = 0.10, alternative = "less")$p_value
  two <- test_threshold(x = 12, n = 100, threshold = 0.10, alternative = "two.sided")
  expect_equal(two$p_value, min(1, 2 * min(g, l)))
  expect_equal(c(two$ci_lower, two$ci_upper), binom.test(12, 100)$conf.int[1:2])
})

test_that("an imperfect test raises the threshold on the measured scale", {
  not_enough <- test_threshold(x = 12, n = 100, threshold = 0.05,
                               sensitivity = 0.9, specificity = 0.95)
  enough <- test_threshold(x = 15, n = 100, threshold = 0.05,
                           sensitivity = 0.9, specificity = 0.95)
  expect_equal(not_enough$p_value, 0.21, tolerance = 0.05)
  expect_false(not_enough$reject)
  expect_equal(enough$p_value, 0.04, tolerance = 0.1)
  expect_true(enough$reject)
  expect_equal(enough$prevalence, (0.15 - 0.05) / 0.85)
})

test_that("clustering makes the same counts less convincing", {
  x <- rep(4, 5); n <- rep(50, 5)
  plain <- test_threshold(x, n, threshold = 0.05, icc = 0)
  clus  <- test_threshold(x, n, threshold = 0.05, icc = 0.05)
  expect_equal(clus$deff, 1 + 49 * 0.05)
  expect_equal(clus$n_eff, 250 / clus$deff)
  expect_gt(clus$p_value, plain$p_value)
})

test_that("Wald test and interval match the hand-check", {
  # se = sqrt(0.08 * 0.92 / 100) = 0.027129, z = 0.03 / 0.027129 = 1.10583,
  # p = pnorm(-1.10583) = 0.13441, L = 0.08 - 1.644854 * 0.027129 = 0.035376
  res <- test_threshold(x = 8, n = 100, threshold = 0.05, method = "wald")
  expect_equal(res$p_value, 0.13441, tolerance = 1e-4)
  expect_equal(res$ci_lower, 0.035376, tolerance = 1e-4)
  expect_equal(res$ci_upper, 1)
  expect_false(res$reject)
})

test_that("reject is TRUE exactly when the threshold is outside the interval", {
  for (method in c("exact", "wald")) {
    for (alt in c("greater", "less", "two.sided")) {
      for (x in list(c(0, 0), c(2, 3), c(5, 6), c(9, 10), c(15, 16))) {
        res <- test_threshold(x = x, n = c(60, 60), threshold = 0.07,
                              alternative = alt, icc = 0.02, method = method)
        outside <- 0.07 < res$ci_lower || 0.07 > res$ci_upper
        expect_identical(res$reject, outside, info = paste(method, alt, x[1]))
      }
    }
  }
})

test_that("invalid inputs are rejected", {
  expect_error(test_threshold(x = 5, n = 100, threshold = 0), "threshold")
  expect_error(test_threshold(x = 5, n = 100, threshold = 0.05,
                              alternative = "up"))
  expect_error(test_threshold(x = 5, n = 100, threshold = 0.05, icc = -1), "icc")
  expect_error(test_threshold(x = 101, n = 100, threshold = 0.05), "between 0")
  expect_error(test_threshold(x = 5, n = 100, threshold = 0.05, method = "score"))
})
