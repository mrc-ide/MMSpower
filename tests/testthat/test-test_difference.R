# Reference cases for test_difference()
#
# Score (default): with no clustering, the pooled z-test is
# prop.test(correct = FALSE). 8/100 vs 20/120: p = 0.05479, score interval
# [-0.174678, 0.001913] (unchanged from the previous implementation).
# Wald, 3/80 vs 12/90: d = -0.095833,
#   se = sqrt(0.0375 * 0.9625 / 80 + 0.13333 * 0.86667 / 90) = 0.041656,
#   z = -2.3006, p = 0.0214, interval d +/- 1.959964 * se = [-0.17748, -0.01419]

test_that("score test matches prop.test with no clustering", {
  res <- test_difference(x1 = 8, n1 = 100, x2 = 20, n2 = 120)
  ref <- prop.test(c(8, 20), c(100, 120), correct = FALSE)
  expect_equal(res$p_value, ref$p.value)
  expect_equal(res$difference, 0.08 - 20 / 120)
  expect_equal(res$ci_lower, -0.174678, tolerance = 1e-4)
  expect_equal(res$ci_upper, 0.001913, tolerance = 1e-3)
  expect_false(res$reject)
})

test_that("Wald test and interval match the hand-check", {
  res <- test_difference(x1 = 3, n1 = 80, x2 = 12, n2 = 90, method = "wald")
  expect_equal(res$p_value, 0.0214, tolerance = 1e-2)
  expect_equal(res$ci_lower, -0.17748, tolerance = 1e-4)
  expect_equal(res$ci_upper, -0.01419, tolerance = 1e-3)
  expect_true(res$reject)
})

test_that("one-sided tests halve the p-value and give one-sided intervals", {
  two  <- test_difference(8, 100, 20, 120)
  less <- test_difference(8, 100, 20, 120, alternative = "less")
  expect_equal(less$p_value, two$p_value / 2)
  expect_equal(less$ci_lower, -1)
  expect_lt(less$ci_upper, 0)
  expect_true(less$reject)
})

test_that("clustering makes the same counts less convincing", {
  x1 <- c(1, 9, 3); n1 <- c(50, 60, 40)
  x2 <- c(4, 18, 9); n2 <- c(55, 70, 45)
  plain <- test_difference(x1, n1, x2, n2, alternative = "less", icc = 0)
  est   <- test_difference(x1, n1, x2, n2, alternative = "less")
  expect_equal(c(est$deff1, est$deff2), c(2.66, 3.36), tolerance = 0.01)
  expect_gt(est$p_value, 10 * plain$p_value)
  two_icc <- test_difference(x1, n1, x2, n2, icc = c(0.02, 0.08))
  expect_equal(two_icc$deff1, 1 + (sum(n1^2) / sum(n1) - 1) * 0.02)
  expect_equal(two_icc$deff2, 1 + (sum(n2^2) / sum(n2) - 1) * 0.08)
})

test_that("an imperfect test widens the difference but not the p-value", {
  perfect <- test_difference(8, 100, 20, 120)
  imperf  <- test_difference(8, 100, 20, 120, sensitivity = 0.9, specificity = 0.97)
  expect_equal(imperf$p_value, perfect$p_value)
  expect_equal(imperf$difference, perfect$difference / 0.87)
  expect_equal(imperf$ci_lower, perfect$ci_lower / 0.87)
})

test_that("reject is TRUE exactly when the interval leaves out 0", {
  for (method in c("score", "wald")) {
    for (alt in c("two.sided", "greater", "less")) {
      for (x2 in c(2, 6, 10, 14, 20)) {
        res <- test_difference(c(5, 6), c(60, 60), c(x2, x2 + 1), c(60, 60),
                               alternative = alt, icc = 0.02, method = method)
        expect_identical(res$reject, res$ci_lower > 0 || res$ci_upper < 0,
                         info = paste(method, alt, x2))
      }
    }
  }
})

test_that("invalid inputs are rejected", {
  expect_error(test_difference(c(1, 2), 10, 3, 10), "same length")
  expect_error(test_difference(1, 10, 11, 10), "x2")
  expect_error(test_difference(1, 10, 3, 10, icc = c(0.1, 0.1, 0.1)), "icc")
  expect_error(test_difference(1, 10, 3, 10, method = "newcombe"))
})
