# Reference cases for design_threshold()
#
# Unclustered answers are the exact binomial search (the same as the
# previous, longer implementation): 5% vs 10% greater -> 191; 10% vs 5%
# less -> 203 (179 is the first n to reach 80%, but 182-190 fall below).

# Brute-force exact power for one n: loop over every count.
brute_power <- function(n, theta, p, alpha, alternative) {
  ks <- 0:n
  if (alternative == "greater") {
    k <- min(c(ks[stats::pbinom(ks - 1, n, theta, lower.tail = FALSE) < alpha], n + 1))
    stats::pbinom(k - 1, n, p, lower.tail = FALSE)
  } else {
    k <- max(c(-1, ks[stats::pbinom(ks, n, theta) < alpha]))
    stats::pbinom(k, n, p)
  }
}

test_that("exact power matches a brute-force count over every outcome", {
  n <- c(1:40, 191, 777)
  expect_equal(.threshold_power(n, 0.05, 0.10, "greater", 0.05),
               sapply(n, brute_power, 0.05, 0.10, 0.05, "greater"))
  expect_equal(.threshold_power(n, 0.10, 0.05, "less", 0.05),
               sapply(n, brute_power, 0.10, 0.05, 0.05, "less"))
})

test_that("unclustered sample sizes match the reference cases", {
  expect_equal(design_threshold(0.05, 0.10)$n, 191)
  expect_equal(design_threshold(0.10, 0.05, alternative = "less")$n, 203)
})

test_that("n is the smallest n from which power never drops below target", {
  pw <- .threshold_power(1:400, 0.10, 0.05, "less", 0.05)
  expect_true(all(pw[203:400] >= 0.8))
  expect_lt(pw[202], 0.8)
  expect_gte(pw[179], 0.8)   # an earlier n reaches it, but is not returned
})

test_that("the design gives a significant test_threshold result at its cut-off", {
  # 191 people, threshold 5%: 16 positives is the smallest count that rejects
  expect_true(test_threshold(16, 191, threshold = 0.05)$reject)
  expect_false(test_threshold(15, 191, threshold = 0.05)$reject)
})

test_that("more power, a two-sided test or an imperfect test need more people", {
  base <- design_threshold(0.05, 0.10)$n
  expect_gt(design_threshold(0.05, 0.10, power = 0.9)$n, base)
  expect_gt(design_threshold(0.05, 0.10, alternative = "two.sided")$n, base)
  expect_gt(design_threshold(0.05, 0.10, sensitivity = 0.9, specificity = 0.95)$n, base)
})

test_that("fixed site size scales n_eff by the design effect", {
  res <- design_threshold(0.05, 0.10, n_per_site = 20, icc = 0.05)
  # 191 * 1.95 = 372.45 -> 19 sites of 20
  expect_equal(res$n_eff, 191)
  expect_equal(res$deff, 1.95)
  expect_equal(res$n_sites, 19)
  expect_equal(res$n, 380)
})

test_that("fixed number of sites uses the closed-form solve", {
  res <- design_threshold(0.05, 0.10, n_sites = 20, icc = 0.02)
  # 191 * 20 * 0.98 / (20 - 191 * 0.02) = 231.4 -> 12 per site
  expect_equal(res$n_per_site, 12)
  expect_equal(res$n, 240)
  expect_error(design_threshold(0.05, 0.10, n_sites = 5, icc = 0.05), "unreachable")
})

test_that("invalid inputs are rejected", {
  expect_error(design_threshold(0.05, 0.05), "differ")
  expect_error(design_threshold(0.10, 0.05), "above")
  expect_error(design_threshold(0.05, 0.10, alternative = "less"), "below")
  expect_error(design_threshold(0.05, 0.051), "10,000")
  expect_error(design_threshold(0.05, 0.10, icc = 0.05), "needs")
})
