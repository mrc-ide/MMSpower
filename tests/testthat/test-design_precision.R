# Hand-checked reference cases for design_precision()
#
# n0 = z^2 * p_app * (1 - p_app) / (moe^2 * (se + sp - 1)^2),
#   p_app = p*se + (1-p)*(1-sp), z = qnorm(0.975) ~= 1.959964
#
# p = 0.3, moe = 0.05, perfect test: n0 = 322.68 -> 323
# se = 0.9, sp = 0.95: p_app = 0.305, n0 = 450.82 -> 451
# n_per_site = 10, icc = 0.05: deff = 1.45, n = 467.89 -> 47 sites x 10 = 470
# n_sites = 50, icc = 0.05: n = 322.68 * 50 * 0.95 / (50 - 322.68 * 0.05)
#   = 452.6 -> 10 per site x 50 = 500
# n_sites = 10, icc = 0.05: 10 < 322.68 * 0.05, unachievable;
#   min moe = z * sqrt(0.21 * 0.05 / 10) = 6.35%

test_that("unclustered, perfect test", {
  res <- design_precision(0.3, 0.05)
  expect_equal(res$n, 323)
  expect_equal(res$n_eff, 323)
  expect_equal(res$deff, 1)
  expect_null(res$n_sites)
  expect_null(res$n_per_site)
})

test_that("imperfect test increases n", {
  res <- design_precision(0.3, 0.05, sensitivity = 0.9, specificity = 0.95)
  expect_equal(res$n, 451)
})

test_that("n is largest at prevalence 0.5 and grows with conf_level", {
  expect_equal(design_precision(0.5, 0.05)$n, 385)
  expect_gt(design_precision(0.3, 0.05, conf_level = 0.99)$n,
            design_precision(0.3, 0.05)$n)
})

test_that("fixed site size solves for the number of sites", {
  res <- design_precision(0.3, 0.05, n_per_site = 10, icc = 0.05)
  expect_equal(res$n_sites, 47)
  expect_equal(res$n_per_site, 10)
  expect_equal(res$n, 470)
  expect_equal(res$deff, 1.45)
  expect_equal(res$n_eff, 323)
})

test_that("fixed number of sites solves for people per site", {
  res <- design_precision(0.3, 0.05, n_sites = 50, icc = 0.05)
  expect_equal(res$n_sites, 50)
  expect_equal(res$n_per_site, 10)
  expect_equal(res$n, 500)
  expect_equal(res$deff, 1.45)
})

test_that("icc = 0 with sites just spreads the unclustered n", {
  res <- design_precision(0.3, 0.05, n_sites = 20)
  expect_equal(res$n_per_site, 17)
  expect_equal(res$n, 340)
  expect_equal(res$deff, 1)
})

test_that("too few sites for the target MOE stops, giving the minimum MOE", {
  expect_error(design_precision(0.3, 0.05, n_sites = 10, icc = 0.05), "6.4%")
})

test_that("invalid inputs are rejected", {
  expect_error(design_precision(0, 0.05), "prevalence")
  expect_error(design_precision(0.3, 5), "moe")
  expect_error(design_precision(0.3, 0.05, sensitivity = 0.5, specificity = 0.5),
               "greater than 1")
  expect_error(design_precision(0.3, 0.05, icc = 0.05), "needs")
  expect_error(design_precision(0.3, 0.05, n_sites = 10, n_per_site = 10), "at most one")
  expect_error(design_precision(0.3, 0.05, n_sites = 2.5, icc = 0.05), "whole number")
})
