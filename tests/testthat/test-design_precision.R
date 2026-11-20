# Hand-checked reference cases for design_precision()
#
# Wald: n0 = z^2 * p_app * (1 - p_app) / (moe^2 * (se + sp - 1)^2),
#   p_app = p*se + (1-p)*(1-sp), z = qnorm(0.975) ~= 1.959964
#
# p = 0.3, moe = 0.05, perfect test: n0 = 322.68 -> 323
# se = 0.9, sp = 0.95: p_app = 0.305, n0 = 450.82 -> 451
# n_per_site = 10, icc = 0.05: deff = 1.45, n = 467.89 -> 47 sites x 10 = 470
# n_sites = 50, icc = 0.05: n = 322.68 * 50 * 0.95 / (50 - 322.68 * 0.05)
#   = 452.6 -> 10 per site x 50 = 500
# n_sites = 10, icc = 0.05: 10 < 322.68 * 0.05, unachievable;
#   min moe = z * sqrt(0.21 * 0.05 / 10) = 6.35%

# the hand-checked cases below use the Wald closed form
dp_wald <- function(...) design_precision(..., ci_method = "wald")

test_that("unclustered, perfect test", {
  res <- dp_wald(0.3, 0.05)
  expect_equal(res$n, 323)
  expect_equal(res$n_eff, 323)
  expect_equal(res$deff, 1)
  expect_null(res$n_sites)
  expect_null(res$n_per_site)
})

test_that("imperfect test increases n", {
  res <- dp_wald(0.3, 0.05, sensitivity = 0.9, specificity = 0.95)
  expect_equal(res$n, 451)
})

test_that("n is largest at prevalence 0.5 and grows with conf_level", {
  expect_equal(dp_wald(0.5, 0.05)$n, 385)
  expect_gt(dp_wald(0.3, 0.05, conf_level = 0.99)$n,
            dp_wald(0.3, 0.05)$n)
})

test_that("fixed site size solves for the number of sites", {
  res <- dp_wald(0.3, 0.05, n_per_site = 10, icc = 0.05)
  expect_equal(res$n_sites, 47)
  expect_equal(res$n_per_site, 10)
  expect_equal(res$n, 470)
  expect_equal(res$deff, 1.45)
  expect_equal(res$n_eff, 323)
})

test_that("fixed number of sites solves for people per site", {
  res <- dp_wald(0.3, 0.05, n_sites = 50, icc = 0.05)
  expect_equal(res$n_sites, 50)
  expect_equal(res$n_per_site, 10)
  expect_equal(res$n, 500)
  expect_equal(res$deff, 1.45)
})

test_that("icc = 0 with sites just spreads the unclustered n", {
  res <- dp_wald(0.3, 0.05, n_sites = 20)
  expect_equal(res$n_per_site, 17)
  expect_equal(res$n, 340)
  expect_equal(res$deff, 1)
})

test_that("too few sites for the target MOE stops, giving the minimum MOE", {
  expect_error(dp_wald(0.3, 0.05, n_sites = 10, icc = 0.05), "6.4%")
})

test_that("invalid inputs are rejected", {
  expect_error(design_precision(0, 0.05), "prevalence")
  expect_error(design_precision(0.3, 5), "moe")
  expect_error(design_precision(0.3, 0.05, sensitivity = 0.5, specificity = 0.5),
               "greater than 1")
  expect_error(design_precision(0.3, 0.05, icc = 0.05), "needs")
  expect_error(design_precision(0.3, 0.05, n_sites = 10, n_per_site = 10), "not both")
  expect_error(design_precision(0.3, 0.05, n_sites = 2.5, icc = 0.05), "whole number")
  expect_error(design_precision(0.3), "Give `moe`, or a full design")
  expect_error(design_precision(0.3, n = 100, n_sites = 5), "without sites")
  expect_error(design_precision(0.3, n_sites = 2, n_per_site = c(10, 20)), "not both")
  expect_error(design_precision(0.3, n_per_site = c(10, 0)), "at least 1")
})

test_that("a full design returns its margin of error", {
  # DRpower::get_margin(rep(40, 6), 0.2, 0.05): 11.31 to 28.69 percent
  res <- dp_wald(0.2, n_sites = 6, n_per_site = 40, icc = 0.05)
  expect_equal(c(res$n, res$n_sites, res$n_per_site), c(240, 6, 40))
  expect_equal(res$deff, 1 + 39 * 0.05)
  expect_equal(res$n_eff, 240 / res$deff)
  expect_equal(res$moe, 0.0869, tolerance = 1e-3)
  # no sites: 323 people give the 5 percent the design was for
  expect_equal(dp_wald(0.3, n = 323)$moe, 0.05, tolerance = 1e-3)
})

test_that("sites of different sizes use the size-weighted mean site size", {
  s   <- c(250, 250, 60, 60, 40, 40)
  res <- dp_wald(0.2, n_per_site = s, icc = 0.05)
  n_w <- sum(s^2) / sum(s)
  expect_equal(res$deff, 1 + (n_w - 1) * 0.05)
  expect_equal(res$moe, qnorm(0.975) * sqrt(0.2 * 0.8 * res$deff / 700))
  expect_equal(res$n_per_site, s)
  # equal sizes given as a vector are the same as n_sites and n_per_site
  expect_equal(design_precision(0.2, n_per_site = rep(40, 6), icc = 0.05),
               design_precision(0.2, n_sites = 6, n_per_site = 40, icc = 0.05))
})

test_that("a solved design reports its own margin, within the target", {
  for (args in list(list(), list(n_per_site = 10, icc = 0.05),
                    list(n_sites = 50, icc = 0.05))) {
    res  <- do.call(design_precision, c(list(0.3, 0.05), args))
    back <- if (is.null(res$n_sites)) design_precision(0.3, n = res$n)
            else design_precision(0.3, n_sites = res$n_sites,
                                  n_per_site = res$n_per_site, icc = args$icc)
    expect_lte(res$moe, 0.05)
    expect_equal(back$moe, res$moe)
  }
})

test_that("Clopper-Pearson is the default, as in estimate_prevalence()", {
  expect_identical(formals(design_precision)$ci_method,
                   formals(estimate_prevalence)$ci_method)
  # moe of a full design: half-width of the Clopper-Pearson interval at the
  # expected prevalence, on n / deff people
  res <- design_precision(0.2, n_sites = 6, n_per_site = 40, icc = 0.05)
  m   <- 240 / res$deff
  ci  <- c(qbeta(0.025, 0.2 * m, 0.8 * m + 1), qbeta(0.975, 0.2 * m + 1, 0.8 * m))
  expect_equal(res$moe, diff(ci) / 2)
})

test_that("each method sizes for the interval estimate_prevalence() gives", {
  for (method in c("clopper-pearson", "wald", "agresti-coull")) {
    for (p in c(0.05, 0.3)) {
      res <- design_precision(p, 0.03, ci_method = method)
      expect_lte(res$moe, 0.03)
      # one fewer person misses the target
      expect_gt(design_precision(p, n = res$n - 1, ci_method = method)$moe, 0.03)
      # analysing the expected counts gives back the margin (exactly, when
      # p * n is a whole number)
      n <- 1000
      expect_equal(design_precision(p, n = n, ci_method = method)$moe,
                   estimate_prevalence(p * n, n, ci_method = method)$moe,
                   info = paste(method, p))
    }
  }
  # Clopper-Pearson is wider than Wald, so needs more people
  expect_gt(design_precision(0.3, 0.05)$n, dp_wald(0.3, 0.05)$n)
})

test_that("an imperfect test is allowed for in every method", {
  for (method in c("clopper-pearson", "agresti-coull")) {
    perfect <- design_precision(0.3, 0.05, ci_method = method)$n
    expect_gt(design_precision(0.3, 0.05, sensitivity = 0.9, specificity = 0.95,
                               ci_method = method)$n, perfect)
  }
})

test_that("too few sites reports the minimum margin for the chosen method", {
  expect_error(design_precision(0.3, 0.05, n_sites = 10, icc = 0.05), "6.6%")
  expect_error(design_precision(0.3, 0.05, ci_method = "exact"))
})
