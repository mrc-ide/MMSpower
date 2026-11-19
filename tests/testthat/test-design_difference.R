# Reference cases for design_difference()
#
# No clustering, 10% vs 20%, 80% power, two-sided: power.prop.test gives
# 198.96 -> 199 per group.
# n_per_site = 30, icc = 0.02: deff = 1.58, sites x with 30x / 1.58 >= 198.96
#   -> x = 10.48 -> 11 sites per group.
# n_sites = 30, icc = 0.02: equal groups, so the closed form applies:
#   n = 198.96 * 30 * 0.98 / (30 - 198.96 * 0.02) = 224.8 -> 8 per site.
# n_sites = 8, icc = 0.05: at most 160 independent people per group, which
#   gives 70.9% power.

test_that("no clustering matches power.prop.test", {
  res <- design_difference(0.10, 0.20)
  ref <- ceiling(power.prop.test(p1 = 0.10, p2 = 0.20, power = 0.8)$n)
  expect_equal(res$n_per_group, c(ref, ref))
  expect_equal(res$n, 2 * ref)
  expect_null(res$n_sites)
  for (pr in list(c(0.05, 0.15), c(0.3, 0.45), c(0.02, 0.08))) {
    ref <- ceiling(power.prop.test(p1 = pr[1], p2 = pr[2], power = 0.9)$n)
    expect_equal(design_difference(pr[1], pr[2], power = 0.9)$n_per_group[1], ref)
  }
})

test_that("one-sided, imperfect test and higher power change n as expected", {
  base <- design_difference(0.10, 0.20)$n
  expect_lt(design_difference(0.10, 0.20, alternative = "less")$n, base)
  expect_gt(design_difference(0.10, 0.20, sensitivity = 0.9, specificity = 0.95)$n, base)
  expect_gt(design_difference(0.10, 0.20, power = 0.9)$n, base)
})

test_that("fixed site size solves for sites per group", {
  res <- design_difference(0.10, 0.20, n_per_site = 30, icc = 0.02)
  expect_equal(res$n_sites, c(11, 11))
  expect_equal(res$n_per_group, c(330, 330))
  expect_equal(res$deff, c(1.58, 1.58))
  expect_equal(res$n_eff, c(199, 199))
})

test_that("fixed number of sites solves for people per site", {
  res <- design_difference(0.10, 0.20, n_sites = 30, icc = 0.02)
  expect_equal(res$n_per_site, c(8, 8))
  expect_equal(res$n_per_group, c(240, 240))
})

test_that("groups can differ in sites and ICC", {
  res <- design_difference(0.10, 0.20, n_sites = c(20, 40), icc = c(0.02, 0.05))
  expect_equal(res$n_sites, c(20, 40))
  expect_equal(res$n_per_site, c(10, 10))
  expect_equal(res$n_per_group, c(200, 400))
  expect_equal(res$deff, 1 + 9 * c(0.02, 0.05))
  # one fewer person per site would miss the target
  m <- c(20, 40) * 9 / (1 + 8 * c(0.02, 0.05))
  z_a <- qnorm(0.975); z_b <- qnorm(0.8)
  expect_gt(.diff_gap(m, c(0.1, 0.2), z_a, z_b), 0)
})

test_that("more sites than needed: 1 person at each of just enough sites", {
  res <- design_difference(0.10, 0.20, n_sites = 500, icc = 0.05)
  expect_equal(res$n_per_site, c(1, 1))
  expect_equal(res$n_sites, c(199, 199))
  expect_equal(res$n_per_group, c(199, 199))
})

test_that("too few sites stops, giving the highest power they can reach", {
  expect_error(design_difference(0.10, 0.20, n_sites = 8, icc = 0.05), "70.9%")
})

test_that("invalid inputs are rejected", {
  expect_error(design_difference(0.2, 0.2), "differ")
  expect_error(design_difference(0.1, 0.2, alternative = "greater"), "above")
  expect_error(design_difference(0.1, 0.2, power = 0.4), "power")
  expect_error(design_difference(0.1, 0.2, conf_level = 0.3), "conf_level")
  expect_error(design_difference(0.1, 0.2, icc = 0.05), "needs")
  expect_error(design_difference(0.1, 0.2, n_sites = c(10, 20, 30)), "one value")
  expect_error(design_difference(0.1, 0.2, n_sites = 10, n_per_site = 10), "at most one")
})
