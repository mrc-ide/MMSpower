# Reference cases for design_difference()
#
# No clustering, 10% vs 20%, 80% power, two-sided: power.prop.test gives
# 198.96 -> 199 per group.
# n_per_site1 = 30, icc = 0.02: deff = 1.58, sites x with 30x / 1.58 >= 198.96
#   -> x = 10.48 -> 11 sites per group.
# n_sites1 = 30, icc = 0.02: equal groups, so the closed form applies:
#   n = 198.96 * 30 * 0.98 / (30 - 198.96 * 0.02) = 224.8 -> 8 per site.
# n_sites1 = 8, icc = 0.05: at most 160 independent people per group, which
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
  res <- design_difference(0.10, 0.20, n_per_site1 = 30, icc = 0.02)
  expect_equal(res$n_sites, c(11, 11))
  expect_equal(res$n_per_group, c(330, 330))
  expect_equal(res$deff, c(1.58, 1.58))
  expect_equal(res$n_eff, c(199, 199))
})

test_that("fixed number of sites solves for people per site", {
  res <- design_difference(0.10, 0.20, n_sites1 = 30, icc = 0.02)
  expect_equal(res$n_per_site, c(8, 8))
  expect_equal(res$n_per_group, c(240, 240))
})

test_that("groups can differ in sites and ICC", {
  res <- design_difference(0.10, 0.20, n_sites1 = 20, n_sites2 = 40, icc = c(0.02, 0.05))
  expect_equal(res$n_sites, c(20, 40))
  expect_equal(res$n_per_site, c(10, 10))
  expect_equal(res$n_per_group, c(200, 400))
  expect_equal(res$deff, 1 + 9 * c(0.02, 0.05))
  # one fewer person per site would miss the target
  m <- c(20, 40) * 9 / (1 + 8 * c(0.02, 0.05))
  z_a <- qnorm(0.975); z_b <- qnorm(0.8)
  expect_gt(.diff_gap(m, c(0.1, 0.2), z_a, z_b), 0)
})

test_that("more sites than needed: every site is used, with 1 person each", {
  res <- design_difference(0.10, 0.20, n_sites1 = 500, icc = 0.05)
  expect_equal(res$n_per_site, c(1, 1))
  expect_equal(res$n_sites, c(500, 500))
  expect_equal(res$n_per_group, c(500, 500))
  expect_equal(res$deff, c(1, 1))
  expect_equal(res$n_eff, c(199, 199))
})

test_that("too few sites stops, giving the highest power they can reach", {
  expect_error(design_difference(0.10, 0.20, n_sites1 = 8, icc = 0.05), "70.9%")
})

test_that("invalid inputs are rejected", {
  expect_error(design_difference(0.2, 0.2), "differ")
  expect_error(design_difference(0.1, 0.2, alternative = "greater"), "above")
  expect_error(design_difference(0.1, 0.2, power = 0.4), "power")
  expect_error(design_difference(0.1, 0.2, conf_level = 0.3), "conf_level")
  expect_error(design_difference(0.1, 0.2, icc = 0.05), "needs")
  expect_error(design_difference(0.1, 0.2, icc = c(0.1, 0.1, 0.1)), "one value")
  expect_error(design_difference(0.1, 0.2, power = 0.9, n1 = 100), "not both")
  expect_error(design_difference(0.1, 0.2, n1 = 100, n2 = NULL), "both groups, or for neither")
  expect_error(design_difference(0.1, 0.2, n_sites1 = 10, n_sites2 = NULL,
                                 n_per_site2 = 20, icc = 0.02), "leave out the same thing")
  expect_error(design_difference(0.1, 0.2, n1 = 100, n2 = NULL,
                                 n_sites2 = 5, n_per_site2 = 20), "sites for both groups")
})

test_that("group-2 arguments default to group 1", {
  expect_equal(design_difference(0.1, 0.2, n_sites1 = 30, icc = 0.02),
               design_difference(0.1, 0.2, n_sites1 = 30, n_sites2 = 30, icc = 0.02))
})

test_that("a full design returns its power", {
  # no sites: power.prop.test() at the same n
  for (n in c(100, 199, 300))
    expect_equal(design_difference(0.10, 0.20, n1 = n)$power,
                 power.prop.test(n = n, p1 = 0.10, p2 = 0.20)$power, tolerance = 1e-4)
  # different group sizes: .diff_power() at those sizes
  res <- design_difference(0.10, 0.20, n1 = 150, n2 = 300)
  expect_equal(res$n_per_group, c(150, 300))
  expect_equal(res$power, .diff_power(c(150, 300), c(0.1, 0.2), qnorm(0.975)))
  # with sites: effective sizes n / deff
  res <- design_difference(0.10, 0.20, n_sites1 = 20, n_per_site1 = 10, icc = 0.02)
  expect_equal(res$deff, c(1.18, 1.18))
  expect_equal(res$power, .diff_power(c(200, 200) / 1.18, c(0.1, 0.2), qnorm(0.975)))
})

test_that("sites of different sizes within each group", {
  s1 <- c(80, 60, 40, 40); s2 <- c(100, 50, 50)
  res <- design_difference(0.10, 0.20, n_per_site1 = s1, n_per_site2 = s2, icc = 0.02)
  d   <- 1 + (c(sum(s1^2) / sum(s1), sum(s2^2) / sum(s2)) - 1) * 0.02
  expect_equal(res$deff, d)
  expect_equal(res$n_per_group, c(220, 200))
  expect_equal(res$n_sites, c(4, 3))
  expect_equal(res$n_per_site, list(s1, s2))
  expect_equal(res$power, .diff_power(c(220, 200) / d, c(0.1, 0.2), qnorm(0.975)))
  # equal sizes as vectors are the same as n_sites and n_per_site
  expect_equal(design_difference(0.1, 0.2, n_per_site1 = rep(10, 20), icc = 0.02)$power,
               design_difference(0.1, 0.2, n_sites1 = 20, n_per_site1 = 10, icc = 0.02)$power)
})

test_that("a solved design reports its own power, at least the target", {
  for (args in list(list(), list(n_per_site1 = 30, icc = 0.02),
                    list(n_sites1 = 20, n_sites2 = 40, icc = c(0.02, 0.05)))) {
    res <- do.call(design_difference, c(list(0.10, 0.20), args))
    expect_gte(res$power, 0.8)
    back <- if (is.null(res$n_sites)) design_difference(0.1, 0.2, n1 = res$n_per_group[1])
            else design_difference(0.1, 0.2, n_sites1 = res$n_sites[1],
                                   n_sites2 = res$n_sites[2],
                                   n_per_site1 = res$n_per_site[1],
                                   n_per_site2 = res$n_per_site[2], icc = args$icc)
    expect_equal(back$power, res$power)
  }
})
