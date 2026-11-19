# Hand-checked reference cases for design_detection()
#
# q = prevalence * sensitivity, n0 = log(1 - detection_prob) / log(1 - q)
#
# 5 percent, 80 percent (Module 7 slide): log(0.2) / log(0.95) = 31.38 -> 32
# 139 people (Module 7 slide): 1 - 0.99^139 = 0.7527, 1 - 0.98^139 = 0.9397
#
# With sites (beta model, as DRpower's get_power_presence()):
#   a = q (1/icc - 1), b = (1 - q) (1/icc - 1)
#   P0(m) = prod over j = 0..m-1 of (b + j) / (a + b + j)
#   chance at k sites of m people: 1 - P0(m)^k
# n_per_site = 10, icc = 0.05: 1 - P0(10)^3 = 71.57, 1 - P0(10)^4 = 81.31
#   percent -> 4 sites, 40 people
# n_sites = 10, icc = 0.05: 1 - P0(3)^10 = 76.83, 1 - P0(4)^10 = 85.10
#   percent -> 4 per site, 40 people

p0_direct <- function(m, q, icc) {
  a <- q * (1 / icc - 1); b <- (1 - q) * (1 / icc - 1)
  prod((b + 0:(m - 1)) / (a + b + 0:(m - 1)))
}

test_that("unclustered: the Module 7 example gives the slide's 32", {
  res <- design_detection(0.05)
  expect_equal(res$n, 32)
  expect_equal(res$n_eff, 32)
  expect_equal(res$deff, 1)
  expect_equal(res$detection_prob, 0.8)
  expect_null(res$n_sites)
  expect_null(res$n_per_site)
  expect_named(res, c("n", "n_sites", "n_per_site", "deff", "n_eff", "detection_prob"))
})

test_that("a higher target, a rarer case or an imperfect test needs more people", {
  expect_equal(design_detection(0.05, 0.95)$n, 59)
  expect_equal(design_detection(0.01)$n, 161)
  expect_equal(design_detection(0.05, sensitivity = 0.8)$n, 40)
})

test_that("reverse: 139 people give the slide's 75 and 94 percent", {
  expect_equal(design_detection(0.01, n = 139)$detection_prob, 1 - 0.99^139)
  expect_equal(design_detection(0.02, n = 139)$detection_prob, 1 - 0.98^139)
})

test_that("forward and reverse agree, also when the target is exactly hit", {
  expect_gte(design_detection(0.05, n = 32)$detection_prob, 0.8)
  expect_lt(design_detection(0.05, n = 31)$detection_prob, 0.8)
})

test_that("fixed site size solves for the number of sites", {
  res <- design_detection(0.05, n_per_site = 10, icc = 0.05)
  expect_equal(c(res$n, res$n_sites, res$n_per_site, res$n_eff), c(40, 4, 10, 32))
  # 40 people at 4 sites do what 40 / 1.2234 = 32.7 people without sites do
  expect_equal(res$deff, 40 * log(0.95) / (4 * log(p0_direct(10, 0.05, 0.05))))
})

test_that("fixed number of sites solves for people per site", {
  res <- design_detection(0.05, n_sites = 10, icc = 0.05)
  expect_equal(c(res$n, res$n_sites, res$n_per_site), c(40, 10, 4))
  expect_equal(res$deff, 40 * log(0.95) / (10 * log(p0_direct(4, 0.05, 0.05))))
})

test_that("sites match DRpower", {
  # people per site: DRpower::get_sample_size_presence(n_clust,
  # target_power = 0.8, prevalence, ICC)
  for (s in list(c(5, 0.01, 0.10, 294), c(5, 0.01, 0.05, 82), c(10, 0.02, 0.10, 13),
                 c(3, 0.05, 0.20, 47),  c(20, 0.01, 0.05, 10), c(10, 0.05, 0.30, 6),
                 c(4, 0.10, 0.60, 101)))
    expect_equal(design_detection(s[2], n_sites = s[1], icc = s[3])$n_per_site, s[4])
  # reverse: DRpower::get_power_presence(rep(4, 10), 0.05, 0.05) = 85.1019
  expect_equal(round(100 * design_detection(0.05, n = 40, n_sites = 10,
                                            icc = 0.05)$detection_prob, 4), 85.1019)
})

test_that("reverse with sites spreads the people as evenly as possible", {
  # 43 people at 10 sites: 3 sites of 5 and 7 of 4
  expect_equal(design_detection(0.05, n = 43, n_sites = 10, icc = 0.05)$detection_prob,
               1 - p0_direct(5, 0.05, 0.05)^3 * p0_direct(4, 0.05, 0.05)^7)
  # 45 people in sites of 10: 4 full sites and one of 5
  res <- design_detection(0.05, n = 45, n_per_site = 10, icc = 0.05)
  expect_equal(res$detection_prob,
               1 - p0_direct(10, 0.05, 0.05)^4 * p0_direct(5, 0.05, 0.05))
  expect_equal(res$n_sites, 5)
  expect_equal(res$deff, 45 / (log(1 - res$detection_prob) / log(0.95)))
})

test_that("icc 0 with sites spreads the unclustered n", {
  # 32 people over 5 sites: 7 per site, 35 sampled
  res <- design_detection(0.05, n_sites = 5)
  expect_equal(c(res$n, res$n_per_site, res$deff), c(35, 7, 1))
})

test_that("more sites than needed: every site is used, with 1 person each", {
  for (icc in c(0, 0.05)) {
    res <- design_detection(0.05, n_sites = 100, icc = icc)
    expect_equal(c(res$n, res$n_sites, res$n_per_site, res$deff), c(100, 100, 1, 1))
    expect_equal(res$n_eff, 32)
  }
})

test_that("one site larger than needed gives 1 site of that size (as design_precision)", {
  # icc 0: 32 needed, so a site of 5000 is 1 site and 5000 people
  res <- design_detection(0.05, n_per_site = 5000)
  expect_equal(c(res$n, res$n_sites, res$n_per_site, res$deff), c(5000, 1, 5000, 1))
  # icc 0.05: one site of 81 is enough, 80 is not; a site of 200 is 1 site
  expect_lt(1 - p0_direct(80, 0.05, 0.05), 0.8)
  expect_equal(design_detection(0.05, n_per_site = 81, icc = 0.05)$n, 81)
  expect_equal(design_detection(0.05, n_per_site = 80, icc = 0.05)$n_sites, 2)
  res <- design_detection(0.05, n_per_site = 200, icc = 0.05)
  expect_equal(c(res$n, res$n_sites), c(200, 1))
})

test_that("a target out of reach of the sites stops", {
  # ICC 0.5, prevalence 1e-5: 5 sites stay far below 80 percent
  expect_error(design_detection(1e-5, n_sites = 5, icc = 0.5), "10^15", fixed = TRUE)
})

test_that("bad inputs give a clear error", {
  expect_error(design_detection(0), "`prevalence` must be a single number in (0, 1)", fixed = TRUE)
  expect_error(design_detection(0.05, 1), "`detection_prob` must be")
  expect_error(design_detection(0.05, sensitivity = 0), "`sensitivity` must be")
  expect_error(design_detection(0.05, n = 2.5), "`n` must be a single whole number")
  expect_error(design_detection(0.05, n_sites = 3, n_per_site = 3), "at most one")
  expect_error(design_detection(0.05, icc = 0.1), "needs `n_sites` or `n_per_site`")
  expect_error(design_detection(0.05, n_sites = 5, icc = 1), "`icc` must be")
  expect_error(design_detection(0.05, n = 2, n_sites = 5, icc = 0.1),
               "`n` must be at least `n_sites`")
  # detection_prob is not used, so not checked, in reverse
  expect_no_error(design_detection(0.05, detection_prob = "x", n = 10))
})
