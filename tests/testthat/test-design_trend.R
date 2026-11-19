# Hand-checked reference cases for design_trend()
#
#   S_xx   = sum((t - mean(t))^2)
#   sigma2 = pbar_app (1 - pbar_app), pbar_app = (p_start_app + p_end_app) / 2
#   n0     = (z_a + z_b)^2 sigma2 / (S_xx slope_app^2)        (per timepoint)
#   z_a = qnorm(0.975) = 1.959964 (two-sided, 95 percent); z_b = qnorm(0.8) = 0.841621
#
# 10 to 20 percent over 5 rounds (t = 0 to 4): S_xx = 10, slope = 0.025,
#   sigma2 = 0.15 * 0.85 = 0.1275, n0 = 2.801585^2 * 0.1275 / (10 * 0.025^2)
#   = 160.12 -> 161 per round
# Reverse, 150 per round: ncp = 0.025 * sqrt(150 * 10 / 0.1275) = 2.7116,
#   power = pnorm(2.7116 - 1.96) + pnorm(-2.7116 - 1.96) = 0.7739
# n_per_site = 20, icc = 0.03: deff = 1.57, n = 160.12 * 1.57 = 251.4 ->
#   13 sites x 20 = 260
# n_sites = 10, icc = 0.03: n = 160.12 * 10 * 0.97 / (10 - 160.12 * 0.03)
#   = 298.8 -> 30 per site x 10 = 300

test_that("10 to 20 percent over 5 rounds gives the hand-checked 161 per round", {
  res <- design_trend(0.10, 0.20, n_timepoints = 5)
  expect_equal(c(res$n_per_timepoint, res$n, res$n_eff, res$deff), c(161, 805, 161, 1))
  expect_null(res$n_sites)
  expect_named(res, c("n_per_timepoint", "n", "n_sites", "n_per_site", "deff",
                      "n_eff", "power"))
  # the same trend given as a slope, or as times 0 to 4
  expect_equal(design_trend(0.10, slope = 0.025, n_timepoints = 5), res)
  expect_equal(design_trend(0.10, 0.20, times = 0:4), res)
})

test_that("power, one-sided tests, uneven times and an imperfect test", {
  expect_equal(design_trend(0.10, 0.20, n_timepoints = 5, power = 0.9)$n_per_timepoint, 215)
  expect_equal(design_trend(0.10, 0.20, n_timepoints = 5,
                            alternative = "greater")$n_per_timepoint, 127)
  expect_equal(design_trend(0.30, 0.20, n_timepoints = 4,
                            alternative = "less")$n_per_timepoint, 209)
  # 2018, 2019, 2021, 2024: S_xx = 21
  expect_equal(design_trend(0.10, 0.20, times = c(2018, 2019, 2021, 2024))$n_per_timepoint, 172)
  # measured slope 0.025 * 0.85
  expect_equal(design_trend(0.10, 0.20, n_timepoints = 5, sensitivity = 0.9,
                            specificity = 0.95)$n_per_timepoint, 254)
})

test_that("reverse: hand-checked power, and the forward n gives back the target", {
  res <- design_trend(0.10, 0.20, n_timepoints = 5, n_per_timepoint = 150)
  expect_equal(res$power, 0.7739, tolerance = 1e-3)
  expect_equal(res$n, 750)
  # 1 person per round: both tails count, 0.041056 + 0.014578
  expect_equal(design_trend(0.10, 0.20, n_timepoints = 5, n_per_timepoint = 1)$power,
               0.055634, tolerance = 1e-4)
  for (s in list(list(), list(alternative = "greater"),
                 list(sensitivity = 0.9, specificity = 0.95),
                 list(n_per_site = 20, icc = 0.03))) {
    fwd <- do.call(design_trend, c(list(0.10, 0.20, n_timepoints = 5), s))
    at  <- do.call(design_trend, c(list(0.10, 0.20, n_timepoints = 5,
                                        n_per_timepoint = fwd$n_per_timepoint), s))
    expect_gte(at$power, 0.8)
  }
})

test_that("fixed site size solves for the number of sites", {
  res <- design_trend(0.10, 0.20, n_timepoints = 5, n_per_site = 20, icc = 0.03)
  expect_equal(c(res$n_per_timepoint, res$n_sites, res$n_per_site, res$deff, res$n_eff),
               c(260, 13, 20, 1.57, 161))
})

test_that("fixed number of sites solves for people per site", {
  res <- design_trend(0.10, 0.20, n_timepoints = 5, n_sites = 10, icc = 0.03)
  expect_equal(c(res$n_per_timepoint, res$n_sites, res$n_per_site, res$deff),
               c(300, 10, 30, 1 + 29 * 0.03))
  # icc 0: 161 spread over 10 sites, 17 each
  res <- design_trend(0.10, 0.20, n_timepoints = 5, n_sites = 10)
  expect_equal(c(res$n_per_timepoint, res$n_per_site, res$deff), c(170, 17, 1))
})

test_that("too few sites stops, giving the highest power they can give", {
  expect_error(design_trend(0.10, 0.20, n_timepoints = 5, n_sites = 3, icc = 0.05),
               "the highest power these sites can give is 40.3 percent")
  # 5 sites can give 60.03 percent
  expect_no_error(design_trend(0.10, 0.20, n_timepoints = 5, n_sites = 5, icc = 0.05,
                               power = 0.6))
  expect_error(design_trend(0.10, 0.20, n_timepoints = 5, n_sites = 5, icc = 0.05,
                            power = 0.601), "60 percent")
})

test_that("more sites than people needed: 1 person at just enough sites", {
  for (icc in c(0, 0.05)) {
    expect_warning(res <- design_trend(0.10, 0.20, n_timepoints = 5, n_sites = 200, icc = icc),
                   "needs 161 people, fewer than the 200 sites")
    expect_equal(c(res$n_per_timepoint, res$n_sites, res$n_per_site, res$deff),
                 c(161, 161, 1, 1))
  }
  # exactly 161 sites: no warning
  expect_no_warning(res <- design_trend(0.10, 0.20, n_timepoints = 5, n_sites = 161,
                                        icc = 0.05))
  expect_equal(res$n_sites, 161)
})

test_that("reverse with sites uses the design effect of that design", {
  pw <- function(n_eff) pnorm(0.025 * sqrt(n_eff * 10 / 0.1275) - qnorm(0.975)) +
    pnorm(-0.025 * sqrt(n_eff * 10 / 0.1275) - qnorm(0.975))
  # 260 in sites of 20, ICC 0.03: deff 1.57
  res <- design_trend(0.10, 0.20, n_timepoints = 5, n_per_timepoint = 260,
                      n_per_site = 20, icc = 0.03)
  expect_equal(c(res$deff, res$n_eff, res$n_sites), c(1.57, 260 / 1.57, 13))
  expect_equal(res$power, pw(260 / 1.57))
  # 300 at 10 sites, ICC 0.03: 30 per site, deff 1.87
  res <- design_trend(0.10, 0.20, n_timepoints = 5, n_per_timepoint = 300,
                      n_sites = 10, icc = 0.03)
  expect_equal(c(res$deff, res$n_per_site), c(1.87, 30))
  expect_equal(res$power, pw(300 / 1.87))
})

test_that("bad inputs give a clear error", {
  expect_error(design_trend(0, 0.2, n_timepoints = 5), "`prevalence_start` must be")
  expect_error(design_trend(0.1, n_timepoints = 5), "exactly one of `prevalence_end` and `slope`")
  expect_error(design_trend(0.1, 0.2, slope = 0.1, n_timepoints = 5), "exactly one of")
  expect_error(design_trend(0.1, slope = 0.3, n_timepoints = 5), "outside (0, 1)", fixed = TRUE)
  expect_error(design_trend(0.1, 0.1, n_timepoints = 5), "must change")
  expect_error(design_trend(0.1, 0.2), "exactly one of `n_timepoints` and `times`")
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, times = 0:4), "exactly one of")
  expect_error(design_trend(0.1, 0.2, n_timepoints = 1), "`n_timepoints` must be")
  expect_error(design_trend(0.1, 0.2, times = c(1, 1)), "`times` must be")
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, power = 0.3), "`power` must be")
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, alternative = "less"), "must fall")
  expect_error(design_trend(0.2, 0.1, n_timepoints = 5, alternative = "greater"), "must rise")
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, icc = 0.1), "needs `n_sites`")
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, n_per_timepoint = 2, n_sites = 5,
                            icc = 0.1), "at least `n_sites`")
  # power is not used, so not checked, in reverse
  expect_no_error(design_trend(0.1, 0.2, n_timepoints = 5, n_per_timepoint = 50, power = 2))
})
