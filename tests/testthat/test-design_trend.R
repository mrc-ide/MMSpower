# Hand-checked reference cases for design_trend()
#
#   S_xx   = sum((t - mean(t))^2)
#   sigma2 = pbar_app * (1 - pbar_app), pbar_app = (p_start_app + p_end_app) / 2
#   n0     = (z_a + z_b)^2 * sigma2 / (S_xx * slope_app^2)       (per timepoint)
#   z_a    = qnorm(0.975) = 1.959964 (two-sided, 95%); z_b = qnorm(0.80) = 0.841621
#
# DTR-1 (10 to 20 percent over 5 equally spaced rounds, t = 0, 1, 2, 3, 4):
#   S_xx = 10; slope_app = 0.1 / 4 = 0.025; pbar = 0.15; sigma2 = 0.1275
#   n0 = 2.801585^2 * 0.1275 / (10 * 0.025^2) = 160.12 -> 161 per round
#
# DTR-5 (reverse mode, 150 per round):
#   ncp = 0.025 * sqrt(150 * 10 / 0.1275) = 2.7116
#   power = pnorm(2.7116 - 1.959964) + pnorm(-2.7116 - 1.959964) = 0.7739


# ---------------------------------------------------------------------------
# Core scenarios (DTR-)
# ---------------------------------------------------------------------------

test_that("DTR-1: 10 to 20 percent over 5 rounds gives the hand-checked 161 per round", {
  res <- design_trend(0.10, 0.20, n_timepoints = 5)
  expect_equal(res$n_per_timepoint, 161)
  expect_equal(res$n,               805)
  expect_equal(res$n_eff,           161)
  expect_equal(res$sxx,             10)
  expect_equal(res$slope,           0.025)
  expect_equal(res$slope_app,       0.025)
  expect_equal(res$sigma2_app,      0.1275)
  expect_equal(res$times,           0:4)
  expect_equal(res$mode,            "solve_n")
  expect_null(res$n_sites)
})

test_that("DTR-2: giving the slope instead of the last prevalence gives the same design", {
  a <- design_trend(0.10, slope = 0.025, n_timepoints = 5)
  b <- design_trend(0.10, 0.20,          n_timepoints = 5)
  expect_equal(a, b)
  expect_equal(a$prevalence_end, 0.20)
})

test_that("DTR-3: more power needs more people; a one-sided test needs fewer", {
  expect_equal(design_trend(0.10, 0.20, n_timepoints = 5, power = 0.90)$n_per_timepoint, 215)
  expect_equal(design_trend(0.10, 0.20, n_timepoints = 5,
                            alternative = "greater")$n_per_timepoint, 127)
  # A falling trend with "less": 30 to 20 percent over 4 rounds (S_xx = 5)
  expect_equal(design_trend(0.30, 0.20, n_timepoints = 4,
                            alternative = "less")$n_per_timepoint, 209)
})

test_that("DTR-4: custom times change S_xx; times 0 to 4 equal n_timepoints = 5", {
  # Years 2018, 2019, 2021, 2024: mean 2020.5, S_xx = 6.25 + 2.25 + 0.25 + 12.25 = 21
  res <- design_trend(0.10, 0.20, times = c(2018, 2019, 2021, 2024))
  expect_equal(res$sxx,             21)
  expect_equal(res$slope,           0.1 / 6)
  expect_equal(res$n_per_timepoint, 172)
  expect_equal(design_trend(0.10, 0.20, times = 0:4)$n_per_timepoint, 161)
})

test_that("DTR-5: reverse mode reports the hand-checked power for 150 per round", {
  res <- design_trend(0.10, 0.20, n_timepoints = 5, n_per_timepoint = 150)
  expect_equal(res$power, 0.7739, tolerance = 1e-3)
  expect_equal(res$n,     750)
  expect_equal(res$mode,  "solve_power")
  # With 1 person per round the power is close to the 5 percent false-alarm
  # rate, and both tails count: ncp = 0.025 * sqrt(10 / 0.1275) = 0.2214,
  # pnorm(0.2214 - 1.96) + pnorm(-0.2214 - 1.96) = 0.041056 + 0.014578
  expect_equal(design_trend(0.10, 0.20, n_timepoints = 5, n_per_timepoint = 1)$power,
               0.055634, tolerance = 1e-4)
})

test_that("DTR-6: reverse mode gives back the target power for every design", {
  # The forward n reaches the target power; one person fewer falls short
  for (s in list(list(), list(power = 0.9), list(alternative = "greater"),
                 list(sensitivity = 0.9, specificity = 0.95),
                 list(n_per_site = 20, icc = 0.03), list(n_sites = 10, icc = 0.03),
                 list(fpc_N = 1000), list(n_per_site = 20, icc = 0.03, fpc_N = 1000),
                 list(n_sites = 10, icc = 0.03, fpc_N = 1000))) {
    fwd <- do.call(design_trend, c(list(0.10, 0.20, n_timepoints = 5), s))
    rev_args <- c(list(0.10, 0.20, n_timepoints = 5), s[setdiff(names(s), "power")])
    at  <- do.call(design_trend, c(rev_args, list(n_per_timepoint = fwd$n_per_timepoint)))
    one <- do.call(design_trend, c(rev_args, list(n_per_timepoint = fwd$n_per_timepoint - 1)))
    expect_gte(at$power, fwd$power - 1e-9)
    expect_lt(one$power, fwd$power)
  }
})

test_that("DTR-7: an imperfect test shrinks the measured slope and needs more people", {
  # slope_app = 0.025 * (0.9 + 0.95 - 1) = 0.02125
  res <- design_trend(0.10, 0.20, n_timepoints = 5, sensitivity = 0.9, specificity = 0.95)
  expect_equal(res$slope_app,       0.02125)
  expect_equal(res$n_per_timepoint, 254)
})


# ---------------------------------------------------------------------------
# Clustering (DTR-C-)
# ---------------------------------------------------------------------------

test_that("DTR-C-1: a fixed number of people per site multiplies n by the design effect", {
  # deff = 1 + 19 * 0.03 = 1.57; n = 160.12 * 1.57 = 251.4 -> 252; 13 sites
  res <- design_trend(0.10, 0.20, n_timepoints = 5, n_per_site = 20, icc = 0.03)
  expect_equal(res$deff,            1.57)
  expect_equal(res$n_per_timepoint, 252)
  expect_equal(res$n_sites,         13)
  expect_equal(res$n_eff,           161)
})

test_that("DTR-C-2: a fixed number of sites solves for people per site, and deff matches it", {
  # n = 160.12 * 10 * 0.97 / (10 - 160.12 * 0.03) = 298.8 -> 299; 30 per site
  res <- design_trend(0.10, 0.20, n_timepoints = 5, n_sites = 10, icc = 0.03)
  expect_equal(res$n_per_timepoint, 299)
  expect_equal(res$n_per_site,      30)
  expect_equal(res$deff,            1 + (30 - 1) * 0.03)
})

test_that("DTR-C-3: too few sites for the target power is rejected, giving the maximum power", {
  expect_error(
    design_trend(0.10, 0.20, n_timepoints = 5, n_sites = 3, icc = 0.05),
    "Target power of 80.0 percent is unachievable with 3 sites per timepoint and ICC = 0.050.\nMaximum achievable power with these settings: 40.3 percent.\nIncrease `n_sites`, lower the ICC assumption, or relax the power target.",
    fixed = TRUE)
  # With 5 sites the maximum is 60.03 percent: asking for 0.600 works (just),
  # asking for 0.601 does not
  expect_no_error(design_trend(0.10, 0.20, n_timepoints = 5, n_sites = 5, icc = 0.05, power = 0.600))
  expect_error(design_trend(0.10, 0.20, n_timepoints = 5, n_sites = 5, icc = 0.05, power = 0.601),
               "Maximum achievable power with these settings: 60.0 percent.", fixed = TRUE)
})

test_that("DTR-C-4: more sites than people needed with no clustering is rejected", {
  expect_error(
    design_trend(0.10, 0.20, n_timepoints = 5, n_sites = 200, icc = 0.05),
    "n_sites = 200 is >= the number of people needed per timepoint with no clustering (about 161), so each site would receive fewer than 1 person on average, which is not a valid cluster design. Use n_sites < 161, or supply `n_per_site` to fix the cluster size and solve for the number of sites instead.",
    fixed = TRUE)
})

test_that("DTR-C-5: icc = 0, or a tiny icc, with sites is the unclustered result spread over the sites", {
  for (icc in c(0, 1e-12)) {
    res <- design_trend(0.10, 0.20, n_timepoints = 5, n_sites = 10, icc = icc)
    expect_equal(res$n_per_timepoint, 161)
    expect_equal(res$deff,            1)
    expect_equal(res$n_per_site,      17)
  }
  expect_equal(design_trend(0.10, 0.20, n_timepoints = 5, icc = 1e-12)$n_per_timepoint, 161)
})


# ---------------------------------------------------------------------------
# Finite population correction (DTR-F-)
# ---------------------------------------------------------------------------

test_that("DTR-F-1: a small population reduces n, and leaves n_eff alone", {
  # 160.12 * 1000 / (160.12 + 999) = 138.1 -> 139
  res <- design_trend(0.10, 0.20, n_timepoints = 5, fpc_N = 1000)
  expect_equal(res$n_per_timepoint, 139)
  expect_equal(res$n_eff,           161)
})

test_that("DTR-F-2: with n_sites fixed, the FPC and the design effect are solved together", {
  n0 <- (stats::qnorm(0.975) + stats::qnorm(0.80))^2 * 0.1275 / (10 * 0.025^2)
  for (case in list(c(sites = 10, icc = 0.03, N = 1000), c(sites = 3, icc = 0.05, N = 2000),
                    c(sites = 25, icc = 0.10, N = 5000))) {
    k <- case[["sites"]]; r <- case[["icc"]]; N <- case[["N"]]
    eq   <- function(n) n0 * (1 + (n / k - 1) * r) * (N - n) / (N - 1) - n
    root <- stats::uniroot(eq, c(1e-9, N), tol = 1e-10)$root
    expect_equal(design_trend(0.10, 0.20, n_timepoints = 5, n_sites = k, icc = r,
                              fpc_N = N)$n_per_timepoint, ceiling(root))
  }
  expect_equal(design_trend(0.10, 0.20, n_timepoints = 5, n_sites = 10, icc = 0.03,
                            fpc_N = 1e9)$n_per_timepoint, 299)
})

test_that("DTR-F-3: a target unachievable with n_sites alone becomes achievable with a small population", {
  # 3 sites: unachievable without fpc_N (DTR-C-3); with 2000 people per round
  # the joint solve gives 1283 (428 per site)
  res <- design_trend(0.10, 0.20, n_timepoints = 5, n_sites = 3, icc = 0.05, fpc_N = 2000)
  expect_equal(res$n_per_timepoint, 1283)
  expect_equal(res$n_per_site,      428)
})

test_that("DTR-F-4: a rounded design larger than the population warns", {
  # Population 50, 13 sites: n = 40, 40 / 13 = 3.08 -> 4 per site,
  # so 13 x 4 = 52 people per round, more than the 50 who exist
  expect_warning(
    res <- design_trend(0.10, 0.20, n_timepoints = 5, n_sites = 13, icc = 0.05, fpc_N = 50),
    "The rounded design (13 sites x 4 = 52 per timepoint) exceeds `fpc_N` = 50: the required n of 40 per timepoint is close to a census of the population. Consider surveying the whole population instead.",
    fixed = TRUE)
  expect_equal(res$n_per_timepoint, 40)
})


# ---------------------------------------------------------------------------
# Validation (DTR-V-)
# ---------------------------------------------------------------------------

test_that("DTR-V-1: prevalences of 0, 1, below 0 or above 1 are rejected", {
  expect_error(design_trend(0, 0.2, n_timepoints = 5),
               "`prevalence_start` must be in (0, 1) (got 0). Use a value from a pilot study, historical data, or a conservative guess.",
               fixed = TRUE)
  expect_error(design_trend(0.1, 1.2, n_timepoints = 5),
               "`prevalence_end` must be in (0, 1) (got 1.2). `prevalence_end` is a fraction of the population and cannot be greater than 1.",
               fixed = TRUE)
})

test_that("DTR-V-2: exactly one of prevalence_end or slope, and the slope must keep prevalence inside (0, 1)", {
  expect_error(design_trend(0.1, n_timepoints = 5),
               "Supply exactly one of `prevalence_end` or `slope` to define the trend.", fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, slope = 0.02, n_timepoints = 5),
               "Supply exactly one of `prevalence_end` or `slope`, not both.", fixed = TRUE)
  expect_error(design_trend(0.1, slope = TRUE, n_timepoints = 5),
               "`slope` must be a single number (got class `logical`, length 1).", fixed = TRUE)
  expect_error(design_trend(0.1, slope = Inf, n_timepoints = 5),
               "`slope` must be a single finite number (got Inf).", fixed = TRUE)
  expect_error(design_trend(0.1, slope = 0.3, n_timepoints = 5),
               "`slope` = 0.3 over a time span of 4 takes the prevalence to 1.3, outside (0, 1). Reduce the slope or the time span.",
               fixed = TRUE)
})

test_that("DTR-V-3: n_timepoints must be an integer of 2 or more", {
  expect_error(design_trend(0.1, 0.2),
               "Supply `n_timepoints` (equally spaced timepoints) or `times`.", fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, n_timepoints = 1),
               "`n_timepoints` must be a single integer of 2 or more (got 1), or supply `times`.", fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, n_timepoints = "5"),
               "`n_timepoints` must be a single integer of 2 or more (got class `character`), or supply `times`.",
               fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, n_timepoints = 4.5), "(got 4.5)", fixed = TRUE)
  # Valid: 2 timepoints works
  expect_no_error(design_trend(0.1, 0.2, n_timepoints = 2))
})

test_that("DTR-V-4: times must be at least two finite, different values; n_timepoints is then ignored with a warning", {
  expect_error(design_trend(0.1, 0.2, times = 2020),
               "`times` must be a numeric vector of at least two timepoints (got class `numeric`, length 1).",
               fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, times = c(2020, NA, 2022)),
               "`times` contains a missing or infinite value (found times[2] = NA).", fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, times = c(2020, 2020)),
               "`times` must contain at least two different values (got all 2020).", fixed = TRUE)
  expect_warning(res <- design_trend(0.1, 0.2, times = 0:4, n_timepoints = 3),
                 "Both `times` and `n_timepoints` were supplied; `n_timepoints` is ignored (the timepoints come from `times`).",
                 fixed = TRUE)
  expect_equal(res$n_timepoints, 5)
})

test_that("DTR-V-5: n_per_timepoint must be a positive integer that fits the sites and the population", {
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, n_per_timepoint = "100"),
               "`n_per_timepoint` must be a single finite positive integer (got class `character`).", fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, n_per_timepoint = 0),
               "`n_per_timepoint` must be a single finite positive integer (got 0).", fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, n_per_timepoint = 5, n_sites = 10, icc = 0.03),
               "`n_per_timepoint` (5) is smaller than `n_sites` (10): there would be fewer than one person per site. Increase `n_per_timepoint`, or lower `n_sites`.",
               fixed = TRUE)
  # With n_per_site fixed, fewer people than one full site is allowed: forward
  # mode can return it (a large trend needs fewer than 20 people per round)
  fwd <- design_trend(0.05, 0.60, n_timepoints = 5, n_per_site = 20, icc = 0.03)
  expect_lt(fwd$n_per_timepoint, 20)
  expect_gte(design_trend(0.05, 0.60, n_timepoints = 5, n_per_site = 20, icc = 0.03,
                          n_per_timepoint = fwd$n_per_timepoint)$power, 0.80)
  # With n_sites fixed and icc = 0, fewer people than sites is fine too
  expect_no_error(design_trend(0.1, 0.2, n_timepoints = 5, n_per_timepoint = 5, n_sites = 10))
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, n_per_timepoint = 101, fpc_N = 100),
               "`n_per_timepoint` (101) is larger than `fpc_N` (100): you cannot test more people than there are in the population.",
               fixed = TRUE)
  # Testing the whole population is allowed: no sampling uncertainty, power 1.
  # Forward mode can return it when the required n rounds up to fpc_N.
  census <- design_trend(0.1, 0.2, n_timepoints = 5, n_per_timepoint = 100, fpc_N = 100)
  expect_equal(census$power, 1)
  expect_equal(census$n_eff, Inf)
  # In reverse mode power is not used, so a bad power is not checked
  expect_no_error(design_trend(0.1, 0.2, n_timepoints = 5, n_per_timepoint = 150, power = 2))
})

test_that("DTR-V-6: power and conf_level outside their ranges are rejected", {
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, power = 0.3),
               "`power` must be at least 0.5 (got 0.3).", fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, power = 1),
               "`power` must be in [0.5, 1) (got 1). For example, use 0.80 for 80 percent power.", fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, power = TRUE),
               "`power` must be a single number in [0.5, 1) (got class `logical`, length 1).", fixed = TRUE)
  # Below 0.5 a one-sided design would come out absurdly small (3 per round)
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, alternative = "greater", conf_level = 0.3),
               "`conf_level` must be at least 0.5 (got 0.3). Below 0.5 the sample-size formula stops increasing with the confidence level, so it no longer gives a sensible answer. Pass a confidence level such as 0.95, not a significance level such as 0.05.",
               fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, conf_level = 1),
               "`conf_level` must be in [0.5, 1) (got 1).", fixed = TRUE)
  expect_no_error(design_trend(0.1, 0.2, n_timepoints = 5, conf_level = 0.5))
})

test_that("DTR-V-7: alternative must be a valid name and match the direction of the trend", {
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, alternative = "up"),
               "`alternative` must be one of 'two.sided', 'greater', or 'less' (got 'up').", fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, alternative = "less"),
               "alternative = 'less' needs a decreasing trend, but the slope is 0.025.", fixed = TRUE)
  expect_error(design_trend(0.2, 0.1, n_timepoints = 5, alternative = "greater"),
               "alternative = 'greater' needs an increasing trend, but the slope is -0.025.", fixed = TRUE)
  expect_error(design_trend(0.1, 0.1, n_timepoints = 5),
               "The slope is zero, or nearly so, so there is no trend to detect.", fixed = TRUE)
  # A flat trend is not "decreasing" (or "increasing") either
  expect_error(design_trend(0.1, 0.1, n_timepoints = 5, alternative = "less"),
               "alternative = 'less' needs a decreasing trend, but the slope is 0.", fixed = TRUE)
  # A slope that is not zero but far too small to detect
  expect_error(design_trend(0.1, 0.1 + 2e-10, n_timepoints = 5),
               "so the required sample size would be astronomically large. Use a larger trend, or a more accurate test.",
               fixed = TRUE)
})

test_that("DTR-V-8: sensitivity and specificity: type, range, sum of 1 or less, and the warning near 1", {
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, sensitivity = NA),
               "`sensitivity` must be a single number in (0, 1] (got class `logical`, length 1). Note: `TRUE`, `FALSE` and a plain `NA` are logical, not numeric. Pass 1 if the test never misses a case.",
               fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, specificity = 0),
               "`specificity` must be in (0, 1] (got 0). A specificity of 0 means everyone who truly does not have the condition tests positive.",
               fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, sensitivity = 0.3,
                            specificity = 0.7 + 2 * .Machine$double.eps),
               "must exceed 1 for the Rogan-Gladen correction", fixed = TRUE)
  expect_warning(design_trend(0.1, 0.2, n_timepoints = 5, sensitivity = 0.55, specificity = 0.5),
                 "`sensitivity` + `specificity` = 1.05, which is very close to 1.", fixed = TRUE)
})

test_that("DTR-V-9: a vector, text, TRUE/FALSE, a bare NA or NULL in a single-number argument is rejected, naming the class", {
  expect_error(design_trend(c(0.1, 0.2), 0.3, n_timepoints = 5),
               "`prevalence_start` must be a single number in (0, 1) (got class `numeric`, length 2).", fixed = TRUE)
  expect_error(design_trend(0.1, "0.2", n_timepoints = 5),
               "`prevalence_end` must be a single number in (0, 1) (got class `character`, length 1).", fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, specificity = TRUE),
               "`specificity` must be a single number in (0, 1] (got class `logical`, length 1). Note: `TRUE`, `FALSE` and a plain `NA` are logical, not numeric. Pass 1 if the test never gives a false positive.",
               fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, conf_level = NULL),
               "`conf_level` must be a single number in [0.5, 1) (got class `NULL`, length 0).", fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, icc = list(0.05)),
               "`icc` must be a single number in [0, 1] (got class `list`, length 1).", fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, alternative = NA),
               "`alternative` must be one of 'two.sided', 'greater', or 'less' (got class `logical`, length 1).",
               fixed = TRUE)
})

test_that("DTR-V-10: NaN, Inf or -Inf in a single-number argument gives the finite-value error", {
  expect_error(design_trend(NaN, 0.2, n_timepoints = 5),
               "`prevalence_start` must be a single finite number (got NaN).", fixed = TRUE)
  expect_error(design_trend(0.1, Inf, n_timepoints = 5),
               "`prevalence_end` must be a single finite number (got Inf).", fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, power = -Inf),
               "`power` must be a single finite number (got -Inf).", fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, sensitivity = NaN),
               "`sensitivity` must be a single finite number (got NaN).", fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, specificity = Inf),
               "`specificity` must be a single finite number (got Inf).", fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, conf_level = NaN),
               "`conf_level` must be a single finite number (got NaN).", fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, icc = Inf),
               "`icc` must be a single finite number (got Inf).", fixed = TRUE)
})

test_that("DTR-V-11: sensitivity and specificity range checks give the right message for each bad value", {
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, sensitivity = 1.1),
               "`sensitivity` must be in (0, 1] (got 1.1). `sensitivity` is a diagnostic probability and cannot be greater than 1.",
               fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, sensitivity = 0),
               "`sensitivity` must be in (0, 1] (got 0). A sensitivity of 0 means everyone who truly has the condition tests negative.",
               fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, specificity = -0.5),
               "`specificity` must be in (0, 1] (got -0.5). `specificity` is a diagnostic probability and cannot be negative.",
               fixed = TRUE)
})

test_that("DTR-V-12: icc, n_sites, n_per_site and fpc_N are checked, naming the problem", {
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, icc = 0.05),
               "`icc` is greater than 0 but there is no cluster structure", fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, n_sites = 10, icc = 1.5),
               "`icc` must be in [0, 1] (got 1.5). `icc` is a correlation and cannot be greater than 1.", fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, n_sites = "10", icc = 0.05),
               "`n_sites` must be a single finite positive integer (got class `character`). `n_sites` is the number of sampling clusters per timepoint.",
               fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, n_per_site = 2.5, icc = 0.05),
               "`n_per_site` must be a single finite positive integer (got 2.5).", fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, fpc_N = TRUE),
               "`fpc_N` must be a single finite positive integer (got class `logical`). `fpc_N` is the population size per timepoint.",
               fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, n_sites = 10, n_per_site = 20, icc = 0.05),
               "Supply at most one of `n_sites` or `n_per_site`, not both.", fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, n_sites = 100, icc = 0.05, fpc_N = 50),
               "`n_sites` (100) exceeds `fpc_N` (50)", fixed = TRUE)
  expect_error(design_trend(0.1, 0.2, n_timepoints = 5, n_per_site = 100, icc = 0.05, fpc_N = 50),
               "`n_per_site` (100) exceeds `fpc_N` (50): a cluster cannot be larger than the whole population.",
               fixed = TRUE)
})


# ---------------------------------------------------------------------------
# Return value (DTR-R-)
# ---------------------------------------------------------------------------

test_that("DTR-R-1: forward and reverse mode return the same fields, in order", {
  fields <- c("n_per_timepoint", "n", "n_eff", "n_timepoints", "times", "power",
              "prevalence_start", "prevalence_end", "slope", "slope_app", "sxx",
              "sigma2_app", "alternative", "conf_level", "sensitivity", "specificity",
              "icc", "deff", "fpc_N", "mode", "n_sites", "n_per_site")
  expect_named(design_trend(0.1, 0.2, n_timepoints = 5), fields, ignore.order = FALSE)
  expect_named(design_trend(0.1, 0.2, n_timepoints = 5, n_per_timepoint = 150), fields,
               ignore.order = FALSE)
})
