# Reference cases for test_trend()
#
#   p_j   = X_j / N_j                         (measured prevalence)
#   w_j   = N_j / (deff * f_j)                (f_j: FPC factor, 1 without fpc_N)
#   slope = sum w (t - t_w) p / sum w (t - t_w)^2, the slope lm(p ~ t, weights = w) gives
#   z     = slope / sqrt(p_w (1 - p_w) / sum w (t - t_w)^2)   (Cochran-Armitage;
#           the same p-value as prop.trend.test() without sites or fpc_N)
#   interval: the slopes b with (slope - b)^2 <= z^2 V(b), V(b) the slope's
#           variance if the line through p_w at t_w with slope b were true.
#           Expected ends below were found by scanning slopes in steps of 1e-6.
#
# TT-1 (12, 15, 21, 24, 30 of 150 at t = 0 to 4): p = 0.08, 0.10, 0.14, 0.16, 0.20,
#   slope = 0.03, p_w = 102 / 750 = 0.136, z = 0.03 / sqrt(0.136 * 0.864 / 1500)
#   = 3.389538, two-sided p = 0.000700, interval [0.012693, 0.046789]


# ---------------------------------------------------------------------------
# Core scenarios (TT-)
# ---------------------------------------------------------------------------

test_that("TT-1: five rounds give the hand-checked slope, test and interval", {
  res <- test_trend(c(12, 15, 21, 24, 30), rep(150, 5), time = 0:4)
  expect_equal(res$slope,     0.03)
  expect_equal(res$statistic, 3.389538, tolerance = 1e-6)
  expect_equal(res$p_value,   0.000700, tolerance = 1e-3)
  expect_true(res$reject)
  expect_equal(c(res$ci_lower, res$ci_upper), c(0.012693, 0.046789), tolerance = 1e-4)
  expect_equal(res$times,     0:4)
  expect_equal(res$n_total,   750)
  expect_equal(res$method,    "cochran-armitage")
})

test_that("TT-2: the test is R's prop.trend.test and the slope is R's regression weighted by number tested", {
  x <- c(12, 15, 21, 24, 30); n <- c(150, 120, 200, 90, 160); t <- c(0, 1, 3, 4, 7)
  res <- test_trend(x, n, time = t)
  expect_equal(res$p_value, stats::prop.trend.test(x, n, t)$p.value)
  expect_equal(res$statistic^2, unname(stats::prop.trend.test(x, n, t)$statistic))
  expect_equal(res$slope_app, unname(stats::coef(stats::lm(I(x / n) ~ t, weights = n))[2]))
})

test_that("TT-3: the reject decision and interval follow alternative and conf_level", {
  res_g <- test_trend(c(12, 15, 21, 24, 30), rep(150, 5), time = 0:4, alternative = "greater")
  expect_equal(res_g$ci_upper, Inf)
  expect_equal(res_g$ci_lower, 0.015493, tolerance = 1e-4)
  expect_equal(res_g$p_value, 0.000350, tolerance = 1e-3)
  expect_true(res_g$reject)
  # The one-sided end uses z at 0.95 (1.644854), not at 0.975: it is the
  # lower end of the two-sided 90 percent interval
  expect_equal(res_g$ci_lower,
               test_trend(c(12, 15, 21, 24, 30), rep(150, 5), time = 0:4, conf_level = 0.90)$ci_lower)
  res_l <- test_trend(c(12, 15, 21, 24, 30), rep(150, 5), time = 0:4, alternative = "less")
  expect_equal(res_l$ci_lower, -Inf)
  expect_equal(res_l$ci_upper, 0.044142, tolerance = 1e-4)
  expect_false(res_l$reject)
  # A weak trend (two-sided p = 0.112): not significant at 0.05, significant at 0.20
  weak <- c(12, 14, 17, 17, 20)
  r95 <- test_trend(weak, rep(150, 5), time = 0:4)
  r80 <- test_trend(weak, rep(150, 5), time = 0:4, conf_level = 0.80)
  expect_equal(r95$p_value, 0.112008, tolerance = 1e-5)
  expect_false(r95$reject)
  expect_true(r80$reject)
  expect_equal(c(r95$ci_lower, r95$ci_upper), c(-0.002952, 0.028067), tolerance = 1e-3)
  expect_equal(c(r80$ci_lower, r80$ci_upper), c(0.002453, 0.022786), tolerance = 1e-3)
})

test_that("TT-4: an imperfect test leaves the test unchanged and divides the slope and interval by Se + Sp - 1", {
  base <- test_trend(c(12, 15, 21, 24, 30), rep(150, 5), time = 0:4)
  res  <- test_trend(c(12, 15, 21, 24, 30), rep(150, 5), time = 0:4,
                     sensitivity = 0.9, specificity = 0.95)
  expect_equal(res$statistic, base$statistic)
  expect_equal(res$slope, base$slope / 0.85)
  expect_equal(c(res$ci_lower, res$ci_upper), c(base$ci_lower, base$ci_upper) / 0.85)
})

test_that("TT-5: rows at the same time are added up; time defaults to 1, 2, ...", {
  # Two half-surveys per round give the same result as one survey per round
  split <- test_trend(c(6, 6, 7, 8, 10, 11, 12, 12, 15, 15), rep(75, 10),
                      time = rep(0:4, each = 2), icc = 0)
  whole <- test_trend(c(12, 15, 21, 24, 30), rep(150, 5), time = 0:4)
  expect_equal(split$slope,   whole$slope)
  expect_equal(split$x_by_time, c(12, 15, 21, 24, 30))
  expect_equal(test_trend(c(12, 15, 21, 24, 30), rep(150, 5))$times, 1:5)
  # Order of the rows does not matter
  expect_equal(test_trend(c(30, 12, 24, 15, 21), rep(150, 5), time = c(4, 0, 3, 1, 2)), whole)
})

test_that("TT-6: a timepoint with 0 positives needs no special handling", {
  res <- test_trend(c(0, 3, 6), c(50, 50, 50), time = 0:2)
  expect_equal(res$slope,     0.06)
  expect_equal(res$statistic, 2.526456, tolerance = 1e-6)
  # (prop.trend.test warns that these three points lie exactly on a line)
  expect_equal(res$p_value, suppressWarnings(stats::prop.trend.test(c(0, 3, 6), c(50, 50, 50), 0:2)$p.value))
  expect_equal(c(res$ci_lower, res$ci_upper), c(0.013530, 0.112342), tolerance = 1e-4)
})

test_that("TT-7: fitted prevalence at the first and last timepoint", {
  # The line through 0.136 at t = 2 with slope 0.03
  res <- test_trend(c(12, 15, 21, 24, 30), rep(150, 5), time = 0:4)
  expect_equal(c(res$prevalence_start_est, res$prevalence_end_est), c(0.076, 0.196))
  expect_equal(res$prevalence_end_est - res$prevalence_start_est, 4 * res$slope)
  # Unequal sizes and uneven times: the ends of R's line weighted by number tested
  x <- c(3, 20, 9); n <- c(40, 200, 60); t <- c(2018, 2021, 2024)
  une <- test_trend(x, n, time = t)
  fit <- stats::lm(I(x / n) ~ t, weights = n)
  expect_equal(c(une$prevalence_start_est, une$prevalence_end_est),
               unname(stats::predict(fit, data.frame(t = c(2018, 2024)))))
})

test_that("TT-8: all negative or all positive everywhere is an error: there is no trend to test", {
  expect_error(test_trend(c(0, 0, 0), rep(50, 3)),
               "Every person tested was negative (`x` is 0 everywhere), so the measured prevalence is 0 at every timepoint and there is no trend to test.",
               fixed = TRUE)
  expect_error(test_trend(c(50, 40, 30), c(50, 40, 30)),
               "Every person tested was positive (`x` equals `n` everywhere), so the measured prevalence is 1 at every timepoint and there is no trend to test.",
               fixed = TRUE)
  # One positive is enough
  expect_no_error(test_trend(c(0, 0, 1), rep(50, 3)))
  expect_no_error(test_trend(c(50, 40, 29), c(50, 40, 30)))
})

test_that("TT-9: interval ends match a scan of slopes, also where the line reaches 0", {
  # Expected ends found by scanning slopes in steps of 1e-6 (see the header)
  # Prevalence falls to 0 and stays there: slopes near the estimate put the
  # line below 0 at the last rounds, where 0 is used
  dec <- test_trend(c(40, 25, 10, 0, 0), rep(100, 5), time = 0:4)
  expect_equal(dec$slope, -0.105)
  expect_equal(c(dec$ci_lower, dec$ci_upper), c(-0.126337, -0.084357), tolerance = 1e-4)
  # Unequal sizes and uneven times
  une <- test_trend(c(3, 20, 9), c(40, 200, 60), time = c(2019, 2021, 2024))
  expect_equal(une$slope, 0.0155)
  expect_equal(c(une$ci_lower, une$ci_upper), c(-0.007228, 0.041219), tolerance = 1e-3)
  expect_false(une$reject)
})

test_that("TT-10: the interval leaves out 0 exactly when the test rejects", {
  set.seed(10)
  for (i in 1:150) {
    k <- sample(2:6, 1); n <- sample(c(1, 2, 5, 20, 100, 500), k, TRUE)
    x <- stats::rbinom(k, n, sample(c(0, 0.02, 0.1, 0.5, 0.9, 1), k, TRUE))
    if (sum(x) == 0 || sum(x) == sum(n)) next
    t <- sort(sample(1:30, k))
    for (alt in c("two.sided", "greater", "less")) for (cl in c(0.8, 0.95)) {
      r <- test_trend(x, n, time = t, alternative = alt, conf_level = cl)
      expect_identical(r$reject, r$ci_lower > 0 || r$ci_upper < 0)
      expect_lte(r$ci_lower, r$ci_upper)
    }
  }
  # A case where the interval's spread rises again between the estimate and
  # 0: the end must be searched for between the estimate and 0, or the
  # interval would contain 0 although the test rejects
  r <- test_trend(c(0, 0, 0, 0, 1), c(2, 2, 2, 50, 1), time = c(3, 6, 7, 8, 10),
                  alternative = "greater")
  expect_true(r$reject)
  expect_gt(r$ci_lower, 0)
  # With sites and fpc_N
  for (i in 1:30) {
    x <- stats::rbinom(8, 30, stats::runif(1, 0.05, 0.4)); tt <- rep(1:4, each = 2)
    r <- suppressWarnings(test_trend(x, rep(30, 8), time = tt, fpc_N = 500))
    expect_identical(r$reject, r$ci_lower > 0 || r$ci_upper < 0)
  }
})


# ---------------------------------------------------------------------------
# Clustering (TT-C-)
# ---------------------------------------------------------------------------

test_that("TT-C-1: a supplied icc gives a design effect from the mean site size", {
  # mean site size 49.75: deff = 1 + 48.75 * 0.05 = 3.4375
  res <- test_trend(c(2, 3, 4, 5, 6, 8, 7, 10), c(50, 45, 55, 50, 48, 52, 47, 51),
                    time = c(1, 1, 2, 2, 3, 3, 4, 4), icc = 0.05)
  expect_equal(res$deff,     3.4375)
  expect_equal(res$icc_used, 0.05)
  expect_equal(res$p_value,  0.117076, tolerance = 1e-5)
  # Same as prop.trend.test on the counts added up within each round, with
  # its chi-square divided by the design effect
  chi <- unname(stats::prop.trend.test(c(5, 9, 14, 17), c(95, 105, 100, 98), 1:4)$statistic)
  expect_equal(res$p_value, stats::pchisq(chi / 3.4375, 1, lower.tail = FALSE))
  expect_equal(c(res$ci_lower, res$ci_upper), c(-0.010399, 0.092616), tolerance = 1e-3)
  # Without clustering the same data give a smaller p-value
  expect_lt(test_trend(c(2, 3, 4, 5, 6, 8, 7, 10), c(50, 45, 55, 50, 48, 52, 47, 51),
                       time = c(1, 1, 2, 2, 3, 3, 4, 4), icc = 0)$p_value, res$p_value)
})

test_that("TT-C-2: with icc = NULL the ICC is estimated within each timepoint and averaged", {
  res <- test_trend(c(1, 9, 2, 10, 3, 12), rep(40, 6), time = c(0, 0, 1, 1, 2, 2))
  expect_equal(res$icc_used, 0.147307, tolerance = 1e-4)
  expect_equal(res$deff,     1 + 39 * res$icc_used)
  # One row per timepoint: nothing to estimate from, no clustering
  one <- test_trend(c(12, 15, 21, 24, 30), rep(150, 5), time = 0:4)
  expect_equal(c(one$icc_used, one$deff), c(0, 1))
  # A timepoint whose sites are all negative, or whose sites have one person
  # each, carries no ICC information and is left out of the average
  ref <- test_trend(c(2, 10, 3, 12), rep(40, 4), time = c(1, 1, 2, 2))$icc_used
  expect_equal(test_trend(c(0, 0, 2, 10, 3, 12), rep(40, 6),
                          time = c(0, 0, 1, 1, 2, 2))$icc_used, ref)
  expect_equal(test_trend(c(0, 1, 2, 10, 3, 12), c(1, 1, 40, 40, 40, 40),
                          time = c(0, 0, 1, 1, 2, 2))$icc_used, ref)
})

test_that("TT-C-3: a supplied icc with one row per timepoint is ignored, with a warning", {
  expect_warning(
    res <- test_trend(c(12, 15, 21, 24, 30), rep(150, 5), time = 0:4, icc = 0.05),
    "`icc` = 0.05 was ignored: no timepoint has two or more sites (or every site has one person), so there is no clustering to adjust for. `icc_used` is reported as 0.",
    fixed = TRUE)
  expect_equal(c(res$icc_used, res$deff), c(0, 1))
})

test_that("TT-C-4: icc = 0 or a tiny icc means no clustering, without a warning", {
  for (icc in c(0, 1e-12)) {
    expect_silent(res <- test_trend(c(1, 2, 3, 4), rep(10, 4), time = c(0, 0, 1, 1), icc = icc))
    expect_equal(c(res$icc_used, res$deff), c(0, 1))
  }
  expect_silent(test_trend(c(12, 15, 21, 24, 30), rep(150, 5), time = 0:4, icc = 1e-12))
})


# ---------------------------------------------------------------------------
# Finite population correction (TT-F-)
# ---------------------------------------------------------------------------

test_that("TT-F-1: fpc_N narrows the interval and strengthens the test", {
  # Each round tests 150 of 600: weights are divided by (600 - 150) / 599
  res  <- test_trend(c(12, 15, 21, 24, 30), rep(150, 5), time = 0:4, fpc_N = 600)
  base <- test_trend(c(12, 15, 21, 24, 30), rep(150, 5), time = 0:4)
  expect_equal(res$statistic, base$statistic / sqrt((600 - 150) / 599))
  expect_equal(res$statistic, 3.910639, tolerance = 1e-6)
  expect_equal(c(res$ci_lower, res$ci_upper), c(0.015013, 0.044597), tolerance = 1e-4)
})

test_that("TT-F-2: fpc_N smaller than, or equal to, the number tested at a timepoint is rejected", {
  expect_error(test_trend(c(12, 15, 21), rep(150, 3), time = 0:2, fpc_N = 100),
               "`fpc_N` (100) is smaller than the number tested at a timepoint (150): you cannot test more people than there are in the population.",
               fixed = TRUE)
  expect_error(test_trend(c(12, 15, 21), rep(150, 3), time = 0:2, fpc_N = 150),
               "`fpc_N` (150) equals the number tested at a timepoint: the whole population was tested there (a census), so that timepoint has no sampling uncertainty. Set `fpc_N = NULL` if no finite-population correction is needed.",
               fixed = TRUE)
  expect_no_error(test_trend(c(12, 15, 21), rep(150, 3), time = 0:2, fpc_N = 151))
})


# ---------------------------------------------------------------------------
# Validation (TT-V-)
# ---------------------------------------------------------------------------

test_that("TT-V-1: x and n are checked as in estimate_prevalence", {
  expect_error(test_trend(NA, 10),
               "`x` and `n` must be numeric, not logical (got class `logical` for x, `numeric` for n). Note: a plain `NA` is logical in R. Remove missing observations before calling.",
               fixed = TRUE)
  expect_error(test_trend("3", 10), "`x` and `n` must be numeric (got class `character` for x, `numeric` for n).", fixed = TRUE)
  expect_error(test_trend(numeric(0), numeric(0)), "`x` and `n` must each contain at least one value.", fixed = TRUE)
  expect_error(test_trend(c(1, NA_real_), c(10, 10)),
               "`x` contains a missing or infinite value (found x[2] = NA).", fixed = TRUE)
  expect_error(test_trend(c(1, 2), c(10, Inf)),
               "`n` contains a missing or infinite value (found n[2] = Inf).", fixed = TRUE)
  expect_error(test_trend(c(1, 2), 10), "`x` and `n` must have the same length (got length(x) = 2, length(n) = 1).", fixed = TRUE)
  expect_error(test_trend(c(-0.5, 2), c(10, 10)), "`x` must be non-negative (found x[1] = -0.5).", fixed = TRUE)
  expect_error(test_trend(c(1, 2), c(10, 0)), "`n` must be positive everywhere (found n[2] = 0).", fixed = TRUE)
  expect_error(test_trend(c(1, 2.5), c(10, 10)), "`x` must contain integers, not decimals (found x[2] = 2.5).", fixed = TRUE)
  expect_error(test_trend(c(1, 2), c(10.5, 10)), "`n` must contain integers, not decimals (found n[1] = 10.5).", fixed = TRUE)
  expect_error(test_trend(c(1, 11), c(10, 10)), "`x` cannot exceed `n` (found x[2] = 11 > n[2] = 10).", fixed = TRUE)
})

test_that("TT-V-2: time must be finite, the same length as x, with at least two different values", {
  expect_error(test_trend(c(3, 4), c(10, 10), time = c(1, NA)),
               "`time` contains a missing or infinite value (found time[2] = NA).", fixed = TRUE)
  expect_error(test_trend(c(3, 4), c(10, 10), time = 1:3),
               "`time` must be a numeric vector the same length as `x` (got class `integer`, length 3, for 2 values of `x`).",
               fixed = TRUE)
  expect_error(test_trend(c(3, 4), c(10, 10), time = c(1, 1)),
               "At least two different timepoints are needed to fit a slope (got 1).", fixed = TRUE)
})

test_that("TT-V-3: single-number arguments: type, NaN/Inf and range", {
  expect_error(test_trend(c(3, 4, 5), rep(10, 3), sensitivity = TRUE),
               "`sensitivity` must be a single number in (0, 1] (got class `logical`, length 1). Note: `TRUE`, `FALSE` and a plain `NA` are logical, not numeric. Pass 1 if the test never misses a case.",
               fixed = TRUE)
  expect_error(test_trend(c(3, 4, 5), rep(10, 3), specificity = NULL),
               "`specificity` must be a single number in (0, 1] (got class `NULL`, length 0).", fixed = TRUE)
  expect_error(test_trend(c(3, 4, 5), rep(10, 3), conf_level = "0.95"),
               "`conf_level` must be a single number in (0, 1) (got class `character`, length 1).", fixed = TRUE)
  expect_error(test_trend(c(3, 4, 5), rep(10, 3), icc = c(0.1, 0.2)),
               "`icc` must be a single number in [0, 1] (got class `numeric`, length 2). To estimate ICC from the data, leave `icc = NULL`.",
               fixed = TRUE)
  expect_error(test_trend(c(3, 4, 5), rep(10, 3), sensitivity = NaN),
               "`sensitivity` must be a single finite number (got NaN).", fixed = TRUE)
  expect_error(test_trend(c(3, 4, 5), rep(10, 3), specificity = -Inf),
               "`specificity` must be a single finite number (got -Inf).", fixed = TRUE)
  expect_error(test_trend(c(3, 4, 5), rep(10, 3), conf_level = Inf),
               "`conf_level` must be a single finite number (got Inf).", fixed = TRUE)
  expect_error(test_trend(c(3, 4, 5), rep(10, 3), icc = NaN),
               "`icc` must be a single finite number (got NaN).", fixed = TRUE)
  expect_error(test_trend(c(3, 4, 5), rep(10, 3), sensitivity = 1.1),
               "`sensitivity` must be in (0, 1] (got 1.1). `sensitivity` is a diagnostic probability and cannot be greater than 1.",
               fixed = TRUE)
  expect_error(test_trend(c(3, 4, 5), rep(10, 3), sensitivity = 0),
               "`sensitivity` must be in (0, 1] (got 0). A sensitivity of 0 means everyone who truly has the condition tests negative.",
               fixed = TRUE)
  expect_error(test_trend(c(3, 4, 5), rep(10, 3), specificity = -0.5),
               "`specificity` must be in (0, 1] (got -0.5). `specificity` is a diagnostic probability and cannot be negative.",
               fixed = TRUE)
  expect_error(test_trend(c(3, 4, 5), rep(10, 3), specificity = 0),
               "`specificity` must be in (0, 1] (got 0). A specificity of 0 means everyone who truly does not have the condition tests positive.",
               fixed = TRUE)
  expect_error(test_trend(c(3, 4, 5), rep(10, 3), icc = 1.5),
               "`icc` must be in [0, 1] (got 1.5). `icc` is a correlation and cannot be greater than 1.", fixed = TRUE)
  for (bad in c(0, 1)) {
    expect_error(test_trend(c(3, 4, 5), rep(10, 3), conf_level = bad),
                 paste0("`conf_level` must be in (0, 1) (got ", bad, "). For example, use 0.95 for a 95% confidence interval."),
                 fixed = TRUE)
  }
  expect_error(test_trend(c(3, 4, 5), rep(10, 3), conf_level = 0.49),
               "`conf_level` must be at least 0.5 (got 0.49). Below 0.5, when prevalence is not changing, the test would wrongly report a trend more than half the time.",
               fixed = TRUE)
  expect_no_error(test_trend(c(3, 4, 5), rep(10, 3), conf_level = 0.5))
})

test_that("TT-V-4: sensitivity + specificity of 1 or less is rejected; below 1.1 warns with the full sum", {
  expect_error(test_trend(c(3, 4, 5), rep(10, 3), sensitivity = 0.5, specificity = 0.5),
               "(got 0.5 + 0.5 = 1). The correction divides by `sensitivity` + `specificity` - 1, so a sum of 1 or less cannot be corrected.",
               fixed = TRUE)
  # Before, a sum just above 1 by rounding gave a slope of about 2e14
  expect_error(test_trend(c(3, 4, 5), rep(10, 3), sensitivity = 0.3,
                          specificity = 0.7 + 2 * .Machine$double.eps),
               "must exceed 1 for the Rogan-Gladen correction", fixed = TRUE)
  expect_warning(test_trend(c(3, 4, 5), rep(10, 3), sensitivity = 0.55, specificity = 0.5),
                 "`sensitivity` + `specificity` = 1.05, which is very close to 1. Correcting for this much test error makes the estimated slope and its confidence interval very wide, and small changes in the assumed sensitivity or specificity will change them a lot.",
                 fixed = TRUE)
})

test_that("TT-V-5: alternative, method and fpc_N are checked, naming the problem", {
  expect_error(test_trend(c(3, 4, 5), rep(10, 3), alternative = "up"),
               "`alternative` must be one of 'two.sided', 'greater', or 'less' (got 'up').", fixed = TRUE)
  expect_error(test_trend(c(3, 4, 5), rep(10, 3), alternative = NA),
               "(got class `logical`, length 1).", fixed = TRUE)
  expect_error(test_trend(c(3, 4, 5), rep(10, 3), method = "wls"),
               "`method` must be 'cochran-armitage', the only method available (got 'wls').", fixed = TRUE)
  expect_error(test_trend(c(3, 4, 5), rep(10, 3), method = 1),
               "(got class `numeric`, length 1).", fixed = TRUE)
  expect_error(test_trend(c(3, 4, 5), rep(10, 3), fpc_N = 100.5),
               "`fpc_N` must be a single finite positive integer (got 100.5).", fixed = TRUE)
  expect_error(test_trend(c(3, 4, 5), rep(10, 3), fpc_N = "100"),
               "(got class `character`). `fpc_N` is the population size per timepoint.", fixed = TRUE)
  expect_error(test_trend(c(3, 4, 5), rep(10, 3), fpc_N = c(100, 200)), "(got length = 2)", fixed = TRUE)
})


# ---------------------------------------------------------------------------
# Return value (TT-R-)
# ---------------------------------------------------------------------------

test_that("TT-R-1: return list contains all expected fields, in order", {
  expect_named(test_trend(c(3, 4, 5), rep(10, 3)),
               c("slope", "slope_app", "statistic", "p_value", "reject", "alternative",
                 "ci_lower", "ci_upper", "prevalence_start_est", "prevalence_end_est",
                 "times", "apparent_prev", "x_by_time", "n_by_time", "n_timepoints",
                 "n_total", "conf_level", "sensitivity", "specificity", "icc_used",
                 "deff", "fpc_N", "method"),
               ignore.order = FALSE)
})


# ---------------------------------------------------------------------------
# Round trips with design_trend() (TT-T-)
# ---------------------------------------------------------------------------

test_that("TT-T-1: data on the planned line at design_trend's n give a significant trend", {
  # design_trend(0.10, 0.20, n_timepoints = 5) -> 161 per round
  n <- design_trend(0.10, 0.20, n_timepoints = 5)$n_per_timepoint
  expect_equal(n, 161)
  p <- seq(0.10, 0.20, length.out = 5)
  expect_true(test_trend(round(p * n), rep(n, 5), time = 0:4)$reject)
})

test_that("TT-T-2: the clustered design from design_trend gives a significant trend when icc is supplied", {
  d <- design_trend(0.10, 0.20, n_timepoints = 5, n_per_site = 20, icc = 0.03)
  expect_equal(c(d$n_sites, d$n_per_site), c(13, 20))
  p <- rep(seq(0.10, 0.20, length.out = 5), each = 13)
  res <- test_trend(round(p * 20), rep(20, 65), time = rep(0:4, each = 13), icc = 0.03)
  expect_equal(res$deff, 1 + 19 * 0.03)
  expect_true(res$reject)
})
