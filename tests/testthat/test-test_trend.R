# Hand-checked reference cases for test_trend()
#
#   p_j   = X_j / N_j, w_j = N_j / deff
#   slope = sum w (t - t_w) p / sum w (t - t_w)^2 (lm(p ~ t, weights = w))
#   z     = slope / sqrt(p_w (1 - p_w) / sum w (t - t_w)^2)  (Cochran-Armitage;
#           the same p-value as prop.trend.test() without sites)
#   interval ends below were found by scanning slopes in steps of 1e-6
#
# 12, 15, 21, 24, 30 of 150 at t = 0 to 4: p = 0.08 ... 0.20, slope = 0.03,
#   p_w = 102 / 750 = 0.136, z = 0.03 / sqrt(0.136 * 0.864 / 1500) = 3.389538,
#   two-sided p = 0.000700, interval [0.012693, 0.046789]

x5 <- c(12, 15, 21, 24, 30)

test_that("five rounds give the hand-checked slope, test and interval", {
  res <- test_trend(x5, rep(150, 5), time = 0:4)
  expect_equal(res$slope, 0.03)
  expect_equal(res$p_value, 2 * pnorm(-0.03 / sqrt(0.136 * 0.864 / 1500)))
  expect_true(res$reject)
  expect_equal(c(res$ci_lower, res$ci_upper), c(0.012693, 0.046789), tolerance = 1e-4)
  expect_equal(c(res$n_eff, res$icc, res$deff), c(750, 0, 1))
  expect_named(res, c("slope", "ci_lower", "ci_upper", "p_value", "reject",
                      "n_eff", "icc", "deff"))
})

test_that("the test is R's prop.trend.test and the slope R's weighted regression", {
  x <- c(12, 15, 21, 24, 30); n <- c(150, 120, 200, 90, 160); t <- c(0, 1, 3, 4, 7)
  res <- test_trend(x, n, time = t)
  expect_equal(res$p_value, stats::prop.trend.test(x, n, t)$p.value)
  expect_equal(res$slope, unname(stats::coef(stats::lm(I(x / n) ~ t, weights = n))[2]))
})

test_that("one-sided tests give one-sided intervals; conf_level sets the level", {
  g <- test_trend(x5, rep(150, 5), time = 0:4, alternative = "greater")
  expect_equal(c(g$ci_lower, g$ci_upper), c(0.015493, Inf), tolerance = 1e-4)
  # the one-sided end is the lower end of the two-sided 90 percent interval
  expect_equal(g$ci_lower, test_trend(x5, rep(150, 5), time = 0:4, conf_level = 0.9)$ci_lower)
  l <- test_trend(x5, rep(150, 5), time = 0:4, alternative = "less")
  expect_equal(c(l$ci_lower, l$ci_upper), c(-Inf, 0.044142), tolerance = 1e-4)
  expect_false(l$reject)
  # a weak trend: not significant at 0.05, significant at 0.20
  weak <- c(12, 14, 17, 17, 20)
  expect_false(test_trend(weak, rep(150, 5), time = 0:4)$reject)
  expect_true(test_trend(weak, rep(150, 5), time = 0:4, conf_level = 0.8)$reject)
})

test_that("an imperfect test divides the slope and interval by Se + Sp - 1", {
  base <- test_trend(x5, rep(150, 5), time = 0:4)
  res  <- test_trend(x5, rep(150, 5), time = 0:4, sensitivity = 0.9, specificity = 0.95)
  expect_equal(res$p_value, base$p_value)
  expect_equal(c(res$slope, res$ci_lower, res$ci_upper),
               c(base$slope, base$ci_lower, base$ci_upper) / 0.85)
})

test_that("rows at the same time are added up; row order does not matter", {
  whole <- test_trend(x5, rep(150, 5), time = 0:4)
  split <- test_trend(c(6, 6, 7, 8, 10, 11, 12, 12, 15, 15), rep(75, 10),
                      time = rep(0:4, each = 2), icc = 0)
  expect_equal(split$slope, whole$slope)
  expect_equal(test_trend(c(30, 12, 24, 15, 21), rep(150, 5), time = c(4, 0, 3, 1, 2)), whole)
  # time defaults to 1, 2, ...
  expect_equal(test_trend(x5, rep(150, 5)), test_trend(x5, rep(150, 5), time = 1:5))
})

test_that("interval ends match a scan of slopes, also where the line reaches 0", {
  dec <- test_trend(c(40, 25, 10, 0, 0), rep(100, 5), time = 0:4)
  expect_equal(c(dec$slope, dec$ci_lower, dec$ci_upper), c(-0.105, -0.126337, -0.084357),
               tolerance = 1e-4)
  une <- test_trend(c(3, 20, 9), c(40, 200, 60), time = c(2019, 2021, 2024))
  expect_equal(c(une$slope, une$ci_lower, une$ci_upper), c(0.0155, -0.007228, 0.041219),
               tolerance = 1e-3)
  # a round with 0 positives needs no special handling
  expect_equal(test_trend(c(0, 3, 6), rep(50, 3), time = 0:2)$slope, 0.06)
  # 0 then 100 percent: the line through the two points is 0 and 1 there,
  # so its variance is 0 and the interval is the single slope 1
  res <- test_trend(c(0, 50), c(50, 50), time = 0:1)
  expect_equal(c(res$slope, res$ci_lower, res$ci_upper), c(1, 1, 1))
})

test_that("the interval leaves out 0 exactly when the test rejects", {
  set.seed(1)
  for (i in 1:200) {
    n <- sample(20:200, 4); x <- rbinom(4, n, runif(1, 0.05, 0.4))
    if (sum(x) == 0) next
    alt <- sample(c("two.sided", "greater", "less"), 1)
    r <- test_trend(x, n, time = 1:4, alternative = alt)
    expect_equal(r$reject, r$ci_lower > 0 || r$ci_upper < 0)
  }
})

test_that("a supplied icc gives a design effect from the size-weighted site size", {
  x <- c(2, 3, 4, 5, 6, 8, 7, 10); n <- c(50, 45, 55, 50, 48, 52, 47, 51)
  t <- c(1, 1, 2, 2, 3, 3, 4, 4)
  res  <- test_trend(x, n, time = t, icc = 0.05)
  deff <- 1 + (sum(n^2) / sum(n) - 1) * 0.05
  expect_equal(c(res$icc, res$deff, res$n_eff), c(0.05, deff, sum(n) / deff))
  # only timepoints with two or more sites count: rounds 2 and 3 have one
  # survey of 500, so n_w is that of the two sites of 50 in round 1
  mixed <- test_trend(c(3, 5, 60, 90), c(50, 50, 500, 500), time = c(1, 1, 2, 3),
                      icc = 0.05)
  expect_equal(mixed$deff, 1 + 49 * 0.05)
  # prop.trend.test on the counts added up within each round, with its
  # chi-square divided by the design effect
  chi <- unname(stats::prop.trend.test(c(5, 9, 14, 17), c(95, 105, 100, 98), 1:4)$statistic)
  expect_equal(res$p_value, stats::pchisq(chi / deff, 1, lower.tail = FALSE))
  # very unequal sites: n_w = (10^2 + 190^2) / 200 = 181, not the mean 100
  r2 <- test_trend(c(1, 20, 2, 25), c(10, 190, 10, 190), time = c(1, 1, 2, 2), icc = 0.1)
  expect_equal(r2$deff, 1 + 180 * 0.1)
})

test_that("icc = NULL: estimated within each timepoint and averaged", {
  res <- test_trend(c(1, 9, 2, 10, 3, 12), rep(40, 6), time = c(0, 0, 1, 1, 2, 2))
  expect_equal(res$icc, 0.147307, tolerance = 1e-4)
  expect_equal(res$deff, 1 + 39 * res$icc)
  # timepoints that are all negative, or have one person per site, are skipped
  ref <- test_trend(c(2, 10, 3, 12), rep(40, 4), time = c(1, 1, 2, 2))$icc
  expect_equal(test_trend(c(0, 0, 2, 10, 3, 12), rep(40, 6),
                          time = c(0, 0, 1, 1, 2, 2))$icc, ref)
  expect_equal(test_trend(c(0, 1, 2, 10, 3, 12), c(1, 1, 40, 40, 40, 40),
                          time = c(0, 0, 1, 1, 2, 2))$icc, ref)
})

test_that("one site per timepoint: no clustering, a supplied icc is ignored with a warning", {
  expect_equal(test_trend(x5, rep(150, 5), time = 0:4)$deff, 1)
  expect_warning(res <- test_trend(x5, rep(150, 5), time = 0:4, icc = 0.05),
                 "`icc` was ignored")
  expect_equal(c(res$icc, res$deff), c(0, 1))
  expect_no_warning(test_trend(x5, rep(150, 5), time = 0:4, icc = 0))
})

test_that("design_trend's design gives a significant trend on data on the planned line", {
  n <- design_trend(0.10, 0.20, n_timepoints = 5)$n_per_timepoint
  p <- seq(0.10, 0.20, length.out = 5)
  expect_true(test_trend(round(p * n), rep(n, 5), time = 0:4)$reject)
  d <- design_trend(0.10, 0.20, n_timepoints = 5, n_per_site = 20, icc = 0.03)
  p <- rep(p, each = d$n_sites)
  res <- test_trend(round(p * 20), rep(20, 5 * d$n_sites), time = rep(0:4, each = d$n_sites),
                    icc = 0.03)
  expect_equal(res$deff, d$deff)
  expect_true(res$reject)
})

test_that("bad inputs give a clear error", {
  expect_error(test_trend(c(1, 2), c(10)), "same length")
  expect_error(test_trend(c(1, 20), c(10, 10)), "between 0 and `n`")
  expect_error(test_trend(c(1, 2), c(10, 10), time = c(1, NA)), "`time` must be")
  expect_error(test_trend(c(1, 2), c(10, 10), time = c(1, 1)), "two different timepoints")
  expect_error(test_trend(c(0, 0, 0), rep(50, 3)), "negative, so there is no trend")
  expect_error(test_trend(c(50, 40), c(50, 40)), "positive, so there is no trend")
  expect_error(test_trend(x5, rep(150, 5), sensitivity = 0.5, specificity = 0.5),
               "greater than 1")
  expect_error(test_trend(x5, rep(150, 5), conf_level = 1), "`conf_level` must be")
  expect_error(test_trend(x5, rep(150, 5), icc = 2), "`icc` must be")
  expect_error(test_trend(x5, rep(150, 5), alternative = "up"), "should be one of")
})
