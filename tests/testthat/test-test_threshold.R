# Reference cases for test_threshold() (exact binomial test)
#
# With one site, no fpc_N and sensitivity = specificity = 1 the one-sided
# p-values are those of binom.test(), and the Clopper-Pearson interval is
# binom.test()'s interval. Numbers below were checked with binom.test():
#   8 of 100 above 0.05:  p = 0.127960478620 -> not rejected at 0.05
#   15 of 100 above 0.05: p = 0.000135854239 -> rejected
#   8 of 100, threshold 0.08: greater p = 0.5528902488, less p = 0.5926282257
#   0 of 100 below 0.05:  p = 0.95^100 = 0.00592052922 -> rejected
#   30 of 100: two-sided interval 0.2124064205 to 0.3998146762;
#              one-sided lower 0.2249232244, one-sided upper 0.3842206128
#
# Clustered (10 sites of 10, 3 positives each, icc 0.05, threshold 0.2):
#   deff = 1.45, effective size m = 100 / 1.45 = 68.96551724138,
#   effective positives 0.3 * m = 20.68965517241,
#   p = pbeta(0.2, 20.689..., 68.965... - 20.689... + 1) = 0.03139052718
#   (without clustering: binom.test(30, 100, 0.2, "greater") = 0.01124897872)


# ---------------------------------------------------------------------------
# Core scenarios (TH-)
# ---------------------------------------------------------------------------

test_that("TH-1: 8 of 100 is not significantly above 0.05 (p = 0.128, as binom.test)", {
  res <- test_threshold(x = 8, n = 100, threshold = 0.05)
  expect_equal(res$threshold,   0.05)
  expect_equal(res$alternative, "greater")
  expect_equal(res$prevalence,  0.08)
  expect_equal(res$n_total,     100)
  expect_equal(res$deff,        1)
  expect_equal(res$statistic,   8)   # effective positives = x with one site
  expect_equal(res$p_value, 0.127960478620, tolerance = 1e-10)
  expect_equal(res$p_value,
               binom.test(8, 100, 0.05, alternative = "greater")$p.value,
               tolerance = 1e-12)
  expect_false(res$reject)
})

test_that("TH-2: 15 of 100 is significantly above 0.05; one-sided p-values are binom.test's", {
  res <- test_threshold(x = 15, n = 100, threshold = 0.05)
  expect_equal(res$p_value, 0.000135854239, tolerance = 1e-10)
  expect_true(res$reject)
  for (alt in c("greater", "less")) {
    expect_equal(test_threshold(x = 12, n = 150, threshold = 0.10, alternative = alt)$p_value,
                 binom.test(12, 150, 0.10, alternative = alt)$p.value,
                 tolerance = 1e-12)
  }
  # The same share in ten times as many people gives a smaller p-value
  big <- test_threshold(x = 80, n = 1000, threshold = 0.05)
  expect_lt(big$p_value, test_threshold(x = 8, n = 100, threshold = 0.05)$p_value)
  expect_equal(big$p_value, binom.test(80, 1000, 0.05, alternative = "greater")$p.value,
               tolerance = 1e-12)
})

test_that("TH-3: a measured prevalence exactly at the threshold gives p-values above 0.5 both ways", {
  gt <- test_threshold(x = 8, n = 100, threshold = 0.08)
  lt <- test_threshold(x = 8, n = 100, threshold = 0.08, alternative = "less")
  ts <- test_threshold(x = 8, n = 100, threshold = 0.08, alternative = "two.sided")
  expect_equal(gt$p_value, 0.5528902488, tolerance = 1e-9)
  expect_equal(lt$p_value, 0.5926282257, tolerance = 1e-9)
  # two-sided = twice the smaller one-sided p-value, at most 1
  expect_equal(ts$p_value, 1)
  expect_false(gt$reject); expect_false(lt$reject); expect_false(ts$reject)
  expect_equal(test_threshold(x = 10, n = 100, threshold = 0.10)$p_value,
               0.5487098346, tolerance = 1e-9)
})

test_that("TH-4: an imperfect test moves the threshold to the measured scale and corrects the prevalence", {
  res <- test_threshold(x = 8, n = 100, threshold = 0.05,
                        sensitivity = 0.9, specificity = 0.95)
  # threshold on the measured scale = 0.05 * 0.9 + 0.95 * 0.05 = 0.0925
  expect_equal(res$p_value,
               binom.test(8, 100, 0.0925, alternative = "greater")$p.value,
               tolerance = 1e-12)
  expect_equal(res$p_value, 0.7166645524, tolerance = 1e-9)
  # Rogan-Gladen: (0.3 - (1 - 0.95)) / (0.9 + 0.95 - 1)
  expect_equal(test_threshold(x = 30, n = 100, threshold = 0.2,
                              sensitivity = 0.9, specificity = 0.95)$prevalence,
               (0.30 - 0.05) / 0.85)
  # 20 of 100 with sensitivity 0.8, specificity 0.9: (0.2 - 0.1) / 0.7, below 0.2
  expect_equal(test_threshold(x = 20, n = 100, threshold = 0.10,
                              sensitivity = 0.80, specificity = 0.90)$prevalence,
               (0.20 - 0.10) / 0.70)
})

test_that("TH-5: alternative = 'less' tests for a prevalence below the threshold", {
  # 0 of 100 below 0.05: p = 0.95^100 = 0.0059
  res_less <- test_threshold(x = 0, n = 100, threshold = 0.05, alternative = "less")
  expect_equal(res_less$p_value, 0.95^100, tolerance = 1e-12)
  expect_true(res_less$reject)
  # 0 of 100 below 0.20: p = 0.8^100
  expect_equal(test_threshold(x = 0, n = 100, threshold = 0.20, alternative = "less")$p_value,
               0.8^100, tolerance = 1e-12)
  # With 0 positives the 'greater' p-value is 1 (0 or more always happens)
  res_greater <- test_threshold(x = 0, n = 100, threshold = 0.05)
  expect_equal(res_greater$p_value, 1)
  expect_false(res_greater$reject)
})

test_that("TH-6: the two-sided p-value is twice the smaller one-sided p-value, at most 1", {
  for (x in c(12, 15)) {
    expect_equal(test_threshold(x = x, n = 100, threshold = 0.05, alternative = "two.sided")$p_value,
                 2 * test_threshold(x = x, n = 100, threshold = 0.05)$p_value,
                 tolerance = 1e-12)
  }
  # This differs from binom.test()'s two-sided rule: 8 of 100 at 0.05
  expect_equal(test_threshold(x = 8, n = 100, threshold = 0.05, alternative = "two.sided")$p_value,
               0.2559209572, tolerance = 1e-9)
})

test_that("TH-7: the two-sided interval is binom.test()'s Clopper-Pearson interval, whatever the threshold", {
  res <- test_threshold(x = 30, n = 100, threshold = 0.20, alternative = "two.sided")
  expect_equal(res$prevalence, 0.3)
  expect_equal(c(res$ci_lower, res$ci_upper), c(0.2124064205, 0.3998146762), tolerance = 1e-9)
  expect_equal(c(res$ci_lower, res$ci_upper),
               as.numeric(binom.test(30, 100)$conf.int), tolerance = 1e-12)
  res2 <- test_threshold(x = 30, n = 100, threshold = 0.40, alternative = "two.sided")
  expect_equal(c(res2$ci_lower, res2$ci_upper), c(res$ci_lower, res$ci_upper))
})

test_that("TH-8: the interval is matched to alternative and reject flips at its ends", {
  gt <- test_threshold(x = 30, n = 100, threshold = 0.20, alternative = "greater")
  lt <- test_threshold(x = 30, n = 100, threshold = 0.20, alternative = "less")
  expect_equal(gt$ci_upper, 1)
  expect_equal(lt$ci_lower, 0)
  expect_equal(gt$ci_lower, 0.2249232244, tolerance = 1e-9)
  expect_equal(lt$ci_upper, 0.3842206128, tolerance = 1e-9)
  # Just below and just above each end
  L <- gt$ci_lower
  expect_true(test_threshold(x = 30, n = 100, threshold = L - 1e-6)$reject)
  expect_false(test_threshold(x = 30, n = 100, threshold = L + 1e-6)$reject)
  U <- lt$ci_upper
  expect_true(test_threshold(x = 30, n = 100, threshold = U + 1e-6, alternative = "less")$reject)
  expect_false(test_threshold(x = 30, n = 100, threshold = U - 1e-6, alternative = "less")$reject)
})

test_that("TH-9: reject agrees with the interval in random designs (sites, fpc_N, imperfect tests)", {
  set.seed(11)
  for (i in 1:300) {
    k   <- sample(1:6, 1)
    n   <- sample(5:150, k, replace = TRUE)
    x   <- vapply(n, function(m) rbinom(1, m, runif(1)), numeric(1))
    th  <- runif(1, 0.01, 0.99)
    alt <- c("greater", "less", "two.sided")[i %% 3 + 1]
    fpc <- if (i %% 4 == 0) sum(n) + sample(1:500, 1) else NULL
    r <- test_threshold(x, n, th, alternative = alt, sensitivity = 0.9,
                        specificity = 0.97, icc = if (i %% 2) NULL else 0.05,
                        fpc_N = fpc)
    expect_identical(r$reject, th < r$ci_lower || th > r$ci_upper)
  }
})

test_that("TH-10: conf_level sets the significance level and the interval width", {
  # A p-value equal to 1 - conf_level is not significant. 1 of 1 against a
  # threshold of 0.05 gives p = 0.05 (stored as 0.050000000000000010), and
  # 1 - 0.95 is stored as 0.050000000000000044: without rounding the
  # significance level, this p-value would count as below it. The
  # interval's lower end is 0.05, the threshold itself.
  res <- test_threshold(x = 1, n = 1, threshold = 0.05)
  expect_equal(res$p_value, 0.05)
  expect_equal(res$p_value, binom.test(1, 1, 0.05, alternative = "greater")$p.value)
  expect_false(res$reject)
  expect_equal(res$ci_lower, 0.05)
  # p-values exactly equal to the significance level, stored exactly:
  # 1 of 1 and 0 of 1 against 0.5 give p = 0.5; 0 of 2 against 0.5 gives
  # p = 0.25. None is significant.
  res <- test_threshold(x = 1, n = 1, threshold = 0.5, conf_level = 0.5)
  expect_identical(res$p_value, 0.5)
  expect_false(res$reject)
  res <- test_threshold(x = 0, n = 1, threshold = 0.5, alternative = "less", conf_level = 0.5)
  expect_identical(res$p_value, 0.5)
  expect_false(res$reject)
  res <- test_threshold(x = 0, n = 2, threshold = 0.5, alternative = "less", conf_level = 0.75)
  expect_identical(res$p_value, 0.25)
  expect_false(res$reject)
  # 9 of 100 above 0.05: exact p = 0.0631, significant at 0.10 only
  expect_equal(binom.test(9, 100, 0.05, alternative = "greater")$p.value, 0.0630896, tolerance = 1e-6)
  expect_true(test_threshold(x = 9, n = 100, threshold = 0.05, conf_level = 0.90)$reject)
  expect_false(test_threshold(x = 9, n = 100, threshold = 0.05)$reject)
  expect_false(test_threshold(x = 9, n = 100, threshold = 0.05, conf_level = 0.99)$reject)
  # conf_level 0.5: 8 of 100 (p = 0.128) is significant at 0.5
  expect_true(test_threshold(x = 8, n = 100, threshold = 0.05, conf_level = 0.5)$reject)
  # conf_level 0.9999: a wide interval
  res <- test_threshold(x = 8, n = 100, threshold = 0.05, conf_level = 0.9999)
  expect_lt(res$ci_lower, 0.02)
  expect_gt(res$ci_upper, 0.12)
})

test_that("TH-11: no positives, or everyone positive, at one or many sites", {
  zero <- test_threshold(x = 0, n = 100, threshold = 0.05)
  expect_equal(c(zero$prevalence, zero$statistic, zero$p_value), c(0, 0, 1))
  expect_false(zero$reject)
  zeros <- test_threshold(x = rep(0, 5), n = rep(20, 5), threshold = 0.05)
  expect_equal(c(zeros$prevalence, zeros$statistic, zeros$p_value), c(0, 0, 1))
  all1 <- test_threshold(x = 100, n = 100, threshold = 0.05)
  expect_equal(all1$prevalence, 1)
  expect_equal(all1$p_value, 0.05^100, tolerance = 1e-12)
  expect_true(all1$reject)
  alls <- test_threshold(x = rep(10, 5), n = rep(10, 5), threshold = 0.50)
  expect_equal(alls$prevalence, 1)
  expect_true(alls$reject)
  # One person tested
  expect_equal(test_threshold(x = 0, n = 1, threshold = 0.3)$p_value, 1)
  expect_equal(test_threshold(x = 1, n = 1, threshold = 0.3)$p_value, 0.3)
})

test_that("TH-12: thresholds near 0 or 1 give valid results", {
  near1 <- test_threshold(x = 80, n = 100, threshold = 0.999)
  expect_equal(near1$p_value, binom.test(80, 100, 0.999, alternative = "greater")$p.value,
               tolerance = 1e-12)
  near0 <- test_threshold(x = 1, n = 100, threshold = 0.001)
  expect_equal(near0$p_value, binom.test(1, 100, 0.001, alternative = "greater")$p.value,
               tolerance = 1e-12)
  # 100 of 100 above 0.999: p = 0.999^100 = 0.905
  expect_equal(test_threshold(x = 100, n = 100, threshold = 0.999)$p_value, 0.999^100,
               tolerance = 1e-12)
})

test_that("TH-13: a corrected prevalence or interval end below 0 or above 1 is set to 0 or 1", {
  # Sensitivity and specificity 0.9: 1 of 100 corrects to (0.01 - 0.1) / 0.8 < 0
  expect_equal(test_threshold(x = 1, n = 100, threshold = 0.30,
                              sensitivity = 0.9, specificity = 0.9)$prevalence, 0)
  # 99 of 100 corrects to (0.99 - 0.1) / 0.8 > 1
  expect_equal(test_threshold(x = 99, n = 100, threshold = 0.50,
                              sensitivity = 0.9, specificity = 0.9)$prevalence, 1)
  # The interval ends are set to 0 or 1 in the same way, for every ci_method.
  for (m in c("clopper-pearson", "agresti-coull", "wald")) {
    # "less" with specificity 0.95: the lower end 0 corrects to
    # (0 - 0.05) / 0.95 < 0, and is returned as 0
    res <- test_threshold(8, 100, 0.2, alternative = "less", specificity = 0.95, ci_method = m)
    expect_identical(res$ci_lower, 0)
    expect_gt(res$ci_upper, 0)
    # 0 of 100 with specificity 0.9: the upper end is below 0.1 on the
    # measured scale (Clopper-Pearson: 0.0295), which corrects to below 0.
    # Both ends are returned as 0.
    res <- suppressWarnings(
      test_threshold(0, 100, 0.2, alternative = "less", specificity = 0.9, ci_method = m))
    expect_identical(c(res$ci_lower, res$ci_upper), c(0, 0))
    # 100 of 100 with sensitivity 0.9: the lower end is above 0.9 on the
    # measured scale (Clopper-Pearson: 0.9705), which corrects to above 1.
    # Both ends are returned as 1.
    res <- suppressWarnings(
      test_threshold(100, 100, 0.5, sensitivity = 0.9, ci_method = m))
    expect_identical(c(res$ci_lower, res$ci_upper), c(1, 1))
  }
  # An end inside 0 to 1 is only corrected: 8 of 100, "less", specificity
  # 0.95: upper end 0.1397171 on the measured scale,
  # (0.1397171 - 0.05) / 0.95 = 0.09443907
  expect_equal(test_threshold(8, 100, 0.2, alternative = "less", specificity = 0.95)$ci_upper,
               0.09443907, tolerance = 1e-7)
})

test_that("TH-14: p-values are between 0 and 1 and the prevalence is inside its interval", {
  for (alt in c("greater", "less", "two.sided")) {
    for (x in c(0, 1, 5, 50, 95, 99, 100)) {
      res <- test_threshold(x = x, n = 100, threshold = 0.30, alternative = alt)
      expect_gte(res$p_value, 0)
      expect_lte(res$p_value, 1)
      expect_lte(res$ci_lower, res$prevalence + 1e-10)
      expect_gte(res$ci_upper, res$prevalence - 1e-10)
    }
  }
})

test_that("TH-15: counts stored as R integers give the same result as the same numbers", {
  expect_equal(test_threshold(8L, 100L, 0.05), test_threshold(8, 100, 0.05))
  expect_equal(test_threshold(c(5L, 3L, 8L), c(40L, 50L, 60L), 0.05, icc = 0.05),
               test_threshold(c(5, 3, 8), c(40, 50, 60), 0.05, icc = 0.05))
})

test_that("TH-16: counts adding up past R's integer limit (about 2.1 billion) still work", {
  # Two sites of 2e9 people: the total of 4e9 does not fit in an R integer
  res <- expect_no_warning(test_threshold(c(1e7L, 1e7L), c(2e9L, 2e9L), 0.005))
  expect_equal(res$n_total, 4e9)
  expect_equal(res, test_threshold(c(1e7, 1e7), c(2e9, 2e9), 0.005))
  # A measured prevalence of exactly the threshold: p-value close to 0.5
  expect_equal(res$p_value, 0.5, tolerance = 1e-3)
})

test_that("TH-17: sensitivity = specificity = 0.999 gives almost the same result as 1", {
  one  <- test_threshold(15, 100, 0.05)
  near <- test_threshold(15, 100, 0.05, sensitivity = 0.999, specificity = 0.999)
  # Within 0.002 on the prevalence scale
  expect_lt(abs(near$prevalence - one$prevalence), 0.002)
  expect_lt(abs(near$ci_lower - one$ci_lower), 0.002)
  expect_equal(near$ci_upper, one$ci_upper)
  expect_equal(near$reject, one$reject)
  expect_lt(near$p_value, 1e-3)
})


# ---------------------------------------------------------------------------
# Clustering (TH-C-)
# ---------------------------------------------------------------------------

test_that("TH-C-1: clustering (icc 0.05) shrinks the effective size and raises the p-value", {
  res_none <- test_threshold(x = rep(3, 10), n = rep(10, 10), threshold = 0.2, icc = 0)
  res      <- test_threshold(x = rep(3, 10), n = rep(10, 10), threshold = 0.2, icc = 0.05)
  expect_equal(res_none$p_value, 0.01124897872, tolerance = 1e-9)
  expect_equal(res$deff, 1.45)
  expect_equal(res$n_eff_adj, 100 / 1.45)
  expect_equal(res$statistic, 0.3 * 100 / 1.45)
  expect_equal(res$p_value, 0.03139052718, tolerance = 1e-9)
})

test_that("TH-C-2: a supplied icc gives deff = 1 + (mean site size - 1) * icc, up to icc = 1", {
  res <- test_threshold(x = c(3, 7, 5), n = c(30, 30, 30), threshold = 0.15, icc = 0.3)
  expect_equal(res$deff, 1 + 29 * 0.3)
  expect_equal(res$icc_used, 0.3)
  # icc = 1: each site counts as one person, deff = mean site size
  one <- test_threshold(x = c(5, 5), n = c(10, 10), threshold = 0.3, icc = 1)
  expect_equal(c(one$deff, one$n_eff), c(10, 2))
})

test_that("TH-C-3: with icc = NULL the ICC is estimated from the sites", {
  # Identical sites: no extra variation, deff 1
  same <- test_threshold(x = c(5, 5, 5, 5), n = rep(50, 4), threshold = 0.05)
  expect_equal(same$deff, 1)
  # Very different sites: deff above 1, smaller effective size
  diff <- test_threshold(x = c(0, 0, 20, 20), n = rep(50, 4), threshold = 0.05)
  expect_gt(diff$deff, 1)
  expect_lt(diff$n_eff, same$n_eff)
  # Half the sites all negative, half all positive: ICC kept at 1, deff = 10
  split <- test_threshold(x = c(rep(0, 5), rep(10, 5)), n = rep(10, 10), threshold = 0.30)
  expect_equal(c(split$icc_used, split$deff), c(1, 10))
  # Two sites are enough to estimate the ICC. 2 of 100 and 20 of 100:
  # measured prevalence 0.11, site prevalences 0.02 and 0.20.
  # Variance between the sites = (0.09^2 + 0.09^2) / 1 = 0.0162.
  # Variance expected if everyone were independent = 0.11 * 0.89 / 100 = 0.000979.
  # Ratio 16.5475, ICC = (16.5475 - 1) / (100 - 1) = 0.1570454,
  # deff = 1 + 99 * 0.1570454 = 16.5475, n_eff = 200 / 16.5475 = 12.08642
  two <- test_threshold(x = c(2, 20), n = c(100, 100), threshold = 0.05)
  expect_equal(two$icc_used, 0.1570454, tolerance = 1e-6)
  expect_equal(two$deff, 16.5475, tolerance = 1e-6)
  expect_equal(two$n_eff, 12.08642, tolerance = 1e-6)
  expect_equal(two$p_value, 0.3130447, tolerance = 1e-6)
  # Unequal site sizes: 21 of 160, measured prevalence 0.13125, site
  # prevalences 0.05, 0.125, 0, 0.2, 0.12 (mean 0.099).
  # Variance between the sites = 0.00588.
  # Variance expected if everyone were independent = 0.13125 * 0.86875 *
  # mean(1 / n) = 0.1140234 * 0.03966667 = 0.004522930.
  # Ratio 1.300042. The ICC divides by the plain mean site size minus 1:
  # (1.300042 - 1) / (32 - 1) = 0.009678784.
  # deff then uses the size-weighted mean site size, 6450 / 160 = 40.3125:
  # 1 + 39.3125 * 0.009678784 = 1.380497
  n <- c(20, 40, 15, 60, 25)
  une <- test_threshold(x = c(1, 5, 0, 12, 3), n = n, threshold = 0.10)
  expect_equal(une$icc_used, 0.009678784, tolerance = 1e-6)
  expect_equal(une$deff, 1.380497, tolerance = 1e-6)
  expect_equal(une$deff, 1 + (sum(n^2) / sum(n) - 1) * une$icc_used)
  # A supplied ICC: the same rule, so unequal sizes give a larger deff than
  # the plain mean would
  sup <- test_threshold(x = c(1, 5, 0, 12, 3), n = n, threshold = 0.10, icc = 0.1)
  expect_equal(sup$deff, 1 + (sum(n^2) / sum(n) - 1) * 0.1)
  expect_gt(sup$deff, 1 + (mean(n) - 1) * 0.1)
})

test_that("TH-C-4: icc = 0 means no clustering, even with very different sites", {
  forced <- test_threshold(x = c(0, 20, 0, 20), n = rep(20, 4), threshold = 0.30, icc = 0)
  est    <- test_threshold(x = c(0, 20, 0, 20), n = rep(20, 4), threshold = 0.30)
  expect_equal(forced$deff, 1)
  expect_gt(est$deff, 1)
  expect_lt(forced$p_value, est$p_value)
  # With one site, icc = 0 and icc = NULL agree
  expect_equal(test_threshold(x = 30, n = 100, threshold = 0.2, icc = 0),
               test_threshold(x = 30, n = 100, threshold = 0.2, icc = NULL))
})

test_that("TH-C-5: a supplied icc is used with one site; one person per site gives deff 1", {
  # With icc = NULL the ICC cannot be estimated from one site, or from one
  # person per site: deff 1, icc_used 0, no warning
  for (args in list(list(x = 15, n = 200), list(x = c(1, 0, 1, 0, 1), n = rep(1, 5)))) {
    expect_silent(res <- test_threshold(args$x, args$n, threshold = 0.05))
    expect_equal(c(res$icc_used, res$deff), c(0, 1))
    expect_equal(res$n_eff, sum(args$n))
  }
  # One site of 50 with icc 0.05: deff = 1 + (50 - 1) * 0.05 = 3.45, so the
  # 50 people count as m = 50 / 3.45 = 14.49 and the 8 positives as
  # 0.16 * m = 2.32. No warning.
  expect_silent(res <- test_threshold(x = 8, n = 50, threshold = 0.05, icc = 0.05))
  expect_identical(res$icc_used, 0.05)
  expect_equal(res$deff, 3.45)
  expect_equal(res$n_eff, 50 / 3.45)
  # It changes the answer: p = 0.1015 (0.0032 without icc)
  m <- 50 / 3.45
  expect_equal(res$p_value, pbeta(0.05, 0.16 * m, m - 0.16 * m + 1))
  expect_equal(res$p_value, 0.1015206, tolerance = 1e-6)
  expect_equal(test_threshold(x = 8, n = 50, threshold = 0.05)$p_value, 0.003188343, tolerance = 1e-6)
  # One person per site: deff = 1 + (1 - 1) * 0.2 = 1, icc reported as given
  expect_silent(res <- test_threshold(x = c(0, 1, 1), n = c(1, 1, 1), threshold = 0.3, icc = 0.2))
  expect_identical(c(res$icc_used, res$deff), c(0.2, 1))
  # A tiny icc is used as given: 1 + 49 * 1e-12
  expect_silent(res <- test_threshold(x = 8, n = 50, threshold = 0.05, icc = 1e-12))
  expect_identical(res$deff, 1 + 49 * 1e-12)
})

test_that("TH-C-6: the order of the sites and their names do not change the result", {
  x <- c(2, 9, 4); n <- c(30, 50, 40)
  # icc estimated from the sites, and icc supplied
  for (icc in list(NULL, 0.1)) {
    a <- test_threshold(x, n, 0.05, icc = icc)
    b <- test_threshold(rev(x), rev(n), 0.05, icc = icc)
    expect_equal(a, b)
    nm <- test_threshold(c(first = 2, second = 9, third = 4), c(a = 30, b = 50, c = 40), 0.05, icc = icc)
    expect_equal(a, nm)
  }
})


# ---------------------------------------------------------------------------
# Finite population correction (TH-F-)
# ---------------------------------------------------------------------------

test_that("TH-F-1: fpc_N raises the effective size and lowers the p-value", {
  res_none <- test_threshold(x = 8, n = 100, threshold = 0.05)
  res      <- test_threshold(x = 8, n = 100, threshold = 0.05, fpc_N = 200)
  # effective size 100 * 199 / 100 = 199, effective positives 0.08 * 199
  expect_equal(res$n_eff, 100)
  expect_equal(res$n_eff_adj, 199)
  expect_equal(res$p_value, 0.04480149153, tolerance = 1e-9)
  expect_lt(res$p_value, res_none$p_value)
  expect_equal(res_none$n_eff_adj, res_none$n_eff)
})

test_that("TH-F-2: fpc_N smaller than, or equal to, the number tested is rejected", {
  expect_error(test_threshold(x = 8, n = 100, threshold = 0.05, fpc_N = 50),
               "`fpc_N` (50) is smaller than the total sample size (100): you cannot test more people than there are in the population.",
               fixed = TRUE)
  # One fewer than the number tested
  expect_error(test_threshold(x = 8, n = 100, threshold = 0.05, fpc_N = 99),
               "`fpc_N` (99) is smaller than the total sample size (100): you cannot test more people than there are in the population.",
               fixed = TRUE)
  expect_error(test_threshold(x = c(5, 5), n = c(50, 50), threshold = 0.05, icc = 0.5, fpc_N = 100),
               "`fpc_N` (100) equals the total sample size: the whole population was tested (a census), so there is no sampling uncertainty and no confidence interval to compute. Set `fpc_N = NULL` if no finite-population correction is needed.",
               fixed = TRUE)
  # One more than the number tested works
  res <- test_threshold(x = 8, n = 100, threshold = 0.05, fpc_N = 101)
  expect_equal(res$n_eff_adj, 100 * 100 / 1)
  expect_true(is.finite(res$p_value))
})

test_that("TH-F-3: a very large population makes the correction negligible", {
  x <- c(5, 3, 8); n <- c(40, 50, 60)
  none <- test_threshold(x, n, 0.05, icc = 0.05)
  big  <- test_threshold(x, n, 0.05, icc = 0.05, fpc_N = 1e12)
  expect_equal(big$n_eff_adj, none$n_eff_adj, tolerance = 1e-8)
  expect_equal(big$p_value, none$p_value, tolerance = 1e-8)
})

test_that("TH-F-4: testing almost the whole population makes the interval very narrow", {
  # 8 of 100 with fpc_N = 101: effective size 10000 and 800 effective positives,
  # so the interval is binom.test(800, 10000)'s: 0.07475563 to 0.08549047
  res <- test_threshold(8, 100, 0.05, alternative = "two.sided", fpc_N = 101)
  b   <- binom.test(800, 10000)$conf.int
  expect_equal(c(res$ci_lower, res$ci_upper), as.numeric(b), tolerance = 1e-9)
  expect_equal(c(res$ci_lower, res$ci_upper), c(0.07475563, 0.08549047), tolerance = 1e-7)
})

test_that("TH-F-5: the correction uses the number of people tested, not the effective size", {
  # 150 people, icc 0.05, fpc_N 1000: n_w = (40^2 + 50^2 + 60^2) / 150 = 50.67,
  # deff = 1 + 49.67 * 0.05 = 3.516667, n_eff = 150 / 3.516667 = 42.65403,
  # factor (1000 - 150) / 999, n_eff_adj = 42.65403 / 0.8508509 = 50.13103
  res <- test_threshold(c(5, 3, 8), c(40, 50, 60), 0.05, icc = 0.05, fpc_N = 1000)
  expect_equal(res$n_eff, 42.65403, tolerance = 1e-7)
  expect_equal(res$n_eff_adj, 50.13103, tolerance = 1e-7)
  # Not (1000 - n_eff) / 999, which would give 44.50990
  expect_false(isTRUE(all.equal(res$n_eff_adj, 44.5099, tolerance = 1e-4)))
})

test_that("TH-F-6: clustering and fpc_N together: both are applied", {
  # 16 positives of 150; the effective positives are 16 / 150 * 50.13103 = 5.347310
  # and the greater p-value is pbeta(0.05, 5.347310, 50.13103 - 5.347310 + 1) = 0.0750499
  res <- test_threshold(c(5, 3, 8), c(40, 50, 60), 0.05, icc = 0.05, fpc_N = 1000)
  expect_equal(res$statistic, 5.347310, tolerance = 1e-6)
  expect_equal(res$p_value, 0.0750499, tolerance = 1e-6)
  clustered_only <- test_threshold(c(5, 3, 8), c(40, 50, 60), 0.05, icc = 0.05)
  fpc_only       <- test_threshold(c(5, 3, 8), c(40, 50, 60), 0.05, fpc_N = 1000)
  expect_lt(res$p_value, clustered_only$p_value)
  expect_gt(res$p_value, fpc_only$p_value)
})


# ---------------------------------------------------------------------------
# Validation (TH-V-)
# ---------------------------------------------------------------------------

test_that("TH-V-1: x and n are checked as in test_difference", {
  e <- function(expr, msg) expect_error(expr, msg, fixed = TRUE)
  e(test_threshold(NA, 10, 0.05),
    "`x` and `n` must be numeric, not logical (got class `logical` for x, `numeric` for n). Note: a plain `NA` is logical in R. Remove missing observations before calling.")
  e(test_threshold("8", 100, 0.05),
    "`x` and `n` must be numeric (got class `character` for x, `numeric` for n).")
  e(test_threshold(numeric(0), numeric(0), 0.05),
    "`x` and `n` must each contain at least one value. Supply at least one site's count and total.")
  # Only one of the two empty: the same message
  e(test_threshold(numeric(0), 10, 0.05),
    "`x` and `n` must each contain at least one value. Supply at least one site's count and total.")
  e(test_threshold(8, numeric(0), 0.05),
    "`x` and `n` must each contain at least one value. Supply at least one site's count and total.")
  e(test_threshold(c(1, NA_real_), c(10, 10), 0.05),
    "`x` contains a missing or infinite value (found x[2] = NA). The number of positives at each site must be a finite non-negative integer.")
  e(test_threshold(1, Inf, 0.05),
    "`n` contains a missing or infinite value (found n[1] = Inf). The number tested at each site must be a finite positive integer.")
  e(test_threshold(1:2, 10, 0.05),
    "`x` and `n` must have the same length (got length(x) = 2, length(n) = 1). Each value of `x` is the number of positives at one site, and each value of `n` is the number tested at that site.")
  e(test_threshold(-1, 10, 0.05), "`x` must be non-negative (found x[1] = -1). Counts cannot be negative.")
  e(test_threshold(c(1, 0), c(10, 0), 0.05),
    "`n` must be positive for every site (found n[2] = 0). Each site must have tested at least one person.")
  e(test_threshold(1.5, 10, 0.05), "`x` must contain integers, not decimals (found x[1] = 1.5).")
  e(test_threshold(1, 10.5, 0.05), "`n` must contain integers, not decimals (found n[1] = 10.5).")
  # A value within rounding of a whole number gets a hint to use round()
  e(test_threshold(3, 100 * 1.1, 0.05),
    "`n` must contain integers, not decimals (found n[1] = 110, which is not exactly an integer; use round(n)).")
  e(test_threshold(0.57 * 100, 100, 0.05),
    "`x` must contain integers, not decimals (found x[1] = 57, which is not exactly an integer; use round(x)).")
  # A real fraction gets no hint, however large the number
  e(test_threshold(3, 1e10 + 0.5, 0.05),
    "`n` must contain integers, not decimals (found n[1] = 10000000000.5).")
  e(test_threshold(60, 50, 0.05),
    "`x` cannot exceed `n` (found x[1] = 60 > n[1] = 50). The number of positives cannot be more than the number tested at that site.")
  # One more positive than people tested, at the second site; x equal to n is fine
  e(test_threshold(c(5, 101), c(100, 100), 0.05),
    "`x` cannot exceed `n` (found x[2] = 101 > n[2] = 100). The number of positives cannot be more than the number tested at that site.")
  expect_no_error(test_threshold(c(5, 100), c(100, 100), 0.05))
  e(test_threshold(1, 2^53 + 2, 0.05),
    "`n` contains a value larger than 2^53 = 9007199254740992, the largest whole number R can store exactly (found n[1] = 9007199254740994). No real site tests this many people.")
  e(test_threshold(1, 1e200, 0.05), "(found n[1] = 1e+200)")
  # Numbers are printed in full whatever the scipen option
  old <- options(scipen = -9); on.exit(options(old), add = TRUE)
  e(test_threshold(1, 100000.5, 0.05), "(found n[1] = 100000.5)")
  options(old)
  expect_no_error(test_threshold(1, 2^53, 0.05))
})

test_that("TH-V-2: threshold must be a single finite number strictly between 0 and 1", {
  e <- function(expr, msg) expect_error(expr, msg, fixed = TRUE)
  e(test_threshold(8, 100, threshold = c(0.05, 0.10)),
    "`threshold` must be a single number in (0, 1) (got class `numeric`, length 2).")
  e(test_threshold(8, 100, threshold = "0.05"),
    "`threshold` must be a single number in (0, 1) (got class `character`, length 1).")
  e(test_threshold(8, 100, threshold = NaN), "`threshold` must be a single finite number (got NaN).")
  e(test_threshold(8, 100, threshold = Inf), "`threshold` must be a single finite number (got Inf).")
  e(test_threshold(8, 100, threshold = -0.1),
    "`threshold` must be in (0, 1) (got -0.1). `threshold` is a prevalence and cannot be negative.")
  e(test_threshold(8, 100, threshold = 1.5),
    "`threshold` must be in (0, 1) (got 1.5). `threshold` is a prevalence and cannot be greater than 1.")
  for (bad in c(0, 1)) {
    e(test_threshold(8, 100, threshold = bad),
      paste0("`threshold` must be in (0, 1) (got ", bad, "). No prevalence can be below 0 or above 1, so there is nothing to test. Use the prevalence that would lead to action, for example 0.05."))
  }
  expect_no_error(test_threshold(8, 100, threshold = 0.5))
})

test_that("TH-V-3: alternative and ci_method must be one of the allowed words", {
  e <- function(expr, msg) expect_error(expr, msg, fixed = TRUE)
  e(test_threshold(8, 100, 0.05, alternative = "up"),
    "`alternative` must be one of 'greater', 'less', or 'two.sided' (got 'up').")
  e(test_threshold(8, 100, 0.05, alternative = 1),
    "`alternative` must be one of 'greater', 'less', or 'two.sided' (got class `numeric`, length 1).")
  e(test_threshold(8, 100, 0.05, ci_method = "wilson"),
    "`ci_method` must be one of 'clopper-pearson', 'agresti-coull', or 'wald' (got 'wilson').")
  e(test_threshold(8, 100, 0.05, ci_method = "Wald"),
    "`ci_method` must be one of 'clopper-pearson', 'agresti-coull', or 'wald' (got 'Wald').")
  e(test_threshold(8, 100, 0.05, ci_method = c("wald", "wald")),
    "`ci_method` must be one of 'clopper-pearson', 'agresti-coull', or 'wald' (got class `character`, length 2).")
  # Valid values: see TH-R-2 (clopper-pearson, agresti-coull, wald)
})

test_that("TH-V-4: sensitivity and specificity: type, NaN and Inf, and range", {
  e <- function(expr, msg) expect_error(expr, msg, fixed = TRUE)
  e(test_threshold(30, 100, 0.2, sensitivity = TRUE),
    "`sensitivity` must be a single number in (0, 1] (got class `logical`, length 1). Note: `TRUE`, `FALSE` and a plain `NA` are logical, not numeric. Pass 1 if the test never misses a case.")
  e(test_threshold(30, 100, 0.2, sensitivity = list(0.9)),
    "`sensitivity` must be a single number in (0, 1] (got class `list`, length 1).")
  e(test_threshold(30, 100, 0.2, specificity = "a"),
    "`specificity` must be a single number in (0, 1] (got class `character`, length 1). Note: `TRUE`, `FALSE` and a plain `NA` are logical, not numeric. Pass 1 if the test never gives a false positive.")
  e(test_threshold(30, 100, 0.2, specificity = NULL),
    "`specificity` must be a single number in (0, 1] (got class `NULL`, length 0).")
  e(test_threshold(30, 100, 0.2, sensitivity = NaN), "`sensitivity` must be a single finite number (got NaN).")
  e(test_threshold(30, 100, 0.2, specificity = Inf), "`specificity` must be a single finite number (got Inf).")
  e(test_threshold(30, 100, 0.2, sensitivity = 1.1),
    "`sensitivity` must be in (0, 1] (got 1.1). `sensitivity` is a diagnostic probability and cannot be greater than 1.")
  e(test_threshold(30, 100, 0.2, sensitivity = -0.5),
    "`sensitivity` must be in (0, 1] (got -0.5). `sensitivity` is a diagnostic probability and cannot be negative.")
  e(test_threshold(30, 100, 0.2, sensitivity = 0),
    "`sensitivity` must be in (0, 1] (got 0). A sensitivity of 0 means everyone who truly has the condition tests negative.")
  e(test_threshold(30, 100, 0.2, specificity = 1.1),
    "`specificity` must be in (0, 1] (got 1.1). `specificity` is a diagnostic probability and cannot be greater than 1.")
  e(test_threshold(30, 100, 0.2, specificity = 0),
    "`specificity` must be in (0, 1] (got 0). A specificity of 0 means everyone who truly does not have the condition tests positive.")
})

test_that("TH-V-5: sensitivity + specificity of 1 or less is rejected; below 1.1 warns, after all checks", {
  e <- function(expr, msg) expect_error(expr, msg, fixed = TRUE)
  e(test_threshold(8, 100, 0.05, sensitivity = 0.4, specificity = 0.4),
    "`sensitivity` + `specificity` must exceed 1 for the Rogan-Gladen correction (got 0.4 + 0.4 = 0.8). The correction divides by `sensitivity` + `specificity` - 1, so a sum of 1 or less cannot be corrected.")
  e(test_threshold(8, 100, 0.05, sensitivity = 0.5, specificity = 0.5),
    "(got 0.5 + 0.5 = 1).")
  # A sum above 1 only by rounding counts as 1 and shows as 1
  e(test_threshold(8, 100, 0.05, sensitivity = 0.3, specificity = 0.7 + 2 * .Machine$double.eps),
    "(got 0.3 + 0.7 = 1).")
  # "Only by rounding" means within 0.000000000001 (1e-12) of 1. A sum of
  # 1.0000000000005 is inside that and is rejected, shown as 1; a sum of
  # 1.00000000001 is outside it and runs (with the below-1.1 warning).
  e(test_threshold(60, 100, 0.5, sensitivity = 0.5, specificity = 0.5 + 5e-13),
    "(got 0.5 + 0.5000000000005 = 1).")
  expect_warning(test_threshold(60, 100, 0.5, sensitivity = 0.5, specificity = 0.5 + 1e-11),
                 "`sensitivity` + `specificity` = 1.00000000001, which is below 1.1.", fixed = TRUE)
  # A sum just below 1.1 warns
  expect_warning(test_threshold(60, 100, 0.5, sensitivity = 0.595, specificity = 0.5),
                 "`sensitivity` + `specificity` = 1.095, which is below 1.1.", fixed = TRUE)
  expect_warning(
    res <- test_threshold(30, 100, 0.05, sensitivity = 0.51, specificity = 0.50),
    "`sensitivity` + `specificity` = 1.01, which is below 1.1. Correcting for this much test error makes the estimated prevalence and its confidence interval very wide, and small changes in the assumed sensitivity or specificity will change them a lot.",
    fixed = TRUE)
  expect_true(is.finite(res$p_value))
  # 1.1 that comes out a hair below by rounding does not warn
  expect_silent(test_threshold(30, 100, 0.05, sensitivity = 0.15, specificity = 0.95))
  # A bad input gives just its error, without the warning first
  expect_no_warning(expect_error(
    test_threshold(30, 100, 0.05, sensitivity = 0.51, specificity = 0.50, fpc_N = 50),
    "`fpc_N` (50) is smaller", fixed = TRUE))
  expect_no_warning(expect_error(
    test_threshold(8, 50, 0.05, icc = 0.05, fpc_N = 50), "`fpc_N` (50) equals", fixed = TRUE))
})

test_that("TH-V-6: conf_level must be a single finite number in [0.5, 1)", {
  e <- function(expr, msg) expect_error(expr, msg, fixed = TRUE)
  e(test_threshold(30, 100, 0.2, conf_level = "0.95"),
    "`conf_level` must be a single number in [0.5, 1) (got class `character`, length 1).")
  e(test_threshold(30, 100, 0.2, conf_level = NaN), "`conf_level` must be a single finite number (got NaN).")
  e(test_threshold(30, 100, 0.2, conf_level = 1),
    "`conf_level` must be in [0.5, 1) (got 1). For example, use 0.95 for a 95% confidence interval.")
  for (bad in c(0, 0.3, 0.4999999)) {
    e(test_threshold(30, 100, 0.2, conf_level = bad),
      paste0("`conf_level` must be at least 0.5 (got ", bad, "). Below 0.5, when the prevalence is exactly at the threshold, the test can wrongly find it above or below the threshold more than half the time."))
  }
  expect_no_error(test_threshold(30, 100, 0.2, conf_level = 0.5))
})

test_that("TH-V-7: icc must be NULL or a single finite number in [0, 1]", {
  e <- function(expr, msg) expect_error(expr, msg, fixed = TRUE)
  e(test_threshold(30, 100, 0.2, icc = FALSE),
    "`icc` must be a single number in [0, 1] (got class `logical`, length 1). To estimate ICC from the data, leave `icc = NULL`.")
  e(test_threshold(30, 100, 0.2, icc = c(0.1, 0.2)),
    "`icc` must be a single number in [0, 1] (got class `numeric`, length 2).")
  e(test_threshold(30, 100, 0.2, icc = NaN),
    "`icc` must be a single finite number (got NaN). To estimate ICC from the data, leave `icc = NULL`.")
  e(test_threshold(30, 100, 0.2, icc = -0.1),
    "`icc` must be in [0, 1] (got -0.1). `icc` is a correlation and cannot be negative. To estimate ICC from the data, leave `icc = NULL`.")
  e(test_threshold(30, 100, 0.2, icc = 1.5),
    "`icc` must be in [0, 1] (got 1.5). `icc` is a correlation and cannot be greater than 1. To estimate ICC from the data, leave `icc = NULL`.")
  expect_no_error(test_threshold(c(5, 5), c(10, 10), 0.3, icc = 1))
})

test_that("TH-V-8: fpc_N must be NULL or a single finite positive integer", {
  e <- function(expr, msg) expect_error(expr, msg, fixed = TRUE)
  e(test_threshold(1, 10, 0.05, fpc_N = "100"),
    "`fpc_N` must be a single finite positive integer (got class `character`). `fpc_N` is the total population size. Set `fpc_N = NULL` to skip the finite-population correction.")
  e(test_threshold(1, 10, 0.05, fpc_N = c(100, 200)),
    "`fpc_N` must be a single finite positive integer (got length = 2).")
  for (bad in c(100.5, -5, 0, Inf)) {
    e(test_threshold(1, 10, 0.05, fpc_N = bad),
      paste0("`fpc_N` must be a single finite positive integer (got ", bad, ")."))
  }
  # Within rounding of a whole number: a hint to use round()
  e(test_threshold(1, 10, 0.05, fpc_N = 100 * 1.1),
    "`fpc_N` must be a single finite positive integer (got 110, which is not exactly an integer; use round(fpc_N)). `fpc_N` is the total population size. Set `fpc_N = NULL` to skip the finite-population correction.")
  expect_no_error(test_threshold(1, 10, 0.05, fpc_N = 11))
  # Above 2^53, the largest whole number R can store exactly (as in
  # design_threshold); 2^53 itself is allowed
  e(test_threshold(1, 10, 0.05, fpc_N = 2^53 + 2),
    "`fpc_N` (9007199254740994) is larger than 2^53 = 9007199254740992, the largest whole number R can store exactly. No real population is this large: set `fpc_N = NULL`.")
  e(test_threshold(1, 10, 0.05, fpc_N = 1e20),
    "`fpc_N` (1e+20) is larger than 2^53 = 9007199254740992, the largest whole number R can store exactly. No real population is this large: set `fpc_N = NULL`.")
  expect_no_error(test_threshold(1, 10, 0.05, fpc_N = 2^53))
})

test_that("TH-V-9: names and 1 x 1 matrices (also in fpc_N and icc) do not leak into the result", {
  res <- test_threshold(8, 100, threshold = c(t = 0.05), sensitivity = matrix(0.9),
                        specificity = c(s = 0.95), conf_level = c(c = 0.95))
  expect_null(names(res$threshold))
  expect_null(dim(res$sensitivity))
  expect_null(names(res$specificity))
  expect_null(names(res$conf_level))
  expect_null(names(res$p_value))
  # fpc_N: the same result as the plain number, with no names or dimensions
  plain <- test_threshold(10, 100, threshold = 0.05, fpc_N = 1000)
  for (f in list(c(a = 1000), matrix(1000))) {
    res <- test_threshold(10, 100, threshold = 0.05, fpc_N = f)
    expect_identical(res, plain)
  }
  # icc: likewise
  plain <- test_threshold(10, 100, threshold = 0.05, icc = 0.05)
  for (i in list(c(a = 0.05), matrix(0.05))) {
    res <- test_threshold(10, 100, threshold = 0.05, icc = i)
    expect_identical(res, plain)
  }
})


# ---------------------------------------------------------------------------
# Return value (TH-R-)
# ---------------------------------------------------------------------------

test_that("TH-R-1: return list contains all expected fields, in order, in every mode", {
  expected <- c("statistic", "p_value", "reject", "threshold", "alternative",
                "prevalence", "ci_lower", "ci_upper", "ci_method", "n_total",
                "n_eff", "n_eff_adj", "conf_level", "sensitivity", "specificity",
                "icc_used", "deff", "fpc_N")
  expect_named(test_threshold(x = 5, n = 100, threshold = 0.10), expected, ignore.order = FALSE)
  expect_named(test_threshold(x = c(5, 10), n = c(50, 50), threshold = 0.10, icc = 0.05),
               expected, ignore.order = FALSE)
  expect_named(test_threshold(x = 5, n = 100, threshold = 0.10, fpc_N = 1000),
               expected, ignore.order = FALSE)
})

test_that("TH-R-2: ci_method changes only the interval: Clopper-Pearson (default), Agresti-Coull or Wald", {
  base <- test_threshold(x = 8, n = 100, threshold = 0.05)
  expect_equal(base$ci_method, "clopper-pearson")
  ac <- test_threshold(x = 8, n = 100, threshold = 0.05, ci_method = "agresti-coull",
                       alternative = "two.sided")
  expect_equal(ac$ci_method, "agresti-coull")
  z  <- qnorm(0.975); nt <- 100 + z^2; pt <- (8 + z^2 / 2) / nt
  expect_equal(c(ac$ci_lower, ac$ci_upper),
               c(pt - z * sqrt(pt * (1 - pt) / nt), pt + z * sqrt(pt * (1 - pt) / nt)),
               tolerance = 1e-12)
  expect_equal(test_threshold(x = 8, n = 100, threshold = 0.05, ci_method = "agresti-coull",
                              alternative = "less")$ci_lower, 0)
  # The test itself is the same
  ac_g <- test_threshold(x = 8, n = 100, threshold = 0.05, ci_method = "agresti-coull")
  expect_equal(c(ac_g$statistic, ac_g$p_value), c(base$statistic, base$p_value))
  # One-sided Agresti-Coull adds z^2 / 2 cases with the one-sided z = 1.644854:
  # n_tilde = 100 + 2.705543 = 102.705543, p_tilde = (8 + 1.352772) / 102.705543
  # = 0.091064, lower = 0.091064 - 1.644854 * sqrt(0.091064 * 0.908936 / 102.705543)
  # = 0.044369
  expect_equal(ac_g$ci_lower, 0.04436896, tolerance = 1e-6)
  expect_equal(ac_g$ci_upper, 1)

  # Wald: 0.08 +/- 1.959964 * sqrt(0.08 * 0.92 / 100) = 0.08 +/- 0.0531725
  # = [0.0268275, 0.1331725]
  wd <- test_threshold(x = 8, n = 100, threshold = 0.05, ci_method = "wald",
                       alternative = "two.sided")
  expect_equal(wd$ci_method, "wald")
  expect_equal(c(wd$ci_lower, wd$ci_upper), c(0.0268275, 0.1331725), tolerance = 1e-6)
  expect_equal(c(wd$statistic, wd$p_value),
               c(base$statistic, test_threshold(8, 100, 0.05, alternative = "two.sided")$p_value))
  # One-sided ends use z = 1.644854: 0.08 - 1.644854 * 0.0271291 = 0.03537624
  expect_equal(test_threshold(8, 100, 0.05, ci_method = "wald")$ci_lower, 0.03537624,
               tolerance = 1e-6)
  expect_equal(test_threshold(8, 100, 0.05, ci_method = "wald")$ci_upper, 1)
  expect_equal(test_threshold(8, 100, 0.05, ci_method = "wald", alternative = "less")$ci_lower, 0)
  # With clustering and fpc_N it uses the effective size m = n_eff_adj
  wc <- test_threshold(x = c(0, 6, 1, 11, 4), n = c(50, 60, 40, 80, 70), threshold = 0.05,
                       icc = 0.03, fpc_N = 2000, ci_method = "wald", alternative = "two.sided")
  p <- 22 / 300; z <- qnorm(0.975); se <- sqrt(p * (1 - p) / wc$n_eff_adj)
  expect_equal(c(wc$ci_lower, wc$ci_upper), c(p - z * se, p + z * se), tolerance = 1e-12)
})

test_that("TH-R-3: Wald and Agresti-Coull with no positives or everyone positive; Wald warns when its interval has zero width", {
  # n = 50, two-sided 95 percent. Agresti-Coull at 0 of 50: z = 1.959964,
  # n_tilde = 53.8415, p_tilde = 1.92073 / 53.8415 = 0.035674,
  # upper = 0.035674 + 1.959964 * sqrt(0.035674 * 0.964326 / 53.8415) = 0.085216
  ac0 <- test_threshold(0, 50, 0.05, ci_method = "agresti-coull", alternative = "two.sided")
  expect_equal(c(ac0$ci_lower, ac0$ci_upper), c(0, 0.08521612), tolerance = 1e-7)
  ac1 <- test_threshold(50, 50, 0.05, ci_method = "agresti-coull", alternative = "two.sided")
  expect_equal(c(ac1$ci_lower, ac1$ci_upper), c(0.9147839, 1), tolerance = 1e-7)
  # A measured prevalence of 0 or 1 gives zero width, and a warning
  for (k in c(0, 50)) {
    expect_warning(res <- test_threshold(k, 50, 0.05, ci_method = "wald", alternative = "two.sided"),
                   paste0("The measured prevalence is ", k / 50, ", so the Wald interval has zero width. Use `ci_method = \"clopper-pearson\"` for an interval that allows for the uncertainty."),
                   fixed = TRUE)
    expect_equal(res$ci_lower, res$ci_upper)
  }
  expect_no_warning(test_threshold(1, 50, 0.05, ci_method = "wald"))
  expect_no_warning(test_threshold(0, 50, 0.05))
  # The warning is for Wald only: the other two methods give no warning at
  # 0 of 50 or 50 of 50, whatever the alternative
  for (m in c("clopper-pearson", "agresti-coull")) {
    for (alt in c("greater", "less", "two.sided")) {
      expect_no_warning(test_threshold(0, 50, 0.05, alternative = alt, ci_method = m))
      expect_no_warning(test_threshold(50, 50, 0.05, alternative = alt, ci_method = m))
    }
  }
  # One-sided: no warning when the interval is [0, 1], warning when it has
  # zero width
  expect_no_warning(res <- test_threshold(0, 50, 0.05, ci_method = "wald", alternative = "greater"))
  expect_equal(c(res$ci_lower, res$ci_upper), c(0, 1))
  expect_no_warning(res <- test_threshold(50, 50, 0.05, ci_method = "wald", alternative = "less"))
  expect_equal(c(res$ci_lower, res$ci_upper), c(0, 1))
  expect_warning(res <- test_threshold(50, 50, 0.05, ci_method = "wald", alternative = "greater"),
                 "has zero width", fixed = TRUE)
  expect_equal(c(res$ci_lower, res$ci_upper), c(1, 1))
  expect_warning(res <- test_threshold(0, 50, 0.05, ci_method = "wald", alternative = "less"),
                 "has zero width", fixed = TRUE)
  expect_equal(c(res$ci_lower, res$ci_upper), c(0, 0))
})
