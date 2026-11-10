# Hand-checked reference cases for test_difference()
#
#   p_g     = x_g / n_g                                (measured prevalence)
#   m_g     = n_g / deff_g / f_g                       (effective sample size;
#             f_g = (N_g - n_g) / (N_g - 1) with fpc_N, else 1)
#   p_pool  = (p_1 m_1 + p_2 m_2) / (m_1 + m_2)
#   z       = (p_1 - p_2) / sqrt(p_pool (1 - p_pool) (1 / m_1 + 1 / m_2))
#   score CI = every difference d with |z(d)| <= 1.959964, where z(d) is the
#             same test with the prevalences estimated again under
#             p_1 - p_2 = d; then divided by sensitivity + specificity - 1.
#             z(0) is the z above, so the CI excludes 0 exactly when the
#             test rejects. It has no closed form: TD-11 checks its ends
#             by a separate numerical calculation.
#
# TD-1  (8 of 100 versus 20 of 120, no clustering):
#   p_1 = 0.08, p_2 = 0.166667, difference = -0.086667
#   p_pool = 28 / 220 = 0.127273; se = 0.045126; z = -1.920545; p = 0.054789
#   score CI = [-0.174678, 0.001913]: it includes 0, as p is above 0.05


# ---------------------------------------------------------------------------
# Core scenarios (TD-)
# ---------------------------------------------------------------------------

test_that("TD-1: two single sites give the hand-checked test and interval", {
  res <- test_difference(8, 100, 20, 120)
  expect_equal(res$difference,  -0.086667, tolerance = 1e-4)
  expect_equal(res$statistic,   -1.920545, tolerance = 1e-6)
  expect_equal(res$p_value,      0.054789, tolerance = 1e-4)
  expect_false(res$reject)
  expect_equal(c(res$ci_lower, res$ci_upper), c(-0.174678, 0.001913), tolerance = 1e-5)
  # (A Wald interval, -0.086667 +/- 1.96 * 0.043513 = [-0.171951, -0.001382],
  # would exclude 0 although the test is not significant.)
})

test_that("TD-2: with no clustering, the test matches R's prop.test(correct = FALSE)", {
  for (s in list(c(8, 100, 20, 120), c(30, 200, 45, 210), c(0, 50, 6, 60),
                 c(15, 40, 15, 40), c(99, 100, 90, 100))) {
    for (alt in c("two.sided", "greater", "less")) {
      res <- suppressWarnings(test_difference(s[1], s[2], s[3], s[4], alternative = alt))
      ref <- suppressWarnings(stats::prop.test(c(s[1], s[3]), c(s[2], s[4]),
                                               correct = FALSE, alternative = alt))
      expect_equal(res$p_value, ref$p.value)
    }
  }
})

test_that("TD-3: the reject decision uses 1 - conf_level as the significance level", {
  # p = 0.0548 (TD-1): not significant at 0.05, significant at 0.10
  expect_false(test_difference(8, 100, 20, 120, conf_level = 0.95)$reject)
  expect_true(test_difference(8, 100, 20, 120,  conf_level = 0.90)$reject)
})

test_that("TD-4: the interval matches the direction of the test", {
  # "less": [-1, U]; "greater": [L, 1]
  less <- test_difference(8, 100, 20, 120, alternative = "less")
  expect_equal(less$ci_lower, -1)
  expect_lt(less$ci_upper, 0)
  expect_true(less$reject)                 # p = 0.0274
  greater <- test_difference(8, 100, 20, 120, alternative = "greater")
  expect_equal(greater$ci_upper, 1)
  expect_false(greater$reject)
  # A one-sided end uses z at 1 - alpha (1.644854), not 1.959964. 30/100
  # versus 15/100: Wald standard error sqrt(0.3 * 0.7 / 100 + 0.15 * 0.85 /
  # 100) = 0.058095, so 0.15 - 1.644854 * 0.058095 = 0.054443 ("greater")
  # and 0.15 + 1.644854 * 0.058095 = 0.245557 ("less")
  z1 <- qnorm(0.95)
  se <- sqrt(0.3 * 0.7 / 100 + 0.15 * 0.85 / 100)
  res <- test_difference(30, 100, 15, 100, alternative = "greater", ci_method = "wald")
  expect_equal(res$ci_lower, 0.15 - z1 * se)
  expect_equal(res$ci_lower, 0.054443, tolerance = 1e-5)
  res <- test_difference(30, 100, 15, 100, alternative = "less", ci_method = "wald")
  expect_equal(res$ci_upper, 0.15 + z1 * se)
  expect_equal(res$ci_upper, 0.245557, tolerance = 1e-5)
  # The same one-sided ends for the other two methods
  expect_equal(test_difference(30, 100, 15, 100, alternative = "greater")$ci_lower,
               0.053718, tolerance = 1e-5)
  expect_equal(test_difference(30, 100, 15, 100, alternative = "less")$ci_upper,
               0.245574, tolerance = 1e-5)
  expect_equal(test_difference(30, 100, 15, 100, alternative = "greater",
                               ci_method = "newcombe")$ci_lower, 0.052981, tolerance = 1e-5)
  expect_equal(test_difference(30, 100, 15, 100, alternative = "less",
                               ci_method = "newcombe")$ci_upper, 0.243914, tolerance = 1e-5)
})

test_that("TD-5: an imperfect test leaves the test unchanged and divides the difference and interval by Se + Sp - 1", {
  # The test is on the measured prevalences, so z and p are as in TD-1;
  # the difference and interval are divided by 0.9 + 0.95 - 1 = 0.85
  res  <- test_difference(8, 100, 20, 120, sensitivity = 0.9, specificity = 0.95)
  base <- test_difference(8, 100, 20, 120)
  expect_equal(res$statistic, base$statistic)
  expect_equal(res$p_value,   base$p_value)
  expect_equal(res$difference, base$difference / 0.85)
  expect_equal(c(res$ci_lower, res$ci_upper), c(base$ci_lower, base$ci_upper) / 0.85)
  # Each group: (p - 0.05) / 0.85: (0.08 - 0.05) / 0.85 = 0.035294
  expect_equal(res$prevalence1, 0.035294, tolerance = 1e-4)
  expect_equal(res$prevalence2, 0.137255, tolerance = 1e-4)
})

test_that("TD-6: difference is the corrected measured difference, even when a group prevalence is set to 0", {
  # Group 1: (0.02 - 0.05) / 0.85 is below 0, so prevalence1 is set to 0.
  # difference = (0.02 - 0.166667) / 0.85 = -0.172549, not 0 - 0.137255.
  res <- test_difference(2, 100, 20, 120, sensitivity = 0.9, specificity = 0.95)
  expect_equal(res$prevalence1, 0)
  expect_equal(res$difference,  -0.172549, tolerance = 1e-4)
  expect_false(isTRUE(all.equal(res$difference, res$prevalence1 - res$prevalence2)))
})

test_that("TD-7: the Newcombe interval matches Newcombe (1998)", {
  # Values from Newcombe (1998), Table II (method 10)
  for (case in list(list(x = c(56, 70, 48, 80), ci = c(0.0524, 0.3339)),
                    list(x = c(9, 10, 3, 10),   ci = c(0.1705, 0.8090)),
                    list(x = c(6, 7, 2, 7),     ci = c(0.0582, 0.8062)),
                    list(x = c(5, 56, 0, 29),   ci = c(-0.0381, 0.1926)),
                    list(x = c(10, 10, 0, 20),  ci = c(0.6791, 1)))) {
    s   <- case$x
    res <- test_difference(s[1], s[2], s[3], s[4], ci_method = "newcombe")
    expect_equal(round(c(res$ci_lower, res$ci_upper), 4), case$ci)
  }
  # The method changes only the interval, not the test
  expect_equal(test_difference(8, 100, 20, 120, ci_method = "newcombe")$p_value,
               test_difference(8, 100, 20, 120)$p_value)
})

test_that("TD-8: swapping the groups flips the sign of the difference and keeps the two-sided p-value", {
  a <- test_difference(8, 100, 20, 120)
  b <- test_difference(20, 120, 8, 100)
  expect_equal(b$difference, -a$difference)
  expect_equal(b$p_value,    a$p_value)
  expect_equal(c(b$ci_lower, b$ci_upper), -c(a$ci_upper, a$ci_lower))
})

test_that("TD-9: both groups at 0 positives (or all positive) stops; 0 versus all positive gives an interval near -1", {
  expect_error(test_difference(0, 50, 0, 60),
               "Both groups have a measured prevalence of 0, so there is no variation and the test statistic cannot be calculated.",
               fixed = TRUE)
  expect_error(test_difference(50, 50, 60, 60),
               "Both groups have a measured prevalence of 1, so there is no variation and the test statistic cannot be calculated.",
               fixed = TRUE)
  # 0 of 50 versus 60 of 60: the difference is -1, so the lower end is -1,
  # but the interval still has real width
  expect_silent(res <- test_difference(0, 50, 60, 60))
  expect_equal(c(res$ci_lower, res$ci_upper), c(-1, -0.928652), tolerance = 1e-5)
  expect_lt(test_difference(0, 50, 60, 60, ci_method = "newcombe")$ci_upper, -0.9)
})

test_that("TD-10: counts stored as R integers give the same result", {
  expect_equal(test_difference(8L, 100L, 20L, 120L), test_difference(8, 100, 20, 120))
  # A named number or a 1 x 1 matrix works like a plain number, and nothing
  # of it leaks into the result
  expect_equal(test_difference(c(2, 5, 3), c(50, 60, 40), c(8, 11, 9), c(55, 70, 45),
                               icc = c(a = 0.05), conf_level = matrix(0.95)),
               test_difference(c(2, 5, 3), c(50, 60, 40), c(8, 11, 9), c(55, 70, 45),
                               icc = 0.05, conf_level = 0.95))
  # Values returned "as supplied" come back as plain numbers, without names
  res <- test_difference(8, 100, 20, 120, sensitivity = c(s = 0.9),
                         conf_level = c(a = 0.95), fpc_N = c(a = 500, b = 600))
  expect_identical(res[c("sensitivity", "conf_level", "fpc_N")],
                   list(sensitivity = 0.9, conf_level = 0.95, fpc_N = c(500, 600)))
})

test_that("TD-11: the score interval matches Mee (1984) and a separate calculation, including ends at -1 or 1", {
  # Newcombe (1998), Table II, 56/70 versus 48/80: Mee's score interval
  res <- test_difference(56, 70, 48, 80)
  expect_equal(round(c(res$ci_lower, res$ci_upper), 4), c(0.0533, 0.3377))
  # At each end the score statistic equals the z value. Checked here by
  # finding the prevalences under each difference with optimize() instead
  # of the closed form the function uses.
  stat <- function(r, d) {
    p1 <- r$apparent_prev1; p2 <- r$apparent_prev2; m1 <- r$n_eff_adj1; m2 <- r$n_eff_adj2
    ll <- function(q2) {
      q1 <- q2 + d
      p1 * m1 * log(q1) + (1 - p1) * m1 * log(1 - q1) +
        (if (p2 > 0) p2 * m2 * log(q2) else 0) + (1 - p2) * m2 * log(1 - q2)
    }
    q2 <- optimize(ll, c(max(0, -d) + 1e-12, min(1, 1 - d) - 1e-12),
                   maximum = TRUE, tol = 1e-13)$maximum
    q1 <- q2 + d
    (p1 - p2 - d) / sqrt(q1 * (1 - q1) / m1 + q2 * (1 - q2) / m2)
  }
  z <- qnorm(0.975)
  # An end at 1: 10 of 10 versus 0 of 20 (difference 1)
  res <- test_difference(10, 10, 0, 20)
  expect_equal(c(res$ci_lower, res$ci_upper), c(0.722467, 1), tolerance = 1e-5)
  # An end at -1 although the difference is not -1: with 5 billion people
  # per group (1 positive versus all positive) the statistic never reaches
  # the z value before -1
  res <- test_difference(1, 5e9, 5e9, 5e9)
  expect_equal(res$ci_lower, -1)
  expect_lt(res$ci_upper, -0.99999)
  # A difference of 1 with no variation left gives an infinite statistic
  expect_equal(.score_stat(0, 50, 1, 60, 1), -Inf)
  # A prevalence of 1 in 10 million: the ends are still exact (checked
  # against a separate high-precision calculation: upper end 5.66493e-07)
  res <- test_difference(1, 1e7, 0, 3)
  expect_equal(res$ci_upper, 5.66493e-07, tolerance = 1e-5)
  # At a difference of 0 the score statistic is the test's own z, even for
  # a tiny prevalence in a huge group (1 of 100 million against 1 of 3)
  for (r in list(test_difference(1, 1e8, 1, 3), test_difference(20, 1e7, 3, 3),
                 test_difference(2, 1e6, 1, 3)))
    expect_equal(.score_stat(r$apparent_prev1, r$n_eff_adj1, r$apparent_prev2, r$n_eff_adj2, 0),
                 r$statistic, tolerance = 1e-9)
  # One group 10 million times larger than the other (via the population
  # correction) and prevalences of exactly 0 and 1: the right solution is
  # found, and swapping the groups mirrors the interval
  res <- test_difference(0, 1e7, 1, 1, sensitivity = 0.8, specificity = 0.999, fpc_N = c(1e7 + 1, 2))
  expect_equal(c(res$ci_lower, res$ci_upper), c(-1, -0.2585098), tolerance = 1e-6)
  sw  <- test_difference(1, 1, 0, 1e7, sensitivity = 0.8, specificity = 0.999, fpc_N = c(2, 1e7 + 1))
  expect_equal(c(sw$ci_lower, sw$ci_upper), -c(res$ci_upper, res$ci_lower))
  # Huge samples make the interval very narrow (here about 1e-12 wide, with
  # an effective size of 1e24): its ends are still exact, so it still agrees
  # with the test
  res <- test_difference(1e11, 1e12, 1e11 + 1, 1e12, fpc_N = 1e12 + 1)
  expect_true(res$reject)
  expect_lt(res$ci_upper, 0)
  expect_equal(.score_stat(res$apparent_prev1, res$n_eff_adj1, res$apparent_prev2, res$n_eff_adj2,
                           res$ci_lower), z, tolerance = 1e-6)
  expect_equal(.score_stat(res$apparent_prev1, res$n_eff_adj1, res$apparent_prev2, res$n_eff_adj2,
                           res$ci_upper), -z, tolerance = 1e-6)
  for (r in list(test_difference(56, 70, 48, 80), test_difference(8, 100, 20, 120),
                 test_difference(c(2, 5, 3), c(50, 60, 40), c(8, 11, 9), c(55, 70, 45), icc = 0.05),
                 test_difference(8, 100, 20, 120, fpc_N = 500))) {
    expect_equal(stat(r, r$ci_lower), z, tolerance = 1e-6)
    expect_equal(stat(r, r$ci_upper), -z, tolerance = 1e-6)
  }
})

test_that("TD-12: the score interval excludes 0 exactly when the test rejects", {
  # Over many data sets, each alternative and two confidence levels,
  # including clustered data (fractional effective sizes)
  set.seed(1)
  for (i in 1:300) {
    n1 <- sample(10:120, 1); n2 <- sample(10:120, 1)
    x1 <- rbinom(1, n1, runif(1, 0, 0.4)); x2 <- rbinom(1, n2, runif(1, 0, 0.4))
    if (x1 + x2 == 0) next
    for (alt in c("two.sided", "greater", "less")) for (cl in c(0.9, 0.95)) {
      r <- test_difference(x1, n1, x2, n2, alternative = alt, conf_level = cl)
      expect_identical(r$ci_lower > 0 || r$ci_upper < 0, r$reject)
    }
  }
  for (i in 1:100) {
    x1 <- rbinom(4, 30, runif(1, 0.05, 0.3)); x2 <- rbinom(4, 30, runif(1, 0.05, 0.3))
    if (sum(x1) + sum(x2) == 0) next
    r <- test_difference(x1, rep(30, 4), x2, rep(30, 4), icc = 0.05)
    expect_identical(r$ci_lower > 0 || r$ci_upper < 0, r$reject)
  }
})

test_that("TD-13: one group nearly a census and all positive, against a group below 1, is tested, not stopped", {
  # Group 2 tests 1e9 of 1e9 + 1 people (effective size about 1e18),
  # all positive; group 1 is 98 of 100. The pooled prevalence would round to
  # exactly 1 and wrongly stop with 'Both groups have a measured prevalence
  # of 1'; the difference of -0.02 is real and clearly significant.
  res <- test_difference(98, 100, 1e9, 1e9, fpc_N = c(1e7, 1e9 + 1))
  expect_equal(res$difference, -0.02)
  expect_equal(c(res$ci_lower, res$ci_upper), c(-0.0700114, -0.005502), tolerance = 1e-5)
  expect_true(res$reject)
  # Swapping the groups mirrors the result
  sw <- test_difference(1e9, 1e9, 98, 100, fpc_N = c(1e9 + 1, 1e7))
  expect_equal(c(sw$ci_lower, sw$ci_upper), -c(res$ci_upper, res$ci_lower))
  expect_true(sw$reject)
})


# ---------------------------------------------------------------------------
# Clustering (TD-C-)
# ---------------------------------------------------------------------------

test_that("TD-14: the Wald interval is the measured difference plus or minus z standard errors", {
  # 8/100 versus 20/120: 0.08 - 0.166667 = -0.086667, standard error
  # sqrt(0.08 * 0.92 / 100 + 0.166667 * 0.833333 / 120) = 0.043513, so
  # -0.086667 +/- 1.959964 * 0.043513 = [-0.171951, -0.001382]
  res <- test_difference(8, 100, 20, 120, ci_method = "wald")
  expect_equal(c(res$ci_lower, res$ci_upper), c(-0.171951, -0.001382), tolerance = 1e-5)
  # The test is the same for every method; here the Wald interval leaves
  # out 0 although the test does not find a difference (p = 0.0548)
  expect_equal(res$p_value, test_difference(8, 100, 20, 120)$p_value)
  expect_equal(res$p_value, 0.054789, tolerance = 1e-5)
  expect_false(res$reject)
  expect_identical(res$ci_method, "wald")
  # One-sided "less" with an imperfect test: upper end only, divided by
  # Se + Sp - 1 = 0.85: (-0.086667 + 1.644854 * 0.043513) / 0.85 = -0.017757
  res <- test_difference(8, 100, 20, 120, alternative = "less", ci_method = "wald",
                         sensitivity = 0.9, specificity = 0.95)
  expect_equal(res$ci_lower, -1)
  expect_equal(res$ci_upper, -0.017757, tolerance = 1e-4)
  # Each group at 0 or 1: zero width, with a warning
  expect_warning(res <- test_difference(0, 50, 50, 50, ci_method = "wald"),
                 "Each group's measured prevalence is 0 or 1 (0 and 1), so the Wald interval has zero width. Use `ci_method = \"score\"` for an interval that allows for the uncertainty.",
                 fixed = TRUE)
  expect_equal(c(res$ci_lower, res$ci_upper), c(-1, -1))
  # The other methods do not warn here
  expect_no_warning(test_difference(0, 50, 50, 50))
  expect_no_warning(test_difference(0, 50, 50, 50, ci_method = "newcombe"))
})

test_that("TD-15: the Wald and Newcombe intervals use each group's effective size, with clustering and fpc_N", {
  # Three sites per group, icc 0.05, fpc_N 500 and 600.
  # Group 1: 10 of 150, deff 1 + (50 - 1) * 0.05 = 3.45,
  #   m1 = 150 / 3.45 / ((500 - 150) / 499) = 61.987578
  # Group 2: 28 of 170, deff 1 + (170 / 3 - 1) * 0.05 = 3.783333,
  #   m2 = 170 / 3.783333 / ((600 - 170) / 599) = 62.593997
  x1 <- c(2, 5, 3); n1 <- c(50, 60, 40); x2 <- c(8, 11, 9); n2 <- c(55, 70, 45)
  m1 <- 150 / 3.45 / ((500 - 150) / 499)
  m2 <- 170 / (1 + (170 / 3 - 1) * 0.05) / ((600 - 170) / 599)
  p1 <- 10 / 150; p2 <- 28 / 170; z <- qnorm(0.975)
  # Wald: (p1 - p2) +/- z * sqrt(p1 (1 - p1) / m1 + p2 (1 - p2) / m2)
  res <- test_difference(x1, n1, x2, n2, icc = 0.05, fpc_N = c(500, 600), ci_method = "wald")
  expect_equal(c(res$n_eff_adj1, res$n_eff_adj2), c(m1, m2))
  expect_equal(c(res$ci_lower, res$ci_upper),
               (p1 - p2) + c(-1, 1) * z * sqrt(p1 * (1 - p1) / m1 + p2 * (1 - p2) / m2))
  expect_equal(c(res$ci_lower, res$ci_upper), c(-0.208941, 0.012863), tolerance = 1e-5)
  # Newcombe: a Wilson interval for each group with m1 and m2 people,
  # written out here: centre (p + z^2 / (2 m)) / (1 + z^2 / m), half-width
  # z / (1 + z^2 / m) * sqrt(p (1 - p) / m + z^2 / (4 m^2))
  wilson <- function(p, m) {
    centre <- (p + z^2 / (2 * m)) / (1 + z^2 / m)
    half   <- z / (1 + z^2 / m) * sqrt(p * (1 - p) / m + z^2 / (4 * m^2))
    c(centre - half, centre + half)
  }
  w1 <- wilson(p1, m1); w2 <- wilson(p2, m2)
  res <- test_difference(x1, n1, x2, n2, icc = 0.05, fpc_N = c(500, 600), ci_method = "newcombe")
  expect_equal(c(res$ci_lower, res$ci_upper),
               c((p1 - p2) - sqrt((p1 - w1[1])^2 + (w2[2] - p2)^2),
                 (p1 - p2) + sqrt((w1[2] - p1)^2 + (p2 - w2[1])^2)))
  expect_equal(c(res$ci_lower, res$ci_upper), c(-0.215729, 0.017644), tolerance = 1e-5)
})

test_that("TD-16: a higher conf_level gives a wider interval, for each method", {
  # 8/100 versus 20/120 at 90, 95 and 99 percent
  for (m in c("score", "wald", "newcombe")) {
    ci <- sapply(c(0.90, 0.95, 0.99), function(cl) {
      r <- test_difference(8, 100, 20, 120, conf_level = cl, ci_method = m)
      c(r$ci_lower, r$ci_upper)
    })
    expect_true(all(diff(ci[1, ]) < 0))   # lower end moves down
    expect_true(all(diff(ci[2, ]) > 0))   # upper end moves up
  }
  # Score: [-0.159959, -0.013171], [-0.174678, 0.001913], [-0.204214, 0.032869]
  expect_equal(c(test_difference(8, 100, 20, 120, conf_level = 0.90)$ci_lower,
                 test_difference(8, 100, 20, 120, conf_level = 0.99)$ci_upper),
               c(-0.159959, 0.032869), tolerance = 1e-5)
})

test_that("TD-17: sensitivity = specificity = 0.999 gives almost the same result as 1", {
  # Se + Sp - 1 = 0.998: the difference and interval are divided by 0.998
  # instead of 1, and the test itself is unchanged
  a <- test_difference(8, 100, 20, 120, sensitivity = 0.999, specificity = 0.999)
  b <- test_difference(8, 100, 20, 120)
  expect_identical(a$p_value, b$p_value)
  expect_equal(a$difference, b$difference / 0.998)
  expect_equal(c(a$ci_lower, a$ci_upper), c(b$ci_lower, b$ci_upper) / 0.998)
  # Each group's prevalence: (0.08 - 0.001) / 0.998 = 0.079158 (0.08 with a
  # test that never errs)
  expect_equal(a$prevalence1, (0.08 - 0.001) / 0.998)
  expect_equal(c(a$difference, a$prevalence1, a$prevalence2),
               c(b$difference, b$prevalence1, b$prevalence2), tolerance = 0.02)
})

test_that("TD-C-1: a supplied icc gives each group its own design effect", {
  # Group 1: mean site size 50, deff = 1 + 49 * 0.05 = 3.45, n_eff = 150 / 3.45 = 43.48
  # Group 2: mean site size 56.67, deff = 1 + 55.67 * 0.05 = 3.78, n_eff = 170 / 3.78 = 44.93
  res <- test_difference(c(2, 5, 3), c(50, 60, 40), c(8, 11, 9), c(55, 70, 45), icc = 0.05)
  expect_equal(c(res$deff1, res$deff2),   c(3.45, 3.783333), tolerance = 1e-6)
  expect_equal(c(res$n_eff1, res$n_eff2), c(43.478261, 44.933921), tolerance = 1e-6)
  expect_equal(c(res$icc_used1, res$icc_used2), c(0.05, 0.05))
  expect_equal(res$p_value, 0.150855, tolerance = 1e-4)
  # Without clustering the same data are significant (p = 0.0068)
  expect_true(test_difference(c(2, 5, 3), c(50, 60, 40), c(8, 11, 9), c(55, 70, 45), icc = 0)$reject)
})

test_that("TD-C-2: with icc = NULL each group's design effect is estimated from its own sites", {
  # Group 1: 2, 9, 1 and 8 positives of 40 each -> p = 20 / 160 = 0.125
  #   var_obs = var(0.05, 0.225, 0.025, 0.2) = 0.010417
  #   var_srs = 0.125 * 0.875 / 40 = 0.002734; deff = 3.809524
  #   icc = (3.809524 - 1) / (40 - 1) = 0.072039
  res <- test_difference(c(2, 9, 1, 8), rep(40, 4), c(10, 3, 12, 15), rep(50, 4))
  expect_equal(res$deff1,     3.809524, tolerance = 1e-6)
  expect_equal(res$icc_used1, 0.072039, tolerance = 1e-4)
  expect_equal(res$deff2,     3.25,     tolerance = 1e-6)
  expect_equal(res$p_value,   0.318006, tolerance = 1e-4)

  # Extreme site differences (two sites all negative, two all positive):
  # the raw icc would be above 1, so it is capped at 1 and
  # deff = 1 + (10 - 1) * 1 = 10, as in estimate_prevalence
  ext <- test_difference(c(0, 0, 10, 10), rep(10, 4), c(3, 4, 5, 4), rep(10, 4))
  expect_equal(c(ext$icc_used1, ext$deff1), c(1, 10))
})

test_that("TD-C-3: a supplied icc is used for a group with one site; one person per site gives deff 1", {
  # One site of 100 with icc 0.05: deff = 1 + (100 - 1) * 0.05 = 5.95, so
  # the 100 people count as 100 / 5.95 = 16.8. No warning.
  expect_silent(res <- test_difference(8, 100, c(5, 7, 8), c(40, 40, 40), icc = 0.05))
  expect_identical(res$icc_used1, 0.05)
  expect_equal(res$deff1, 5.95)
  expect_equal(res$n_eff1, 100 / 5.95)
  expect_equal(res$deff2, 1 + (40 - 1) * 0.05)
  # It changes the answer: one site per group, 8 of 100 v 20 of 100. With
  # icc 0.05, m = 100 / 5.95 in each group and pooled prevalence 0.14, so
  # z = -0.12 / sqrt(0.14 * 0.86 * 2 / m) and p = 0.316 (0.014 without icc)
  res <- test_difference(8, 100, 20, 100, icc = 0.05)
  m <- 100 / 5.95
  expect_equal(res$p_value, 2 * pnorm(-0.12 / sqrt(0.14 * 0.86 * 2 / m)))
  expect_equal(res$p_value, 0.3160909, tolerance = 1e-6)
  expect_equal(test_difference(8, 100, 20, 100)$p_value, 0.01446846, tolerance = 1e-6)
  # One person per site: deff = 1 + (1 - 1) * 0.05 = 1, icc reported as given
  expect_silent(res <- test_difference(c(1, 0, 1), c(1, 1, 1), c(5, 7, 8), c(40, 40, 40), icc = 0.05))
  expect_identical(c(res$icc_used1, res$deff1), c(0.05, 1))
  # With icc = NULL the ICC cannot be estimated from one site: deff 1,
  # icc_used 0, no warning (a single count is also how a sample without
  # sites is given)
  expect_silent(res <- test_difference(8, 100, c(5, 7, 8), c(40, 40, 40)))
  expect_identical(c(res$icc_used1, res$deff1), c(0, 1))
  # With icc = NULL and one person per site there is no clustering to
  # estimate: icc_used 0, deff 1, also when every one of them is negative
  expect_silent(res <- test_difference(c(1, 0, 1, 0, 1), rep(1, 5), c(5, 7, 8), c(40, 40, 40)))
  expect_identical(c(res$icc_used1, res$deff1), c(0, 1))
  expect_silent(res <- test_difference(c(0, 0, 0), rep(1, 3), c(5, 7, 8), c(40, 40, 40)))
  expect_identical(c(res$icc_used1, res$deff1), c(0, 1))
})

test_that("TD-C-4: icc = 0 means no clustering; a tiny icc is used as given, without a warning", {
  base <- test_difference(c(2, 5, 3), c(50, 60, 40), c(8, 11, 9), c(55, 70, 45), icc = 0)
  expect_equal(c(base$deff1, base$deff2, base$icc_used1, base$icc_used2), c(1, 1, 0, 0))
  # 1e-12 is used, not set to 0: deff = 1 + (50 - 1) * 1e-12, which changes
  # nothing visible
  expect_silent(res <- test_difference(c(2, 5, 3), c(50, 60, 40), c(8, 11, 9), c(55, 70, 45), icc = 1e-12))
  expect_identical(res$icc_used1, 1e-12)
  expect_identical(res$deff1, 1 + 49 * 1e-12)
  expect_equal(res$p_value, base$p_value)
  # With very large sites a tiny icc matters: sites of 10 million with
  # icc 1.4e-8 have deff 1 + (1e7 - 1) * 1.4e-8 = 1.14
  res <- test_difference(c(5e5, 5e5), c(1e7, 1e7), c(5.1e5, 5.1e5), c(1e7, 1e7), icc = 1.4e-8)
  expect_equal(res$deff1, 1 + (1e7 - 1) * 1.4e-8)
  # One site with icc 1e-12: also used as given, without a warning
  expect_silent(res <- test_difference(8, 100, 20, 120, icc = 1e-12))
  expect_identical(res$icc_used1, 1e-12)
  expect_identical(res$deff1, 1 + 99 * 1e-12)
})

test_that("TD-C-5: the estimated design effect with sites of different sizes", {
  # Group 1: 1 of 20, 20 of 50, 3 of 30. Overall prevalence 24 / 100 = 0.24;
  # observed variance of 0.05, 0.4, 0.1 is 0.0358333; expected variance
  # mean(0.24 * 0.76 / c(20, 50, 30)) = 0.0062827; deff = 5.70352,
  # icc = (5.70352 - 1) / (33.333 - 1) = 0.145470
  res <- test_difference(c(1, 20, 3), c(20, 50, 30), c(5, 6, 4), c(40, 40, 40))
  expect_equal(res$deff1, 5.70352292, tolerance = 1e-8)
  expect_equal(res$icc_used1, 0.14546978, tolerance = 1e-7)
  # Group 2: site prevalences 0.125, 0.15, 0.1 vary less than chance, so
  # the design effect is floored at 1
  expect_equal(c(res$deff2, res$icc_used2), c(1, 0))
})


# ---------------------------------------------------------------------------
# Finite population correction (TD-F-)
# ---------------------------------------------------------------------------

test_that("TD-C-6: two icc values give each group its own ICC", {
  x1 <- c(2, 5, 3); n1 <- c(50, 60, 40); x2 <- c(8, 11, 9); n2 <- c(55, 70, 45)
  # Group 1: 50 per site, 1 + 49 * 0.02 = 1.98. Group 2: 170 / 3 per site,
  # 1 + (170 / 3 - 1) * 0.08 = 5.453333
  res <- test_difference(x1, n1, x2, n2, icc = c(0.02, 0.08))
  expect_equal(c(res$deff1, res$deff2), c(1.98, 5.453333), tolerance = 1e-6)
  expect_identical(c(res$icc_used1, res$icc_used2), c(0.02, 0.08))
  expect_equal(res$p_value, 0.116531, tolerance = 1e-5)
  # Each group's design effect is the one a single value would give it
  expect_identical(res$deff1, test_difference(x1, n1, x2, n2, icc = 0.02)$deff1)
  expect_identical(res$deff2, test_difference(x1, n1, x2, n2, icc = 0.08)$deff2)
  # Two equal values are the same as one
  expect_identical(test_difference(x1, n1, x2, n2, icc = c(0.05, 0.05)),
                   test_difference(x1, n1, x2, n2, icc = 0.05))
  # A group with one site uses its own value too: 1 + (100 - 1) * 0.05 = 5.95
  expect_silent(res <- test_difference(8, 100, x2, n2, icc = c(0.05, 0.1)))
  expect_identical(c(res$icc_used1, res$icc_used2), c(0.05, 0.1))
  expect_equal(res$deff1, 5.95)
  expect_equal(res$deff2, 1 + (170 / 3 - 1) * 0.1)
})

test_that("TD-C-7: the order of the sites and their names do not change the result", {
  # icc estimated from the sites (deff 4.34 and 3.69), so the sites matter
  x1 <- c(1, 12, 3); n1 <- c(50, 60, 40); x2 <- c(8, 20, 4); n2 <- c(55, 70, 45)
  res <- test_difference(x1, n1, x2, n2)
  expect_equal(c(res$deff1, res$deff2), c(4.343832, 3.686405), tolerance = 1e-6)
  expect_identical(test_difference(rev(x1), rev(n1), x2[c(2, 3, 1)], n2[c(2, 3, 1)]), res)
  expect_identical(test_difference(c(A = 1, B = 12, C = 3), c(A = 50, B = 60, C = 40), x2, n2), res)
})

test_that("TD-F-1: fpc_N gives each group a larger effective sample size", {
  # One value for both groups: m1 = 100 / ((500 - 100) / 499) = 124.75,
  # m2 = 120 / ((500 - 120) / 499) = 157.58
  res <- test_difference(8, 100, 20, 120, fpc_N = 500)
  expect_equal(c(res$n_eff_adj1, res$n_eff_adj2), c(124.75, 157.578947), tolerance = 1e-6)
  expect_equal(c(res$n_eff1, res$n_eff2), c(100, 120))
  expect_equal(res$p_value, 0.030623, tolerance = 1e-4)
  # Two values, one per group: m1 = 100 / ((200 - 100) / 199) = 199
  res2 <- test_difference(8, 100, 20, 120, fpc_N = c(200, 1e6))
  expect_equal(res2$n_eff_adj1, 199)
  expect_equal(res2$n_eff_adj2, 120 / ((1e6 - 120) / (1e6 - 1)))
})

test_that("TD-F-2: fpc_N smaller than, or equal to, a group's number tested is rejected", {
  expect_error(test_difference(8, 100, 20, 120, fpc_N = 110),
               "`fpc_N` for group 2 (110) is smaller than the number tested in that group (120): you cannot test more people than there are in the population.",
               fixed = TRUE)
  expect_error(test_difference(8, 100, 20, 120, fpc_N = c(100, 500)),
               "`fpc_N` for group 1 (100) equals the number tested in that group: the whole group was tested (a census), so its prevalence is known exactly, and this test needs a sample from each group. Set `fpc_N = NULL` if no finite-population correction is needed.",
               fixed = TRUE)
  # One more than the number tested works
  expect_no_error(test_difference(8, 100, 20, 120, fpc_N = c(101, 121)))
  # Large numbers are printed in full, not as 1e+05
  expect_error(test_difference(5, 1e5, 3, 10, fpc_N = 5e4),
               "`fpc_N` for group 1 (50000) is smaller than the number tested in that group (100000): you cannot test more people than there are in the population.",
               fixed = TRUE)
  expect_error(test_difference(5, 1e5, 3, 10, fpc_N = 1e5),
               "`fpc_N` for group 1 (100000) equals the number tested in that group: the whole group was tested (a census), so its prevalence is known exactly, and this test needs a sample from each group. Set `fpc_N = NULL` if no finite-population correction is needed.",
               fixed = TRUE)
})


# ---------------------------------------------------------------------------
# Validation (TD-V-)
# ---------------------------------------------------------------------------

test_that("TD-V-1: counts of the wrong type, empty, missing or infinite, or of different lengths are rejected, naming the group's arguments", {
  # Valid counts: every test from TD-1 on
  expect_error(test_difference(NA, 10, 3, 10),
               "`x1` and `n1` must be numeric, not logical (got class `logical` for x1, `numeric` for n1). Note: a plain `NA` is logical in R. Remove missing observations before calling.",
               fixed = TRUE)
  expect_error(test_difference(3, 10, "3", 10),
               "`x2` and `n2` must be numeric (got class `character` for x2, `numeric` for n2).", fixed = TRUE)
  expect_error(test_difference(numeric(0), numeric(0), 3, 10),
               "`x1` and `n1` must each contain at least one value. Supply at least one site's count and total.",
               fixed = TRUE)
  expect_error(test_difference(c(1, NA_real_), c(10, 10), 3, 10),
               "`x1` contains a missing or infinite value (found x1[2] = NA). The number of positives at each site must be a finite non-negative integer.",
               fixed = TRUE)
  expect_error(test_difference(3, 10, 3, c(10, Inf)),
               "`n2` contains a missing or infinite value (found n2[2] = Inf). The number tested at each site must be a finite positive integer.",
               fixed = TRUE)
  expect_error(test_difference(c(1, 2), 10, 3, 10),
               "`x1` and `n1` must have the same length (got length(x1) = 2, length(n1) = 1). Each value of `x1` is the number of positives at one site, and each value of `n1` is the number tested at that site.",
               fixed = TRUE)
})

test_that("TD-V-2: negative, zero, decimal or too-large counts are rejected, naming the position", {
  expect_error(test_difference(-0.5, 10, 3, 10),
               "`x1` must be non-negative (found x1[1] = -0.5). Counts cannot be negative.", fixed = TRUE)
  expect_error(test_difference(3, 10, 0, 0),
               "`n2` must be positive for every site (found n2[1] = 0). Each site must have tested at least one person.",
               fixed = TRUE)
  expect_error(test_difference(3, -10, 3, 10),
               "`n1` must be positive for every site (found n1[1] = -10). Each site must have tested at least one person.",
               fixed = TRUE)
  expect_error(test_difference(c(1, 2.5), c(10, 10), 3, 10),
               "`x1` must contain integers, not decimals (found x1[2] = 2.5).", fixed = TRUE)
  expect_error(test_difference(3, 10.5, 3, 10),
               "`n1` must contain integers, not decimals (found n1[1] = 10.5).", fixed = TRUE)
  # A value that only rounding keeps from being a whole number prints as that
  # number, so the message says why: 100 * 1.1 is 110.00000000000001 and
  # 0.57 * 100 is 56.99999999999999
  expect_error(test_difference(5, 100 * 1.1, 3, 10),
               "`n1` must contain integers, not decimals (found n1[1] = 110, which is not exactly an integer; use round(n1)).",
               fixed = TRUE)
  expect_error(test_difference(5, 100, c(1, 0.57 * 100), c(10, 100)),
               "`x2` must contain integers, not decimals (found x2[2] = 57, which is not exactly an integer; use round(x2)).",
               fixed = TRUE)
  # A real decimal, even a large one, gets no such hint
  expect_error(test_difference(5, 1e10 + 0.5, 3, 10),
               "`n1` must contain integers, not decimals (found n1[1] = 10000000000.5).",
               fixed = TRUE)
  expect_error(test_difference(3, 10, 11, 10),
               "`x2` cannot exceed `n2` (found x2[1] = 11 > n2[1] = 10). The number of positives cannot be more than the number tested at that site.",
               fixed = TRUE)
  # Large counts are printed in full, not as 1e+05
  expect_error(test_difference(c(5, 100001), c(10, 1e5), 3, 10),
               "`x1` cannot exceed `n1` (found x1[2] = 100001 > n1[2] = 100000). The number of positives cannot be more than the number tested at that site.",
               fixed = TRUE)
  # Counts above 2^53 are rejected; 2^53 itself is allowed. (Far larger
  # counts, like 1e150, used to crash the score interval.)
  big <- "contains a value larger than 2^53 = 9007199254740992, the largest whole number R can store exactly"
  expect_error(test_difference(1, 1e150, 5, 10),
               paste0("`n1` ", big, " (found n1[1] = 1e+150). No real site tests this many people."),
               fixed = TRUE)
  expect_error(test_difference(5, 10, c(1, 1), c(10, 2^53 + 2)),
               paste0("`n2` ", big, " (found n2[2] = 9007199254740994). No real site tests this many people."),
               fixed = TRUE)
  expect_no_error(test_difference(1, 2^53, 5, 10))
  # x = n is allowed: see TD-9
})

test_that("TD-V-3: a vector, text, TRUE/FALSE, a bare NA or NULL in a single-number argument is rejected, naming the class", {
  expect_error(test_difference(8, 100, 20, 120, sensitivity = TRUE),
               "`sensitivity` must be a single number in (0, 1] (got class `logical`, length 1). Note: `TRUE`, `FALSE` and a plain `NA` are logical, not numeric. Pass 1 if the test never misses a case.",
               fixed = TRUE)
  expect_error(test_difference(8, 100, 20, 120, specificity = NULL),
               "`specificity` must be a single number in (0, 1] (got class `NULL`, length 0). Note: `TRUE`, `FALSE` and a plain `NA` are logical, not numeric. Pass 1 if the test never gives a false positive.",
               fixed = TRUE)
  expect_error(test_difference(8, 100, 20, 120, conf_level = c(0.9, 0.95)),
               "`conf_level` must be a single number in [0.5, 1) (got class `numeric`, length 2).", fixed = TRUE)
  expect_error(test_difference(8, 100, 20, 120, icc = NA),
               "`icc` must be one or two numbers in [0, 1]: one for both groups, or c(group 1, group 2) (got class `logical`, length 1). To estimate ICC from the data, leave `icc = NULL`.",
               fixed = TRUE)
  expect_error(test_difference(8, 100, 20, 120, icc = "0.05"),
               "`icc` must be one or two numbers in [0, 1]: one for both groups, or c(group 1, group 2) (got class `character`, length 1). To estimate ICC from the data, leave `icc = NULL`.",
               fixed = TRUE)
  # icc takes one or two values (one per group), not three
  expect_error(test_difference(8, 100, 20, 120, icc = c(0.1, 0.2, 0.3)),
               "`icc` must be one or two numbers in [0, 1]: one for both groups, or c(group 1, group 2) (got class `numeric`, length 3). To estimate ICC from the data, leave `icc = NULL`.",
               fixed = TRUE)
})

test_that("TD-V-4: alternative and ci_method must be one of their names", {
  expect_error(test_difference(8, 100, 20, 120, alternative = "bigger"),
               "`alternative` must be one of 'two.sided', 'greater', or 'less' (got 'bigger').", fixed = TRUE)
  expect_error(test_difference(8, 100, 20, 120, alternative = NA),
               "`alternative` must be one of 'two.sided', 'greater', or 'less' (got class `logical`, length 1).", fixed = TRUE)
  expect_error(test_difference(8, 100, 20, 120, alternative = c("less", "greater")),
               "`alternative` must be one of 'two.sided', 'greater', or 'less' (got class `character`, length 2).", fixed = TRUE)
  expect_error(test_difference(8, 100, 20, 120, ci_method = "Newcombe"),
               "`ci_method` must be one of 'score', 'newcombe' or 'wald' (got 'Newcombe').", fixed = TRUE)
  expect_error(test_difference(8, 100, 20, 120, ci_method = 1),
               "`ci_method` must be one of 'score', 'newcombe' or 'wald' (got class `numeric`, length 1).", fixed = TRUE)
  expect_error(test_difference(8, 100, 20, 120, ci_method = c("score", "newcombe")),
               "`ci_method` must be one of 'score', 'newcombe' or 'wald' (got class `character`, length 2).", fixed = TRUE)
  # Valid values: see TD-1 (score), TD-7 (newcombe), TD-14 (wald)
})

test_that("TD-V-5: NaN, Inf or -Inf in a single-number argument gives the finite-value error", {
  expect_error(test_difference(8, 100, 20, 120, sensitivity = NaN),
               "`sensitivity` must be a single finite number (got NaN).", fixed = TRUE)
  expect_error(test_difference(8, 100, 20, 120, specificity = -Inf),
               "`specificity` must be a single finite number (got -Inf).", fixed = TRUE)
  expect_error(test_difference(8, 100, 20, 120, conf_level = Inf),
               "`conf_level` must be a single finite number (got Inf).", fixed = TRUE)
  expect_error(test_difference(8, 100, 20, 120, icc = NaN),
               "`icc` must be finite (got NaN). To estimate ICC from the data, leave `icc = NULL`.",
               fixed = TRUE)
  # With two values, the first bad one is named
  expect_error(test_difference(8, 100, 20, 120, icc = c(0.05, NA)),
               "`icc` must be finite (found icc[2] = NA). To estimate ICC from the data, leave `icc = NULL`.",
               fixed = TRUE)
  expect_error(test_difference(8, 100, 20, 120, icc = c(Inf, 0.1)),
               "`icc` must be finite (found icc[1] = Inf). To estimate ICC from the data, leave `icc = NULL`.",
               fixed = TRUE)
})

test_that("TD-V-6: sensitivity, specificity, icc and conf_level outside their ranges are rejected", {
  expect_error(test_difference(8, 100, 20, 120, sensitivity = 0),
               "`sensitivity` must be in (0, 1] (got 0). A sensitivity of 0 means everyone who truly has the condition tests negative.",
               fixed = TRUE)
  expect_error(test_difference(8, 100, 20, 120, sensitivity = 1.1),
               "`sensitivity` must be in (0, 1] (got 1.1). `sensitivity` is a diagnostic probability and cannot be greater than 1.",
               fixed = TRUE)
  # Large values are printed in full, not as 1e+05
  expect_error(test_difference(8, 100, 20, 120, sensitivity = 1e5),
               "`sensitivity` must be in (0, 1] (got 100000). `sensitivity` is a diagnostic probability and cannot be greater than 1.",
               fixed = TRUE)
  expect_error(test_difference(8, 100, 20, 120, sensitivity = -0.5),
               "`sensitivity` must be in (0, 1] (got -0.5). `sensitivity` is a diagnostic probability and cannot be negative.",
               fixed = TRUE)
  expect_error(test_difference(8, 100, 20, 120, specificity = 0),
               "`specificity` must be in (0, 1] (got 0). A specificity of 0 means everyone who truly does not have the condition tests positive.",
               fixed = TRUE)
  expect_error(test_difference(8, 100, 20, 120, specificity = -0.2),
               "`specificity` must be in (0, 1] (got -0.2). `specificity` is a diagnostic probability and cannot be negative.",
               fixed = TRUE)
  expect_error(test_difference(8, 100, 20, 120, specificity = 1.1),
               "`specificity` must be in (0, 1] (got 1.1). `specificity` is a diagnostic probability and cannot be greater than 1.",
               fixed = TRUE)
  expect_error(test_difference(8, 100, 20, 120, icc = -0.1),
               "`icc` must be in [0, 1] (got -0.1). `icc` is a correlation and cannot be negative. To estimate ICC from the data, leave `icc = NULL`.",
               fixed = TRUE)
  expect_error(test_difference(8, 100, 20, 120, icc = 1.5),
               "`icc` must be in [0, 1] (got 1.5). `icc` is a correlation and cannot be greater than 1. To estimate ICC from the data, leave `icc = NULL`.",
               fixed = TRUE)
  expect_error(test_difference(8, 100, 20, 120, icc = c(0.05, 1.5)),
               "`icc` must be in [0, 1] (found icc[2] = 1.5). `icc` is a correlation and cannot be greater than 1. To estimate ICC from the data, leave `icc = NULL`.",
               fixed = TRUE)
  expect_error(test_difference(8, 100, 20, 120, icc = c(-0.1, 0.5)),
               "`icc` must be in [0, 1] (found icc[1] = -0.1). `icc` is a correlation and cannot be negative. To estimate ICC from the data, leave `icc = NULL`.",
               fixed = TRUE)
  expect_error(test_difference(8, 100, 20, 120, icc = c(0.5, -0.2)),
               "`icc` must be in [0, 1] (found icc[2] = -0.2). `icc` is a correlation and cannot be negative. To estimate ICC from the data, leave `icc = NULL`.",
               fixed = TRUE)
  # icc = 1 is allowed: with sites of 50 and 56.67 the design effects are
  # 1 + 49 = 50 and 1 + 55.67 = 56.67
  res <- test_difference(c(2, 5, 3), c(50, 60, 40), c(8, 11, 9), c(55, 70, 45), icc = 1)
  expect_equal(c(res$deff1, res$deff2), c(50, 56.666667), tolerance = 1e-7)
  # A conf_level one rounding step below 1 still gives proper intervals,
  # for both methods (z is about 8.3, not infinite)
  for (cm in c("score", "newcombe")) {
    res <- test_difference(0, 1000, 1000, 1000, conf_level = 1 - 1e-16, ci_method = cm)
    expect_true(all(is.finite(c(res$ci_lower, res$ci_upper))))
    expect_gt(res$ci_upper, -1)
    expect_lt(res$ci_upper, -0.9)
  }
  # Wald too (0 versus 1000 of 1000 would have zero width: see TD-14):
  # -0.086667 +/- 8.2 * 0.043513 (1 - 1e-16 is stored as 1 - 1.1e-16, so
  # each tail is 5.55e-17)
  res <- test_difference(8, 100, 20, 120, conf_level = 1 - 1e-16, ci_method = "wald")
  expect_equal(c(res$ci_lower, res$ci_upper),
               -0.086667 + c(-1, 1) * stats::qnorm((1 - (1 - 1e-16)) / 2, lower.tail = FALSE) * 0.043513,
               tolerance = 1e-4)
  for (bad in c(1, 1.05)) {
    expect_error(test_difference(8, 100, 20, 120, conf_level = bad),
                 paste0("`conf_level` must be in [0.5, 1) (got ", bad, "). For example, use 0.95 for a 95% confidence interval."),
                 fixed = TRUE)
  }
  # Below 0.5 a one-sided interval would start above the estimate itself
  # (0.3: lower end 0.126 for a difference of 0.1), and with no real
  # difference the test would report one more than half the time.
  # 0 and negative values get the same message.
  for (bad in c(-0.1, 0, 0.3, 0.49, 0.4999999)) {
    expect_error(test_difference(20, 100, 10, 100, alternative = "greater", conf_level = bad),
                 paste0("`conf_level` must be at least 0.5 (got ", bad, "). Below 0.5, when there is ",
                        "no real difference, the test would wrongly report one more than half the time."),
                 fixed = TRUE)
  }
  # 0.5 itself works: z = qnorm(0.75) = 0.6745, a narrower interval than
  # the 95 percent one in TD-1
  res <- test_difference(8, 100, 20, 120, conf_level = 0.5)
  expect_equal(c(res$ci_lower, res$ci_upper), c(-0.116168, -0.057207), tolerance = 1e-5)
})

test_that("TD-V-7: sensitivity + specificity of 1 or less is rejected; below 1.1 warns with the full sum", {
  sum1 <- function(se, sp) paste0(
    "`sensitivity` + `specificity` must exceed 1 for the Rogan-Gladen correction (got ",
    se, " + ", sp, " = 1). The correction divides by `sensitivity` + `specificity` - 1, ",
    "so a sum of 1 or less cannot be corrected.")
  expect_error(test_difference(8, 100, 20, 120, sensitivity = 0.5, specificity = 0.5),
               sum1("0.5", "0.5"), fixed = TRUE)
  # A sum above 1 only by rounding counts as 1, and is shown as 1
  expect_error(test_difference(8, 100, 20, 120, sensitivity = 0.3,
                               specificity = 0.7 + 2 * .Machine$double.eps),
               sum1("0.3", "0.7"), fixed = TRUE)
  expect_error(test_difference(8, 100, 20, 120, sensitivity = 0.5, specificity = 0.5 + 1e-13),
               sum1("0.5", "0.5000000000001"), fixed = TRUE)
  expect_warning(test_difference(8, 100, 20, 120, sensitivity = 0.55, specificity = 0.5),
                 "`sensitivity` + `specificity` = 1.05, which is very close to 1. Correcting for this much test error makes the estimated difference and its confidence interval very wide, and small changes in the assumed sensitivity or specificity will change them a lot.",
                 fixed = TRUE)
  expect_warning(test_difference(8, 100, 20, 120, sensitivity = 0.5000001, specificity = 0.5),
                 "= 1.0000001, which is very close to 1.", fixed = TRUE)
  expect_no_warning(test_difference(8, 100, 20, 120, sensitivity = 0.61, specificity = 0.5))
  # A sum truly above 1, however little, is not refused (only rounding counts)
  expect_warning(res <- test_difference(8, 100, 20, 120, sensitivity = 0.5, specificity = 0.500000001),
                 "= 1.000000001, which is very close to 1.", fixed = TRUE)
  expect_equal(c(res$ci_lower, res$ci_upper), c(-1, 1))
  # A sum of 1.1 that rounding puts a hair below it (0.15 + 0.95) does not warn
  expect_no_warning(test_difference(8, 100, 20, 120, sensitivity = 0.15, specificity = 0.95))
  # The warning comes only after every check has passed
  expect_no_warning(expect_error(
    test_difference(8, 100, 20, 120, sensitivity = 0.55, specificity = 0.5, conf_level = 2),
    "`conf_level` must be in [0.5, 1) (got 2). For example, use 0.95 for a 95% confidence interval.",
    fixed = TRUE))
  expect_no_warning(expect_error(
    test_difference(0, 50, 0, 60, sensitivity = 0.55, specificity = 0.5),
    "Both groups have a measured prevalence of 0, so there is no variation and the test statistic cannot be calculated.",
    fixed = TRUE))
})

test_that("TD-V-8: fpc_N must be one or two positive integers", {
  fpc_msg <- function(got) paste0(
    "`fpc_N` must be one or two finite positive integers ", got, ". `fpc_N` is ",
    "the population size of each group: one value for both groups, or c(group 1, ",
    "group 2). Set `fpc_N = NULL` to skip the finite-population correction.")
  bad <- list("500", TRUE, c(500, 600, 700), 100.5, 0, -5, Inf, NA_real_,
              c(1e6, 1e7 + 0.5), c(0, 1e7 + 0.5), c(1e6, NaN))
  got <- c("(got class `character`)", "(got class `logical`)", "(got length = 3)",
           "(got 100.5)", "(got 0)", "(got -5)", "(got Inf)", "(got NA)",
           "(found fpc_N[2] = 10000000.5)", "(found fpc_N[1] = 0)",
           "(found fpc_N[2] = NaN)")
  for (i in seq_along(bad))
    expect_error(test_difference(8, 100, 20, 120, fpc_N = bad[[i]]), fpc_msg(got[i]), fixed = TRUE)
  # A value that only rounding keeps from being a whole number says why:
  # 100 * 1.1 is 110.00000000000001
  expect_error(test_difference(8, 100, 20, 120, fpc_N = 100 * 1.1),
               fpc_msg("(got 110, which is not exactly an integer; use round(fpc_N))"),
               fixed = TRUE)
  expect_error(test_difference(8, 100, 20, 120, fpc_N = c(1000, 100 * 1.1)),
               fpc_msg("(found fpc_N[2] = 110, which is not exactly an integer; use round(fpc_N))"),
               fixed = TRUE)
  # R's printing options do not change the value shown
  op <- options(scipen = -9, digits = 3)
  expect_error(test_difference(5, 10, 3, 10, fpc_N = 100000.5), fpc_msg("(got 100000.5)"), fixed = TRUE)
  options(op)
  # Valid values: 500 and c(200, 1e6) (TD-F-1)
})


# ---------------------------------------------------------------------------
# Return value (TD-R-)
# ---------------------------------------------------------------------------

test_that("TD-R-1: return list contains all expected fields, in order", {
  expect_named(test_difference(8, 100, 20, 120),
               c("difference", "difference_app", "statistic", "p_value", "reject",
                 "alternative", "prevalence1", "prevalence2", "apparent_prev1",
                 "apparent_prev2", "ci_lower", "ci_upper", "ci_method", "n_total1",
                 "n_total2", "n_eff1", "n_eff2", "n_eff_adj1", "n_eff_adj2",
                 "conf_level", "sensitivity", "specificity", "icc_used1", "icc_used2",
                 "deff1", "deff2", "fpc_N"),
               ignore.order = FALSE)
})


# ---------------------------------------------------------------------------
# Round trips with design_difference() (TD-T-)
# ---------------------------------------------------------------------------

test_that("TD-T-1: data at the planned rates and design_difference's n give a significant result", {
  # design_difference(0.10, 0.20) -> 199 per group; 20 and 40 positives
  n <- design_difference(0.10, 0.20)$n_per_group1
  expect_equal(n, 199)
  res <- test_difference(round(0.1 * n), n, round(0.2 * n), n)
  expect_true(res$reject)
  expect_equal(res$p_value, 0.005082, tolerance = 1e-3)
})

test_that("TD-T-2: the clustered design from design_difference gives a significant result when icc is supplied", {
  # design_difference(0.10, 0.20, n_sites = 30, icc = 0.05) -> 30 sites of 10;
  # 1 and 2 positives per site
  d <- design_difference(0.10, 0.20, n_sites = 30, icc = 0.05)
  expect_equal(c(d$n_sites1, d$n_per_site1), c(30, 10))
  res <- test_difference(rep(1, 30), rep(10, 30), rep(2, 30), rep(10, 30), icc = 0.05)
  expect_equal(res$deff1, 1.45)
  expect_true(res$reject)
})

test_that("TD-T-3: a design with different sites and ICCs in each group gives a significant result with those ICCs", {
  # design_difference(0.10, 0.20, n_sites = c(20, 40), icc = c(0.02, 0.05))
  # -> 20 and 40 sites of 10; 1 and 2 positives per site
  d <- design_difference(0.10, 0.20, n_sites = c(20, 40), icc = c(0.02, 0.05))
  expect_equal(c(d$n_sites1, d$n_sites2, d$n_per_site1, d$n_per_site2), c(20, 40, 10, 10))
  res <- test_difference(rep(1, 20), rep(10, 20), rep(2, 40), rep(10, 40), icc = c(0.02, 0.05))
  expect_equal(c(res$deff1, res$deff2), c(d$deff1, d$deff2))
  expect_true(res$reject)
})

test_that("TD-T-4: a one-site design from design_difference is analysed with its icc", {
  # design_difference(0.10, 0.20, n_sites = 1, icc = 0.002) -> one site of
  # 330 per group (deff 1 + 329 * 0.002 = 1.658). Data at exactly those
  # prevalences, analysed with the same icc, give z close to the design's
  # z_a + z_b = qnorm(0.975) + qnorm(0.8) = 2.80 (the design uses the
  # unpooled variance under the alternative, so not exactly)
  d <- design_difference(0.10, 0.20, n_sites = 1, icc = 0.002)
  expect_identical(c(d$n_sites1, d$n_per_site1), c(1, 330))
  res <- test_difference(33, 330, 66, 330, icc = 0.002)
  expect_equal(res$deff1, d$deff1)
  expect_equal(res$statistic, -2.793796, tolerance = 1e-6)
  expect_true(res$reject)
})
