# Hand-checked reference cases for estimate_prevalence()
#
# Method: Wald CI on the apparent (measured) prevalence, then the
# Rogan-Gladen correction applied to the estimate and both CI endpoints.
#
#   p_hat   = sum(x) / sum(n)                          (apparent prevalence)
#   n_eff   = n_total / deff
#   se      = sqrt(p_hat * (1 - p_hat) / n_eff)        (standard error)
#   moe_app = z * se,  z = qnorm(0.975) = 1.959964
#   CI_app  = [p_hat - moe_app, p_hat + moe_app]
#   p_true  = (p_hat - (1 - spec)) / (sens + spec - 1) (Rogan-Gladen;
#             sens = sensitivity, spec = specificity)
#   CI_true = Rogan-Gladen applied to both CI_app endpoints, clamped to [0, 1]
#
# EP-1  (one site, perfect test, x = 30, n = 100):
#   se  = sqrt(0.3 * 0.7 / 100) = 0.045826
#   moe = 1.959964 * 0.045826 = 0.089817
#   CI  = [0.210183, 0.389817]
#
# EP-2  (imperfect test, sensitivity 0.9, specificity 0.95, x = 30, n = 100):
#   correction = 0.9 + 0.95 - 1 = 0.85
#   p_true = (0.3 - 0.05) / 0.85 = 0.294118
#   moe    = 0.089817 / 0.85 = 0.105667
#   CI     = [(0.210183 - 0.05) / 0.85, (0.389817 - 0.05) / 0.85]
#          = [0.188451, 0.399785]
#
# EP-C-1  (icc = 0.05 supplied, 10 sites of 10, 3 positives each):
#   deff  = 1 + (10 - 1) * 0.05 = 1.45, n_eff = 100 / 1.45 = 68.9655
#   se    = sqrt(0.3 * 0.7 / 68.9655) = 0.055182
#   moe   = 1.959964 * 0.055182 = 0.108154
#   CI    = [0.191846, 0.408154]
#
# Values written to 6 decimal places are checked with tolerance = 1e-4
# (relative), which still catches any real mistake.


# ---------------------------------------------------------------------------
# Core scenarios (EP-)
# ---------------------------------------------------------------------------

test_that("EP-1: one site, perfect test -- Wald CI matches the hand-check", {
  res <- estimate_prevalence(x = 30, n = 100)
  expect_equal(res$prevalence, 0.3)
  expect_equal(res$ci_lower,   0.210183, tolerance = 1e-4)
  expect_equal(res$ci_upper,   0.389817, tolerance = 1e-4)
  expect_equal(res$moe,        0.089817, tolerance = 1e-4)
  expect_equal(res$n_total,    100)
  expect_equal(res$n_eff,      100)
  expect_equal(res$icc_used,   0)
  expect_equal(res$deff,       1)
  # Typing the default sensitivity = specificity = 1 changes nothing
  res_explicit <- estimate_prevalence(x = 30, n = 100,
                                      sensitivity = 1, specificity = 1)
  expect_equal(res_explicit, res)
})

test_that("EP-2: imperfect test -- Rogan-Gladen applied to the estimate and both CI endpoints", {
  res <- estimate_prevalence(x = 30, n = 100,
                             sensitivity = 0.9, specificity = 0.95)
  expect_equal(res$prevalence, 0.294118, tolerance = 1e-4)
  expect_equal(res$ci_lower,   0.188451, tolerance = 1e-4)
  expect_equal(res$ci_upper,   0.399785, tolerance = 1e-4)
  # The interval is the perfect-test interval divided by the correction (0.85)
  expect_equal(res$moe, estimate_prevalence(x = 30, n = 100)$moe / 0.85)
})

test_that("EP-3: sensitivity 0.5 with a perfect specificity doubles the estimate", {
  # correction = 0.5 + 1 - 1 = 0.5, so p_true = p_hat / 0.5
  # p_hat = 15 / 100 = 0.15 -> p_true = 0.30
  # apparent moe = 1.959964 * sqrt(0.15 * 0.85 / 100) = 0.069985 -> / 0.5 = 0.139969
  res <- estimate_prevalence(x = 15, n = 100, sensitivity = 0.5)
  expect_equal(res$prevalence, 0.3)
  expect_equal(res$moe,        0.139969, tolerance = 1e-4)
})

test_that("EP-4: a higher conf_level gives a wider interval", {
  # 99%: z = 2.575829 -> moe = 2.575829 * 0.045826 = 0.118039
  # 90%: z = 1.644854 -> moe = 1.644854 * 0.045826 = 0.075377
  expect_equal(estimate_prevalence(x = 30, n = 100, conf_level = 0.99)$moe,
               0.118039, tolerance = 1e-4)
  expect_equal(estimate_prevalence(x = 30, n = 100, conf_level = 0.90)$moe,
               0.075377, tolerance = 1e-4)
})

test_that("EP-5: no positives, or all positive, collapses the interval to a point", {
  # p_hat = 0 or 1 -> p_hat * (1 - p_hat) = 0 -> the Wald interval has no width
  one_site <- estimate_prevalence(x = 0, n = 100)
  expect_equal(c(one_site$prevalence, one_site$ci_lower, one_site$ci_upper, one_site$moe),
               c(0, 0, 0, 0))

  # Three sites, all zero: no difference between sites, so icc = 0, deff = 1
  zeros <- estimate_prevalence(x = c(0, 0, 0), n = c(10, 10, 10))
  expect_equal(c(zeros$prevalence, zeros$ci_lower, zeros$ci_upper, zeros$moe),
               c(0, 0, 0, 0))
  expect_equal(zeros$icc_used, 0)
  expect_equal(zeros$deff,     1)

  all_pos <- estimate_prevalence(x = c(10, 10, 10), n = c(10, 10, 10))
  expect_equal(c(all_pos$prevalence, all_pos$ci_lower, all_pos$ci_upper, all_pos$moe),
               c(1, 1, 1, 0))
})

test_that("EP-6: whole-number inputs typed as integers (100L) give the same result as doubles", {
  # moe = 1.959964 * sqrt(0.1 * 0.9 / 1000) = 0.018594
  r_int <- estimate_prevalence(x = 100L, n = 1000L)
  r_dbl <- estimate_prevalence(x = 100,  n = 1000)
  expect_equal(r_int, r_dbl)
  expect_equal(r_int$moe, 0.018594, tolerance = 1e-4)
})

test_that("EP-7: a large single site gives a narrow interval", {
  # moe = 1.959964 * sqrt(0.3 * 0.7 / 1000) = 0.028403
  res <- estimate_prevalence(x = 300, n = 1000)
  expect_equal(res$prevalence, 0.3)
  expect_equal(res$moe,        0.028403, tolerance = 1e-4)
})

test_that("EP-8: a Wald interval that runs past 1 is cut at 1, and the message says it is lopsided", {
  # p_hat = 0.99, moe_app = 1.959964 * sqrt(0.99 * 0.01 / 100) = 0.019501
  # CI = [0.970499, 1.009501] -> upper end cut to 1
  expect_message(
    res <- estimate_prevalence(x = 99, n = 100),
    "wald CI is asymmetric: moe_lower = 0.0195, moe_upper = 0.01. moe = 0.0148 is the average half-width; report moe_lower and moe_upper separately.",
    fixed = TRUE
  )
  expect_equal(res$ci_lower,  0.970499, tolerance = 1e-4)
  expect_equal(res$ci_upper,  1)
  expect_equal(res$moe_lower, 0.019501, tolerance = 1e-4)
  expect_equal(res$moe_upper, 0.01,     tolerance = 1e-4)
  expect_equal(res$moe,       0.014751, tolerance = 1e-4)  # average of the two
})

test_that("EP-9: very low prevalence -- the lower end of the interval is cut at 0", {
  # 1 positive out of 50 over 5 sites: p_hat = 0.02. The icc estimated from
  # the data is small (0.0023, deff 1.02), so n_eff = 49 rather than 50.
  res <- suppressMessages(
    estimate_prevalence(x = c(0, 0, 0, 0, 1), n = rep(10, 5))
  )
  expect_equal(res$prevalence, 0.02)
  expect_equal(res$ci_lower,   0)
  expect_equal(res$ci_upper,   0.059199, tolerance = 1e-4)
})

test_that("EP-10: Rogan-Gladen overshoot warns, except when the measured interval already has no width", {
  # p_hat = 0.98 with sensitivity 0.8, specificity 0.9: the whole measured
  # interval maps above 1, so the estimate and CI are all pinned to 1.
  # moe = 0 here is not real precision, so the function warns.
  expect_warning(
    res <- estimate_prevalence(x = 49, n = 50, sensitivity = 0.8, specificity = 0.9),
    "Rogan-Gladen overshoot: after correcting for test error, the estimate and both ends of the confidence interval are all cut to 1, because the corrected interval lies entirely outside 0 to 1 for this sensitivity and specificity.",
    fixed = TRUE
  )
  expect_equal(res$prevalence, 1)
  expect_equal(res$moe,        0)

  # 1 positive out of 1: the measured interval is already a single point
  # ([1, 1]) before correcting, so there is nothing to overshoot -- no warning.
  expect_silent(res1 <- estimate_prevalence(x = 1, n = 1,
                                            sensitivity = 0.8, specificity = 0.9))
  expect_equal(res1$prevalence, 1)   # (1 - 0.1) / 0.7 = 1.29, cut to 1

  # A perfect test with x = 0 is the ordinary [0, 0] case -- no warning
  expect_silent(estimate_prevalence(x = 0, n = 50))
})


# ---------------------------------------------------------------------------
# Clustering: icc and design effect (EP-C-)
# ---------------------------------------------------------------------------

test_that("EP-C-1: a supplied icc widens the interval through the design effect", {
  res <- estimate_prevalence(x = rep(3, 10), n = rep(10, 10), icc = 0.05)
  expect_equal(res$prevalence, 0.3)
  expect_equal(res$deff,       1.45)
  expect_equal(res$n_eff,      100 / 1.45)
  expect_equal(res$ci_lower,   0.191846, tolerance = 1e-4)
  expect_equal(res$ci_upper,   0.408154, tolerance = 1e-4)
  # wider than the same 100 people with no clustering (EP-1: 0.089817)
  expect_gt(res$moe, estimate_prevalence(x = 30, n = 100)$moe)
})

test_that("EP-C-2: icc estimated from the data (the workshop's 8-site example)", {
  # Same sites as Module 5 "The Design Effect - worked example".
  # Observed variance of the site prevalences = 0.042917.
  # Expected variance if people were independent, using the overall
  # prevalence 87 / 550 = 0.158182: mean(0.158182 * 0.841818 / n) = 0.002101.
  # deff = 0.042917 / 0.002101 = 20.43; icc = (20.43 - 1) / (68.75 - 1) = 0.2867
  # NOTE: the slide gets 17.73 because it uses the average of the site
  # prevalences (0.186) instead of the overall prevalence -- see the
  # TODO(review) in R/estimate_prevalence.R. This test pins the current method.
  res <- estimate_prevalence(
    x = c(0, 4, 0, 22, 25, 16, 12, 8),
    n = c(60, 80, 70, 100, 40, 60, 50, 90)
  )
  expect_equal(res$n_total,    550)
  expect_equal(res$prevalence, 0.158182, tolerance = 1e-4)
  expect_equal(res$deff,       20.425830, tolerance = 1e-4)
  expect_equal(res$icc_used,   0.286728,  tolerance = 1e-4)
  expect_equal(res$n_eff,      26.926691, tolerance = 1e-4)
  expect_equal(res$ci_lower,   0.020352,  tolerance = 1e-4)
  expect_equal(res$ci_upper,   0.296012,  tolerance = 1e-4)
})

test_that("EP-C-3: icc = 0, or icc = NULL with one site, gives the no-clustering result", {
  # With one site there is nothing to estimate icc from, so NULL falls back to 0
  srs <- estimate_prevalence(x = 30, n = 100)
  for (res in list(estimate_prevalence(x = 30, n = 100, icc = 0),
                   estimate_prevalence(x = 30, n = 100, icc = NULL))) {
    expect_equal(res$icc_used, 0)
    expect_equal(res$deff,     1)
    expect_equal(res$moe,      srs$moe)   # 0.089817, as in EP-1
  }
})

test_that("EP-C-4: sites of 1 person each (or a single person) give icc = 0, deff = 1", {
  # icc = (deff - 1) / (n_bar - 1) would divide by 0 when every site has
  # 1 person, so the function falls back to no clustering.
  res <- suppressMessages(
    estimate_prevalence(x = c(0, 1, 0, 1, 1), n = c(1, 1, 1, 1, 1))
  )
  expect_equal(res$icc_used,   0)
  expect_equal(res$deff,       1)
  expect_equal(res$prevalence, 0.6)
  expect_equal(res$ci_lower,   0.170593, tolerance = 1e-4)
  expect_equal(res$ci_upper,   1)

  one <- estimate_prevalence(x = 1, n = 1)
  expect_equal(one$icc_used,   0)
  expect_equal(one$deff,       1)
  expect_equal(one$prevalence, 1)
  expect_equal(one$moe,        0)
})

test_that("EP-C-5: extreme site differences cap icc at 1, and deff is recomputed to match", {
  # Two sites all negative, two all positive: the raw icc would be above 1,
  # so it is capped at 1 and deff = 1 + (10 - 1) * 1 = 10.
  res <- estimate_prevalence(x = c(0, 0, 10, 10), n = rep(10, 4))
  expect_equal(res$icc_used, 1)
  expect_equal(res$deff,     10)
  expect_equal(res$n_eff,    4)
  expect_equal(res$ci_lower, 0.010009, tolerance = 1e-4)
  expect_equal(res$ci_upper, 0.989991, tolerance = 1e-4)

  # Two sites, one all negative, one all positive: same cap
  res2 <- estimate_prevalence(x = c(0, 10), n = c(10, 10))
  expect_equal(res2$icc_used,   1)
  expect_equal(res2$deff,       10)
  expect_equal(res2$prevalence, 0.5)
  expect_equal(c(res2$ci_lower, res2$ci_upper), c(0, 1))
})

test_that("EP-C-6: a supplied icc = 1 makes deff equal to the site size", {
  # deff = 1 + (10 - 1) * 1 = 10; n_eff = 100 / 10 = 10 (one per site)
  # moe = 1.959964 * sqrt(0.3 * 0.7 / 10) = 0.284026
  res <- estimate_prevalence(x = rep(3, 10), n = rep(10, 10), icc = 1)
  expect_equal(res$deff,  10)
  expect_equal(res$n_eff, 10)
  expect_equal(res$moe,   0.284026, tolerance = 1e-4)
})

test_that("EP-C-7: sites with identical prevalence give an estimated icc of 0", {
  # No difference between sites -> observed variance 0 -> deff = 1
  # moe = 1.959964 * sqrt(0.3 * 0.7 / 80) = 0.100418
  res <- estimate_prevalence(x = rep(3, 8), n = rep(10, 8))
  expect_equal(res$icc_used, 0)
  expect_equal(res$deff,     1)
  expect_equal(res$moe,      0.100418, tolerance = 1e-4)
})

test_that("EP-C-8: two sites are enough to estimate icc", {
  # Site prevalences 0.2 and 0.4: observed variance 0.02 is below the
  # 0.021 expected for independent people, so deff is set to 1, icc to 0.
  # moe = 1.959964 * sqrt(0.3 * 0.7 / 20) = 0.200837
  res <- estimate_prevalence(x = c(2, 4), n = c(10, 10))
  expect_equal(res$icc_used, 0)
  expect_equal(res$deff,     1)
  expect_equal(res$moe,      0.200837, tolerance = 1e-4)
})

test_that("EP-C-9: sites of different sizes are handled", {
  # Sizes 5, 10, 15, 20: p_hat = 10 / 50 = 0.2
  res <- estimate_prevalence(x = c(1, 2, 3, 4), n = c(5, 10, 15, 20))
  expect_equal(res$n_total,    50)
  expect_equal(res$prevalence, 0.2)
  expect_equal(res$moe,        0.110872, tolerance = 1e-4)

  # Very different sizes (1 and 1000)
  res2 <- estimate_prevalence(x = c(0, 300), n = c(1, 1000))
  expect_equal(res2$n_total,    1001)
  expect_equal(res2$prevalence, 300 / 1001)
  expect_equal(res2$moe,        0.028380, tolerance = 1e-4)
})

test_that("EP-C-10: 100 sites (simulated, fixed seed) give the expected icc and deff", {
  set.seed(42)
  x100 <- rbinom(100, 20, 0.2)          # 414 positives in total
  res  <- estimate_prevalence(x = x100, n = rep(20, 100))
  expect_equal(res$n_total,    2000)
  expect_equal(res$prevalence, 0.207)
  expect_equal(res$icc_used,   0.008910, tolerance = 1e-4)
  expect_equal(res$deff,       1.169286, tolerance = 1e-4)
  expect_equal(res$moe,        0.019201, tolerance = 1e-4)
})

test_that("EP-C-11: a supplied icc with one site is ignored (with a warning); a tiny icc counts as 0", {
  srs <- estimate_prevalence(x = 8, n = 50)

  # One site: there is no cluster structure, so icc = 0.05 cannot be used
  expect_warning(
    res <- estimate_prevalence(x = 8, n = 50, icc = 0.05),
    "`icc` = 0.05 was ignored: the design effect needs a cluster structure (>= 2 clusters, mean size > 1). `icc_used` is reported as 0.",
    fixed = TRUE
  )
  expect_equal(res$icc_used, 0)
  expect_equal(res$deff,     1)
  expect_equal(res$moe,      srs$moe)

  # icc = 0 is just "no clustering" -- nothing to warn about
  expect_silent(estimate_prevalence(x = 8, n = 50, icc = 0))

  # Any icc below about 0.000000015 counts as 0. A value like 1e-12 is
  # almost always leftover rounding from another calculation, not real
  # clustering: no warning with one site, and deff stays exactly 1 with ten.
  expect_silent(estimate_prevalence(x = 8, n = 50, icc = 1e-12))
  tiny <- estimate_prevalence(x = rep(3, 10), n = rep(10, 10), icc = 1e-12)
  expect_equal(tiny$icc_used, 0)
  expect_equal(tiny$deff,     1)
})


# ---------------------------------------------------------------------------
# Finite population correction (EP-F-)
# ---------------------------------------------------------------------------

test_that("EP-F-1: a very large population makes the FPC negligible", {
  # FPC factor = (1e8 - 100) / (1e8 - 1) = 0.99999901 -> n_eff_adj = 100.000099
  # (checked against the exact formula: 100.000099 is too close to 100 for a
  # rounded value to tell "FPC applied" from "FPC not applied")
  res <- estimate_prevalence(x = 30, n = 100, fpc_N = 1e8)
  expect_equal(res$n_eff_adj, 100 / ((1e8 - 100) / (1e8 - 1)))
  expect_gt(res$n_eff_adj, 100)
  expect_equal(res$moe,       0.089817, tolerance = 1e-4)   # same as EP-1
})

test_that("EP-F-2: testing almost the whole population makes the interval very narrow", {
  # 100 of 101 people tested: FPC factor = (101 - 100) / (101 - 1) = 0.01
  # n_eff_adj = 100 / 0.01 = 10000, so moe = 0.089817 * sqrt(0.01) = 0.008982
  res <- estimate_prevalence(x = 30, n = 100, fpc_N = 101)
  expect_equal(res$n_eff_adj, 10000)
  expect_equal(res$moe,       0.008982, tolerance = 1e-4)
})

test_that("EP-F-3: the FPC uses the number of people actually tested, not n_eff", {
  # 5 sites of 60 (300 people), icc 0.1: deff = 1 + (60 - 1) * 0.1 = 6.9,
  # n_eff = 300 / 6.9 = 43.48. Population 310, so 300 / 310 were tested --
  # nearly a census. FPC factor = (310 - 300) / (310 - 1) = 0.03236.
  # n_eff_adj = 43.48 / 0.03236 = 1343.48; moe falls from 0.140122 to 0.025207.
  x <- c(5, 30, 10, 40, 15)
  n <- rep(60, 5)
  r_no  <- estimate_prevalence(x, n, icc = 0.1)
  r_fpc <- estimate_prevalence(x, n, icc = 0.1, fpc_N = 310)
  expect_equal(r_fpc$n_eff,     r_no$n_eff)          # n_eff itself unchanged
  expect_equal(r_fpc$n_eff_adj, 1343.478261, tolerance = 1e-6)
  expect_equal(r_no$moe,        0.140122, tolerance = 1e-4)
  expect_equal(r_fpc$moe,       0.025207, tolerance = 1e-4)
})

test_that("EP-F-4: clustering and FPC together -- both are applied", {
  # deff = 1.45 (as EP-C-1), n_eff = 68.9655; population 500:
  # FPC factor = (500 - 100) / (500 - 1) = 0.8016 -> n_eff_adj = 86.0345
  # moe = 1.959964 * sqrt(0.3 * 0.7 / 86.0345) = 0.096833 (vs 0.108154 without FPC)
  res <- estimate_prevalence(x = rep(3, 10), n = rep(10, 10), icc = 0.05, fpc_N = 500)
  expect_equal(res$deff,      1.45)
  expect_equal(res$n_eff_adj, 86.034483, tolerance = 1e-6)
  expect_equal(res$moe,       0.096833,  tolerance = 1e-4)
})


# ---------------------------------------------------------------------------
# CI methods (EP-M-)
# ---------------------------------------------------------------------------

test_that("EP-M-1: wald (the default) gives a symmetric interval", {
  res <- estimate_prevalence(x = 30, n = 100)
  expect_equal(res$method,    "wald")
  expect_equal(res$moe_lower, res$moe)
  expect_equal(res$moe_upper, res$moe)
})

test_that("EP-M-2: clopper-pearson matches R's own binom.test, and is lopsided", {
  res <- suppressMessages(
    estimate_prevalence(x = 30, n = 100, method = "clopper-pearson")
  )
  ref <- binom.test(30, 100)$conf.int                 # 0.212406, 0.399815
  expect_equal(res$method,   "clopper-pearson")
  expect_equal(c(res$ci_lower, res$ci_upper), c(ref[1], ref[2]))
  # moe is the average half-width; moe_lower / moe_upper are the two sides
  expect_equal(res$moe,       0.093704, tolerance = 1e-4)
  expect_equal(res$moe_lower, 0.087594, tolerance = 1e-4)
  expect_equal(res$moe_upper, 0.099815, tolerance = 1e-4)
})

test_that("EP-M-3: agresti-coull matches the hand-check, and is lopsided", {
  # n~ = 100 + 1.959964^2 = 103.84; p~ = (30 + 1.959964^2 / 2) / 103.84 = 0.30740
  # moe = 1.959964 * sqrt(0.30740 * 0.69260 / 103.84) = 0.088747
  # CI  = [0.218651, 0.396146] (centred on p~, not on 0.3)
  res <- suppressMessages(
    estimate_prevalence(x = 30, n = 100, method = "agresti-coull")
  )
  expect_equal(res$method,    "agresti-coull")
  expect_equal(res$ci_lower,  0.218651, tolerance = 1e-4)
  expect_equal(res$ci_upper,  0.396146, tolerance = 1e-4)
  expect_equal(res$moe_lower, 0.081349, tolerance = 1e-4)
  expect_equal(res$moe_upper, 0.096146, tolerance = 1e-4)
})

test_that("EP-M-4: the method changes the interval, not the estimate", {
  # 15 / 80 = 0.1875 for all three
  wald <- estimate_prevalence(x = 15, n = 80)
  cp   <- suppressMessages(estimate_prevalence(x = 15, n = 80, method = "clopper-pearson"))
  ac   <- suppressMessages(estimate_prevalence(x = 15, n = 80, method = "agresti-coull"))
  expect_equal(c(wald$prevalence, cp$prevalence, ac$prevalence), rep(0.1875, 3))
  expect_equal(c(wald$ci_lower, wald$ci_upper), c(0.101971, 0.273029), tolerance = 1e-4)
  expect_equal(c(cp$ci_lower,   cp$ci_upper),   c(0.108914, 0.290328), tolerance = 1e-4)
  expect_equal(c(ac$ci_lower,   ac$ci_upper),   c(0.115907, 0.287729), tolerance = 1e-4)
})

test_that("EP-M-5: clopper-pearson and agresti-coull at 0 and at all-positive", {
  cp0 <- suppressMessages(estimate_prevalence(x = 0,  n = 50, method = "clopper-pearson"))
  expect_equal(c(cp0$ci_lower, cp0$ci_upper), c(0, 0.071122), tolerance = 1e-4)

  cpn <- suppressMessages(estimate_prevalence(x = 50, n = 50, method = "clopper-pearson"))
  expect_equal(c(cpn$ci_lower, cpn$ci_upper), c(0.928878, 1), tolerance = 1e-4)

  ac0 <- suppressMessages(estimate_prevalence(x = 0,  n = 50, method = "agresti-coull"))
  expect_equal(c(ac0$ci_lower, ac0$ci_upper), c(0, 0.085216), tolerance = 1e-4)
})

test_that("EP-M-6: clopper-pearson with an imperfect test corrects both endpoints", {
  # (0.212406 - 0.05) / 0.85 = 0.191066; (0.399815 - 0.05) / 0.85 = 0.411547
  res <- suppressMessages(
    estimate_prevalence(x = 30, n = 100, sensitivity = 0.9, specificity = 0.95,
                        method = "clopper-pearson")
  )
  expect_equal(res$prevalence, 0.294118, tolerance = 1e-4)
  expect_equal(res$ci_lower,   0.191066, tolerance = 1e-4)
  expect_equal(res$ci_upper,   0.411547, tolerance = 1e-4)
})

test_that("EP-M-7: clopper-pearson with clustering gives a wider interval", {
  # icc 0.05 -> deff 1.45 -> CP computed on n_eff = 68.97 instead of 100
  srs <- suppressMessages(estimate_prevalence(x = rep(3, 10), n = rep(10, 10),
                                              icc = 0,    method = "clopper-pearson"))
  cl  <- suppressMessages(estimate_prevalence(x = rep(3, 10), n = rep(10, 10),
                                              icc = 0.05, method = "clopper-pearson"))
  expect_equal(c(srs$ci_lower, srs$ci_upper), c(0.212406, 0.399815), tolerance = 1e-4)
  expect_equal(c(cl$ci_lower,  cl$ci_upper),  c(0.195497, 0.422331), tolerance = 1e-4)
})

test_that("EP-M-8: the lopsided-interval message appears only when the two sides differ by > 10% of moe", {
  expect_message(
    estimate_prevalence(x = 30, n = 100, method = "clopper-pearson"),
    "clopper-pearson CI is asymmetric: moe_lower = 0.0876, moe_upper = 0.0998. moe = 0.0937 is the average half-width; report moe_lower and moe_upper separately.",
    fixed = TRUE
  )
  expect_message(
    estimate_prevalence(x = 30, n = 100, method = "agresti-coull"),
    "agresti-coull CI is asymmetric: moe_lower = 0.0813, moe_upper = 0.0961. moe = 0.0887 is the average half-width; report moe_lower and moe_upper separately.",
    fixed = TRUE
  )
  # No message: a wald interval away from 0 and 1, and a clopper-pearson
  # interval at p = 0.5 with n = 1000 (sides 0.031451 each)
  expect_no_message(estimate_prevalence(x = 30,  n = 100))
  expect_no_message(estimate_prevalence(x = 500, n = 1000, method = "clopper-pearson"))
})


# ---------------------------------------------------------------------------
# Validation (EP-V-)
# ---------------------------------------------------------------------------

test_that("EP-V-1: TRUE/FALSE or a bare NA in x or n is rejected", {
  # Valid x and n: every test from EP-1 on (e.g. x = 30, n = 100).
  # A bare NA is logical in R, so it is caught here, and the message says so
  expect_error(estimate_prevalence(x = NA, n = 100),
               "`x` and `n` must be numeric, not logical (got class `logical` for x, `numeric` for n). Note: a plain `NA` is logical in R -- remove missing observations before calling.",
               fixed = TRUE)
  expect_error(estimate_prevalence(x = 30, n = TRUE),
               "`x` and `n` must be numeric, not logical (got class `numeric` for x, `logical` for n).",
               fixed = TRUE)
})

test_that("EP-V-2: text in x or n is rejected, naming the class", {
  # Valid x and n: every test from EP-1 on (e.g. x = 30, n = 100).
  expect_error(estimate_prevalence(x = "30", n = 100),
               "`x` and `n` must be numeric vectors (got class `character` for x, `numeric` for n).",
               fixed = TRUE)
  expect_error(estimate_prevalence(x = 30, n = "100"),
               "`x` and `n` must be numeric vectors (got class `numeric` for x, `character` for n).",
               fixed = TRUE)
})

test_that("EP-V-3: empty x and n are rejected", {
  # The smallest valid input is one site: x = 1, n = 1 (EP-C-4).
  expect_error(estimate_prevalence(x = numeric(0), n = numeric(0)),
               "`x` and `n` must be non-empty vectors.", fixed = TRUE)
})

test_that("EP-V-4: NA, NaN or Inf inside x or n is rejected", {
  # Valid x and n: every test from EP-1 on (e.g. x = 30, n = 100).
  expect_error(estimate_prevalence(x = c(1, NA_real_), n = c(10, 10)),
               "`x` contains NA, NaN, or infinite values.", fixed = TRUE)
  expect_error(estimate_prevalence(x = c(1, 2, 3), n = c(10, NA, 10)),
               "`n` contains NA, NaN, or infinite values.", fixed = TRUE)
  expect_error(estimate_prevalence(x = c(1, 2, 3), n = c(10, Inf, 10)),
               "`n` contains NA, NaN, or infinite values.", fixed = TRUE)
})

test_that("EP-V-5: x and n of different lengths are rejected", {
  # Equal lengths are valid: e.g. 10 sites in EP-C-1.
  expect_error(estimate_prevalence(x = c(3, 5), n = 10),
               "`x` and `n` must have the same length (got length(x) = 2, length(n) = 1).",
               fixed = TRUE)
  expect_error(estimate_prevalence(x = c(1, 2, 3), n = c(10, 10)),
               "`x` and `n` must have the same length (got length(x) = 3, length(n) = 2).",
               fixed = TRUE)
})

test_that("EP-V-6: a negative count in x is rejected, naming its position", {
  # x = 0 is allowed (no positives): see EP-5.
  expect_error(estimate_prevalence(x = c(1, -1, 2), n = c(10, 10, 10)),
               "`x` must be non-negative (found x[2] = -1).", fixed = TRUE)
  # -0.5 is both negative and not whole: the negative check comes first
  expect_error(estimate_prevalence(x = -0.5, n = 10),
               "`x` must be non-negative (found x[1] = -0.5).", fixed = TRUE)
})

test_that("EP-V-7: a zero or negative n is rejected, naming its position", {
  # n = 1 is allowed (a site of one person): see EP-C-4.
  expect_error(estimate_prevalence(x = c(1, 0, 2), n = c(10, 0, 10)),
               "`n` must be positive for every cluster (found n[2] = 0).", fixed = TRUE)
  expect_error(estimate_prevalence(x = 1, n = -5),
               "`n` must be positive for every cluster (found n[1] = -5).", fixed = TRUE)
})

test_that("EP-V-8: a fraction in x or n is rejected, naming its position", {
  # Whole numbers stored as decimals (100) or integers (100L) both work: see EP-6.
  expect_error(estimate_prevalence(x = c(1.5, 2, 3), n = c(10, 10, 10)),
               "`x` must contain whole numbers -- counts cannot be fractional (found x[1] = 1.5).",
               fixed = TRUE)
  expect_error(estimate_prevalence(x = c(1, 2, 3), n = c(10, 10.5, 10)),
               "`n` must contain whole numbers -- sample sizes cannot be fractional (found n[2] = 10.5).",
               fixed = TRUE)
})

test_that("EP-V-9: more positives than people tested is rejected, naming the position", {
  expect_error(estimate_prevalence(x = 60, n = 50),
               "`x` cannot exceed `n` (found x[1] = 60 > n[1] = 50).", fixed = TRUE)
  expect_error(estimate_prevalence(x = c(3, 12), n = c(10, 10)),
               "`x` cannot exceed `n` (found x[2] = 12 > n[2] = 10).", fixed = TRUE)
  # x = n is allowed (everyone positive): see EP-5
})

test_that("EP-V-10: a vector in a single-number argument is rejected, giving its length", {
  expect_error(estimate_prevalence(x = 30, n = 100, sensitivity = c(0.9, 0.8)),
               "`sensitivity` must be a single number in (0, 1] (got class `numeric`, length 2).",
               fixed = TRUE)
  expect_error(estimate_prevalence(x = 30, n = 100, conf_level = c(0.9, 0.95)),
               "`conf_level` must be a single number in (0, 1) (got class `numeric`, length 2).",
               fixed = TRUE)
  expect_error(estimate_prevalence(x = c(3, 3), n = c(10, 10), icc = c(0.1, 0.2)),
               "`icc` must be a single number in [0, 1] (got class `numeric`, length 2). To estimate ICC from the data, leave `icc = NULL`.",
               fixed = TRUE)
})

test_that("EP-V-11: text, TRUE/FALSE, a bare NA, a list or NULL is rejected, naming the class", {
  # None of these is a number, so each is stopped by the type check and the
  # message says what was passed. A bare NA counts as TRUE/FALSE-type
  # (logical) in R, so it gets the same "class `logical`" message.
  expect_error(estimate_prevalence(x = 30, n = 100, sensitivity = "0.9"),
               "`sensitivity` must be a single number in (0, 1] (got class `character`, length 1).",
               fixed = TRUE)
  # TRUE must not be quietly read as a perfect test (1); the message says so.
  expect_error(estimate_prevalence(x = 30, n = 100, sensitivity = TRUE),
               "`sensitivity` must be a single number in (0, 1] (got class `logical`, length 1). Note: `TRUE`/`FALSE` is logical, not numeric -- pass 1 for a perfect test.",
               fixed = TRUE)
  expect_error(estimate_prevalence(x = 30, n = 100, sensitivity = NULL),
               "`sensitivity` must be a single number in (0, 1] (got class `NULL`, length 0).",
               fixed = TRUE)
  expect_error(estimate_prevalence(x = 30, n = 100, sensitivity = list(0.9)),
               "`sensitivity` must be a single number in (0, 1] (got class `list`, length 1).",
               fixed = TRUE)
  expect_error(estimate_prevalence(x = 30, n = 100, specificity = FALSE),
               "`specificity` must be a single number in (0, 1] (got class `logical`, length 1).",
               fixed = TRUE)
  expect_error(estimate_prevalence(x = 30, n = 100, specificity = NA),
               "`specificity` must be a single number in (0, 1] (got class `logical`, length 1).",
               fixed = TRUE)
  expect_error(estimate_prevalence(x = c(3, 3), n = c(10, 10), icc = FALSE),
               "`icc` must be a single number in [0, 1] (got class `logical`, length 1).",
               fixed = TRUE)
  expect_error(estimate_prevalence(x = 30, n = 100, conf_level = NA),
               "`conf_level` must be a single number in (0, 1) (got class `logical`, length 1).",
               fixed = TRUE)
  expect_error(estimate_prevalence(x = 30, n = 100, conf_level = list(0.95)),
               "`conf_level` must be a single number in (0, 1) (got class `list`, length 1).",
               fixed = TRUE)
})

test_that("EP-V-12: NaN or Inf in a single-number argument gives the finite-value error", {
  # NaN and Inf are numbers, so they get through the type check and reach this one
  expect_error(estimate_prevalence(x = 30, n = 100, sensitivity = NaN),
               "`sensitivity` must be a single finite number (got NaN).", fixed = TRUE)
  expect_error(estimate_prevalence(x = 30, n = 100, specificity = Inf),
               "`specificity` must be a single finite number (got Inf).", fixed = TRUE)
  expect_error(estimate_prevalence(x = 30, n = 100, conf_level = NaN),
               "`conf_level` must be a single finite number (got NaN).", fixed = TRUE)
  expect_error(estimate_prevalence(x = c(3, 3), n = c(10, 10), icc = Inf),
               "`icc` must be a single finite number (got Inf). To estimate ICC from the data, leave `icc = NULL`.",
               fixed = TRUE)
})

test_that("EP-V-13: sensitivity's own range check (0, 1] gives the right message for each bad value", {
  # Valid values: 1 (EP-1), 0.9 (EP-2), 0.5 (EP-3).
  # sensitivity = 0 with the default specificity = 1 is the "0 + 1 = 1"
  # case: it is stopped here, by sensitivity's own check, before the
  # combined check runs.
  expect_error(estimate_prevalence(x = 30, n = 100, sensitivity = 0),
               "A sensitivity of 0 means", fixed = TRUE)
  expect_error(estimate_prevalence(x = 30, n = 100, sensitivity = 1.1),
               "`sensitivity` is a diagnostic probability and cannot be greater than 1",
               fixed = TRUE)
  expect_error(estimate_prevalence(x = 30, n = 100, sensitivity = -0.5),
               "`sensitivity` is a diagnostic probability and cannot be negative",
               fixed = TRUE)
})

test_that("EP-V-14: specificity's own range check (0, 1] gives the right message for each bad value", {
  # Valid values: 1 (EP-1), 0.95 (EP-2).
  # Mirror of the sensitivity test above.
  expect_error(estimate_prevalence(x = 30, n = 100, specificity = 0),
               "A specificity of 0 means", fixed = TRUE)
  expect_error(estimate_prevalence(x = 30, n = 100, specificity = 1.1),
               "`specificity` is a diagnostic probability and cannot be greater than 1",
               fixed = TRUE)
  expect_error(estimate_prevalence(x = 30, n = 100, specificity = -0.5),
               "`specificity` is a diagnostic probability and cannot be negative",
               fixed = TRUE)
})

test_that("EP-V-15: sensitivity + specificity <= 1 is rejected (Rogan-Gladen correction)", {
  # A sum just above 1 is accepted, with a warning: see EP-V-16.
  # The correction divides by (sensitivity + specificity - 1), so that
  # value must be above 0.
  # Sum below 1 (0.3 + 0.3 = 0.6): it would divide by a negative number.
  expect_error(estimate_prevalence(x = 10, n = 50, sensitivity = 0.3, specificity = 0.3),
               "(got 0.3 + 0.3 = 0.6). The correction divides by `sensitivity` + `specificity` - 1, so a sum of 1 or less cannot be corrected.",
               fixed = TRUE)
  # Sum exactly 1 (0.5 + 0.5): it would divide by 0.
  expect_error(estimate_prevalence(x = 10, n = 50, sensitivity = 0.5, specificity = 0.5),
               "(got 0.5 + 0.5 = 1).", fixed = TRUE)
})

test_that("EP-V-16: sensitivity + specificity above 1 but below 1.1 warns, 1.1 and above does not", {
  # A sum just above 1 can be corrected, but the correction divides by a
  # very small number, so the interval becomes very wide: 55 / 100 positive
  # with sensitivity 0.59, specificity 0.5 gives a CI of [0, 1].
  expect_warning(
    res <- suppressMessages(
      estimate_prevalence(x = 55, n = 100, sensitivity = 0.59, specificity = 0.5)
    ),
    "`sensitivity` + `specificity` = 1.09, which is very close to 1. Correcting for this much test error makes the prevalence estimate and its confidence interval very wide",
    fixed = TRUE
  )
  expect_equal(c(res$ci_lower, res$ci_upper), c(0, 1))
  # 1.10 does not warn: 0.6 + 0.5 - 1 comes out as 0.10000000000000009 in
  # computer arithmetic, just above the 0.1 cut-off. 1.11 does not warn.
  expect_no_warning(suppressMessages(
    estimate_prevalence(x = 55, n = 100, sensitivity = 0.60, specificity = 0.5)))
  expect_no_warning(suppressMessages(
    estimate_prevalence(x = 55, n = 100, sensitivity = 0.61, specificity = 0.5)))
})

test_that("EP-V-17: icc must be in [0, 1] -- below 0 and above 1 are rejected, 1 works", {
  expect_error(estimate_prevalence(x = c(3, 3), n = c(10, 10), icc = -0.1),
               "`icc` must be in [0, 1] (got -0.1). `icc` is a correlation and cannot be negative. To estimate ICC from the data, leave `icc = NULL`.",
               fixed = TRUE)
  expect_error(estimate_prevalence(x = c(3, 3), n = c(10, 10), icc = 2),
               "`icc` must be in [0, 1] (got 2). `icc` is a correlation and cannot be greater than 1. To estimate ICC from the data, leave `icc = NULL`.",
               fixed = TRUE)
  # icc = 1 itself is allowed: see EP-C-6 (deff = 10, moe = 0.284026)
  expect_equal(estimate_prevalence(x = rep(3, 10), n = rep(10, 10), icc = 1)$deff, 10)
})

test_that("EP-V-18: conf_level must be in (0, 1) -- boundaries and outside values rejected, inside works", {
  for (bad in c(0, 1, -0.05, 1.05)) {
    expect_error(estimate_prevalence(x = 30, n = 100, conf_level = bad),
                 paste0("`conf_level` must be in (0, 1) (got ", bad, "). For example, use 0.95 for a 95% confidence interval."),
                 fixed = TRUE)
  }
  # A valid value gives the hand-checked moe (see EP-4): 0.075377
  expect_equal(estimate_prevalence(x = 30, n = 100, conf_level = 0.90)$moe,
               0.075377, tolerance = 1e-4)
})

test_that("EP-V-19: fpc_N must be a single positive whole number", {
  for (bad in list(0, -10, 100.5, Inf)) {
    expect_error(estimate_prevalence(x = 30, n = 100, fpc_N = bad),
                 paste0("`fpc_N` must be a single finite positive integer (got ", bad, ")."),
                 fixed = TRUE)
  }
  # A string or TRUE is reported by its class; a vector by its length.
  # TRUE must not be quietly read as fpc_N = 1.
  expect_error(estimate_prevalence(x = 30, n = 100, fpc_N = "500"),
               "`fpc_N` must be a single finite positive integer (got class `character`).",
               fixed = TRUE)
  expect_error(estimate_prevalence(x = 30, n = 100, fpc_N = TRUE),
               "`fpc_N` must be a single finite positive integer (got class `logical`).",
               fixed = TRUE)
  expect_error(estimate_prevalence(x = 30, n = 100, fpc_N = c(500, 600)),
               "`fpc_N` must be a single finite positive integer (got length = 2).",
               fixed = TRUE)
  # A valid value: FPC factor = (400 - 100) / (400 - 1) -> n_eff_adj = 133
  expect_equal(estimate_prevalence(x = 30, n = 100, fpc_N = 400)$n_eff_adj, 133)
})

test_that("EP-V-20: fpc_N smaller than, or equal to, the number tested is rejected", {
  # Smaller: more people tested than exist
  expect_error(estimate_prevalence(x = c(3, 3), n = c(10, 10), fpc_N = 5),
               "`fpc_N` (5) is smaller than the total sample size (20): you cannot test more people than there are in the population.",
               fixed = TRUE)
  # Equal: everyone was tested (a census), so there is no interval to compute.
  # Holds for one site and for several (the FPC uses the total tested).
  census_msg <- "equals the total sample size: the whole population was tested (a census), so there is no sampling uncertainty and no confidence interval to compute."
  expect_error(estimate_prevalence(x = 30, n = 100, fpc_N = 100),
               paste("`fpc_N` (100)", census_msg), fixed = TRUE)
  expect_error(estimate_prevalence(x = c(10, 10, 10), n = c(40, 40, 40), fpc_N = 120),
               paste("`fpc_N` (120)", census_msg), fixed = TRUE)
  # One more than the number tested works: see EP-F-2
})

test_that("EP-V-21: method must be one of the three names", {
  expect_error(estimate_prevalence(x = 30, n = 100, method = "exact"),
               "`method` must be one of 'wald', 'clopper-pearson', or 'agresti-coull' (got 'exact').",
               fixed = TRUE)
  expect_error(estimate_prevalence(x = 30, n = 100, method = c("wald", "clopper-pearson")),
               "`method` must be a single character string: 'wald', 'clopper-pearson', or 'agresti-coull' (got class `character`, length 2).",
               fixed = TRUE)
  expect_error(estimate_prevalence(x = 30, n = 100, method = 1),
               "`method` must be a single character string: 'wald', 'clopper-pearson', or 'agresti-coull' (got class `numeric`, length 1).",
               fixed = TRUE)
})


# ---------------------------------------------------------------------------
# Return value (EP-R-)
# ---------------------------------------------------------------------------

test_that("EP-R-1: return list contains all expected fields, in order", {
  res <- estimate_prevalence(x = 30, n = 100)
  expect_named(res, c("prevalence", "ci_lower", "ci_upper", "moe",
                      "moe_lower", "moe_upper", "method",
                      "n_total", "n_eff", "n_eff_adj", "conf_level",
                      "sensitivity", "specificity", "icc_used", "deff", "fpc_N"),
               ignore.order = FALSE)
})

test_that("EP-R-2: n_eff_adj equals n_eff without an FPC, and is larger with one", {
  # n_eff is the equivalent simple-random-sample size before the FPC;
  # n_eff_adj is the size the interval actually uses after it.
  no_fpc <- estimate_prevalence(x = 30, n = 100)
  expect_equal(no_fpc$n_eff_adj, no_fpc$n_eff)          # both 100

  # Population 400: FPC factor = (400 - 100) / (400 - 1) = 0.7519
  # n_eff_adj = 100 / 0.7519 = 133; moe = 1.959964 * sqrt(0.21 / 133) = 0.077881
  with_fpc <- estimate_prevalence(x = 30, n = 100, fpc_N = 400)
  expect_equal(with_fpc$n_eff,     100)
  expect_equal(with_fpc$n_eff_adj, 133)
  expect_equal(with_fpc$moe,       0.077881, tolerance = 1e-4)
})


# ---------------------------------------------------------------------------
# Round trips with design_precision() (EP-T-)
# ---------------------------------------------------------------------------

test_that("EP-T-1: the n from design_precision gives back the target moe (no clustering)", {
  # design_precision(0.3, 0.05) -> n = 323; observing 30% -> x = round(0.3 * 323) = 97
  res <- estimate_prevalence(x = 97, n = 323)
  expect_equal(res$moe, 0.049990, tolerance = 1e-4)   # target 0.05
})

test_that("EP-T-2: the clustered design from design_precision gives back the target moe when icc is supplied", {
  # design_precision(0.3, 0.05, n_per_site = 10, icc = 0.05) -> 47 sites of 10;
  # 3 positives per site
  res <- estimate_prevalence(x = rep(3, 47), n = rep(10, 47), icc = 0.05)
  expect_equal(res$moe, 0.049888, tolerance = 1e-4)   # target 0.05
})

test_that("EP-T-3: without icc supplied, identical sites estimate icc = 0 and the round trip misses", {
  # Identical site prevalences -> observed variance 0 -> estimated icc = 0,
  # deff = 1 -> the interval is narrower than design_precision planned for.
  # The clustered round trip only holds when icc is supplied (EP-T-2).
  res <- estimate_prevalence(x = rep(3, 47), n = rep(10, 47))
  expect_equal(res$icc_used, 0)
  expect_equal(res$deff,     1)
  expect_equal(res$moe,      0.041429, tolerance = 1e-4)   # not 0.05
})
