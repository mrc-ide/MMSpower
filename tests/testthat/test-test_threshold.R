# Hand-checked reference cases for test_threshold()
#
# Formula: z = (p_hat - theta_app) / sqrt(theta_app*(1-theta_app) / n_eff_adj)
#   theta_app = threshold * se + (1-threshold)*(1-sp)   [threshold on apparent scale]
#   n_eff_adj = n_total / deff / fpc^2                  [variance-equivalent n]
#   z = qnorm(0.975) = 1.959964 for 95% conf_level
#
# Case 1  (single site, x=8, n=100, threshold=0.05, alternative="greater"):
#   p_hat = 0.08, theta_app = 0.05 (perfect test)
#   n_eff_adj = 100, se_null = sqrt(0.05*0.95/100) = sqrt(0.000475) = 0.02179
#   z = (0.08 - 0.05) / 0.02179 = 1.376
#   p_value = 1 - pnorm(1.376) = 0.0844 → fail to reject at alpha=0.05
#
# Case 2  (x=15, n=100, threshold=0.05, alternative="greater"):
#   p_hat = 0.15, theta_app = 0.05
#   z = (0.15 - 0.05) / 0.02179 = 4.588
#   p_value = 1 - pnorm(4.588) ≈ 2.24e-6 → reject
#
# Case 3  (two.sided, x=8, n=100, threshold=0.08):
#   p_hat = 0.08 = threshold → z = 0 → p_value = 1 → fail to reject
#
# Case 4  (imperfect test, se=0.9, sp=0.95, x=8, n=100, threshold=0.05):
#   theta_app = 0.05*0.9 + 0.95*0.05 = 0.045 + 0.0475 = 0.0925
#   p_hat = 0.08
#   se_null = sqrt(0.0925*0.9075/100) = sqrt(0.0000839) = 0.009158... hmm wait
#   se_null = sqrt(0.0925*(1-0.0925)/100) = sqrt(0.0925*0.9075/100)
#           = sqrt(0.083944/100) = sqrt(0.00083944) = 0.028973
#   z = (0.08 - 0.0925) / 0.028973 = -0.01250 / 0.028973 = -0.4314
#   p_value (greater) = 1 - pnorm(-0.4314) = 0.667 → fail to reject
#
# Case 5  (clustered, x=rep(3,10), n=rep(10,10), icc=0.05, threshold=0.2):
#   p_hat = 30/100 = 0.3, n_bar=10, deff = 1+(10-1)*0.05 = 1.45
#   n_eff = 100/1.45 = 68.97, n_eff_adj = 68.97 (no FPC)
#   theta_app = 0.2 (perfect test)
#   se_null = sqrt(0.2*0.8/68.97) = sqrt(0.002320) = 0.04817
#   z = (0.3 - 0.2) / 0.04817 = 2.076
#   p_value (greater) = 1 - pnorm(2.076) = 0.0190 → reject at alpha=0.05

test_that("Case 1: x=8, n=100, threshold=0.05 -- fail to reject (p=0.084)", {
  res <- test_threshold(x = 8, n = 100, threshold = 0.05)

  expect_equal(res$threshold,   0.05)
  expect_equal(res$alternative, "greater")
  expect_equal(res$prevalence,  0.08, tolerance = 1e-6)
  expect_equal(res$n_total,     100)
  expect_equal(res$deff,        1,    tolerance = 1e-6)

  # z = (0.08 - 0.05) / sqrt(0.05*0.95/100)
  z_expected <- (0.08 - 0.05) / sqrt(0.05 * 0.95 / 100)
  expect_equal(res$statistic, z_expected, tolerance = 1e-6)
  expect_equal(res$p_value, pnorm(z_expected, lower.tail = FALSE), tolerance = 1e-6)
  expect_false(res$reject)   # p > 0.05
})

test_that("Case 2: x=15, n=100, threshold=0.05 -- reject (large z)", {
  res <- test_threshold(x = 15, n = 100, threshold = 0.05)

  z_expected <- (0.15 - 0.05) / sqrt(0.05 * 0.95 / 100)
  expect_equal(res$statistic, z_expected, tolerance = 1e-6)
  expect_lt(res$p_value, 0.05)
  expect_true(res$reject)
})

test_that("Case 3: p_hat = threshold exactly → z=0, p_value=0.5 (greater)", {
  res <- test_threshold(x = 8, n = 100, threshold = 0.08)

  expect_equal(res$statistic, 0, tolerance = 1e-6)
  expect_equal(res$p_value,   0.5, tolerance = 1e-6)
  expect_false(res$reject)
})

test_that("Case 3b: two.sided with p_hat = threshold → z=0, p_value=1", {
  res <- test_threshold(x = 8, n = 100, threshold = 0.08, alternative = "two.sided")
  expect_equal(res$statistic, 0,   tolerance = 1e-6)
  expect_equal(res$p_value,   1.0, tolerance = 1e-6)
  expect_false(res$reject)
})

test_that("Case 4: imperfect test shifts theta_app, changes z and p-value", {
  res_perfect  <- test_threshold(x = 8, n = 100, threshold = 0.05)
  res_imperfect <- test_threshold(x = 8, n = 100, threshold = 0.05,
                                  sensitivity = 0.9, specificity = 0.95)

  # theta_app with imperfect test = 0.05*0.9 + 0.95*0.05 = 0.0925
  expect_equal(res_imperfect$statistic,
               (0.08 - 0.0925) / sqrt(0.0925 * 0.9075 / 100),
               tolerance = 1e-5)
  # RG-corrected prevalence differs from perfect test
  expect_false(isTRUE(all.equal(res_imperfect$prevalence, res_perfect$prevalence)))
})

test_that("Case 5: clustering (icc=0.05) widens null SE and reduces z", {
  res_srs     <- test_threshold(x = rep(3, 10), n = rep(10, 10),
                                 threshold = 0.2, icc = 0)
  res_clustered <- test_threshold(x = rep(3, 10), n = rep(10, 10),
                                   threshold = 0.2, icc = 0.05)

  # Clustering inflates SE → smaller |z| → larger p-value
  expect_lt(abs(res_clustered$statistic), abs(res_srs$statistic))
  expect_gt(res_clustered$p_value, res_srs$p_value)
  expect_equal(res_clustered$deff, 1.45, tolerance = 1e-6)
})

test_that("alternative='less' reverses the test direction", {
  # x=0, n=100, threshold=0.05: p_hat=0, z strongly negative → reject 'less'
  # z = (0 - 0.05)/sqrt(0.05*0.95/100) = -2.294; p(less) = 0.011 < 0.05
  res_less    <- test_threshold(x = 0, n = 100, threshold = 0.05, alternative = "less")
  res_greater <- test_threshold(x = 0, n = 100, threshold = 0.05, alternative = "greater")

  expect_true(res_less$reject)
  expect_false(res_greater$reject)
  # p-values sum to 1 for one-sided (z<0): pnorm(z) + pnorm(-z) = 1
  expect_equal(res_less$p_value + res_greater$p_value, 1, tolerance = 1e-10)
})

test_that("two.sided p-value = 2 * min(one-sided p-values)", {
  res_ts <- test_threshold(x = 15, n = 100, threshold = 0.05, alternative = "two.sided")
  res_gt <- test_threshold(x = 15, n = 100, threshold = 0.05, alternative = "greater")
  expect_equal(res_ts$p_value, 2 * res_gt$p_value, tolerance = 1e-10)
})

test_that("CI is the two-sided Wald interval (hand-checked, threshold-independent)", {
  # threshold does not affect the CI -- only x, n, conf_level matter
  res <- test_threshold(x = 30, n = 100, threshold = 0.20)

  # Wald CI: p_hat=0.3, n_eff_adj=100, z=1.96
  # moe = 1.96 * sqrt(0.3*0.7/100) = 0.08982
  z_ci      <- qnorm(0.975)
  moe_check <- z_ci * sqrt(0.3 * 0.7 / 100)
  expect_equal(res$prevalence, 0.3,              tolerance = 1e-6)
  expect_equal(res$ci_lower,   0.3 - moe_check,  tolerance = 1e-6)
  expect_equal(res$ci_upper,   0.3 + moe_check,  tolerance = 1e-6)

  # CI does not change when threshold changes (threshold only affects statistic)
  res2 <- test_threshold(x = 30, n = 100, threshold = 0.40)
  expect_equal(res$ci_lower, res2$ci_lower, tolerance = 1e-10)
  expect_equal(res$ci_upper, res2$ci_upper, tolerance = 1e-10)
})

test_that("FPC tightens n_eff and increases |z|", {
  res_nofpc <- test_threshold(x = 8, n = 100, threshold = 0.05)
  res_fpc   <- test_threshold(x = 8, n = 100, threshold = 0.05, fpc_N = 200)

  # FPC makes n_eff_adj larger → narrower null SE → larger |z|
  expect_gt(abs(res_fpc$statistic), abs(res_nofpc$statistic))
  expect_lt(res_fpc$p_value, res_nofpc$p_value)
})

test_that("return list has all expected fields in order", {
  res <- test_threshold(x = 8, n = 100, threshold = 0.05)
  expect_named(res, c("statistic", "p_value", "reject", "threshold", "alternative",
                       "prevalence", "ci_lower", "ci_upper",
                       "n_total", "n_eff", "conf_level",
                       "sensitivity", "specificity", "icc_used", "deff", "fpc_N"),
               ignore.order = FALSE)
})

# ---- input validation ----

test_that("validation: non-numeric x gives informative error", {
  expect_error(test_threshold(x = "8", n = 100, threshold = 0.05), "numeric")
})

test_that("validation: x > n gives informative error", {
  expect_error(test_threshold(x = 60, n = 50, threshold = 0.05), "cannot exceed")
})

test_that("validation: threshold=0 is rejected", {
  expect_error(test_threshold(x = 8, n = 100, threshold = 0), "`threshold`")
})

test_that("validation: threshold=1 is rejected", {
  expect_error(test_threshold(x = 8, n = 100, threshold = 1), "`threshold`")
})

test_that("validation: threshold outside (0,1) is rejected", {
  expect_error(test_threshold(x = 8, n = 100, threshold = 1.5), "`threshold`")
  expect_error(test_threshold(x = 8, n = 100, threshold = -0.1), "`threshold`")
})

test_that("validation: invalid alternative is rejected", {
  expect_error(test_threshold(x = 8, n = 100, threshold = 0.05, alternative = "up"),
               "'greater', 'less', or 'two.sided'")
})

test_that("validation: se+sp <= 1 is rejected", {
  expect_error(
    test_threshold(x = 8, n = 100, threshold = 0.05,
                   sensitivity = 0.4, specificity = 0.4),
    "must exceed 1"
  )
})

test_that("validation: icc out of range is rejected", {
  expect_error(test_threshold(x = 8, n = 100, threshold = 0.05, icc = 1.5), "`icc`")
  expect_error(test_threshold(x = 8, n = 100, threshold = 0.05, icc = -0.1), "`icc`")
})

test_that("validation: fpc_N < n_total is rejected", {
  expect_error(
    test_threshold(x = 8, n = 100, threshold = 0.05, fpc_N = 50),
    "at least as large"
  )
})

test_that("validation: conf_level out of range is rejected", {
  expect_error(test_threshold(x = 8, n = 100, threshold = 0.05, conf_level = 0),
               "`conf_level`")
  expect_error(test_threshold(x = 8, n = 100, threshold = 0.05, conf_level = 1),
               "`conf_level`")
})

test_that("validation: threshold as vector is rejected", {
  expect_error(
    test_threshold(x = 8, n = 100, threshold = c(0.05, 0.10)),
    "`threshold`"
  )
})

# ---- edge cases ----

test_that("x=0 (no positives): p_hat=0, z negative (less than threshold)", {
  res <- test_threshold(x = 0, n = 100, threshold = 0.05)
  expect_equal(res$prevalence, 0,    tolerance = 1e-10)
  expect_lt(res$statistic,     0)
  expect_false(res$reject)
})

test_that("x=n (all positive): p_hat=1, z large positive → reject 'greater'", {
  res <- test_threshold(x = 100, n = 100, threshold = 0.05)
  expect_equal(res$prevalence, 1, tolerance = 1e-10)
  expect_gt(res$statistic, 0)
  expect_true(res$reject)
})

test_that("icc=0 explicit is same as icc=NULL with one cluster", {
  res_null <- test_threshold(x = 30, n = 100, threshold = 0.2, icc = NULL)
  res_zero <- test_threshold(x = 30, n = 100, threshold = 0.2, icc = 0)
  expect_equal(res_null$statistic, res_zero$statistic, tolerance = 1e-10)
  expect_equal(res_null$p_value,   res_zero$p_value,   tolerance = 1e-10)
})

test_that("conf_level=0.99 raises the bar for rejection", {
  # p ≈ 0.084 for Case 1 -- rejects at 0.10 but not at 0.05 or 0.01
  res_90 <- test_threshold(x = 8, n = 100, threshold = 0.05, conf_level = 0.90)
  res_99 <- test_threshold(x = 8, n = 100, threshold = 0.05, conf_level = 0.99)
  expect_true(res_90$reject)   # alpha=0.10: 0.084 < 0.10
  expect_false(res_99$reject)  # alpha=0.01: 0.084 > 0.01
})

test_that("n_bar=1 guard: all clusters of size 1 → icc_used=0, deff=1", {
  res <- test_threshold(x = c(0, 1, 0, 1), n = c(1, 1, 1, 1), threshold = 0.3)
  expect_equal(res$icc_used, 0)
  expect_equal(res$deff,     1)
  expect_true(is.finite(res$statistic))
})

test_that("ICC estimated from heterogeneous clusters inflates deff and reduces n_eff", {
  # Identical clusters → zero between-cluster variance → deff clamped to 1
  # Highly heterogeneous clusters → large ICC → deff >> 1 → much smaller n_eff
  res_identical <- test_threshold(
    x = c(5, 5, 5, 5), n = c(50, 50, 50, 50), threshold = 0.05
  )
  res_heterog <- test_threshold(
    x = c(0, 0, 20, 20), n = c(50, 50, 50, 50), threshold = 0.05
  )
  # Identical clusters: var_obs = 0 → deff = max(0, 1) = 1
  expect_equal(res_identical$deff, 1, tolerance = 1e-10)
  # Heterogeneous: strong between-cluster signal → deff >> 1
  expect_gt(res_heterog$deff, 1)
  expect_lt(res_heterog$n_eff, res_identical$n_eff)
})
