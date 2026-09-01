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
#   p_value = 1 - pnorm(1.376) = 0.0844 -> fail to reject at alpha=0.05
#
# Case 2  (x=15, n=100, threshold=0.05, alternative="greater"):
#   p_hat = 0.15, theta_app = 0.05
#   z = (0.15 - 0.05) / 0.02179 = 4.588
#   p_value = 1 - pnorm(4.588) ~= 2.24e-6 -> reject
#
# Case 3  (two.sided, x=8, n=100, threshold=0.08):
#   p_hat = 0.08 = threshold -> z = 0 -> p_value = 1 -> fail to reject
#
# Case 4  (imperfect test, se=0.9, sp=0.95, x=8, n=100, threshold=0.05):
#   theta_app = 0.05*0.9 + 0.95*0.05 = 0.045 + 0.0475 = 0.0925
#   p_hat = 0.08
#   se_null = sqrt(0.0925*0.9075/100) = sqrt(0.0000839) = 0.009158... hmm wait
#   se_null = sqrt(0.0925*(1-0.0925)/100) = sqrt(0.0925*0.9075/100)
#           = sqrt(0.083944/100) = sqrt(0.00083944) = 0.028973
#   z = (0.08 - 0.0925) / 0.028973 = -0.01250 / 0.028973 = -0.4314
#   p_value (greater) = 1 - pnorm(-0.4314) = 0.667 -> fail to reject
#
# Case 5  (clustered, x=rep(3,10), n=rep(10,10), icc=0.05, threshold=0.2):
#   p_hat = 30/100 = 0.3, n_bar=10, deff = 1+(10-1)*0.05 = 1.45
#   n_eff = 100/1.45 = 68.97, n_eff_adj = 68.97 (no FPC)
#   theta_app = 0.2 (perfect test)
#   se_null = sqrt(0.2*0.8/68.97) = sqrt(0.002320) = 0.04817
#   z = (0.3 - 0.2) / 0.04817 = 2.076
#   p_value (greater) = 1 - pnorm(2.076) = 0.0190 -> reject at alpha=0.05

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

test_that("Case 3: p_hat = threshold exactly -> z=0, p_value=0.5 (greater)", {
  res <- test_threshold(x = 8, n = 100, threshold = 0.08)

  expect_equal(res$statistic, 0, tolerance = 1e-6)
  expect_equal(res$p_value,   0.5, tolerance = 1e-6)
  expect_false(res$reject)
})

test_that("Case 3b: two.sided with p_hat = threshold -> z=0, p_value=1", {
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

  # Clustering inflates SE -> smaller |z| -> larger p-value
  expect_lt(abs(res_clustered$statistic), abs(res_srs$statistic))
  expect_gt(res_clustered$p_value, res_srs$p_value)
  expect_equal(res_clustered$deff, 1.45, tolerance = 1e-6)
})

test_that("alternative='less' reverses the test direction", {
  # x=0, n=100, threshold=0.05: p_hat=0, z strongly negative -> reject 'less'
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

  # FPC makes n_eff_adj larger -> narrower null SE -> larger |z|
  expect_gt(abs(res_fpc$statistic), abs(res_nofpc$statistic))
  expect_lt(res_fpc$p_value, res_nofpc$p_value)
})

test_that("return list has all expected fields in order", {
  res <- test_threshold(x = 8, n = 100, threshold = 0.05)
  expect_named(res, c("statistic", "p_value", "reject", "threshold", "alternative",
                       "prevalence", "ci_lower", "ci_upper", "ci_method",
                       "n_total", "n_eff", "n_eff_adj", "conf_level",
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

test_that("x=n (all positive): p_hat=1, z large positive -> reject 'greater'", {
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
  # p ~= 0.084 for Case 1 -- rejects at 0.10 but not at 0.05 or 0.01
  res_90 <- test_threshold(x = 8, n = 100, threshold = 0.05, conf_level = 0.90)
  res_99 <- test_threshold(x = 8, n = 100, threshold = 0.05, conf_level = 0.99)
  expect_true(res_90$reject)   # alpha=0.10: 0.084 < 0.10
  expect_false(res_99$reject)  # alpha=0.01: 0.084 > 0.01
})

test_that("n_bar=1 guard: all clusters of size 1 -> icc_used=0, deff=1", {
  res <- test_threshold(x = c(0, 1, 0, 1), n = c(1, 1, 1, 1), threshold = 0.3)
  expect_equal(res$icc_used, 0)
  expect_equal(res$deff,     1)
  expect_true(is.finite(res$statistic))
})

test_that("ICC estimated from heterogeneous clusters inflates deff and reduces n_eff", {
  # Identical clusters -> zero between-cluster variance -> deff clamped to 1
  # Highly heterogeneous clusters -> large ICC -> deff >> 1 -> much smaller n_eff
  res_identical <- test_threshold(
    x = c(5, 5, 5, 5), n = c(50, 50, 50, 50), threshold = 0.05
  )
  res_heterog <- test_threshold(
    x = c(0, 0, 20, 20), n = c(50, 50, 50, 50), threshold = 0.05
  )
  # Identical clusters: var_obs = 0 -> deff = max(0, 1) = 1
  expect_equal(res_identical$deff, 1, tolerance = 1e-10)
  # Heterogeneous: strong between-cluster signal -> deff >> 1
  expect_gt(res_heterog$deff, 1)
  expect_lt(res_heterog$n_eff, res_identical$n_eff)
})

# ---- stress tests: 15 edge cases to break test_threshold ----

test_that("TT-S1: n_total=1 (single obs): p_hat=0 or 1, z finite", {
  res0 <- test_threshold(x = 0, n = 1, threshold = 0.3)
  res1 <- test_threshold(x = 1, n = 1, threshold = 0.3)
  expect_true(is.finite(res0$statistic))
  expect_true(is.finite(res1$statistic))
  expect_equal(res0$prevalence, 0)
  expect_equal(res1$prevalence, 1)
})

test_that("TT-S2: threshold=0.999 (near boundary) produces valid statistic", {
  res <- test_threshold(x = 80, n = 100, threshold = 0.999)
  expect_true(is.finite(res$statistic))
  expect_true(is.finite(res$p_value))
  expect_gte(res$p_value, 0)
  expect_lte(res$p_value, 1)
})

test_that("TT-S3: threshold=0.001 (near boundary) produces valid statistic", {
  res <- test_threshold(x = 1, n = 100, threshold = 0.001)
  expect_true(is.finite(res$statistic))
  expect_gte(res$p_value, 0)
})

test_that("TT-S4: imperfect test with se+sp just above 1 (correction=0.01): warning emitted", {
  expect_warning(
    test_threshold(x = 30, n = 100, threshold = 0.05,
                   sensitivity = 0.51, specificity = 0.50),
    "numerically unstable"
  )
})

test_that("TT-S5: p_hat clamped to 0 via RG (p_true would be negative)", {
  # se=0.9, sp=0.9: threshold_app = 0.3*0.9 + 0.7*0.1 = 0.27+0.07 = 0.34
  # But with x=1, n=100: p_hat = 0.01; rg(0.01) = (0.01 - 0.1)/0.8 = -0.1125 -> clamped to 0
  res <- test_threshold(x = 1, n = 100, threshold = 0.30,
                        sensitivity = 0.9, specificity = 0.9)
  expect_equal(res$prevalence, 0)
  expect_true(is.finite(res$statistic))
})

test_that("TT-S6: p_hat clamped to 1 via RG (p_true would exceed 1)", {
  # se=0.9, sp=0.9: rg(0.99) = (0.99 - 0.1)/0.8 = 1.1125 -> clamped to 1
  res <- test_threshold(x = 99, n = 100, threshold = 0.50,
                        sensitivity = 0.9, specificity = 0.9)
  expect_equal(res$prevalence, 1)
  expect_true(is.finite(res$statistic))
})

test_that("TT-S7: icc=1 with multi-cluster data -> deff = n_bar (max possible)", {
  # ICC=1: every cluster is perfectly correlated; deff = 1 + (n_bar-1)*1 = n_bar
  res <- test_threshold(x = c(5, 5), n = c(10, 10), threshold = 0.3, icc = 1)
  expect_equal(res$deff, 10, tolerance = 1e-10)  # n_bar=10
  expect_equal(res$n_eff, 20 / 10, tolerance = 1e-10)  # n_total/deff = 2
})

test_that("TT-S8: large n (10000 obs) gives extremely small p-value when far from threshold", {
  res <- test_threshold(x = 500, n = 10000, threshold = 0.10)
  # p_hat=0.05, threshold=0.10: testing greater -> z strongly negative -> p near 1
  expect_gt(res$p_value, 0.99)
  expect_false(res$reject)
})

test_that("TT-S9: x=n=100, p_hat=1.0, threshold=0.999 -> barely reject 'greater'", {
  res <- test_threshold(x = 100, n = 100, threshold = 0.999)
  # p_hat=1; threshold_app=0.999 (perfect test): z=(1-0.999)/se_null, very small z
  expect_true(is.finite(res$statistic))
})

test_that("TT-S10: conf_level=0.5 (alpha=0.5) -> most tests reject", {
  res <- test_threshold(x = 8, n = 100, threshold = 0.05, conf_level = 0.5)
  # p_value=0.084, alpha=0.5 -> 0.084 < 0.5 -> reject
  expect_true(res$reject)
})

test_that("TT-S11: fpc_N = n_total + 1 (population barely larger than sample)", {
  # Near-census: fpc nearly 0 -> n_eff_adj very large -> se_null tiny -> |z| huge
  res <- suppressWarnings(
    test_threshold(x = 8, n = 100, threshold = 0.05, fpc_N = 101)
  )
  expect_true(is.finite(res$statistic))
  expect_true(is.finite(res$p_value))
})

test_that("TT-S12: all x=0 across many clusters -> p_hat=0 consistently", {
  res <- test_threshold(x = rep(0, 20), n = rep(50, 20), threshold = 0.05)
  expect_equal(res$prevalence, 0)
  expect_lt(res$statistic, 0)   # z negative (p_hat < threshold)
  expect_false(res$reject)      # testing 'greater'
})

test_that("TT-S13: x=n for all clusters -> p_hat=1 -> reject 'greater' for any threshold", {
  res <- test_threshold(x = rep(10, 5), n = rep(10, 5), threshold = 0.50)
  expect_equal(res$prevalence, 1)
  expect_gt(res$statistic, 0)
  expect_true(res$reject)
})

test_that("TT-S14: sensitivity=1, specificity=1 -> statistic equals naive formula", {
  # Perfect test: theta_app = threshold, p_hat_app = p_hat
  # z = (p_hat - threshold) / sqrt(threshold*(1-threshold)/n_eff_adj)
  res <- test_threshold(x = 15, n = 100, threshold = 0.10)
  z_expected <- (0.15 - 0.10) / sqrt(0.10 * 0.90 / 100)
  expect_equal(res$statistic, z_expected, tolerance = 1e-8)
})

test_that("TT-S15: unequal cluster sizes estimated ICC is non-negative and consistent", {
  res <- test_threshold(
    x = c(1, 5, 0, 12, 3),
    n = c(20, 40, 15, 60, 25),
    threshold = 0.10
  )
  expect_gte(res$icc_used, 0)
  expect_lte(res$icc_used, 1)
  expect_equal(res$deff, 1 + (mean(c(20,40,15,60,25)) - 1) * res$icc_used,
               tolerance = 1e-10)
  expect_true(is.finite(res$statistic))
})

# -- Second batch of stress tests: TT-S16 through TT-S30 ----------------------

test_that("TT-S16: fpc_N = n_total with clustering -> valid result (fpc > 0 because n_eff < n_total)", {
  # With clustering n_eff << n_total, so fpc_N = n_total does NOT collapse fpc to 0.
  # Old code wrongly warned here; fixed to check fpc_N <= n_eff instead.
  res <- test_threshold(x = c(5, 5), n = c(50, 50), threshold = 0.05,
                        icc = 0.5, fpc_N = 100)
  expect_true(is.finite(res$statistic))
  expect_true(is.finite(res$p_value))
})

test_that("TT-S17: all-singleton clusters (n=1 each) -> icc skipped, finite result", {
  # n_bar = 1 triggers the short-circuit in ICC estimation
  res <- test_threshold(x = c(1, 0, 1, 0, 1), n = c(1, 1, 1, 1, 1),
                        threshold = 0.30)
  expect_equal(res$icc_used, 0)
  expect_equal(res$deff, 1)
  expect_true(is.finite(res$statistic))
})

test_that("TT-S18: icc=0 forces deff=1 even with maximally heterogeneous clusters", {
  res_forced <- test_threshold(x = c(0, 20, 0, 20), n = rep(20, 4),
                               threshold = 0.30, icc = 0)
  expect_equal(res_forced$deff, 1)
  # Estimated ICC on the same data must give deff > 1
  res_est <- test_threshold(x = c(0, 20, 0, 20), n = rep(20, 4),
                            threshold = 0.30)
  expect_gt(res_est$deff, 1)
  # And a larger z magnitude under SRS (n_eff = n_total, not deflated)
  expect_gt(abs(res_forced$statistic), abs(res_est$statistic))
})

test_that("TT-S19: perfectly bimodal clusters -> estimated ICC clamped to 1, deff = n_bar", {
  # Half clusters all-zero, half all-positive: maximum within-cluster homogeneity
  res <- test_threshold(
    x = c(rep(0, 5), rep(10, 5)),
    n = rep(10, 10),
    threshold = 0.30
  )
  expect_equal(res$icc_used, 1, tolerance = 1e-10)  # clamped to 1
  expect_equal(res$deff, 10, tolerance = 1e-10)       # n_bar = 10
})

test_that("TT-S20: p_hat exactly at null -> z = 0, p_value = 0.5 for 'greater'", {
  # Perfect test: theta_app = threshold = 0.10; x/n = 0.10 exactly
  res <- test_threshold(x = 10, n = 100, threshold = 0.10)
  expect_equal(res$statistic, 0, tolerance = 1e-10)
  expect_equal(res$p_value, 0.5, tolerance = 1e-10)
  expect_false(res$reject)
})

test_that("TT-S21: conf_level=0.9999 -> CI nearly spans [0, 1] for middling prevalence", {
  res <- test_threshold(x = 8, n = 100, threshold = 0.05, conf_level = 0.9999)
  expect_lt(res$ci_lower, 0.02)
  expect_gt(res$ci_upper, 0.12)
  expect_true(res$ci_upper > res$ci_lower)
})

test_that("TT-S22: two-sided p_value = 2 x one-sided p_value when z > 0", {
  r_g  <- test_threshold(x = 12, n = 100, threshold = 0.05, alternative = "greater")
  r_2s <- test_threshold(x = 12, n = 100, threshold = 0.05, alternative = "two.sided")
  expect_gt(r_g$statistic, 0)  # verify z > 0
  expect_equal(r_2s$p_value, 2 * r_g$p_value, tolerance = 1e-12)
})

test_that("TT-S23: larger n with same rate -> strictly larger |z|", {
  r_small <- test_threshold(x = 8,  n = 100,  threshold = 0.05)
  r_large <- test_threshold(x = 80, n = 1000, threshold = 0.05)
  expect_equal(r_small$statistic / r_large$statistic, 1 / sqrt(10), tolerance = 1e-6)
})

test_that("TT-S24: prevalence is always inside its own Wald CI", {
  for (x_val in c(0, 1, 5, 50, 95, 99, 100)) {
    res <- test_threshold(x = x_val, n = 100, threshold = 0.30)
    expect_lte(res$ci_lower, res$prevalence + 1e-10)
    expect_gte(res$ci_upper, res$prevalence - 1e-10)
  }
})

test_that("TT-S25: supplied icc -> deff = 1 + (n_bar - 1) * icc exactly", {
  res <- test_threshold(x = c(3, 7, 5), n = c(30, 30, 30),
                        threshold = 0.15, icc = 0.3)
  expect_equal(res$deff, 1 + (30 - 1) * 0.3, tolerance = 1e-10)
  expect_equal(res$icc_used, 0.3)
})

test_that("TT-S26: single cluster with icc=NULL -> icc_used=0, deff=1", {
  res <- test_threshold(x = 15, n = 200, threshold = 0.05)
  expect_equal(res$icc_used, 0)
  expect_equal(res$deff, 1)
  expect_equal(res$n_eff, 200)
})

test_that("TT-S27: imperfect test RG-corrects prevalence below raw apparent rate", {
  # Raw rate = 20/100 = 0.20; se=0.80, sp=0.90 -> correction = 0.70
  # True prev = (0.20 - 0.10) / 0.70 ~= 0.1429
  res <- test_threshold(x = 20, n = 100, threshold = 0.10,
                        sensitivity = 0.80, specificity = 0.90)
  expect_lt(res$prevalence, 0.20)
  expect_equal(res$prevalence, (0.20 - 0.10) / 0.70, tolerance = 1e-6)
})

test_that("TT-S28: 'less' alternative with x=0 (p=0 far below threshold) -> reject", {
  # z = (0 - theta_app) / se_null << 0 -> one-sided p << 0.05
  res <- test_threshold(x = 0, n = 100, threshold = 0.20, alternative = "less")
  expect_lt(res$statistic, -4)
  expect_true(res$reject)
})

test_that("TT-S29: return list has all documented names regardless of clustering mode", {
  expected <- c("statistic", "p_value", "reject", "threshold", "alternative",
                "prevalence", "ci_lower", "ci_upper", "ci_method", "n_total",
                "n_eff", "n_eff_adj", "conf_level", "sensitivity", "specificity",
                "icc_used", "deff", "fpc_N")
  # SRS
  expect_named(test_threshold(x = 5, n = 100, threshold = 0.10), expected)
  # Supplied ICC
  expect_named(test_threshold(x = c(5, 10), n = c(50, 50), threshold = 0.10, icc = 0.05),
               expected)
  # FPC
  expect_named(test_threshold(x = 5, n = 100, threshold = 0.10, fpc_N = 1000), expected)
})

test_that("TT-S30: p_value in [0, 1] for all alternatives and extreme x values", {
  for (alt in c("greater", "less", "two.sided")) {
    for (x_val in c(0, 1, 50, 99, 100)) {
      res <- test_threshold(x = x_val, n = 100, threshold = 0.50, alternative = alt)
      expect_gte(res$p_value, 0)
      expect_lte(res$p_value, 1 + 1e-12)
      expect_true(is.logical(res$reject) && !is.na(res$reject))
    }
  }
})

test_that("TT-R7-1: non-numeric sensitivity/conf_level/icc give a friendly class error", {
  expect_error(test_threshold(x = 30, n = 100, threshold = 0.2, sensitivity = TRUE), "`sensitivity`")
  expect_error(test_threshold(x = 30, n = 100, threshold = 0.2, conf_level = "0.95"), "`conf_level`")
  expect_error(test_threshold(x = 30, n = 100, threshold = 0.2, icc = FALSE), "`icc`")
  expect_error(test_threshold(x = 30, n = 100, threshold = 0.2, sensitivity = list(0.9)), "`sensitivity`")
})

test_that("TT-R7-2: RG back-transform via .rogan_gladen matches the closed form", {
  res <- test_threshold(x = 30, n = 100, threshold = 0.2,
                        sensitivity = 0.9, specificity = 0.95)
  expect_equal(res$prevalence,
               (30 / 100 - (1 - 0.95)) / (0.9 + 0.95 - 1),
               tolerance = 1e-9)
})

test_that("TT-R8-1: ci_method offers wald, clopper-pearson, agresti-coull", {
  base <- test_threshold(x = 8, n = 100, threshold = 0.05)
  expect_equal(base$ci_method, "wald")

  for (m in c("wald", "clopper-pearson", "agresti-coull")) {
    res <- test_threshold(x = 8, n = 100, threshold = 0.05, ci_method = m)
    expect_equal(res$ci_method, m)
    expect_gte(res$ci_lower, 0)
    expect_lte(res$ci_upper, 1)
    expect_lte(res$ci_lower, res$ci_upper)
    # the hypothesis test is unaffected by ci_method
    expect_equal(res$statistic, base$statistic, tolerance = 1e-12)
    expect_equal(res$p_value,   base$p_value,   tolerance = 1e-12)
  }
})

test_that("TT-R8-2: invalid ci_method is rejected", {
  expect_error(test_threshold(x = 8, n = 100, threshold = 0.05, ci_method = "wilson"),
               "ci_method")
  expect_error(test_threshold(x = 8, n = 100, threshold = 0.05,
                              ci_method = c("wald", "wald")),
               "ci_method")
})

test_that("TT-R8-3: a supplied icc has no effect with a single cluster", {
  a <- test_threshold(x = 8, n = 50, threshold = 0.05)
  b <- test_threshold(x = 8, n = 50, threshold = 0.05, icc = 0.05)
  expect_equal(b$deff, 1)
  expect_equal(b$icc_used, 0)
  expect_equal(a$statistic, b$statistic, tolerance = 1e-12)
})

test_that("TT-R8-4: n_eff_adj equals n_eff without FPC, exceeds it with FPC", {
  no_fpc <- test_threshold(x = 30, n = 100, threshold = 0.2)
  expect_equal(no_fpc$n_eff_adj, no_fpc$n_eff)

  with_fpc <- test_threshold(x = 30, n = 100, threshold = 0.2, fpc_N = 400)
  expect_gt(with_fpc$n_eff_adj, with_fpc$n_eff)
})
