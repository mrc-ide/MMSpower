# Hand-checked reference cases for estimate_prevalence()
#
# Formula: Wald CI on apparent prevalence, then Rogan-Gladen to true scale.
#
#   p_hat   = sum(x) / sum(n)
#   n_eff   = n_total / deff
#   se      = sqrt(p_hat*(1-p_hat) / n_eff)
#   moe_app = z * se          where z = qnorm(0.975) = 1.959964
#   CI_app  = [p_hat - moe_app, p_hat + moe_app]
#   p_true  = (p_hat - (1-sp)) / (se_test + sp - 1)    [Rogan-Gladen]
#   CI_true = Rogan-Gladen applied to CI_app endpoints
#
# Case 1  (single site, perfect test, x=30, n=100):
#   p_hat = 0.3, deff = 1, n_eff = 100
#   se = sqrt(0.3*0.7/100) = 0.04583
#   moe = 1.959964 * 0.04583 = 0.08982
#   CI = [0.2102, 0.3898]
#
# Case 2  (supplied icc=0.05, 10 uniform sites of 10 each):
#   deff = 1 + 9*0.05 = 1.45, n_eff = 100/1.45 = 68.97
#   se = sqrt(0.3*0.7/68.97) = 0.05519
#   moe = 1.959964 * 0.05519 = 0.10817
#   CI = [0.1918, 0.4082]
#
# Case 3  (imperfect test, se=0.9, sp=0.95, x=30, n=100):
#   p_apparent = 0.3, correction = 0.85
#   p_true = (0.3 - 0.05) / 0.85 = 0.2941
#   moe_apparent = 0.08982  (same as Case 1)
#   moe_true = 0.08982 / 0.85 = 0.10567
#   CI_true = [(0.2102-0.05)/0.85, (0.3898-0.05)/0.85]
#           = [0.1885, 0.3998]

test_that("Case 1: single site, perfect test -- Wald CI", {
  res <- estimate_prevalence(x = 30, n = 100)

  expect_equal(res$prevalence,     0.3,    tolerance = 1e-6)
  expect_equal(res$deff,           1,      tolerance = 1e-6)
  expect_equal(res$icc_used,       0,      tolerance = 1e-6)
  expect_equal(res$n_total,        100)
  expect_equal(res$n_eff,          100,    tolerance = 1e-6)

  # Wald 95% CI (hand-checked)
  expect_equal(res$ci_lower,       0.2102, tolerance = 5e-4)
  expect_equal(res$ci_upper,       0.3898, tolerance = 5e-4)
  expect_equal(res$margin_of_error, 0.08982, tolerance = 5e-4)
})

test_that("Case 2: supplied ICC widens CI via design effect", {
  res <- estimate_prevalence(
    x   = rep(3, 10),
    n   = rep(10, 10),
    icc = 0.05
  )

  expect_equal(res$prevalence, 0.3,  tolerance = 1e-6)
  expect_equal(res$deff,       1.45, tolerance = 1e-6)
  expect_equal(res$n_eff,      100 / 1.45, tolerance = 1e-4)

  # CI should be wider than the no-clustering case
  res_plain <- estimate_prevalence(x = 30, n = 100, icc = 0)
  expect_gt(res$ci_upper - res$ci_lower,
            res_plain$ci_upper - res_plain$ci_lower)

  expect_equal(res$ci_lower, 0.1918, tolerance = 5e-4)
  expect_equal(res$ci_upper, 0.4082, tolerance = 5e-4)
})

test_that("Case 3: imperfect test applies Rogan-Gladen to estimate and CI", {
  res <- estimate_prevalence(x = 30, n = 100,
                             sensitivity = 0.9, specificity = 0.95)

  # p_true = (0.3 - 0.05) / 0.85
  expect_equal(res$prevalence, 0.25 / 0.85, tolerance = 1e-6)

  # MOE inflated by 1/correction
  res_perfect <- estimate_prevalence(x = 30, n = 100)
  expect_equal(res$margin_of_error,
               res_perfect$margin_of_error / 0.85,
               tolerance = 1e-6)

  expect_equal(res$ci_lower, 0.1885, tolerance = 5e-4)
  expect_equal(res$ci_upper, 0.3998, tolerance = 5e-4)
})

test_that("ICC estimated from data is non-negative and finite", {
  res <- estimate_prevalence(
    x = c(0, 4, 0, 22, 25, 16, 12, 8),
    n = c(60, 80, 70, 100, 40, 60, 50, 90)
  )
  expect_true(is.finite(res$icc_used))
  expect_gte(res$icc_used, 0)
  expect_lte(res$icc_used, 1)
  expect_gte(res$deff, 1)
  expect_equal(res$n_total, 550)
})

test_that("icc = 0 forces no clustering adjustment", {
  res_null <- estimate_prevalence(x = 30, n = 100, icc = 0)
  res_no   <- estimate_prevalence(x = 30, n = 100, icc = NULL)  # single site

  expect_equal(res_null$deff, 1, tolerance = 1e-6)
  # Both single-site paths should agree
  expect_equal(res_null$ci_lower, res_no$ci_lower, tolerance = 1e-10)
})

test_that("perfect test is identity: se=sp=1 changes nothing", {
  res_default <- estimate_prevalence(x = 30, n = 100)
  res_explicit <- estimate_prevalence(x = 30, n = 100,
                                      sensitivity = 1, specificity = 1)
  expect_equal(res_default$prevalence, res_explicit$prevalence)
  expect_equal(res_default$ci_lower,   res_explicit$ci_lower)
  expect_equal(res_default$ci_upper,   res_explicit$ci_upper)
})

test_that("input validation catches bad arguments", {
  expect_error(estimate_prevalence(x = c(3, 5), n = c(10)),
               "length")
  expect_error(estimate_prevalence(x = 60, n = 50),
               "cannot exceed")
  expect_error(estimate_prevalence(x = 10, n = 50,
                                   sensitivity = 0.3, specificity = 0.3),
               "must exceed 1")
  # negative counts
  expect_error(estimate_prevalence(x = -1, n = 100),
               "non-negative")
  # conf_level out of range
  expect_error(estimate_prevalence(x = 30, n = 100, conf_level = 0), "`conf_level`")
  expect_error(estimate_prevalence(x = 30, n = 100, conf_level = 1), "`conf_level`")
  # icc out of range when supplied
  expect_error(estimate_prevalence(x = 30, n = 100, icc =  2),  "`icc`")
  expect_error(estimate_prevalence(x = 30, n = 100, icc = -0.1), "`icc`")
  # fpc_N must be positive
  expect_error(estimate_prevalence(x = 30, n = 100, fpc_N =   0), "`fpc_N`")
  expect_error(estimate_prevalence(x = 30, n = 100, fpc_N = -10), "`fpc_N`")
  # fpc_N smaller than sample
  expect_error(estimate_prevalence(x = c(3,3), n = c(10,10), fpc_N = 5),
               "at least as large")
})

test_that("n_bar=1 (clusters of 1) does not produce Inf/NaN ICC", {
  # Kish formula has (n_bar-1) in denominator; n_bar=1 means division by zero.
  # Should fall back to icc=0, deff=1, not crash or return NaN.
  res <- estimate_prevalence(x = c(0, 1, 0, 1, 1), n = c(1, 1, 1, 1, 1))
  expect_equal(res$icc_used, 0)
  expect_equal(res$deff,     1)
  expect_true(is.finite(res$prevalence))
  expect_true(is.finite(res$margin_of_error))
})

test_that("deff and icc_used are mutually consistent after clamping", {
  # Extreme clustering: all variation between sites, none within.
  # Raw deff could imply icc > 1; after clamping, deff must be recomputed.
  res <- estimate_prevalence(
    x = c(0, 0, 10, 10),
    n = c(10, 10, 10, 10)
  )
  expect_equal(res$deff, 1 + (10 - 1) * res$icc_used, tolerance = 1e-10)
  expect_lte(res$icc_used, 1)
  expect_gte(res$icc_used, 0)
})

test_that("round-trip: estimate_prevalence recovers design_precision MOE (SRS)", {
  # design_precision(0.3, 0.05) -> n = 323
  # observe exactly 30% -> x = round(0.3 * 323) = 97
  res <- estimate_prevalence(x = 97, n = 323)
  expect_equal(res$margin_of_error, 0.05, tolerance = 5e-3)
})

test_that("round-trip: estimate_prevalence recovers design_precision MOE (clustered)", {
  # design_precision(0.3, 0.05, n_per_site=10, icc=0.05) -> 47 sites of 10
  # observe 30% per site -> x = 3 per site
  res <- estimate_prevalence(x = rep(3, 47), n = rep(10, 47), icc = 0.05)
  expect_equal(res$margin_of_error, 0.05, tolerance = 5e-3)
})

# ---- Round 4: 15 new edge cases ----

test_that("EP-1: all-zero prevalence collapses CI to [0, 0]", {
  # p=0 → se=0 → Wald CI=[0,0], moe=0. All outputs finite.
  res <- estimate_prevalence(x = c(0, 0, 0), n = c(10, 10, 10))
  expect_equal(res$prevalence,      0)
  expect_equal(res$ci_lower,        0)
  expect_equal(res$ci_upper,        0)
  expect_equal(res$margin_of_error, 0)
  expect_true(is.finite(res$deff))
})

test_that("EP-2: all-100% prevalence collapses CI to [1, 1]", {
  res <- estimate_prevalence(x = c(10, 10, 10), n = c(10, 10, 10))
  expect_equal(res$prevalence,      1)
  expect_equal(res$ci_lower,        1)
  expect_equal(res$ci_upper,        1)
  expect_equal(res$margin_of_error, 0)
})

test_that("EP-3: single site x=0 returns finite output, moe=0", {
  res <- estimate_prevalence(x = 0, n = 100)
  expect_equal(res$prevalence, 0)
  expect_equal(res$margin_of_error, 0)
  expect_true(is.finite(res$n_eff))
})

test_that("EP-4: conf_level=0.99 produces wider CI than 0.95", {
  r99 <- estimate_prevalence(x = 30, n = 100, conf_level = 0.99)
  r95 <- estimate_prevalence(x = 30, n = 100)
  width99 <- r99$ci_upper - r99$ci_lower
  width95 <- r95$ci_upper - r95$ci_lower
  expect_gt(width99, width95)
  # z_0.99=2.576 vs z_0.95=1.960: ratio should be ~1.315
  expect_equal(width99 / width95, 2.576 / 1.960, tolerance = 0.01)
})

test_that("EP-5: supplied icc=1 (maximum) → deff=n_bar, n_eff=n_clusters", {
  # deff = 1 + (n_bar - 1)*1 = n_bar = 10
  # n_eff = 100 / 10 = 10 = n_clusters
  res <- estimate_prevalence(x = rep(3, 10), n = rep(10, 10), icc = 1)
  expect_equal(res$deff,  10,  tolerance = 1e-10)
  expect_equal(res$n_eff, 10,  tolerance = 1e-10)
})

test_that("EP-6: identical cluster prevalences → estimated ICC=0, deff=1", {
  # No between-cluster variation: var_obs=0, so deff=max(0,1)=1, icc=0
  res <- estimate_prevalence(x = rep(3, 8), n = rep(10, 8))
  expect_equal(res$icc_used, 0, tolerance = 1e-10)
  expect_equal(res$deff,     1, tolerance = 1e-10)
})

test_that("EP-7: character x gives informative type error (not NA/NaN/Inf message)", {
  expect_error(
    estimate_prevalence(x = "30", n = 100),
    "numeric vectors"
  )
})

test_that("EP-8: fpc_N = n_total (census) gives warning and moe=0", {
  expect_warning(
    res <- estimate_prevalence(x = 30, n = 100, fpc_N = 100),
    "entire population"
  )
  expect_equal(res$margin_of_error, 0)
})

test_that("EP-9: sensitivity=0.5, specificity=1 doubles the prevalence estimate", {
  # correction = 0.5; p_true = p_apparent / 0.5 = 2 * p_apparent
  # x=15, n=100 → p_apparent=0.15 → p_true=0.30
  res <- estimate_prevalence(x = 15, n = 100, sensitivity = 0.5, specificity = 1)
  expect_equal(res$prevalence, 0.30, tolerance = 1e-6)
})

test_that("EP-10: integer inputs 100L, 1000L work identically to double", {
  r_int <- estimate_prevalence(x = 100L, n = 1000L)
  r_dbl <- estimate_prevalence(x = 100,  n = 1000)
  expect_equal(r_int$prevalence,      r_dbl$prevalence)
  expect_equal(r_int$margin_of_error, r_dbl$margin_of_error)
})

test_that("EP-11: maximum between-cluster heterogeneity → icc=1, deff=n_bar", {
  # Two clusters: one all-negative, one all-positive → ICC clamped to 1
  res <- estimate_prevalence(x = c(0, 10), n = c(10, 10))
  expect_equal(res$icc_used, 1,  tolerance = 1e-10)
  expect_equal(res$deff,     10, tolerance = 1e-10)
  expect_equal(res$prevalence, 0.5, tolerance = 1e-6)
})

test_that("EP-12: fpc_N=Inf is rejected", {
  expect_error(
    estimate_prevalence(x = 30, n = 100, fpc_N = Inf),
    "`fpc_N`"
  )
})

test_that("EP-13: very large fpc_N has negligible effect on moe", {
  r_fpc <- estimate_prevalence(x = 30, n = 100, fpc_N = 1e8)
  r_no  <- estimate_prevalence(x = 30, n = 100)
  # FPC factor ≈ sqrt((1e8 - 100) / (1e8 - 1)) ≈ 1 − 5e-7
  expect_equal(r_fpc$margin_of_error, r_no$margin_of_error, tolerance = 1e-5)
})

test_that("EP-14: all-zero multi-cluster: var_obs=0, estimated icc=0, deff=1", {
  res <- estimate_prevalence(x = c(0, 0), n = c(10, 10))
  expect_equal(res$prevalence, 0)
  expect_equal(res$icc_used,   0, tolerance = 1e-10)
  expect_equal(res$deff,       1, tolerance = 1e-10)
})

test_that("EP-15: single observation (x=1, n=1): n_bar=1 guard, icc=0, deff=1", {
  res <- estimate_prevalence(x = 1, n = 1)
  expect_equal(res$icc_used, 0)
  expect_equal(res$deff,     1)
  expect_equal(res$prevalence, 1)
  expect_true(is.finite(res$margin_of_error))
})
