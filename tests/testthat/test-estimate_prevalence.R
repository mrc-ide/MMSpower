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
  expect_equal(res$moe, 0.08982, tolerance = 5e-4)
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
  expect_equal(res$moe,
               res_perfect$moe / 0.85,
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
  expect_true(is.finite(res$moe))
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
  expect_equal(res$moe, 0.05, tolerance = 5e-3)
})

test_that("round-trip: estimate_prevalence recovers design_precision MOE (clustered, icc supplied)", {
  # design_precision(0.3, 0.05, n_per_site=10, icc=0.05) -> 47 sites of 10
  # observe 30% per site -> x = 3 per site
  # icc MUST be supplied explicitly: identical clusters give estimated icc=0 (no between-cluster
  # variance to estimate from), which would return moe ≈ 0.041, not 0.05.
  res <- estimate_prevalence(x = rep(3, 47), n = rep(10, 47), icc = 0.05)
  expect_equal(res$moe, 0.05, tolerance = 5e-3)
})

test_that("round-trip caveat: without icc supplied, identical clusters estimate icc=0 and round-trip fails", {
  # This test documents the limitation: the clustered round-trip only holds when icc is supplied.
  # Identical cluster prevalences → var_obs = 0 → estimated icc = 0 → deff = 1 → moe ≈ SRS moe.
  # design_precision assumed icc=0.05 and returned n=468 (47 sites x 10); EP without icc gives
  # moe ≈ 0.041 (SRS-like), not the 0.05 that design_precision targeted.
  res_no_icc  <- estimate_prevalence(x = rep(3, 47), n = rep(10, 47))
  res_with_icc <- estimate_prevalence(x = rep(3, 47), n = rep(10, 47), icc = 0.05)

  expect_equal(res_no_icc$icc_used, 0, tolerance = 1e-10)   # estimated from identical clusters
  expect_equal(res_no_icc$deff,     1, tolerance = 1e-10)
  expect_lt(res_no_icc$moe, 0.05)                # narrower than intended
  expect_equal(res_with_icc$moe, 0.05, tolerance = 5e-3)  # correct only with icc
})

# ---- Round 4: 15 new edge cases ----

test_that("EP-1: all-zero prevalence collapses CI to [0, 0]", {
  # p=0 → se=0 → Wald CI=[0,0], moe=0. All outputs finite.
  res <- estimate_prevalence(x = c(0, 0, 0), n = c(10, 10, 10))
  expect_equal(res$prevalence,      0)
  expect_equal(res$ci_lower,        0)
  expect_equal(res$ci_upper,        0)
  expect_equal(res$moe, 0)
  expect_true(is.finite(res$deff))
})

test_that("EP-2: all-100% prevalence collapses CI to [1, 1]", {
  res <- estimate_prevalence(x = c(10, 10, 10), n = c(10, 10, 10))
  expect_equal(res$prevalence,      1)
  expect_equal(res$ci_lower,        1)
  expect_equal(res$ci_upper,        1)
  expect_equal(res$moe, 0)
})

test_that("EP-3: single site x=0 returns finite output, moe=0", {
  res <- estimate_prevalence(x = 0, n = 100)
  expect_equal(res$prevalence, 0)
  expect_equal(res$moe, 0)
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

test_that("EP-8: fpc_N = n_total (SRS census) → error: n_eff = fpc_N, variance undefined", {
  # SRS: n_eff = n_total = fpc_N → fpc = 0 → CI undefined. Clustered case is valid.
  expect_error(
    estimate_prevalence(x = 30, n = 100, fpc_N = 100),
    "n_eff"
  )
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
  expect_equal(r_int$moe, r_dbl$moe)
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
  expect_equal(r_fpc$moe, r_no$moe, tolerance = 1e-5)
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
  expect_true(is.finite(res$moe))
})

# ---- Round 5: 8 new edge cases ----

test_that("EP-R5-1: n as character gives informative type error", {
  expect_error(
    estimate_prevalence(x = 30, n = "100"),
    "numeric vectors"
  )
})

test_that("EP-R5-2: length mismatch x=length 3, n=length 2 → informative error", {
  expect_error(
    estimate_prevalence(x = c(1, 2, 3), n = c(10, 10)),
    "same length"
  )
})

test_that("EP-R5-3: single large cluster (icc=NULL) gives correct narrow CI", {
  # n=1000: moe = 1.96 * sqrt(0.3*0.7/1000) = 0.02837
  res <- estimate_prevalence(x = 300, n = 1000)
  expect_equal(res$prevalence, 0.3, tolerance = 1e-6)
  expect_equal(res$deff, 1, tolerance = 1e-6)
  expect_equal(res$moe, 0.02840, tolerance = 5e-4)
})

test_that("EP-R5-4: fpc_N = n_total+1 (near-census) dramatically reduces moe", {
  # sampling fraction = 100/101 = 99.0%; FPC = sqrt(1/100) = 0.1
  r_fpc <- estimate_prevalence(x = 30, n = 100, fpc_N = 101)
  r_no  <- estimate_prevalence(x = 30, n = 100)
  expect_lt(r_fpc$moe, r_no$moe / 5)
  expect_equal(r_fpc$moe, r_no$moe * sqrt(1/100),
               tolerance = 1e-5)
})

test_that("EP-R5-5: specificity=0 is rejected (boundary, not in (0,1])", {
  expect_error(estimate_prevalence(x = 30, n = 100, specificity = 0), "`specificity`")
})

test_that("EP-R5-6: conf_level as vector is rejected with length error", {
  expect_error(
    estimate_prevalence(x = 30, n = 100, conf_level = c(0.9, 0.95)),
    "single number"
  )
})

test_that("EP-R5-7: unequal cluster sizes return finite, correct n_total", {
  res <- estimate_prevalence(x = c(1, 2, 3, 4), n = c(5, 10, 15, 20))
  expect_equal(res$n_total, 50)
  expect_equal(res$prevalence, 10/50, tolerance = 1e-6)
  expect_true(is.finite(res$moe))
  expect_true(is.finite(res$icc_used))
})

test_that("EP-R5-8: RG-corrected p > 1 is clamped to 1 without error", {
  # p_apparent=1, se=0.8, sp=0.9 → p_true = (1-0.1)/0.7 = 1.286 → clamped to 1
  res <- estimate_prevalence(x = 1, n = 1, sensitivity = 0.8, specificity = 0.9)
  expect_equal(res$prevalence, 1)
  expect_true(is.finite(res$moe))
})

# ---- Round 6: 15 + 1 new edge cases ----

test_that("EP-R6-1: x=NA (logical NA) gives informative class error mentioning NA", {
  expect_error(
    estimate_prevalence(x = NA, n = 100),
    "logical"
  )
  expect_error(
    estimate_prevalence(x = NA, n = 100),
    "NA"
  )
})

test_that("EP-R6-2: n vector with NA gives is.finite error at position", {
  expect_error(
    estimate_prevalence(x = c(1, 2, 3), n = c(10, NA, 10)),
    "NA, NaN, or infinite"
  )
})

test_that("EP-R6-3: x negative at position 2 gives position-aware error", {
  expect_error(
    estimate_prevalence(x = c(1, -1, 2), n = c(10, 10, 10)),
    "non-negative"
  )
})

test_that("EP-R6-4: n=0 at position 2 gives position-aware error", {
  expect_error(
    estimate_prevalence(x = c(1, 0, 2), n = c(10, 0, 10)),
    "positive for every cluster"
  )
})

test_that("EP-R6-5: x=99, n=100 → CI upper clamped to 1", {
  res <- estimate_prevalence(x = 99, n = 100)
  expect_equal(res$prevalence, 0.99, tolerance = 1e-6)
  expect_equal(res$ci_upper,   1.0,  tolerance = 1e-10)
})

test_that("EP-R6-6: clustering + FPC together: deff and FPC both applied", {
  r_both <- estimate_prevalence(x = rep(3, 10), n = rep(10, 10),
                                icc = 0.05, fpc_N = 500)
  r_noFPC <- estimate_prevalence(x = rep(3, 10), n = rep(10, 10), icc = 0.05)
  # FPC reduces moe when sampling fraction is non-negligible
  expect_lt(r_both$moe, r_noFPC$moe)
  expect_equal(r_both$deff, 1.45, tolerance = 1e-6)
})

test_that("EP-R6-7: fractional x at position 1 gives position-aware error", {
  expect_error(
    estimate_prevalence(x = c(1.5, 2, 3), n = c(10, 10, 10)),
    "whole numbers"
  )
})

test_that("EP-R6-8: exactly 2 clusters allows ICC estimation without error", {
  res <- estimate_prevalence(x = c(2, 4), n = c(10, 10))
  expect_true(is.finite(res$icc_used))
  expect_true(is.finite(res$moe))
})

test_that("EP-R6-9: sensitivity=1.1 is rejected", {
  expect_error(estimate_prevalence(x = 30, n = 100, sensitivity = 1.1), "`sensitivity`")
})

test_that("EP-R6-10: very low prevalence (1 positive / 50) returns valid output", {
  res <- estimate_prevalence(x = c(0, 0, 0, 0, 1), n = rep(10, 5))
  expect_equal(res$prevalence, 0.02, tolerance = 1e-6)
  expect_true(is.finite(res$moe))
  expect_gte(res$ci_lower, 0)
})

test_that("EP-R6-11: 100 clusters returns finite icc, deff, moe", {
  set.seed(42)
  x100 <- rbinom(100, 20, 0.2)
  res  <- estimate_prevalence(x = x100, n = rep(20, 100))
  expect_equal(res$n_total, 2000)
  expect_true(is.finite(res$icc_used))
  expect_true(is.finite(res$moe))
  expect_gte(res$icc_used, 0)
  expect_lte(res$icc_used, 1)
})

test_that("EP-R6-12: fractional n at position 2 gives position-aware error", {
  expect_error(
    estimate_prevalence(x = c(1, 2, 3), n = c(10, 10.5, 10)),
    "whole numbers"
  )
})

test_that("EP-R6-13: non-integer fpc_N (e.g. 100.5) is accepted", {
  # No integer requirement on fpc_N; population size can be approximated
  expect_no_error(
    suppressWarnings(estimate_prevalence(x = 30, n = 100, fpc_N = 100.5))
  )
})

test_that("EP-R6-14: icc=NULL with single cluster falls back to icc=0, deff=1", {
  res <- estimate_prevalence(x = 30, n = 100, icc = NULL)
  expect_equal(res$icc_used, 0, tolerance = 1e-10)
  expect_equal(res$deff,     1, tolerance = 1e-10)
})

test_that("EP-R6-15: very unequal clusters (1 vs 1000) return finite output", {
  res <- estimate_prevalence(x = c(0, 300), n = c(1, 1000))
  expect_equal(res$n_total, 1001)
  expect_equal(res$prevalence, 300/1001, tolerance = 1e-6)
  expect_true(is.finite(res$moe))
})

test_that("EP-R6-bonus: multi-cluster n_bar=1 (n=c(1,1,1,1)) uses guard, not ICC formula", {
  # n_bar = mean(c(1,1,1,1)) = 1 → Kish denominator (n_bar-1) = 0 → div/0 without guard
  # Guard: n_bar==1 → icc_used=0, deff=1 (same path as single-cluster guard)
  res <- estimate_prevalence(x = c(0, 1, 0, 1), n = c(1, 1, 1, 1))
  expect_equal(res$icc_used, 0)
  expect_equal(res$deff,     1)
  expect_true(is.finite(res$prevalence))
  expect_true(is.finite(res$moe))
})

# ---- CI methods: return list structure, validation, and method behaviour ----

test_that("return list contains all expected fields (wald default)", {
  res <- estimate_prevalence(x = 30, n = 100)
  expect_named(res, c("prevalence", "ci_lower", "ci_upper", "moe",
                       "moe_lower", "moe_upper", "method",
                       "n_total", "n_eff", "conf_level",
                       "sensitivity", "specificity", "icc_used", "deff", "fpc_N"),
               ignore.order = FALSE)
})

test_that("wald: moe_lower = moe_upper = moe (symmetric interval)", {
  res <- estimate_prevalence(x = 30, n = 100)
  expect_equal(res$method, "wald")
  expect_equal(res$moe_lower, res$moe, tolerance = 1e-10)
  expect_equal(res$moe_upper, res$moe, tolerance = 1e-10)
})

test_that("clopper-pearson: method field set, interval valid and finite", {
  res <- estimate_prevalence(x = 30, n = 100, method = "clopper-pearson")
  expect_equal(res$method, "clopper-pearson")
  expect_true(is.finite(res$ci_lower))
  expect_true(is.finite(res$ci_upper))
  expect_gte(res$ci_lower, 0)
  expect_lte(res$ci_upper, 1)
  expect_lte(res$ci_lower, res$prevalence)
  expect_gte(res$ci_upper, res$prevalence)
  # moe is average half-width
  expect_equal(res$moe, (res$ci_upper - res$ci_lower) / 2, tolerance = 1e-10)
})

test_that("clopper-pearson: moe_lower and moe_upper are one-sided distances", {
  res <- estimate_prevalence(x = 30, n = 100, method = "clopper-pearson")
  expect_equal(res$moe_lower, res$prevalence - res$ci_lower, tolerance = 1e-10)
  expect_equal(res$moe_upper, res$ci_upper  - res$prevalence, tolerance = 1e-10)
  # CP is asymmetric for interior p_hat
  expect_false(isTRUE(all.equal(res$moe_lower, res$moe_upper, tolerance = 1e-6)))
})

test_that("agresti-coull: method field set, interval valid and finite", {
  res <- estimate_prevalence(x = 30, n = 100, method = "agresti-coull")
  expect_equal(res$method, "agresti-coull")
  expect_true(is.finite(res$ci_lower))
  expect_true(is.finite(res$ci_upper))
  expect_gte(res$ci_lower, 0)
  expect_lte(res$ci_upper, 1)
  expect_equal(res$moe, (res$ci_upper - res$ci_lower) / 2, tolerance = 1e-10)
})

test_that("agresti-coull: moe_lower and moe_upper are one-sided distances", {
  res <- estimate_prevalence(x = 30, n = 100, method = "agresti-coull")
  expect_equal(res$moe_lower, res$prevalence - res$ci_lower, tolerance = 1e-10)
  expect_equal(res$moe_upper, res$ci_upper  - res$prevalence, tolerance = 1e-10)
  # AC centres on p_tilde != p_hat, so moe_lower != moe_upper
  expect_false(isTRUE(all.equal(res$moe_lower, res$moe_upper, tolerance = 1e-6)))
})

test_that("all three methods agree on point estimate (prevalence unchanged)", {
  # CI method changes the interval, not the point estimate
  wald <- estimate_prevalence(x = 15, n = 80)
  cp   <- estimate_prevalence(x = 15, n = 80, method = "clopper-pearson")
  ac   <- estimate_prevalence(x = 15, n = 80, method = "agresti-coull")
  expect_equal(wald$prevalence, cp$prevalence,  tolerance = 1e-10)
  expect_equal(wald$prevalence, ac$prevalence,  tolerance = 1e-10)
})

test_that("clopper-pearson: x=0 gives ci_lower=0 without error", {
  res <- estimate_prevalence(x = 0, n = 50, method = "clopper-pearson")
  expect_equal(res$ci_lower, 0)
  expect_gt(res$ci_upper,    0)
  expect_true(is.finite(res$ci_upper))
})

test_that("clopper-pearson: x=n gives ci_upper=1 without error", {
  res <- estimate_prevalence(x = 50, n = 50, method = "clopper-pearson")
  expect_equal(res$ci_upper, 1)
  expect_lt(res$ci_lower,    1)
  expect_true(is.finite(res$ci_lower))
})

test_that("agresti-coull: x=0 gives ci_lower >= 0 without error", {
  res <- estimate_prevalence(x = 0, n = 50, method = "agresti-coull")
  expect_gte(res$ci_lower, 0)
  expect_true(is.finite(res$ci_upper))
})

test_that("clopper-pearson with imperfect test applies Rogan-Gladen to CI endpoints", {
  res_perfect  <- estimate_prevalence(x = 30, n = 100, method = "clopper-pearson")
  res_imperfect <- estimate_prevalence(x = 30, n = 100,
                                        sensitivity = 0.9, specificity = 0.95,
                                        method = "clopper-pearson")
  # Imperfect test: correction = 0.85 < 1 → wider CI (moe inflated by 1/correction)
  expect_gt(res_imperfect$moe, res_perfect$moe)
  # Point estimate shifted by RG
  expect_false(isTRUE(all.equal(res_imperfect$prevalence, res_perfect$prevalence)))
})

test_that("clopper-pearson with clustering widens CI via deff", {
  res_srs      <- estimate_prevalence(x = rep(3, 10), n = rep(10, 10),
                                       icc = 0,    method = "clopper-pearson")
  res_clustered <- estimate_prevalence(x = rep(3, 10), n = rep(10, 10),
                                        icc = 0.05, method = "clopper-pearson")
  # Clustering inflates n_eff (deff > 1) → wider interval
  expect_gt(res_clustered$ci_upper - res_clustered$ci_lower,
            res_srs$ci_upper       - res_srs$ci_lower)
})

test_that("invalid method gives informative error", {
  expect_error(
    estimate_prevalence(x = 30, n = 100, method = "exact"),
    "'wald', 'clopper-pearson', or 'agresti-coull'"
  )
})

test_that("method as vector is rejected with informative error", {
  expect_error(
    estimate_prevalence(x = 30, n = 100, method = c("wald", "clopper-pearson")),
    "single character string"
  )
})

test_that("asymmetric method emits a message when moe_lower != moe_upper", {
  expect_message(
    estimate_prevalence(x = 30, n = 100, method = "clopper-pearson"),
    "asymmetric"
  )
  expect_message(
    estimate_prevalence(x = 30, n = 100, method = "agresti-coull"),
    "asymmetric"
  )
})

test_that("wald emits no message about asymmetry", {
  expect_no_message(
    estimate_prevalence(x = 30, n = 100, method = "wald")
  )
})

test_that("EP-R7-1: logical sensitivity/specificity are rejected, not coerced", {
  expect_error(estimate_prevalence(x = 30, n = 100, sensitivity = TRUE), "`sensitivity`")
  expect_error(estimate_prevalence(x = 30, n = 100, specificity = FALSE), "`specificity`")
})

test_that("EP-R7-2: logical icc is rejected with a class error, not returned as icc_used", {
  expect_error(
    estimate_prevalence(x = c(3, 3), n = c(10, 10), icc = FALSE),
    "`icc`"
  )
})

test_that("EP-R7-3: list-valued scalar params give a friendly class error, not a raw R error", {
  expect_error(estimate_prevalence(x = 30, n = 100, sensitivity = list(0.9)), "`sensitivity`")
  expect_error(estimate_prevalence(x = 30, n = 100, conf_level  = list(0.95)), "`conf_level`")
})

test_that("EP-R7-4: near-symmetric clopper-pearson interval emits no asymmetry message", {
  # p_hat = 0.5, large n: the two half-widths differ by < 10% of moe
  expect_no_message(estimate_prevalence(x = 500, n = 1000, method = "clopper-pearson"))
})
