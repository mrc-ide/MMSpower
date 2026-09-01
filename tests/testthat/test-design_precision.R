# Hand-checked reference cases for design_precision()
#
# Formula: n = z^2 * p_app * (1 - p_app) / (moe^2 * correction^2)
#   where p_app = p*se + (1-p)*(1-sp)   [apparent prevalence]
#         correction = se + sp - 1       [Rogan-Gladen denominator]
#         z = qnorm(0.975) ≈ 1.959964
#
# Case 1  (perfect test, p=0.30, moe=0.05):
#   p_app = 0.30, correction = 1
#   n_cont = 1.959964^2 * 0.30 * 0.70 / 0.05^2 = 322.68 → ceiling → 323
#
# Case 2  (se=0.90, sp=0.95, p=0.30, moe=0.05):
#   p_app = 0.30*0.90 + 0.70*0.05 = 0.305, correction = 0.85
#   n_cont = 1.959964^2 * 0.305 * 0.695 / (0.0025 * 0.7225) = 450.82 → ceiling → 451
#
# Case 3  (n_per_site=10, icc=0.05, p=0.30, moe=0.05):
#   deff = 1 + (10-1)*0.05 = 1.45
#   n_cont = 322.68 * 1.45 = 467.89 → ceiling → 468
#   n_sites = ceiling(468 / 10) = 47
#
# Case 4  (n_sites=50, icc=0.05, p=0.30, moe=0.05):
#   denom = 50 - 322.68*0.05 = 33.87
#   n_cont = 322.68 * 50 * 0.95 / 33.87 = 452.6 → ceiling → 453
#   n_per_site = ceiling(453 / 50) = 10
#
# Case 5  (n_sites=10, icc=0.05, p=0.30, moe=0.05):
#   denom = 10 - 322.68*0.05 = -6.13 < 0 → infeasible
#   min_moe = z * sqrt(0.3*0.7*0.05 / (10*1)) = 1.959964*sqrt(0.00105) ≈ 6.35%
#
# Case 6  (FPC, fpc_N=500, SRS):
#   n_cont = 322.68 (SRS base)
#   n_adj  = 322.68 * 500 / (322.68 + 500 - 1) = 161340 / 821.68 = 196.35 → ceiling → 197

test_that("Case 1: perfect test, SRS returns correct n and deff", {
  res <- design_precision(prevalence = 0.3, moe = 0.05)

  expect_equal(res$n,            323)
  expect_equal(res$deff,         1)
  expect_equal(res$apparent_prev, 0.3, tolerance = 1e-6)
  expect_null(res$n_sites)
  expect_null(res$n_per_site)
})

test_that("Case 2: imperfect test inflates n via Rogan-Gladen variance", {
  res <- design_precision(0.3, 0.05, sensitivity = 0.9, specificity = 0.95)

  # apparent prevalence = 0.30*0.90 + 0.70*0.05 = 0.305
  expect_equal(res$apparent_prev, 0.305, tolerance = 1e-6)
  expect_equal(res$n, 451)
  expect_gt(res$n, 323)   # imperfect test always needs more samples
})

test_that("Case 3: fixed n_per_site inflates n by deff = 1.45", {
  res <- design_precision(0.3, 0.05, n_per_site = 10, icc = 0.05)

  expect_equal(res$deff, 1.45, tolerance = 1e-6)
  expect_equal(res$n, 468)
  expect_gt(res$n, 323)

  # sites needed = ceiling(468 / 10) = 47
  expect_equal(res$n_sites,    47)
  expect_equal(res$n_per_site, 10)
})

test_that("Case 4: fixed n_sites resolves circularity via closed-form", {
  res <- design_precision(0.3, 0.05, n_sites = 50, icc = 0.05)

  expect_equal(res$n,         453)
  expect_equal(res$n_sites,    50)
  expect_equal(res$n_per_site, 10)   # ceiling(453/50) = 10
})

test_that("Case 5: infeasible n_sites produces informative error", {
  # denom = 10 - 322.68*0.05 < 0 → unachievable
  expect_error(
    design_precision(0.3, 0.05, n_sites = 10, icc = 0.05),
    "unachievable"
  )
})

test_that("Case 6: FPC reduces required n for a small population", {
  res_fpc  <- design_precision(0.3, 0.05, fpc_N = 500)
  res_nofpc <- design_precision(0.3, 0.05)

  # FPC-adjusted n should be smaller
  expect_lt(res_fpc$n, res_nofpc$n)
  # Hand-checked: n_adj = 322.68*500/(322.68+500-1) = 196.35 → 197
  expect_equal(res_fpc$n, 197)
  expect_equal(res_fpc$fpc_N, 500)
})

test_that("return list contains all expected fields", {
  res <- design_precision(0.2, 0.05)
  expect_named(res, c("n", "n_eff", "n_sites", "n_per_site", "prevalence",
                       "apparent_prev", "moe", "conf_level", "sensitivity",
                       "specificity", "icc", "deff", "fpc_N"),
               ignore.order = FALSE)
})

test_that("input validation catches bad arguments", {
  expect_error(design_precision(0,   0.05),  "`prevalence`")
  expect_error(design_precision(0.3, 0),     "`moe`")
  expect_error(design_precision(0.3, 0.05, sensitivity = 0.2, specificity = 0.2),
               "must exceed 1")
  expect_error(design_precision(0.3, 0.05, icc = 0.05),
               "cluster structure")
  expect_error(design_precision(0.3, 0.05, n_sites = 50, n_per_site = 10, icc = 0.05),
               "at most one")
  # conf_level out of range
  expect_error(design_precision(0.3, 0.05, conf_level = 0),  "`conf_level`")
  expect_error(design_precision(0.3, 0.05, conf_level = 1),  "`conf_level`")
  # fpc_N must be positive
  expect_error(design_precision(0.3, 0.05, fpc_N =  0),  "`fpc_N`")
  expect_error(design_precision(0.3, 0.05, fpc_N = -50), "`fpc_N`")
  # n_per_site must be a positive integer
  expect_error(design_precision(0.3, 0.05, n_per_site = 0,    icc = 0.05), "`n_per_site`")
  expect_error(design_precision(0.3, 0.05, n_per_site = 10.7, icc = 0.05), "`n_per_site`")
  # n_sites must be a positive integer
  expect_error(design_precision(0.3, 0.05, n_sites = 0,    icc = 0.05), "`n_sites`")
  expect_error(design_precision(0.3, 0.05, n_sites = 50.5, icc = 0.05), "`n_sites`")
})

# ---- Round 4: 15 new edge cases ----

test_that("DP-1: prevalence=0.5 (max variance) gives largest n", {
  # n = ceiling(1.96^2 * 0.25 / 0.05^2) = ceiling(384.16) = 385
  r50 <- design_precision(0.5, 0.05)
  r30 <- design_precision(0.3, 0.05)
  expect_equal(r50$n, 385)
  expect_gt(r50$n, r30$n)
})

test_that("DP-2: moe=0.499 (nearly maximum) → tiny n", {
  # Extremely wide CI: n should be very small (< 10)
  res <- design_precision(0.3, 0.499)
  expect_lt(res$n, 10)
  expect_gte(res$n, 1)
})

test_that("DP-3: conf_level=0.99 requires more samples than 0.95", {
  r99 <- design_precision(0.3, 0.05, conf_level = 0.99)
  r95 <- design_precision(0.3, 0.05)
  expect_gt(r99$n, r95$n)
  # n scales with z^2: ratio ≈ (2.576/1.960)^2 = 1.727
  expect_equal(r99$n / r95$n, (qnorm(0.995) / qnorm(0.975))^2,
               tolerance = 0.02)
})

test_that("DP-4: n_per_site=1 with icc>0 → deff=1 (no clustering when cluster size=1)", {
  # deff = 1 + (1-1)*icc = 1 regardless of icc
  res <- design_precision(0.3, 0.05, n_per_site = 1, icc = 0.05)
  expect_equal(res$deff, 1, tolerance = 1e-10)
  expect_equal(res$n, design_precision(0.3, 0.05)$n)  # same as SRS
})

test_that("DP-5: n_sites > n_base produces deff<=1 error (impossible design)", {
  # n_base ≈ 323; n_sites=400 → closed-form gives deff < 1 → error
  expect_error(
    design_precision(0.3, 0.05, n_sites = 400, icc = 0.05),
    "SRS sample size"
  )
})

test_that("DP-6: very large fpc_N has negligible effect on n", {
  r_fpc <- design_precision(0.3, 0.05, fpc_N = 1e8)
  r_srs <- design_precision(0.3, 0.05)
  # FPC factor ≈ 1 − n/(2N) → barely changes n
  expect_equal(r_fpc$n, r_srs$n)
})

test_that("DP-7: near-perfect se=sp=0.999 barely inflates n above perfect test", {
  # correction = 0.998 vs 1.000: n ≈ 323 / 0.998^2 ≈ 324
  r_perfect <- design_precision(0.3, 0.05)
  r_near    <- design_precision(0.3, 0.05, sensitivity = 0.999, specificity = 0.999)
  expect_gt(r_near$n, r_perfect$n)
  expect_lt(r_near$n, r_perfect$n + 5)  # inflation should be tiny
})

test_that("DP-8: moe=NA → informative finite-value error", {
  expect_error(design_precision(0.3, NA), "`moe`")
})

test_that("DP-9: prevalence=Inf → informative finite-value error", {
  expect_error(design_precision(Inf, 0.05), "`prevalence`")
})

test_that("DP-10: n_sites=1 with icc=0 → SRS branch returns n=323, n_sites=1", {
  # icc=0 triggers SRS branch; distribution step still assigns n_sites=1
  res <- design_precision(0.3, 0.05, n_sites = 1, icc = 0)
  expect_equal(res$n,        323)
  expect_equal(res$deff,     1,   tolerance = 1e-10)
  expect_equal(res$n_sites,  1)
})

test_that("DP-11: icc=0 with n_sites=50 → SRS n distributed across 50 sites", {
  # icc=0 dominates: deff=1, n=323; but n_sites=50 is used in distribution
  res <- design_precision(0.3, 0.05, n_sites = 50, icc = 0)
  expect_equal(res$n,         323)
  expect_equal(res$deff,      1,   tolerance = 1e-10)
  expect_equal(res$n_sites,   50)
  expect_equal(res$n_per_site, ceiling(323 / 50))
})

test_that("DP-12: n_per_site=1L (integer type) works identically to double", {
  r_int <- design_precision(0.3, 0.05, n_per_site = 1L,  icc = 0.05)
  r_dbl <- design_precision(0.3, 0.05, n_per_site = 1.0, icc = 0.05)
  expect_equal(r_int$n,    r_dbl$n)
  expect_equal(r_int$deff, r_dbl$deff)
})

test_that("DP-13: very small moe=0.001 → very large n, scales correctly", {
  # n ∝ 1/moe^2; ratio vs moe=0.05: (0.05/0.001)^2 = 2500
  r_large <- design_precision(0.3, 0.001)
  r_small <- design_precision(0.3, 0.05)
  expect_equal(r_large$n / r_small$n, 2500, tolerance = 0.01)
})

test_that("DP-14: clustering + FPC combined: n_per_site + fpc_N both reduce final n", {
  # Clustered alone
  r_cluster <- design_precision(0.3, 0.05, n_per_site = 10, icc = 0.05)
  # Clustered + FPC (population = 1000)
  r_both    <- design_precision(0.3, 0.05, n_per_site = 10, icc = 0.05, fpc_N = 1000)
  # FPC should reduce n below the clustering-only case
  expect_lt(r_both$n, r_cluster$n)
  expect_equal(r_both$fpc_N, 1000)
})

test_that("DP-15: very low prevalence (p=0.001) returns valid n and apparent_prev", {
  res <- design_precision(0.001, 0.001)
  expect_true(is.finite(res$n))
  expect_gt(res$n, 0)
  expect_equal(res$apparent_prev, 0.001, tolerance = 1e-6)  # perfect test
})

# ---- Round 5: 7 new edge cases ----

test_that("DP-R5-1: sensitivity=0 is rejected (boundary, not in (0,1])", {
  expect_error(design_precision(0.3, 0.05, sensitivity = 0), "`sensitivity`")
})

test_that("DP-R5-2: n_sites as vector gives informative length error", {
  expect_error(
    design_precision(0.3, 0.05, n_sites = c(10, 20), icc = 0.05),
    "single finite positive integer"
  )
})

test_that("DP-R5-3: conf_level=0.5 produces a small n (low confidence threshold)", {
  # z_0.5 = qnorm(0.75) ≈ 0.674; n scales with z^2
  r_50 <- design_precision(0.3, 0.05, conf_level = 0.50)
  r_95 <- design_precision(0.3, 0.05)
  expect_lt(r_50$n, r_95$n)
  expect_equal(r_50$n / r_95$n, (qnorm(0.75) / qnorm(0.975))^2, tolerance = 0.02)
})

test_that("DP-R5-4: icc=0.999, n_per_site=2 → deff≈2, n≈2*n_base", {
  res <- design_precision(0.3, 0.05, n_per_site = 2, icc = 0.999)
  expect_equal(res$deff, 1.999, tolerance = 1e-3)
  expect_equal(res$n, ceiling(design_precision(0.3, 0.05)$n * 1.999), tolerance = 1)
})

test_that("DP-R5-5: moe=0.5 is rejected (boundary, must be strictly < 0.5)", {
  expect_error(design_precision(0.3, 0.5), "`moe`")
})

test_that("DP-R5-6: prevalence vector is rejected with length error", {
  expect_error(design_precision(c(0.2, 0.3), 0.05), "`prevalence`")
})

test_that("DP-R5-7: n_per_site given but icc=0 → deff=1, same n as SRS", {
  # deff = 1 + (n_per_site - 1)*0 = 1 regardless of cluster size
  res <- design_precision(0.3, 0.05, n_per_site = 50, icc = 0)
  expect_equal(res$deff, 1, tolerance = 1e-10)
  expect_equal(res$n, design_precision(0.3, 0.05)$n)
})

# ---- Round 6: 15 new edge cases ----

test_that("DP-R6-1: moe=-0.05 is rejected (must be > 0)", {
  expect_error(design_precision(0.3, -0.05), "`moe`")
})

test_that("DP-R6-2: sensitivity=1.001 is rejected (must be <= 1)", {
  expect_error(design_precision(0.3, 0.05, sensitivity = 1.001), "`sensitivity`")
})

test_that("DP-R6-3: specificity=NA is rejected by is.finite check", {
  expect_error(design_precision(0.3, 0.05, specificity = NA), "`specificity`")
})

test_that("DP-R6-4: n_sites=1 with icc=0.5 is infeasible (denom < 0)", {
  # denom = 1 - n_base*0.5 ≈ 1 - 161 = -160 < 0 → unachievable
  expect_error(
    design_precision(0.3, 0.05, n_sites = 1, icc = 0.5),
    "unachievable"
  )
})

test_that("DP-R6-5: large n_per_site with high icc → very large n", {
  # deff = 1 + (10000-1)*0.3 = 3000.7; n ≈ 323 * 3000.7 ≈ 969k
  res <- design_precision(0.3, 0.05, n_per_site = 10000, icc = 0.3)
  expect_equal(res$deff, 3000.7, tolerance = 0.1)
  expect_gt(res$n, 900000)
})

test_that("DP-R6-6: icc as vector is rejected with length error", {
  expect_error(
    design_precision(0.3, 0.05, n_per_site = 10, icc = c(0.05, 0.1)),
    "`icc`"
  )
})

test_that("DP-R6-7: fpc_N=1 (population of 1) → n=1 regardless of SRS n", {
  # n_adj = n_cont * 1 / (n_cont + 1 - 1) = 1.0 exactly
  res <- design_precision(0.3, 0.05, fpc_N = 1)
  expect_equal(res$n, 1)
})

test_that("DP-R6-8: prevalence=-0.01 is rejected", {
  expect_error(design_precision(-0.01, 0.05), "`prevalence`")
})

test_that("DP-R6-9: wide moe=0.3 gives very small n", {
  # n = ceil(1.96^2 * 0.3*0.7 / 0.3^2) = ceil(3.84*0.21/0.09) = ceil(8.96) = 9
  res <- design_precision(0.3, 0.3)
  expect_equal(res$n, 9)
})

test_that("DP-R6-10: n_sites=323 (= ceiling n_base) triggers deff<=1 error", {
  # n_base_cont ≈ 322.68; n_sites=323 >= n_base → deff < 1 → impossible design
  expect_error(
    design_precision(0.3, 0.05, n_sites = 323, icc = 0.05),
    "SRS sample size"
  )
})

test_that("DP-R6-11: fpc_N=2 → n=2 (census of a 2-person population)", {
  res <- design_precision(0.3, 0.05, fpc_N = 2)
  expect_equal(res$n, 2)
})

test_that("DP-R6-12: conf_level=NA is rejected by is.finite check", {
  expect_error(design_precision(0.3, 0.05, conf_level = NA), "`conf_level`")
})

test_that("DP-R6-13: n_per_site as vector is rejected with length error", {
  expect_error(
    design_precision(0.3, 0.05, n_per_site = c(5, 10), icc = 0.05),
    "single finite positive integer"
  )
})

test_that("DP-R6-14: n_sites + fpc_N together is valid (not mutually exclusive)", {
  res <- design_precision(0.3, 0.05, n_sites = 50, icc = 0.05, fpc_N = 1000)
  expect_true(is.finite(res$n))
  expect_equal(res$n_sites, 50)
  expect_equal(res$fpc_N,   1000)
  # FPC reduces n below the non-FPC clustered case
  res_no_fpc <- design_precision(0.3, 0.05, n_sites = 50, icc = 0.05)
  expect_lt(res$n, res_no_fpc$n)
})

test_that("DP-R6-15: prevalence=0.9999 (near-boundary) returns finite n", {
  # p_app ≈ 0.9999; p*(1-p) ≈ 0.0001 → tiny variance → very small n
  res <- design_precision(0.9999, 0.05)
  expect_true(is.finite(res$n))
  expect_gte(res$n, 1)
  expect_equal(res$apparent_prev, 0.9999, tolerance = 1e-6)
})

test_that("DP-R7-1: character prevalence/moe give a friendly class error", {
  expect_error(design_precision(prevalence = "0.3", moe = 0.05), "`prevalence`")
  expect_error(design_precision(prevalence = 0.3, moe = "0.05"), "`moe`")
})

test_that("DP-R7-2: logical params are rejected, not coerced", {
  expect_error(design_precision(0.3, 0.05, sensitivity = TRUE), "`sensitivity`")
  expect_error(design_precision(0.3, 0.05, icc = FALSE), "`icc`")
})

test_that("DP-R7-3: list-valued conf_level gives a friendly class error", {
  expect_error(design_precision(0.3, 0.05, conf_level = list(0.95)), "`conf_level`")
})
