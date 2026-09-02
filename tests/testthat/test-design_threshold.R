# Hand-checked reference cases for design_threshold()
#
# Power formula: n = ((z_a * SE0 + z_b * SE1) / delta)^2
#   SE0   = sqrt(theta_app * (1-theta_app))   [null SE per sqrt(n)]
#   SE1   = sqrt(p1_app * (1-p1_app))         [alternative SE]
#   delta = |p1_app - theta_app|
#   z_a   = qnorm(1-alpha) for one-sided, qnorm(1-alpha/2) for two-sided
#   z_b   = qnorm(power)
#
# Case 1 (threshold=0.05, prevalence=0.10, power=0.80, one-sided "greater"):
#   alpha=0.05, z_a=qnorm(0.95)=1.6449, z_b=qnorm(0.80)=0.8416
#   theta_app=0.05, p1_app=0.10 (perfect test)
#   SE0 = sqrt(0.05*0.95) = 0.2179, SE1 = sqrt(0.10*0.90) = 0.3000
#   delta = 0.05
#   n = ((1.6449*0.2179 + 0.8416*0.3000) / 0.05)^2
#     = ((0.3583 + 0.2525) / 0.05)^2 = (0.6108/0.05)^2 = 12.216^2 = 149.2 -> 150
#
# Case 2 (same, power=0.90):
#   z_b = qnorm(0.90) = 1.2816
#   n = ((1.6449*0.2179 + 1.2816*0.3000) / 0.05)^2
#     = ((0.3583 + 0.3845) / 0.05)^2 = (0.7428/0.05)^2 = 14.856^2 = 220.7 -> 221
#
# Case 3 (two-sided, threshold=0.05, prevalence=0.10, power=0.80):
#   z_a = qnorm(0.975) = 1.9600
#   n = ((1.9600*0.2179 + 0.8416*0.3000) / 0.05)^2
#     = ((0.4271 + 0.2525) / 0.05)^2 = (0.6796/0.05)^2 = 13.592^2 = 184.7 -> 185

test_that("Case 1: SRS, one-sided 'greater', 80% power -- hand-checked n=150", {
  res <- design_threshold(threshold = 0.05, prevalence = 0.10)

  expect_equal(res$threshold,     0.05)
  expect_equal(res$prevalence,    0.10)
  expect_equal(res$power,         0.80)
  expect_equal(res$alternative,   "greater")
  expect_equal(res$deff,          1,    tolerance = 1e-6)
  expect_null(res$n_sites)
  expect_null(res$n_per_site)

  # Hand-computed
  z_a <- qnorm(0.95); z_b <- qnorm(0.80)
  n_expected <- ceiling(((z_a * sqrt(0.05*0.95) + z_b * sqrt(0.10*0.90)) / 0.05)^2)
  expect_equal(res$n, n_expected)
})

test_that("Case 2: power=0.90 requires more samples than power=0.80", {
  r80 <- design_threshold(threshold = 0.05, prevalence = 0.10, power = 0.80)
  r90 <- design_threshold(threshold = 0.05, prevalence = 0.10, power = 0.90)
  expect_gt(r90$n, r80$n)

  z_a <- qnorm(0.95); z_b90 <- qnorm(0.90)
  n_expected <- ceiling(((z_a*sqrt(0.05*0.95) + z_b90*sqrt(0.10*0.90)) / 0.05)^2)
  expect_equal(r90$n, n_expected)
})

test_that("Case 3: two-sided test uses z_{alpha/2} and requires more than one-sided", {
  r_1s <- design_threshold(threshold = 0.05, prevalence = 0.10, alternative = "greater")
  r_2s <- design_threshold(threshold = 0.05, prevalence = 0.10, alternative = "two.sided")
  expect_gt(r_2s$n, r_1s$n)

  z_a2 <- qnorm(0.975); z_b <- qnorm(0.80)
  n_expected <- ceiling(((z_a2*sqrt(0.05*0.95) + z_b*sqrt(0.10*0.90)) / 0.05)^2)
  expect_equal(r_2s$n, n_expected)
})

test_that("n scales correctly with effect size: larger gap needs fewer samples", {
  r_small <- design_threshold(threshold = 0.05, prevalence = 0.07)   # gap=0.02
  r_large <- design_threshold(threshold = 0.05, prevalence = 0.15)   # gap=0.10
  expect_gt(r_small$n, r_large$n)
})

test_that("larger conf_level (smaller alpha) requires more samples", {
  r95 <- design_threshold(threshold = 0.05, prevalence = 0.10)
  r99 <- design_threshold(threshold = 0.05, prevalence = 0.10, conf_level = 0.99)
  expect_gt(r99$n, r95$n)
})

test_that("perfect test: apparent_prev = prevalence, threshold_app = threshold", {
  res <- design_threshold(threshold = 0.05, prevalence = 0.10)
  expect_equal(res$apparent_prev, 0.10, tolerance = 1e-10)
  expect_equal(res$threshold_app, 0.05, tolerance = 1e-10)
})

test_that("imperfect test shifts apparent prevalences and changes n", {
  r_perfect  <- design_threshold(threshold = 0.05, prevalence = 0.10)
  r_imperfect <- design_threshold(threshold = 0.05, prevalence = 0.10,
                                   sensitivity = 0.9, specificity = 0.95)
  # apparent values differ from true
  expect_false(isTRUE(all.equal(r_imperfect$apparent_prev, 0.10)))
  expect_false(isTRUE(all.equal(r_imperfect$threshold_app, 0.05)))
  # n may increase or decrease depending on the apparent delta
  expect_true(is.finite(r_imperfect$n))
  expect_gte(r_imperfect$n, 1)
})

test_that("alternative='less': prevalence < threshold, hand-checked n", {
  # threshold=0.10, prevalence=0.05, alternative="less"
  # SE0=sqrt(0.10*0.90)=0.3000, SE1=sqrt(0.05*0.95)=0.2179, delta=0.05
  res <- design_threshold(threshold = 0.10, prevalence = 0.05, alternative = "less")
  z_a <- qnorm(0.95); z_b <- qnorm(0.80)
  n_expected <- ceiling(((z_a*sqrt(0.10*0.90) + z_b*sqrt(0.05*0.95)) / 0.05)^2)
  expect_equal(res$n, n_expected)
  expect_equal(res$alternative, "less")
})

test_that("clustering with n_per_site inflates n via deff", {
  r_srs     <- design_threshold(threshold = 0.05, prevalence = 0.10)
  r_cluster <- design_threshold(threshold = 0.05, prevalence = 0.10,
                                 n_per_site = 10, icc = 0.05)
  expect_equal(r_cluster$deff, 1.45, tolerance = 1e-6)
  expect_gt(r_cluster$n, r_srs$n)
  expect_equal(r_cluster$n_sites, ceiling(r_cluster$n / 10))
})

test_that("clustering with n_sites: closed-form produces consistent deff", {
  res <- design_threshold(threshold = 0.05, prevalence = 0.10,
                           n_sites = 50, icc = 0.05)
  expect_equal(res$n_sites,    50)
  expect_gt(res$deff,          1)
  expect_equal(res$n_per_site, ceiling(res$n / 50))
})

test_that("FPC reduces required n", {
  r_nofpc <- design_threshold(threshold = 0.05, prevalence = 0.10)
  r_fpc   <- design_threshold(threshold = 0.05, prevalence = 0.10, fpc_N = 500)
  expect_lt(r_fpc$n, r_nofpc$n)
  expect_equal(r_fpc$fpc_N, 500)
})

test_that("return list has all expected fields in order (SRS)", {
  res <- design_threshold(threshold = 0.05, prevalence = 0.10)
  expect_named(res, c("n", "n_eff", "n_sites", "n_per_site",
                       "threshold", "prevalence", "apparent_prev", "threshold_app",
                       "power", "alternative", "conf_level",
                       "sensitivity", "specificity", "icc", "deff", "fpc_N"),
               ignore.order = FALSE)
})

test_that("n_eff is the SRS-equivalent size: <= n, and n/deff up to rounding", {
  res <- design_threshold(threshold = 0.05, prevalence = 0.10,
                           n_per_site = 10, icc = 0.05)
  expect_lte(res$n_eff, res$n)
  # n_eff = ceiling(n_base_cont); n = ceiling(n_base_cont * deff); they agree
  # to within one ceiling step per side
  expect_equal(res$n_eff, res$n / res$deff, tolerance = 1)
  # SRS: n_eff == n exactly
  srs <- design_threshold(threshold = 0.05, prevalence = 0.10)
  expect_equal(srs$n_eff, srs$n)
})

# ---- validation ----

test_that("validation: threshold=0 is rejected", {
  expect_error(design_threshold(threshold = 0, prevalence = 0.10), "`threshold`")
})

test_that("validation: threshold=1 is rejected", {
  expect_error(design_threshold(threshold = 1, prevalence = 0.10), "`threshold`")
})

test_that("validation: prevalence=0 is rejected", {
  expect_error(design_threshold(threshold = 0.05, prevalence = 0), "`prevalence`")
})

test_that("validation: prevalence <= threshold for 'greater' is rejected", {
  expect_error(
    design_threshold(threshold = 0.10, prevalence = 0.10, alternative = "greater"),
    "must be >"
  )
  expect_error(
    design_threshold(threshold = 0.10, prevalence = 0.05, alternative = "greater"),
    "must be >"
  )
})

test_that("validation: prevalence >= threshold for 'less' is rejected", {
  expect_error(
    design_threshold(threshold = 0.10, prevalence = 0.10, alternative = "less"),
    "must be <"
  )
  expect_error(
    design_threshold(threshold = 0.10, prevalence = 0.15, alternative = "less"),
    "must be <"
  )
})

test_that("validation: power=0 is rejected", {
  expect_error(design_threshold(threshold = 0.05, prevalence = 0.10, power = 0),
               "`power`")
})

test_that("validation: power=1 is rejected", {
  expect_error(design_threshold(threshold = 0.05, prevalence = 0.10, power = 1),
               "`power`")
})

test_that("validation: invalid alternative is rejected", {
  expect_error(
    design_threshold(threshold = 0.05, prevalence = 0.10, alternative = "up"),
    "'greater', 'less', or 'two.sided'"
  )
})

test_that("validation: se+sp <= 1 is rejected", {
  expect_error(
    design_threshold(threshold = 0.05, prevalence = 0.10,
                     sensitivity = 0.4, specificity = 0.4),
    "must exceed 1"
  )
})

test_that("validation: icc > 0 without clustering args is rejected", {
  expect_error(
    design_threshold(threshold = 0.05, prevalence = 0.10, icc = 0.05),
    "cluster structure"
  )
})

test_that("validation: n_sites and n_per_site both supplied is rejected", {
  expect_error(
    design_threshold(threshold = 0.05, prevalence = 0.10,
                     n_sites = 20, n_per_site = 10, icc = 0.05),
    "at most one"
  )
})

test_that("validation: threshold as vector is rejected", {
  expect_error(
    design_threshold(threshold = c(0.05, 0.10), prevalence = 0.15),
    "`threshold`"
  )
})

test_that("validation: prevalence as vector is rejected", {
  expect_error(
    design_threshold(threshold = 0.05, prevalence = c(0.10, 0.15)),
    "`prevalence`"
  )
})

test_that("validation: power as vector is rejected", {
  expect_error(
    design_threshold(threshold = 0.05, prevalence = 0.10, power = c(0.8, 0.9)),
    "`power`"
  )
})

test_that("validation: conf_level out of range is rejected", {
  expect_error(
    design_threshold(threshold = 0.05, prevalence = 0.10, conf_level = 0),
    "`conf_level`"
  )
  expect_error(
    design_threshold(threshold = 0.05, prevalence = 0.10, conf_level = 1),
    "`conf_level`"
  )
})

test_that("validation: icc out of range is rejected", {
  expect_error(
    design_threshold(threshold = 0.05, prevalence = 0.10,
                     n_per_site = 10, icc = 1.5),
    "`icc`"
  )
})

test_that("validation: fpc_N as vector is rejected", {
  expect_error(
    design_threshold(threshold = 0.05, prevalence = 0.10, fpc_N = c(100, 200)),
    "`fpc_N`"
  )
})

# ---- edge cases ----

test_that("prevalence very close to threshold: n grows very large", {
  r_far   <- design_threshold(threshold = 0.05, prevalence = 0.15)
  r_close <- design_threshold(threshold = 0.05, prevalence = 0.06)
  expect_gt(r_close$n, r_far$n)
  expect_gt(r_close$n, 1000)
})

test_that("power=0.999 (very high) requires much more than power=0.80", {
  r80   <- design_threshold(threshold = 0.05, prevalence = 0.10)
  r999  <- design_threshold(threshold = 0.05, prevalence = 0.10, power = 0.999)
  expect_gt(r999$n, r80$n * 2)
})

test_that("n_per_site=1: deff=1, same as SRS", {
  r_srs <- design_threshold(threshold = 0.05, prevalence = 0.10)
  r_1ps <- design_threshold(threshold = 0.05, prevalence = 0.10,
                              n_per_site = 1, icc = 0.05)
  expect_equal(r_1ps$deff, 1, tolerance = 1e-10)
  expect_equal(r_1ps$n,    r_srs$n)
})

test_that("infeasible n_sites with high ICC gives informative error", {
  # Need to find an n_sites that's too small. n_base for threshold=0.05, prev=0.10 ~= 150
  # n_sites=5 with icc=0.5: denom = 5 - 150*0.5 = -70 < 0
  expect_error(
    design_threshold(threshold = 0.05, prevalence = 0.10,
                     n_sites = 5, icc = 0.5),
    "unachievable"
  )
})

test_that("very large fpc_N barely changes n", {
  r_srs <- design_threshold(threshold = 0.05, prevalence = 0.10)
  r_fpc <- design_threshold(threshold = 0.05, prevalence = 0.10, fpc_N = 1e8)
  expect_equal(r_fpc$n, r_srs$n)
})

test_that("fpc_N=1: census of 1 person -> n=1", {
  res <- design_threshold(threshold = 0.05, prevalence = 0.10, fpc_N = 1)
  expect_equal(res$n, 1)
})

test_that("'less' and 'greater' symmetric around 0.5 give same n", {
  # threshold=0.50: p*(1-p) is symmetric, so p=0.45 and p=0.55 give same SE1
  # SE1(0.55) = sqrt(0.55*0.45) = SE1(0.45) = sqrt(0.45*0.55)
  r_gt <- design_threshold(threshold = 0.50, prevalence = 0.55, alternative = "greater")
  r_lt <- design_threshold(threshold = 0.50, prevalence = 0.45, alternative = "less")
  expect_equal(r_gt$n, r_lt$n)
})

test_that("two.sided n > one-sided n for same effect size", {
  r_1s <- design_threshold(threshold = 0.05, prevalence = 0.10, alternative = "greater")
  r_2s <- design_threshold(threshold = 0.05, prevalence = 0.10, alternative = "two.sided")
  expect_gt(r_2s$n, r_1s$n)
})

# ---- stress tests: 15 edge cases to break design_threshold ----

test_that("DD-S1: threshold=0.001 (near-zero boundary): valid finite n", {
  res <- design_threshold(threshold = 0.001, prevalence = 0.01)
  expect_true(is.finite(res$n))
  expect_gte(res$n, 1)
})

test_that("DD-S2: threshold=0.999 (near-one boundary): valid finite n (less)", {
  res <- design_threshold(threshold = 0.999, prevalence = 0.99, alternative = "less")
  expect_true(is.finite(res$n))
  expect_gte(res$n, 1)
})

test_that("DD-S3: power=0.999: very large n, still finite", {
  res <- design_threshold(threshold = 0.05, prevalence = 0.10, power = 0.999)
  expect_true(is.finite(res$n))
  r80 <- design_threshold(threshold = 0.05, prevalence = 0.10, power = 0.80)
  expect_gt(res$n, r80$n * 2)
})

test_that("DD-S4: power < 0.5 is rejected (formula breaks down below 0.5)", {
  expect_error(
    design_threshold(threshold = 0.05, prevalence = 0.10, power = 0.001),
    "at least 0.5"
  )
  expect_error(
    design_threshold(threshold = 0.05, prevalence = 0.10, power = 0.49),
    "at least 0.5"
  )
  # exactly 0.5 is allowed (z_beta = 0)
  res <- design_threshold(threshold = 0.05, prevalence = 0.10, power = 0.5)
  expect_true(is.finite(res$n) && res$n >= 1)
})

test_that("DD-S5: prevalence very close to threshold: large n", {
  r_far   <- design_threshold(threshold = 0.10, prevalence = 0.20)
  r_close <- design_threshold(threshold = 0.10, prevalence = 0.101)
  expect_gt(r_close$n, r_far$n)
  expect_gt(r_close$n, 10000)
})

test_that("DD-S6: se+sp just above 1: warning emitted, n finite", {
  expect_warning(
    res <- design_threshold(threshold = 0.05, prevalence = 0.10,
                             sensitivity = 0.51, specificity = 0.50),
    "numerically unstable"
  )
  expect_true(is.finite(res$n))
})

test_that("DD-S7: n_per_site=1, icc>0 -> deff=1, same as SRS", {
  r_srs <- design_threshold(threshold = 0.05, prevalence = 0.10)
  r_1ps <- design_threshold(threshold = 0.05, prevalence = 0.10,
                              n_per_site = 1, icc = 0.5)
  expect_equal(r_1ps$deff, 1, tolerance = 1e-10)
  expect_equal(r_1ps$n, r_srs$n)
})

test_that("DD-S8: icc=1, n_per_site=10 -> deff=10 (maximum)", {
  res <- design_threshold(threshold = 0.05, prevalence = 0.10,
                           n_per_site = 10, icc = 1)
  expect_equal(res$deff, 10, tolerance = 1e-6)
})

test_that("DD-S9: imperfect test shifts apparent_prev and threshold_app", {
  # se=0.9, sp=0.8
  res <- design_threshold(threshold = 0.05, prevalence = 0.20,
                           sensitivity = 0.9, specificity = 0.8)
  expect_true(is.finite(res$n))
  expect_equal(res$threshold_app, 0.05*0.9 + 0.95*0.2, tolerance = 1e-10)
  expect_equal(res$apparent_prev, 0.20*0.9 + 0.80*0.2, tolerance = 1e-10)
})

test_that("DD-S10: n_eff <= n and ~ n/deff in all three design modes", {
  r_srs <- design_threshold(threshold = 0.05, prevalence = 0.10)
  r_nps <- design_threshold(threshold = 0.05, prevalence = 0.10,
                              n_per_site = 10, icc = 0.05)
  r_ns  <- design_threshold(threshold = 0.05, prevalence = 0.10,
                              n_sites = 50, icc = 0.05)
  expect_equal(r_srs$n_eff, r_srs$n)          # SRS: exact
  expect_lte(r_nps$n_eff, r_nps$n)
  expect_lte(r_ns$n_eff,  r_ns$n)
  expect_equal(r_nps$n_eff, r_nps$n / r_nps$deff, tolerance = 1)
  expect_equal(r_ns$n_eff,  r_ns$n  / r_ns$deff,  tolerance = 1)
})

test_that("DD-S11: conf_level=0.50 (alpha=0.50) requires fewer samples than 0.95", {
  r95 <- design_threshold(threshold = 0.05, prevalence = 0.10)
  r50 <- design_threshold(threshold = 0.05, prevalence = 0.10, conf_level = 0.50)
  expect_lt(r50$n, r95$n)
})

test_that("DD-S12: fpc_N=2: n=1 or 2", {
  res <- design_threshold(threshold = 0.05, prevalence = 0.10, fpc_N = 2)
  expect_lte(res$n, 2)
  expect_gte(res$n, 1)
})

test_that("DD-S13: NA threshold is rejected", {
  expect_error(design_threshold(threshold = NA, prevalence = 0.10), "`threshold`")
})

test_that("DD-S14: NA prevalence is rejected", {
  expect_error(design_threshold(threshold = 0.05, prevalence = NA), "`prevalence`")
})

test_that("DD-S15: n_sites fixed at n_base or above triggers deff<=1 error", {
  r_srs <- design_threshold(threshold = 0.05, prevalence = 0.10)
  expect_error(
    design_threshold(threshold = 0.05, prevalence = 0.10,
                     n_sites = r_srs$n, icc = 0.05),
    "SRS sample size|< 1 person"
  )
})

# -- Second batch of stress tests: DD-S16 through DD-S30 ----------------------

test_that("DD-S16: prevalence barely above threshold -> n much larger than farther case", {
  r_far  <- design_threshold(threshold = 0.10, prevalence = 0.20)
  r_near <- design_threshold(threshold = 0.10, prevalence = 0.101)
  expect_gt(r_near$n, r_far$n * 50)  # tiny delta -> enormous n
})

test_that("DD-S17: power=0.50 (z_b=0) -> n governed by type-I term only", {
  res <- design_threshold(threshold = 0.10, prevalence = 0.20, power = 0.50)
  z_a <- qnorm(0.95)
  se0 <- sqrt(0.10 * 0.90)
  n_expected <- ceiling((z_a * se0 / 0.10)^2)  # z_b=0, so only z_a*se0 term
  expect_equal(res$n, n_expected)
})

test_that("DD-S18: two-sided uses qnorm(1-alpha/2) -> n strictly larger than one-sided", {
  r_one <- design_threshold(threshold = 0.10, prevalence = 0.20, alternative = "greater")
  r_two <- design_threshold(threshold = 0.10, prevalence = 0.20, alternative = "two.sided")
  expect_gt(r_two$n, r_one$n)
})

test_that("DD-S19: perfect test: n matches closed-form power formula exactly", {
  res <- design_threshold(threshold = 0.10, prevalence = 0.20)
  z_a <- qnorm(0.95); z_b <- qnorm(0.80)
  se0 <- sqrt(0.10 * 0.90); se1 <- sqrt(0.20 * 0.80)
  n_expected <- ceiling(((z_a * se0 + z_b * se1) / 0.10)^2)
  expect_equal(res$n, n_expected)
})

test_that("DD-S20: alternative='less' with prevalence < threshold -> valid finite n", {
  res <- design_threshold(threshold = 0.20, prevalence = 0.10, alternative = "less")
  expect_gte(res$n, 1)
  expect_true(is.finite(res$n))
  expect_equal(res$alternative, "less")
})

test_that("DD-S21: n_sites far too small for target power -> error (denom <= 0)", {
  # SRS n ~= 69; n_sites=3 with icc=0.5 -> denom = 3 - 69*0.5 << 0
  expect_error(
    design_threshold(threshold = 0.10, prevalence = 0.20, n_sites = 3, icc = 0.5),
    "unachievable"
  )
})

test_that("DD-S22: n_eff is SRS-equivalent (== n for SRS, <= n clustered)", {
  r_srs <- design_threshold(threshold = 0.10, prevalence = 0.20)
  expect_equal(r_srs$n_eff, r_srs$n)

  r_cl <- design_threshold(threshold = 0.10, prevalence = 0.20,
                           n_per_site = 10, icc = 0.05)
  expect_lt(r_cl$n_eff, r_cl$n)
  expect_equal(r_cl$n_eff, r_cl$n / r_cl$deff, tolerance = 1)
})

test_that("DD-S23: FPC reduces n relative to infinite-population case", {
  r_inf <- design_threshold(threshold = 0.10, prevalence = 0.20)
  r_fpc <- design_threshold(threshold = 0.10, prevalence = 0.20, fpc_N = 200)
  expect_lt(r_fpc$n, r_inf$n)
})

test_that("DD-S24: higher target power monotonically increases required n", {
  r80 <- design_threshold(threshold = 0.10, prevalence = 0.20, power = 0.80)
  r90 <- design_threshold(threshold = 0.10, prevalence = 0.20, power = 0.90)
  r99 <- design_threshold(threshold = 0.10, prevalence = 0.20, power = 0.99)
  expect_lt(r80$n, r90$n)
  expect_lt(r90$n, r99$n)
})

test_that("DD-S25: imperfect test with small apparent-scale delta -> very large n", {
  # se=sp=0.6 compresses delta: p1_app-theta_app = (prev-thresh)*correction = 0.10*0.2 = 0.02
  res <- design_threshold(threshold = 0.50, prevalence = 0.60,
                          sensitivity = 0.6, specificity = 0.6,
                          alternative = "greater")
  expect_gt(res$n, 1000)
})

test_that("DD-S26: return list has all documented names for all design modes", {
  expected <- c("n", "n_eff", "n_sites", "n_per_site", "threshold", "prevalence",
                "apparent_prev", "threshold_app", "power", "alternative",
                "conf_level", "sensitivity", "specificity", "icc", "deff", "fpc_N")
  # SRS
  expect_named(design_threshold(threshold = 0.10, prevalence = 0.20), expected)
  # Fixed n_per_site
  expect_named(design_threshold(threshold = 0.10, prevalence = 0.20,
                                n_per_site = 10, icc = 0.05), expected)
})

test_that("DD-S27: apparent_prev and threshold_app use Rogan-Gladen forward transform", {
  res <- design_threshold(threshold = 0.10, prevalence = 0.20,
                          sensitivity = 0.80, specificity = 0.90)
  expect_equal(res$threshold_app, 0.10 * 0.80 + 0.90 * 0.10, tolerance = 1e-10)
  expect_equal(res$apparent_prev, 0.20 * 0.80 + 0.80 * 0.10, tolerance = 1e-10)
  # apparent_prev > threshold_app for "greater" (monotone RG)
  expect_gt(res$apparent_prev, res$threshold_app)
})

test_that("DD-S28: fixed n_sites path -> n_per_site = ceiling(n / n_sites)", {
  k   <- 10
  res <- design_threshold(threshold = 0.10, prevalence = 0.20,
                          n_sites = k, icc = 0.05)
  expect_equal(res$n_sites, k)
  expect_equal(res$n_per_site, ceiling(res$n / k))
})

test_that("DD-S29: fixed n_per_site path -> n_sites = ceiling(n / n_per_site)", {
  m   <- 20
  res <- design_threshold(threshold = 0.10, prevalence = 0.20,
                          n_per_site = m, icc = 0.05)
  expect_equal(res$n_per_site, m)
  expect_equal(res$n_sites, ceiling(res$n / m))
})

test_that("DD-S30: conf_level=0.90 (alpha=0.10) uses qnorm(0.90) for one-sided z_a", {
  res <- design_threshold(threshold = 0.10, prevalence = 0.20, conf_level = 0.90)
  z_a <- qnorm(0.90); z_b <- qnorm(0.80)
  se0 <- sqrt(0.10 * 0.90); se1 <- sqrt(0.20 * 0.80)
  n_expected <- ceiling(((z_a * se0 + z_b * se1) / 0.10)^2)
  expect_equal(res$n, n_expected)
})

test_that("DT-R7-1: non-numeric power/conf_level/icc give a friendly class error", {
  expect_error(design_threshold(0.05, 0.10, power = "0.8"), "`power`")
  expect_error(design_threshold(0.05, 0.10, conf_level = "0.95"), "`conf_level`")
  expect_error(design_threshold(0.05, 0.10, icc = "0"), "`icc`")
})

test_that("DT-R7-2: logical params are rejected, not coerced", {
  expect_error(design_threshold(0.05, 0.10, power = TRUE), "`power`")
  expect_error(design_threshold(0.05, 0.10, sensitivity = TRUE), "`sensitivity`")
  expect_error(design_threshold(0.05, 0.10, icc = FALSE), "`icc`")
})


# ---------------------------------------------------------------------------
# Round 8 -- fresh code-review fixes (2026-09-02)
# ---------------------------------------------------------------------------

test_that("DT-R8-1: conf_level below 0.5 is rejected (mirrors the power rule)", {
  expect_error(design_threshold(0.05, 0.10, conf_level = 0.30), "at least 0.5")
  expect_error(design_threshold(0.05, 0.10, conf_level = 0.05), "at least 0.5")
  # exactly 0.5 is still allowed (z_alpha = 0), as for power
  expect_silent(design_threshold(0.05, 0.10, conf_level = 0.50))
})

test_that("DT-R8-2: fpc_N must be a whole number", {
  expect_error(design_threshold(0.05, 0.10, fpc_N = 500.5), "integer")
  expect_silent(design_threshold(0.05, 0.10, fpc_N = 500))
})

test_that("DT-R8-3: reported deff is consistent with the returned n_per_site", {
  a <- design_threshold(0.05, 0.10, n_sites = 50, icc = 0.05)
  expect_equal(a$deff, 1 + (a$n_per_site - 1) * 0.05, tolerance = 1e-9)

  b <- design_threshold(0.05, 0.10, n_sites = 50, icc = 0.05, fpc_N = 5000)
  expect_equal(b$deff, 1 + (b$n_per_site - 1) * 0.05, tolerance = 1e-9)
})
