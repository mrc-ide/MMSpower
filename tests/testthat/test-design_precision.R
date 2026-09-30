# Hand-checked reference cases for design_precision()
#
# Formula: n = z^2 * p_app * (1 - p_app) / (moe^2 * correction^2)
#   where p_app = p*se + (1-p)*(1-sp)   [apparent prevalence]
#         correction = se + sp - 1       [Rogan-Gladen denominator]
#         z = qnorm(0.975) ~= 1.959964
#
# DP-1  (perfect test, p=0.30, moe=0.05):
#   p_app = 0.30, correction = 1
#   n_cont = 1.959964^2 * 0.30 * 0.70 / 0.05^2 = 322.68 -> ceiling -> 323
#
# DP-2  (se=0.90, sp=0.95, p=0.30, moe=0.05):
#   p_app = 0.30*0.90 + 0.70*0.05 = 0.305, correction = 0.85
#   n_cont = 1.959964^2 * 0.305 * 0.695 / (0.0025 * 0.7225) = 450.82 -> ceiling -> 451
#
# DP-C-1  (n_per_site=10, icc=0.05, p=0.30, moe=0.05):
#   deff = 1 + (10-1)*0.05 = 1.45
#   n_cont = 322.68 * 1.45 = 467.89 -> ceiling -> 468
#   n_sites = ceiling(468 / 10) = 47
#
# DP-C-2  (n_sites=50, icc=0.05, p=0.30, moe=0.05):
#   denom = 50 - 322.68*0.05 = 33.87
#   n_cont = 322.68 * 50 * 0.95 / 33.87 = 452.6 -> ceiling -> 453
#   n_per_site = ceiling(453 / 50) = 10
#
# DP-C-3  (n_sites=10, icc=0.05, p=0.30, moe=0.05):
#   denom = 10 - 322.68*0.05 = -6.13 < 0 -> infeasible
#   min_moe = z * sqrt(0.3*0.7*0.05 / (10*1)) = 1.959964*sqrt(0.00105) ~= 6.35%
#
# DP-F-1  (FPC, fpc_N=500, SRS):
#   n_cont = 322.68 (SRS base)
#   n_adj  = 322.68 * 500 / (322.68 + 500 - 1) = 161340 / 821.68 = 196.35 -> ceiling -> 197

test_that("DP-1: perfect test, SRS returns correct n and deff", {
  res <- design_precision(prevalence = 0.3, moe = 0.05)

  expect_equal(res$n,            323)
  expect_equal(res$deff,         1)
  expect_equal(res$apparent_prev, 0.3, tolerance = 1e-6)
  expect_null(res$n_sites)
  expect_null(res$n_per_site)
})

test_that("DP-2: imperfect test inflates n via Rogan-Gladen variance", {
  res <- design_precision(0.3, 0.05, sensitivity = 0.9, specificity = 0.95)

  # apparent prevalence = 0.30*0.90 + 0.70*0.05 = 0.305
  expect_equal(res$apparent_prev, 0.305, tolerance = 1e-6)
  expect_equal(res$n, 451)
  expect_gt(res$n, 323)   # imperfect test always needs more samples
})

test_that("DP-3: prevalence=0.5 (max variance) gives largest n", {
  # n = ceiling(1.959964^2 * 0.25 / 0.05^2) = ceiling(384.15) = 385
  r50 <- design_precision(0.5, 0.05)
  r30 <- design_precision(0.3, 0.05)
  expect_equal(r50$n, 385)
  expect_gt(r50$n, r30$n)
})

test_that("DP-4: moe=0.499 (nearly maximum) -> tiny n", {
  # Extremely wide CI, so n is very small:
  # n = ceiling(1.959964^2 * 0.3 * 0.7 / 0.499^2) = ceiling(3.24) = 4
  res <- design_precision(0.3, 0.499)
  expect_equal(res$n, 4)
})

test_that("DP-5: conf_level=0.99 requires more samples than 0.95", {
  r99 <- design_precision(0.3, 0.05, conf_level = 0.99)
  r95 <- design_precision(0.3, 0.05)
  expect_gt(r99$n, r95$n)
  # n scales with z^2: ratio ~= (2.576/1.960)^2 = 1.727
  expect_equal(r99$n / r95$n, (qnorm(0.995) / qnorm(0.975))^2,
               tolerance = 0.02)
  # n = ceiling(qnorm(0.995)^2 * 0.3 * 0.7 / 0.05^2) = ceiling(557.33) = 558
  expect_equal(r99$n, 558)
})

test_that("DP-6: near-perfect se=sp=0.999 barely inflates n above perfect test", {
  # p_app = 0.3*0.999 + 0.7*0.001 = 0.3004, correction = 0.998
  # n = ceiling(1.959964^2 * 0.3004 * 0.6996 / (0.05^2 * 0.998^2))
  #   = ceiling(324.22) = 325
  r_perfect <- design_precision(0.3, 0.05)
  r_near    <- design_precision(0.3, 0.05, sensitivity = 0.999, specificity = 0.999)
  expect_equal(r_near$n, 325)
  expect_lt(r_near$n, r_perfect$n + 5)  # inflation should be tiny
})

test_that("DP-7: very low prevalence (p=0.001) returns valid n and apparent_prev", {
  # n = ceiling(1.959964^2 * 0.001 * 0.999 / 0.001^2) = ceiling(3837.62) = 3838
  res <- design_precision(0.001, 0.001)
  expect_equal(res$n, 3838)
  expect_equal(res$apparent_prev, 0.001, tolerance = 1e-6)  # perfect test
})

test_that("DP-8: conf_level=0.5 produces a small n (low confidence threshold)", {
  # z for 50% confidence = qnorm(0.75) = 0.6745
  # n = ceiling(0.6745^2 * 0.3 * 0.7 / 0.05^2) = ceiling(38.21) = 39
  res <- design_precision(0.3, 0.05, conf_level = 0.50)
  expect_equal(res$n, 39)
})

test_that("DP-9: wide moe=0.3 gives very small n", {
  # n = ceiling(1.959964^2 * 0.3 * 0.7 / 0.3^2) = ceiling(8.96) = 9
  res <- design_precision(0.3, 0.3)
  expect_equal(res$n, 9)
})

test_that("DP-10: prevalence=0.9999 (near 1) gives the smallest possible n", {
  # n = ceiling(1.959964^2 * 0.9999 * 0.0001 / 0.05^2) = ceiling(0.15) = 1
  res <- design_precision(0.9999, 0.05)
  expect_equal(res$n, 1)
  # sensitivity and specificity are 1 (the defaults): the test never misses
  # a case and never gives a false positive, so the prevalence it would
  # measure (apparent_prev) is the same as the true prevalence.
  expect_equal(res$apparent_prev, 0.9999, tolerance = 1e-6)
})

test_that("DP-C-1: fixed n_per_site inflates n by deff = 1.45", {
  res <- design_precision(0.3, 0.05, n_per_site = 10, icc = 0.05)

  expect_equal(res$deff, 1.45, tolerance = 1e-6)
  expect_equal(res$n, 468)
  expect_gt(res$n, 323)

  # sites needed = ceiling(468 / 10) = 47
  expect_equal(res$n_sites,    47)
  expect_equal(res$n_per_site, 10)
})

test_that("DP-C-2: fixed n_sites resolves circularity via closed-form", {
  res <- design_precision(0.3, 0.05, n_sites = 50, icc = 0.05)

  expect_equal(res$n,         453)
  expect_equal(res$n_sites,    50)
  expect_equal(res$n_per_site, 10)   # ceiling(453/50) = 10
})

test_that("DP-C-3: infeasible n_sites produces informative error", {
  # denom = 10 - 322.68*0.05 < 0 -> unachievable
  expect_error(
    design_precision(0.3, 0.05, n_sites = 10, icc = 0.05),
    "unachievable"
  )
})

test_that("DP-C-4: the reported minimum achievable MOE matches the min_moe formula", {
  # DP-C-3 only checks that the error fires. This checks the minimum-MOE
  # number the message reports, for both a perfect and an imperfect test.
  # With a perfect test p_app = p and correction = 1, so a mistake in
  # either term would not change the number; the imperfect case covers that.
  z <- qnorm(0.975)

  # Perfect test: p_app = 0.3, correction = 1 -> 6.4%
  min_moe <- z * sqrt(0.3 * 0.7 * 0.05 / (10 * 1^2))
  expect_error(
    design_precision(0.3, 0.05, n_sites = 10, icc = 0.05),
    sprintf("%.1f%%", 100 * min_moe),
    fixed = TRUE
  )

  # Imperfect test: p_app = 0.305, correction = 0.85 -> 7.5%
  p_app   <- 0.3 * 0.9 + 0.7 * 0.05
  min_moe <- z * sqrt(p_app * (1 - p_app) * 0.05 / (10 * 0.85^2))
  expect_error(
    design_precision(0.3, 0.05, sensitivity = 0.9, specificity = 0.95,
                     n_sites = 10, icc = 0.05),
    sprintf("%.1f%%", 100 * min_moe),
    fixed = TRUE
  )
})

test_that("DP-C-5: n_per_site=1 with icc>0 -> deff=1 (no clustering when cluster size=1)", {
  # deff = 1 + (1-1)*icc = 1 regardless of icc
  res <- design_precision(0.3, 0.05, n_per_site = 1, icc = 0.05)
  expect_equal(res$deff, 1, tolerance = 1e-10)
  expect_equal(res$n, design_precision(0.3, 0.05)$n)  # same as SRS
})

test_that("DP-C-6: n_sites=1 with icc=0 -> SRS branch returns n=323, n_sites=1", {
  # icc=0 triggers SRS branch; distribution step still assigns n_sites=1
  res <- design_precision(0.3, 0.05, n_sites = 1, icc = 0)
  expect_equal(res$n,        323)
  expect_equal(res$deff,     1,   tolerance = 1e-10)
  expect_equal(res$n_sites,  1)
})

test_that("DP-C-7: icc=0 with n_sites=50 -> SRS n distributed across 50 sites", {
  # icc=0 dominates: deff=1, n=323; but n_sites=50 is used in distribution
  res <- design_precision(0.3, 0.05, n_sites = 50, icc = 0)
  expect_equal(res$n,         323)
  expect_equal(res$deff,      1,   tolerance = 1e-10)
  expect_equal(res$n_sites,   50)
  expect_equal(res$n_per_site, ceiling(323 / 50))
})

test_that("DP-C-8: n_per_site=10L (integer type) works identically to double", {
  # Uses 10, not 1: with n_per_site = 1 the design effect is 1 whatever
  # happens, so the test could not tell integer and double apart.
  r_int <- design_precision(0.3, 0.05, n_per_site = 10L, icc = 0.05)
  r_dbl <- design_precision(0.3, 0.05, n_per_site = 10,  icc = 0.05)
  expect_equal(r_int$n,    r_dbl$n)
  expect_equal(r_int$deff, r_dbl$deff)
  expect_equal(r_int$n,    468)   # same as DP-C-1
  expect_equal(r_int$deff, 1.45)
})

test_that("DP-C-9: n_per_site given but icc=0 -> deff=1, same n as SRS", {
  # deff = 1 + (n_per_site - 1)*0 = 1 regardless of cluster size
  res <- design_precision(0.3, 0.05, n_per_site = 50, icc = 0)
  expect_equal(res$deff, 1, tolerance = 1e-10)
  expect_equal(res$n, design_precision(0.3, 0.05)$n)
})

test_that("DP-C-10: large n_per_site with high icc -> very large n", {
  # deff = 1 + (10000 - 1) * 0.3 = 3000.7
  # n = ceiling(322.68 * 3000.7) = ceiling(968273.5) = 968274
  res <- design_precision(0.3, 0.05, n_per_site = 10000, icc = 0.3)
  expect_equal(res$deff, 3000.7)
  expect_equal(res$n, 968274)
})

test_that("DP-C-11: more sites than the study needs people is rejected; the same sites work for a bigger study", {
  # Whether a number of sites is too many depends on how many people the
  # study needs. Every site must get at least 1 person.
  #
  # prevalence 0.3, moe 0.05: the study needs only 322.68 people (before
  # clustering). Spread over 323 sites, that is 322.68 / 323 = 0.999 people
  # per site -- less than 1 -- so it is rejected.
  expect_error(design_precision(0.3, 0.05, n_sites = 323, icc = 0.05),
               "n_sites = 323 is >= the SRS sample size (n_base ~= 323)",
               fixed = TRUE)
  # prevalence 0.3, moe 0.02: the study needs 2016.77 people (before
  # clustering), so 323 sites is fine. With clustering: n = 2786,
  # 9 people per site, deff = 1 + (9 - 1) * 0.05 = 1.4.
  res <- design_precision(0.3, 0.02, n_sites = 323, icc = 0.05)
  expect_equal(res$n, 2786)
  expect_equal(res$n_per_site, 9)
  expect_equal(res$deff, 1.4)
})

test_that("DP-C-12: the reported deff matches the reported people per site", {
  # With n_sites fixed, the function rounds people per site up to a whole
  # number, then works deff out again from that rounded number, so that
  # deff = 1 + (n_per_site - 1) * icc describes the design you would
  # actually carry out.

  # 50 sites: the function needs n = 452.6 -> 453 people in total.
  # 453 / 50 = 9.06 per site, rounded up to 10, so deff is worked out
  # again from 10: 1 + (10 - 1) * 0.05 = 1.45.
  # (The plan then collects 50 * 10 = 500, more than n = 453 -- see the
  # note on n in ?design_precision.)
  a <- design_precision(0.3, 0.05, n_sites = 50, icc = 0.05)
  expect_equal(a$n, 453)
  expect_equal(a$n_per_site, 10)
  expect_equal(a$deff, 1.45)

  # 50 sites in a population of 1000: the FPC lowers n to 312.
  # 312 / 50 = 6.24 per site, rounded up to 7: deff = 1 + (7 - 1) * 0.05 = 1.3
  # (collects 50 * 7 = 350)
  b <-design_precision(0.3, 0.05, n_sites = 50, icc = 0.05, fpc_N = 1000)
  expect_equal(b$n, 312)
  expect_equal(b$n_per_site, 7)
  expect_equal(b$deff, 1.3)
})

test_that("DP-C-13: a tiny icc (1e-12) is treated as exactly 0, with or without sites", {
  # Any icc below about 0.000000015 counts as 0. A value like 1e-12 is
  # almost always leftover rounding from another calculation, not real
  # clustering. Without this rule, icc = 1e-12 with no sites would be
  # rejected like icc = 0.05 with no sites is (DP-V-6).

  # No sites: accepted, and gives the same n as no clustering (323)
  res <- design_precision(0.3, 0.05, icc = 1e-12)
  expect_equal(res$n, 323)
  expect_equal(res$deff, 1)

  # 50 sites: n = 323 split over 50 sites (7 each), deff exactly 1
  res <- design_precision(0.3, 0.05, n_sites = 50, icc = 1e-12)
  expect_equal(res$n, 323)
  expect_equal(res$n_per_site, 7)
  expect_equal(res$deff, 1)
})

# FPC matters once the sampling fraction f = n/N is large (Cochran 1977's
# standard threshold: f > ~5-10%). fpc_N=500 here gives f = 197/500 ~= 39%
# of the population sampled, well into the range where FPC matters.
test_that("DP-F-1: FPC reduces required n for a small population", {
  res_fpc  <- design_precision(0.3, 0.05, fpc_N = 500)
  res_nofpc <- design_precision(0.3, 0.05)

  # FPC-adjusted n should be smaller
  expect_lt(res_fpc$n, res_nofpc$n)
  # Hand-checked: n_adj = 322.68*500/(322.68+500-1) = 196.35 -> 197
  expect_equal(res_fpc$n, 197)
  expect_equal(res_fpc$fpc_N, 500)
})

test_that("DP-F-2: very large fpc_N has negligible effect on n", {
  r_fpc <- design_precision(0.3, 0.05, fpc_N = 1e8)
  r_srs <- design_precision(0.3, 0.05)
  # f = 323/1e8 ~= 0.0003% of the population sampled -> FPC negligible
  # FPC factor = N/(n + N - 1) ~= 1 - n/N ~= 0.999997 -> barely changes n
  # Before rounding up: 322.6815 with FPC vs 322.6825 without -- both -> 323.
  expect_equal(r_fpc$n, r_srs$n)
  expect_equal(r_fpc$n, 323)
})

test_that("DP-F-3: clustering + FPC combined: n_per_site + fpc_N both reduce final n", {
  # Clustered alone
  r_cluster <- design_precision(0.3, 0.05, n_per_site = 10, icc = 0.05)
  # Clustered + FPC (population = 1000)
  r_both    <- design_precision(0.3, 0.05, n_per_site = 10, icc = 0.05, fpc_N = 1000)
  # FPC should reduce n below the clustering-only case
  expect_lt(r_both$n, r_cluster$n)
  # n = ceiling(467.89 * 1000 / (467.89 + 1000 - 1)) = ceiling(318.97) = 319
  expect_equal(r_both$n, 319)
  expect_equal(r_both$n_sites, 32)   # ceiling(319 / 10)
  expect_equal(r_both$fpc_N, 1000)
})

test_that("DP-F-4: a population of 1 or 2 -> n is the whole population", {
  # FPC: n = n0 * N / (n0 + N - 1), with n0 = 322.68 (the SRS n before rounding)
  # N = 1: n = 322.68 / 322.68 = 1
  # N = 2: n = 645.36 / 323.68 = 1.99 -> rounds up to 2
  expect_equal(design_precision(0.3, 0.05, fpc_N = 1)$n, 1)
  expect_equal(design_precision(0.3, 0.05, fpc_N = 2)$n, 2)
})

test_that("DP-V-1: prevalence = 0 and prevalence = 1 are rejected as degenerate", {
  # Both are valid probabilities, but there is nothing uncertain to estimate.
  expect_error(design_precision(0, 0.05), "no uncertainty to estimate")
  expect_error(design_precision(1, 0.05), "no uncertainty to estimate")
})

test_that("DP-V-2: prevalence below 0 or above 1 is rejected as an invalid probability", {
  # prevalence has two different error messages: one for 0 or 1 (DP-V-1:
  # a valid probability, but nothing to estimate) and one for values
  # outside 0-1 (not a probability at all). These should get the second.
  expect_error(design_precision(-0.01, 0.05),
               "`prevalence` is a fraction of the population and cannot be negative",
               fixed = TRUE)
  expect_error(design_precision(1.5, 0.05),
               "`prevalence` is a fraction of the population and cannot be greater than 1",
               fixed = TRUE)
})

test_that("DP-V-3: moe must be > 0 -- 0 and negative are rejected, a small positive value works", {
  # 0 and negative get different explanations in the message.
  expect_error(design_precision(0.3, 0),
               "a margin of error of 0 would demand infinite precision", fixed = TRUE)
  expect_error(design_precision(0.3, -0.05),
               "it is a width and cannot be negative", fixed = TRUE)
  # Just above the boundary is accepted and gives the hand-checked n:
  # n = ceiling(1.959964^2 * 0.3 * 0.7 / 0.001^2) = 806707
  expect_equal(design_precision(0.3, 0.001)$n, 806707)
})

test_that("DP-V-4: sensitivity + specificity <= 1 is rejected (Rogan-Gladen correction)", {
  # The correction divides by (sensitivity + specificity - 1), so that
  # value must be above 0.
  # Sum below 1 (0.2 + 0.2 = 0.4): it would divide by a negative number.
  expect_error(design_precision(0.3, 0.05, sensitivity = 0.2, specificity = 0.2),
               "must exceed 1")
  # Sum exactly 1 (0.5 + 0.5): it would divide by 0.
  expect_error(design_precision(0.3, 0.05, sensitivity = 0.5, specificity = 0.5),
               "must exceed 1")
})

test_that("DP-V-5: sensitivity + specificity above 1 but below 1.1 warns, 1.1 and above does not", {
  # Sums between 1 and 1.1 are allowed (no error), but n gets very large
  # (20,760 at 1.05 vs 323 for a perfect test), so the function warns.
  expect_warning(
    design_precision(0.3, 0.05, sensitivity = 1, specificity = 0.05),   # 1.05
    "which is very close to 1", fixed = TRUE
  )
  # Either side of the 1.1 cutoff:
  expect_warning(
    design_precision(0.3, 0.05, sensitivity = 1, specificity = 0.09),   # 1.09
    "which is very close to 1", fixed = TRUE
  )
  expect_no_warning(design_precision(0.3, 0.05, sensitivity = 1, specificity = 0.10))  # 1.10
  expect_no_warning(design_precision(0.3, 0.05, sensitivity = 1, specificity = 0.11))  # 1.11
})

test_that("DP-V-6: icc > 0 with no cluster structure (no n_sites or n_per_site) is rejected", {
  expect_error(design_precision(0.3, 0.05, icc = 0.05),
               "there is no cluster structure", fixed = TRUE)
  # The reverse is allowed: n_sites or n_per_site without icc uses the
  # default icc = 0 (no clustering effect), so n is the SRS n of 323.
  expect_equal(design_precision(0.3, 0.05, n_sites = 50)$n,    323)
  expect_equal(design_precision(0.3, 0.05, n_per_site = 10)$n, 323)
})

test_that("DP-V-7: supplying both n_sites and n_per_site is rejected", {
  expect_error(design_precision(0.3, 0.05, n_sites = 50, n_per_site = 10, icc = 0.05),
               "Supply at most one of `n_sites` or `n_per_site`, not both", fixed = TRUE)
})

test_that("DP-V-8: conf_level must be in (0, 1) -- boundaries and outside values rejected, inside works", {
  for (cl in c(0, 1, -0.05, 1.05)) {
    expect_error(design_precision(0.3, 0.05, conf_level = cl),
                 "`conf_level` must be in (0, 1)", fixed = TRUE)
  }
  # A valid value inside the range gives the hand-checked n:
  # n = ceiling(qnorm(0.95)^2 * 0.3 * 0.7 / 0.05^2) = ceiling(227.27) = 228
  expect_equal(design_precision(0.3, 0.05, conf_level = 0.90)$n, 228)
})

test_that("DP-V-9: fpc_N must be a positive whole number", {
  expect_error(design_precision(0.3, 0.05, fpc_N =  0),
               "`fpc_N` must be a single finite positive integer", fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, fpc_N = -50),
               "`fpc_N` must be a single finite positive integer", fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, fpc_N = 500.5),
               "`fpc_N` must be a single finite positive integer", fixed = TRUE)
  # A valid value gives the hand-checked n (see DP-F-1 in the header): 197
  expect_equal(design_precision(0.3, 0.05, fpc_N = 500)$n, 197)
})

test_that("DP-V-10: n_per_site must be a positive whole number", {
  expect_error(design_precision(0.3, 0.05, n_per_site = 0,    icc = 0.05),
               "`n_per_site` must be a single finite positive integer", fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, n_per_site = -5,   icc = 0.05),
               "`n_per_site` must be a single finite positive integer", fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, n_per_site = 10.7, icc = 0.05),
               "`n_per_site` must be a single finite positive integer", fixed = TRUE)
  # A valid value gives the hand-checked n (see DP-C-1 in the header): 468
  expect_equal(design_precision(0.3, 0.05, n_per_site = 10, icc = 0.05)$n, 468)
})

test_that("DP-V-11: n_sites must be a positive whole number", {
  expect_error(design_precision(0.3, 0.05, n_sites = 0,    icc = 0.05),
               "`n_sites` must be a single finite positive integer", fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, n_sites = -5,   icc = 0.05),
               "`n_sites` must be a single finite positive integer", fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, n_sites = 50.5, icc = 0.05),
               "`n_sites` must be a single finite positive integer", fixed = TRUE)
  # A valid value gives the hand-checked n (see DP-C-2 in the header): 453
  expect_equal(design_precision(0.3, 0.05, n_sites = 50, icc = 0.05)$n, 453)
})

test_that("DP-V-12: specificity's own range check (0, 1] gives the right message for each bad value", {
  # These can't use DP-V-4's values (0.2 + 0.2): each of those is fine on
  # its own, so they are stopped by the combined "must exceed 1" check
  # instead of by specificity's own check.
  # specificity = 0 with the default sensitivity = 1 is the "1 + 0 = 1"
  # case: it is stopped here, by specificity's own check, before the
  # combined check runs. The mirror case (sensitivity = 0, specificity = 1)
  # is in DP-V-13.
  expect_error(design_precision(0.3, 0.05, specificity = 0),
               "A specificity of 0 means", fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, specificity = 1.001),
               "`specificity` is a diagnostic probability and cannot be greater than 1",
               fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, specificity = -0.5),
               "`specificity` is a diagnostic probability and cannot be negative",
               fixed = TRUE)
})

test_that("DP-V-13: sensitivity's own range check (0, 1] gives the right message for each bad value", {
  # Mirror of the specificity test above. sensitivity = 0 with the default
  # specificity = 1 is the "0 + 1 = 1" case: it is stopped here, by
  # sensitivity's own check, before the combined check runs.
  expect_error(design_precision(0.3, 0.05, sensitivity = 0),
               "A sensitivity of 0 means", fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, sensitivity = 1.001),
               "`sensitivity` is a diagnostic probability and cannot be greater than 1",
               fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, sensitivity = -0.5),
               "`sensitivity` is a diagnostic probability and cannot be negative",
               fixed = TRUE)
})

test_that("DP-V-14: NaN, Inf and -Inf in any single-number argument give the finite-value error", {
  # Note: a bare NA is a logical value in R, so it is caught earlier by the
  # type check ("got class `logical`"), not by this check. NaN and Inf are
  # numbers, so they get through the type check and reach this one.
  expect_error(design_precision(0.3, NaN),
               "`moe` must be a single finite number (got NaN)", fixed = TRUE)
  expect_error(design_precision(0.3, Inf),
               "`moe` must be a single finite number (got Inf)", fixed = TRUE)
  expect_error(design_precision(0.3, -Inf),
               "`moe` must be a single finite number (got -Inf)", fixed = TRUE)
  expect_error(design_precision(NaN, 0.05),
               "`prevalence` must be a single finite number (got NaN)", fixed = TRUE)
  expect_error(design_precision(Inf, 0.05),
               "`prevalence` must be a single finite number (got Inf)", fixed = TRUE)
  expect_error(design_precision(-Inf, 0.05),
               "`prevalence` must be a single finite number (got -Inf)", fixed = TRUE)
  # The other four single-number arguments each have their own check line
  expect_error(design_precision(0.3, 0.05, sensitivity = NaN),
               "`sensitivity` must be a single finite number (got NaN)", fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, specificity = Inf),
               "`specificity` must be a single finite number (got Inf)", fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, conf_level = NaN),
               "`conf_level` must be a single finite number (got NaN)", fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, icc = Inf),
               "`icc` must be a single finite number (got Inf)", fixed = TRUE)
})

test_that("DP-V-15: n_sites as a string, vector, TRUE or Inf is rejected, naming the problem", {
  # A string or TRUE is reported by its class; a vector by its length.
  # TRUE must not be quietly read as n_sites = 1, and Inf must not slip
  # through as a "whole number". A valid n_sites is checked in DP-V-11.
  expect_error(design_precision(0.3, 0.05, n_sites = "50", icc = 0.05),
               "`n_sites` must be a single finite positive integer (got class `character`)",
               fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, n_sites = c(10, 20), icc = 0.05),
               "`n_sites` must be a single finite positive integer (got length = 2)",
               fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, n_sites = TRUE, icc = 0.05),
               "`n_sites` must be a single finite positive integer (got class `logical`)",
               fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, n_sites = Inf, icc = 0.05),
               "`n_sites` must be a single finite positive integer (got Inf)",
               fixed = TRUE)
})

test_that("DP-V-16: n_per_site as a string, vector, TRUE or Inf is rejected, naming the problem", {
  # Same four cases as DP-V-15. A valid n_per_site is checked in DP-V-10.
  expect_error(design_precision(0.3, 0.05, n_per_site = "10", icc = 0.05),
               "`n_per_site` must be a single finite positive integer (got class `character`)",
               fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, n_per_site = c(5, 10), icc = 0.05),
               "`n_per_site` must be a single finite positive integer (got length = 2)",
               fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, n_per_site = TRUE, icc = 0.05),
               "`n_per_site` must be a single finite positive integer (got class `logical`)",
               fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, n_per_site = Inf, icc = 0.05),
               "`n_per_site` must be a single finite positive integer (got Inf)",
               fixed = TRUE)
})

test_that("DP-V-17: fpc_N as a string, vector, TRUE or Inf is rejected, naming the problem", {
  # Same four cases as DP-V-15. A valid fpc_N is checked in DP-V-9.
  expect_error(design_precision(0.3, 0.05, fpc_N = "500"),
               "`fpc_N` must be a single finite positive integer (got class `character`)",
               fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, fpc_N = c(500, 600)),
               "`fpc_N` must be a single finite positive integer (got length = 2)",
               fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, fpc_N = TRUE),
               "`fpc_N` must be a single finite positive integer (got class `logical`)",
               fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, fpc_N = Inf),
               "`fpc_N` must be a single finite positive integer (got Inf)",
               fixed = TRUE)
})

test_that("DP-V-18: icc must be in [0, 1] -- below 0 and above 1 are rejected, 1 works", {
  expect_error(design_precision(0.3, 0.05, icc = -0.1),
               "`icc` must be in [0, 1] (got -0.1). `icc` is a correlation and cannot be negative.",
               fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, icc = 1.5),
               "`icc` must be in [0, 1] (got 1.5). `icc` is a correlation and cannot be greater than 1.",
               fixed = TRUE)
  # icc = 1 itself is allowed: deff = 1 + (2 - 1) * 1 = 2, so
  # n = ceiling(322.68 * 2) = ceiling(645.37) = 646
  res <- design_precision(0.3, 0.05, n_per_site = 2, icc = 1)
  expect_equal(res$deff, 2)
  expect_equal(res$n, 646)
})

test_that("DP-V-19: moe must be < 0.5 -- 0.5 and 0.6 are rejected, 0.499 works", {
  expect_error(design_precision(0.3, 0.5),
               "`moe` must be less than 0.5 (got 0.5)", fixed = TRUE)
  expect_error(design_precision(0.3, 0.6),
               "`moe` must be less than 0.5 (got 0.6)", fixed = TRUE)
  # Just below the boundary is accepted (same case as DP-4): n = 4
  expect_equal(design_precision(0.3, 0.499)$n, 4)
})

test_that("DP-V-20: a vector in a single-number argument is rejected, giving its length", {
  expect_error(design_precision(c(0.2, 0.3), 0.05),
               "`prevalence` must be a single number in (0, 1) (got class `numeric`, length 2)",
               fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, n_per_site = 10, icc = c(0.05, 0.1)),
               "`icc` must be a single number in [0, 1] (got class `numeric`, length 2)",
               fixed = TRUE)
})

test_that("DP-V-21: text, TRUE/FALSE, a bare NA or a list is rejected, naming the class", {
  # None of these is a number, so each is stopped by the type check and the
  # message says what was passed. A bare NA counts as TRUE/FALSE-type
  # (logical) in R, so it gets the same "class `logical`" message.
  expect_error(design_precision("0.3", 0.05),
               "`prevalence` must be a single number in (0, 1) (got class `character`, length 1)",
               fixed = TRUE)
  expect_error(design_precision(0.3, "0.05"),
               "`moe` must be a single number in (0, 0.5) (got class `character`, length 1)",
               fixed = TRUE)
  # TRUE must not be quietly read as a perfect test (1); the message says so.
  expect_error(design_precision(0.3, 0.05, sensitivity = TRUE),
               "`sensitivity` must be a single number in (0, 1] (got class `logical`, length 1). Note: `TRUE`/`FALSE` is logical, not numeric -- pass 1 for a perfect test.",
               fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, specificity = FALSE),
               "`specificity` must be a single number in (0, 1] (got class `logical`, length 1)",
               fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, specificity = NA),
               "`specificity` must be a single number in (0, 1] (got class `logical`, length 1)",
               fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, icc = FALSE),
               "`icc` must be a single number in [0, 1] (got class `logical`, length 1)",
               fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, conf_level = NA),
               "`conf_level` must be a single number in (0, 1) (got class `logical`, length 1)",
               fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, conf_level = list(0.95)),
               "`conf_level` must be a single number in (0, 1) (got class `list`, length 1)",
               fixed = TRUE)
})

test_that("DP-V-22: more sites, or more people per site, than the whole population is rejected", {
  expect_error(design_precision(0.3, 0.05, n_sites = 100, fpc_N = 50, icc = 0.05),
               "`n_sites` (100) exceeds `fpc_N` (50): you cannot have more clusters than individuals in the population.",
               fixed = TRUE)
  expect_error(design_precision(0.3, 0.05, n_per_site = 100, fpc_N = 50, icc = 0.01),
               "`n_per_site` (100) exceeds `fpc_N` (50): a cluster cannot be larger than the whole population.",
               fixed = TRUE)
  # 20 sites in a population of 5000 is fine: n = 1205, 61 per site,
  # deff = 1 + (61 - 1) * 0.05 = 4
  res <- design_precision(0.3, 0.05, n_sites = 20, fpc_N = 5000, icc = 0.05)
  expect_equal(res$n, 1205)
  expect_equal(res$n_per_site, 61)
  expect_equal(res$deff, 4)
})

test_that("DP-R-1: return list contains all expected fields", {
  res <- design_precision(0.2, 0.05)
  expect_named(res, c("n", "n_eff", "n_sites", "n_per_site", "prevalence",
                       "apparent_prev", "moe", "conf_level", "sensitivity",
                       "specificity", "icc", "deff", "fpc_N"),
               ignore.order = FALSE)
  # prevalence is passed back exactly as given
  expect_equal(res$prevalence, 0.2)
})

test_that("DP-R-2: n_eff is the equivalent simple-random-sample size (323), even with FPC or clustering", {
  # n_eff answers: "how many people picked completely at random would give
  # this precision?" Here that is always 323. Only n, the number you
  # actually collect, changes with the design.

  # No FPC, no clustering: collect 323, worth 323
  srs <- design_precision(0.3, 0.05)
  expect_equal(srs$n_eff, srs$n)

  # Population of 500: collect only 197, still worth 323
  fpc <- design_precision(0.3, 0.05, fpc_N = 500)
  expect_lt(fpc$n, srs$n)
  expect_equal(fpc$n_eff, srs$n_eff)

  # 20 people per site, icc 0.05 (deff 1.95): people in a site are alike,
  # so collect 630 to be worth 323
  clus <- design_precision(0.3, 0.05, n_per_site = 20, icc = 0.05)
  expect_equal(clus$n_eff, srs$n_eff)
  expect_gt(clus$n, clus$n_eff)
})
