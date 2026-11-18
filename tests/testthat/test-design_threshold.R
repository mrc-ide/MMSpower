# Reference cases for design_threshold() (exact binomial power)
#
# The expected numbers come from a separate brute-force calculation written
# here (bf_power / bf_first_safe): for every possible count of positives
# 0..n it asks whether binom.test() would reject (one-sided p-values from
# pbinom(); two-sided = twice the smaller one), and adds up dbinom() for
# those counts. It shares no code with the package.
#
#   threshold 0.05, prevalence 0.10, greater, 80 percent:  first 169, safe 191
#   threshold 0.05, prevalence 0.10, 90 percent:           checked in DT-2
#   threshold 0.05, prevalence 0.10, two.sided:            first 213, safe 234
#   threshold 0.10, prevalence 0.05, less (saw-tooth):     first 179, safe 203
#     (179 people give 81.3 percent, 182 to 190 all below 80 percent)

bf_power <- function(n, th, p, alt = "greater", a = 0.05) {
  x  <- 0:n
  up <- pbinom(x - 1, n, th, lower.tail = FALSE)   # P(X >= x)
  dn <- pbinom(x, n, th)                           # P(X <= x)
  pv <- switch(alt, greater = up, less = dn, two.sided = pmin(1, 2 * pmin(up, dn)))
  sum(dbinom(x, n, p)[pv < a])
}
bf_first_safe <- function(th, p, alt = "greater", power = 0.8, a = 0.05, n_max = 2000) {
  pw <- vapply(1:n_max, bf_power, numeric(1), th = th, p = p, alt = alt, a = a)
  c(first = which(pw >= power)[1], safe = max(which(pw < power)) + 1)
}


# ---------------------------------------------------------------------------
# Core scenarios (DT-)
# ---------------------------------------------------------------------------

test_that("DT-1: 5 versus 10 percent, 'greater', 80 percent power gives n = 191 (first 169)", {
  res <- design_threshold(threshold = 0.05, prevalence = 0.10)
  expect_equal(c(first = res$n_first, safe = res$n), bf_first_safe(0.05, 0.10))
  expect_equal(res$n, 191)
  expect_equal(res$n_first, 169)
  expect_equal(res$n_eff, 191)
  expect_equal(res$deff, 1)
  expect_null(res$n_sites)
  expect_null(res$n_per_site)
  # Every n from 191 to 600 reaches 80 percent; 190 does not
  expect_lt(bf_power(190, 0.05, 0.10), 0.8)
  expect_true(all(vapply(191:600, bf_power, numeric(1), th = 0.05, p = 0.10) >= 0.8))
})

test_that("DT-2: more power needs more people; 0.5, 0.9 and 0.999 power match or pass the brute-force search", {
  r50  <- design_threshold(threshold = 0.10, prevalence = 0.20, power = 0.50)
  r80  <- design_threshold(threshold = 0.10, prevalence = 0.20, power = 0.80)
  r90  <- design_threshold(threshold = 0.10, prevalence = 0.20, power = 0.90)
  r99  <- design_threshold(threshold = 0.10, prevalence = 0.20, power = 0.99)
  expect_lt(r50$n, r80$n)
  expect_lt(r80$n, r90$n)
  expect_lt(r90$n, r99$n)
  expect_equal(c(first = r50$n_first, safe = r50$n), bf_first_safe(0.10, 0.20, power = 0.5))
  expect_equal(c(first = r80$n_first, safe = r80$n), bf_first_safe(0.10, 0.20))
  g90 <- design_threshold(threshold = 0.05, prevalence = 0.10, power = 0.90)
  expect_equal(c(first = g90$n_first, safe = g90$n), bf_first_safe(0.05, 0.10, power = 0.90))
  # A very high power gives a large but finite n
  r999 <- design_threshold(threshold = 0.05, prevalence = 0.10, power = 0.999)
  expect_true(is.finite(r999$n))
  expect_gt(r999$n, 2 * 191)
})

test_that("DT-3: a two-sided test needs more people than a one-sided one", {
  r_1s <- design_threshold(threshold = 0.05, prevalence = 0.10, alternative = "greater")
  r_2s <- design_threshold(threshold = 0.05, prevalence = 0.10, alternative = "two.sided")
  expect_gt(r_2s$n, r_1s$n)
  expect_equal(c(r_2s$n_first, r_2s$n), c(213, 234))
  expect_equal(c(first = r_2s$n_first, safe = r_2s$n),
               bf_first_safe(0.05, 0.10, alt = "two.sided"))
  expect_gt(design_threshold(threshold = 0.10, prevalence = 0.20, alternative = "two.sided")$n,
            design_threshold(threshold = 0.10, prevalence = 0.20, alternative = "greater")$n)
})

test_that("DT-4: saw-tooth below the threshold: 'less', 10 versus 5 percent gives n 203, first 179", {
  res <- design_threshold(threshold = 0.10, prevalence = 0.05, alternative = "less")
  expect_equal(res$alternative, "less")
  expect_equal(res$n, 203)
  expect_equal(res$n_first, 179)
  expect_equal(c(first = res$n_first, safe = res$n),
               bf_first_safe(0.10, 0.05, alt = "less"))
  expect_equal(bf_power(179, 0.10, 0.05, "less"), 0.8128, tolerance = 1e-3)
  # n_first far below the normal formula (about 70): the search starts at 1,
  # so it is still found
  r_low <- design_threshold(threshold = 0.067, prevalence = 0.005, alternative = "less")
  expect_equal(c(r_low$n_first, r_low$n), c(44, 69))
  expect_equal(c(first = r_low$n_first, safe = r_low$n),
               bf_first_safe(0.067, 0.005, alt = "less"))
  expect_true(all(vapply(182:190, bf_power, numeric(1), th = 0.10, p = 0.05,
                         alt = "less") < 0.8))
  # Another 'less' design
  r <- design_threshold(threshold = 0.20, prevalence = 0.10, alternative = "less")
  expect_equal(c(first = r$n_first, safe = r$n), bf_first_safe(0.20, 0.10, alt = "less"))
})

test_that("DT-5: conf_level 0.5, 0.90, 0.95 and 0.99: a higher level needs more people, matching the brute-force search", {
  for (cl in c(0.90, 0.99)) {
    r <- design_threshold(threshold = 0.10, prevalence = 0.20, conf_level = cl)
    expect_equal(c(first = r$n_first, safe = r$n), bf_first_safe(0.10, 0.20, a = 1 - cl))
  }
  r50 <- design_threshold(threshold = 0.05, prevalence = 0.10, conf_level = 0.50)
  r95 <- design_threshold(threshold = 0.05, prevalence = 0.10)
  r99 <- design_threshold(threshold = 0.05, prevalence = 0.10, conf_level = 0.99)
  expect_lt(r50$n, r95$n)
  expect_lt(r95$n, r99$n)
})

test_that("DT-6: random designs without sites match the brute-force search", {
  set.seed(5)
  for (i in 1:12) {
    th  <- runif(1, 0.03, 0.6)
    p   <- th + sample(c(-1, 1), 1) * runif(1, 0.08, 0.2)
    if (p <= 0.01 || p >= 0.99) next
    alt <- if (i %% 3 == 0) "two.sided" else if (p > th) "greater" else "less"
    pw  <- sample(c(0.8, 0.9), 1)
    r <- design_threshold(threshold = th, prevalence = p, alternative = alt, power = pw)
    expect_equal(c(first = r$n_first, safe = r$n),
                 bf_first_safe(th, p, alt = alt, power = pw))
  }
})

test_that("DT-7: a prevalence closer to the threshold needs more people", {
  r_far   <- design_threshold(threshold = 0.05, prevalence = 0.15)
  r_mid   <- design_threshold(threshold = 0.05, prevalence = 0.07)
  r_close <- design_threshold(threshold = 0.05, prevalence = 0.06)
  expect_gt(r_mid$n, r_far$n)
  expect_gt(r_close$n, r_mid$n)
  expect_gt(r_close$n, 1000)
  # 10 versus 10.1 percent needs more than 50 times as many as 10 versus 20
  expect_gt(design_threshold(threshold = 0.10, prevalence = 0.101)$n,
            50 * design_threshold(threshold = 0.10, prevalence = 0.20)$n)
})

test_that("DT-8: an imperfect test moves both values to the measured scale and needs more people", {
  # Sensitivity = specificity = 1: measured = true
  res <- design_threshold(threshold = 0.05, prevalence = 0.10)
  expect_equal(res$apparent_prev, 0.10)
  expect_equal(res$threshold_app, 0.05)
  # Measured value = p * Se + (1 - p) * (1 - Sp)
  imp <- design_threshold(threshold = 0.10, prevalence = 0.20,
                          sensitivity = 0.80, specificity = 0.90)
  expect_equal(imp$threshold_app, 0.10 * 0.80 + 0.90 * 0.10)
  expect_equal(imp$apparent_prev, 0.20 * 0.80 + 0.80 * 0.10)
  expect_gt(imp$n, design_threshold(threshold = 0.10, prevalence = 0.20)$n)
  # Sensitivity = specificity = 0.6 shrinks a 10-point gap to 2 points
  expect_gt(design_threshold(threshold = 0.50, prevalence = 0.60,
                             sensitivity = 0.6, specificity = 0.6)$n, 1000)
})

test_that("DT-9: 'less' and 'greater' around a threshold of 0.5 give the same n", {
  r_gt <- design_threshold(threshold = 0.50, prevalence = 0.55, alternative = "greater")
  r_lt <- design_threshold(threshold = 0.50, prevalence = 0.45, alternative = "less")
  expect_equal(r_gt$n, r_lt$n)
})

test_that("DT-10: thresholds close to 0 or 1 give a finite n", {
  lo <- design_threshold(threshold = 0.001, prevalence = 0.01)
  hi <- design_threshold(threshold = 0.999, prevalence = 0.99, alternative = "less")
  expect_true(is.finite(lo$n) && lo$n >= 1)
  expect_true(is.finite(hi$n) && hi$n >= 1)
})

test_that("DT-11: a design needing more than 2 million people is too large to work out exactly", {
  expect_error(design_threshold(threshold = 0.10, prevalence = 0.1001),
               "These settings need more than 2 million people (about 55,700,000), too many to work out exactly. Use a prevalence further from the threshold, or a more accurate test.",
               fixed = TRUE)
})


# ---------------------------------------------------------------------------
# Clustering (DT-C-)
# ---------------------------------------------------------------------------

test_that("DT-C-1: with sites, power matches a separate effective-size calculation", {
  # The test with sites uses effective size m = n / deff and effective
  # positives X * m / n; power is worked out on that scale.
  eff_power <- function(n, d, th, p, a = 0.05) {
    m <- n / d; X <- 0:n; xe <- X * m / n
    pv <- ifelse(xe <= 0, 1, pbeta(th, pmax(xe, 1e-300), m - xe + 1))
    k  <- X[which(pv < a)[1]]
    pbeta(p, k * m / n, m - k * m / n + 1)
  }
  r <- design_threshold(threshold = 0.05, prevalence = 0.10, n_per_site = 20, icc = 0.05)
  expect_equal(r$deff, 1.95)
  expect_equal(r$n_sites, 17)
  expect_equal(r$n, r$n_sites * 20)
  expect_gte(eff_power(r$n, 1.95, 0.05, 0.10), 0.8)
  expect_lt(eff_power(r$n - 20, 1.95, 0.05, 0.10), 0.8)

  r2 <- design_threshold(threshold = 0.05, prevalence = 0.10, n_sites = 10, icc = 0.05)
  expect_equal(r2$n, 10 * r2$n_per_site)
  j <- r2$n_per_site
  expect_gte(eff_power(10 * j, 1 + (j - 1) * 0.05, 0.05, 0.10), 0.8)
  expect_lt(eff_power(10 * (j - 1), 1 + (j - 2) * 0.05, 0.05, 0.10), 0.8)
  expect_equal(c(r2$n_first, r2$n), c(760, 910))
})

test_that("DT-C-2: a fixed number of people per site multiplies up through the design effect", {
  res <- design_threshold(threshold = 0.05, prevalence = 0.10, n_per_site = 10, icc = 0.05)
  expect_equal(res$deff, 1.45)
  expect_equal(res$n_per_site, 10)
  expect_equal(res$n_sites, ceiling(res$n / 10))
  expect_gt(res$n, 191)
  r20 <- design_threshold(threshold = 0.10, prevalence = 0.20, n_per_site = 20, icc = 0.05)
  expect_equal(r20$n_sites, ceiling(r20$n / 20))
  # icc = 1: each site counts as one person, deff = n_per_site
  expect_equal(design_threshold(threshold = 0.05, prevalence = 0.10,
                                n_per_site = 10, icc = 1)$deff, 10)
})

test_that("DT-C-3: a fixed number of sites solves for people per site, and deff matches it", {
  res <- design_threshold(threshold = 0.05, prevalence = 0.10, n_sites = 50, icc = 0.05)
  expect_equal(res$n_sites, 50)
  expect_equal(res$n_per_site, ceiling(res$n / 50))
  expect_equal(res$deff, 1 + (res$n_per_site - 1) * 0.05)
  expect_gt(res$deff, 1)
  r10 <- design_threshold(threshold = 0.10, prevalence = 0.20, n_sites = 10, icc = 0.05)
  expect_equal(r10$n_per_site, ceiling(r10$n / 10))
})

test_that("DT-C-4: one person per site gives deff = 1 and the same n as no sites", {
  for (icc in c(0.05, 0.5)) {
    r <- design_threshold(threshold = 0.05, prevalence = 0.10, n_per_site = 1, icc = icc)
    expect_equal(r$deff, 1)
    expect_equal(r$n, 191)
  }
})

test_that("DT-C-5: too few sites for the target power is rejected; the exact search decides", {
  # The effective size can never pass n_sites / icc; when even that is below
  # the target, no number per site is enough
  expect_error(design_threshold(threshold = 0.05, prevalence = 0.10, n_sites = 5, icc = 0.5),
               "Target power of 80 percent cannot be reached with 5 sites and ICC = 0.5. Increase `n_sites`, lower the ICC, or lower the target power.",
               fixed = TRUE)
  expect_error(design_threshold(threshold = 0.05, prevalence = 0.10, n_sites = 1, icc = 0.5),
               "Target power of 80 percent cannot be reached with 1 site and ICC = 0.5.",
               fixed = TRUE)
  # 8 sites, icc 0.05: the effective size can never pass 8 / 0.05 = 160, but
  # 191 people are needed without sites; the exact search finds no answer
  expect_error(design_threshold(threshold = 0.05, prevalence = 0.10, n_sites = 8, icc = 0.05),
               "Target power of 80 percent cannot be reached with 8 sites and ICC = 0.05. Increase `n_sites`, lower the ICC, or lower the target power.",
               fixed = TRUE)
  # 2 sites, icc 0.05, 'less', 10 versus 1 percent: the normal formula needs
  # an effective size of 41.1, above the limit 2 / 0.05 = 40, but the exact
  # test reaches 80 percent from 110 per site (79 per site is the first;
  # checked by brute force over 1 to 20000 per site)
  r2 <- design_threshold(threshold = 0.10, prevalence = 0.01, alternative = "less",
                         n_sites = 2, icc = 0.05)
  expect_equal(c(r2$n_first, r2$n, r2$n_per_site), c(158, 220, 110))
  expect_gte(MMSpower:::.thr_power(220, r2$deff, NULL, 0.10, 0.01, "less", 0.05), 0.8)
  # Enough sites: a valid design
  expect_equal(design_threshold(threshold = 0.05, prevalence = 0.10, n_sites = 10, icc = 0.05)$n, 910)
})

test_that("DT-C-6: more sites than people needed gives a warning and uses fewer sites, 1 person at each", {
  # 191 people are needed with no sites. 1000 sites, with or without an ICC:
  # 191 sites with 1 person at each (no clustering with 1 per site)
  msg <- "With this prevalence, threshold and power, the study needs 191 people, fewer than the 1000 sites given (`n_sites`). The design uses 191 sites with 1 person at each."
  for (icc in c(0, 0.05)) {
    expect_warning(r <- design_threshold(threshold = 0.05, prevalence = 0.10, n_sites = 1000, icc = icc),
                   msg, fixed = TRUE)
    expect_equal(c(r$n, r$n_first, r$n_sites, r$n_per_site, r$deff), c(191, 169, 191, 1, 1))
  }
  # One site more than needed already warns
  expect_warning(design_threshold(threshold = 0.05, prevalence = 0.10, n_sites = 192),
                 "the study needs 191 people, fewer than the 192 sites given", fixed = TRUE)
  # Exactly as many sites as needed: no warning, 1 per site
  expect_no_warning(r <- design_threshold(threshold = 0.05, prevalence = 0.10, n_sites = 191, icc = 0.05))
  expect_equal(c(r$n, r$n_sites, r$n_per_site), c(191, 191, 1))
  # One site fewer than needed: 2 per site
  r <- design_threshold(threshold = 0.05, prevalence = 0.10, n_sites = 190, icc = 0.05)
  expect_equal(r$n_per_site, 2)
  # 2 people needed: "2 sites"; this used to hang at 2^53 sites
  expect_warning(r <- design_threshold(threshold = 0.05, prevalence = 0.99, n_sites = 2^53, icc = 0.05),
                 "the study needs 2 people, fewer than the 9007199254740992 sites given (`n_sites`). The design uses 2 sites with 1 person at each.",
                 fixed = TRUE)
  expect_equal(r$n, 2)
  # 1 person needed (1 positive in 1 is below 0.05 at a 1 percent
  # threshold, and has a 99 percent chance): "1 person", "1 site"
  expect_warning(design_threshold(threshold = 0.01, prevalence = 0.99, n_sites = 5),
                 "the study needs 1 person, fewer than the 5 sites given (`n_sites`). The design uses 1 site with 1 person at each.",
                 fixed = TRUE)
  # With fpc_N the number needed is the one with the correction: 143 of 500
  # (the same as with no sites), even when n_sites equals the population
  expect_warning(r <- design_threshold(threshold = 0.05, prevalence = 0.10, n_sites = 500, fpc_N = 500),
                 "the study needs 143 people, fewer than the 500 sites given", fixed = TRUE)
  expect_equal(r$n, design_threshold(threshold = 0.05, prevalence = 0.10, fpc_N = 500)$n)
})

test_that("DT-C-6b: one site larger than the number of people needed stops", {
  # 191 people needed: a site of 192 or more is more than enough
  msg <- function(k, need) paste0("`n_per_site` (", k, ") is larger than the number of people needed (", need, "), so one site would be more than enough. Use a smaller `n_per_site`, or set `n_sites = 1` to find how many people one site needs.")
  expect_error(design_threshold(threshold = 0.05, prevalence = 0.10, n_per_site = 500), msg(500, 191), fixed = TRUE)
  expect_error(design_threshold(threshold = 0.05, prevalence = 0.10, n_per_site = 192), msg(192, 191), fixed = TRUE)
  # Exactly the number needed: 1 site
  r <- design_threshold(threshold = 0.05, prevalence = 0.10, n_per_site = 191)
  expect_equal(c(r$n, r$n_sites), c(191, 1))
  # ...the same as n_sites = 1, as the message suggests
  expect_equal(design_threshold(threshold = 0.05, prevalence = 0.10, n_sites = 1)$n, 191)
  # A tiny ICC barely changes the number needed: still too big
  expect_error(design_threshold(threshold = 0.05, prevalence = 0.10, n_per_site = 500, icc = 1e-6),
               msg(500, 191), fixed = TRUE)
  # The number needed includes the design effect: with ICC 0.002 one site
  # of 241 is just enough (power 0.8010, design effect 1.48), so no error,
  # although 241 is more than the 191 needed without sites
  r <- design_threshold(threshold = 0.05, prevalence = 0.10, n_per_site = 241, icc = 0.002)
  expect_equal(c(r$n, r$n_sites, r$deff), c(241, 1, 1.48))
  # With ICC 0.05 one site of 500 is not enough: 9 sites
  expect_equal(design_threshold(threshold = 0.05, prevalence = 0.10, n_per_site = 500, icc = 0.05)$n_sites, 9)
  # With fpc_N the number needed is the one with the correction (131 of 400,
  # the same as with no sites)
  expect_error(design_threshold(threshold = 0.05, prevalence = 0.10, n_per_site = 300, fpc_N = 400),
               msg(300, 131), fixed = TRUE)
  expect_equal(design_threshold(threshold = 0.05, prevalence = 0.10, fpc_N = 400)$n, 131)
})

test_that("DT-C-7: n_eff is the number needed with no sites and no fpc_N", {
  srs <- design_threshold(threshold = 0.05, prevalence = 0.10)
  expect_equal(srs$n_eff, srs$n)
  for (args in list(list(n_per_site = 10, icc = 0.05), list(n_sites = 50, icc = 0.05))) {
    r <- do.call(design_threshold, c(list(threshold = 0.05, prevalence = 0.10), args))
    expect_equal(r$n_eff, 191)
    expect_lte(r$n_eff, r$n)
    expect_equal(r$n_eff, r$n / r$deff, tolerance = 0.1)
  }
  cl <- design_threshold(threshold = 0.10, prevalence = 0.20, n_per_site = 10, icc = 0.05)
  expect_lt(cl$n_eff, cl$n)
  # Only whole sites are checked: with icc 0 and 20 per site, n can come out
  # a little below n_eff (180 people give 80 percent; 190 is not a multiple of 20)
  r0 <- design_threshold(threshold = 0.05, prevalence = 0.10, n_per_site = 20, icc = 0)
  expect_equal(c(r0$n, r0$n_eff), c(180, 191))
  expect_gte(bf_power(180, 0.05, 0.10), 0.8)
})

test_that("DT-C-8: a tiny icc without sites counts as 0; with sites it is used as given", {
  # No sites: 1e-12 means no clustering, with no error or warning
  expect_silent(res <- design_threshold(threshold = 0.05, prevalence = 0.10, icc = 1e-12))
  expect_equal(res$deff, 1)
  expect_equal(res$n, 191)
  # With 20 people per site, 1e-12 is used: deff is exactly 1 + 19 * 1e-12,
  # not 1 (expect_equal would treat the two as the same)
  r <- design_threshold(threshold = 0.05, prevalence = 0.10, n_per_site = 20, icc = 1e-12)
  expect_identical(r$deff, 1 + (20 - 1) * 1e-12)
  expect_gt(r$deff, 1)
  expect_equal(r$icc, 1e-12)
})


# ---------------------------------------------------------------------------
# Finite population correction (DT-F-)
# ---------------------------------------------------------------------------

test_that("DT-F-1: a small population reduces n and leaves n_eff alone; a huge one changes nothing", {
  r_fpc <- design_threshold(threshold = 0.05, prevalence = 0.10, fpc_N = 500)
  expect_lt(r_fpc$n, 191)
  expect_equal(r_fpc$fpc_N, 500)
  fpc300 <- design_threshold(threshold = 0.05, prevalence = 0.10, fpc_N = 300)
  expect_equal(fpc300$n_eff, 191)
  expect_lt(fpc300$n, fpc300$n_eff)
  expect_lt(design_threshold(threshold = 0.10, prevalence = 0.20, fpc_N = 200)$n,
            design_threshold(threshold = 0.10, prevalence = 0.20)$n)
  expect_equal(design_threshold(threshold = 0.05, prevalence = 0.10, fpc_N = 1e8)$n, 191)
})

test_that("DT-F-2: a population of 1 or 2 leaves too little room", {
  # The test needs at least one person left untested
  expect_error(design_threshold(threshold = 0.05, prevalence = 0.10, fpc_N = 1),
               "`fpc_N` (1) leaves no room for a design: the test needs at least one person in the population who is not tested.",
               fixed = TRUE)
  # With 2, at most 1 person can be tested; that gives at most 10 percent power
  expect_error(design_threshold(threshold = 0.05, prevalence = 0.10, fpc_N = 2),
               "Target power of 80 percent cannot be reached within the population (`fpc_N`) with these settings. Use a prevalence further from the threshold, or lower the target power.",
               fixed = TRUE)
  # 50 versus 95 percent in a population of 20 is reachable
  expect_lte(design_threshold(threshold = 0.5, prevalence = 0.95, fpc_N = 20)$n, 19)
})

test_that("DT-F-3: with sites and a population, deff matches the reported people per site", {
  b <- design_threshold(threshold = 0.05, prevalence = 0.10, n_sites = 50, icc = 0.05, fpc_N = 5000)
  expect_equal(b$deff, 1 + (b$n_per_site - 1) * 0.05)
  expect_silent(design_threshold(threshold = 0.05, prevalence = 0.10, n_sites = 20,
                                 fpc_N = 5000, icc = 0.05))
})


# ---------------------------------------------------------------------------
# Validation (DT-V-)
# ---------------------------------------------------------------------------

test_that("DT-V-1: threshold and prevalence must be strictly between 0 and 1", {
  # Valid values: 0.05 and 0.10 (DT-1); 0.001 and 0.999 (DT-10)
  expect_error(design_threshold(0, 0.10),
               "`threshold` must be in (0, 1) (got 0). No prevalence can be below 0 or above 1, so there is nothing to test. Use the prevalence that would lead to action, for example 0.05.", fixed = TRUE)
  expect_error(design_threshold(1, 0.10),
               "`threshold` must be in (0, 1) (got 1). No prevalence can be below 0 or above 1, so there is nothing to test. Use the prevalence that would lead to action, for example 0.05.", fixed = TRUE)
  expect_error(design_threshold(-0.1, 0.10),
               "`threshold` must be in (0, 1) (got -0.1). `threshold` is a prevalence and cannot be negative.",
               fixed = TRUE)
  expect_error(design_threshold(0.05, 0),
               "`prevalence` must be in (0, 1) (got 0). Use a value strictly between 0 and 1.", fixed = TRUE)
  expect_error(design_threshold(0.05, 1.2),
               "`prevalence` must be in (0, 1) (got 1.2). `prevalence` is a fraction of the population and cannot be greater than 1.",
               fixed = TRUE)
  expect_error(design_threshold(0.05, -0.1),
               "`prevalence` must be in (0, 1) (got -0.1). `prevalence` is a fraction of the population and cannot be negative.",
               fixed = TRUE)
  # Large values are printed in full, not as 1e+05
  expect_error(design_threshold(1e5, 0.10), "`threshold` must be in (0, 1) (got 100000).", fixed = TRUE)
})

test_that("DT-V-2: power at or below the false-alarm rate is rejected; below 0.5 warns", {
  for (bad in c(1, 1.2)) {
    expect_error(design_threshold(0.05, 0.10, power = bad),
                 paste0("`power` must be in (0, 1) (got ", bad, "). For example, use 0.80 for 80 percent power."),
                 fixed = TRUE)
  }
  # At 95 percent confidence the false-alarm rate is 0.05: 0.05 and below,
  # including 0 and negative values, are rejected
  for (bad in c(-0.2, 0, 0.03, 0.05)) {
    expect_error(design_threshold(0.05, 0.10, power = bad),
                 paste0("`power` must be greater than 1 - `conf_level` = 0.05 (got ", bad, "). Even when the ",
                        "prevalence is exactly at the threshold, the test reports a difference 5 percent of the ",
                        "time, so a lower power is not meaningful."),
                 fixed = TRUE)
  }
  # The false-alarm rate follows conf_level; equal only up to rounding
  # (1 - 0.8 = 0.19999999999999996) counts as equal
  expect_error(design_threshold(0.05, 0.10, power = 0.2, conf_level = 0.8),
               "`power` must be greater than 1 - `conf_level` = 0.2 (got 0.2).", fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, power = 1 - 0.8, conf_level = 0.8),
               "`power` must be greater than 1 - `conf_level` = 0.2 (got 0.2).", fixed = TRUE)
  # A power truly above the rate is allowed, however close
  expect_no_error(suppressWarnings(design_threshold(0.05, 0.10, power = 0.05 + 1e-8)))
  # Between the false-alarm rate and 0.5: warns, and the result is still
  # exact: 55 people and every larger number give at least 30 percent power
  # (54 gives less), and 28 is the first to reach it
  expect_warning(res <- design_threshold(0.05, 0.10, power = 0.3),
                 "`power` = 0.3 is below 0.5, so the study would be more likely to miss a real difference than to find it.",
                 fixed = TRUE)
  expect_equal(c(res$n_first, res$n), c(28, 55))
  expect_warning(design_threshold(0.05, 0.10, power = 0.45), "is below 0.5", fixed = TRUE)
  # A p-value exactly 0.05 does not count as significant, so 1 person (whose
  # one positive gives p = 0.05 exactly) has no power: the first size to
  # reach 10 percent is 6, as a direct binomial calculation gives
  expect_equal(suppressWarnings(design_threshold(0.05, 0.10, power = 0.1))$n_first, 6)
  # 0.5 itself works with no warning (DT-2)
  expect_no_warning(res <- design_threshold(0.05, 0.10, power = 0.5))
  expect_equal(res$n, 97)
})

test_that("DT-V-3: alternative must be one of the three names", {
  expect_error(design_threshold(0.05, 0.10, alternative = "up"),
               "`alternative` must be one of 'greater', 'less', or 'two.sided' (got 'up').", fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, alternative = NA),
               "`alternative` must be one of 'greater', 'less', or 'two.sided' (got class `logical`, length 1).",
               fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, alternative = c("a", "b")),
               "(got class `character`, length 2).", fixed = TRUE)
})

test_that("DT-V-4: the prevalence must be on the side the alternative says, and must differ from the threshold", {
  expect_error(design_threshold(0.10, 0.05, alternative = "greater"),
               "`prevalence` (0.05) must be greater than `threshold` (0.1) when alternative = 'greater'. Use alternative = 'less'.",
               fixed = TRUE)
  expect_error(design_threshold(0.10, 0.15, alternative = "less"),
               "`prevalence` (0.15) must be less than `threshold` (0.1) when alternative = 'less'. Use alternative = 'greater'.",
               fixed = TRUE)
  # Equal, or within 1e-8, for every alternative
  for (alt in c("greater", "less", "two.sided")) {
    expect_error(design_threshold(0.10, 0.10, alternative = alt),
                 "`prevalence` (0.1) and `threshold` (0.1) are equal, or nearly so, so there is no difference to detect. Use a prevalence further from the threshold.",
                 fixed = TRUE)
  }
  expect_error(design_threshold(0.05, 0.05 + 1e-9, alternative = "two.sided"),
               "`prevalence` (0.050000001) and `threshold` (0.05) are equal, or nearly so",
               fixed = TRUE)
})

test_that("DT-V-5: a vector, text, TRUE/FALSE, a bare NA or NULL in a single-number argument is rejected, naming the class", {
  expect_error(design_threshold(c(0.05, 0.10), 0.15),
               "`threshold` must be a single number in (0, 1) (got class `numeric`, length 2).", fixed = TRUE)
  expect_error(design_threshold("0.05", 0.10),
               "`threshold` must be a single number in (0, 1) (got class `character`, length 1).", fixed = TRUE)
  expect_error(design_threshold(0.05, c(0.10, 0.15)),
               "`prevalence` must be a single number in (0, 1) (got class `numeric`, length 2).", fixed = TRUE)
  expect_error(design_threshold(0.05, NA),
               "`prevalence` must be a single number in (0, 1) (got class `logical`, length 1).", fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, power = TRUE),
               "`power` must be a single number in (0, 1) (got class `logical`, length 1).", fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, power = c(0.8, 0.9)),
               "`power` must be a single number in (0, 1) (got class `numeric`, length 2).", fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, power = "0.8"),
               "`power` must be a single number in (0, 1) (got class `character`, length 1).", fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, sensitivity = TRUE),
               "`sensitivity` must be a single number in (0, 1] (got class `logical`, length 1). Note: `TRUE`, `FALSE` and a plain `NA` are logical, not numeric. Pass 1 if the test never misses a case.",
               fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, specificity = NULL),
               "`specificity` must be a single number in (0, 1] (got class `NULL`, length 0). Note: `TRUE`, `FALSE` and a plain `NA` are logical, not numeric. Pass 1 if the test never gives a false positive.",
               fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, specificity = "x"),
               "`specificity` must be a single number in (0, 1] (got class `character`, length 1).", fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, conf_level = "0.95"),
               "`conf_level` must be a single number in [0.5, 1) (got class `character`, length 1).", fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, icc = FALSE),
               "`icc` must be a single number in [0, 1] (got class `logical`, length 1).", fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, icc = "0"),
               "`icc` must be a single number in [0, 1] (got class `character`, length 1).", fixed = TRUE)
})

test_that("DT-V-6: NaN, Inf or -Inf in a single-number argument gives the finite-value error", {
  expect_error(design_threshold(NaN, 0.10), "`threshold` must be a single finite number (got NaN).", fixed = TRUE)
  expect_error(design_threshold(NA_real_, 0.10), "`threshold` must be a single finite number (got NA).", fixed = TRUE)
  expect_error(design_threshold(0.05, Inf), "`prevalence` must be a single finite number (got Inf).", fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, power = NaN), "`power` must be a single finite number (got NaN).", fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, sensitivity = NaN), "`sensitivity` must be a single finite number (got NaN).", fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, specificity = -Inf), "`specificity` must be a single finite number (got -Inf).", fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, conf_level = Inf), "`conf_level` must be a single finite number (got Inf).", fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, icc = NaN), "`icc` must be a single finite number (got NaN).", fixed = TRUE)
})

test_that("DT-V-7: sensitivity and specificity range checks give the right message for each bad value", {
  # Valid values: 0.9 and 0.95 (DT-8)
  expect_error(design_threshold(0.05, 0.10, sensitivity = 1.2),
               "`sensitivity` must be in (0, 1] (got 1.2). `sensitivity` is a diagnostic probability and cannot be greater than 1.",
               fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, sensitivity = 0),
               "`sensitivity` must be in (0, 1] (got 0). A sensitivity of 0 means everyone who truly has the condition tests negative.",
               fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, sensitivity = -0.5),
               "`sensitivity` must be in (0, 1] (got -0.5). `sensitivity` is a diagnostic probability and cannot be negative.",
               fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, specificity = 1.5),
               "`specificity` must be in (0, 1] (got 1.5). `specificity` is a diagnostic probability and cannot be greater than 1.",
               fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, specificity = -0.5),
               "`specificity` must be in (0, 1] (got -0.5). `specificity` is a diagnostic probability and cannot be negative.",
               fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, specificity = 0),
               "`specificity` must be in (0, 1] (got 0). A specificity of 0 means everyone who truly does not have the condition tests positive.",
               fixed = TRUE)
})

test_that("DT-V-8: sensitivity + specificity of 1 or less is rejected, including 1 plus rounding", {
  expect_error(design_threshold(0.05, 0.10, sensitivity = 0.4, specificity = 0.4),
               "`sensitivity` + `specificity` must exceed 1 (got 0.4 + 0.4 = 0.8). The difference the test can see is the true difference times `sensitivity` + `specificity` - 1, so with a sum of 1 or less the test cannot see the difference, or sees it reversed.",
               fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, sensitivity = 0.5, specificity = 0.5),
               "(got 0.5 + 0.5 = 1).", fixed = TRUE)
  # A sum above 1 only by rounding counts as 1
  expect_error(design_threshold(0.05, 0.10, sensitivity = 0.3,
                                specificity = 0.7 + 2 * .Machine$double.eps),
               "(got 0.3 + 0.7 = 1).", fixed = TRUE)
})

test_that("DT-V-9: sensitivity + specificity above 1 but below 1.1 warns, showing the full sum", {
  expect_warning(res <- design_threshold(0.05, 0.10, sensitivity = 0.53, specificity = 0.56),
                 "`sensitivity` + `specificity` = 1.09, which is below 1.1. The test then barely tells apart people with and without the condition, so the required n is very large, and small changes in the assumed sensitivity or specificity will change n a lot.",
                 fixed = TRUE)
  expect_equal(res$n, 75887)
  # 0.15 + 0.95 is 1.1 up to rounding: no warning
  expect_silent(design_threshold(0.05, 0.10, sensitivity = 0.15, specificity = 0.95))
  # 0.51 + 0.50 squeezes the gap to 0.0005 on the measured scale: millions
  # of people, beyond the 2 million limit
  expect_error(suppressWarnings(design_threshold(0.05, 0.10, sensitivity = 0.51, specificity = 0.50)),
               "These settings need more than 2 million people (about 6,180,000)", fixed = TRUE)
})

test_that("DT-V-10: icc and conf_level outside their ranges are rejected", {
  expect_error(design_threshold(0.05, 0.10, n_per_site = 10, icc = 1.5),
               "`icc` must be in [0, 1] (got 1.5). `icc` is a correlation and cannot be greater than 1.", fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, n_per_site = 10, icc = -0.1),
               "`icc` must be in [0, 1] (got -0.1). `icc` is a correlation and cannot be negative.", fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, conf_level = 1),
               "`conf_level` must be in [0.5, 1) (got 1). For example, use 0.95 for a 95% confidence interval.", fixed = TRUE)
  for (bad in c(0, 0.05, 0.3)) {
    expect_error(design_threshold(0.05, 0.10, conf_level = bad),
                 paste0("`conf_level` must be at least 0.5 (got ", bad, "). Below 0.5, when the prevalence is exactly at the threshold, the test can wrongly find it above or below the threshold more than half the time."),
                 fixed = TRUE)
  }
  # 0.5 itself is allowed
  expect_silent(design_threshold(0.05, 0.10, conf_level = 0.50))
})

test_that("DT-V-11: icc above 0 with no sites, or both n_sites and n_per_site, is rejected", {
  expect_error(design_threshold(0.05, 0.10, icc = 0.05),
               "`icc` is greater than 0 but no sites were given (`n_sites` and `n_per_site` are both `NULL`), so the design effect cannot be computed. Supply `n_sites` (fix the number of sites) or `n_per_site` (fix the number of people per site), or set `icc = 0` for a design without clustering.",
               fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, n_sites = 20, n_per_site = 10, icc = 0.05),
               "Supply at most one of `n_sites` or `n_per_site`, not both. `n_sites` fixes the number of sites and solves for people per site; `n_per_site` fixes the number of people per site and solves for the number of sites.",
               fixed = TRUE)
})

test_that("DT-V-12: n_sites, n_per_site and fpc_N must be single positive integers, naming the problem", {
  expect_error(design_threshold(0.05, 0.10, n_sites = 2.5, icc = 0.05),
               "`n_sites` must be a single finite positive integer (got 2.5). `n_sites` is the number of sites.", fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, n_sites = c(2, 3), icc = 0.05),
               "`n_sites` must be a single finite positive integer (got length = 2).", fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, n_per_site = c(10, 20)),
               "`n_per_site` must be a single finite positive integer (got length = 2). `n_per_site` is the fixed number of people sampled per site.",
               fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, n_sites = "10", icc = 0.05),
               "`n_sites` must be a single finite positive integer (got class `character`).", fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, n_per_site = 0, icc = 0.05),
               "`n_per_site` must be a single finite positive integer (got 0). `n_per_site` is the fixed number of people sampled per site.",
               fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, n_per_site = TRUE, icc = 0.05),
               "`n_per_site` must be a single finite positive integer (got class `logical`).", fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, fpc_N = 500.5),
               "`fpc_N` must be a single finite positive integer (got 500.5). `fpc_N` is the total population size. Set `fpc_N = NULL` to skip the finite-population correction.",
               fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, fpc_N = c(100, 200)),
               "`fpc_N` must be a single finite positive integer (got length = 2).", fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, fpc_N = "500"),
               "`fpc_N` must be a single finite positive integer (got class `character`).", fixed = TRUE)
  # A value within rounding of a whole number gets a hint to use round()
  expect_error(design_threshold(0.05, 0.15, n_sites = 100 * 1.1),
               "`n_sites` must be a single finite positive integer (got 110, which is not exactly an integer; use round(n_sites)). `n_sites` is the number of sites.",
               fixed = TRUE)
  expect_error(design_threshold(0.05, 0.15, n_per_site = 100 * 1.1),
               "`n_per_site` must be a single finite positive integer (got 110, which is not exactly an integer; use round(n_per_site)). `n_per_site` is the fixed number of people sampled per site.",
               fixed = TRUE)
  expect_error(design_threshold(0.05, 0.15, fpc_N = 100 * 1.1),
               "`fpc_N` must be a single finite positive integer (got 110, which is not exactly an integer; use round(fpc_N)). `fpc_N` is the total population size. Set `fpc_N = NULL` to skip the finite-population correction.",
               fixed = TRUE)
  # A real fraction gets no hint
  expect_error(design_threshold(0.05, 0.15, fpc_N = 1e10 + 0.5),
               "`fpc_N` must be a single finite positive integer (got 10000000000.5). `fpc_N` is the total population size. Set `fpc_N = NULL` to skip the finite-population correction.",
               fixed = TRUE)
  # Above 2^53, the largest whole number R can store exactly
  expect_error(design_threshold(0.05, 0.10, fpc_N = 1e20),
               "`fpc_N` (1e+20) is larger than 2^53 = 9007199254740992, the largest whole number R can store exactly. No real population is this large: set `fpc_N = NULL`.",
               fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, n_sites = 1e17, icc = 0.05),
               "`n_sites` (1e+17) is larger than 2^53 = 9007199254740992, the largest whole number R can store exactly. No real study has this many sites.",
               fixed = TRUE)
  # Just above 2^53 but below 1e16: printed in full
  expect_error(design_threshold(0.05, 0.10, n_sites = 9.1e15, icc = 0.05),
               "`n_sites` (9100000000000000) is larger than 2^53 = 9007199254740992, the largest whole number R can store exactly. No real study has this many sites.",
               fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, n_per_site = 1e17, icc = 0.05),
               "`n_per_site` (1e+17) is larger than 2^53 = 9007199254740992, the largest whole number R can store exactly. No real site has this many people.",
               fixed = TRUE)
  # Valid whole numbers
  expect_equal(design_threshold(0.05, 0.10, fpc_N = 1e5)$n, 191)
  expect_silent(design_threshold(0.05, 0.10, fpc_N = 500))
})

test_that("DT-V-13: more sites, or more people per site, than the population is rejected", {
  expect_error(design_threshold(0.05, 0.10, n_sites = 100, fpc_N = 50, icc = 0.05),
               "`n_sites` (100) exceeds `fpc_N` (50): you cannot have more sites than people in the population.",
               fixed = TRUE)
  expect_error(design_threshold(0.05, 0.10, n_per_site = 100, fpc_N = 50, icc = 0.05),
               "`n_per_site` (100) exceeds `fpc_N` (50): a site cannot have more people than there are in the whole population.",
               fixed = TRUE)
  # Valid: 20 sites from a population of 5000 (DT-F-3)
})

test_that("DT-V-14: an imperfect test that brings the two values almost together is rejected", {
  expect_error(suppressWarnings(design_threshold(0.05, 0.10, sensitivity = 0.5,
                                                 specificity = 0.5 + 1e-7)),
               "The measured prevalence and threshold differ by only 0.000000005 (0.49999991 and 0.499999905), so the required sample size would be extremely large. Use a prevalence further from the threshold, or a more accurate test.",
               fixed = TRUE)
})


# ---------------------------------------------------------------------------
# Return value (DT-R-)
# ---------------------------------------------------------------------------

test_that("DT-R-1: return list contains all expected fields, in order, for all three designs", {
  expected <- c("n", "n_first", "n_eff", "n_sites", "n_per_site", "threshold",
                "prevalence", "apparent_prev", "threshold_app", "power",
                "alternative", "conf_level", "sensitivity", "specificity", "icc",
                "deff", "fpc_N")
  expect_named(design_threshold(0.10, 0.20), expected, ignore.order = FALSE)
  expect_named(design_threshold(0.10, 0.20, n_per_site = 10, icc = 0.05), expected,
               ignore.order = FALSE)
  expect_named(design_threshold(0.10, 0.20, n_sites = 10, icc = 0.05), expected,
               ignore.order = FALSE)
})

test_that("DT-R-2: inputs are echoed back as plain numbers", {
  res <- design_threshold(c(a = 0.05), 0.10, power = 0.9, alternative = "greater",
                          sensitivity = 0.9, specificity = 0.95, conf_level = 0.9,
                          n_per_site = 10L, icc = 0.05, fpc_N = 5000)
  expect_identical(res$threshold, 0.05)
  expect_equal(c(res$prevalence, res$power, res$conf_level, res$sensitivity,
                 res$specificity, res$icc, res$n_per_site, res$fpc_N),
               c(0.10, 0.9, 0.9, 0.9, 0.95, 0.05, 10, 5000))
  expect_equal(res$alternative, "greater")
})


# ---------------------------------------------------------------------------
# Round trip with test_threshold() (DT-T-)
# ---------------------------------------------------------------------------

test_that("DT-T-1: simulated power of test_threshold() at the design n reaches the target", {
  set.seed(8)
  d <- design_threshold(threshold = 0.05, prevalence = 0.10)
  rej <- vapply(1:3000, function(i)
    test_threshold(rbinom(1, d$n, 0.10), d$n, 0.05)$reject, logical(1))
  expect_gt(mean(rej), 0.78)   # exact power at 191 is at least 0.80; 3000 draws
})

test_that("DT-T-2: the clustered design gives a significant result at the planned prevalence when icc is supplied", {
  # 10 sites of 91 people, icc 0.05: deff 5.5; planned prevalence 0.10 is 9 of 91 per site
  d <- design_threshold(threshold = 0.05, prevalence = 0.10, n_sites = 10, icc = 0.05)
  expect_equal(c(d$n_sites, d$n_per_site, d$deff), c(10, 91, 5.5))
  res <- test_threshold(rep(9, 10), rep(d$n_per_site, 10), threshold = 0.05, icc = 0.05)
  expect_equal(res$deff, d$deff)
  expect_true(res$reject)
})

test_that("DT-T-4: a one-site design is analysed with its icc and gives a significant result", {
  d <- design_threshold(threshold = 0.05, prevalence = 0.10, n_sites = 1, icc = 0.005)
  expect_equal(c(d$n_sites, d$n, d$deff), c(1, 962, 5.805))
  res <- test_threshold(round(0.10 * d$n), d$n, threshold = 0.05, icc = 0.005)
  expect_equal(res$deff, d$deff)
  expect_true(res$reject)
})

test_that("DT-T-5: a design for an imperfect test gives a significant result at the planned apparent prevalence", {
  d <- design_threshold(threshold = 0.05, prevalence = 0.10, sensitivity = 0.95,
                        specificity = 0.98, n_sites = 20, icc = 0.05)
  expect_equal(c(d$n_sites, d$n_per_site, d$deff), c(20, 28, 2.35))
  # apparent prevalence 0.113: 3 of 28 per site
  res <- test_threshold(rep(3, 20), rep(d$n_per_site, 20), threshold = 0.05,
                        sensitivity = 0.95, specificity = 0.98, icc = 0.05)
  expect_equal(res$deff, d$deff)
  expect_true(res$reject)
})

test_that("DT-P-1: the exact power at its edges, checked by hand", {
  P <- MMSpower:::.thr_power
  # A p-value exactly at the level is not below it: 5 positives of 5 at a
  # threshold of 0.5 has p-value 1/32 exactly, so at level 1/32 no count is
  # significant (0.9^5 = 0.59049 if it were)
  expect_equal(P(5, 1, NULL, 0.5, 0.9, "greater", 1 / 32), 0)
  # and 0 positives of 5 on the "less" side
  expect_equal(P(5, 1, NULL, 0.5, 0.1, "less", 1 / 32), 0)
  # No count significant at all: with 2 people and a threshold of 0.5 the
  # smallest p-value is 0.25, on either side
  expect_equal(P(2, 1, NULL, 0.5, 0.9, "greater", 0.05), 0)
  expect_equal(P(2, 1, NULL, 0.5, 0.1, "less", 0.05), 0)
  # The same with sites (4 people, design effect 2: effective size 2), where
  # the counts are fractions
  expect_equal(P(4, 2, NULL, 0.5, 0.9, "greater", 0.05), 0)
  expect_equal(P(4, 2, NULL, 0.5, 0.1, "less", 0.05), 0)
  # fpc_N: 100 tested of 150 is an effective size of 100 * 149 / 50 = 298;
  # the first significant count is 8, and its chance at 8 percent is 0.5306583
  # (0.530556 with 150 / 50 in place of 149 / 50)
  expect_equal(P(100, 1, 150, 0.05, 0.08, "greater", 0.05), 0.5306583, tolerance = 1e-6)
})
