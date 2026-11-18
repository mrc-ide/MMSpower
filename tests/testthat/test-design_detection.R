# Hand-checked reference cases for design_detection()
#
#   q  = prevalence * sensitivity        (chance one person is a case and tests positive)
#   n0 = log(1 - detection_prob) / log(1 - q)
#
# DET-1 (Module 7 worked example: one health facility, 5 percent, 80 percent):
#   n0 = log(0.2) / log(0.95) = -1.609438 / -0.051293 = 31.38 -> 32
#   (the slide gives 32 samples)
#
# DET-2 (Module 7: 139 samples):
#   1 - 0.98^139 = 0.9397 (2 percent);  1 - 0.99^139 = 0.7527 (1 percent,
#   the slide's "75% powered to detect a 1% prevalence")
#
# With sites (beta-binomial, as in DRpower's get_power_presence()):
#   a = q (1/icc - 1), b = (1 - q) (1/icc - 1)
#   chance of no positive among m people at one site:
#   P0(m) = B(a, b + m) / B(a, b) = prod over j = 0..m-1 of (b + j) / (a + b + j)
#   chance of at least one positive at k sites: 1 - P0(m)^k
#
# DET-C-1 (n_per_site = 10, icc = 0.05): a = 0.95, b = 18.05
#   1 - P0(10)^3 = 71.57 percent, 1 - P0(10)^4 = 81.31 percent -> 4 sites, 40 people
#
# DET-C-2 (n_sites = 10, icc = 0.05):
#   1 - P0(3)^10 = 76.83 percent, 1 - P0(4)^10 = 85.10 percent -> 4 per site, 40 people
#
# DET-F-1 (population 200, 5 percent, 95 percent):
#   10 cases. 1 - choose(190, 51) / choose(200, 51) = 0.9514 (enough),
#   1 - choose(190, 50) / choose(200, 50) = 0.9479 (not enough) -> 51


# ---------------------------------------------------------------------------
# Core scenarios (DET-)
# ---------------------------------------------------------------------------

test_that("DET-1: the Module 7 example (5 percent, 80 percent) gives the slide's 32", {
  res <- design_detection(prevalence = 0.05)
  expect_equal(res$n,              32)
  expect_equal(res$n_eff,          32)
  expect_equal(res$q,              0.05)
  expect_equal(res$detection_prob, 0.80)
  expect_equal(res$deff,           1)
  expect_equal(res$mode,           "solve_n")
  expect_null(res$n_sites)
  expect_null(res$n_per_site)
  expect_null(res$n_cases)
})

test_that("DET-2: in reverse, 139 people give the slide's 94 and 75 percent", {
  r2 <- design_detection(prevalence = 0.02, n = 139)
  r1 <- design_detection(prevalence = 0.01, n = 139)
  expect_equal(r2$detection_prob, 1 - 0.98^139)
  expect_equal(r1$detection_prob, 1 - 0.99^139)
  expect_equal(round(r1$detection_prob, 3), 0.753)
  expect_equal(r1$mode, "solve_detection_prob")
})

test_that("DET-3: a higher target or a rarer case needs more people", {
  expect_equal(design_detection(0.05, detection_prob = 0.95)$n, 59)
  expect_equal(design_detection(0.05, detection_prob = 0.99)$n, 90)
  expect_equal(design_detection(0.01)$n,                        161)
})

test_that("DET-4: an imperfect test needs more people", {
  # q = 0.05 * 0.8 = 0.04; log(0.2) / log(0.96) = 39.43 -> 40
  res <- design_detection(0.05, sensitivity = 0.8)
  expect_equal(res$q, 0.04)
  expect_equal(res$n, 40)
})

test_that("DET-5: forward and reverse agree: n reaches the target and n - 1 does not", {
  expect_gte(design_detection(0.05, n = 32)$detection_prob, 0.80)
  expect_lt(design_detection(0.05, n = 31)$detection_prob,  0.80)
})

test_that("DET-6: integers stored as R integers (30L) give the same result", {
  expect_equal(design_detection(0.05, n_per_site = 10L, icc = 0.05),
               design_detection(0.05, n_per_site = 10,  icc = 0.05))
  expect_equal(design_detection(0.05, n = 139L), design_detection(0.05, n = 139))
  expect_equal(design_detection(0.05, fpc_N = 200L), design_detection(0.05, fpc_N = 200))
})


# ---------------------------------------------------------------------------
# Clustering (DET-C-)
# ---------------------------------------------------------------------------

# Chance of no positive at one site of m people, worked out as a product
# (a second way, not using lbeta as the function does)
p0_direct <- function(m, q, icc) {
  a <- q * (1 / icc - 1); b <- (1 - q) * (1 / icc - 1)
  prod((b + 0:(m - 1)) / (a + b + 0:(m - 1)))
}

test_that("DET-C-1: a fixed number of people per site solves for the number of sites", {
  res <- design_detection(0.05, n_per_site = 10, icc = 0.05)
  expect_equal(res$n,       40)
  expect_equal(res$n_eff,   32)
  expect_equal(res$n_sites, 4)
  # 4 sites reach 80 percent, 3 do not
  expect_gte(1 - p0_direct(10, 0.05, 0.05)^4, 0.80)
  expect_lt(1 - p0_direct(10, 0.05, 0.05)^3,  0.80)
  # deff: 40 people with sites do what 40 / 1.2234 = 32.7 people without do
  expect_equal(res$deff, 40 * log(0.95) / (4 * log(p0_direct(10, 0.05, 0.05))))
  expect_equal(round(res$deff, 6), 1.223414)
  # A target exactly at the chance of 4 sites needs 4 sites, not 5
  expect_equal(design_detection(0.05, 1 - p0_direct(10, 0.05, 0.05)^4,
                                n_per_site = 10, icc = 0.05)$n_sites, 4)
})

test_that("DET-C-2: a fixed number of sites solves for people per site", {
  res <- design_detection(0.05, n_sites = 10, icc = 0.05)
  expect_equal(res$n,          40)
  expect_equal(res$n_per_site, 4)
  expect_gte(1 - p0_direct(4, 0.05, 0.05)^10, 0.80)
  expect_lt(1 - p0_direct(3, 0.05, 0.05)^10,  0.80)
  expect_equal(round(res$deff, 6), 1.077626)
  # A target exactly at the chance of 4 per site needs 4 per site, not 5
  expect_equal(design_detection(0.05, 1 - p0_direct(4, 0.05, 0.05)^10,
                                n_sites = 10, icc = 0.05)$n_per_site, 4)
})

test_that("DET-C-3: people per site and number of sites match DRpower", {
  # Expected values from DRpower::get_sample_size_presence(n_clust,
  # target_power = 0.8, prevalence, ICC)
  for (s in list(c(5, 0.01, 0.10, 294), c(5, 0.01, 0.05, 82), c(10, 0.02, 0.10, 13),
                 c(3, 0.05, 0.20, 47),  c(20, 0.01, 0.05, 10), c(10, 0.05, 0.30, 6),
                 c(4, 0.10, 0.60, 101))) {
    expect_equal(design_detection(s[2], n_sites = s[1], icc = s[3])$n_per_site, s[4])
  }
  # Number of sites, target 90 percent: the smallest k whose chance in
  # DRpower::get_power_presence(rep(n_per_site, k), prevalence, ICC) is above
  # the target (6 sites: 91.92 percent, 5 sites: 87.71; 14: 90.52, 13: 88.79;
  # 30: 90.32, 29: 89.54)
  for (s in list(c(10, 0.05, 0.05, 6), c(25, 0.02, 0.20, 14), c(50, 0.01, 0.30, 30))) {
    expect_equal(design_detection(s[2], 0.90, n_per_site = s[1], icc = s[3])$n_sites, s[4])
  }
})

test_that("DET-C-4: a high ICC needs many people per site; at ICC 1 too few sites is rejected", {
  # The old design-effect method said this could not be reached; it can,
  # with 223 per site
  expect_equal(design_detection(0.05, 0.95, n_sites = 2, icc = 0.1)$n_per_site, 223)
  # ICC 1: each site counts as one person, so 5 sites give at most
  # 1 - 0.95^5 = 22.62 percent
  expect_error(
    design_detection(0.05, n_sites = 5, icc = 1),
    "A detection probability of 80 percent cannot be reached with 5 sites and ICC = 1: with an ICC of 1 everyone at a site has the same result, so these sites can give at most 22.62191 percent, however many people each site has. Increase `n_sites`, lower the ICC, or lower `detection_prob`.",
    fixed = TRUE)
  # Exactly that maximum is reached with 1 person per site
  res <- design_detection(0.05, 1 - 0.95^5, n_sites = 5, icc = 1)
  expect_equal(res$n_per_site, 1)
  expect_equal(res$deff, 1)
  # ICC below 1 but so high that more than a trillion people per site would
  # be needed (beyond that, the answer could not be worked out exactly)
  expect_error(
    design_detection(0.001, n_sites = 3, icc = 0.5),
    "A detection probability of 80 percent cannot be reached with 3 sites and ICC = 0.5: it would need more than a trillion (1000000000000) people per site (that many give 8.114521 percent). Increase `n_sites`, lower the ICC, or lower `detection_prob`.",
    fixed = TRUE)
})

test_that("DET-C-5: more sites than people needed with no clustering is rejected", {
  expect_error(
    design_detection(0.05, n_sites = 40, icc = 0.05),
    "`n_sites` (40) is more than the number of people needed with no clustering (32), so one person per site would already be more than enough. Use `n_sites` of 32 or fewer, or supply `n_per_site` instead.",
    fixed = TRUE)
  # 32 sites of 1 person is the unclustered design
  res <- design_detection(0.05, n_sites = 32, icc = 0.05)
  expect_equal(res$n_per_site, 1)
  expect_equal(res$n,          32)
  expect_equal(res$deff,       1)
  # 31 sites need 2 people each
  expect_equal(design_detection(0.05, n_sites = 31, icc = 0.05)$n, 62)
})

test_that("DET-C-10: one site with more people than needed is rejected", {
  msg <- function(nps, needed)
    paste0("`n_per_site` (", nps, ") is larger than the number of people needed (",
           needed, "), so one site would be more than enough. Use a smaller ",
           "`n_per_site`, or set `n_sites = 1` to find how many people one site needs.")
  # icc 0: 32 people are needed (DET-1), so a site of 33 is too big; 32 is fine
  expect_error(design_detection(0.05, n_per_site = 33), msg(33, 32), fixed = TRUE)
  res <- design_detection(0.05, n_per_site = 32)
  expect_equal(c(res$n, res$n_sites, res$n_per_site), c(32, 1, 32))
  # icc 0.05: one site of 81 reaches 80 percent, one site of 80 does not
  expect_gte(1 - p0_direct(81, 0.05, 0.05), 0.80)
  expect_lt(1 - p0_direct(80, 0.05, 0.05),  0.80)
  expect_error(design_detection(0.05, n_per_site = 82, icc = 0.05), msg(82, 81),
               fixed = TRUE)
  res <- design_detection(0.05, n_per_site = 81, icc = 0.05)
  expect_equal(c(res$n, res$n_sites, res$n_per_site), c(81, 1, 81))
  # A target exactly at the chance of one site of 81 still needs 81, not 82
  expect_error(design_detection(0.05, 1 - p0_direct(81, 0.05, 0.05),
                                n_per_site = 82, icc = 0.05), msg(82, 81),
               fixed = TRUE)
  # The number in the message is what n_sites = 1 gives
  expect_equal(design_detection(0.05, n_sites = 1, icc = 0.05)$n_per_site, 81)
  # ICC 1: one person per site already gives q = 0.5, above a 30 percent target
  expect_error(design_detection(0.5, 0.3, n_per_site = 2, icc = 1), msg(2, 1),
               fixed = TRUE)
  # Population of 200 at 80 percent: 30 people are needed
  expect_equal(design_detection(0.05, fpc_N = 200)$n, 30)
  expect_error(design_detection(0.05, n_per_site = 31, fpc_N = 200), msg(31, 30),
               fixed = TRUE)
  expect_equal(design_detection(0.05, n_per_site = 30, fpc_N = 200)$n_sites, 1)
  # In reverse, a large site is fine
  expect_equal(design_detection(0.05, n = 30, n_per_site = 100)$n_sites, 1)
})

test_that("DET-C-6: icc = 0, or a tiny icc, with sites is the unclustered result spread over the sites", {
  for (icc in c(0, 1e-12)) {
    res <- design_detection(0.05, n_sites = 10, icc = icc)
    expect_equal(res$n,          32)
    expect_equal(res$deff,       1)
    expect_equal(res$n_per_site, 4)   # 32 / 10 = 3.2, rounded up
  }
  # A tiny icc without sites is not an error: it counts as 0
  expect_equal(design_detection(0.05, icc = 1e-12)$n, 32)
})

test_that("DET-C-7: in reverse, the chance matches DRpower, with the people spread as evenly as possible", {
  # 40 people at 10 sites: 4 each (DRpower: 85.1019 percent)
  r1 <- design_detection(0.05, n = 40, n_sites = 10, icc = 0.05)
  expect_equal(r1$detection_prob, 1 - p0_direct(4, 0.05, 0.05)^10)
  expect_equal(round(100 * r1$detection_prob, 4), 85.1019)
  expect_equal(r1$n_per_site, 4)
  # 43 people at 10 sites: 3 sites of 5 and 7 of 4
  r2 <- design_detection(0.05, n = 43, n_sites = 10, icc = 0.05)
  expect_equal(r2$detection_prob,
               1 - p0_direct(5, 0.05, 0.05)^3 * p0_direct(4, 0.05, 0.05)^7)
  expect_equal(r2$n_per_site, 5)
  # 45 people in sites of 10: 4 full sites and one of 5
  r3 <- design_detection(0.05, n = 45, n_per_site = 10, icc = 0.05)
  expect_equal(r3$detection_prob,
               1 - p0_direct(10, 0.05, 0.05)^4 * p0_direct(5, 0.05, 0.05))
  expect_equal(r3$n_sites, 5)
  # n_eff: the people with no clustering who give the same chance,
  # rounded up; deff = n / that number
  n_same <- log(1 - r3$detection_prob) / log(0.95)
  expect_equal(r3$n_eff, ceiling(n_same))
  expect_equal(r3$deff,  45 / n_same)
  # One person per site is the same as no clustering
  r4 <- design_detection(0.05, n = 30, n_per_site = 1, icc = 0.5)
  expect_equal(r4$detection_prob, 1 - 0.95^30)
  expect_equal(r4$n_eff, 30)
  expect_equal(r4$deff,  1)
})

test_that("DET-C-8: forward n reaches the target and one person per site (or one site) fewer does not", {
  for (s in list(c(0.05, 0.80, 0.05), c(0.05, 0.95, 0.20), c(0.02, 0.90, 0.10),
                 c(0.01, 0.80, 0.50))) {
    res <- design_detection(s[1], s[2], n_per_site = 10, icc = s[3])
    expect_gte(1 - p0_direct(10, s[1], s[3])^res$n_sites,       s[2])
    expect_lt(1 - p0_direct(10, s[1], s[3])^(res$n_sites - 1),  s[2])
    res <- design_detection(s[1], s[2], n_sites = 15, icc = s[3])
    expect_gte(1 - p0_direct(res$n_per_site, s[1], s[3])^15,     s[2])
    expect_lt(1 - p0_direct(res$n_per_site - 1, s[1], s[3])^15,  s[2])
    # and reverse mode gives the same chance
    expect_equal(design_detection(s[1], n = 15 * res$n_per_site, n_sites = 15, icc = s[3])$detection_prob,
                 1 - p0_direct(res$n_per_site, s[1], s[3])^15)
  }
})

test_that("DET-C-9: with an imperfect test, sites use q = prevalence x sensitivity", {
  # Same as DRpower with prevalence 0.05 x 0.8 = 0.04: 5 per site
  res <- design_detection(0.05, sensitivity = 0.8, n_sites = 10, icc = 0.05)
  expect_equal(res$n_per_site, 5)
  expect_equal(res$q, 0.04)
  expect_equal(res, modifyList(design_detection(0.04, n_sites = 10, icc = 0.05),
                               list(prevalence = 0.05, sensitivity = 0.8)))
})


# ---------------------------------------------------------------------------
# Finite population (DET-F-)
# ---------------------------------------------------------------------------

test_that("DET-F-1: a population of 200 gives the exact 51, and leaves n_eff alone", {
  res <- design_detection(0.05, 0.95, fpc_N = 200)
  expect_equal(res$n,       51)
  expect_equal(res$n_eff,   59)
  expect_equal(res$n_cases, 10)
  expect_gte(1 - choose(190, 51) / choose(200, 51), 0.95)
  expect_lt(1 - choose(190, 50) / choose(200, 50),  0.95)
})

test_that("DET-F-2: matches a brute-force search with R's dhyper, including an imperfect test", {
  brute <- function(p, d, N, se) {
    D <- floor(N * p + 1e-8)
    for (n in 1:N) {
      k <- 0:min(n, D)
      if (1 - sum(stats::dhyper(k, D, N - D, n) * (1 - se)^k) >= d) return(n)
    }
  }
  for (s in list(c(0.05, 0.80, 200, 1), c(0.02, 0.95, 500, 1), c(0.05, 0.95, 200, 0.8),
                 c(0.02, 0.80, 100, 1), c(0.10, 0.90, 80, 0.6))) {
    expect_equal(design_detection(s[1], s[2], sensitivity = s[4], fpc_N = s[3])$n,
                 brute(s[1], s[2], s[3], s[4]))
  }
})

test_that("DET-F-3: a target above what testing everyone can give is rejected", {
  # 5 cases, each missed with probability 0.7: testing all 100 people finds
  # a case with probability 1 - 0.7^5 = 83.2 percent
  expect_error(
    design_detection(0.05, 0.99, sensitivity = 0.3, fpc_N = 100),
    "A detection probability of 99 percent is not achievable in a population of 100 with 5 cases and `sensitivity` = 0.3: even testing everyone finds a case with probability 83.2 percent. Lower `detection_prob`, or use a more sensitive test.",
    fixed = TRUE)
  # With sensitivity 1, testing everyone always finds a case: 2 cases in 40
  expect_equal(design_detection(0.05, 0.95, fpc_N = 40)$n, 31)
})

test_that("DET-F-4: the number of cases is rounded down, without rounding errors", {
  # 0.29 * 100 is stored as 28.999999999999996 and must count as 29;
  # 0.07 * 100 is stored as 7.000000000000001 and must count as 7
  expect_equal(design_detection(0.29, fpc_N = 100)$n_cases, 29)
  expect_equal(design_detection(0.07, fpc_N = 100)$n_cases, 7)
  # 1.5 cases counts as 1
  expect_equal(design_detection(0.05, fpc_N = 30)$n_cases,  1)
  expect_error(
    design_detection(0.005, fpc_N = 100),
    "`fpc_N` (100) x `prevalence` (0.005) is less than 1: a population of 100 at this prevalence contains no case to find. Use a larger population or a higher prevalence.",
    fixed = TRUE)
})

test_that("DET-F-5: fpc_N cannot be used together with an icc above 0", {
  msg <- "`fpc_N` cannot be used together with an `icc` above 0 (got icc = 0.05): the calculation with sites assumes a very large population. Set `fpc_N = NULL`, or set `icc = 0`."
  expect_error(design_detection(0.05, 0.95, n_per_site = 10, icc = 0.05, fpc_N = 200),
               msg, fixed = TRUE)
  expect_error(design_detection(0.05, 0.95, n_sites = 10, icc = 0.05, fpc_N = 200),
               msg, fixed = TRUE)
  expect_error(design_detection(0.05, n = 40, n_sites = 10, icc = 0.05, fpc_N = 200),
               msg, fixed = TRUE)
  # With icc 0, sites and fpc_N together are fine: 51 people (DET-F-1),
  # 6 sites of 10
  res <- design_detection(0.05, 0.95, n_per_site = 10, fpc_N = 200)
  expect_equal(res$n,       51)
  expect_equal(res$n_sites, 6)
})

test_that("DET-F-6: a rounded design larger than the population warns", {
  # Population 40, 21 sites, icc 0: n = 22, 22 / 21 = 1.05 -> 2 per site,
  # so 42 people
  expect_warning(
    res <- design_detection(0.05, n_sites = 21, fpc_N = 40),
    "The rounded design (21 sites x 2 = 42) exceeds `fpc_N` = 40: the required n of 22 is close to a census of the population. Consider surveying the whole population instead.",
    fixed = TRUE)
  expect_equal(res$n, 22)
})

test_that("DET-F-7: in reverse, a finite population raises the chance, and testing everyone gives the maximum", {
  expect_gt(design_detection(0.05, n = 100, fpc_N = 200)$detection_prob,
            design_detection(0.05, n = 100)$detection_prob)
  expect_equal(design_detection(0.05, n = 200, fpc_N = 200)$detection_prob, 1)
  expect_equal(design_detection(0.05, n = 200, fpc_N = 200, sensitivity = 0.5)$detection_prob,
               1 - 0.5^10)
})


# ---------------------------------------------------------------------------
# Validation (DET-V-)
# ---------------------------------------------------------------------------

test_that("DET-V-1: a prevalence of 0, 1, or outside [0, 1] is rejected", {
  # Valid value: 0.05 (DET-1)
  expect_error(design_detection(0),
               "`prevalence` must be in (0, 1) (got 0). With a prevalence of 0 there is no case to find.",
               fixed = TRUE)
  expect_error(design_detection(1),
               "`prevalence` must be in (0, 1) (got 1). With a prevalence of 1 everyone is a case, so there is no sample size to plan.",
               fixed = TRUE)
  expect_error(design_detection(-0.1),
               "`prevalence` must be in (0, 1) (got -0.1). `prevalence` is a fraction of the population and cannot be negative.",
               fixed = TRUE)
  expect_error(design_detection(1.2),
               "`prevalence` must be in (0, 1) (got 1.2). `prevalence` is a fraction of the population and cannot be greater than 1.",
               fixed = TRUE)
})

test_that("DET-V-2: detection_prob outside (0, 1) is rejected, but only when solving for n", {
  # Valid values: 0.80 (DET-1), 0.95 and 0.99 (DET-3)
  for (bad in c(0, 1, -0.2, 80)) {
    expect_error(design_detection(0.05, detection_prob = bad),
                 paste0("`detection_prob` must be in (0, 1) (got ", bad, "). For example, use 0.80 for an 80 percent chance of finding at least one case."),
                 fixed = TRUE)
  }
  # In reverse detection_prob is not used, so it is not checked
  expect_no_error(design_detection(0.05, detection_prob = 80, n = 32))
  expect_no_error(design_detection(0.05, detection_prob = "x", n = 32))
})

test_that("DET-V-3: a vector, text, TRUE/FALSE, a bare NA or NULL in a single-number argument is rejected, naming the class", {
  expect_error(design_detection(c(0.05, 0.1)),
               "`prevalence` must be a single number in (0, 1) (got class `numeric`, length 2).", fixed = TRUE)
  expect_error(design_detection("0.05"),
               "`prevalence` must be a single number in (0, 1) (got class `character`, length 1).", fixed = TRUE)
  expect_error(design_detection(0.05, detection_prob = TRUE),
               "`detection_prob` must be a single number in (0, 1) (got class `logical`, length 1).", fixed = TRUE)
  expect_error(design_detection(0.05, detection_prob = NULL),
               "`detection_prob` must be a single number in (0, 1) (got class `NULL`, length 0).", fixed = TRUE)
  expect_error(design_detection(0.05, sensitivity = NA),
               "`sensitivity` must be a single number in (0, 1] (got class `logical`, length 1). Note: `TRUE`, `FALSE` and a plain `NA` are logical, not numeric. Pass 1 if the test never misses a case.",
               fixed = TRUE)
  expect_error(design_detection(0.05, icc = list(0.05)),
               "`icc` must be a single number in [0, 1] (got class `list`, length 1).", fixed = TRUE)
})

test_that("DET-V-4: NaN, Inf or -Inf in a single-number argument gives the finite-value error", {
  expect_error(design_detection(NaN), "`prevalence` must be a single finite number (got NaN).", fixed = TRUE)
  expect_error(design_detection(0.05, detection_prob = Inf),
               "`detection_prob` must be a single finite number (got Inf).", fixed = TRUE)
  expect_error(design_detection(0.05, sensitivity = -Inf),
               "`sensitivity` must be a single finite number (got -Inf).", fixed = TRUE)
  expect_error(design_detection(0.05, icc = NaN),
               "`icc` must be a single finite number (got NaN).", fixed = TRUE)
})

test_that("DET-V-5: sensitivity and icc range checks give the right message for each bad value", {
  # Valid values: sensitivity 1 (DET-1) and 0.8 (DET-4); icc 0.05 (DET-C-1)
  expect_error(design_detection(0.05, sensitivity = 0),
               "`sensitivity` must be in (0, 1] (got 0). A sensitivity of 0 means everyone who truly has the condition tests negative.",
               fixed = TRUE)
  expect_error(design_detection(0.05, sensitivity = -0.1),
               "`sensitivity` must be in (0, 1] (got -0.1). `sensitivity` is a diagnostic probability and cannot be negative.",
               fixed = TRUE)
  expect_error(design_detection(0.05, sensitivity = 1.1),
               "`sensitivity` must be in (0, 1] (got 1.1). `sensitivity` is a diagnostic probability and cannot be greater than 1.",
               fixed = TRUE)
  expect_error(design_detection(0.05, n_per_site = 10, icc = -0.1),
               "`icc` must be in [0, 1] (got -0.1). `icc` is a correlation and cannot be negative.", fixed = TRUE)
  expect_error(design_detection(0.05, n_per_site = 10, icc = 1.5),
               "`icc` must be in [0, 1] (got 1.5). `icc` is a correlation and cannot be greater than 1.", fixed = TRUE)
})

test_that("DET-V-6: icc above 0 with no sites, or both n_sites and n_per_site, is rejected", {
  expect_error(design_detection(0.05, icc = 0.05),
               "`icc` is greater than 0 but no sites were given (`n_sites` and `n_per_site` are both `NULL`), so the design effect cannot be computed. Supply `n_sites` (fix the number of sites) or `n_per_site` (fix the number of people per site), or set `icc = 0` for a design without clustering.",
               fixed = TRUE)
  expect_error(design_detection(0.05, n_sites = 10, n_per_site = 10, icc = 0.05),
               "Supply at most one of `n_sites` or `n_per_site`, not both. `n_sites` fixes the number of sites and solves for people per site; `n_per_site` fixes the number of people per site and solves for the number of sites.", fixed = TRUE)
})

test_that("DET-V-7: n, n_sites, n_per_site and fpc_N must be single positive integers, naming the problem", {
  for (arg in c("n", "n_sites", "n_per_site")) {
    call_with <- function(v) {
      a <- list(0.05); if (arg != "n") a$icc <- 0.05; a[[arg]] <- v
      do.call(design_detection, a)
    }
    expect_error(call_with("10"),      paste0("`", arg, "` must be a single finite positive integer (got class `character`)."), fixed = TRUE)
    expect_error(call_with(TRUE),      paste0("`", arg, "` must be a single finite positive integer (got class `logical`)."), fixed = TRUE)
    expect_error(call_with(c(10, 20)), paste0("`", arg, "` must be a single finite positive integer (got length = 2)."), fixed = TRUE)
    for (bad in c(0, 2.5, Inf)) {
      expect_error(call_with(bad), paste0("`", arg, "` must be a single finite positive integer (got ", bad, ")."), fixed = TRUE)
    }
  }
  expect_error(design_detection(0.05, fpc_N = "500"),
               "`fpc_N` must be a single finite positive integer (got class `character`). `fpc_N` is the population size. Set `fpc_N = NULL` to treat the population as very large.",
               fixed = TRUE)
  for (bad in c(0, 100.5, Inf)) {
    expect_error(design_detection(0.05, fpc_N = bad),
                 paste0("`fpc_N` must be a single finite positive integer (got ", bad, ")."), fixed = TRUE)
  }
  # The whole messages, with what each argument is
  expect_error(design_detection(0.05, n = 2.5),
               "`n` must be a single finite positive integer (got 2.5). `n` is the number of people tested.", fixed = TRUE)
  expect_error(design_detection(0.05, n_sites = 2.5, icc = 0.05),
               "`n_sites` must be a single finite positive integer (got 2.5). `n_sites` is the number of sites.", fixed = TRUE)
  expect_error(design_detection(0.05, n_per_site = 2.5, icc = 0.05),
               "`n_per_site` must be a single finite positive integer (got 2.5). `n_per_site` is the fixed number of people sampled per site.", fixed = TRUE)
  # A value within rounding of a whole number gets a hint to use round()
  expect_error(design_detection(0.05, n = 100 * 1.1),
               "`n` must be a single finite positive integer (got 110, which is not exactly an integer; use round(n)). `n` is the number of people tested.", fixed = TRUE)
  expect_error(design_detection(0.05, n_sites = 100 * 1.1, icc = 0.05),
               "`n_sites` must be a single finite positive integer (got 110, which is not exactly an integer; use round(n_sites)). `n_sites` is the number of sites.", fixed = TRUE)
  expect_error(design_detection(0.05, n_per_site = 100 * 1.1, icc = 0.05),
               "`n_per_site` must be a single finite positive integer (got 110, which is not exactly an integer; use round(n_per_site)). `n_per_site` is the fixed number of people sampled per site.", fixed = TRUE)
  expect_error(design_detection(0.05, fpc_N = 100 * 1.1),
               "`fpc_N` must be a single finite positive integer (got 110, which is not exactly an integer; use round(fpc_N)). `fpc_N` is the population size. Set `fpc_N = NULL` to treat the population as very large.", fixed = TRUE)
  # A real fraction gets no hint, however large; numbers are printed in full
  expect_error(design_detection(0.05, fpc_N = 1e10 + 0.5),
               "`fpc_N` must be a single finite positive integer (got 10000000000.5).", fixed = TRUE)
  old <- options(scipen = -9); on.exit(options(old), add = TRUE)
  expect_error(design_detection(0.05, n = 100000.5), "(got 100000.5)", fixed = TRUE)
  options(old)
  # Above 2^53 (the largest whole number R stores exactly) is rejected
  expect_error(design_detection(0.05, n = 2^53 + 2),
               "`n` (9007199254740994) is larger than 2^53 = 9007199254740992, the largest whole number R can store exactly. No real study tests this many people.", fixed = TRUE)
  expect_error(design_detection(0.05, n_sites = 1e200, icc = 0.05),
               "`n_sites` (1e+200) is larger than 2^53 = 9007199254740992, the largest whole number R can store exactly. No real study has this many sites.", fixed = TRUE)
  expect_error(design_detection(0.05, n_per_site = 2^53 + 2, icc = 0.05),
               "`n_per_site` (9007199254740994) is larger than 2^53 = 9007199254740992, the largest whole number R can store exactly. No real site has this many people.", fixed = TRUE)
  expect_error(design_detection(0.05, fpc_N = 2^53 + 2),
               "`fpc_N` (9007199254740994) is larger than 2^53 = 9007199254740992, the largest whole number R can store exactly. No real population is this large: set `fpc_N = NULL`.", fixed = TRUE)
  # Valid values: n = 139 (DET-2), n_sites = 10 (DET-C-2), n_per_site = 10
  # (DET-C-1), fpc_N = 200 (DET-F-1)
})

test_that("DET-V-8: more sites, or more people per site, than the population is rejected", {
  expect_error(design_detection(0.05, n_sites = 100, icc = 0.05, fpc_N = 50),
               "`n_sites` (100) exceeds `fpc_N` (50): you cannot have more sites than people in the population.",
               fixed = TRUE)
  expect_error(design_detection(0.05, n_per_site = 100, icc = 0.05, fpc_N = 50),
               "`n_per_site` (100) exceeds `fpc_N` (50): a site cannot have more people than there are in the whole population.",
               fixed = TRUE)
})

test_that("DET-V-9: in reverse, n must fit the sites and the population", {
  expect_error(design_detection(0.05, n = 5, n_sites = 10, icc = 0.05),
               "`n` (5) is smaller than `n_sites` (10): there would be fewer than one person per site. Increase `n`, or lower `n_sites`.",
               fixed = TRUE)
  expect_error(design_detection(0.05, n = 300, fpc_N = 200),
               "`n` (300) is larger than `fpc_N` (200): you cannot test more people than there are in the population.",
               fixed = TRUE)
  # Fewer people than one full site is fine with n_per_site (forward mode
  # can return it), and so is fewer people than sites when icc is 0
  expect_no_error(design_detection(0.05, n = 5, n_per_site = 10, icc = 0.05))
  expect_no_error(design_detection(0.05, n = 5, n_sites = 10))
})


# ---------------------------------------------------------------------------
# Return value (DET-R-)
# ---------------------------------------------------------------------------

test_that("DET-R-1: return list contains all expected fields, in order, in both directions", {
  fields <- c("n", "n_eff", "n_sites", "n_per_site", "detection_prob", "prevalence",
              "q", "n_cases", "sensitivity", "icc", "deff", "fpc_N", "mode")
  expect_named(design_detection(0.05), fields, ignore.order = FALSE)
  expect_named(design_detection(0.05, n = 32), fields, ignore.order = FALSE)
})

test_that("DET-R-2: inputs are echoed back as supplied", {
  res <- design_detection(0.05, 0.9, sensitivity = 0.9, n_per_site = 10, icc = 0.05)
  expect_equal(res[c("detection_prob", "prevalence", "sensitivity", "icc", "fpc_N")],
               list(detection_prob = 0.9, prevalence = 0.05, sensitivity = 0.9,
                    icc = 0.05, fpc_N = NULL))
  expect_equal(design_detection(0.05, fpc_N = 5000)$fpc_N, 5000)
  expect_equal(design_detection(0.05, n = 139)$n, 139)
})
