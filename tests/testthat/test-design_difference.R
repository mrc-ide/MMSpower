# Hand-checked reference cases for design_difference()
#
#   p_app = p * Se + (1 - p) * (1 - Sp)              (measured prevalence)
#   pbar  = (p1_app + p2_app) / 2
#   n0    = ((z_a * sqrt(2 * pbar * (1 - pbar)) +
#             z_b * sqrt(p1_app * (1 - p1_app) + p2_app * (1 - p2_app)))
#            / |p1_app - p2_app|)^2                   (people per group)
#   z_a   = qnorm(0.975) = 1.959964 (two-sided, 95%); z_b = qnorm(0.80) = 0.841621
#
# DD-1  (10 versus 20 percent, sensitivity = specificity = 1):
#   pbar = 0.15; sqrt(2 * 0.15 * 0.85) = 0.504975; sqrt(0.09 + 0.16) = 0.5
#   n0 = ((1.959964 * 0.504975 + 0.841621 * 0.5) / 0.1)^2 = 198.96 -> 199
#
# DD-C-1 (n_per_site = 20, icc = 0.05):
#   deff = 1 + (20 - 1) * 0.05 = 1.95; n = 198.96 * 1.95 = 387.98 -> 388,
#   20 sites (388 / 20 = 19.4, rounded up)
#
# DD-C-2 (n_sites = 30, icc = 0.05):
#   n = 198.96 * 30 * 0.95 / (30 - 198.96 * 0.05) = 282.79 -> 283;
#   283 / 30 = 9.43 -> 10 per site; deff = 1 + (10 - 1) * 0.05 = 1.45


# ---------------------------------------------------------------------------
# Core scenarios (DD-)
# ---------------------------------------------------------------------------

test_that("DD-1: 10 versus 20 percent with 80 percent power gives the hand-checked 199 per group", {
  res <- design_difference(0.10, 0.20)
  expect_equal(res$n_per_group, 199)
  expect_equal(res$n,           398)
  expect_equal(res$n_eff,       199)
  expect_equal(res$deff,        1)
  expect_null(res$n_sites)
  expect_null(res$n_per_site)
})

test_that("DD-2: matches R's own power.prop.test for several designs", {
  for (s in list(c(0.10, 0.20, 0.80), c(0.10, 0.20, 0.90), c(0.30, 0.50, 0.80),
                 c(0.05, 0.10, 0.80), c(0.60, 0.40, 0.80), c(0.01, 0.02, 0.80))) {
    ref <- ceiling(stats::power.prop.test(p1 = s[1], p2 = s[2], power = s[3])$n)
    expect_equal(design_difference(s[1], s[2], power = s[3])$n_per_group, ref)
  }
  # One-sided: 157 per group
  ref <- ceiling(stats::power.prop.test(p1 = 0.10, p2 = 0.20, power = 0.80,
                                        alternative = "one.sided")$n)
  expect_equal(design_difference(0.10, 0.20, alternative = "less")$n_per_group, ref)
  expect_equal(ref, 157)
})

test_that("DD-3: more power, or a higher conf_level, needs more people", {
  expect_equal(design_difference(0.10, 0.20, power = 0.90)$n_per_group,      266)
  expect_equal(design_difference(0.10, 0.20, conf_level = 0.99)$n_per_group, 297)
})

test_that("DD-4: a one-sided test needs fewer people than a two-sided one", {
  # 157 (one-sided, DD-2) versus 199 (two-sided, DD-1); "greater" and
  # "less" give the same n for the mirrored prevalences
  expect_equal(design_difference(0.20, 0.10, alternative = "greater")$n_per_group, 157)
  expect_equal(design_difference(0.10, 0.20, alternative = "less")$n_per_group,    157)
})

test_that("DD-5: a two-sided design does not depend on which group is listed first", {
  expect_equal(design_difference(0.20, 0.10)$n_per_group,
               design_difference(0.10, 0.20)$n_per_group)
})

test_that("DD-6: an imperfect test shrinks the measured difference and needs more people", {
  # p1_app = 0.1 * 0.9 + 0.9 * 0.05 = 0.135; p2_app = 0.2 * 0.9 + 0.8 * 0.05 = 0.22
  # difference 0.085 = 0.1 * (0.9 + 0.95 - 1); pbar = 0.1775
  # n0 = ((1.959964 * 0.540359 +0.841621 * 0.537006) / 0.085)^2 = 316.0 -> 317
  res <- design_difference(0.10, 0.20, sensitivity = 0.9, specificity = 0.95)
  expect_equal(res$apparent_prev1, 0.135)
  expect_equal(res$apparent_prev2, 0.22)
  expect_equal(res$difference,     -0.1)
  expect_equal(res$difference_app, -0.085)
  expect_equal(res$n_per_group,    317)
})

test_that("DD-7: sensitivity = specificity = 0.999 gives almost the same n as 1", {
  # 201 versus 199 (DD-1)
  expect_equal(design_difference(0.10, 0.20, sensitivity = 0.999,
                                 specificity = 0.999)$n_per_group, 201)
})

test_that("DD-8: a prevalence of 0 or 1 is allowed, as long as the two differ", {
  # For example, 0 versus 5 percent: has the condition appeared? Matches
  # R's power.prop.test, including at the edges.
  for (s in list(c(0, 0.05), c(0, 0.20), c(0.95, 1), c(0, 1), c(1, 0.5))) {
    ref <- ceiling(stats::power.prop.test(p1 = s[1], p2 = s[2], power = 0.80)$n)
    expect_equal(design_difference(s[1], s[2])$n_per_group, ref)
  }
  expect_equal(design_difference(0, 0.05)$n_per_group, 152)
  # 0 versus 1, one-sided at conf_level 0.5: the formula gives 0, but a group
  # always needs at least 1 person
  res <- design_difference(0, 1, conf_level = 0.5, alternative = "less")
  expect_equal(res$n_per_group, 1)
  expect_equal(res$n_eff, 1)
  # Equal prevalences of 0 are still "no difference" (see DD-V-4)
  expect_error(design_difference(0, 0), "are equal, or nearly so", fixed = TRUE)
})

test_that("DD-9: integers stored as R integers (30L) give the same result", {
  expect_equal(design_difference(0.10, 0.20, n_sites = 30L, icc = 0.05),
               design_difference(0.10, 0.20, n_sites = 30,  icc = 0.05))
  expect_equal(design_difference(0.10, 0.20, n_per_site = 20L, icc = 0.05),
               design_difference(0.10, 0.20, n_per_site = 20,  icc = 0.05))
  expect_equal(design_difference(0.10, 0.20, fpc_N = 2000L),
               design_difference(0.10, 0.20, fpc_N = 2000))
  # A 1 x 1 matrix, a 1-d array or a named number works like a plain number,
  # and nothing of it leaks into the result
  expect_equal(design_difference(0.10, 0.20, n_per_site = matrix(30), icc = 0.02, fpc_N = array(1000)),
               design_difference(0.10, 0.20, n_per_site = 30, icc = 0.02, fpc_N = 1000))
  expect_identical(design_difference(0.10, 0.20, n_sites = c(s = 30), icc = 0.05)$deff, 1.45)
})

test_that("DD-10: a tiny difference gives a huge n that is still a finite whole number", {
  # 50 versus 50.001 percent: about 39 billion per group. The formula gives
  # 39244398666.9993, just below a whole number, so 39244398667: rounding
  # must not drop any whole people
  res <- design_difference(0.5, 0.50001)
  expect_equal(res$n_per_group, 39244398667, tolerance = 0)
  expect_equal(res$n,           2 * 39244398667, tolerance = 0)
  expect_equal(res$n_per_group %% 1, 0)
  # The same with more sites than people needed (1 person per site)
  expect_warning(res <- design_difference(0.5, 0.50001, n_sites = 1e11, icc = 0.05),
                 "the study needs 39244398667 people per group", fixed = TRUE)
  expect_equal(res$n_per_group, 39244398667, tolerance = 0)
  # icc just above the zero cutoff: deff is only a hair below 1, which must
  # not lower n either
  expect_warning(res <- design_difference(0.5, 0.50001, n_sites = 1e11, icc = 2e-8),
                 "the study needs 39244398667 people per group", fixed = TRUE)
  expect_equal(res$n_per_group, 39244398667, tolerance = 0)
  # Too large for R to compute: an error, not Inf. Only a conf_level this
  # close to 1 can do it (z is then infinite), in every design
  for (a in list(list(), list(n_sites = 10, icc = 0.05), list(n_sites = 10, icc = 0.05, fpc_N = 1000),
                 list(n_per_site = 10, icc = 0.05)))
    expect_error(do.call(design_difference, c(list(0.1, 0.2, conf_level = 1 - 1e-16), a)),
                 "The required sample size is too large for R to compute. Use a lower `conf_level`.",
                 fixed = TRUE)
  # Above a trillion people per group the answer cannot be exact (rounding
  # error in the calculation is then more than one person), so it stops:
  # 50 versus 50.0001 percent needs about 3.9 trillion
  expect_error(design_difference(0.5, 0.5 + 1e-6),
               "The required sample size is more than a trillion (1000000000000) people per group, too large to work out to the nearest person. Use prevalences further apart, or a lower `power` or `conf_level`.",
               fixed = TRUE)
  # The advice fits the design: with large sites, a smaller n_per_site;
  # with an imperfect test, a more accurate one
  expect_error(design_difference(0.1, 0.2, n_per_site = 1e13, icc = 0.5),
               "Use prevalences further apart, a lower `power` or `conf_level`, or a smaller `n_per_site`.",
               fixed = TRUE)
  expect_error(design_difference(0.5, 0.5 + 1e-5, sensitivity = 0.55, specificity = 0.55),
               "Use prevalences further apart, a lower `power` or `conf_level`, or a more accurate test.",
               fixed = TRUE)
  # With a fixed number of sites just enough to reach the target, more
  # sites or a lower ICC (icc a hair below n_sites / n0, here 10 / 198.96)
  n0 <- ((qnorm(0.975) * sqrt(2 * 0.15 * 0.85) + qnorm(0.8) * 0.5) / 0.1)^2
  expect_error(design_difference(0.1, 0.2, n_sites = 10, icc = 10 / n0 * (1 - 1e-12)),
               "Use prevalences further apart, a lower `power` or `conf_level`, or more sites or a lower ICC.",
               fixed = TRUE)
  # Just under a trillion is still given; and a population brings it down
  expect_equal(design_difference(0.5, 0.5 + 3e-6)$n_per_group, 436048874133, tolerance = 0)
  expect_equal(design_difference(0.5, 0.5 + 1e-6, fpc_N = 1e6)$n_per_group, 1e6)
})


# ---------------------------------------------------------------------------
# Clustering (DD-C-)
# ---------------------------------------------------------------------------

test_that("DD-C-1: a fixed number of people per site multiplies n by the design effect", {
  res <- design_difference(0.10, 0.20, n_per_site = 20, icc = 0.05)
  expect_equal(res$deff,        1.95)
  expect_equal(res$n_per_group, 388)
  expect_equal(res$n_sites,     20)
  expect_equal(res$n_per_site,  20)
  expect_equal(res$n_eff,       199)   # the unclustered n, unchanged
  # 1 person per site: no two people share a site, so deff = 1 whatever the icc
  res1 <- design_difference(0.10, 0.20, n_per_site = 1, icc = 0.05)
  expect_equal(res1$deff,        1)
  expect_equal(res1$n_per_group, 199)
  expect_equal(res1$n_sites,     199)
})

test_that("DD-C-2: a fixed number of sites solves for people per site, and deff matches it", {
  res <- design_difference(0.10, 0.20, n_sites = 30, icc = 0.05)
  expect_equal(res$n_per_group, 283)
  expect_equal(res$n_sites,     30)
  expect_equal(res$n_per_site,  10)
  # deff is worked out again from the rounded 10 per site
  expect_equal(res$deff, 1 + (10 - 1) * 0.05)
})

test_that("DD-C-3: too few sites for the target power is rejected, giving the maximum power", {
  # 198.96 * 0.05 = 9.95 > 9 sites, so 80 percent cannot be reached
  expect_error(
    design_difference(0.10, 0.20, n_sites = 9, icc = 0.05),
    "Target power of 80 percent cannot be reached with 9 sites per group and ICC = 0.05. The highest power these sites can give is 75.9226 percent. Increase `n_sites`, lower the ICC, or lower the target power.",
    fixed = TRUE)
  # One site says "site", and a small ICC is shown as given, not rounded to 0
  expect_error(
    design_difference(0.10, 0.12, n_sites = 1, icc = 0.0004),
    "Target power of 80 percent cannot be reached with 1 site per group and ICC = 0.0004. The highest power this site can give is 61.79543 percent.",
    fixed = TRUE)
  # The maximum power is printed like the target, so one just below 80
  # (79.986) is not shown as 80.0
  expect_error(
    design_difference(0.10, 0.20, n_sites = 9, icc = 0.04525),
    "The highest power these sites can give is 79.98642 percent.",
    fixed = TRUE)
  # A maximum close to 0 (2e-16 percent) is not printed as 2.011586e-16
  expect_error(design_difference(0.001, 0.999, n_sites = 1, icc = 1),
               "The highest power this site can give is less than 0.0001 percent.", fixed = TRUE)
  # Closer still: 7 digits would show 80, so more are shown (79.999999)
  expect_error(
    design_difference(0.10, 0.20, n_sites = 9, icc = 0.045234448280423757),
    "Target power of 80 percent cannot be reached with 9 sites per group and ICC = 0.0452344482804238. The highest power these sites can give is 79.999999 percent.",
    fixed = TRUE)
  # The highest power 9 sites can give is 75.9 percent: asking for 75.5
  # works and asking for 76.0 does not
  expect_error(design_difference(0.10, 0.20, power = 0.760, n_sites = 9, icc = 0.05),
               "cannot be reached", fixed = TRUE)
  expect_no_error(design_difference(0.10, 0.20, power = 0.755, n_sites = 9, icc = 0.05))
  # With an imperfect test the maximum power uses the measured prevalences
  # (0.135 and 0.22, see DD-6): 12 sites are worth at most 12 / 0.05 = 240
  # people, which gives 68.43652 percent power
  p1 <- 0.135; p2 <- 0.22; pbar <- (p1 + p2) / 2
  max_power <- pnorm((0.085 * sqrt(12 / 0.05) - qnorm(0.975) * sqrt(2 * pbar * (1 - pbar))) /
                       sqrt(p1 * (1 - p1) + p2 * (1 - p2)))
  expect_equal(format(100 * max_power), "68.43652")
  expect_error(
    design_difference(0.10, 0.20, sensitivity = 0.9, specificity = 0.95, n_sites = 12, icc = 0.05),
    "Target power of 80 percent cannot be reached with 12 sites per group and ICC = 0.05. The highest power these sites can give is 68.43652 percent.",
    fixed = TRUE)
})

test_that("DD-C-4: more sites than people needed uses only as many sites as people, 1 each, with a warning", {
  # The study needs 199 people per group (198.96 rounded up). 250 sites is
  # more than that, so 199 of them are used with 1 person at each.
  expect_warning(
    res <- design_difference(0.10, 0.20, n_sites = 250, icc = 0.05),
    "With these prevalences and this power, the study needs 199 people per group, fewer than the 250 sites given (`n_sites`). The design uses 199 sites with 1 person at each.",
    fixed = TRUE)
  expect_equal(res$n_per_group, 199)
  expect_equal(res$n_sites,     199)
  expect_equal(res$n_per_site,  1)
  expect_equal(res$deff,        1)
  # The same with icc = 0: 199 sites are used, not all 250
  expect_warning(res0 <- design_difference(0.10, 0.20, n_sites = 250), "uses 199 sites", fixed = TRUE)
  expect_equal(res0$n_sites, 199)
  # Large site numbers are printed in full, not as 1e+05
  expect_warning(design_difference(0.10, 0.20, n_sites = 1e5, icc = 0.05),
                 "fewer than the 100000 sites given", fixed = TRUE)
  # 199 sites is exactly enough: 1 person each, no warning
  expect_no_warning(res <- design_difference(0.10, 0.20, n_sites = 199, icc = 0.05))
  expect_equal(res$n_per_group, 199)
  expect_equal(res$n_per_site,  1)
  expect_equal(res$deff,        1)
  # One site fewer than that: 198 sites need 200 people, so 2 per site and
  # 396 collected, nearly double n_per_group (see the help page for n)
  res <- design_difference(0.10, 0.20, n_sites = 198, icc = 0.05)
  expect_equal(c(res$n_per_group, res$n_per_site), c(200, 2))
  # A study needing 1 person says "1 person" and "1 site"
  expect_warning(design_difference(0, 1, conf_level = 0.5, alternative = "less", n_sites = 5, icc = 0.05),
                 "the study needs 1 person per group, fewer than the 5 sites given (`n_sites`). The design uses 1 site with 1 person at each.",
                 fixed = TRUE)
  # Power chosen so that exactly 150 / 230 / 300 / 500 people are needed with
  # no clustering: that many sites gives 1 person per site and deff = 1, a
  # real design. Rounding puts deff a hair below 1 for 500, and n a hair
  # above 300 for 300 (which must not round up to 301 and 2 per site).
  for (k in c(150, 230, 300, 500)) {
    z_a <- qnorm(0.975)
    se0 <- sqrt(2 * 0.15 * 0.85)
    se1 <- sqrt(0.10 * 0.90 + 0.20 * 0.80)
    pw  <- pnorm((0.10 * sqrt(k) - z_a * se0) / se1)
    res <- design_difference(0.10, 0.20, power = pw, n_sites = k, icc = 0.05)
    expect_equal(res$n_per_group, k)
    expect_equal(res$n_eff,       k)
    expect_equal(res$n_per_site,  1)
    expect_equal(res$deff,        1)
    # 50 more sites than needed: still exactly k people, with a warning
    expect_warning(res <- design_difference(0.10, 0.20, power = pw, n_sites = k + 50, icc = 0.05),
                   paste0("the study needs ", k, " people per group"), fixed = TRUE)
    expect_equal(res$n_per_group, k)
    expect_equal(res$n_sites,     k)
  }
})

test_that("DD-C-5: icc = 0 with sites is the unclustered result spread over the sites; a tiny icc is used as given", {
  res <- design_difference(0.10, 0.20, n_sites = 30, icc = 0)
  expect_equal(res$n_per_group, 199)
  expect_equal(res$n_per_site,  7)    # 199 / 30 = 6.6, rounded up
  expect_identical(res$deff, 1)
  # A tiny icc with sites is used, not set to 0: here it changes deff only
  # in the 12th decimal place (1 + 6 * 1e-12), and n not at all
  res <- design_difference(0.10, 0.20, n_sites = 30, icc = 1e-12)
  expect_equal(c(res$n_per_group, res$n_per_site), c(199, 7))
  expect_identical(res$deff, 1 + 6 * 1e-12)
  expect_identical(res$icc,  1e-12)   # echoed as given
  # With a very large site even a tiny icc matters: 15.7 million people at
  # one site with icc 1.4e-8 have deff 1.28, so more people are needed
  res <- design_difference(0.5, 0.5005, n_sites = 1, icc = 1.4e-8)
  expect_equal(res$n_per_group, 20119357)
  expect_equal(res$deff, 1 + (20119357 - 1) * 1.4e-8)
  # A tiny icc without sites is not an error: it counts as 0
  expect_equal(design_difference(0.10, 0.20, icc = 1e-12)$n_per_group, 199)
  # One site with icc 0: all 199 at that site
  res1 <- design_difference(0.10, 0.20, n_sites = 1)
  expect_equal(res1$n_per_group, 199)
  expect_equal(res1$n_sites,     1)
  expect_equal(res1$n_per_site,  199)
  # n_per_site with icc 0: the same 199, on 10 sites of 20 (199 / 20 = 9.95)
  res2 <- design_difference(0.10, 0.20, n_per_site = 20)
  expect_equal(res2$n_per_group, 199)
  expect_equal(res2$deff,        1)
  expect_equal(res2$n_sites,     10)
})

test_that("DD-C-6: a larger icc needs more people", {
  expect_gt(design_difference(0.10, 0.20, n_per_site = 20, icc = 0.10)$n_per_group,
            design_difference(0.10, 0.20, n_per_site = 20, icc = 0.05)$n_per_group)
})


# ---------------------------------------------------------------------------
# Finite population correction (DD-F-)
# ---------------------------------------------------------------------------

test_that("DD-F-1: a small population reduces n, and leaves n_eff alone", {
  # 198.96 * 2000 / (198.96 + 2000 - 1) = 181.04 -> 182
  res <- design_difference(0.10, 0.20, fpc_N = 2000)
  expect_equal(res$n_per_group, 182)
  expect_equal(res$n_eff,       199)
  # With n_per_site: 387.98 * 1000 / (387.98 + 999) = 279.73 -> 280, 14 sites
  res2 <- design_difference(0.10, 0.20, n_per_site = 20, icc = 0.05, fpc_N = 1000)
  expect_equal(res2$n_per_group, 280)
  expect_equal(res2$n_sites,     14)
})

test_that("DD-F-2: with n_sites fixed, the FPC and the design effect are solved together", {
  # n solves n0 * (1 + (n / n_sites - 1) * icc) * (N - n) / (N - 1) = n;
  # checked against a numerical solution of that equation (uniroot)
  n0 <- ((stats::qnorm(0.975) * sqrt(2 * 0.15 * 0.85) +
            stats::qnorm(0.80) * sqrt(0.09 + 0.16)) / 0.1)^2
  for (case in list(c(sites = 30, icc = 0.05, N = 1000), c(sites = 10, icc = 0.05, N = 2000),
                    c(sites = 50, icc = 0.10, N = 5000), c(sites = 9, icc = 0.05, N = 2000),
                    c(sites = 150, icc = 1, N = 2000))) {
    k <- case[["sites"]]; r <- case[["icc"]]; N <- case[["N"]]
    eq   <- function(n) n0 * (1 + (n / k - 1) * r) * (N - n) / (N - 1) - n
    root <- stats::uniroot(eq, c(1e-9, N), tol = 1e-10)$root
    # (The 9-site and icc = 1 cases also warn: see DD-F-3.)
    res  <- suppressWarnings(design_difference(0.10, 0.20, n_sites = k, icc = r, fpc_N = N))
    expect_equal(res$n_per_group, ceiling(root))
  }
  # A very large population gives the same n as no population (283)
  expect_equal(design_difference(0.10, 0.20, n_sites = 30, icc = 0.05, fpc_N = 1e9)$n_per_group, 283)
  # With icc = 1 one root of the equation is 0. Here the other is 0 too
  # (n0 = 1.6, 2 sites, population 5: 1.6 * (n / 2) * (5 - n) / 4 = n only
  # at n = 0), so the unclustered answer is used: 2 people, 1 per site
  res <- design_difference(1, 0, alternative = "greater", conf_level = 0.96318086493984867,
                           n_sites = 2, icc = 1, fpc_N = 5)
  expect_equal(c(res$n_per_group, res$n_sites, res$n_per_site), c(2, 2, 1))
})

test_that("DD-F-3: a target unachievable with n_sites alone becomes achievable with a population, with a warning", {
  # 10 sites, icc 0.05: 198.96 * 0.05 = 9.95 is just under 10, so the target is
  # only just achievable and needs 36469 per group. With a population of
  # 2000 the joint solve gives 525 (53 per site, deff 1 + 52 * 0.05 = 3.6).
  expect_equal(design_difference(0.10, 0.20, n_sites = 10, icc = 0.05)$n_per_group, 36469)
  res <- design_difference(0.10, 0.20, n_sites = 10, icc = 0.05, fpc_N = 2000)
  expect_equal(res$n_per_group, 525)
  expect_equal(res$n_per_site,  53)
  expect_equal(res$deff,        3.6)
  # 10 sites can reach the target on their own, so no warning
  expect_no_warning(design_difference(0.10, 0.20, n_sites = 10, icc = 0.05, fpc_N = 2000))
  # 9 sites: unachievable without a population (DD-C-3). With a population
  # of 2000 it works, but only by testing 596 of the 2000, so it warns
  expect_warning(
    res <- design_difference(0.10, 0.20, n_sites = 9, icc = 0.05, fpc_N = 2000),
    "With 9 sites per group and ICC = 0.05, the target power can only be reached because of the finite-population correction: the design tests 596 of the 2000 people in each group (67 per site). Check that `fpc_N` is right, or use more sites.",
    fixed = TRUE)
  expect_equal(c(res$n_per_group, res$n_per_site), c(596, 67))
  # A huge population can give a huge design: 95312682 people per group
  expect_warning(design_difference(0.10, 0.20, n_sites = 9, icc = 0.05, fpc_N = 1e9),
                 "the design tests 95312682 of the 1000000000 people in each group (10590298 per site)",
                 fixed = TRUE)
  # One site says "site"
  expect_warning(design_difference(0.10, 0.20, n_sites = 1, icc = 0.5, fpc_N = 500),
                 "With 1 site per group and ICC = 0.5,", fixed = TRUE)
  # When some sites are dropped as too many (the other warning), this one
  # does not also say "use more sites"
  w <- character(0)
  withCallingHandlers(design_difference(0.1, 0.5, n_sites = 10, icc = 1, fpc_N = 10),
                      warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_length(w, 2)
  expect_true(grepl("The design uses 7 sites with 1 person at each.", w[1], fixed = TRUE))
  expect_true(endsWith(w[2], "the design tests 7 of the 10 people in each group (1 per site). Check that `fpc_N` is right."))
})

test_that("DD-F-4: a rounded design larger than the population warns; too many sites is judged after the population correction", {
  # Population 60, 23 sites: n = 47, 47 / 23 = 2.04 -> 3 per site,
  # so 23 x 3 = 69 people per group, more than the 60 who exist
  expect_warning(
    res <- design_difference(0.10, 0.20, n_sites = 23, icc = 0.05, fpc_N = 60),
    "The rounded design (23 sites x 3 = 69 per group) exceeds `fpc_N` = 60: the required n of 47 per group is close to a census of the population. Consider surveying the whole population instead.",
    fixed = TRUE)
  expect_equal(res$n_per_group, 47)
  # Large numbers are printed in full, not as 3e+05 or 1e+06
  expect_warning(
    design_difference(0.5, 0.5005, n_sites = 300000, icc = 1e-7, fpc_N = 1e6),
    "The rounded design (300000 sites x 4 = 1200000 per group) exceeds `fpc_N` = 1000000: the required n of 940112 per group",
    fixed = TRUE)
  expect_warning(design_difference(0.5, 0.50001, n_sites = 3, fpc_N = 1e5),
                 "the required n of 100000 per group", fixed = TRUE)
  # Too many sites is judged against the population-corrected n (166.08, so
  # 167 whole people): 167 sites is exactly enough, 250 is cut to 167
  expect_no_warning(res <- design_difference(0.10, 0.20, n_sites = 167, icc = 0.05, fpc_N = 1000))
  expect_equal(res$n_per_site, 1)
  expect_warning(res <- design_difference(0.10, 0.20, n_sites = 250, icc = 0.05, fpc_N = 1000),
                 "the study needs 167 people per group, fewer than the 250 sites given", fixed = TRUE)
  expect_equal(res$n_sites,    167)
  expect_equal(res$n_per_site, 1)
})

test_that("DD-F-5: a population of 1 or 2 per group means testing everyone", {
  for (N in c(1, 2)) {
    res <- design_difference(0.10, 0.20, fpc_N = N)
    expect_equal(res$n_per_group, N)
    expect_equal(res$n_eff,       199)   # without the population, unchanged
  }
  # Never more people than the population, even when rounding error in the
  # correction comes out a hair above 1 (here 1.0000000000088)
  res <- suppressWarnings(design_difference(0.17, 0.163828627001494, power = 0.15060057245877,
                                            alternative = "greater", conf_level = 0.84940042754123,
                                            fpc_N = 1))
  expect_equal(res$n_per_group, 1)
  # 0 versus 1 needs 0 people by the formula (see DD-8): with a population
  # of 1 that is still 1 person, with or without a site
  expect_equal(design_difference(0, 1, conf_level = 0.5, alternative = "less", fpc_N = 1)$n_per_group, 1)
  res <- design_difference(0, 1, conf_level = 0.5, alternative = "less", n_sites = 1,
                           icc = 0.05, fpc_N = 1)
  expect_equal(c(res$n_per_group, res$n_sites, res$n_per_site), c(1, 1, 1))
})


# ---------------------------------------------------------------------------
# Validation (DD-V-)
# ---------------------------------------------------------------------------

test_that("DD-V-1: prevalences below 0 or above 1 are rejected", {
  # Valid values: 0.10 and 0.20 (DD-1); 0 and 1 themselves (DD-8)
  expect_error(design_difference(-0.1, 0.2),
               "`prevalence1` must be in [0, 1] (got -0.1). `prevalence1` is a fraction of the population and cannot be negative.",
               fixed = TRUE)
  expect_error(design_difference(0.1, 1.2),
               "`prevalence2` must be in [0, 1] (got 1.2). `prevalence2` is a fraction of the population and cannot be greater than 1.",
               fixed = TRUE)
  # Large values are printed in full, not as 1e+05
  expect_error(design_difference(1e5, 0.2), "`prevalence1` must be in [0, 1] (got 100000).", fixed = TRUE)
})

test_that("DD-V-2: power at or below the false-alarm rate is rejected; below 0.5 warns", {
  # Valid values: 0.80 (DD-1), 0.90 (DD-3); 0.5 itself works with no warning
  for (bad in c(1, 1.2)) {
    expect_error(design_difference(0.1, 0.2, power = bad),
                 paste0("`power` must be in (0, 1) (got ", bad, "). For example, use 0.80 for 80 percent power."),
                 fixed = TRUE)
  }
  # At 95 percent confidence the false-alarm rate is 0.05: 0.05 and below,
  # including 0 and negative values, are rejected
  for (bad in c(-0.2, 0, 0.03, 0.05)) {
    expect_error(design_difference(0.1, 0.2, power = bad),
                 paste0("`power` must be greater than 1 - `conf_level` = 0.05 (got ", bad, "). Even when there is ",
                        "no real difference, the test reports one 5 percent of the time, so a lower power is not meaningful."),
                 fixed = TRUE)
  }
  # The false-alarm rate follows conf_level: at 0.80 it is 0.2
  expect_error(design_difference(0.1, 0.2, power = 0.2, conf_level = 0.8),
               "`power` must be greater than 1 - `conf_level` = 0.2 (got 0.2).", fixed = TRUE)
  # Equal only up to rounding (1 - 0.8 = 0.19999999999999996) counts as equal,
  # but a power truly above the rate is allowed, however close
  expect_error(design_difference(0.1, 0.2, power = 1 - 0.8, conf_level = 0.8),
               "`power` must be greater than 1 - `conf_level` = 0.2", fixed = TRUE)
  expect_equal(suppressWarnings(design_difference(0.1, 0.2, power = 0.05 + 1e-8))$n_per_group, 3)
  expect_error(design_difference(0.1, 0.2, power = 1e5), "`power` must be in (0, 1) (got 100000).", fixed = TRUE)
  # Between the false-alarm rate and 0.5: warns, and the result is still correct
  # (53 per group really gives 30 percent power, see power.prop.test)
  expect_warning(res <- design_difference(0.1, 0.2, power = 0.3),
                 "`power` = 0.3 is below 0.5, so the study would be more likely to miss a real difference than to find it.",
                 fixed = TRUE)
  expect_equal(res$n_per_group, 53)
  expect_equal(res$n_per_group, ceiling(stats::power.prop.test(p1 = 0.1, p2 = 0.2, power = 0.3)$n))
  expect_warning(design_difference(0.1, 0.2, power = 0.45), "is below 0.5", fixed = TRUE)
  expect_no_warning(design_difference(0.1, 0.2, power = 0.5))
})

test_that("DD-V-3: alternative must be one of the three names", {
  expect_error(design_difference(0.1, 0.2, alternative = "bigger"),
               "`alternative` must be one of 'two.sided', 'greater', or 'less' (got 'bigger').",
               fixed = TRUE)
  expect_error(design_difference(0.1, 0.2, alternative = NA),
               "`alternative` must be one of 'two.sided', 'greater', or 'less' (got class `logical`, length 1).",
               fixed = TRUE)
  expect_error(design_difference(0.1, 0.2, alternative = c("less", "greater")),
               "(got class `character`, length 2).", fixed = TRUE)
})

test_that("DD-V-4: the prevalences must point the way the alternative says, and must differ", {
  expect_error(design_difference(0.1, 0.2, alternative = "greater"),
               "`prevalence1` (0.1) must be greater than `prevalence2` (0.2) when alternative = 'greater'. Use alternative = 'less'.",
               fixed = TRUE)
  expect_error(design_difference(0.2, 0.1, alternative = "less"),
               "`prevalence1` (0.2) must be less than `prevalence2` (0.1) when alternative = 'less'. Use alternative = 'greater'.",
               fixed = TRUE)
  # Equal prevalences get the "equal" message for every alternative, so the
  # direction message only suggests switching when switching would work
  for (alt in c("greater", "less"))
    expect_error(design_difference(0.1, 0.1, alternative = alt),
                 "`prevalence1` (0.1) and `prevalence2` (0.1) are equal, or nearly so", fixed = TRUE)
  expect_error(design_difference(0.2, 0.2 + 1e-12, alternative = "greater"),
               "`prevalence1` (0.2) and `prevalence2` (0.200000000001) are equal, or nearly so", fixed = TRUE)
  expect_error(design_difference(0.1, 0.1),
               "`prevalence1` (0.1) and `prevalence2` (0.1) are equal, or nearly so, so there is no difference to detect. Use two different prevalences.",
               fixed = TRUE)
  # Less than 1e-8 apart counts as equal, for every alternative
  expect_error(design_difference(0.1, 0.1 + 1e-9),
               "`prevalence1` (0.1) and `prevalence2` (0.100000001) are equal, or nearly so", fixed = TRUE)
  expect_error(design_difference(0.1, 0.1 - 1e-9, alternative = "greater"),
               "`prevalence1` (0.1) and `prevalence2` (0.099999999) are equal, or nearly so", fixed = TRUE)
  expect_error(design_difference(0.1, 0.1 + 1e-9, alternative = "less"),
               "are equal, or nearly so", fixed = TRUE)
  # An imperfect test can bring different prevalences almost together. The
  # two measured prevalences are shown with enough digits to differ
  expect_error(design_difference(0.1, 0.10000005, sensitivity = 0.55, specificity = 0.55),
               "The measured prevalences differ by only 0.000000005 (0.46 and 0.460000005), so the required sample size would be extremely large. Use prevalences further apart, or a more accurate test.",
               fixed = TRUE)
  # Very small numbers are printed in full too, not as 8.25e-09
  expect_error(design_difference(0, 1.5e-8, sensitivity = 0.55, specificity = 1),
               "differ by only 0.00000000825 (0 and 0.00000000825)", fixed = TRUE)
})

test_that("DD-V-5: a vector, text, TRUE/FALSE, a bare NA or NULL in a single-number argument is rejected, naming the class", {
  expect_error(design_difference(c(0.1, 0.2), 0.3),
               "`prevalence1` must be a single number in [0, 1] (got class `numeric`, length 2).", fixed = TRUE)
  expect_error(design_difference(TRUE, 0.2),
               "`prevalence1` must be a single number in [0, 1] (got class `logical`, length 1).", fixed = TRUE)
  expect_error(design_difference(0.1, "0.2"),
               "`prevalence2` must be a single number in [0, 1] (got class `character`, length 1).", fixed = TRUE)
  expect_error(design_difference(0.1, 0.2, power = TRUE),
               "`power` must be a single number in (0, 1) (got class `logical`, length 1).", fixed = TRUE)
  expect_error(design_difference(0.1, 0.2, sensitivity = NA),
               "`sensitivity` must be a single number in (0, 1] (got class `logical`, length 1). Note: `TRUE`, `FALSE` and a plain `NA` are logical, not numeric. Pass 1 if the test never misses a case.",
               fixed = TRUE)
  expect_error(design_difference(0.1, 0.2, specificity = NULL),
               "`specificity` must be a single number in (0, 1] (got class `NULL`, length 0). Note: `TRUE`, `FALSE` and a plain `NA` are logical, not numeric. Pass 1 if the test never gives a false positive.",
               fixed = TRUE)
  expect_error(design_difference(0.1, 0.2, specificity = TRUE),
               "`specificity` must be a single number in (0, 1] (got class `logical`, length 1).", fixed = TRUE)
  expect_error(design_difference(0.1, 0.2, conf_level = list(0.95)),
               "`conf_level` must be a single number in [0.5, 1) (got class `list`, length 1).", fixed = TRUE)
  expect_error(design_difference(0.1, 0.2, icc = FALSE),
               "`icc` must be a single number in [0, 1] (got class `logical`, length 1).", fixed = TRUE)
})

test_that("DD-V-6: NaN, Inf or -Inf in a single-number argument gives the finite-value error", {
  expect_error(design_difference(NaN, 0.2),  "`prevalence1` must be a single finite number (got NaN).", fixed = TRUE)
  expect_error(design_difference(0.1, Inf),  "`prevalence2` must be a single finite number (got Inf).", fixed = TRUE)
  expect_error(design_difference(0.1, 0.2, power = -Inf),
               "`power` must be a single finite number (got -Inf).", fixed = TRUE)
  expect_error(design_difference(0.1, 0.2, sensitivity = NaN),
               "`sensitivity` must be a single finite number (got NaN).", fixed = TRUE)
  expect_error(design_difference(0.1, 0.2, specificity = Inf),
               "`specificity` must be a single finite number (got Inf).", fixed = TRUE)
  expect_error(design_difference(0.1, 0.2, conf_level = NaN),
               "`conf_level` must be a single finite number (got NaN).", fixed = TRUE)
  expect_error(design_difference(0.1, 0.2, icc = Inf),
               "`icc` must be a single finite number (got Inf).", fixed = TRUE)
})

test_that("DD-V-7: sensitivity and specificity range checks give the right message for each bad value", {
  # Valid values: 1 (DD-1), 0.9 and 0.95 (DD-6)
  expect_error(design_difference(0.1, 0.2, sensitivity = 0),
               "`sensitivity` must be in (0, 1] (got 0). A sensitivity of 0 means everyone who truly has the condition tests negative.",
               fixed = TRUE)
  expect_error(design_difference(0.1, 0.2, sensitivity = 1.1),
               "`sensitivity` must be in (0, 1] (got 1.1). `sensitivity` is a diagnostic probability and cannot be greater than 1.",
               fixed = TRUE)
  expect_error(design_difference(0.1, 0.2, sensitivity = -0.5),
               "`sensitivity` must be in (0, 1] (got -0.5). `sensitivity` is a diagnostic probability and cannot be negative.",
               fixed = TRUE)
  expect_error(design_difference(0.1, 0.2, specificity = 1.1),
               "`specificity` must be in (0, 1] (got 1.1). `specificity` is a diagnostic probability and cannot be greater than 1.",
               fixed = TRUE)
  expect_error(design_difference(0.1, 0.2, specificity = 0),
               "`specificity` must be in (0, 1] (got 0). A specificity of 0 means everyone who truly does not have the condition tests positive.",
               fixed = TRUE)
  expect_error(design_difference(0.1, 0.2, specificity = -0.5),
               "`specificity` must be in (0, 1] (got -0.5). `specificity` is a diagnostic probability and cannot be negative.",
               fixed = TRUE)
})

test_that("DD-V-8: sensitivity + specificity of 1 or less is rejected, including 1 plus rounding", {
  expect_error(design_difference(0.1, 0.2, sensitivity = 0.5, specificity = 0.5),
               paste0("`sensitivity` + `specificity` must exceed 1 (got 0.5 + 0.5 = 1). The difference the test can ",
                      "see is the true difference times `sensitivity` + `specificity` - 1, so with a sum of 1 or less ",
                      "the test cannot see the difference, or sees it reversed."),
               fixed = TRUE)
  # A sum truly above 1, however little, is not refused here (it then fails
  # as a measured difference too small to detect, see DD-V-4)
  expect_error(design_difference(0.1, 0.2, sensitivity = 0.5, specificity = 0.50000001),
               "The measured prevalences differ by only 0.000000001", fixed = TRUE)
  # Small values are printed in full, not as 1e-05
  expect_error(design_difference(0.1, 0.2, sensitivity = 1e-5, specificity = 0.99),
               "(got 0.00001 + 0.99 = 0.99001)", fixed = TRUE)
  expect_error(design_difference(0.1, 0.2, sensitivity = 0.3, specificity = 0.3),
               "(got 0.3 + 0.3 = 0.6)", fixed = TRUE)
  # A calculated specificity of 0.7000000000000004 counts as 0.7
  expect_error(design_difference(0.1, 0.2, sensitivity = 0.3,
                                 specificity = 0.7 + 2 * .Machine$double.eps),
               "`sensitivity` + `specificity` must exceed 1 (got 0.3 + 0.7", fixed = TRUE)
})

test_that("DD-V-9: sensitivity + specificity above 1 but below 1.1 warns, showing the full sum", {
  expect_warning(
    res <- design_difference(0.1, 0.5, sensitivity = 0.55, specificity = 0.5),
    paste0("`sensitivity` + `specificity` = 1.05, which is very close to 1. The test then barely tells apart people ",
           "with and without the condition, so the required n is very large, and small changes in the assumed ",
           "sensitivity or specificity will change n a lot."),
    fixed = TRUE)
  expect_equal(res$n_per_group, 9802)
  # A sum of 1.0000001 warns, showing all the digits; the study would then
  # need more than a trillion people per group, which stops (see DD-10)
  expect_error(
    expect_warning(design_difference(0.1, 0.5, sensitivity = 0.5000001, specificity = 0.5),
                   "= 1.0000001, which is very close to 1.", fixed = TRUE),
    "more than a trillion", fixed = TRUE)
  # The warning comes only after every check has passed: not before an
  # unrelated error
  expect_no_warning(expect_error(
    design_difference(0.1, 0.2, sensitivity = 0.5, specificity = 0.55, icc = 2),
    "`icc` must be in [0, 1] (got 2)", fixed = TRUE))
  expect_no_warning(design_difference(0.1, 0.5, sensitivity = 0.61, specificity = 0.5))
  # The cutoff: 1.09 warns, 1.10 does not
  expect_warning(design_difference(0.1, 0.5, sensitivity = 0.59, specificity = 0.5),
                 "= 1.09, which is very close to 1.", fixed = TRUE)
  expect_no_warning(design_difference(0.1, 0.5, sensitivity = 0.6, specificity = 0.5))
  # A sum of 1.1 that rounding puts a hair below it (0.15 + 0.95) does not warn
  expect_no_warning(design_difference(0.1, 0.5, sensitivity = 0.15, specificity = 0.95))
})

test_that("DD-V-10: icc and conf_level outside their ranges are rejected", {
  expect_error(design_difference(0.1, 0.2, n_sites = 30, icc = -0.1),
               "`icc` must be in [0, 1] (got -0.1). `icc` is a correlation and cannot be negative.", fixed = TRUE)
  expect_error(design_difference(0.1, 0.2, n_sites = 30, icc = 1.5),
               "`icc` must be in [0, 1] (got 1.5). `icc` is a correlation and cannot be greater than 1.", fixed = TRUE)
  for (bad in c(1, 1.05)) {
    expect_error(design_difference(0.1, 0.2, conf_level = bad),
                 paste0("`conf_level` must be in [0.5, 1) (got ", bad, "). For example, use 0.95 for a 95% confidence interval."),
                 fixed = TRUE)
  }
  # Below 0.5 the one-sided formula gives nonsense n (0.3 -> 3, 0.1 -> 6);
  # 0 and negative values get the same message
  for (bad in c(-0.1, 0, 0.1, 0.3, 0.49)) {
    expect_error(design_difference(0.2, 0.1, conf_level = bad, alternative = "greater"),
                 paste0("`conf_level` must be at least 0.5 (got ", bad, "). Below 0.5, when there is ",
                        "no real difference, the test would wrongly report one more than half the time."),
                 fixed = TRUE)
  }
  expect_equal(design_difference(0.2, 0.1, conf_level = 0.5, alternative = "greater")$n_per_group, 18)
  # icc = 1 itself works: everyone at a site is alike, so deff = site size
  # (20) and n = 198.96 * 20 = 3979.3 -> 3980
  res <- design_difference(0.1, 0.2, n_per_site = 20, icc = 1)
  expect_equal(res$deff,        20)
  expect_equal(res$n_per_group, 3980)
  # With n_sites and icc = 1 each site is worth one person, so 199 sites of
  # 1 person is the only design that works: 198 cannot reach the target
  expect_equal(design_difference(0.1, 0.2, n_sites = 199, icc = 1)$n_per_site, 1)
  expect_error(design_difference(0.1, 0.2, n_sites = 198, icc = 1),
               "cannot be reached", fixed = TRUE)
  # Valid values: icc = 0.05 (DD-C-1), conf_level = 0.99 (DD-3)
})

test_that("DD-V-11: icc above 0 with no sites, or both n_sites and n_per_site, is rejected", {
  expect_error(design_difference(0.1, 0.2, icc = 0.05),
               paste0("`icc` is greater than 0 but no sites were given (`n_sites` and `n_per_site` are both `NULL`), ",
                      "so the design effect cannot be computed. Supply `n_sites` (fix the number of sites) or ",
                      "`n_per_site` (fix the number of people per site), or set `icc = 0` for a design without clustering."),
               fixed = TRUE)
  expect_error(design_difference(0.1, 0.2, n_sites = 10, n_per_site = 10, icc = 0.05),
               "Supply at most one of `n_sites` or `n_per_site`, not both.", fixed = TRUE)
})

test_that("DD-V-12: n_sites, n_per_site and fpc_N must be single positive integers, naming the problem", {
  for (arg in c("n_sites", "n_per_site")) {
    call_with <- function(v) {
      a <- list(0.1, 0.2, icc = 0.05); a[[arg]] <- v
      do.call(design_difference, a)
    }
    expect_error(call_with("10"),     paste0("`", arg, "` must be a single finite positive integer (got class `character`)."), fixed = TRUE)
    expect_error(call_with(TRUE),     paste0("`", arg, "` must be a single finite positive integer (got class `logical`)."), fixed = TRUE)
    expect_error(call_with(c(10, 20)), paste0("`", arg, "` must be a single finite positive integer (got length = 2)."), fixed = TRUE)
    for (bad in c(-5, 0, 2.5, Inf)) {
      expect_error(call_with(bad), paste0("`", arg, "` must be a single finite positive integer (got ", bad, ")."), fixed = TRUE)
    }
    # Large values are printed in full, not as -1e+05
    expect_error(call_with(-1e5), paste0("`", arg, "` must be a single finite positive integer (got -100000)."), fixed = TRUE)
  }
  expect_error(design_difference(0.1, 0.2, fpc_N = "500"),
               "`fpc_N` must be a single finite positive integer (got class `character`). `fpc_N` is the population size per group. Set `fpc_N = NULL` to skip the finite-population correction.",
               fixed = TRUE)
  expect_error(design_difference(0.1, 0.2, fpc_N = TRUE),
               "`fpc_N` must be a single finite positive integer (got class `logical`).", fixed = TRUE)
  expect_error(design_difference(0.1, 0.2, fpc_N = c(100, 200)),
               "`fpc_N` must be a single finite positive integer (got length = 2).", fixed = TRUE)
  for (bad in c(-5, 0, 100.5, Inf)) {
    expect_error(design_difference(0.1, 0.2, fpc_N = bad),
                 paste0("`fpc_N` must be a single finite positive integer (got ", bad, ")."), fixed = TRUE)
  }
  # Above 2^53 R cannot tell whole numbers apart (and the arithmetic breaks
  # above about 1e150); 2^53 itself works
  expect_error(design_difference(0.1, 0.2, n_sites = 30, icc = 0.05, fpc_N = 1e200),
               "`fpc_N` (1e+200) is larger than 2^53 = 9007199254740992, the largest whole number R can store exactly. No real population is this large: set `fpc_N = NULL`.",
               fixed = TRUE)
  expect_error(design_difference(0.1, 0.2, fpc_N = 2^53 + 2), "is larger than 2^53", fixed = TRUE)
  expect_equal(design_difference(0.1, 0.2, n_sites = 30, icc = 0.05, fpc_N = 2^53)$n_per_group, 283)
  # Just above 2^53 is printed in full, not as 9.007199e+15
  expect_error(design_difference(0.1, 0.2, n_sites = 2^53 + 2, icc = 0.05),
               "`n_sites` (9007199254740994) is larger than 2^53", fixed = TRUE)
  expect_error(design_difference(0.1, 0.2, n_sites = 9.1e15, icc = 0.05),
               "`n_sites` (9100000000000000) is larger than 2^53", fixed = TRUE)
  # The same limit for n_sites and n_per_site; 2^53 itself works
  expect_error(design_difference(0.1, 0.2, n_sites = 1e16, icc = 0.05),
               "`n_sites` (1e+16) is larger than 2^53 = 9007199254740992, the largest whole number R can store exactly. No real study has this many sites.",
               fixed = TRUE)
  expect_error(design_difference(0.1, 0.2, n_per_site = 1e16, icc = 0.05),
               "`n_per_site` (1e+16) is larger than 2^53 = 9007199254740992, the largest whole number R can store exactly. No real site has this many people.",
               fixed = TRUE)
  expect_warning(design_difference(0.1, 0.2, n_sites = 2^53, icc = 0.05),
                 "fewer than the 9007199254740992 sites given", fixed = TRUE)
  # Valid values: n_sites = 30 (DD-C-2), n_per_site = 20 (DD-C-1), fpc_N = 2000 (DD-F-1)
})

test_that("DD-V-13: more sites, or more people per site, than the population is rejected", {
  expect_error(design_difference(0.1, 0.2, n_sites = 100, icc = 0.05, fpc_N = 50),
               "`n_sites` (100) exceeds `fpc_N` (50): you cannot have more sites than people in the population.",
               fixed = TRUE)
  expect_error(design_difference(0.1, 0.2, n_per_site = 100, icc = 0.05, fpc_N = 50),
               "`n_per_site` (100) exceeds `fpc_N` (50): a site cannot have more people than there are in the whole population.",
               fixed = TRUE)
  # Large numbers are printed in full, not as 1e+06 or 1e+05
  expect_error(design_difference(0.1, 0.2, n_sites = 1e6, icc = 0.05, fpc_N = 1e5),
               "`n_sites` (1000000) exceeds `fpc_N` (100000)", fixed = TRUE)
  expect_error(design_difference(0.1, 0.2, n_per_site = 1e6, icc = 0.05, fpc_N = 1e5),
               "`n_per_site` (1000000) exceeds `fpc_N` (100000)", fixed = TRUE)
  # As many sites, or people per site, as the population is allowed. With a
  # population of 2 the design is everyone (2 people), which is not more than
  # the population. With 2 sites that only works because of the population
  # (2 sites alone are worth at most 40 people), so it warns (see DD-F-3).
  expect_warning(res <- design_difference(0.1, 0.2, n_sites = 2, icc = 0.05, fpc_N = 2),
                 "the design tests 2 of the 2 people in each group (1 per site)", fixed = TRUE)
  expect_equal(c(res$n_per_group, res$n_sites, res$n_per_site), c(2, 2, 1))
  expect_no_warning(res <- design_difference(0.1, 0.2, n_per_site = 2, icc = 0.05, fpc_N = 2))
  expect_equal(c(res$n_per_group, res$n_sites, res$n_per_site), c(2, 1, 2))
  # Other valid values: 23 sites with a population of 60 (DD-F-4), 20 per site
  # with a population of 1000 (DD-F-1)
})

test_that("DD-V-14: a site larger than a whole group is rejected", {
  # icc 0.001: deff at 1000 per site is 1.999, so 198.96 x 1.999 = 397.7 -> 398 per group,
  # but one site of 249 would do (n_sites = 1)
  expect_error(design_difference(0.1, 0.2, n_per_site = 1000, icc = 0.001),
               paste0("`n_per_site` (1000) is larger than the number of people needed per group (398), ",
                      "so one site would be more than enough. Use a smaller `n_per_site`, or set ",
                      "`n_sites = 1` to find how many people one site needs."),
               fixed = TRUE)
  expect_equal(design_difference(0.1, 0.2, n_sites = 1, icc = 0.001)$n_per_site, 249)
  # The check is against the clustered need, not the 199 without clustering:
  # 220 per site is above 199 but below the 243 needed, so it is allowed
  res <- design_difference(0.1, 0.2, n_per_site = 220, icc = 0.001)
  expect_equal(res$n_per_group, 243)
  expect_equal(res$n_sites, 2)
  # With icc 0 the group needs 199: a site of 199 is exactly one site, 200 is too big
  expect_equal(design_difference(0.1, 0.2, n_per_site = 199)$n_sites, 1)
  expect_error(design_difference(0.1, 0.2, n_per_site = 200),
               "`n_per_site` (200) is larger than the number of people needed per group (199)",
               fixed = TRUE)
  # Large numbers are printed in full, not as 1e+05
  expect_error(design_difference(0.1, 0.2, n_per_site = 1e5),
               "`n_per_site` (100000) is larger than", fixed = TRUE)
  # Power chosen so that exactly 100000 people are needed (50 versus 51 percent)
  pw <- pnorm((0.01 * sqrt(1e5) - qnorm(0.975) * sqrt(2 * 0.505 * 0.495)) / sqrt(0.25 + 0.51 * 0.49))
  expect_error(design_difference(0.5, 0.51, power = pw, n_per_site = 2e5),
               "`n_per_site` (200000) is larger than the number of people needed per group (100000)", fixed = TRUE)
})


# ---------------------------------------------------------------------------
# Return value (DD-R-)
# ---------------------------------------------------------------------------

test_that("DD-R-1: return list contains all expected fields, in order", {
  expect_named(design_difference(0.1, 0.2),
               c("n_per_group", "n", "n_eff", "n_sites", "n_per_site",
                 "prevalence1", "prevalence2", "apparent_prev1", "apparent_prev2",
                 "difference", "difference_app", "power", "alternative", "conf_level", "sensitivity",
                 "specificity", "icc", "deff", "fpc_N"),
               ignore.order = FALSE)
})

test_that("DD-R-2: inputs are echoed back as supplied", {
  res <- design_difference(0.2, 0.1, power = 0.9, alternative = "greater",
                           sensitivity = 0.9, specificity = 0.95, conf_level = 0.9,
                           n_per_site = 20, icc = 0.05, fpc_N = 5000)
  expect_equal(res[c("prevalence1", "prevalence2", "power", "alternative", "conf_level",
                     "sensitivity", "specificity", "icc", "fpc_N")],
               list(prevalence1 = 0.2, prevalence2 = 0.1, power = 0.9, alternative = "greater",
                    conf_level = 0.9, sensitivity = 0.9, specificity = 0.95, icc = 0.05,
                    fpc_N = 5000))
  expect_equal(res$n, 2 * res$n_per_group)
})
