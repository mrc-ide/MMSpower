# Cost model (MMS-SD budget-officer activity):
#   Total Variable = n * cost_per_sample
#   Total Fixed    = sum over sites i of (fixed_i + transport_i)
#                    (= n_sites * (fixed + transport) when costs are the same)
#   Total          = Variable + Fixed


# ---------------------------------------------------------------------------
# Core arithmetic
# ---------------------------------------------------------------------------

test_that("EC-1: variable cost only (no sites)", {
  res <- estimate_cost(n = 323, cost_per_sample = 50)
  expect_equal(res$total_variable_cost, 323 * 50)
  expect_equal(res$total_fixed_cost, 0)
  expect_equal(res$total_cost, 323 * 50)
  expect_equal(res$n_sites, 0)
  expect_null(res$by_site)
  expect_null(res$budget)
})

test_that("EC-2: same costs at every site -- fixed = n_sites * (fixed + transport)", {
  res <- estimate_cost(n = 600, cost_per_sample = 50, n_sites = 12,
                       fixed_cost_per_site = 5000, transport_cost_per_site = 1000)
  expect_equal(res$total_variable_cost, 600 * 50)
  expect_equal(res$total_fixed_cost, 12 * (5000 + 1000))
  expect_equal(res$total_cost, 600 * 50 + 12 * 6000)
  expect_equal(res$n_sites, 12)
  expect_equal(nrow(res$by_site), 12)
  expect_equal(res$by_site$site, 1:12)
  expect_equal(res$by_site$site_total, rep(6000, 12))
})

test_that("EC-3: a different fixed cost at each site -- hand-checked", {
  res <- estimate_cost(n = 50, cost_per_sample = 10, n_sites = 3,
                       fixed_cost_per_site = c(100, 200, 300))
  expect_equal(res$total_variable_cost, 500)
  expect_equal(res$total_fixed_cost, 600)                 # 100 + 200 + 300
  expect_equal(res$total_cost, 1100)
  expect_equal(res$n_sites, 3)
  bs <- res$by_site
  expect_equal(bs$site, 1:3)
  expect_equal(bs$fixed_cost_per_site, c(100, 200, 300))
  expect_equal(bs$transport_cost_per_site, c(0, 0, 0))
  expect_equal(bs$site_total, c(100, 200, 300))
})

test_that("EC-4: both costs differ by site", {
  res <- estimate_cost(n = 10, cost_per_sample = 1, n_sites = 3,
                       fixed_cost_per_site = c(100, 200, 300),
                       transport_cost_per_site = c(10, 20, 30))
  expect_equal(res$by_site$site_total, c(110, 220, 330))
  expect_equal(res$total_fixed_cost, 660)
  # total is always variable + fixed, and fixed is the sum over sites
  expect_equal(res$total_cost, res$total_variable_cost + res$total_fixed_cost)
  expect_equal(sum(res$by_site$site_total), res$total_fixed_cost)
})

test_that("EC-5: a single cost is applied to every site alongside a vector", {
  res <- estimate_cost(n = 10, cost_per_sample = 1, n_sites = 3,
                       fixed_cost_per_site = c(100, 200, 300),
                       transport_cost_per_site = 50)
  expect_equal(res$by_site$transport_cost_per_site, c(50, 50, 50))
  expect_equal(res$by_site$site_total, c(150, 250, 350))
  expect_equal(res$total_fixed_cost, 750)
})

test_that("EC-6: per-site costs are matched by position; names are ignored", {
  res <- estimate_cost(n = 10, cost_per_sample = 1, n_sites = 2,
                       fixed_cost_per_site = c(B = 200, A = 100))
  expect_equal(res$by_site$fixed_cost_per_site, c(200, 100))
  expect_identical(rownames(res$by_site), c("1", "2"))
  # n_sites given as an integer works the same
  expect_equal(estimate_cost(n = 10, cost_per_sample = 1, n_sites = 2L,
                             fixed_cost_per_site = 100)$total_fixed_cost, 200)
})

test_that("EC-7: non-integer costs, n = 1, cost_per_sample = 0", {
  expect_equal(estimate_cost(n = 200, cost_per_sample = 12.50)$total_cost, 2500)
  expect_equal(estimate_cost(n = 1, cost_per_sample = 99.99)$total_cost, 99.99)
  expect_equal(estimate_cost(n = 500, cost_per_sample = 0)$total_cost, 0)
})

test_that("EC-8: large n * unit cost does not overflow to NA", {
  res <- estimate_cost(n = 100000L, cost_per_sample = 100000L)
  expect_equal(res$total_cost, 1e10)
})


# ---------------------------------------------------------------------------
# Budget
# ---------------------------------------------------------------------------

test_that("EC-B1: budget within -- remaining and flags", {
  res <- estimate_cost(n = 800, cost_per_sample = 50, n_sites = 30,
                       fixed_cost_per_site = 5000, transport_cost_per_site = 1000,
                       budget = 300000)
  # variable 40,000 + fixed 30 * 6000 = 180,000 -> total 220,000
  expect_equal(res$budget, 300000)
  expect_equal(res$budget_remaining, 80000)
  expect_false(res$over_budget)
})

test_that("EC-B2: budget exceeded -- negative remaining, over_budget TRUE", {
  res <- estimate_cost(n = 800, cost_per_sample = 50, n_sites = 20,
                       fixed_cost_per_site = 5000, transport_cost_per_site = 2000,
                       budget = 100000)
  # variable 40,000 + fixed 20*7000 = 140,000 -> total 180,000
  expect_equal(res$total_cost, 180000)
  expect_equal(res$budget_remaining, 100000 - 180000)
  expect_true(res$over_budget)
})

test_that("EC-B3: no budget -> budget fields are NULL", {
  res <- estimate_cost(n = 100, cost_per_sample = 50, n_sites = 5,
                       fixed_cost_per_site = 5000)
  expect_null(res$budget)
  expect_null(res$budget_remaining)
  expect_null(res$over_budget)
})

test_that("EC-B4: money is compared in whole cents, so decimal costs exactly on budget are on budget, up to about 9e13", {
  # 7 * 1.1 is stored a hair above 7.7, so a plain comparison said
  # "over budget by 0.00"
  for (d in list(c(7, 1.1, 7.7), c(3, 0.1, 0.3), c(323, 0.1, 32.3))) {
    res <- estimate_cost(n = d[1], cost_per_sample = d[2], budget = d[3])
    expect_false(res$over_budget)
    expect_identical(res$budget_remaining, 0)
  }
  out <- capture.output(print(estimate_cost(n = 7, cost_per_sample = 1.1, budget = 7.7)))
  expect_true("  Budget 7.70 -> under budget by 0.00" %in% out)
  # One cent short is over budget
  res <- estimate_cost(n = 7, cost_per_sample = 1.1, budget = 7.69)
  expect_true(res$over_budget)
  expect_equal(res$budget_remaining, -0.01)
  expect_true("  Budget 7.69 -> over budget by 0.01" %in% capture.output(print(res)))
  # Money is compared to the cent at any size: one cent over 100 million
  # or 20 billion, or 9 over 10 trillion, is over budget
  res <- estimate_cost(n = 1, cost_per_sample = 1e8 + 0.01, budget = 1e8)
  expect_true(res$over_budget)
  expect_equal(res$budget_remaining, -0.01, tolerance = 1e-6)
  res <- estimate_cost(n = 1, cost_per_sample = 2e10 + 0.01, budget = 2e10)
  expect_true(res$over_budget)
  expect_equal(res$budget_remaining, -0.01, tolerance = 1e-3)
  expect_true("  Budget 20,000,000,000.00 -> over budget by 0.01" %in% capture.output(print(res)))
  res <- estimate_cost(n = 1, cost_per_sample = 1e13 + 9, budget = 1e13)
  expect_true(res$over_budget)
  expect_equal(res$budget_remaining, -9)
  res <- estimate_cost(n = 1, cost_per_sample = 1e13 - 9, budget = 1e13)
  expect_false(res$over_budget)
  expect_equal(res$budget_remaining, 9)
  res <- estimate_cost(n = 1, cost_per_sample = 1e12 + 900, budget = 1e12)
  expect_true(res$over_budget)
  expect_equal(res$budget_remaining, -900)
  res <- estimate_cost(n = 1, cost_per_sample = 1e12 - 900, budget = 1e12)
  expect_false(res$over_budget)
  expect_equal(res$budget_remaining, 900)
  # Over by less than a cent rounds to the same cent, so it is on budget
  res <- estimate_cost(n = 1, cost_per_sample = 100.004, budget = 100)
  expect_false(res$over_budget)
  expect_identical(res$budget_remaining, 0)
  # A budget with a fraction of a cent counts to the nearest cent
  for (b in c(99.996, 100.004)) {
    res <- estimate_cost(n = 1, cost_per_sample = 100, budget = b)
    expect_false(res$over_budget)
    expect_identical(res$budget_remaining, 0)
  }
  # Total and budget are each rounded to the cent first: 100.00 and 100.01
  # (rounding only the difference, 0.002, would give 0)
  res <- estimate_cost(n = 1, cost_per_sample = 100.004, budget = 100.006)
  expect_false(res$over_budget)
  expect_equal(res$budget_remaining, 0.01)
  expect_true("  Budget 100.01 -> under budget by 0.01" %in% capture.output(print(res)))
  # A total of 9 trillion exactly equal to the budget is on budget
  # (round(x, 2) alone gives two different numbers for the same cent here)
  res <- estimate_cost(n = 6105655, cost_per_sample = 1479073.34, budget = 9030711533737.70)
  expect_false(res$over_budget)
  expect_identical(res$budget_remaining, 0)
  # ...while a large total that equals the budget is still on budget
  res <- estimate_cost(n = 1e6, cost_per_sample = 1234.56, budget = 1234560000)
  expect_false(res$over_budget)
  expect_identical(res$budget_remaining, 0)
})


# ---------------------------------------------------------------------------
# Structure / class
# ---------------------------------------------------------------------------

test_that("EC-S1: result class and top-level names", {
  res <- estimate_cost(n = 50, cost_per_sample = 5)
  expect_s3_class(res, "mms_cost")
  expect_setequal(names(res),
                  c("total_cost", "total_variable_cost", "total_fixed_cost",
                    "by_site", "n", "n_sites", "cost_per_sample",
                    "fixed_cost_per_site", "transport_cost_per_site", "budget",
                    "budget_remaining", "over_budget"))
})

test_that("EC-S2: both per-site costs are echoed back as supplied", {
  res <- estimate_cost(n = 10, cost_per_sample = 1, n_sites = 2,
                       fixed_cost_per_site = c(6000, 4500),
                       transport_cost_per_site = 1200)
  expect_equal(res$fixed_cost_per_site, c(6000, 4500))
  expect_equal(res$transport_cost_per_site, 1200)
})

test_that("EC-S3: by_site columns", {
  res <- estimate_cost(n = 10, cost_per_sample = 1, n_sites = 2,
                       fixed_cost_per_site = 100, transport_cost_per_site = 10)
  expect_named(res$by_site,
               c("site", "fixed_cost_per_site", "transport_cost_per_site",
                 "site_total"))
})


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

test_that("EC-V1: n must be a single finite positive integer", {
  type_msg <- function(cls, len) paste0(
    "`n` must be a single positive integer (got class `", cls, "`, length ", len, ").")
  # TRUE would otherwise be read as 1 sample
  expect_error(estimate_cost(n = TRUE, cost_per_sample = 10), type_msg("logical", 1), fixed = TRUE)
  expect_error(estimate_cost(n = "100", cost_per_sample = 10), type_msg("character", 1), fixed = TRUE)
  expect_error(estimate_cost(n = c(1, 2), cost_per_sample = 10), type_msg("numeric", 2), fixed = TRUE)
  for (bad in c(0, -5, 10.5, Inf, -Inf, NaN, NA_real_))
    expect_error(estimate_cost(n = bad, cost_per_sample = 10),
                 paste0("`n` must be a finite positive integer (got ", bad, ")."), fixed = TRUE)
  # A value a hair off an integer says so (100 * 1.1 = 110.00000000000001)
  expect_error(estimate_cost(n = 100 * 1.1, cost_per_sample = 10),
               "`n` must be a finite positive integer (got 110, which is not exactly an integer; use round(n)).",
               fixed = TRUE)
  # Large values too (1e8 * 1.1 is 1.5e-8 above 110000000)
  expect_error(estimate_cost(n = 1e8 * 1.1, cost_per_sample = 10),
               "`n` must be a finite positive integer (got 110000000, which is not exactly an integer; use round(n)).",
               fixed = TRUE)
  # No hint for values that are not a rounding slip, or round to 0
  for (bad in list(c(1e10 + 3.7, "10000000003.7"), c(1e6 + 0.0005, "1000000.0005"),
                   c(1e-12, "0.000000000001"), c(0.4, "0.4")))
    expect_error(estimate_cost(n = as.numeric(bad[1]), cost_per_sample = 10),
                 paste0("`n` must be a finite positive integer (got ", bad[2], ")."),
                 fixed = TRUE)
  # valid value: n = 1, see EC-7
})

test_that("EC-V2: cost_per_sample must be a single finite non-negative number", {
  type_msg <- function(cls, len) paste0(
    "`cost_per_sample` must be a single non-negative number (got class `", cls,
    "`, length ", len, ").")
  expect_error(estimate_cost(n = 100, cost_per_sample = TRUE), type_msg("logical", 1), fixed = TRUE)
  expect_error(estimate_cost(n = 100, cost_per_sample = "5"), type_msg("character", 1), fixed = TRUE)
  expect_error(estimate_cost(n = 100, cost_per_sample = c(5, 10)), type_msg("numeric", 2), fixed = TRUE)
  for (bad in c(-5, -0.01, Inf, -Inf, NaN, NA_real_))
    expect_error(estimate_cost(n = 100, cost_per_sample = bad),
                 paste0("`cost_per_sample` must be a finite non-negative number (got ", bad, ")."),
                 fixed = TRUE)
  # valid value: cost_per_sample = 0, see EC-7
})

test_that("EC-V3: n_sites must be a single number; the region hint only for several values", {
  msg <- function(cls, len) paste0(
    "`n_sites` must be a single positive integer or NULL (got class `", cls,
    "`, length ", len, ").",
    if (len > 1) " To cost several regions, call estimate_cost() once per region and add the totals.")
  expect_error(estimate_cost(n = 100, cost_per_sample = 5, n_sites = c(3, 4),
                             fixed_cost_per_site = 10),
               msg("numeric", 2), fixed = TRUE)
  expect_error(estimate_cost(n = 100, cost_per_sample = 5,
                             n_sites = c(North = 10, South = 8),
                             fixed_cost_per_site = 10),
               msg("numeric", 2), fixed = TRUE)
  # A single wrong value gets the message without the region hint (the
  # whole message is compared, so extra text would be caught). TRUE would
  # otherwise be read as 1 site.
  err <- function(expr) tryCatch(expr, error = conditionMessage)
  expect_identical(err(estimate_cost(n = 100, cost_per_sample = 5, n_sites = TRUE,
                                     fixed_cost_per_site = 10)), msg("logical", 1))
  expect_identical(err(estimate_cost(n = 100, cost_per_sample = 5, n_sites = "3",
                                     fixed_cost_per_site = 10)), msg("character", 1))
  expect_identical(err(estimate_cost(n = 100, cost_per_sample = 5, n_sites = numeric(0),
                                     fixed_cost_per_site = 10)), msg("numeric", 0))
  expect_identical(err(estimate_cost(n = 100, cost_per_sample = 5, n_sites = c(3, 4),
                                     fixed_cost_per_site = 10)), msg("numeric", 2))
})

test_that("EC-V4: n_sites must be a finite positive integer", {
  for (bad in c(0, -3, 2.5, Inf, -Inf, NaN, NA_real_))
    expect_error(estimate_cost(n = 100, cost_per_sample = 5, n_sites = bad,
                               fixed_cost_per_site = 10),
                 paste0("`n_sites` must be a finite positive integer (got ",
                        bad, ")."), fixed = TRUE)
  # A value a hair off an integer says so (0.1 * 3 * 10 = 3.0000000000000004)
  expect_error(estimate_cost(n = 100, cost_per_sample = 5, n_sites = 0.1 * 3 * 10,
                             fixed_cost_per_site = 10),
               "`n_sites` must be a finite positive integer (got 3, which is not exactly an integer; use round(n_sites)).",
               fixed = TRUE)
  # valid value: see EC-2
})

test_that("EC-V5: a per-site cost vector must have 1 or n_sites values", {
  expect_error(
    estimate_cost(n = 10, cost_per_sample = 1, n_sites = 3,
                  fixed_cost_per_site = c(100, 200)),
    paste0("`fixed_cost_per_site` has 2 values but `n_sites` is 3. Give a ",
           "single number for every site, or one value per site ",
           "(3 values)."), fixed = TRUE)
  expect_error(
    estimate_cost(n = 10, cost_per_sample = 1, n_sites = 3,
                  fixed_cost_per_site = 100,
                  transport_cost_per_site = c(10, 20, 30, 40)),
    paste0("`transport_cost_per_site` has 4 values but `n_sites` is 3. Give a ",
           "single number for every site, or one value per site ",
           "(3 values)."), fixed = TRUE)
  # valid values: see EC-3 (vector) and EC-5 (vector + single number)
})

test_that("EC-V6: per-site cost values must be finite and non-negative", {
  start <- function(nm) paste0(
    "`", nm, "` must be a finite non-negative number, or a vector of ",
    "non-negative costs with one value per site (")
  msg <- function(nm, cls, len) paste0(start(nm), "got class `", cls, "`, length ", len, ").")
  # Wrong type or empty: the class and length are given
  expect_error(estimate_cost(n = 100, cost_per_sample = 5, n_sites = 3,
                             fixed_cost_per_site = "100"),
               msg("fixed_cost_per_site", "character", 1), fixed = TRUE)
  expect_error(estimate_cost(n = 100, cost_per_sample = 5, n_sites = 3,
                             transport_cost_per_site = "10"),
               msg("transport_cost_per_site", "character", 1), fixed = TRUE)
  expect_error(estimate_cost(n = 100, cost_per_sample = 5, n_sites = 3,
                             fixed_cost_per_site = numeric(0)),
               msg("fixed_cost_per_site", "numeric", 0), fixed = TRUE)
  # A single bad value is shown
  expect_error(estimate_cost(n = 100, cost_per_sample = 5, n_sites = 3,
                             fixed_cost_per_site = -1),
               paste0(start("fixed_cost_per_site"), "got -1)."), fixed = TRUE)
  expect_error(estimate_cost(n = 100, cost_per_sample = 5, n_sites = 5,
                             fixed_cost_per_site = 100,
                             transport_cost_per_site = -50),
               paste0(start("transport_cost_per_site"), "got -50)."), fixed = TRUE)
  # In a vector, the first bad position is named (here 2, not 3)
  expect_error(estimate_cost(n = 100, cost_per_sample = 5, n_sites = 3,
                             fixed_cost_per_site = c(100, NA, -300)),
               paste0(start("fixed_cost_per_site"), "found fixed_cost_per_site[2] = NA)."),
               fixed = TRUE)
  expect_error(estimate_cost(n = 100, cost_per_sample = 5, n_sites = 3,
                             transport_cost_per_site = c(10, 20, Inf)),
               paste0(start("transport_cost_per_site"), "found transport_cost_per_site[3] = Inf)."),
               fixed = TRUE)
  # 0 is allowed
  expect_no_error(estimate_cost(n = 100, cost_per_sample = 5, n_sites = 3,
                                fixed_cost_per_site = c(0, 100, 0)))
})

test_that("EC-V7: budget must be a single finite positive number or NULL", {
  type_msg <- function(cls, len) paste0(
    "`budget` must be a single positive number or NULL (got class `", cls,
    "`, length ", len, ").")
  # Text is reported by its class, not printed as if it were the number 100
  expect_error(estimate_cost(n = 100, cost_per_sample = 5, budget = "100"),
               type_msg("character", 1), fixed = TRUE)
  expect_error(estimate_cost(n = 100, cost_per_sample = 5, budget = TRUE),
               type_msg("logical", 1), fixed = TRUE)
  expect_error(estimate_cost(n = 100, cost_per_sample = 5, budget = c(1, 2)),
               type_msg("numeric", 2), fixed = TRUE)
  for (bad in c(0, -1, Inf, -Inf, NaN, NA_real_))
    expect_error(estimate_cost(n = 100, cost_per_sample = 5, budget = bad),
                 paste0("`budget` must be a finite positive number (got ", bad, ")."),
                 fixed = TRUE)
  # A small positive budget is allowed; no budget (NULL) is EC-B3
  expect_no_error(estimate_cost(n = 100, cost_per_sample = 5, budget = 0.01))
})

test_that("EC-V8: per-site costs with no n_sites -> warn, still compute", {
  expect_warning(
    res <- estimate_cost(n = 100, cost_per_sample = 50, fixed_cost_per_site = 5000),
    "`fixed_cost_per_site` and `transport_cost_per_site` are ignored because `n_sites` is NULL. Supply `n_sites` to include per-site costs.",
    fixed = TRUE
  )
  expect_equal(res$total_cost, 100 * 50)
  # a per-site vector with no n_sites warns the same way
  expect_warning(
    estimate_cost(n = 100, cost_per_sample = 50,
                  transport_cost_per_site = c(10, 20)),
    "`fixed_cost_per_site` and `transport_cost_per_site` are ignored because `n_sites` is NULL. Supply `n_sites` to include per-site costs.",
    fixed = TRUE
  )
  # default (0) per-site costs with no n_sites is silent
  expect_silent(estimate_cost(n = 100, cost_per_sample = 50))
})

test_that("EC-V9: n_sites with no per-site costs -> warn, still compute", {
  expect_warning(
    res <- estimate_cost(n = 100, cost_per_sample = 50, n_sites = 5),
    paste0("`n_sites` supplied but both per-site costs are 0. Supply ",
           "`fixed_cost_per_site` or `transport_cost_per_site` to cost the sites."),
    fixed = TRUE
  )
  expect_equal(res$total_cost, 100 * 50)   # fixed cost is 0
  expect_equal(res$n_sites, 5)
  # supplying either cost silences it
  expect_silent(estimate_cost(n = 100, cost_per_sample = 50, n_sites = 5,
                              transport_cost_per_site = 800))
})

test_that("EC-V10: numbers in messages are printed in full", {
  expect_error(estimate_cost(n = -1e5, cost_per_sample = 10),
               "`n` must be a finite positive integer (got -100000).", fixed = TRUE)
  expect_error(estimate_cost(n = 10.0000001, cost_per_sample = 10),
               "`n` must be a finite positive integer (got 10.0000001).", fixed = TRUE)
  expect_error(estimate_cost(n = 10, cost_per_sample = -1e5),
               "`cost_per_sample` must be a finite non-negative number (got -100000).", fixed = TRUE)
  expect_error(estimate_cost(n = 10, cost_per_sample = 1, budget = -1e6),
               "`budget` must be a finite positive number (got -1000000).", fixed = TRUE)
  expect_error(estimate_cost(n = 10, cost_per_sample = 1, n_sites = -1e5),
               "`n_sites` must be a finite positive integer (got -100000).", fixed = TRUE)
  expect_error(estimate_cost(n = 10, cost_per_sample = 1, n_sites = 3,
                             fixed_cost_per_site = c(1, -1e6, 1)),
               "(found fixed_cost_per_site[2] = -1000000).", fixed = TRUE)
  expect_error(estimate_cost(n = 10, cost_per_sample = 1, n_sites = 1e6,
                             fixed_cost_per_site = c(1, 2)),
               paste0("`fixed_cost_per_site` has 2 values but `n_sites` is 1000000. Give a ",
                      "single number for every site, or one value per site (1000000 values)."),
               fixed = TRUE)
  # One site: "1 value"
  expect_error(estimate_cost(n = 10, cost_per_sample = 1, n_sites = 1,
                             fixed_cost_per_site = c(1, 2)),
               "one value per site (1 value).", fixed = TRUE)
})

test_that("EC-V11: a named number or a 1 x 1 matrix works like a plain number", {
  expect_identical(
    estimate_cost(n = c(a = 50), cost_per_sample = matrix(10), n_sites = c(s = 3),
                  fixed_cost_per_site = 5, budget = c(b = 1000)),
    estimate_cost(n = 50, cost_per_sample = 10, n_sites = 3,
                  fixed_cost_per_site = 5, budget = 1000))
})


# ---------------------------------------------------------------------------
# Print method
# ---------------------------------------------------------------------------

test_that("EC-P1: print returns the object invisibly and shows the split", {
  res <- estimate_cost(n = 100, cost_per_sample = 50, n_sites = 5,
                       fixed_cost_per_site = 5000, transport_cost_per_site = 1000)
  out <- capture.output(ret <- print(res))
  expect_identical(ret, res)
  expect_true(any(grepl("Variable:", out)))
  expect_true(any(grepl("Fixed:",    out)))
  expect_true(any(grepl("Total",     out)))
})

test_that("EC-P2: same costs everywhere -> one 'sites x' line, no site lines", {
  res <- estimate_cost(n = 600, cost_per_sample = 50, n_sites = 12,
                       fixed_cost_per_site = 5000, transport_cost_per_site = 1000)
  out <- capture.output(print(res))
  expect_true(any(grepl("(fixed + transport cost per site)", out, fixed = TRUE)))
  expect_true(any(grepl("12 sites x (5,000.00 + 1,000.00)", out, fixed = TRUE)))
  expect_false(any(grepl("site 1", out, fixed = TRUE)))
})

test_that("EC-P3: costs that differ -> one line per site, plus budget verdict", {
  res <- estimate_cost(n = 50, cost_per_sample = 10, n_sites = 3,
                       fixed_cost_per_site = c(100, 200, 300), budget = 1500)
  out <- capture.output(print(res))
  expect_true(any(grepl("(fixed + transport cost per site)", out, fixed = TRUE)))
  expect_true(any(grepl("site 1    (100.00 + 0.00)   =   100.00", out, fixed = TRUE)))
  expect_true(any(grepl("site 3    (300.00 + 0.00)   =   300.00", out, fixed = TRUE)))
  expect_false(any(grepl("not shown", out)))
  expect_true(any(grepl("under budget by 400.00", out, fixed = TRUE)))
  # leftover 400: 40 more samples at 10, or 2 more sites at the
  # average cost 600 / 3 = 200 -> "one of:" header
  expect_true("      40 more samples, or" %in% out)
  expect_true("      2 more sites (at the average per-site cost)" %in% out)
  expect_true("    leftover budget would cover roughly one of:" %in% out)
})

test_that("EC-P4: more than 10 differing sites -> first 10 + a note", {
  res <- estimate_cost(n = 100, cost_per_sample = 1, n_sites = 25,
                       fixed_cost_per_site = 1:25 * 100)
  out <- capture.output(print(res))
  expect_true(any(grepl("site 10 ", out, fixed = TRUE)))
  expect_false(any(grepl("site 11 ", out, fixed = TRUE)))
  expect_true(any(grepl("... 15 more sites not shown (all 25 are in $by_site)",
                        out, fixed = TRUE)))
  # The dashed line above the total is as wide as the note
  note <- out[grepl("not shown", out, fixed = TRUE)]
  dash <- out[grepl("^  -+$", out)]
  expect_identical(nchar(dash), nchar(note))
  # exactly 10 differing sites -> all shown, no note
  res10 <- estimate_cost(n = 100, cost_per_sample = 1, n_sites = 10,
                         fixed_cost_per_site = 1:10 * 100)
  out10 <- capture.output(print(res10))
  expect_true(any(grepl("site 10 ", out10, fixed = TRUE)))
  expect_false(any(grepl("not shown", out10)))
})

test_that("EC-P5: a single leftover option uses 'roughly:' not 'one of:'", {
  # variable-only, under budget -> only "more samples" is offered
  res <- estimate_cost(n = 100, cost_per_sample = 50, budget = 10000)
  out <- capture.output(print(res))
  expect_true(any(grepl("would cover roughly:", out, fixed = TRUE)))
  expect_false(any(grepl("one of", out)))
  # 5,000 left at 50 a sample: 100 more
  expect_true("      100 more samples" %in% out)
})

test_that("EC-P6: over-budget print says over budget", {
  res <- estimate_cost(n = 800, cost_per_sample = 50, n_sites = 20,
                       fixed_cost_per_site = 5000, transport_cost_per_site = 2000,
                       budget = 100000)
  out <- capture.output(print(res))
  expect_true("  Budget 100,000.00 -> over budget by 80,000.00" %in% out)
  expect_false(any(grepl("would cover", out, fixed = TRUE)))
})

test_that("EC-P7: large totals print grouped, not scientific", {
  res <- estimate_cost(n = 1e6, cost_per_sample = 10)
  out <- capture.output(print(res))
  expect_true(any(grepl("10,000,000", out, fixed = TRUE)))
  expect_false(any(grepl("e\\+0", out)))
})

test_that("EC-P8: exactly on budget -- no 'would cover' block", {
  res <- estimate_cost(n = 1000, cost_per_sample = 100, budget = 100000)
  expect_equal(res$budget_remaining, 0)
  expect_false(res$over_budget)
  out <- capture.output(print(res))
  expect_true(any(grepl("under budget by 0.00", out, fixed = TRUE)))
  expect_false(any(grepl("would cover", out)))
})

test_that("EC-P9: the = signs and amounts line up, however large the numbers", {
  # Bob's example, printed in full
  res <- estimate_cost(n = 50, cost_per_sample = 10, n_sites = 3,
                       fixed_cost_per_site = c(100, 200, 300))
  expect_identical(capture.output(print(res)), c(
    "Study cost estimate",
    "  Variable:  50 samples x 10.00 =   500.00",
    "  Fixed:     3 sites            =   600.00",
    "    (fixed + transport cost per site)",
    "    site 1    (100.00 + 0.00)   =   100.00",
    "    site 2    (200.00 + 0.00)   =   200.00",
    "    site 3    (300.00 + 0.00)   =   300.00",
    "  ----------------------------------------",
    "  Total                         = 1,100.00"))
  # Huge numbers: every line with an amount has its = in the same column
  # and ends at the same width, and the widest line still fits in 80
  # characters
  res <- estimate_cost(n = 1e8, cost_per_sample = 1234.56, n_sites = 3,
                       fixed_cost_per_site = c(2.5e6, 12e6, 750000),
                       transport_cost_per_site = c(40000, 1.1e6, 5000))
  out <- capture.output(print(res))
  with_amount <- out[grepl(" = ", out, fixed = TRUE)]
  expect_length(with_amount, 6)
  expect_length(unique(regexpr(" = ", with_amount, fixed = TRUE)), 1)
  expect_length(unique(nchar(with_amount)), 1)
  expect_lte(max(nchar(out)), 80)
  expect_true(any(grepl("Total                                      = 123,472,395,000.00",
                        out, fixed = TRUE)))
})

test_that("EC-P10: the leftover hint rounds down, and rounding error does not lose one", {
  hint <- function(res) {
    out <- capture.output(print(res))
    i <- grep("would cover", out, fixed = TRUE)
    if (length(i) == 0) character(0) else out[-seq_len(i)]
  }
  # 149.97 left at 49.99 a sample is exactly 3 samples (the division comes
  # out a hair below 3)
  expect_identical(hint(estimate_cost(n = 7, cost_per_sample = 49.99, budget = 499.9)),
                   "      3 more samples")
  # 0.29 left at 0.01 a sample is 29 samples (0.29 / 0.01 comes out as
  # 28.999999999999996)
  expect_identical(hint(estimate_cost(n = 1, cost_per_sample = 0.01, budget = 0.30)),
                   "      29 more samples")
  # Large counts are exact too: 521,373,200.84 left at 93.26 is exactly
  # 5,590,534 samples (the division comes out a hair below), and
  # 30,000,000,002.99 left at 3 is 10,000,000,000 samples, not one more
  expect_identical(hint(estimate_cost(n = 1, cost_per_sample = 93.26, budget = 521373294.10)),
                   "      5,590,534 more samples")
  expect_identical(hint(estimate_cost(n = 1, cost_per_sample = 3, budget = 3e10 + 5.99)),
                   "      10,000,000,000 more samples")
  expect_identical(hint(estimate_cost(n = 1, cost_per_sample = 0, n_sites = 2,
                                      fixed_cost_per_site = 1.5,
                                      budget = 3 + 1.5e10 + 1.49)),
                   "      10,000,000,000 more sites (at the average per-site cost)")
  # 350 left at 200 a sample is 1.75 samples: 1, not 2
  expect_identical(hint(estimate_cost(n = 1, cost_per_sample = 200, budget = 550)),
                   "      1 more sample")
  # 30.30 left at 10.10 a site is exactly 3 sites
  expect_identical(hint(estimate_cost(n = 1, cost_per_sample = 0, n_sites = 3,
                                      fixed_cost_per_site = 10.1, budget = 60.6)),
                   "      3 more sites (at the average per-site cost)")
  # 15 left at 10 a site: 1 site
  expect_identical(hint(estimate_cost(n = 1, cost_per_sample = 0, n_sites = 1,
                                      fixed_cost_per_site = 10, budget = 25)),
                   "      1 more site (at the average per-site cost)")
  # Free samples are not offered; nor are sites that cost nothing
  expect_identical(hint(suppressWarnings(
    estimate_cost(n = 5, cost_per_sample = 10, n_sites = 2, budget = 100))),
    "      5 more samples")
  # Many decimal site costs: 0.01 left (stored as 0.0099999999802) buys 1
  # sample at 0.01, and no site, so only the sample is offered
  expect_identical(hint(estimate_cost(n = 1088, cost_per_sample = 0.01, n_sites = 11,
                                      fixed_cost_per_site = 10331.70,
                                      transport_cost_per_site = 1716.11,
                                      budget = 132536.80)),
                   "      1 more sample")
  # The average site cost includes transport: 300 left at 100 + 50 a site
  # is 2 sites, not 3
  expect_identical(hint(estimate_cost(n = 1, cost_per_sample = 0, n_sites = 2,
                                      fixed_cost_per_site = 100,
                                      transport_cost_per_site = 50, budget = 600)),
                   "      2 more sites (at the average per-site cost)")
  # 5 left at 10 a sample cannot buy one: no hint at all
  out <- capture.output(print(estimate_cost(n = 1, cost_per_sample = 10, budget = 15)))
  expect_true("  Budget 15.00 -> under budget by 5.00" %in% out)
  expect_false(any(grepl("would cover", out, fixed = TRUE)))
  # Less than a cent left: no hint at all
  out <- capture.output(print(estimate_cost(n = 1, cost_per_sample = 10, budget = 10.004)))
  expect_true("  Budget 10.00 -> under budget by 0.00" %in% out)
  expect_false(any(grepl("would cover", out, fixed = TRUE)))
})

test_that("EC-P11: one sample or one site is written in the singular", {
  out <- capture.output(print(estimate_cost(n = 1, cost_per_sample = 5, n_sites = 1,
                                            fixed_cost_per_site = 10)))
  expect_true("  Variable:  1 sample x 5.00 =  5.00" %in% out)
  expect_true("  Fixed:     1 site          = 10.00" %in% out)
  expect_true("    1 site x (10.00 + 0.00)" %in% out)
  out <- capture.output(print(estimate_cost(n = 1, cost_per_sample = 1, n_sites = 11,
                                            fixed_cost_per_site = 1:11)))
  expect_true("    ... 1 more site not shown (all 11 are in $by_site)" %in% out)
})

test_that("EC-P12: sites whose transport costs differ are listed one per line", {
  out <- capture.output(print(estimate_cost(n = 1, cost_per_sample = 1, n_sites = 2,
                                            fixed_cost_per_site = 100,
                                            transport_cost_per_site = c(10, 20))))
  expect_true("    site 1    (100.00 + 10.00) = 110.00" %in% out)
  expect_true("    site 2    (100.00 + 20.00) = 120.00" %in% out)
  expect_false(any(grepl("2 sites x", out, fixed = TRUE)))
})

test_that("EC-P13: printed amounts are not affected by options(digits)", {
  op <- options(digits = 22)
  on.exit(options(op))
  out <- capture.output(print(estimate_cost(n = 7, cost_per_sample = 1.1, budget = 7.7)))
  expect_true("  Variable:  7 samples x 1.10 = 7.70" %in% out)
  expect_true("  Budget 7.70 -> under budget by 0.00" %in% out)
})

test_that("EC-P14: a cost of -0 prints as 0.00, not -0.00", {
  out <- capture.output(print(suppressWarnings(
    estimate_cost(n = 2, cost_per_sample = 10, n_sites = 2, fixed_cost_per_site = -0))))
  expect_true("    2 sites x (0.00 + 0.00)" %in% out)
  expect_false(any(grepl("-0.00", out, fixed = TRUE)))
})
