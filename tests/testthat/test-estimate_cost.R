# Cost model (MMS-SD budget-officer activity):
#   Total Variable = n * cost_per_sample
#   Total Fixed    = sum_r  n_HF_r * (fixed_cost_per_site + transport_r)
#   Total          = Variable + Fixed


# ---------------------------------------------------------------------------
# Round 1 -- core arithmetic
# ---------------------------------------------------------------------------

test_that("EC-1: variable cost only (no facilities)", {
  res <- estimate_cost(n = 323, cost_per_sample = 50)
  expect_equal(res$total_variable_cost, 323 * 50)
  expect_equal(res$total_fixed_cost, 0)
  expect_equal(res$total_cost, 323 * 50)
  expect_equal(res$n_sites, 0)
  expect_null(res$by_region)
  expect_null(res$budget)
})

test_that("EC-2: single region -- fixed = n_HF * (fixed + transport)", {
  res <- estimate_cost(n = 600, cost_per_sample = 50, n_sites = 12,
                       fixed_cost_per_site = 5000, transport_cost_per_site = 1000)
  expect_equal(res$total_variable_cost, 600 * 50)
  expect_equal(res$total_fixed_cost, 12 * (5000 + 1000))
  expect_equal(res$total_cost, 600 * 50 + 12 * 6000)
  expect_equal(res$n_sites, 12)
  expect_equal(nrow(res$by_region), 1)
  expect_equal(res$by_region$region, "(all)")
})

test_that("EC-3: multi-region with per-region transport -- hand-checked", {
  # North 10, South 8, West 12; fixed 5000; transport 1000 / 1500 / 800
  res <- estimate_cost(
    n = 800, cost_per_sample = 50,
    n_sites = c(North = 10, South = 8, West = 12),
    fixed_cost_per_site = 5000,
    transport_cost_per_site = c(North = 1000, South = 1500, West = 800)
  )
  expect_equal(res$total_variable_cost, 800 * 50)                 # 40,000
  expect_equal(res$total_fixed_cost,
               10 * 6000 + 8 * 6500 + 12 * 5800)                  # 60k + 52k + 69.6k = 181,600
  expect_equal(res$total_cost, 40000 + 181600)                    # 221,600
  # total is always variable + fixed, and fixed is the sum of the regions
  # (was EC-8, a random sweep -- folded here as fixed numbers).
  expect_equal(res$total_cost, res$total_variable_cost + res$total_fixed_cost)
  expect_equal(res$n_sites, 30)

  br <- res$by_region
  expect_equal(br$region, c("North", "South", "West"))
  expect_equal(br$region_total, c(60000, 52000, 69600))
  expect_equal(br$transport_subtotal, c(10 * 1000, 8 * 1500, 12 * 800))
  expect_equal(sum(br$region_total), res$total_fixed_cost)
})

test_that("EC-4: scalar transport is applied to every region", {
  res <- estimate_cost(n = 100, cost_per_sample = 50,
                       n_sites = c(A = 3, B = 5), fixed_cost_per_site = 5000,
                       transport_cost_per_site = 900)
  expect_equal(res$by_region$transport_cost_per_site, c(900, 900))
  expect_equal(res$total_fixed_cost, (3 + 5) * (5000 + 900))
})

test_that("EC-5: transport vector may be given in any order; extras warn", {
  # matched by name, order-independent
  a <- estimate_cost(n = 10, cost_per_sample = 1,
                     n_sites = c(North = 2, South = 4), fixed_cost_per_site = 100,
                     transport_cost_per_site = c(South = 50, North = 20))
  expect_equal(a$by_region$transport_cost_per_site, c(20, 50))

  # a name not in n_sites (here "East") is dropped, with a warning
  expect_warning(
    b <- estimate_cost(n = 10, cost_per_sample = 1,
                       n_sites = c(North = 2, South = 4), fixed_cost_per_site = 100,
                       transport_cost_per_site = c(South = 50, North = 20, East = 999)),
    "not in `n_sites`"
  )
  expect_equal(b$by_region$transport_cost_per_site, c(20, 50))
})

test_that("EC-6: non-integer costs, n = 1, cost_per_sample = 0", {
  expect_equal(estimate_cost(n = 200, cost_per_sample = 12.50)$total_cost, 2500)
  expect_equal(estimate_cost(n = 1, cost_per_sample = 99.99)$total_cost, 99.99)
  expect_equal(estimate_cost(n = 500, cost_per_sample = 0)$total_cost, 0)
})

test_that("EC-7: large n * unit cost does not overflow to NA", {
  res <- estimate_cost(n = 100000L, cost_per_sample = 100000L)
  expect_equal(res$total_cost, 1e10)
  expect_false(is.na(res$total_cost))
})

# ---------------------------------------------------------------------------
# Round 2 -- budget
# ---------------------------------------------------------------------------

test_that("EC-B1: budget within -- remaining and flags", {
  res <- estimate_cost(
    n = 800, cost_per_sample = 50,
    n_sites = c(North = 10, South = 8, West = 12),
    fixed_cost_per_site = 5000,
    transport_cost_per_site = c(North = 1000, South = 1500, West = 800),
    budget = 300000
  )
  expect_equal(res$budget, 300000)
  expect_equal(res$budget_remaining, 300000 - 221600)   # 78,400
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


# ---------------------------------------------------------------------------
# Round 3 -- structure / class
# ---------------------------------------------------------------------------

test_that("EC-S1: result class and top-level names", {
  res <- estimate_cost(n = 50, cost_per_sample = 5)
  expect_s3_class(res, "mms_cost")
  expect_true(all(c("total_cost", "total_variable_cost", "total_fixed_cost",
                    "by_region", "n", "n_sites", "cost_per_sample",
                    "fixed_cost_per_site", "transport_cost_per_site", "budget",
                    "budget_remaining", "over_budget") %in% names(res)))
})

test_that("EC-S1b: both per-facility rates are echoed back as supplied", {
  res <- estimate_cost(n = 10, cost_per_sample = 1,
                       n_sites = c(North = 3, South = 4),
                       fixed_cost_per_site = c(North = 6000, South = 4500),
                       transport_cost_per_site = 1200)
  expect_equal(res$fixed_cost_per_site, c(North = 6000, South = 4500))
  expect_equal(res$transport_cost_per_site, 1200)
})

test_that("EC-S2: by_region columns", {
  res <- estimate_cost(n = 10, cost_per_sample = 1, n_sites = c(X = 2, Y = 3),
                       fixed_cost_per_site = 100, transport_cost_per_site = 10)
  expect_named(res$by_region,
               c("region", "n_sites", "fixed_cost_per_site",
                 "transport_cost_per_site", "fixed_subtotal",
                 "transport_subtotal", "region_total"))
})


# ---------------------------------------------------------------------------
# Round 4 -- validation
# ---------------------------------------------------------------------------

test_that("EC-V1: n guards", {
  expect_error(estimate_cost(n = 0, cost_per_sample = 10),    "positive integer")
  expect_error(estimate_cost(n = -5, cost_per_sample = 10),   "positive integer")
  expect_error(estimate_cost(n = 10.5, cost_per_sample = 10), "positive integer")
  expect_error(estimate_cost(n = Inf, cost_per_sample = 10),  "finite")
  expect_error(estimate_cost(n = "100", cost_per_sample = 10), "character")
  expect_error(estimate_cost(n = c(1, 2), cost_per_sample = 10), "single")
})

test_that("EC-V2: cost_per_sample guards", {
  expect_error(estimate_cost(n = 100, cost_per_sample = -5),        "non-negative")
  expect_error(estimate_cost(n = 100, cost_per_sample = NA_real_),  "finite")
  expect_error(estimate_cost(n = 100, cost_per_sample = c(5, 10)),  "single")
})

test_that("EC-V3: fixed_cost_per_site guards", {
  expect_error(estimate_cost(n = 100, cost_per_sample = 5, n_sites = 3,
                             fixed_cost_per_site = -1), "non-negative")
  # The unnamed-multi-value "must be named" path is shared with
  # transport_cost_per_site through .resolve_region_cost(); EC-V5 covers it.
  # per-region fixed cost missing a region
  expect_error(estimate_cost(n = 100, cost_per_sample = 5,
                             n_sites = c(A = 2, B = 3),
                             fixed_cost_per_site = c(A = 10)),
               "missing a value for region")
})

test_that("EC-V3b: fixed_cost_per_site varies by region", {
  res <- estimate_cost(n = 10, cost_per_sample = 15,
                       n_sites = c(North = 3, South = 4),
                       fixed_cost_per_site = c(North = 3, South = 100))
  expect_equal(res$total_fixed_cost, 3 * 3 + 4 * 100)
  expect_equal(res$by_region$fixed_cost_per_site, c(3, 100))
  expect_equal(res$total_cost, 9 + 400 + 10 * 15)
})

test_that("EC-V4: n_sites guards", {
  expect_error(estimate_cost(n = 100, cost_per_sample = 5, n_sites = 0),
               "positive integer")
  expect_error(estimate_cost(n = 100, cost_per_sample = 5, n_sites = 2.5),
               "positive integer")
  expect_error(estimate_cost(n = 100, cost_per_sample = 5, n_sites = c(3, 4)),
               "must name every element")
  # partially named multi-region vector is rejected too
  expect_error(estimate_cost(n = 100, cost_per_sample = 5,
                             n_sites = c(North = 3, 4)),
               "must name every element")
  # duplicate region names are rejected
  expect_error(estimate_cost(n = 100, cost_per_sample = 5,
                             n_sites = c(A = 3, A = 4)),
               "repeated region name")
})

test_that("EC-V5: transport_cost_per_site guards", {
  expect_error(
    estimate_cost(n = 10, cost_per_sample = 1, n_sites = c(A = 2, B = 3),
                  fixed_cost_per_site = 100, transport_cost_per_site = c(A = 10)),
    "missing a value for region"
  )
  expect_error(
    estimate_cost(n = 10, cost_per_sample = 1, n_sites = c(A = 2, B = 3),
                  fixed_cost_per_site = 100, transport_cost_per_site = c(10, 20)),
    "must be named"
  )
  # a region named more than once in the cost vector is rejected
  expect_error(
    estimate_cost(n = 10, cost_per_sample = 1, n_sites = c(A = 2, B = 3),
                  fixed_cost_per_site = 100,
                  transport_cost_per_site = c(A = 10, A = 20, B = 5)),
    "more than one value for region"
  )
  expect_error(
    estimate_cost(n = 10, cost_per_sample = 1, n_sites = 5,
                  fixed_cost_per_site = 100, transport_cost_per_site = -50),
    "non-negative"
  )
})

test_that("EC-V6: budget guard", {
  expect_error(estimate_cost(n = 100, cost_per_sample = 5, budget = -1),
               "positive number")
  expect_error(estimate_cost(n = 100, cost_per_sample = 5, budget = c(1, 2)),
               "length = 2")
})

test_that("EC-V7: per-facility costs with no n_sites -> warn, still compute", {
  expect_warning(
    res <- estimate_cost(n = 100, cost_per_sample = 50, fixed_cost_per_site = 5000),
    "ignored"
  )
  expect_equal(res$total_cost, 100 * 50)
  # default (0) per-facility costs with no n_sites is silent
  expect_silent(estimate_cost(n = 100, cost_per_sample = 50))
})

test_that("EC-V8: n_sites with no per-facility rates -> warn, still compute", {
  expect_warning(
    res <- estimate_cost(n = 100, cost_per_sample = 50, n_sites = 5),
    "both per-facility rates are 0"
  )
  expect_equal(res$total_cost, 100 * 50)   # fixed cost is 0
  expect_equal(res$n_sites, 5)
  # supplying either rate silences it
  expect_silent(estimate_cost(n = 100, cost_per_sample = 50, n_sites = 5,
                              transport_cost_per_site = 800))
})


# ---------------------------------------------------------------------------
# Round 5 -- print method
# ---------------------------------------------------------------------------

test_that("EC-P1: print returns the object invisibly and shows the split", {
  res <- estimate_cost(n = 100, cost_per_sample = 50, n_sites = c(N = 5),
                       fixed_cost_per_site = 5000, transport_cost_per_site = 1000)
  expect_identical(print(res), res)
  out <- capture.output(print(res))
  expect_true(any(grepl("Variable:", out)))
  expect_true(any(grepl("Fixed:",    out)))
  expect_true(any(grepl("Total",     out)))
})

test_that("EC-P2: print shows per-region lines and the budget verdict", {
  res <- estimate_cost(
    n = 800, cost_per_sample = 50,
    n_sites = c(North = 10, South = 8, West = 12),
    fixed_cost_per_site = 5000,
    transport_cost_per_site = c(North = 1000, South = 1500, West = 800),
    budget = 300000
  )
  out <- capture.output(print(res))
  expect_true(any(grepl("North", out)))
  expect_true(any(grepl("South", out)))
  expect_true(any(grepl("under budget by 78,400", out, fixed = TRUE)))
  expect_true(any(grepl("more samples", out)))
})

test_that("EC-P3: over-budget print says over budget", {
  res <- estimate_cost(n = 800, cost_per_sample = 50, n_sites = 20,
                       fixed_cost_per_site = 5000, transport_cost_per_site = 2000,
                       budget = 100000)
  out <- capture.output(print(res))
  expect_true(any(grepl("over budget by 80,000", out, fixed = TRUE)))
})

test_that("EC-P4: large totals print grouped, not scientific", {
  res <- estimate_cost(n = 1e6, cost_per_sample = 10)
  out <- capture.output(print(res))
  expect_true(any(grepl("10,000,000", out, fixed = TRUE)))
  expect_false(any(grepl("e\\+0", out)))
})

test_that("EC-P5: exactly on budget -- no 'would cover' block", {
  res <- estimate_cost(n = 1000, cost_per_sample = 100, budget = 100000)
  expect_equal(res$budget_remaining, 0)
  expect_false(res$over_budget)
  out <- capture.output(print(res))
  expect_true(any(grepl("under budget by 0.00", out, fixed = TRUE)))
  expect_false(any(grepl("would cover", out)))
})
