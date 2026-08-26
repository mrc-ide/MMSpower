test_that("EC-1: basic SRS cost: total = n * cost_per_sample", {
  res <- estimate_cost(n = 100, cost_per_sample = 10)
  expect_equal(res$total_cost, 1000)
  expect_equal(res$sampling_cost, 1000)
  expect_equal(res$site_cost, 0)
  expect_null(res$n_sites)
})

test_that("EC-2: with sites: total = sampling + site fixed costs", {
  res <- estimate_cost(n = 323, cost_per_sample = 15, n_sites = 10, cost_per_site = 200)
  expect_equal(res$sampling_cost, 323 * 15)
  expect_equal(res$site_cost, 10 * 200)
  expect_equal(res$total_cost, 323 * 15 + 10 * 200)
})

test_that("EC-3: n_sites supplied but cost_per_site=0 → site_cost = 0", {
  res <- estimate_cost(n = 100, cost_per_sample = 5, n_sites = 10, cost_per_site = 0)
  expect_equal(res$site_cost, 0)
  expect_equal(res$total_cost, 500)
})

test_that("EC-4: cost_per_sample = 0 (free test) → sampling_cost = 0", {
  res <- estimate_cost(n = 500, cost_per_sample = 0)
  expect_equal(res$sampling_cost, 0)
  expect_equal(res$total_cost, 0)
})

test_that("EC-5: non-integer cost values work (e.g. $12.50 per sample)", {
  res <- estimate_cost(n = 200, cost_per_sample = 12.50)
  expect_equal(res$total_cost, 2500)
})

test_that("EC-6: n = 1 (minimum valid sample)", {
  res <- estimate_cost(n = 1, cost_per_sample = 99.99)
  expect_equal(res$total_cost, 99.99)
})

test_that("EC-7: n_sites = n (one sample per site)", {
  res <- estimate_cost(n = 5, cost_per_sample = 10, n_sites = 5, cost_per_site = 50)
  expect_equal(res$n_sites, 5)
  expect_equal(res$total_cost, 5 * 10 + 5 * 50)
})

test_that("EC-8: large n and costs produce correct arithmetic", {
  res <- estimate_cost(n = 1e6L, cost_per_sample = 0.01)
  expect_equal(res$total_cost, 10000)
})

test_that("EC-9: return list has all seven expected names", {
  res <- estimate_cost(n = 50, cost_per_sample = 5)
  expect_named(res, c("total_cost", "sampling_cost", "site_cost",
                      "n", "n_sites", "cost_per_sample", "cost_per_site"))
})

test_that("EC-10: result has class 'mms_cost'", {
  res <- estimate_cost(n = 50, cost_per_sample = 5)
  expect_s3_class(res, "mms_cost")
})

test_that("EC-11: input values are echoed back unchanged", {
  res <- estimate_cost(n = 77, cost_per_sample = 8.25,
                       n_sites = 7, cost_per_site = 150)
  expect_equal(res$n, 77)
  expect_equal(res$cost_per_sample, 8.25)
  expect_equal(res$n_sites, 7)
  expect_equal(res$cost_per_site, 150)
})

test_that("EC-12: total_cost = sampling_cost + site_cost always", {
  for (k in 1:5) {
    n_val    <- sample(10:500, 1)
    cps      <- runif(1, 0, 100)
    n_s      <- sample(1:min(n_val, 20), 1)
    cpsite   <- runif(1, 0, 500)
    res      <- estimate_cost(n = n_val, cost_per_sample = cps,
                              n_sites = n_s, cost_per_site = cpsite)
    expect_equal(res$total_cost, res$sampling_cost + res$site_cost,
                 tolerance = 1e-10)
  }
})

# ── Validation errors ─────────────────────────────────────────────────────────

test_that("EC-V1: n = 0 → error", {
  expect_error(estimate_cost(n = 0, cost_per_sample = 10), "positive integer")
})

test_that("EC-V2: n negative → error", {
  expect_error(estimate_cost(n = -5, cost_per_sample = 10), "positive integer")
})

test_that("EC-V3: n non-integer → error", {
  expect_error(estimate_cost(n = 10.5, cost_per_sample = 10), "positive integer")
})

test_that("EC-V4: n = Inf → error", {
  expect_error(estimate_cost(n = Inf, cost_per_sample = 10), "finite")
})

test_that("EC-V5: n as character → error", {
  expect_error(estimate_cost(n = "100", cost_per_sample = 10), "character")
})

test_that("EC-V6: n as vector → error", {
  expect_error(estimate_cost(n = c(100, 200), cost_per_sample = 10), "single")
})

test_that("EC-V7: cost_per_sample negative → error", {
  expect_error(estimate_cost(n = 100, cost_per_sample = -5), "non-negative")
})

test_that("EC-V8: cost_per_sample = NA → error", {
  expect_error(estimate_cost(n = 100, cost_per_sample = NA_real_), "finite")
})

test_that("EC-V9: n_sites = 0 → error", {
  expect_error(estimate_cost(n = 100, cost_per_sample = 5, n_sites = 0), "positive integer")
})

test_that("EC-V10: n_sites non-integer → error", {
  expect_error(estimate_cost(n = 100, cost_per_sample = 5, n_sites = 2.5), "positive integer")
})

test_that("EC-V11: n_sites > n → error", {
  expect_error(
    estimate_cost(n = 5, cost_per_sample = 10, n_sites = 10),
    "cannot exceed"
  )
})

test_that("EC-V12: cost_per_site negative → error", {
  expect_error(
    estimate_cost(n = 100, cost_per_sample = 5, n_sites = 10, cost_per_site = -1),
    "non-negative"
  )
})

test_that("EC-V13: cost_per_site as vector → error", {
  expect_error(
    estimate_cost(n = 100, cost_per_sample = 5, n_sites = 10,
                  cost_per_site = c(100, 200)),
    "single"
  )
})

test_that("EC-V14: cost_per_sample as vector → error", {
  expect_error(estimate_cost(n = 100, cost_per_sample = c(5, 10)), "single")
})

test_that("EC-V15: n = NA → error", {
  expect_error(estimate_cost(n = NA_integer_, cost_per_sample = 10), "finite")
})

# ── print method ──────────────────────────────────────────────────────────────

test_that("EC-P1: print returns the object invisibly", {
  res <- estimate_cost(n = 100, cost_per_sample = 10)
  expect_identical(print(res), res)
})

test_that("EC-P2: print output contains 'Total cost'", {
  res <- estimate_cost(n = 100, cost_per_sample = 10)
  expect_output(print(res), "Total cost")
})

test_that("EC-P3: print shows site line only when site_cost > 0", {
  res_no_site <- estimate_cost(n = 100, cost_per_sample = 10,
                               n_sites = 5, cost_per_site = 0)
  out_no_site <- capture.output(print(res_no_site))
  expect_false(any(grepl("Site fixed", out_no_site)))

  res_with_site <- estimate_cost(n = 100, cost_per_sample = 10,
                                 n_sites = 5, cost_per_site = 50)
  expect_output(print(res_with_site), "Site fixed")
})
