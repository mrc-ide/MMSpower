test_that("variable cost only", {
  res <- estimate_cost(n = 323, cost_per_sample = 50)
  expect_equal(res$variable_cost, 16150)
  expect_equal(res$site_cost, 0)
  expect_equal(res$total_cost, 16150)
  expect_null(res$budget_remaining)
})

test_that("the same costs at every site", {
  # 600 * 50 + 12 * (5000 + 1000) = 30000 + 72000
  res <- estimate_cost(n = 600, cost_per_sample = 50, n_sites = 12,
                       fixed_cost_per_site = 5000, transport_cost_per_site = 1000)
  expect_equal(res$site_cost, 72000)
  expect_equal(res$total_cost, 102000)
})

test_that("per-site costs and a budget", {
  # 50 * 10 + (100 + 200 + 300) + 3 * 20 = 500 + 600 + 60
  res <- estimate_cost(n = 50, cost_per_sample = 10, n_sites = 3,
                       fixed_cost_per_site = c(100, 200, 300),
                       transport_cost_per_site = 20, budget = 1500)
  expect_equal(res$total_cost, 1160)
  expect_equal(res$budget_remaining, 340)
  over <- estimate_cost(n = 50, cost_per_sample = 10, budget = 400)
  expect_equal(over$budget_remaining, -100)
})

test_that("a design's output can be costed", {
  des <- design_precision(0.15, 0.05, n_per_site = 30, icc = 0.05)
  res <- estimate_cost(n = des$n, cost_per_sample = 50, n_sites = des$n_sites,
                       fixed_cost_per_site = 5000)
  expect_equal(res$total_cost, des$n * 50 + des$n_sites * 5000)
})

test_that("invalid inputs are rejected", {
  expect_error(estimate_cost(n = 10.5, cost_per_sample = 5), "whole number")
  expect_error(estimate_cost(n = 10, cost_per_sample = -5), "cost_per_sample")
  expect_error(estimate_cost(n = 10, cost_per_sample = 5, fixed_cost_per_site = 100),
               "needs `n_sites`")
  expect_error(estimate_cost(n = 10, cost_per_sample = 5, n_sites = 3,
                             fixed_cost_per_site = c(1, 2)), "one per site")
  expect_error(estimate_cost(n = 10, cost_per_sample = 5, budget = 0), "budget")
})
