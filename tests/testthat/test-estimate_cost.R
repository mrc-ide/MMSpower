test_that("one stage only", {
  res <- estimate_cost(n_sequenced = 323, cost_per_sequenced = 50)
  expect_equal(res$variable_cost, 16150)
  expect_equal(res$stage_cost, c(screened = 0, collected = 0, sequenced = 16150))
  expect_equal(res$site_cost, 0)
  expect_equal(res$total_cost, 16150)
  expect_null(res$budget_remaining)
})

test_that("three stages add up", {
  # 1950 * 2 + 780 * 15 + 700 * 40 = 3900 + 11700 + 28000
  res <- estimate_cost(n_screened = 1950, cost_per_screened = 2,
                       n_collected = 780, cost_per_collected = 15,
                       n_sequenced = 700, cost_per_sequenced = 40)
  expect_equal(res$stage_cost, c(screened = 3900, collected = 11700, sequenced = 28000))
  expect_equal(res$variable_cost, 43600)
  # any stage can be left out
  res2 <- estimate_cost(n_screened = 1950, cost_per_screened = 2,
                        n_sequenced = 700, cost_per_sequenced = 40)
  expect_equal(res2$stage_cost, c(screened = 3900, collected = 0, sequenced = 28000))
  # equal numbers at successive stages are allowed
  expect_no_error(estimate_cost(n_collected = 700, cost_per_collected = 15,
                                n_sequenced = 700, cost_per_sequenced = 40))
})

test_that("the same costs at every site", {
  # 600 * 50 + 12 * (5000 + 1000) = 30000 + 72000
  res <- estimate_cost(n_sequenced = 600, cost_per_sequenced = 50, n_sites = 12,
                       fixed_cost_per_site = 5000, transport_cost_per_site = 1000)
  expect_equal(res$site_cost, 72000)
  expect_equal(res$total_cost, 102000)
})

test_that("per-site costs and a budget", {
  # 50 * 10 + (100 + 200 + 300) + 3 * 20 = 500 + 600 + 60
  res <- estimate_cost(n_sequenced = 50, cost_per_sequenced = 10, n_sites = 3,
                       fixed_cost_per_site = c(100, 200, 300),
                       transport_cost_per_site = 20, budget = 1500)
  expect_equal(res$total_cost, 1160)
  expect_equal(res$budget_remaining, 340)
  over <- estimate_cost(n_sequenced = 50, cost_per_sequenced = 10, budget = 400)
  expect_equal(over$budget_remaining, -100)
})

test_that("a design's output can be costed", {
  des <- design_precision(0.15, 0.05, n_per_site = 30, icc = 0.05)
  res <- estimate_cost(n_sequenced = des$n, cost_per_sequenced = 50,
                       n_sites = des$n_sites, fixed_cost_per_site = 5000)
  expect_equal(res$total_cost, des$n * 50 + des$n_sites * 5000)
})

test_that("invalid inputs are rejected", {
  expect_error(estimate_cost(), "at least one stage")
  expect_error(estimate_cost(n_sequenced = 10), "together, or neither")
  expect_error(estimate_cost(cost_per_screened = 2), "together, or neither")
  expect_error(estimate_cost(n_sequenced = 10.5, cost_per_sequenced = 5), "whole number")
  expect_error(estimate_cost(n_sequenced = 10, cost_per_sequenced = -5),
               "cost_per_sequenced")
  expect_error(estimate_cost(n_screened = 100, cost_per_screened = 1,
                             n_collected = 200, cost_per_collected = 1), "subset")
  expect_error(estimate_cost(n_screened = 100, cost_per_screened = 1,
                             n_sequenced = 200, cost_per_sequenced = 1), "subset")
  expect_error(estimate_cost(n_sequenced = 10, cost_per_sequenced = 5,
                             fixed_cost_per_site = 100), "needs `n_sites`")
  expect_error(estimate_cost(n_sequenced = 10, cost_per_sequenced = 5, n_sites = 3,
                             fixed_cost_per_site = c(1, 2)), "one per site")
  expect_error(estimate_cost(n_sequenced = 10, cost_per_sequenced = 5, budget = 0),
               "budget")
})
