test_that("U-1: .wilson_ci matches R's own prop.test (no continuity correction)", {
  w <- .wilson_ci(0.3, 100, qnorm(0.975))
  ref <- prop.test(30, 100, correct = FALSE)$conf.int
  expect_equal(c(w$lower, w$upper), c(ref[1], ref[2]))   # 0.2189, 0.3958
})

test_that("U-2: .rogan_gladen undoes .apparent_prev", {
  # 0.3 * 0.9 + 0.7 * 0.05 = 0.305
  expect_equal(.apparent_prev(0.3, 0.9, 0.95), 0.305)
  expect_equal(.rogan_gladen(0.305, 0.9, 0.95), 0.3)
})
