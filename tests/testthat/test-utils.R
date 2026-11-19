test_that("U-2: .rogan_gladen undoes .apparent_prev", {
  # 0.3 * 0.9 + 0.7 * 0.05 = 0.305
  expect_equal(.apparent_prev(0.3, 0.9, 0.95), 0.305)
  expect_equal(.rogan_gladen(0.305, 0.9, 0.95), 0.3)
})
