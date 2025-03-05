test_that("Check for N > 25", {
  expect_no_error(ccutN_f(30))
})
#> Test passed 🥇

test_that("Sample fewer than 25 gives cutoff of 0.5", {
  expect_equal(ccutN_f(24), 0.5)
})
#> Test passed 🥇

test_that("N = 0 returns NA", {
  expect_equal(ccutN_f(0), NA)
})
#> Test passed 😸
