dat = read.table(test_path("testdata.tsv"))

test_that("make sure we return a list of CMRs", {
  expect_true(is.list(cmr(dat)))
})
#> Test passed 😸

test_that("We should have a set of CMRs for each chromosome", {
  expect_equal(length(cmr(dat)), 24)
})
#> Test passed 🥇

