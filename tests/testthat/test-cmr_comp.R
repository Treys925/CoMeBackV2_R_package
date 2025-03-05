num_rows <- 5
num_cols <- 24
data_ <- data.frame(matrix(runif(num_rows * num_cols, min = 0, max = 1), nrow = num_rows, ncol = num_cols))
names <- c("cg00000165", "cg11454459", "cg06365538", "cg20836156", "cg03667047", "cg14149503", "cg18093619", "cg07866476", 
           "cg03168497", "cg22306155", "cg18850089", "cg06215035", "cg14206848", "cg08009711", "cg09103083", "cg07138366", 
           "cg20099968", "cg15103180", "cg09149149", "cg27183030", "cg06942649", "cg20401521", "cg20359650", "cg27443332")
colnames(data_) <- names

test_that("make sure we return a matrix", {
  expect_true(is.matrix(cmr_comp(names, data_)))
})
#> Test passed 😸

test_that("the number of samples should not change", {
  expect_equal(nrow(cmr_comp(names, data_)), 5)
})
#> Test passed 🥇

