library(testthat)
library(rust)
library(hyperSpec)
library(Matrix)

test_that("graphSmooth computes successfully", {
  # 1. Create a dummy hyperSpec object
  spc <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10) # 100 pixels, 10 bands
  hy <- new("hyperSpec", spc = spc)

  # 2. Smooth it
  hy_smooth <- graphSmooth(hy, alpha = 0.5, neighbors = 4L)

  # 3. Check the dimensions and layout
  expect_equal(dim(hy_smooth), dim(hy))

  # 4. Check if the variance actually reduced
  orig_var <- mean(apply(hy@data$spc, 2, var))
  smooth_var <- mean(apply(hy_smooth@data$spc, 2, var))
  expect_true(smooth_var < orig_var)
})
