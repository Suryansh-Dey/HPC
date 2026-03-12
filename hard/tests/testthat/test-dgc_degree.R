library(testthat)
library(Matrix)
library(dgcmat)

test_that("3x3 sparse matrix", {
  m1 <- sparseMatrix(
    i = c(1, 2, 3, 1),
    j = c(1, 2, 3, 3),
    x = c(4.0, 5.0, 6.0, 7.0),
    dims = c(3, 3)
  )
  rust_result <- dgc_degree(m1)
  r_result    <- as.numeric(Matrix::rowSums(m1))
  expect_equal(rust_result, r_result)
})

test_that("5x5 identity", {
  m2 <- Diagonal(5)
  m2 <- as(m2, "dgCMatrix")
  rust_result <- dgc_degree(m2)
  r_result    <- as.numeric(Matrix::rowSums(m2))
  expect_equal(rust_result, r_result)
})

test_that("100x50 random sparse (300 nnz)", {
  set.seed(42)
  m3 <- rsparsematrix(100, 50, nnz = 300)
  rust_result <- dgc_degree(m3)
  r_result    <- as.numeric(Matrix::rowSums(m3))
  expect_equal(rust_result, r_result)
})

test_that("4x3 with empty row 2 and 4", {
  m4 <- sparseMatrix(
    i = c(1, 1, 3),
    j = c(1, 2, 3),
    x = c(10.0, 20.0, 30.0),
    dims = c(4, 3)
  )
  rust_result <- dgc_degree(m4)
  r_result    <- as.numeric(Matrix::rowSums(m4))
  expect_equal(rust_result, r_result)
})
