test_that("row_sums_faer matches base R rowSums on a 3x4 matrix", {
  m <- matrix(1:12, nrow = 3, ncol = 4)
  expect_equal(row_sums_faer(m), rowSums(m))
})

test_that("row_sums_faer works on a single-row matrix", {
  m <- matrix(c(1.0, 2.0, 3.0), nrow = 1)
  expect_equal(row_sums_faer(m), rowSums(m))
})

test_that("row_sums_faer works on a single-column matrix", {
  m <- matrix(c(10.0, 20.0, 30.0), ncol = 1)
  expect_equal(row_sums_faer(m), rowSums(m))
})

test_that("row_sums_faer works on an identity matrix", {
  m <- diag(5)
  expect_equal(row_sums_faer(m), rowSums(m))
})
