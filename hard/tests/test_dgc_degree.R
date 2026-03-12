library(Matrix)
library(dgcmat)

# Test 1: simple 3x3 matrix
m1 <- sparseMatrix(
  i = c(1, 2, 3, 1),
  j = c(1, 2, 3, 3),
  x = c(4.0, 5.0, 6.0, 7.0),
  dims = c(3, 3)
)
cat("Test 1: 3x3 sparse matrix\n")
rust_result <- dgc_degree(m1)
r_result    <- as.numeric(rowSums(m1))
cat("  Rust row sums:", rust_result, "\n")
cat("  R    row sums:", r_result, "\n")
stopifnot(all.equal(rust_result, r_result))
cat("  PASS\n\n")

# Test 2: identity matrix
m2 <- Diagonal(5)
m2 <- as(m2, "dgCMatrix")
cat("Test 2: 5x5 identity\n")
rust_result <- dgc_degree(m2)
r_result    <- as.numeric(rowSums(m2))
cat("  Rust row sums:", rust_result, "\n")
cat("  R    row sums:", r_result, "\n")
stopifnot(all.equal(rust_result, r_result))
cat("  PASS\n\n")

# Test 3: larger random sparse matrix
set.seed(42)
m3 <- rsparsematrix(100, 50, nnz = 300)
cat("Test 3: 100x50 random sparse (300 nnz)\n")
rust_result <- dgc_degree(m3)
r_result    <- as.numeric(rowSums(m3))
stopifnot(all.equal(rust_result, r_result))
cat("  PASS\n\n")

# Test 4: matrix with empty rows
m4 <- sparseMatrix(
  i = c(1, 1, 3),
  j = c(1, 2, 3),
  x = c(10.0, 20.0, 30.0),
  dims = c(4, 3)
)
cat("Test 4: 4x3 with empty row 2 and 4\n")
rust_result <- dgc_degree(m4)
r_result    <- as.numeric(rowSums(m4))
cat("  Rust row sums:", rust_result, "\n")
cat("  R    row sums:", r_result, "\n")
stopifnot(all.equal(rust_result, r_result))
cat("  PASS\n\n")

cat("=== All tests passed! ===\n")
