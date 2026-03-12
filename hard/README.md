# Hard Task: Sparse Matrix Row Sums via Rust

This task demonstrates high-performance sparse matrix handling by manually extracting slots from an R `dgCMatrix` and using the `faer` sparse module for computations.

## How it Works

1. Since `extendr` does not provide an automatic wrapper for sparse matrices, this implementation manually extracts the `i`, `p`, and `x` slots (standard Compressed Sparse Column format) using R attributes.
2. The numeric vectors from R are wrapped as Rust slices to avoid copying data where possible.
3. A `faer::sparse::SparseColMat` is constructed using these slices. This requires building a symbolic structure first.
4. The implementation iterates over the columns of the sparse matrix, adding values to the corresponding row indices. This is an $O(NNZ)$ operation.
5. The implementation validates that the slots exist and are of the correct type, providing clear errors if an incompatible object is passed.

## How to Run

### 1. Build and Install the Package

From the `hard/` directory, run:

```bash
R CMD INSTALL .
```

### 2. Run Tests

Verify the implementation against base R `rowSums`:

```bash
Rscript -e "testthat::test_local('.')"
```

### 3. Usage Example

```R
library(Matrix)
library(dgcmat)

# Create a sparse matrix
m <- sparseMatrix(i = 1:3, j = 1:3, x = 10, dims = c(5, 5))

# Compute row sums using Rust
dgc_degree(m)
```
