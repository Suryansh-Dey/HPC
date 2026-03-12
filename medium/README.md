# Medium Task: Dense Matrix Operations via Rust

This task demonstrates how to pass a dense matrix from R to Rust using `extendr` and compute row sums using the `faer` linear algebra library.

## How it Works

1. The `rowSum` package uses `extendr` to create a bridge between R and Rust.
2. In Rust, the `RMatrix<f64>` is converted into a `faer::Mat<f64>`. 
   R uses column-major ordering, so the data is indexed as `j * nrow + i` when constructing the `faer` matrix.
3. The row sums are computed by iterating over the `faer` matrix and returned as a standard Rust `Vec<f64>`, which `extendr` automatically converts back to an R numeric vector.
4. Using `faer` allows for highly optimized linear algebra operations that can outperform pure R for large, dense datasets.

## How to Run

### 1. Build and Install the Package

From the `medium/` directory, run:

```bash
R CMD INSTALL .
```

### 2. Run Tests

To verify the implementation, you can run the included test suite:

```bash
Rscript -e "testthat::test_local('.')"
```

### 3. Usage Example

You can also use the function directly in R:

```R
library(rowSum)
m <- matrix(1:12, nrow = 3, ncol = 4)
row_sums_faer(m)
# Should match rowSums(m)
```
