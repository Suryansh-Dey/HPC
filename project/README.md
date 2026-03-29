# hyperSpec-HPC Project Prototype

This directory contains the final architectural prototype for the Fast Sparse Spectral Kernels project. It bridges R's `dgCMatrix` object to Rust's high-performance `faer` sparse matrix interface via `extendr`, efficiently solving the spatial graph-smoothing objective:

$$(I + \alpha L) x = b$$

## Usage

Build and test using standard R commands:

```bash
R CMD INSTALL .
Rscript -e "testthat::test_local('.')"
```

Example in R:
```R
library(hyperSpec)
library(rust) # This module

data(laser)
# Apply a smoothing alpha = 0.5 to spatial grid neighbors
smoothed_laser <- graphSmooth(laser, alpha = 0.5, neighbors = 4)
```
