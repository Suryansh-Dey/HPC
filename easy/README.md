# Easy Task: Graph Laplacian Smoothing in R

This task demonstrates the construction and application of a spatial graph Laplacian for hyperspectral data smoothing using R.

## How it Works

1. The script loads the `laser` dataset from `hyperSpec`.
2. It constructs a synthetic 2D spatial grid corresponding to the number of spectra.
3. A 4-connectivity adjacency matrix is built using the `Matrix` package to represent spatial relationships between pixels.
4. It computes the combinatorial Laplacian ($L = D - A$) and the normalized Laplacian.
5. It performs spatial smoothing at a specific wavelength by solving the system $(I + \alpha L)x = b$, where $b$ is the original signal and $\alpha$ is the smoothing strength.
6. The result is saved to `results.png`, showing the adjacency pattern, degree distribution, and comparing the original vs. smoothed signal.

## How to Run

Ensure you have the `hyperSpec` and `Matrix` packages installed in R. Then, run the script from the command line:

```bash
Rscript graph_laplacian.R
```

After execution, check `results.png` for the generated plots.
