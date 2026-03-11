library(hyperSpec)
library(Matrix)

data(laser)
cat("Laser dataset loaded successfully.\n")
spc <- laser@data$spc
wl_vec <- wl(laser)
cat(sprintf(
  "  Dimensions: %d spectra × %d wavelengths\n",
  nrow(spc), ncol(spc)
))
cat(sprintf(
  "  Wavelength range: %.1f – %.1f\n",
  min(wl_vec), max(wl_vec)
))

n_spectra <- nrow(spc)

grid_rows <- floor(sqrt(n_spectra))
while (n_spectra %% grid_rows != 0) {
  grid_rows <- grid_rows - 1
}
grid_cols <- n_spectra %/% grid_rows
cat(sprintf(
  "\nSynthetic spatial grid: %d × %d = %d pixels\n",
  grid_rows, grid_cols, grid_rows * grid_cols
))

pixel_coords <- expand.grid(
  row = seq_len(grid_rows),
  col = seq_len(grid_cols)
)

n <- nrow(pixel_coords)
cat(sprintf("\nBuilding 4-connectivity adjacency matrix for %d pixels...\n", n))

# Helper: convert (r, c) → linear index (1-based, column-major order)
rc_to_idx <- function(r, c) (c - 1) * grid_rows + r

edges_i <- integer(0)
edges_j <- integer(0)

for (r in seq_len(grid_rows)) {
  for (c in seq_len(grid_cols)) {
    idx <- rc_to_idx(r, c)
    # Right neighbor
    if (c < grid_cols) {
      nb <- rc_to_idx(r, c + 1)
      edges_i <- c(edges_i, idx, nb)
      edges_j <- c(edges_j, nb, idx)
    }
    # Down neighbor
    if (r < grid_rows) {
      nb <- rc_to_idx(r + 1, c)
      edges_i <- c(edges_i, idx, nb)
      edges_j <- c(edges_j, nb, idx)
    }
  }
}

a <- sparseMatrix(i = edges_i, j = edges_j, x = 1, dims = c(n, n))
cat(sprintf(
  "  Adjacency matrix: %d × %d, %d non-zeros\n",
  nrow(a), ncol(a), nnzero(a)
))

degree <- rowSums(a)
d <- Diagonal(x = degree)
l <- d - a # combinatorial Laplacian

cat(sprintf("  Degree range: %d – %d\n", min(degree), max(degree)))
cat(sprintf(
  "  Laplacian: %d × %d, %d non-zeros\n",
  nrow(l), ncol(l), nnzero(l)
))

d_inv_sqrt <- Diagonal(x = 1 / sqrt(degree))
l_norm <- d_inv_sqrt %*% l %*% d_inv_sqrt

cat("  Normalized Laplacian computed.\n")

wl_idx <- ncol(spc) %/% 2 # middle wavelength
b <- as.numeric(spc[, wl_idx]) # signal vector (one value per pixel)
alpha <- 0.5 # smoothing strength

cat(sprintf(
  "\nSpatial smoothing at wavelength index %d (%.1f nm), α = %.2f\n",
  wl_idx, wl_vec[wl_idx], alpha
))

# Solve (I + αL) x = b using base R's sparse solve
i_mat <- Diagonal(n)
system_mat <- i_mat + alpha * l
x_smooth <- as.numeric(solve(system_mat, b))

cat(sprintf("  Original signal: mean=%.2f, sd=%.4f\n", mean(b), sd(b)))
cat(sprintf(
  "  Smoothed signal: mean=%.2f, sd=%.4f\n",
  mean(x_smooth), sd(x_smooth)
))
cat(sprintf(
  "  Smoothing reduced variance by %.1f%%\n",
  100 * (1 - var(x_smooth) / var(b))
))

output_file <- "results.png"
cat(sprintf("\nSaving plots to %s ...\n", output_file))

png(output_file, width = 1400, height = 1000, res = 150)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

sub_n <- min(n, 40)
image(as.matrix(a[1:sub_n, 1:sub_n]),
  main = sprintf("Adjacency Matrix (first %d pixels)", sub_n),
  xlab = "Pixel index", ylab = "Pixel index",
  col = c("white", "steelblue"),
  axes = FALSE
)
axis(1,
  at = seq(0, 1, length.out = 5),
  labels = round(seq(1, sub_n, length.out = 5))
)
axis(2,
  at = seq(0, 1, length.out = 5),
  labels = round(seq(1, sub_n, length.out = 5))
)

barplot(table(degree),
  main = "Degree Distribution (4-connectivity)",
  xlab = "Degree", ylab = "Count",
  col = "coral", border = "darkred"
)

signal_mat_orig <- matrix(b, nrow = grid_rows, ncol = grid_cols)
image(signal_mat_orig,
  main = sprintf("Original Signal (wl = %.1f nm)", wl(laser)[wl_idx]),
  xlab = "Column", ylab = "Row",
  col = hcl.colors(64, "YlOrRd"),
  axes = FALSE
)

signal_mat_smooth <- matrix(x_smooth, nrow = grid_rows, ncol = grid_cols)
image(signal_mat_smooth,
  main = sprintf("Smoothed Signal (α = %.2f)", alpha),
  xlab = "Column", ylab = "Row",
  col = hcl.colors(64, "YlOrRd"),
  axes = FALSE
)

dev.off()
cat("Plots saved")
