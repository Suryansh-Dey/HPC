#' Spatial Smoothing for hyperSpec Objects
#'
#' @param hy hyperSpec object
#' @param alpha smoothing parameter
#' @param neighbors adjacency connectivity
#' @export
setGeneric("graphSmooth",
           function(hy, alpha, neighbors = 4L) standardGeneric("graphSmooth"))

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(hyperSpec))

# Helper to build combinatorial Laplacian
build_pixel_laplacian <- function(hy, neighbors = 4L) {
  spc <- hy@data$spc
  n_spectra <- nrow(spc)
  grid_rows <- floor(sqrt(n_spectra))

  while (n_spectra %% grid_rows != 0) {
    grid_rows <- grid_rows - 1
  }
  grid_cols <- n_spectra %/% grid_rows

  rc_to_idx <- function(r, c) (c - 1) * grid_rows + r
  edges_i <- integer(0)
  edges_j <- integer(0)

  for (r in seq_len(grid_rows)) {
    for (c in seq_len(grid_cols)) {
      idx <- rc_to_idx(r, c)
      if (c < grid_cols) {
        nb <- rc_to_idx(r, c + 1)
        edges_i <- c(edges_i, idx, nb)
        edges_j <- c(edges_j, nb, idx)
      }
      if (r < grid_rows) {
        nb <- rc_to_idx(r + 1, c)
        edges_i <- c(edges_i, idx, nb)
        edges_j <- c(edges_j, nb, idx)
      }
    }
  }

  a <- sparseMatrix(i = edges_i, j = edges_j, x = 1,
                    dims = c(n_spectra, n_spectra))
  degree <- rowSums(a)
  d <- Diagonal(x = degree)
  l <- d - a
  return(l)
}

#' @export
setMethod("graphSmooth", signature(hy = "hyperSpec"),
          function(hy, alpha, neighbors = 4L) {
  L <- build_pixel_laplacian(hy, neighbors)
  i_slot <- L@i
  p_slot <- L@p
  x_slot <- L@x

            signal <- hy@data$spc

  # R calls the Rust extendr bridge
  smoothed_signal <- rust_graph_smooth(
    signal = signal,
    i_slot = i_slot,
    p_slot = p_slot,
    x_slot = x_slot,
    nrow   = nrow(L),
    ncol   = ncol(L),
    alpha  = alpha
  )

  # Return a new hyperSpec object with smoothed data
  res <- hy
  res@data$spc <- smoothed_signal
  return(res)
})
