#' @docType package
#' @usage NULL
#' @useDynLib rust, .registration = TRUE
NULL

#' @export
rust_graph_smooth <- function(signal, i_slot, p_slot, x_slot, nrow, ncol, alpha) .Call(
    "wrap__rust_graph_smooth",
    signal, i_slot, p_slot, x_slot, nrow, ncol, alpha
)
