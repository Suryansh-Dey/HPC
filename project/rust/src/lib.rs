mod solve_per_band;
mod validate_build;
use extendr_api::prelude::*;
use solve_per_band::solve_per_band;
use validate_build::CscRef;

#[extendr]
fn rust_graph_smooth(
    signal: RMatrix<f64>,
    i_slot: &[i32],
    p_slot: &[i32],
    x_slot: &[f64],
    nrow: i32,
    ncol: i32,
    alpha: f64,
) -> RMatrix<f64> {
    assert_eq!(p_slot.len(), ncol as usize + 1);
    let csc = CscRef::new_checked(nrow as usize, ncol as usize, p_slot, i_slot, x_slot);
    solve_per_band(signal, &csc, alpha)
}

extendr_module! {
    mod rust;
    fn rust_graph_smooth;
}

#[cfg(test)]
mod tests;
