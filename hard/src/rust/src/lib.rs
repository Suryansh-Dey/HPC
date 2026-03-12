use extendr_api::prelude::*;
use faer::sparse::{SparseColMat, SymbolicSparseColMat};

/// Compute the degree vector (row sums) of an R dgCMatrix sparse matrix.
///
/// @param mat An R dgCMatrix object (S4, from the Matrix package).
/// @return A numeric vector of row sums.
/// @export
#[extendr]
fn dgc_degree_rust(mat: Robj) -> Doubles {
    let dim_robj = mat
        .get_attrib("Dim")
        .expect("dgCMatrix must have a 'Dim' slot");
    let dim_slice = dim_robj
        .as_integer_slice()
        .expect("'Dim' slot must be an integer vector");
    let nrows = dim_slice[0] as usize;
    let ncols = dim_slice[1] as usize;

    let i_robj = mat
        .get_attrib("i")
        .expect("dgCMatrix must have an 'i' slot");
    let i_slice: &[i32] = i_robj
        .as_integer_slice()
        .expect("'i' slot must be an integer vector");

    let p_robj = mat
        .get_attrib("p")
        .expect("dgCMatrix must have a 'p' slot");
    let p_slice: &[i32] = p_robj
        .as_integer_slice()
        .expect("'p' slot must be an integer vector");

    let x_robj = mat
        .get_attrib("x")
        .expect("dgCMatrix must have an 'x' slot");
    let x_slice: &[f64] = x_robj
        .as_real_slice()
        .expect("'x' slot must be a numeric vector");

    let col_ptr: Vec<usize> = p_slice.iter().map(|&v| v as usize).collect();
    let row_idx: Vec<usize> = i_slice.iter().map(|&v| v as usize).collect();
    let values: Vec<f64> = x_slice.to_vec();

    let symbolic = SymbolicSparseColMat::new_checked(nrows, ncols, col_ptr, None, row_idx);
    let sparse_mat = SparseColMat::<usize, f64>::new(symbolic, values);

    let mut result = Doubles::new(nrows);

    for i in 0..nrows {
        result[i] = Rfloat::from(0.0);
    }

    let (sym, val) = sparse_mat.parts();
    let (_, _, col_ptrs, _, row_indices) = sym.parts();
    for col in 0..ncols {
        let start = col_ptrs[col];
        let end = col_ptrs[col + 1];
        for idx in start..end {
            let row = row_indices[idx];
            let cur: f64 = result[row].inner();
            result[row] = Rfloat::from(cur + val[idx]);
        }
    }

    result
}

extendr_module! {
    mod dgcmat;
    fn dgc_degree_rust;
}
