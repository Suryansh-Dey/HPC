use extendr_api::prelude::*;
use faer::Mat;

/// Compute row sums of a dense matrix using the faer linear algebra library.
///
/// @param matrix A numeric matrix passed from R.
/// @return A numeric vector of row sums.
/// @export
#[extendr]
fn row_sums_faer(matrix: RMatrix<f64>) -> Vec<f64> {
    let nrow = matrix.nrows();
    let ncol = matrix.ncols();

    let data = matrix.as_real_slice().expect("expected numeric matrix");

    let mat = Mat::<f64>::from_fn(nrow, ncol, |i, j| {
        data[j * nrow + i]
    });

    let mut row_sums = vec![0.0_f64; nrow];
    for i in 0..nrow {
        let mut sum = 0.0_f64;
        for j in 0..ncol {
            sum += mat[(i, j)];
        }
        row_sums[i] = sum;
    }

    row_sums
}

extendr_module! {
    mod rowSum;
    fn row_sums_faer;
}
