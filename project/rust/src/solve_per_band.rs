use rayon::prelude::*;
use crate::validate_build::CscRef;
use extendr_api::prelude::*;

// Sparse matrix-vector product: y = L * x
pub fn spmv(mat: &CscRef, x: &[f64], y: &mut [f64]) {
    for yi in y.iter_mut() {
        *yi = 0.0;
    }

    for col in 0..mat.ncols {
        let x_col = x[col];
        if x_col == 0.0 {
            continue;
        }
        let start = mat.col_ptrs[col];
        let end = mat.col_ptrs[col + 1];
        for idx in start..end {
            let row = mat.row_idx[idx];
            y[row] += mat.values[idx] * x_col;
        }
    }
}

// Conjugate Gradient solver for (I + alpha * L) x = b
pub fn solve_cg(mat: &CscRef, alpha: f64, b: &[f64], x: &mut [f64], max_iter: usize, tol: f64) {
    let n = b.len();
    let mut r = vec![0.0; n];
    let mut p = vec![0.0; n];
    let mut ap = vec![0.0; n];

    // Initial guess: x = b
    x.copy_from_slice(b);

    // r0 = b - (I + alpha L) x
    spmv(mat, x, &mut ap);

    let mut r_dot_r = 0.0;
    for i in 0..n {
        r[i] = b[i] - (x[i] + alpha * ap[i]);
        p[i] = r[i];
        r_dot_r += r[i] * r[i];
    }

    let tol_sq = tol * tol;
    if r_dot_r < tol_sq {
        return;
    }

    for _ in 0..max_iter {
        spmv(mat, &p, &mut ap);

        let mut p_dot_ap = 0.0;
        for i in 0..n {
            let a_pi = p[i] + alpha * ap[i];
            p_dot_ap += p[i] * a_pi;
            ap[i] = a_pi;
        }

        let alpha_step = if p_dot_ap > 0.0 {
            r_dot_r / p_dot_ap
        } else {
            0.0
        };

        let mut r_dot_r_new = 0.0;
        for i in 0..n {
            x[i] += alpha_step * p[i];
            r[i] -= alpha_step * ap[i];
            r_dot_r_new += r[i] * r[i];
        }

        if r_dot_r_new < tol_sq {
            break;
        }

        let beta = r_dot_r_new / r_dot_r;
        for i in 0..n {
            p[i] = r[i] + beta * p[i];
        }

        r_dot_r = r_dot_r_new;
    }
}

pub fn solve_per_band(signal: RMatrix<f64>, csc: &CscRef, alpha: f64) -> RMatrix<f64> {
    let n_pixels = signal.nrows();
    let n_bands = signal.ncols();

    let signal_slice = signal.as_real_slice().expect("Expected numeric matrix");

    let mut out_data = vec![0.0; n_pixels * n_bands];

    out_data
        .par_chunks_mut(n_pixels)
        .enumerate()
        .for_each(|(band_idx, out_band)| {
            let start = band_idx * n_pixels;
            let end = start + n_pixels;
            let b = &signal_slice[start..end];
            solve_cg(csc, alpha, b, out_band, 1000, 1e-6);
        });

    RMatrix::new_matrix(n_pixels, n_bands, |r, c| out_data[c * n_pixels + r])
}
