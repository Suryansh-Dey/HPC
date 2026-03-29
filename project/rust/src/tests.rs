use crate::solve_per_band::*;
use crate::validate_build::*;

// Build the 4x4 Laplacian for a 2x2 grid (4-connectivity):
//   L = [ 2 -1 -1  0]
//       [-1  2  0 -1]
//       [-1  0  2 -1]
//       [ 0 -1 -1  2]
//
// CSC (column-major):
//   col 0: rows {0,1,2}, vals {2,-1,-1}
//   col 1: rows {0,1,3}, vals {-1,2,-1}
//   col 2: rows {0,2,3}, vals {-1,2,-1}
//   col 3: rows {1,2,3}, vals {-1,-1,2}
fn laplacian_2x2() -> CscRef<'static> {
    let p: &[i32] = &[0, 3, 6, 9, 12];
    let i: &[i32] = &[0, 1, 2, 0, 1, 3, 0, 2, 3, 1, 2, 3];
    let x: &[f64] = &[
        2.0, -1.0, -1.0, -1.0, 2.0, -1.0, -1.0, 2.0, -1.0, -1.0, -1.0, 2.0,
    ];
    CscRef::new_checked(4, 4, p, i, x)
}

#[test]
fn test_build_csc_dimensions() {
    let mat = laplacian_2x2();
    assert_eq!(mat.nrows, 4);
    assert_eq!(mat.ncols, 4);
    assert_eq!(mat.values.len(), 12);
    assert_eq!(mat.col_ptrs.len(), 5);
    assert_eq!(mat.row_idx.len(), 12);
}

#[test]
fn test_spmv_laplacian_nullspace() {
    // L * ones = 0 (row sums of a Laplacian are zero)
    let mat = laplacian_2x2();
    let x = vec![1.0, 1.0, 1.0, 1.0];
    let mut y = vec![0.0; 4];
    spmv(&mat, &x, &mut y);
    for val in &y {
        assert!(val.abs() < 1e-12, "L * ones should be zero, got {}", val);
    }
}

#[test]
fn test_spmv_first_column() {
    // L * e0 = column 0 of L = [2, -1, -1, 0]
    let mat = laplacian_2x2();
    let x = vec![1.0, 0.0, 0.0, 0.0];
    let mut y = vec![0.0; 4];
    spmv(&mat, &x, &mut y);
    let expected = vec![2.0, -1.0, -1.0, 0.0];
    for i in 0..4 {
        assert!(
            (y[i] - expected[i]).abs() < 1e-12,
            "spmv mismatch at {}: {} vs {}",
            i,
            y[i],
            expected[i]
        );
    }
}

#[test]
fn test_cg_residual_accuracy() {
    // Solve (I + 0.5 L) x = b and verify the residual is tiny
    let mat = laplacian_2x2();
    let alpha = 0.5;
    let b = vec![5.0, 4.8, 5.1, 4.7];
    let mut x = vec![0.0; 4];

    solve_cg(&mat, alpha, &b, &mut x, 1000, 1e-10);

    // Verify: r = b - (I + alpha*L)*x
    let mut lx = vec![0.0; 4];
    spmv(&mat, &x, &mut lx);
    let mut max_residual = 0.0_f64;
    for i in 0..4 {
        let ax_i = x[i] + alpha * lx[i];
        let residual = (b[i] - ax_i).abs();
        if residual > max_residual {
            max_residual = residual;
        }
    }
    assert!(
        max_residual < 1e-8,
        "CG residual too large: {}",
        max_residual
    );
}

#[test]
fn test_cg_preserves_dc() {
    // Constant signal passes through unchanged: L * ones = 0 => x = b
    let mat = laplacian_2x2();
    let b = vec![3.0, 3.0, 3.0, 3.0];
    let mut x = vec![0.0; 4];

    solve_cg(&mat, 1.0, &b, &mut x, 1000, 1e-12);

    for i in 0..4 {
        assert!(
            (x[i] - 3.0).abs() < 1e-8,
            "DC component not preserved at {}: {}",
            i,
            x[i]
        );
    }
}

#[test]
fn test_cg_reduces_variance() {
    // Checkerboard noise should have reduced variance after smoothing
    let mat = laplacian_2x2();
    let b = vec![10.0, 1.0, 1.0, 10.0];
    let mut x = vec![0.0; 4];

    solve_cg(&mat, 2.0, &b, &mut x, 1000, 1e-10);

    let mean_b: f64 = b.iter().sum::<f64>() / 4.0;
    let var_b: f64 = b.iter().map(|v| (v - mean_b).powi(2)).sum::<f64>() / 4.0;
    let mean_x: f64 = x.iter().sum::<f64>() / 4.0;
    let var_x: f64 = x.iter().map(|v| (v - mean_x).powi(2)).sum::<f64>() / 4.0;

    assert!(
        var_x < var_b,
        "Smoothing should reduce variance: var_b={}, var_x={}",
        var_b,
        var_x
    );
}

#[test]
fn test_zero_copy_x_slot() {
    // Verify that CscRef.values points to the same memory as the input slice
    let x_data: Vec<f64> = vec![1.0, 2.0, 3.0];
    let p: Vec<i32> = vec![0, 1, 2, 3];
    let i: Vec<i32> = vec![0, 1, 2];

    let csc = CscRef::new_checked(3, 3, &p, &i, &x_data);

    // Same pointer address — proves zero copy
    assert!(
        std::ptr::eq(csc.values.as_ptr(), x_data.as_ptr()),
        "x_slot should be borrowed, not copied"
    );
}
