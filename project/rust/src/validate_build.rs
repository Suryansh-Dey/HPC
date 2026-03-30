// Raw CSC layout — borrows directly from R memory, total zero-copy for all slots (p, i, x)
pub struct CscRef<'a> {
    pub nrows: usize,
    pub ncols: usize,
    pub col_ptrs: &'a [i32], // zero copy
    pub row_idx: &'a [i32],  // zero copy
    pub values: &'a [f64],   // zero copy
}
impl<'a> CscRef<'a> {
    pub fn new_checked(
        nrow: usize,
        ncol: usize,
        p_slot: &'a [i32],
        i_slot: &'a [i32],
        x_slot: &'a [f64],
    ) -> CscRef<'a> {
        // Validate CSC structure (mirrors faer::SymbolicSparseColMat::new_checked)
        assert_eq!(
            p_slot.len(),
            ncol + 1,
            "p slot length must be ncols + 1, got {} for {} columns",
            p_slot.len(),
            ncol
        );

        // col_ptrs must be monotonically non-decreasing and non-negative
        for j in 0..ncol {
            assert!(
                p_slot[j] >= 0,
                "negative column pointer at {}: {}",
                j,
                p_slot[j]
            );
            assert!(
                p_slot[j] <= p_slot[j + 1],
                "col_ptrs not monotonic at column {}: {} > {}",
                j,
                p_slot[j],
                p_slot[j + 1]
            );
        }
        assert!(
            p_slot[ncol] >= 0,
            "negative column pointer at {}: {}",
            ncol,
            p_slot[ncol]
        );

        let nnz = p_slot[ncol] as usize;
        assert_eq!(
            i_slot.len(),
            nnz,
            "i slot length ({}) does not match nnz from p slot ({})",
            i_slot.len(),
            nnz
        );
        assert_eq!(
            x_slot.len(),
            nnz,
            "x slot length ({}) does not match nnz from p slot ({})",
            x_slot.len(),
            nnz
        );

        // Validate row indices capacity and strictly ascending order explicitly
        for j in 0..ncol {
            let start = p_slot[j] as usize;
            let end = p_slot[j + 1] as usize;
            if start < end {
                let mut prev = i_slot[start];
                assert!(
                    prev >= 0 && (prev as usize) < nrow,
                    "row index out of bounds: {} (nrows={})",
                    prev,
                    nrow
                );
                for idx in (start + 1)..end {
                    let curr = i_slot[idx];
                    assert!(
                        curr >= 0 && (curr as usize) < nrow,
                        "row index out of bounds: {} (nrows={})",
                        curr,
                        nrow
                    );
                    assert!(
                        prev < curr,
                        "row indices not strictly sorted in column {}: {} >= {}",
                        j,
                        prev,
                        curr
                    );
                    prev = curr;
                }
            }
        }

        CscRef {
            nrows: nrow,
            ncols: ncol,
            col_ptrs: p_slot,
            row_idx: i_slot,
            values: x_slot,
        }
    }
}
