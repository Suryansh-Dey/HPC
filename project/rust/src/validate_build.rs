// Raw CSC layout — borrows directly from R memory, zero allocation for x_slot
pub struct CscRef<'a> {
    pub nrows: usize,
    pub ncols: usize,
    pub col_ptrs: Vec<usize>,
    pub row_idx: Vec<usize>,
    pub values: &'a [f64],
}
impl<'a> CscRef<'a> {
    pub fn new_checked(
        nrow: usize,
        ncol: usize,
        p_slot: &[i32],
        i_slot: &[i32],
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

        // i32 -> usize conversion is unavoidable since R uses 32-bit indices
        let col_ptrs: Vec<usize> = p_slot
            .iter()
            .map(|&v| {
                assert!(v >= 0, "negative column pointer: {}", v);
                v as usize
            })
            .collect();

        // col_ptrs must be monotonically non-decreasing
        for j in 0..ncol {
            assert!(
                col_ptrs[j] <= col_ptrs[j + 1],
                "col_ptrs not monotonic at column {}: {} > {}",
                j,
                col_ptrs[j],
                col_ptrs[j + 1]
            );
        }

        let nnz = col_ptrs[ncol];
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

        let row_idx: Vec<usize> = i_slot
            .iter()
            .map(|&v| {
                assert!(
                    v >= 0 && (v as usize) < nrow,
                    "row index out of bounds: {} (nrows={})",
                    v,
                    nrow
                );
                v as usize
            })
            .collect();

        // Verify row indices are strictly sorted and unique within each column
        // (This is exactly what faer's `check_row_idx` enforces)
        for j in 0..ncol {
            let start = col_ptrs[j];
            let end = col_ptrs[j + 1];
            if start < end {
                let mut prev = row_idx[start];
                for idx in (start + 1)..end {
                    let curr = row_idx[idx];
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
            col_ptrs,
            row_idx,
            values: x_slot, // direct pointer borrow, no .to_vec()
        }
    }
}
