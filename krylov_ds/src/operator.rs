//! Matrix-free linear operator abstraction.
//!
//! Everything in this crate works against `apply()` alone, never against an
//! explicit dense/sparse matrix. That's what lets
//! [`crate::Arnoldi`] run on things like `nbsc`'s `2n x 2n` Bass-reduced
//! Hashimoto linearization without ever materializing a `2n x 2n` (or, worse,
//! `2m x 2m`) matrix.

/// A linear map on `R^dim`, applied one matrix-vector product at a time.
pub trait LinearOperator {
    /// Dimension of the space this operator acts on.
    fn dim(&self) -> usize;

    /// `y = A x`.
    fn apply(&self, x: &[f64]) -> Vec<f64>;
}

/// A plain dense matrix (row-major), for tests / small operators where
/// matrix-free isn't necessary.
pub struct DenseOperator {
    pub n: usize,
    pub rows: Vec<Vec<f64>>,
}

impl LinearOperator for DenseOperator {
    fn dim(&self) -> usize {
        self.n
    }

    fn apply(&self, x: &[f64]) -> Vec<f64> {
        self.rows
            .iter()
            .map(|row| row.iter().zip(x.iter()).map(|(a, b)| a * b).sum())
            .collect()
    }
}
