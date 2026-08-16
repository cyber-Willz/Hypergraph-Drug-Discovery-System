//! Matrix-free Arnoldi iteration.
//!
//! Builds an orthonormal Krylov basis `V = [v_0 .. v_{k-1}]` (each `v_i` a
//! length-`dim` vector) and the `k x k` upper-Hessenberg matrix `H` such
//! that `A V_k ~= V_k H` (modified Gram-Schmidt, with early termination on
//! Krylov-subspace breakdown, i.e. an invariant subspace was found before
//! reaching the requested dimension -- not an error, just fewer than `m`
//! basis vectors).

use crate::operator::LinearOperator;

#[derive(Debug, thiserror::Error)]
pub enum ArnoldiError {
    #[error("start vector length does not match operator dimension")]
    DimMismatch,
    #[error("start vector is (numerically) zero")]
    ZeroStartVector,
}

/// Output of one Arnoldi run: the orthonormal basis and the (square,
/// leading) Hessenberg projection of the operator onto that basis.
#[derive(Debug, Clone)]
pub struct ArnoldiResult {
    /// `k` orthonormal basis vectors, each of length `dim`.
    pub basis: Vec<Vec<f64>>,
    /// `k x k` upper-Hessenberg matrix, row-major (`hessenberg[i][j]`).
    pub hessenberg: Vec<Vec<f64>>,
    /// Actual Krylov dimension reached (`<= max_dim`; can be smaller if an
    /// invariant subspace was found early).
    pub k: usize,
    /// Dimension of the underlying operator's domain/range.
    pub dim: usize,
}

pub struct Arnoldi {
    pub max_dim: usize,
    pub tol: f64,
}

impl Arnoldi {
    pub fn new(max_dim: usize, tol: f64) -> Self {
        Self { max_dim: max_dim.max(1), tol }
    }

    pub fn run(
        &self,
        op: &dyn LinearOperator,
        v0: &[f64],
    ) -> Result<ArnoldiResult, ArnoldiError> {
        let dim = op.dim();
        if v0.len() != dim {
            return Err(ArnoldiError::DimMismatch);
        }
        let n0 = norm(v0);
        if n0 < self.tol {
            return Err(ArnoldiError::ZeroStartVector);
        }
        let m = self.max_dim.min(dim).max(1);

        let mut basis: Vec<Vec<f64>> = Vec::with_capacity(m + 1);
        basis.push(v0.iter().map(|x| x / n0).collect());

        // Column-major accumulation (h_cols[j][i] = H[i][j]); square k x k
        // Hessenberg is assembled from this once we know the final k.
        let mut h_cols: Vec<Vec<f64>> = Vec::with_capacity(m);
        let mut k_actual = 0usize;

        for j in 0..m {
            let mut w = op.apply(&basis[j]);
            let mut hcol = vec![0.0; m];
            for i in 0..=j {
                let hij = dot(&basis[i], &w);
                hcol[i] = hij;
                for (wc, vc) in w.iter_mut().zip(basis[i].iter()) {
                    *wc -= hij * vc;
                }
            }
            // Re-orthogonalize once against the existing basis (modified
            // Gram-Schmidt is only first-order stable; one extra pass keeps
            // orthogonality to close-to-machine precision even for the
            // larger Krylov dimensions this crate is used at).
            for i in 0..=j {
                let corr = dot(&basis[i], &w);
                hcol[i] += corr;
                for (wc, vc) in w.iter_mut().zip(basis[i].iter()) {
                    *wc -= corr * vc;
                }
            }

            let hnext = norm(&w);
            k_actual = j + 1;
            h_cols.push(hcol);
            if j + 1 < m {
                if hnext < self.tol {
                    // Invariant subspace: stop here, no breakdown error.
                    break;
                }
                let vnext: Vec<f64> = w.iter().map(|x| x / hnext).collect();
                h_cols[j][j + 1] = hnext;
                basis.push(vnext);
            }
        }

        basis.truncate(k_actual);
        let mut hessenberg = vec![vec![0.0; k_actual]; k_actual];
        for (j, col) in h_cols.iter().enumerate().take(k_actual) {
            for i in 0..k_actual {
                hessenberg[i][j] = col.get(i).copied().unwrap_or(0.0);
            }
        }

        Ok(ArnoldiResult { basis, hessenberg, k: k_actual, dim })
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::operator::DenseOperator;

    #[test]
    fn arnoldi_reproduces_full_hessenberg_for_symmetric_matrix() {
        // Symmetric 3x3, so H should come out (numerically) symmetric
        // tridiagonal and its trace should equal the matrix trace.
        let op = DenseOperator {
            n: 3,
            rows: vec![
                vec![2.0, 1.0, 0.0],
                vec![1.0, 2.0, 1.0],
                vec![0.0, 1.0, 2.0],
            ],
        };
        let arnoldi = Arnoldi::new(3, 1e-12);
        let v0 = vec![1.0, 0.0, 0.0];
        let result = arnoldi.run(&op, &v0).unwrap();
        assert_eq!(result.k, 3);
        let trace: f64 = (0..3).map(|i| result.hessenberg[i][i]).sum();
        assert!((trace - 6.0).abs() < 1e-8);
    }
}
