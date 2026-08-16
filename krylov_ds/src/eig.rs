//! Extract real Ritz pairs (eigenvalue/eigenvector approximations of the
//! original operator) from an [`ArnoldiResult`]'s small, dense Hessenberg
//! projection.
//!
//! Strategy: compute the real Schur form of the `k x k` Hessenberg matrix
//! `H` (via `nalgebra`, which is happy to do this for a general, possibly
//! non-symmetric, real matrix), read off the real eigenvalues from the
//! (quasi-triangular) Schur form's diagonal, then recover each real
//! eigenvalue's eigenvector by a few steps of shifted inverse iteration on
//! `H` directly. The eigenvector in the Krylov basis is then lifted back
//! into the full operator's domain via `basis` (`V y`).
//!
//! Complex-conjugate eigenvalue pairs (which do occur for non-symmetric
//! operators like a non-backtracking-walk linearization, off the dominant
//! Perron branch) are simply skipped: this crate only ever needs the real
//! part of the spectrum.

use crate::arnoldi::ArnoldiResult;
use nalgebra::{DMatrix, DVector, Schur};

#[derive(Debug, Clone)]
pub struct RitzPair {
    pub value: f64,
    /// Full-dimension vector (`basis * eigenvector_in_krylov_space`),
    /// *not* separately renormalized -- callers that care about scale
    /// (e.g. min-max normalizing a centrality score) should do that
    /// themselves.
    pub vector: Vec<f64>,
}

pub fn arnoldi_real_ritz_pairs(result: &ArnoldiResult) -> Vec<RitzPair> {
    let k = result.k;
    if k == 0 {
        return Vec::new();
    }
    let h = DMatrix::from_fn(k, k, |i, j| result.hessenberg[i][j]);

    if k == 1 {
        let lambda = h[(0, 0)];
        return vec![lift(result, lambda, &DVector::from_element(1, 1.0))];
    }

    let schur = Schur::new(h.clone());
    let complex_eigs = schur.complex_eigenvalues();

    let mut pairs = Vec::new();
    let mut seen: Vec<f64> = Vec::new();
    for c in complex_eigs.iter() {
        let scale = 1.0 + c.re.abs();
        if c.im.abs() > 1e-7 * scale {
            continue; // genuinely complex Ritz value: not part of this crate's contract
        }
        let lambda = c.re;
        // De-duplicate (a real eigenvalue appears once per algebraic
        // multiplicity in complex_eigenvalues(), and we only want one
        // eigenvector per distinct value here).
        if seen.iter().any(|&s| (s - lambda).abs() < 1e-9 * scale) {
            continue;
        }
        seen.push(lambda);
        if let Some(y) = inverse_iterate(&h, lambda) {
            pairs.push(lift(result, lambda, &y));
        }
    }
    pairs
}

fn lift(result: &ArnoldiResult, lambda: f64, y: &DVector<f64>) -> RitzPair {
    let dim = result.dim;
    let mut full = vec![0.0f64; dim];
    for (i, coeff) in y.iter().enumerate() {
        let basis_vec = &result.basis[i];
        for d in 0..dim {
            full[d] += coeff * basis_vec[d];
        }
    }
    RitzPair { value: lambda, vector: full }
}

/// A handful of steps of shifted inverse iteration on the small dense
/// matrix `h`, targeting eigenvalue `lambda`. Returns `None` only if the
/// shifted matrix is (numerically) exactly singular even after nudging the
/// shift, which shouldn't happen in practice for the shift nudges used
/// here.
fn inverse_iterate(h: &DMatrix<f64>, lambda: f64) -> Option<DVector<f64>> {
    let k = h.nrows();
    let scale = 1.0 + lambda.abs();
    for nudge_pow in 0..4 {
        let nudge = scale * 10f64.powi(-8 + nudge_pow);
        let shifted = h - DMatrix::identity(k, k) * (lambda + nudge);
        let Some(lu) = shifted.clone().lu().try_inverse() else { continue };
        let mut x = DVector::from_element(k, 1.0 / (k as f64).sqrt());
        for _ in 0..6 {
            let mut next = &lu * &x;
            let n = next.norm();
            if n < 1e-300 {
                break;
            }
            next /= n;
            x = next;
        }
        // Sanity check: residual of the *unshifted* problem should be small
        // relative to the vector's own scale.
        let resid = (h * &x) - &x * lambda;
        if resid.norm() < 1e-4 * scale {
            return Some(x);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arnoldi::Arnoldi;
    use crate::operator::DenseOperator;

    #[test]
    fn dominant_ritz_pair_matches_known_eigenvalue() {
        // [[2,1],[1,2]] has eigenvalues 1 and 3, eigenvector for 3 is (1,1)/sqrt(2).
        let op = DenseOperator { n: 2, rows: vec![vec![2.0, 1.0], vec![1.0, 2.0]] };
        let arnoldi = Arnoldi::new(2, 1e-12);
        let result = arnoldi.run(&op, &[1.0, 0.0]).unwrap();
        let pairs = arnoldi_real_ritz_pairs(&result);
        let dominant = pairs.iter().max_by(|a, b| a.value.abs().partial_cmp(&b.value.abs()).unwrap()).unwrap();
        assert!((dominant.value - 3.0).abs() < 1e-6, "value={}", dominant.value);
        let ratio = dominant.vector[0] / dominant.vector[1];
        assert!((ratio - 1.0).abs() < 1e-4, "ratio={ratio}");
    }
}
