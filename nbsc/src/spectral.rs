//! Non-backtracking (Hashimoto) spectral machinery.
//!
//! Rather than forming the `2m x 2m` directed-edge non-backtracking matrix
//! `B` explicitly, this module works with the standard Bass/Ihara `2n x 2n`
//! linearization `M`, which has the same nonzero spectrum as `B` (up to the
//! trivial `+-1` eigenvalues Ihara's theorem accounts for separately) and
//! is far cheaper for sparse graphs (`n` vertices vs. `2m` directed edges).
//!
//! Ihara's theorem: `det(I - u B) = (1 - u^2)^(m-n) det(I - u A + u^2 (D - I))`.
//! Substituting `u = 1/lambda` and clearing denominators, a nonzero,
//! non-trivial eigenvalue `lambda` of `B` satisfies, for some vertex
//! vector `v`:
//!
//! ```text
//! lambda^2 v - lambda A v + (D - I) v = 0
//! ```
//!
//! Linearizing with `w = lambda v` gives the `2n x 2n` eigenproblem
//! `M (w, v)^T = lambda (w, v)^T` for
//!
//! ```text
//! M = [ A      (I - D) ]
//!     [ I         0     ]
//! ```
//!
//! i.e. `(Mz)_w = A w + (I - D) v` and `(Mz)_v = w`, applied here
//! matrix-free in [`HashimotoLinearization::apply`].
//!
//! For a connected graph, `B` (and hence, on its dominant branch, `M`) is a
//! non-negative matrix in the Perron-Frobenius sense, so the
//! spectral-radius eigenvalue `rho_B` is real and positive, and its
//! eigenvector is (up to sign) non-negative -- which is what makes reading
//! off `|w_i|` as a per-vertex non-backtracking centrality meaningful.

use crate::graph::Graph;
use krylov_ds::eig::arnoldi_real_ritz_pairs;
use krylov_ds::operator::LinearOperator;
use krylov_ds::Arnoldi;

/// Matrix-free `2n x 2n` Bass/Ihara linearization of a graph's
/// non-backtracking (Hashimoto) operator. See module docs for the exact
/// linear map.
pub struct HashimotoLinearization<'a> {
    graph: &'a Graph,
}

impl<'a> HashimotoLinearization<'a> {
    pub fn new(graph: &'a Graph) -> Self {
        Self { graph }
    }
}

impl<'a> LinearOperator for HashimotoLinearization<'a> {
    fn dim(&self) -> usize {
        2 * self.graph.n
    }

    fn apply(&self, z: &[f64]) -> Vec<f64> {
        let n = self.graph.n;
        let (w, v) = z.split_at(n);
        let mut new_w = vec![0.0; n];
        for i in 0..n {
            let deg = self.graph.degree(i) as f64;
            let mut s = (1.0 - deg) * v[i];
            for &j in &self.graph.neighbors[i] {
                s += w[j];
            }
            new_w[i] = s;
        }
        let mut out = new_w;
        out.extend_from_slice(w); // new_v = w
        out
    }
}

fn norm(a: &[f64]) -> f64 {
    a.iter().map(|x| x * x).sum::<f64>().sqrt()
}

fn seeded_unit_vector(dim: usize, seed: u64) -> Vec<f64> {
    // xorshift64*, matching the seeding style used elsewhere in this
    // workspace for reproducible (but not cryptographic) start vectors.
    let mut state = seed.wrapping_mul(2685821657736338717).wrapping_add(1);
    let mut v = vec![0.0f64; dim];
    for vi in v.iter_mut() {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        *vi = (state as f64 / u64::MAX as f64) - 0.5;
    }
    let n = norm(&v);
    if n > 0.0 {
        for vi in v.iter_mut() {
            *vi /= n;
        }
    } else {
        v[0] = 1.0;
    }
    v
}

/// Estimate `rho_B`, the Perron-Frobenius spectral radius of the
/// non-backtracking (Hashimoto) matrix, as the largest-magnitude real Ritz
/// value of an Arnoldi projection of the matrix-free `2n x 2n`
/// linearization onto a `krylov_dim`-dimensional Krylov subspace.
///
/// Using the same Arnoldi + real-Ritz-pair machinery as
/// [`crate::graph`]-consuming callers that need eigenvectors (not just the
/// eigenvalue) matters here for more than code reuse: plain power
/// iteration's convergence rate is governed by the ratio of the top two
/// eigenvalue magnitudes, which degenerates badly (or fails outright to
/// pick a direction) when those are tied or nearly so -- as happens for
/// e.g. any regular cycle, where the non-backtracking spectrum sits
/// entirely on the unit circle. Reading `rho_B` off the Schur form of the
/// small Hessenberg projection instead handles repeated/near-repeated
/// dominant eigenvalues correctly.
///
/// Graphs with fewer than 2 nodes have no meaningful non-backtracking
/// structure; returns `0.0` for those degenerate cases rather than
/// panicking, since a whole-target "how expander-like is this fold" score
/// of zero is a reasonable answer for "there's no fold here."
pub fn estimate_spectral_radius(graph: &Graph, krylov_dim: usize, seed: u64) -> f64 {
    if graph.n < 2 {
        return 0.0;
    }
    let op = HashimotoLinearization::new(graph);
    let dim = op.dim();
    let v0 = seeded_unit_vector(dim, seed);

    let arnoldi = Arnoldi::new(krylov_dim.max(2), 1e-12);
    let Ok(result) = arnoldi.run(&op, &v0) else { return 0.0 };
    let pairs = arnoldi_real_ritz_pairs(&result);
    pairs.iter().map(|p| p.value.abs()).fold(0.0, f64::max)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cycle_graph_has_known_rho_b() {
        // For a k-regular graph, rho_B = degree - 1 (here: a 6-cycle, all
        // degree 2, so rho_B = 1).
        let mut g = Graph::new(6);
        for i in 0..6 {
            g.add_edge(i, (i + 1) % 6);
        }
        let rho = estimate_spectral_radius(&g, 60, 7);
        assert!((rho - 1.0).abs() < 1e-3, "rho_B={rho}");
    }

    #[test]
    fn complete_graph_has_known_rho_b() {
        // For K_n (n-1 regular), rho_B = n - 2.
        let mut g = Graph::new(5);
        for i in 0..5 {
            for j in (i + 1)..5 {
                g.add_edge(i, j);
            }
        }
        let rho = estimate_spectral_radius(&g, 60, 7);
        assert!((rho - 3.0).abs() < 1e-2, "rho_B={rho}");
    }
}
