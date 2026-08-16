//! Non-backtracking (Hashimoto) spectral scoring of a residue contact
//! network, used as a structural proxy for "this region of the protein is
//! densely, non-locally coupled" -- the same structural signature that
//! shows up at allosteric sites and cryptic pockets in the protein-network
//! literature (network hubs / high betweenness residues correlate with
//! allosteric communication paths; e.g. Vishveshwara-style protein
//! structure network analysis, and elastic-network/normal-mode studies of
//! cryptic pocket opening).
//!
//! This module makes exactly one modeling claim and is explicit about it:
//! high non-backtracking eigenvector centrality identifies *candidate*
//! long-range structurally-coupled residues, worth prioritizing for
//! follow-up (MD simulation, fragment screening, experimental validation)
//! -- it is not a druggability predictor on its own and is not validated
//! against experimental pocket data here. Treat its output as a ranked
//! hypothesis list, not a verdict.
//!
//! ## Why non-backtracking instead of plain adjacency/degree centrality
//! Ordinary degree or eigenvector centrality on a contact graph is
//! dominated by locally dense, sterically packed regions (the protein
//! core), which is mostly just "this residue is buried" -- not
//! informative for pocket-finding. The non-backtracking (Hashimoto)
//! operator suppresses that local-density signal by construction (walks
//! that immediately reverse themselves don't count), which empirically
//! shifts weight toward residues that bridge otherwise-distant parts of
//! the structure -- the graph-theoretic signature of an allosteric
//! communication path. This is the same rho_B / non-backtracking-spectrum
//! machinery already used elsewhere in this workspace (`nbsc`), applied
//! here to a physical contact network instead of a citation/social graph.

use krylov_ds::eig::arnoldi_real_ritz_pairs;
use krylov_ds::operator::LinearOperator;
use krylov_ds::Arnoldi;
use nbsc::graph::Graph;
use nbsc::spectral::{estimate_spectral_radius, HashimotoLinearization};

/// Global structural-coupling score for the whole target: `rho_B`, the
/// Perron-Frobenius spectral radius of the Hashimoto matrix. Higher values
/// indicate a more expander-like, densely non-locally-connected fold;
/// useful for comparing targets against each other, not for localizing
/// pockets within one target (see [`residue_centrality`] for that).
pub fn global_coupling_score(graph: &Graph, seed: u64) -> f64 {
    let krylov_dim = (2 * graph.n).min(40).max(2);
    estimate_spectral_radius(graph, krylov_dim, seed)
}

#[derive(Debug, thiserror::Error)]
pub enum SpectralError {
    #[error("graph has fewer than 2 nodes; cannot run non-backtracking spectral analysis")]
    TooSmall,
    #[error("Arnoldi failed to converge on the Hashimoto linearization: {0}")]
    ArnoldiFailed(String),
    #[error("no real dominant eigenpair found (Hashimoto spectrum was entirely complex at this Krylov dimension -- try increasing krylov_dim)")]
    NoRealDominantEigenpair,
}

/// Per-residue non-backtracking eigenvector centrality.
///
/// Computed via the same Bass-reduced `2n x 2n` linearization `M` that
/// `nbsc::spectral::estimate_spectral_radius` uses for `rho_B`, but here we
/// keep the eigenvector (not just the eigenvalue): the vertex-space block
/// (first `n` of the `2n` components) of the dominant real Ritz pair is,
/// up to normalization, the vertex-projected non-backtracking Perron
/// eigenvector -- the standard non-backtracking centrality of Martin,
/// Zhang & Newman (2014), computed here without ever forming the `2m x 2m`
/// directed-edge Hashimoto matrix explicitly.
///
/// Returns one score per node, index-aligned with `graph`'s node indices
/// (and therefore with `ContactGraph::residues`), min-max normalized to
/// `[0, 1]` for readability.
pub fn residue_centrality(
    graph: &Graph,
    krylov_dim: usize,
    seed: u64,
) -> Result<Vec<f64>, SpectralError> {
    if graph.n < 2 {
        return Err(SpectralError::TooSmall);
    }
    let op = HashimotoLinearization::new(graph);
    let dim = op.dim(); // 2n
    let m = krylov_dim.min(dim).max(2);

    let mut state = seed.wrapping_mul(2685821657736338717).wrapping_add(1);
    let mut v0 = vec![0.0f64; dim];
    for vi in v0.iter_mut() {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        *vi = (state as f64 / u64::MAX as f64) - 0.5;
    }

    let arnoldi = Arnoldi::new(m, 1e-12);
    let result =
        arnoldi.run(&op, &v0).map_err(|e| SpectralError::ArnoldiFailed(format!("{e:?}")))?;

    let pairs = arnoldi_real_ritz_pairs(&result);
    let dominant = pairs
        .into_iter()
        .max_by(|a, b| a.value.abs().partial_cmp(&b.value.abs()).unwrap())
        .ok_or(SpectralError::NoRealDominantEigenpair)?;

    let n = graph.n;
    let mut scores: Vec<f64> = dominant.vector[..n].iter().map(|x| x.abs()).collect();

    let max = scores.iter().cloned().fold(0.0, f64::max);
    if max > 0.0 {
        for s in scores.iter_mut() {
            *s /= max;
        }
    }
    Ok(scores)
}

/// A spatially contiguous cluster of high-scoring residues: a candidate
/// cryptic pocket / allosteric hotspot.
#[derive(Debug, Clone, serde::Serialize)]
pub struct PocketCandidate {
    /// Node indices (into the residue list) making up this cluster.
    pub residue_indices: Vec<usize>,
    pub mean_score: f64,
    pub max_score: f64,
}

/// Cluster the top-`percentile` scoring residues into spatially contiguous
/// groups using connected components of the contact graph restricted to
/// that residue set. `percentile` is e.g. 0.9 for "top 10% of residues by
/// score". Clusters of size 1 are dropped (an isolated high-scoring
/// residue with no high-scoring neighbors is more likely noise than a
/// pocket).
pub fn cluster_pockets(
    graph: &Graph,
    scores: &[f64],
    percentile: f64,
    min_cluster_size: usize,
) -> Vec<PocketCandidate> {
    assert_eq!(graph.n, scores.len());
    let mut sorted = scores.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let idx = ((sorted.len() as f64 - 1.0) * percentile).round().max(0.0) as usize;
    let threshold = sorted.get(idx).copied().unwrap_or(1.0);

    let hot: Vec<bool> = scores.iter().map(|&s| s >= threshold).collect();
    let mut visited = vec![false; graph.n];
    let mut clusters = Vec::new();

    for start in 0..graph.n {
        if !hot[start] || visited[start] {
            continue;
        }
        let mut stack = vec![start];
        let mut component = Vec::new();
        visited[start] = true;
        while let Some(v) = stack.pop() {
            component.push(v);
            for &u in &graph.neighbors[v] {
                if hot[u] && !visited[u] {
                    visited[u] = true;
                    stack.push(u);
                }
            }
        }
        if component.len() >= min_cluster_size {
            let mean_score = component.iter().map(|&i| scores[i]).sum::<f64>() / component.len() as f64;
            let max_score = component.iter().map(|&i| scores[i]).fold(0.0, f64::max);
            clusters.push(PocketCandidate { residue_indices: component, mean_score, max_score });
        }
    }

    clusters.sort_by(|a, b| b.mean_score.partial_cmp(&a.mean_score).unwrap());
    clusters
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two dense cliques bridged by a single path -- the bridge residues
    /// should score highest under non-backtracking centrality, since a
    /// clique's internal high degree is exactly the "locally dense, not
    /// informative" signal this metric is designed to de-emphasize.
    #[test]
    fn bridge_residues_score_higher_than_clique_interior() {
        let mut g = Graph::new(10);
        // clique A: 0..5
        for i in 0..5 {
            for j in (i + 1)..5 {
                g.add_edge(i, j);
            }
        }
        // clique B: 5..10
        for i in 5..10 {
            for j in (i + 1)..10 {
                g.add_edge(i, j);
            }
        }
        // bridge: connect the cliques with one extra edge between 2 and 7
        g.add_edge(2, 7);

        let scores = residue_centrality(&g, 40, 42).unwrap();
        let bridge_score = (scores[2] + scores[7]) / 2.0;
        let interior_score =
            [0usize, 1, 3, 4].iter().map(|&i| scores[i]).sum::<f64>() / 4.0;
        assert!(
            bridge_score > interior_score,
            "expected bridge residues to outscore clique interior: bridge={bridge_score}, interior={interior_score}"
        );
    }

    #[test]
    fn cluster_pockets_groups_contiguous_hotspots() {
        let mut g = Graph::new(5);
        g.add_edge(0, 1);
        g.add_edge(1, 2);
        // node 4 isolated hotspot, node 3 low score
        let scores = vec![0.9, 0.95, 0.85, 0.1, 0.9];
        let clusters = cluster_pockets(&g, &scores, 0.3, 1);
        // {0,1,2} should form one contiguous cluster; {4} isolated singleton
        assert!(clusters.iter().any(|c| c.residue_indices.len() == 3));
    }
}
