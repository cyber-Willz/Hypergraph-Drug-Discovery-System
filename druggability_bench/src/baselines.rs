//! Classical centrality baselines for the residue contact graph, computed
//! independently of `nbsc`/`krylov_ds` so the comparison isn't leaning on
//! the same machinery being validated. All four are the standard measures
//! cited in protein-structure-network literature as candidate pocket/
//! allosteric-site predictors, which is exactly what non-backtracking
//! centrality is being benchmarked against.
//!
//! All scores are min-max normalized to `[0, 1]` so they're on the same
//! footing as `target_druggability::spectral_druggability::residue_centrality`
//! for the AUC comparison (AUC is rank-invariant under monotone rescaling
//! anyway, but normalizing keeps any downstream reporting/plots sane).

use nbsc::graph::Graph;
use std::collections::VecDeque;

fn normalize(mut v: Vec<f64>) -> Vec<f64> {
    let max = v.iter().cloned().fold(0.0, f64::max);
    if max > 0.0 {
        for x in v.iter_mut() {
            *x /= max;
        }
    }
    v
}

/// Plain degree centrality: fraction of the graph a node is connected to.
/// This is exactly the "locally dense / buried residue" signal that
/// non-backtracking centrality is designed to suppress -- the baseline
/// most likely to lose this comparison, included because it's the
/// simplest possible one and a useful floor.
pub fn degree_centrality(graph: &Graph) -> Vec<f64> {
    let raw: Vec<f64> = (0..graph.n).map(|i| graph.degree(i) as f64).collect();
    normalize(raw)
}

/// Closeness centrality (Wasserman-Faust variant, robust to
/// disconnected graphs): for node `i`, `(reachable - 1)^2 / ((n-1) *
/// sum_of_distances_to_reachable_nodes)`. Unweighted BFS shortest paths.
pub fn closeness_centrality(graph: &Graph) -> Vec<f64> {
    let n = graph.n;
    let mut raw = vec![0.0f64; n];
    for s in 0..n {
        let mut dist = vec![-1i64; n];
        dist[s] = 0;
        let mut q = VecDeque::new();
        q.push_back(s);
        let mut sum_dist = 0i64;
        let mut reachable = 0i64;
        while let Some(u) = q.pop_front() {
            for &v in &graph.neighbors[u] {
                if dist[v] == -1 {
                    dist[v] = dist[u] + 1;
                    sum_dist += dist[v];
                    reachable += 1;
                    q.push_back(v);
                }
            }
        }
        raw[s] = if sum_dist > 0 && n > 1 {
            (reachable as f64).powi(2) / ((n as f64 - 1.0) * sum_dist as f64)
        } else {
            0.0
        };
    }
    normalize(raw)
}

/// Betweenness centrality via Brandes' algorithm (unweighted, O(n*m)),
/// the standard measure network-science papers mean by "betweenness" and
/// the closest classical analogue to what non-backtracking centrality is
/// claimed to capture (bridging residues between otherwise-distant
/// regions).
pub fn betweenness_centrality(graph: &Graph) -> Vec<f64> {
    let n = graph.n;
    let mut centrality = vec![0.0f64; n];

    for s in 0..n {
        let mut stack: Vec<usize> = Vec::new();
        let mut preds: Vec<Vec<usize>> = vec![Vec::new(); n];
        let mut sigma = vec![0.0f64; n];
        let mut dist = vec![-1i64; n];
        sigma[s] = 1.0;
        dist[s] = 0;
        let mut queue = VecDeque::new();
        queue.push_back(s);

        while let Some(v) = queue.pop_front() {
            stack.push(v);
            for &w in &graph.neighbors[v] {
                if dist[w] < 0 {
                    dist[w] = dist[v] + 1;
                    queue.push_back(w);
                }
                if dist[w] == dist[v] + 1 {
                    sigma[w] += sigma[v];
                    preds[w].push(v);
                }
            }
        }

        let mut delta = vec![0.0f64; n];
        while let Some(w) = stack.pop() {
            for &v in &preds[w] {
                delta[v] += (sigma[v] / sigma[w]) * (1.0 + delta[w]);
            }
            if w != s {
                centrality[w] += delta[w];
            }
        }
    }

    // Undirected graph: every shortest path is counted from both
    // endpoints' BFS, so divide by 2.
    for c in centrality.iter_mut() {
        *c /= 2.0;
    }
    normalize(centrality)
}

/// Plain (backtracking-allowed) eigenvector centrality via power
/// iteration on the adjacency matrix. This is the direct "what if we
/// didn't suppress backtracking walks" ablation baseline -- the most
/// informative comparison for whether the non-backtracking construction
/// specifically is doing the work, versus eigenvector centrality generally.
pub fn eigenvector_centrality(graph: &Graph, iterations: usize) -> Vec<f64> {
    let n = graph.n;
    if n == 0 {
        return Vec::new();
    }
    let mut x = vec![1.0f64 / (n as f64).sqrt(); n];
    for _ in 0..iterations {
        let mut next = vec![0.0f64; n];
        for u in 0..n {
            for &v in &graph.neighbors[u] {
                next[u] += x[v];
            }
        }
        let norm = next.iter().map(|v| v * v).sum::<f64>().sqrt();
        if norm > 1e-15 {
            for v in next.iter_mut() {
                *v /= norm;
            }
        }
        x = next;
    }
    normalize(x.into_iter().map(f64::abs).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn path_graph(n: usize) -> Graph {
        let mut g = Graph::new(n);
        for i in 0..n - 1 {
            g.add_edge(i, i + 1);
        }
        g
    }

    #[test]
    fn betweenness_peaks_at_path_center() {
        let g = path_graph(5); // 0-1-2-3-4
        let b = betweenness_centrality(&g);
        assert!(b[2] > b[0] && b[2] > b[4], "center of a path should be the betweenness peak");
    }

    #[test]
    fn degree_uniform_on_cycle() {
        let mut g = Graph::new(4);
        g.add_edge(0, 1);
        g.add_edge(1, 2);
        g.add_edge(2, 3);
        g.add_edge(3, 0);
        let d = degree_centrality(&g);
        for &x in &d {
            assert!((x - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn eigenvector_favors_hub() {
        // star graph: node 0 connected to 1..5
        let mut g = Graph::new(6);
        for i in 1..6 {
            g.add_edge(0, i);
        }
        let ev = eigenvector_centrality(&g, 100);
        assert!(ev[0] > ev[1]);
    }

    #[test]
    fn closeness_favors_center_of_star() {
        let mut g = Graph::new(6);
        for i in 1..6 {
            g.add_edge(0, i);
        }
        let c = closeness_centrality(&g);
        assert!(c[0] > c[1]);
    }
}
