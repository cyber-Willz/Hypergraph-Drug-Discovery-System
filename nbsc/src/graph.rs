//! Minimal undirected simple graph, adjacency-list representation.

/// An undirected, unweighted, simple graph on `n` labeled `0..n` nodes.
#[derive(Debug, Clone)]
pub struct Graph {
    pub n: usize,
    /// `neighbors[i]` is the (deduplicated) adjacency list of node `i`.
    pub neighbors: Vec<Vec<usize>>,
}

impl Graph {
    pub fn new(n: usize) -> Self {
        Self { n, neighbors: vec![Vec::new(); n] }
    }

    /// Add an undirected edge `i -- j`. Self-loops and duplicate edges are
    /// silently ignored (idempotent).
    pub fn add_edge(&mut self, i: usize, j: usize) {
        assert!(i < self.n && j < self.n, "node index out of bounds");
        if i == j {
            return;
        }
        if !self.neighbors[i].contains(&j) {
            self.neighbors[i].push(j);
            self.neighbors[j].push(i);
        }
    }

    pub fn degree(&self, i: usize) -> usize {
        self.neighbors[i].len()
    }

    /// Number of (undirected) edges.
    pub fn m(&self) -> usize {
        self.neighbors.iter().map(|v| v.len()).sum::<usize>() / 2
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add_edge_is_undirected_and_deduplicated() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1);
        g.add_edge(1, 0); // duplicate, opposite order
        g.add_edge(1, 1); // self loop, ignored
        assert_eq!(g.neighbors[0], vec![1]);
        assert_eq!(g.neighbors[1], vec![0]);
        assert_eq!(g.m(), 1);
    }
}
