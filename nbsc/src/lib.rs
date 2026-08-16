//! `nbsc`: Non-Backtracking Spectral Convolution primitives -- the
//! `Graph` type and Hashimoto/non-backtracking-operator spectral machinery
//! (`rho_B` estimation, and the matrix-free Bass/Ihara linearization other
//! crates in this workspace build per-vertex centrality on top of).
//!
//! This is the minimal `default-features = false` surface: just
//! [`graph::Graph`] and [`spectral`]. Feature-gated extras that exist
//! elsewhere in this workspace (compliance channel scoring, hypergraph
//! bridges, incremental resync, etc.) are out of scope here.

pub mod graph;
pub mod spectral;
