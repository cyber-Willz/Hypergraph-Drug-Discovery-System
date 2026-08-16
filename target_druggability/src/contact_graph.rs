//! Build a residue-residue contact graph (the standard "residue interaction
//! network", RIN) from C-alpha coordinates, in the representation
//! [`nbsc::graph::Graph`] expects.
//!
//! Two residues are connected if their C-alpha atoms are within
//! `distance_cutoff` Angstroms *and* they are separated by at least
//! `min_seq_separation` positions along the chain. The sequence-separation
//! filter is standard practice in protein contact-network analysis: without
//! it, the graph is dominated by trivial local backbone contacts (i, i+1),
//! (i, i+2) that carry no tertiary-structure information and would drown
//! out the long-range contacts that actually encode allosteric coupling.
//!
//! Default cutoff (8.0 A) is the common choice for Cα-Cα contact maps in
//! the protein structure network literature (e.g. Vishveshwara et al.'s
//! "protein structure networks"; elastic network model cutoffs are usually
//! quoted in the 7-12 A range depending on whether Cα or Cβ is used).

use crate::pdb::Structure;
use nbsc::graph::Graph;

#[derive(Debug, Clone, Copy)]
pub struct ContactParams {
    pub distance_cutoff: f64,
    pub min_seq_separation: i32,
}

impl Default for ContactParams {
    fn default() -> Self {
        Self { distance_cutoff: 8.0, min_seq_separation: 3 }
    }
}

/// A contact graph plus the residue metadata needed to map node indices
/// back to (chain, resSeq, resName) for reporting.
pub struct ContactGraph {
    pub graph: Graph,
    /// `residues[i]` is the residue at graph node `i`.
    pub residues: Vec<crate::pdb::Residue>,
}

fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

pub fn build_contact_graph(structure: &Structure, params: ContactParams) -> ContactGraph {
    let n = structure.residues.len();
    let mut graph = Graph::new(n);

    for i in 0..n {
        for j in (i + 1)..n {
            let ri = &structure.residues[i];
            let rj = &structure.residues[j];
            let same_chain = ri.chain == rj.chain;
            let seq_sep = (ri.res_seq - rj.res_seq).abs();
            if same_chain && seq_sep < params.min_seq_separation {
                continue;
            }
            if dist(ri.ca, rj.ca) <= params.distance_cutoff {
                graph.add_edge(i, j);
            }
        }
    }

    ContactGraph { graph, residues: structure.residues.clone() }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pdb::Residue;

    fn res(chain: char, seq: i32, ca: [f64; 3]) -> Residue {
        Residue { chain, res_seq: seq, res_name: "ALA".into(), ca }
    }

    #[test]
    fn connects_nearby_nonlocal_residues() {
        let structure = Structure {
            id: "t".into(),
            residues: vec![
                res('A', 1, [0.0, 0.0, 0.0]),
                res('A', 2, [3.8, 0.0, 0.0]),
                res('A', 3, [7.6, 0.0, 0.0]),
                res('A', 4, [4.0, 4.0, 0.0]), // ~ within 8A of residue 1 (dist 5.66)
            ],
        };
        let cg = build_contact_graph(&structure, ContactParams::default());
        // residue 0 and 3: seq_sep = 3 (>= min_seq_separation), dist ~5.66 <= 8.0
        assert!(cg.graph.neighbors[0].contains(&3));
        // residue 0 and 1: seq_sep = 1 (< 3), must be excluded even though close in space
        assert!(!cg.graph.neighbors[0].contains(&1));
    }

    #[test]
    fn different_chains_ignore_seq_separation() {
        let structure = Structure {
            id: "t".into(),
            residues: vec![res('A', 1, [0.0, 0.0, 0.0]), res('B', 1, [3.0, 0.0, 0.0])],
        };
        let cg = build_contact_graph(&structure, ContactParams::default());
        assert!(cg.graph.neighbors[0].contains(&1));
    }
}
