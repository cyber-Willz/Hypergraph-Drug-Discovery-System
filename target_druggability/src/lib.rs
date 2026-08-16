//! `target_druggability`: non-backtracking (Hashimoto) spectral analysis of
//! protein residue contact networks, applied to prioritizing NIH
//! Illuminating the Druggable Genome (IDG) Tdark/Tbio ("difficult to drug")
//! targets for follow-up.
//!
//! ## Pipeline
//! 1. [`pdb`] -- parse C-alpha coordinates from a PDB structure file.
//! 2. [`contact_graph`] -- build a residue contact network (an
//!    [`nbsc::graph::Graph`]) from those coordinates.
//! 3. [`spectral_druggability`] -- run the non-backtracking / Hashimoto
//!    spectral machinery already implemented in `nbsc` (`rho_B`) plus a
//!    new per-residue non-backtracking centrality, and cluster
//!    high-scoring residues into candidate pocket/hotspot regions.
//! 4. [`tcrd`] -- load NIH IDG target development level (TDL) tier
//!    metadata (Tclin/Tchem/Tbio/Tdark) to know *which* targets are
//!    actually understudied/difficult.
//! 5. [`report`] -- combine the two into a ranked, serializable report.
//!
//! ## Honest scope statement
//! This is a structural-network heuristic for prioritizing where to look,
//! built on real spectral graph theory machinery, applied to real PDB
//! coordinates. It is **not** a validated druggability predictor: it has
//! not been benchmarked against known cryptic-pocket datasets (e.g.
//! CryptoSite) here, and non-backtracking centrality's correlation with
//! true allosteric/pocket sites is a reasonable structural-biology
//! hypothesis (borne out in the broader protein-network literature for
//! *related* centrality measures) rather than something this codebase
//! itself has validated. Treat output as a ranked hypothesis list to
//! triage against domain expertise and wet-lab/simulation follow-up, not
//! as a final answer. This crate also has no network access: fetching PDB
//! structures and TCRD/Pharos exports is left to the user (see module docs
//! in `pdb` and `tcrd` for exact sources/URLs).

pub mod contact_graph;
pub mod pdb;
pub mod report;
pub mod spectral_druggability;
pub mod tcrd;
