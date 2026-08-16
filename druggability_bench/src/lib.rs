//! `druggability_bench`: validates `target_druggability`'s non-backtracking
//! residue centrality against classical centrality baselines
//! (degree, closeness, betweenness, plain eigenvector) on a labeled
//! cryptic-pocket dataset (e.g. CryptoSite), reporting ROC-AUC per method.
//!
//! See [`labels`] for the expected ground-truth CSV schema and where to
//! get real CryptoSite labels; see `data/example_labels.csv` for a small
//! synthetic smoke-test set that ships with this crate and validates
//! nothing about real-world performance.

pub mod baselines;
pub mod labels;
pub mod roc;

use labels::StructureLabels;
use roc::{roc_auc, AucResult, RocError};
use serde::Serialize;
use target_druggability::contact_graph::{self, ContactGraph, ContactParams};
use target_druggability::pdb::{self, PdbError, Structure};
use target_druggability::spectral_druggability::{self, SpectralError};

pub const METHODS: [&str; 5] =
    ["nonbacktracking", "degree", "closeness", "betweenness", "eigenvector"];

#[derive(Debug, thiserror::Error)]
pub enum BenchError {
    #[error(transparent)]
    Pdb(#[from] PdbError),
    #[error(transparent)]
    Spectral(#[from] SpectralError),
    #[error("no labels found for structure_id {0:?}")]
    NoLabels(String),
    #[error("structure {structure_id:?} has {labeled} labeled residues but 0 are positive, 0 are negative, or its contact graph has no residues at all -- cannot compute AUC")]
    Unusable { structure_id: String, labeled: usize },
}

#[derive(Debug, Serialize)]
pub struct MethodAuc {
    pub method: String,
    pub auc: Option<f64>,
    pub n_pos: usize,
    pub n_neg: usize,
}

#[derive(Debug, Serialize)]
pub struct StructureResult {
    pub structure_id: String,
    pub n_residues: usize,
    pub n_labeled: usize,
    pub methods: Vec<MethodAuc>,
}

fn score_all_methods(cg: &ContactGraph, seed: u64) -> Result<Vec<(&'static str, Vec<f64>)>, SpectralError> {
    let krylov_dim = (2 * cg.graph.n).min(200).max(4);
    let nb = spectral_druggability::residue_centrality(&cg.graph, krylov_dim, seed)?;
    Ok(vec![
        ("nonbacktracking", nb),
        ("degree", baselines::degree_centrality(&cg.graph)),
        ("closeness", baselines::closeness_centrality(&cg.graph)),
        ("betweenness", baselines::betweenness_centrality(&cg.graph)),
        ("eigenvector", baselines::eigenvector_centrality(&cg.graph, 200)),
    ])
}

fn auc_or_none(scores: &[f64], mask_labels: &[bool]) -> (Option<f64>, usize, usize) {
    match roc_auc(scores, mask_labels) {
        Ok(AucResult { auc, n_pos, n_neg }) => (Some(auc), n_pos, n_neg),
        Err(RocError::DegenerateLabels { n_pos, n_neg }) => (None, n_pos, n_neg),
        Err(RocError::LengthMismatch { .. }) => unreachable!("built from the same index set"),
    }
}

/// Per-method labeled-subset scores for one structure, for pooling across
/// structures without re-running the spectral solve.
pub type MethodSubsetScores = Vec<(&'static str, Vec<f64>)>;

/// Score one structure against its labels with every method, returning
/// per-method AUC (or `None` if that structure alone has only one class
/// present, e.g. an all-negative structure -- still reported with counts
/// so it's visible in output, just not misleadingly averaged as 0/1) plus
/// the raw labeled-subset scores/labels so the caller can pool across
/// structures without redoing the (relatively expensive) Arnoldi solve.
pub fn evaluate_structure(
    structure: &Structure,
    contact_params: ContactParams,
    structure_labels: &StructureLabels,
    seed: u64,
) -> Result<(StructureResult, MethodSubsetScores, Vec<bool>), BenchError> {
    let cg = contact_graph::build_contact_graph(structure, contact_params);

    // Join labels onto contact-graph node indices via (chain, res_seq);
    // residues not present in the label file are excluded, not treated
    // as negative (see labels module docs).
    let mut idx: Vec<usize> = Vec::new();
    let mut lbl: Vec<bool> = Vec::new();
    for (i, r) in cg.residues.iter().enumerate() {
        if let Some(&is_pos) = structure_labels.get(&(r.chain, r.res_seq)) {
            idx.push(i);
            lbl.push(is_pos);
        }
    }

    if idx.is_empty() {
        return Err(BenchError::Unusable { structure_id: structure.id.clone(), labeled: 0 });
    }

    let scored = score_all_methods(&cg, seed)?;
    let mut subset_scores: MethodSubsetScores = Vec::with_capacity(scored.len());
    let mut methods = Vec::with_capacity(scored.len());
    for (name, full_scores) in scored {
        let subset: Vec<f64> = idx.iter().map(|&i| full_scores[i]).collect();
        let (auc, n_pos, n_neg) = auc_or_none(&subset, &lbl);
        methods.push(MethodAuc { method: name.to_string(), auc, n_pos, n_neg });
        subset_scores.push((name, subset));
    }

    let result = StructureResult {
        structure_id: structure.id.clone(),
        n_residues: cg.graph.n,
        n_labeled: idx.len(),
        methods,
    };
    Ok((result, subset_scores, lbl))
}

#[derive(Debug, Serialize)]
pub struct PooledAuc {
    pub method: String,
    pub auc: Option<f64>,
    pub n_pos: usize,
    pub n_neg: usize,
}

#[derive(Debug, Serialize)]
pub struct BenchmarkReport {
    pub structures: Vec<StructureResult>,
    /// AUC pooling every labeled residue across every structure together
    /// (as opposed to averaging per-structure AUCs) -- the standard
    /// aggregate for benchmarks where individual structures may have few
    /// labeled residues.
    pub pooled: Vec<PooledAuc>,
}

pub fn build_report(structures: Vec<StructureResult>, per_method_pooled: Vec<(String, Vec<f64>, Vec<bool>)>) -> BenchmarkReport {
    let pooled = per_method_pooled
        .into_iter()
        .map(|(method, scores, labels)| {
            let (auc, n_pos, n_neg) = auc_or_none(&scores, &labels);
            PooledAuc { method, auc, n_pos, n_neg }
        })
        .collect();
    BenchmarkReport { structures, pooled }
}

pub fn load_structure(path: &std::path::Path) -> Result<Structure, PdbError> {
    pdb::parse_ca_atoms_file(path)
}
