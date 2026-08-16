//! Combine per-target spectral pocket candidates with NIH IDG druggability
//! tier metadata into one ranked report: "difficult (Tdark/Tbio) targets
//! whose structure shows strong candidate cryptic-pocket/allosteric
//! signal, ranked by how strong that signal is."
//!
//! This is a prioritization aid, not a discovery pipeline end-to-end. It
//! tells you *where to look next* (which targets, which residues), not
//! "this is a druggable pocket" -- that still needs structural biology
//! judgment, MD simulation / cryptic-pocket-opening simulation, and
//! ultimately experimental validation (fragment screening, SPR, etc).

use crate::contact_graph::ContactGraph;
use crate::spectral_druggability::PocketCandidate;
use crate::tcrd::TcrdTarget;
use serde::Serialize;

#[derive(Debug, Serialize)]
pub struct ResidueReport {
    pub chain: char,
    pub res_seq: i32,
    pub res_name: String,
    pub score: f64,
}

#[derive(Debug, Serialize)]
pub struct PocketReport {
    pub rank: usize,
    pub mean_score: f64,
    pub max_score: f64,
    pub size: usize,
    pub residues: Vec<ResidueReport>,
}

#[derive(Debug, Serialize)]
pub struct TargetReport {
    pub structure_id: String,
    pub n_residues: usize,
    pub n_contacts: usize,
    pub global_coupling_rho_b: f64,
    /// None if no TCRD/Pharos metadata matched this structure (no symbol
    /// or UniProt ID supplied, or not found in the loaded export).
    pub tcrd: Option<TcrdSummary>,
    pub pockets: Vec<PocketReport>,
}

#[derive(Debug, Serialize)]
pub struct TcrdSummary {
    pub symbol: String,
    pub uniprot: String,
    pub tdl: String,
    pub family: String,
    pub is_difficult_target: bool,
}

pub fn build_report(
    structure_id: &str,
    contact_graph: &ContactGraph,
    residue_scores: &[f64],
    global_coupling_rho_b: f64,
    pockets: &[PocketCandidate],
    tcrd: Option<&TcrdTarget>,
) -> TargetReport {
    let pocket_reports = pockets
        .iter()
        .enumerate()
        .map(|(i, p)| PocketReport {
            rank: i + 1,
            mean_score: p.mean_score,
            max_score: p.max_score,
            size: p.residue_indices.len(),
            residues: p
                .residue_indices
                .iter()
                .map(|&idx| {
                    let r = &contact_graph.residues[idx];
                    ResidueReport {
                        chain: r.chain,
                        res_seq: r.res_seq,
                        res_name: r.res_name.clone(),
                        score: residue_scores[idx],
                    }
                })
                .collect(),
        })
        .collect();

    TargetReport {
        structure_id: structure_id.to_string(),
        n_residues: contact_graph.graph.n,
        n_contacts: contact_graph.graph.m(),
        global_coupling_rho_b,
        tcrd: tcrd.map(|t| TcrdSummary {
            symbol: t.symbol.clone(),
            uniprot: t.uniprot.clone(),
            tdl: t.tdl.clone(),
            family: t.family.clone(),
            is_difficult_target: t.is_difficult(),
        }),
        pockets: pocket_reports,
    }
}

/// Rank a batch of target reports for triage: difficult (Tdark/Tbio)
/// targets with strong pocket signal first, then everything else by
/// top-pocket score.
pub fn rank_targets(mut reports: Vec<TargetReport>) -> Vec<TargetReport> {
    reports.sort_by(|a, b| {
        let a_diff = a.tcrd.as_ref().map(|t| t.is_difficult_target).unwrap_or(false);
        let b_diff = b.tcrd.as_ref().map(|t| t.is_difficult_target).unwrap_or(false);
        let a_top = a.pockets.first().map(|p| p.mean_score).unwrap_or(0.0);
        let b_top = b.pockets.first().map(|p| p.mean_score).unwrap_or(0.0);
        b_diff.cmp(&a_diff).then(b_top.partial_cmp(&a_top).unwrap())
    });
    reports
}
