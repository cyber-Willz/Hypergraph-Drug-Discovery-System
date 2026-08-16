//! Loader for NIH Illuminating the Druggable Genome (IDG) target metadata
//! -- specifically the Target Development Level (TDL) tier, which is the
//! actual "how difficult/understudied is this target" signal from the NIH
//! Common Fund's IDG program:
//!
//! - **Tclin**: target of an approved drug with a known mechanism of action.
//! - **Tchem**: has a potent, selective small-molecule/biologic modulator
//!   but no approved drug yet.
//! - **Tbio**: biologically characterized (e.g. has a GO annotation, is
//!   studied in model organisms) but has no known potent chemical modulator.
//! - **Tdark**: minimally characterized -- these are the actual "dark
//!   genome" targets the IDG program exists to illuminate.
//!
//! This crate has no network access, so it cannot query TCRD/Pharos
//! directly. Get the real data yourself:
//! - Pharos UI + CSV export: <https://pharos.nih.gov/targets>
//! - Pharos GraphQL API: <https://pharos-api.ncats.io/graphql>
//! - Full TCRD relational dump: <http://juniper.health.unm.edu/tcrd/download/>
//!
//! and export/convert to the flat CSV schema this module reads (see
//! `data/sample_tcrd_export.csv` in this crate for the exact header and a
//! few real, publicly documented example rows).

use serde::Deserialize;
use std::collections::HashMap;
use std::path::Path;

#[derive(Debug, Clone, Deserialize)]
pub struct TcrdTarget {
    pub symbol: String,
    pub uniprot: String,
    pub name: String,
    /// One of Tclin / Tchem / Tbio / Tdark (case-insensitive on load).
    pub tdl: String,
    /// IDG protein family, e.g. GPCR, Kinase, IonChannel, NR, Other.
    pub family: String,
    /// Optional: IDG "novelty score" or similar continuous under-study
    /// metric, if present in your export. NaN/absent if not provided.
    #[serde(default)]
    pub novelty_score: Option<f64>,
}

impl TcrdTarget {
    pub fn is_difficult(&self) -> bool {
        matches!(self.tdl.to_ascii_uppercase().as_str(), "TDARK" | "TBIO")
    }
}

#[derive(Debug, Default)]
pub struct TcrdIndex {
    by_symbol: HashMap<String, TcrdTarget>,
    by_uniprot: HashMap<String, TcrdTarget>,
}

impl TcrdIndex {
    pub fn lookup_symbol(&self, symbol: &str) -> Option<&TcrdTarget> {
        self.by_symbol.get(&symbol.to_ascii_uppercase())
    }

    pub fn lookup_uniprot(&self, uniprot: &str) -> Option<&TcrdTarget> {
        self.by_uniprot.get(&uniprot.to_ascii_uppercase())
    }
}

#[derive(Debug, thiserror::Error)]
pub enum TcrdError {
    #[error("I/O error reading TCRD export: {0}")]
    Io(#[from] std::io::Error),
    #[error("CSV parse error at line {line}: {msg}")]
    Parse { line: usize, msg: String },
}

/// Parse a flat CSV with header `symbol,uniprot,name,tdl,family,novelty_score`
/// (novelty_score column optional / may be empty per-row).
pub fn load_csv(path: &Path) -> Result<TcrdIndex, TcrdError> {
    let text = std::fs::read_to_string(path)?;
    let mut lines = text.lines();
    let header = lines.next().unwrap_or_default();
    let cols: Vec<&str> = header.split(',').map(|s| s.trim()).collect();
    let col_idx = |name: &str| cols.iter().position(|c| c.eq_ignore_ascii_case(name));
    let (Some(i_sym), Some(i_uni), Some(i_name), Some(i_tdl), Some(i_fam)) = (
        col_idx("symbol"),
        col_idx("uniprot"),
        col_idx("name"),
        col_idx("tdl"),
        col_idx("family"),
    ) else {
        return Err(TcrdError::Parse { line: 1, msg: "missing required column(s); expected symbol,uniprot,name,tdl,family[,novelty_score]".into() });
    };
    let i_novelty = col_idx("novelty_score");

    let mut index = TcrdIndex::default();
    for (n, line) in lines.enumerate() {
        if line.trim().is_empty() || line.trim_start().starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split(',').map(|s| s.trim()).collect();
        let get = |i: usize| fields.get(i).copied().unwrap_or("").to_string();
        let novelty_score = i_novelty
            .and_then(|i| fields.get(i))
            .and_then(|s| if s.is_empty() { None } else { s.parse().ok() });

        let target = TcrdTarget {
            symbol: get(i_sym),
            uniprot: get(i_uni),
            name: get(i_name),
            tdl: get(i_tdl),
            family: get(i_fam),
            novelty_score,
        };
        if target.symbol.is_empty() {
            return Err(TcrdError::Parse { line: n + 2, msg: "empty symbol field".into() });
        }
        index.by_symbol.insert(target.symbol.to_ascii_uppercase(), target.clone());
        if !target.uniprot.is_empty() {
            index.by_uniprot.insert(target.uniprot.to_ascii_uppercase(), target);
        }
    }
    Ok(index)
}
