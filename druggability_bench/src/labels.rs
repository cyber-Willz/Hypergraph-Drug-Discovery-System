//! Ground-truth residue labels for benchmarking, in a CryptoSite-style
//! schema.
//!
//! ## Format
//! CSV, one row per labeled residue:
//! ```text
//! structure_id,chain,res_seq,is_pocket_lining
//! 1a8o,A,154,1
//! 1a8o,A,155,0
//! ```
//! - `structure_id` must match the file stem (or explicit `--id`) used for
//!   the corresponding `--pdb` input to `bench_cryptosite`.
//! - `chain`/`res_seq` must match the PDB numbering in that structure file
//!   (author chain ID and `resSeq`, i.e. exactly [`target_druggability::pdb::Residue`]'s
//!   fields) -- this is what lets labels be joined onto the contact graph's
//!   node indices unambiguously even when there are unresolved-loop gaps.
//! - `is_pocket_lining` is `1` if this residue lines an experimentally
//!   characterized cryptic pocket (e.g. a CryptoSite "positive" residue --
//!   see Cimermancic et al. 2016, *J Mol Biol* 428(4):709-719, and the
//!   curated apo/holo structure pairs at
//!   <https://github.com/bowman-lab/CryptoSite>), `0` otherwise. Residues
//!   present in the structure but absent from the label file are treated
//!   as unlabeled and excluded from AUC (not silently treated as negative
//!   -- CryptoSite negatives are usually a curated "definitely not
//!   pocket" set, not "everything else").
//!
//! This crate does not bundle real CryptoSite data (it isn't network-
//! accessible from here) -- point `--labels` at a CSV built from the
//! actual dataset for a real benchmark run. `data/example_labels.csv`
//! ships a small synthetic label set over the sample 1A8O structure
//! purely so the harness is exercised end-to-end in CI/tests; treat any
//! AUC computed from it as a harness smoke test, not a validation result.

use std::collections::HashMap;
use std::path::Path;

#[derive(Debug, thiserror::Error)]
pub enum LabelError {
    #[error("I/O error reading labels file: {0}")]
    Io(#[from] std::io::Error),
    #[error("labels file line {line}: expected 4 columns (structure_id,chain,res_seq,is_pocket_lining), got {found}")]
    BadColumnCount { line: usize, found: usize },
    #[error("labels file line {line}: could not parse res_seq: {source}")]
    BadResSeq { line: usize, source: std::num::ParseIntError },
    #[error("labels file line {line}: could not parse is_pocket_lining as 0/1: {value:?}")]
    BadLabel { line: usize, value: String },
    #[error("labels file line {line}: chain field must be exactly one character, got {value:?}")]
    BadChain { line: usize, value: String },
}

/// `(chain, res_seq) -> is_pocket_lining`, scoped to one structure.
pub type StructureLabels = HashMap<(char, i32), bool>;

/// All labels loaded from one CSV, keyed by `structure_id`.
#[derive(Debug, Default)]
pub struct LabelSet {
    pub by_structure: HashMap<String, StructureLabels>,
}

impl LabelSet {
    pub fn get(&self, structure_id: &str) -> Option<&StructureLabels> {
        self.by_structure.get(structure_id)
    }
}

pub fn load_csv(path: &Path) -> Result<LabelSet, LabelError> {
    let text = std::fs::read_to_string(path)?;
    let mut set = LabelSet::default();

    for (i, raw_line) in text.lines().enumerate() {
        let line_no = i + 1;
        let line = raw_line.trim();
        if line.is_empty() {
            continue;
        }
        // Skip an optional header row.
        if line_no == 1 && line.to_ascii_lowercase().starts_with("structure_id") {
            continue;
        }
        let cols: Vec<&str> = line.split(',').map(str::trim).collect();
        if cols.len() != 4 {
            return Err(LabelError::BadColumnCount { line: line_no, found: cols.len() });
        }
        let structure_id = cols[0].to_string();
        let chain = {
            let mut chars = cols[1].chars();
            let c = chars.next().ok_or_else(|| LabelError::BadChain {
                line: line_no,
                value: cols[1].to_string(),
            })?;
            if chars.next().is_some() {
                return Err(LabelError::BadChain { line: line_no, value: cols[1].to_string() });
            }
            c
        };
        let res_seq: i32 = cols[2]
            .parse()
            .map_err(|source| LabelError::BadResSeq { line: line_no, source })?;
        let is_pocket_lining = match cols[3] {
            "1" | "true" | "TRUE" | "True" => true,
            "0" | "false" | "FALSE" | "False" => false,
            other => {
                return Err(LabelError::BadLabel { line: line_no, value: other.to_string() })
            }
        };

        set.by_structure
            .entry(structure_id)
            .or_default()
            .insert((chain, res_seq), is_pocket_lining);
    }

    Ok(set)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn parses_with_header_and_without() {
        let mut f = tempfile_with_content(
            "structure_id,chain,res_seq,is_pocket_lining\n1a8o,A,154,1\n1a8o,A,155,0\n1abc,B,10,1\n",
        );
        let set = load_csv(f.path()).unwrap();
        assert_eq!(set.get("1a8o").unwrap().get(&('A', 154)), Some(&true));
        assert_eq!(set.get("1a8o").unwrap().get(&('A', 155)), Some(&false));
        assert_eq!(set.get("1abc").unwrap().get(&('B', 10)), Some(&true));
        f.flush().unwrap();
    }

    fn tempfile_with_content(content: &str) -> tempfile_shim::NamedTempFile {
        let mut f = tempfile_shim::NamedTempFile::new();
        f.write_all(content.as_bytes()).unwrap();
        f
    }

    // Minimal std-only stand-in so this crate doesn't need a `tempfile`
    // dev-dependency just for one test.
    mod tempfile_shim {
        use std::fs::File;
        use std::io::Write;
        use std::path::{Path, PathBuf};

        pub struct NamedTempFile {
            path: PathBuf,
            file: File,
        }
        impl NamedTempFile {
            pub fn new() -> Self {
                let path = std::env::temp_dir()
                    .join(format!("druggability_bench_test_{}.csv", std::process::id()));
                let file = File::create(&path).unwrap();
                Self { path, file }
            }
            pub fn path(&self) -> &Path {
                &self.path
            }
        }
        impl Write for NamedTempFile {
            fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
                self.file.write(buf)
            }
            fn flush(&mut self) -> std::io::Result<()> {
                self.file.flush()
            }
        }
        impl Drop for NamedTempFile {
            fn drop(&mut self) {
                let _ = std::fs::remove_file(&self.path);
            }
        }
    }
}
