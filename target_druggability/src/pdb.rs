//! Minimal PDB parser: pulls C-alpha coordinates out of `ATOM` records.
//!
//! This deliberately does not try to be a general-purpose PDB library (no
//! HETATM/ligand handling, no altloc resolution beyond "first one wins", no
//! mmCIF support). It reads exactly what the contact-graph builder in
//! [`crate::contact_graph`] needs: one 3D point per residue, in chain/seq
//! order, plus enough identity metadata (chain, residue number, residue
//! name) to report results back in terms a biologist can act on.
//!
//! Fetching structures: this crate has no network access baked in. Download
//! PDB files yourself from RCSB, e.g.
//! `https://files.rcsb.org/download/<PDB_ID>.pdb`, or from a local
//! AlphaFold/ColabFold model, and pass the path on the command line.

use std::collections::HashSet;
use std::fmt;

#[derive(Debug, Clone, PartialEq)]
pub struct Residue {
    pub chain: char,
    /// Author-assigned residue sequence number (PDB `resSeq` field), not a
    /// dense 0-based index -- there can be gaps from unresolved loops.
    pub res_seq: i32,
    pub res_name: String,
    pub ca: [f64; 3],
}

#[derive(Debug, Clone, Default)]
pub struct Structure {
    pub id: String,
    /// Residues in file order (chain-major, then sequence order), one entry
    /// per residue that had a resolved C-alpha.
    pub residues: Vec<Residue>,
}

#[derive(Debug)]
pub enum PdbError {
    Empty,
    Io(std::io::Error),
}

impl fmt::Display for PdbError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PdbError::Empty => write!(f, "no C-alpha atoms found in PDB input"),
            PdbError::Io(e) => write!(f, "I/O error reading PDB file: {e}"),
        }
    }
}
impl std::error::Error for PdbError {}
impl From<std::io::Error> for PdbError {
    fn from(e: std::io::Error) -> Self {
        PdbError::Io(e)
    }
}

/// Parse C-alpha residues out of raw PDB text.
///
/// Column layout follows the standard fixed-width PDB `ATOM` record:
/// atom name in columns 13-16, altLoc in 17, resName 18-20, chainID 22,
/// resSeq 23-26, coordinates in 31-54 (8.3f each). We only keep records
/// where the atom name (trimmed) is exactly `CA` and, if an altloc is
/// present, only the first one encountered per residue.
pub fn parse_ca_atoms(text: &str, id: &str) -> Result<Structure, PdbError> {
    let mut residues = Vec::new();
    let mut seen: HashSet<(char, i32)> = HashSet::new();

    for line in text.lines() {
        if !(line.starts_with("ATOM") || line.starts_with("HETATM")) {
            continue;
        }
        if line.len() < 54 {
            continue;
        }
        let atom_name = line.get(12..16).unwrap_or("").trim();
        if atom_name != "CA" {
            continue;
        }
        let res_name = line.get(17..20).unwrap_or("").trim().to_string();
        let chain = line.get(21..22).unwrap_or(" ").chars().next().unwrap_or(' ');
        let res_seq: i32 = match line.get(22..26).unwrap_or("").trim().parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        if !seen.insert((chain, res_seq)) {
            continue; // altloc duplicate: keep the first
        }
        let parse_coord = |s: Option<&str>| -> Option<f64> { s?.trim().parse().ok() };
        let x = parse_coord(line.get(30..38));
        let y = parse_coord(line.get(38..46));
        let z = parse_coord(line.get(46..54));
        let (Some(x), Some(y), Some(z)) = (x, y, z) else { continue };

        residues.push(Residue { chain, res_seq, res_name, ca: [x, y, z] });
    }

    if residues.is_empty() {
        return Err(PdbError::Empty);
    }

    Ok(Structure { id: id.to_string(), residues })
}

pub fn parse_ca_atoms_file(path: &std::path::Path) -> Result<Structure, PdbError> {
    let text = std::fs::read_to_string(path)?;
    let id = path.file_stem().map(|s| s.to_string_lossy().to_string()).unwrap_or_default();
    parse_ca_atoms(&text, &id)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_minimal_ca_record() {
        let line = "ATOM      1  CA  ALA A   1      11.104  13.207   2.500  1.00 20.00           C  ";
        let s = parse_ca_atoms(line, "test").unwrap();
        assert_eq!(s.residues.len(), 1);
        assert_eq!(s.residues[0].chain, 'A');
        assert_eq!(s.residues[0].res_seq, 1);
        assert_eq!(s.residues[0].res_name, "ALA");
        assert!((s.residues[0].ca[0] - 11.104).abs() < 1e-9);
    }

    #[test]
    fn skips_non_ca_atoms() {
        let text = "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00 20.00           N  \n\
                    ATOM      2  CA  ALA A   1      11.104  13.207   2.500  1.00 20.00           C  \n\
                    ATOM      3  CB  ALA A   1      12.000  14.000   3.000  1.00 20.00           C  ";
        let s = parse_ca_atoms(text, "test").unwrap();
        assert_eq!(s.residues.len(), 1);
    }

    #[test]
    fn empty_input_errors() {
        assert!(parse_ca_atoms("", "empty").is_err());
    }
}
