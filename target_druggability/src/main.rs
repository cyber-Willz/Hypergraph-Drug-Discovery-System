//! CLI: analyze one or more PDB structures for candidate cryptic-pocket /
//! allosteric-hotspot residues via non-backtracking spectral centrality,
//! optionally cross-referenced against an NIH TCRD/Pharos export to flag
//! which targets are actually understudied (Tdark/Tbio).
//!
//! ```text
//! target_druggability --pdb structure.pdb [--pdb another.pdb ...] \
//!     [--tcrd tcrd_export.csv] [--symbol GENE_SYMBOL] [--uniprot P12345] \
//!     [--cutoff 8.0] [--min-seq-sep 3] [--top-percentile 0.9] \
//!     [--out report.json]
//! ```
//!
//! `--symbol`/`--uniprot` apply to the *first* `--pdb` given (single-target
//! mode). For batch mode, run once per structure and merge the JSON
//! reports downstream -- keeping this CLI simple rather than inventing a
//! manifest format no one asked for.

use std::path::PathBuf;
use target_druggability::{contact_graph, pdb, report, spectral_druggability, tcrd};

struct Args {
    pdb_paths: Vec<PathBuf>,
    tcrd_path: Option<PathBuf>,
    symbol: Option<String>,
    uniprot: Option<String>,
    cutoff: f64,
    min_seq_sep: i32,
    top_percentile: f64,
    out: Option<PathBuf>,
}

fn parse_args() -> Result<Args, String> {
    let mut a = Args {
        pdb_paths: Vec::new(),
        tcrd_path: None,
        symbol: None,
        uniprot: None,
        cutoff: 8.0,
        min_seq_sep: 3,
        top_percentile: 0.9,
        out: None,
    };
    let mut it = std::env::args().skip(1);
    while let Some(flag) = it.next() {
        let mut need = |name: &str| it.next().ok_or_else(|| format!("{name} requires a value"));
        match flag.as_str() {
            "--pdb" => a.pdb_paths.push(PathBuf::from(need("--pdb")?)),
            "--tcrd" => a.tcrd_path = Some(PathBuf::from(need("--tcrd")?)),
            "--symbol" => a.symbol = Some(need("--symbol")?),
            "--uniprot" => a.uniprot = Some(need("--uniprot")?),
            "--cutoff" => a.cutoff = need("--cutoff")?.parse().map_err(|e| format!("{e}"))?,
            "--min-seq-sep" => {
                a.min_seq_sep = need("--min-seq-sep")?.parse().map_err(|e| format!("{e}"))?
            }
            "--top-percentile" => {
                a.top_percentile = need("--top-percentile")?.parse().map_err(|e| format!("{e}"))?
            }
            "--out" => a.out = Some(PathBuf::from(need("--out")?)),
            other => return Err(format!("unrecognized argument: {other}")),
        }
    }
    if a.pdb_paths.is_empty() {
        return Err("at least one --pdb <path> is required".into());
    }
    Ok(a)
}

fn main() {
    let args = match parse_args() {
        Ok(a) => a,
        Err(e) => {
            eprintln!("error: {e}");
            eprintln!(
                "usage: target_druggability --pdb <file.pdb> [--tcrd <export.csv>] [--symbol SYM] [--uniprot ID] [--cutoff 8.0] [--min-seq-sep 3] [--top-percentile 0.9] [--out report.json]"
            );
            std::process::exit(2);
        }
    };

    let tcrd_index = match &args.tcrd_path {
        Some(p) => match tcrd::load_csv(p) {
            Ok(idx) => Some(idx),
            Err(e) => {
                eprintln!("error loading --tcrd {}: {e}", p.display());
                std::process::exit(1);
            }
        },
        None => None,
    };

    let mut reports = Vec::new();

    for (i, pdb_path) in args.pdb_paths.iter().enumerate() {
        let structure = match pdb::parse_ca_atoms_file(pdb_path) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("error parsing {}: {e}", pdb_path.display());
                std::process::exit(1);
            }
        };
        eprintln!(
            "[{}] parsed {} C-alpha residues from {}",
            structure.id,
            structure.residues.len(),
            pdb_path.display()
        );

        let params = contact_graph::ContactParams {
            distance_cutoff: args.cutoff,
            min_seq_separation: args.min_seq_sep,
        };
        let cg = contact_graph::build_contact_graph(&structure, params);
        eprintln!(
            "[{}] contact graph: {} nodes, {} edges",
            structure.id,
            cg.graph.n,
            cg.graph.m()
        );

        let rho_b = spectral_druggability::global_coupling_score(&cg.graph, 42);

        let krylov_dim = (2 * cg.graph.n).min(200).max(4);
        let scores = match spectral_druggability::residue_centrality(&cg.graph, krylov_dim, 42) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("error computing spectral centrality for {}: {e}", pdb_path.display());
                std::process::exit(1);
            }
        };

        let pockets = spectral_druggability::cluster_pockets(
            &cg.graph,
            &scores,
            args.top_percentile,
            3, // min_cluster_size
        );
        eprintln!("[{}] {} candidate pocket cluster(s) found", structure.id, pockets.len());

        // TCRD lookup only applies unambiguously to the first structure
        // when --symbol/--uniprot are given (see module doc).
        let tcrd_target = if i == 0 {
            args.symbol
                .as_deref()
                .and_then(|s| tcrd_index.as_ref().and_then(|idx| idx.lookup_symbol(s)))
                .or_else(|| {
                    args.uniprot
                        .as_deref()
                        .and_then(|u| tcrd_index.as_ref().and_then(|idx| idx.lookup_uniprot(u)))
                })
        } else {
            None
        };

        let target_report = report::build_report(
            &structure.id,
            &cg,
            &scores,
            rho_b,
            &pockets,
            tcrd_target,
        );
        reports.push(target_report);
    }

    let ranked = report::rank_targets(reports);
    let json = serde_json::to_string_pretty(&ranked).expect("serialization cannot fail here");

    match args.out {
        Some(path) => {
            std::fs::write(&path, &json).unwrap_or_else(|e| {
                eprintln!("error writing {}: {e}", path.display());
                std::process::exit(1);
            });
            eprintln!("report written to {}", path.display());
        }
        None => println!("{json}"),
    }
}
