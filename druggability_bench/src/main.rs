//! `bench_cryptosite`: score every `--pdb` structure with non-backtracking
//! centrality and four classical baselines, join against `--labels`
//! ground truth, and report per-structure + pooled ROC-AUC.
//!
//! ```text
//! bench_cryptosite --pdb-dir <dir of .pdb files> --labels <labels.csv> \
//!     [--cutoff 8.0] [--min-seq-sep 3] [--out results.json]
//! ```
//!
//! Every `.pdb`/`.ent` file directly under `--pdb-dir` is scored; a
//! structure with no rows in `--labels` is skipped with a warning rather
//! than failing the whole run, since real CryptoSite-derived label sets
//! are usually built incrementally and rarely cover every structure you
//! might have PDB files for.

use druggability_bench::{build_report, evaluate_structure, labels as labels_mod, load_structure, METHODS};
use std::path::PathBuf;
use target_druggability::contact_graph::ContactParams;

struct Args {
    pdb_dir: PathBuf,
    labels_path: PathBuf,
    cutoff: f64,
    min_seq_sep: i32,
    out: Option<PathBuf>,
    seed: u64,
}

fn parse_args() -> Result<Args, String> {
    let mut pdb_dir = None;
    let mut labels_path = None;
    let mut cutoff = 8.0;
    let mut min_seq_sep = 3;
    let mut out = None;
    let mut seed = 42u64;

    let mut it = std::env::args().skip(1);
    while let Some(flag) = it.next() {
        let mut need = |name: &str| it.next().ok_or_else(|| format!("{name} requires a value"));
        match flag.as_str() {
            "--pdb-dir" => pdb_dir = Some(PathBuf::from(need("--pdb-dir")?)),
            "--labels" => labels_path = Some(PathBuf::from(need("--labels")?)),
            "--cutoff" => cutoff = need("--cutoff")?.parse().map_err(|e| format!("{e}"))?,
            "--min-seq-sep" => {
                min_seq_sep = need("--min-seq-sep")?.parse().map_err(|e| format!("{e}"))?
            }
            "--out" => out = Some(PathBuf::from(need("--out")?)),
            "--seed" => seed = need("--seed")?.parse().map_err(|e| format!("{e}"))?,
            other => return Err(format!("unrecognized argument: {other}")),
        }
    }

    Ok(Args {
        pdb_dir: pdb_dir.ok_or("--pdb-dir <dir> is required")?,
        labels_path: labels_path.ok_or("--labels <file.csv> is required")?,
        cutoff,
        min_seq_sep,
        out,
        seed,
    })
}

fn usage() -> &'static str {
    "usage: bench_cryptosite --pdb-dir <dir> --labels <labels.csv> [--cutoff 8.0] [--min-seq-sep 3] [--seed 42] [--out results.json]"
}

fn main() {
    let args = match parse_args() {
        Ok(a) => a,
        Err(e) => {
            eprintln!("error: {e}\n{}", usage());
            std::process::exit(2);
        }
    };

    let label_set = match labels_mod::load_csv(&args.labels_path) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("error loading --labels {}: {e}", args.labels_path.display());
            std::process::exit(1);
        }
    };

    let mut pdb_files: Vec<PathBuf> = match std::fs::read_dir(&args.pdb_dir) {
        Ok(rd) => rd
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .filter(|p| {
                p.extension().and_then(|e| e.to_str()).map(|e| e.eq_ignore_ascii_case("pdb") || e.eq_ignore_ascii_case("ent")).unwrap_or(false)
            })
            .collect(),
        Err(e) => {
            eprintln!("error reading --pdb-dir {}: {e}", args.pdb_dir.display());
            std::process::exit(1);
        }
    };
    pdb_files.sort();

    if pdb_files.is_empty() {
        eprintln!("error: no .pdb/.ent files found in {}", args.pdb_dir.display());
        std::process::exit(1);
    }

    let params = ContactParams { distance_cutoff: args.cutoff, min_seq_separation: args.min_seq_sep };

    let mut structure_results = Vec::new();
    // method -> (pooled scores, pooled labels)
    let mut pooled: Vec<(String, Vec<f64>, Vec<bool>)> =
        METHODS.iter().map(|m| (m.to_string(), Vec::new(), Vec::new())).collect();

    for path in &pdb_files {
        let structure = match load_structure(path) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("[skip] {}: {e}", path.display());
                continue;
            }
        };

        let Some(structure_labels) = label_set.get(&structure.id) else {
            eprintln!("[skip] {}: no rows in --labels for structure_id {:?}", path.display(), structure.id);
            continue;
        };

        match evaluate_structure(&structure, params, structure_labels, args.seed) {
            Ok((result, subset_scores, lbl)) => {
                eprintln!(
                    "[{}] {} residues, {} labeled",
                    result.structure_id, result.n_residues, result.n_labeled
                );
                for m in &result.methods {
                    match m.auc {
                        Some(auc) => eprintln!(
                            "    {:<16} AUC={auc:.4}  (n_pos={}, n_neg={})",
                            m.method, m.n_pos, m.n_neg
                        ),
                        None => eprintln!(
                            "    {:<16} AUC=n/a    (n_pos={}, n_neg={} -- single class, excluded from pooling)",
                            m.method, m.n_pos, m.n_neg
                        ),
                    }
                }

                for (name, subset) in subset_scores {
                    if let Some(slot) = pooled.iter_mut().find(|(m, _, _)| m == name) {
                        slot.1.extend(subset);
                        slot.2.extend(lbl.iter().copied());
                    }
                }

                structure_results.push(result);
            }
            Err(e) => eprintln!("[skip] {}: {e}", path.display()),
        }
    }

    if structure_results.is_empty() {
        eprintln!("error: no structure could be evaluated (check --labels structure_id values match PDB file stems)");
        std::process::exit(1);
    }

    let report = build_report(structure_results, pooled);

    eprintln!("\n=== pooled ROC-AUC across {} structure(s) ===", pdb_files.len());
    for p in &report.pooled {
        match p.auc {
            Some(auc) => eprintln!("{:<16} AUC={auc:.4}  (n_pos={}, n_neg={})", p.method, p.n_pos, p.n_neg),
            None => eprintln!("{:<16} AUC=n/a", p.method),
        }
    }

    let json = serde_json::to_string_pretty(&report).expect("serialization cannot fail here");
    match args.out {
        Some(path) => {
            std::fs::write(&path, &json).unwrap_or_else(|e| {
                eprintln!("error writing {}: {e}", path.display());
                std::process::exit(1);
            });
            eprintln!("\nfull report written to {}", path.display());
        }
        None => println!("{json}"),
    }
}
