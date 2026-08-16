# `druggability_bench`

Validates `target_druggability`'s non-backtracking residue centrality
against four classical centrality baselines (degree, closeness,
betweenness, plain/backtracking eigenvector) on a labeled cryptic-pocket
dataset, reporting ROC-AUC per method, per structure and pooled.

## What's real here and what isn't

- The ROC-AUC math (`src/roc.rs`), the baseline centrality implementations
  (`src/baselines.rs`), and the harness plumbing (`src/lib.rs`,
  `src/main.rs`) are complete, tested, and not dataset-specific.
- **`data/example_structure.pdb` and `data/example_labels.csv` are
  synthetic** (two small residue clusters bridged by one interface
  residue, built by hand — not a real protein, not CryptoSite). They
  exist only so `cargo test` and a first `bench_cryptosite` run exercise
  the full pipeline without a network connection. Any AUC number you get
  from them is a smoke test, not a validation result — do not cite it.
- This crate has no bundled real CryptoSite data and cannot fetch it (no
  network access baked in, same constraint as `target_druggability`
  itself). Getting a real number requires the steps below.

## Running a real CryptoSite benchmark

1. **Get the dataset.** CryptoSite (Cimermancic et al. 2016, *J Mol Biol*
   428(4):709–719, PMID 26854760) defines a curated set of apo structures
   with experimentally known cryptic pockets (identified by comparing apo
   vs. holo/ligand-bound structures). The curated structure list and
   pocket-residue annotations are distributed with the paper's
   supplementary material and have been re-packaged by later benchmark
   papers (e.g. PocketMiner, CryptoBench) — search for the specific
   release you want; there is no single canonical machine-readable file
   this crate can point you to without risking a stale/wrong URL.
2. **Get structures.** For each CryptoSite PDB ID, download the apo
   structure, e.g. `https://files.rcsb.org/download/<PDB_ID>.pdb`, into a
   directory (`--pdb-dir`).
3. **Build the label CSV.** For each structure, write one row per residue
   you have ground truth for:
   ```text
   structure_id,chain,res_seq,is_pocket_lining
   1jwp,A,42,1
   1jwp,A,43,0
   ```
   `structure_id` must equal the PDB file's stem (`1jwp.pdb` → `1jwp`).
   `chain`/`res_seq` must match the PDB file's own numbering (author
   chain ID / `resSeq`) — see [`src/labels.rs`] for the exact join
   semantics, including why unlabeled residues are excluded from AUC
   rather than treated as negatives.
4. **Run it:**
   ```bash
   cargo build --release -p druggability_bench
   ./target/release/bench_cryptosite \
       --pdb-dir path/to/cryptosite_structures \
       --labels path/to/cryptosite_labels.csv \
       --cutoff 8.0 --min-seq-sep 3 \
       --out cryptosite_results.json
   ```
   This prints per-structure and pooled AUC for all five methods to
   stderr and writes the full per-structure/per-method breakdown as JSON.

## Interpreting results

- **Pooled AUC** (all labeled residues across all structures, one ROC
  curve) is the headline number — more stable than averaging noisy
  per-structure AUCs when individual structures have few labeled
  residues.
- AUC = 0.5 is chance; AUC = 1.0 is perfect ranking of positives above
  negatives.
- The `eigenvector` baseline (plain, backtracking-allowed) is the most
  informative single comparison: it isolates whether *specifically*
  suppressing backtracking walks helps, versus eigenvector-style
  centrality in general.
- A structure contributes `AUC=n/a` for a method if it happens to have
  only positive or only negative labeled residues (ROC-AUC is undefined
  with one class) — those structures are still shown for transparency,
  and are naturally handled correctly in the pooled AUC since pooling
  happens before the single-class check.
- Non-backtracking centrality outperforming the baselines here would be
  evidence for (not proof of) the crate's design rationale; underperforming
  would be a real finding worth investigating, not something to explain
  away. Report whatever the numbers say.

## Tests

```bash
cargo test --release -p druggability_bench
```

Covers: ROC-AUC against perfect/inverted/tied/random-case brute-force
cross-checks, each baseline centrality on graphs with known structure
(path center has max betweenness, star hub has max eigenvector/closeness,
regular cycle has uniform degree), and the labels CSV parser.
