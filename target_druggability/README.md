# target_druggability

Non-backtracking (Hashimoto) spectral analysis of protein residue contact
networks, for prioritizing NIH **Illuminating the Druggable Genome (IDG)**
Tdark/Tbio ("difficult to drug") targets for structural follow-up.

Builds on this workspace's existing `nbsc` (Non-Backtracking Spectral
Convolution) and `krylov_ds` crates -- the same Hashimoto/`rho_B` machinery
used elsewhere for graph anomaly detection, applied here to a physical
protein contact network instead of a citation/social graph.

## What it actually does

1. **`pdb`** — parses C-alpha coordinates out of a PDB file.
2. **`contact_graph`** — builds a residue contact network (nodes = residues,
   edges = spatially close, sequence-distant Cα pairs) as an `nbsc::graph::Graph`.
3. **`spectral_druggability`** —
   - `global_coupling_score`: `rho_B`, the Hashimoto spectral radius, as a
     whole-target structural-coupling score (reuses `nbsc::spectral::estimate_spectral_radius` directly).
   - `residue_centrality`: **new** per-residue non-backtracking eigenvector
     centrality, extracted from the vertex-space block of the dominant
     Ritz eigenvector of the same `2n x 2n` Bass-reduced linearization
     `nbsc` already uses for `rho_B` — so no new eigensolver was written,
     just a new read of an existing computation.
   - `cluster_pockets`: connected-component clustering of high-scoring
     residues into candidate contiguous "hotspot" regions.
4. **`tcrd`** — loads NIH IDG Target Development Level (TDL) tier data
   (Tclin/Tchem/Tbio/Tdark) from a flat CSV export.
5. **`report`** — combines both into a ranked JSON report.

## Why non-backtracking centrality (the actual hypothesis being tested)

Plain degree/eigenvector centrality on a contact graph mostly just finds
"buried, densely packed core residues" — not useful for pocket-finding.
The non-backtracking operator suppresses that local-density signal by
construction (it discounts walks that immediately reverse), which shifts
weight toward residues that *bridge* otherwise-distant parts of the fold —
the graph-theoretic signature of a long-range allosteric communication
path or cryptic-pocket-lining region. `tests::bridge_residues_score_higher_than_clique_interior`
is a minimal synthetic check of exactly this claim (two dense cliques
joined by one edge; the bridge residues score higher than clique
interiors).

**This is a structural-network heuristic, not a validated druggability
predictor.** It has not been benchmarked here against a real cryptic-pocket
dataset (e.g. CryptoSite). Treat its output as a ranked hypothesis list —
somewhere to point MD simulation, fragment screening, or a structural
biologist's attention — not a verdict.

## Getting real data (this crate has no network access)

- **PDB structures**: `https://files.rcsb.org/download/<PDB_ID>.pdb`, or a
  local AlphaFold/ColabFold model for a target with no experimental structure
  (common for Tdark targets — that's part of why they're understudied).
- **NIH IDG target tiers**: browse/export from Pharos
  (`https://pharos.nih.gov/targets`), query the Pharos GraphQL API
  (`https://pharos-api.ncats.io/graphql`), or pull the full TCRD relational
  dump (`http://juniper.health.unm.edu/tcrd/download/`) and reshape it into
  the flat CSV schema in `data/sample_tcrd_export.csv`
  (`symbol,uniprot,name,tdl,family,novelty_score`).

`data/sample_tcrd_export.csv` ships with a handful of illustrative rows for
format demonstration only — verify current TDL tiers against Pharos before
using this for real prioritization; they change as new probes/annotations land.

## Usage

```bash
# single target, no TCRD metadata
cargo run --release -p target_druggability -- \
  --pdb my_target.pdb --out report.json

# with NIH IDG tier lookup
cargo run --release -p target_druggability -- \
  --pdb my_target.pdb \
  --tcrd tcrd_export.csv --symbol MYGENE \
  --cutoff 8.0 --min-seq-sep 3 --top-percentile 0.9 \
  --out report.json
```

Batch mode: run once per structure (one PDB file per invocation) and merge
the JSON reports downstream — `rank_targets` in `report.rs` is the sort you
want to apply across the merged set (difficult targets first, then by
top-pocket score).

## Suggested next step: local ΠGPX / GPCR / kinase Tdark batch run

Given the rest of this workspace's focus on brane-tiling spectral geometry,
a natural extension (not built here) is scoring every AlphaFold model for
the ~400 Tdark GPCR/kinase/ion-channel/NR targets IDG tracks, ranking by
`(is_difficult_target, top pocket mean_score)`, and handing the top N to a
domain expert — the same "honest next step" pattern already used in the
Hashimoto/Kasteleyn Newton-polygon work: get the pipeline right on a couple
of concrete cases before generalizing to a batch.
