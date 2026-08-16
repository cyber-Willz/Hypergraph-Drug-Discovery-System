# Filter-to-dock pipeline

Implements:

```
Pharos GraphQL API --> AlphaFold PDB --> target_druggability --> OpenMM / DiffDock-Pocket prep
```

## Setup

```bash
cd pipeline
pip install -r requirements.txt
# target_druggability's release binary is built automatically on first
# run (needs `cargo`); or build it yourself once:
cargo build --release -p target_druggability --manifest-path ../Cargo.toml
```

Requires network access to `pharos-api.ncats.io` and `alphafold.ebi.ac.uk`.

## Usage

**Single target** (you already know the gene symbol or UniProt accession):
```bash
python pipeline.py --symbol CA2 --out-dir runs/ca2
```

**Disease-driven discovery + batch triage** (the full roadmap loop —
find targets, keep the understudied ones, score them all):
```bash
python pipeline.py --disease "pancreatic cancer" --top 20 \
    --difficult-only --out-dir runs/panc_cancer
```

Each target gets its own `runs/<label>/<SYMBOL>/` with:
- `<uniprot>.pdb` — the downloaded AlphaFold model
- `tcrd_row.csv` — that target's Pharos TDL/family metadata, in the CSV
  schema `target_druggability --tcrd` expects
- `report.json` — the full `target_druggability` output
- `pocket_grid.json` — top pocket's center/box + residue list (only if a
  pocket was found)
- `openmm_pocket_prep.py` — ready-to-run OpenMM restraint script that
  relaxes the pocket region while holding the rest of the structure fixed
  (needs `pip install openmm pdbfixer` separately — not a dependency of
  this pipeline itself, since most triage runs won't need it)

Plus, at the top level:
- `summary.json` — one row per successfully processed target, sorted
  difficult-targets-with-pockets first
- `diffdock_pocket_batch.csv` — one row per target with a pocket, pocket
  center included; **edit the `ligand_description` column** (it's a
  placeholder) before feeding this to DiffDock-Pocket, and verify the
  header matches your DiffDock-Pocket checkout's expected columns (see
  `docking_prep.py`'s docstring — this is the one place in the pipeline
  where an upstream tool's exact CSV schema wasn't verifiable from here)

## What each stage actually does, and its real caveats

1. **Pharos GraphQL** (`pharos_client.py`) — `target(q:{sym|uniprot})` for
   single lookups, `targets(filter:{associatedDisease})` for discovery.
   Confirmed against Pharos' own published example queries at the time
   this was written; Pharos' schema does evolve, so if a query starts
   failing, check https://pharos.nih.gov/api's live Playground before
   assuming this code is broken.
2. **AlphaFold DB** (`alphafold_client.py`) — `GET /api/prediction/{uniprot}`.
   These are *predicted* structures. Low-pLDDT regions are often
   genuinely disordered, not just "uncertain but real" — a pocket flagged
   in one is weaker evidence than the same call on a crystallographic apo
   structure. This client doesn't currently pull per-residue pLDDT; if
   you're triaging Tdark targets that only have AlphaFold coverage,
   fetch the PAE/confidence file too (`paeDocUrl` in the raw API
   response) before trusting a pocket call in a specific region.
3. **target_druggability** — see its own module docs and
   `../druggability_bench` for what non-backtracking centrality is and
   isn't validated to do. This pipeline doesn't change or bypass any of
   that; it just automates getting structures into it.
4. **Docking/MD prep** (`docking_prep.py`) — pocket center is the mean
   Cα position of the top pocket's residues, not a true cavity centroid
   from probe-based pocket detection (fpocket/CASTp). Good enough to
   point a pocket-aware docking tool or an OpenMM restraint selection at
   the right neighborhood; not a substitute for visually inspecting the
   pocket (PyMOL/ChimeraX) before committing compute to a screen.

## What this pipeline does not do

- It does not run DiffDock-Pocket or OpenMM itself — both are heavy,
  GPU-shaped dependencies out of scope to auto-install here. It produces
  correctly-computed inputs for them and stops.
- It does not validate that a flagged pocket is real (that's what
  `../druggability_bench` is for, against labeled data).
- Batch mode does not paginate beyond `--top`; Pharos' disease-association
  ranking determines which targets you see first.
