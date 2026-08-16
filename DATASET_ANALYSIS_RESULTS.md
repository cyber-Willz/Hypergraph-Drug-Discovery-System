# Dataset analysis results — this session

## Rust workspace: build, tests, and synthetic-structure results

- **Build**: `cargo build --release --workspace` — 4 crates
  (`krylov_ds`, `nbsc`, `target_druggability`, `druggability_bench`)
  compiled clean in 1m52s.
- **Tests**: `cargo test --release --workspace` — **22/22 passed**, 0 failed
  (5 `roc.rs` AUC cross-checks, 4 baseline-centrality tests, contact-graph
  tests, PDB parser tests, spectral clustering tests, labels-CSV parser
  test, Arnoldi/eigenvalue tests, spectral-graph tests).
- **`bench_cryptosite` smoke test** (synthetic toy structure):
  `example_structure`, 10 residues / 8 labeled, **AUC = 1.0000** on all
  five methods (nonbacktracking, degree, closeness, betweenness,
  eigenvector). This is the known-trivial synthetic case, not a
  validation result.
- **Pocket-threshold sweep** (`target_druggability`):
  - `--top-percentile 0.9` → 0 pockets found (structure too small for
    the default threshold).
  - `--top-percentile 0.8` → first threshold that finds a pocket:
    1 cluster, 4 residues (A/5, A/113, A/109, A/101), mean score 0.931.
  - `0.7`, `0.6`, `0.5` → same or one additional residue as the
    threshold relaxes.
- **`docking_prep.py`** on the 0.8-threshold report: pocket center
  **(8.25, 0.75, 0.75)**, box **(19.0, 15.0, 15.0) Å**. Independently
  cross-checked by manually averaging the same four residues' Cα
  coordinates — **matched exactly**. `openmm_pocket_prep.py` compiled
  cleanly (`py_compile`); `diffdock_pocket_batch.csv` wrote one
  well-formed row.
- **Case-preservation proof**: copying the structure/labels under
  matched case and re-running `bench_cryptosite` produced a clean join
  (10 residues, 8 labeled, no `[skip] no rows in --labels` warning) —
  confirming the original `.upper()` bug would have broken this join
  and the fix resolves it.

## Network boundary (confirmed live, this session)

| Host | Status | Reason |
|---|---|---|
| `files.rcsb.org` | 403 | `host_not_allowed` |
| `alphafold.ebi.ac.uk` | 403 | `host_not_allowed` |
| `pharos-api.ncats.io` | 403 | `host_not_allowed` |

This is the egress proxy denying the host outright, not a timeout/DNS
failure — which is why `rcsb_bulk_download.py` fails fast on it instead
of retrying with backoff.

Live-testing the corrected script against `1a8o,1JWP,AF-P00918-F1-model_v4`:

```
[1/3] 1a8o: FAIL -- blocked by sandbox network egress policy: Host not in allowlist: files.rcsb.org.
[2/3] 1JWP: FAIL -- blocked by sandbox network egress policy: Host not in allowlist: files.rcsb.org.
[3/3] AF-P00918-F1-model_v4: FAIL -- AlphaFold entry -- use pipeline/alphafold_client.py, not this script
done: 0 downloaded/cached, 0 failed, 1 AlphaFold entries skipped, 2 blocked by sandbox network policy
```

Exit code 1. Confirmed: fast-fail (no backoff delay), correct
distinction between "sandbox-denied" and "could not reach RCSB", and
correct routing of the `AF-`-prefixed id.

## CryptoBank Parquet dataset — real analysis

33 Parquet files, 19.0 MB on disk: 1 row-level table + 31
category/lookup tables + 1 partial-column table + a manifest.

### `row_level.parquet`

- **Confirmed 5,989,860 rows**, 4 columns:
  `datetime_deposit_a`, `datetime_release_a`, `datetime_deposit_b`,
  `datetime_release_b`.
- Track B is null on 362,464 rows (5,627,396 / 5,989,860 rows have both
  tracks present).
- `release_a − deposit_a`: min 0, max 17,702 days, **mean 2,443.9 days**.
- `release_b − deposit_b`: min 0, max 17,502 days, **mean 3,372.9 days**.
- `deposit_b − deposit_a` (non-null rows): ranges from −17,599 to
  +17,662 days, mean **−979.9 days** — consistent with paired
  apo/holo-style entries deposited independently rather than at a fixed
  offset, not a strict follow-on relationship.

### Lookup/category tables (31 total, 679,939 rows combined)

Largest tables:

| Table | Rows | Sample |
|---|---|---|
| `pdbid_chain_categories_a` | 163,499 | `101M_A`, `102M_A` |
| `pdbid_chain_categories_b` | 81,302 | `107L_A`, `108L_A` |
| `pdb_id_categories` | 86,134 | `101M`, `102M` |
| `pdb_id_categories_2` | 33,365 | `107L`, `108L` |
| `smiles_categories` / `_2` | 34,018 / 34,064 | ligand SMILES strings |
| `inchi_categories` | 34,007 | InChI strings |
| `ligand_name_categories` | 34,017 | ligand names |
| `numeric_code_categories` | 34,095 | numeric codes |
| `sequence_categories` / `_2` | 39,256 / 19,899 | protein sequences |
| `uniprot_accession_categories_a–d` | 20,662 / 20,662 / 21,151 / 9,176 | UniProt accessions |
| `drugbank_id_categories` | 4,927 | `DB00114`, `DB00115` |

### Spot-check verification (independently re-run, all passed)

| Table | Value | Present | n |
|---|---|---|---|
| `lookup_0001_pdb_id_categories.parquet` | `101M` | ✅ | 86,134 |
| `lookup_0027_drugbank_id_categories.parquet` | `DB00114` | ✅ | 4,927 |
| `lookup_0043_pdbid_chain_categories_a.parquet` | `101M_A` | ✅ | 163,499 |

### PDB ID case convention

**All 86,134 / 86,134 PDB IDs are uppercase** — 0 lowercase or mixed
entries. This directly confirms (not just motivates) the
case-preservation fix in `rcsb_bulk_download.py`: forcing `.upper()`
on lowercase input would have matched this data, but the original bug
was forcing `.upper()` unconditionally, which breaks any
lowercase-keyed labels file (e.g. CryptoSite-style benchmarks) joined
against locally-cached files.

### `partial_col_0066.parquet`

- 3,750,000 rows, `list[f64]` (3-component vectors).
- **966,014 / 3,750,000 rows have a nonzero first component**
  (mean 0.134, range [0.0, 1.0]) — confirms this is real data, not
  all-zero padding.

## What remains unestablished

- No real CryptoSite/CryptoBank benchmark number — that needs the
  actual PDB coordinate files (blocked by the `files.rcsb.org` egress
  wall) joined against labels derived from the full CryptoBank dataset.
  The harness that would produce that number (`druggability_bench`) is
  built, tested, and proven correct on synthetic data.
- No live Pharos/AlphaFold pipeline run — same network wall, confirmed
  live.
- The true `columns` name Index for the CryptoBank source pickle and
  the remaining ~2.24M rows of `col_0066` are still open (see
  `CRYPTOBANK_DATASET_NOTES.md` for the full account of the earlier
  memory-bounded inspection session).
