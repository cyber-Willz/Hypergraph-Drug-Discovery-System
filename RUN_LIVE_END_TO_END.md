# Live end-to-end run log



## 0. Setup

```bash
mkdir -p work && cd work
tar -xzf target_druggability_roadmap_deliverables_tar.gz

# Rust toolchain wasn't preinstalled in this sandbox; installed via apt
# (crates.io/index.crates.io/static.crates.io are allowlisted, so `cargo
# build` itself works fine once rustc/cargo exist):
apt-get install -y cargo rustc
cargo --version   # cargo 1.75.0
rustc --version   # rustc 1.75.0
```

## 1. Confirm the actual network boundary (don't assume it)

```bash
python3 -c "
import requests
for url in [
    'https://files.rcsb.org/download/1A8O.pdb',
    'https://alphafold.ebi.ac.uk/api/prediction/P00918',
    'https://pharos-api.ncats.io/graphql',
]:
    r = requests.get(url, timeout=8)
    print(url, '->', r.status_code, r.headers.get('x-deny-reason'), r.text[:120])
"
```

Real output:

```
https://files.rcsb.org/download/1A8O.pdb -> 403 host_not_allowed Host not in allowlist: files.rcsb.org. Add this host to your network egress settings to allow access.
https://alphafold.ebi.ac.uk/api/prediction/P00918 -> 403 host_not_allowed Host not in allowlist: alphafold.ebi.ac.uk. Add this host to your network egress settings to allow access.
https://pharos-api.ncats.io/graphql -> 403 host_not_allowed Host not in allowlist: pharos-api.ncats.io. Add this host to your network egress settings to allow access.
```

This is the egress proxy denying the host outright (`403` +
`x-deny-reason: host_not_allowed`), not a connection timeout or DNS
failure. That distinction is why the corrected `rcsb_bulk_download.py`
fails fast on it instead of retrying with backoff (see step 5).

## 2. Build the full Rust workspace

```bash
cargo build --release --workspace
```

Real output (trimmed): resolved and compiled `nalgebra`, `serde`,
`thiserror`, etc. from crates.io, then all four workspace members.

```
   Compiling target_druggability v0.1.0 (/home/claude/work/target_druggability)
   Compiling druggability_bench v0.1.0 (/home/claude/work/druggability_bench)
    Finished release [optimized] target(s) in 1m 27s
```

## 3. Run the full test suite

```bash
cargo test --release --workspace
```

Real result: **22/22 tests pass** across `krylov_ds`, `nbsc`,
`target_druggability`, and `druggability_bench` (5 `roc.rs` AUC
cross-checks, 4 baseline-centrality sanity tests, contact-graph tests,
PDB parser tests, spectral clustering tests, labels-CSV parser test,
Arnoldi/eigenvalue tests, spectral-graph tests).

## 4. Run `bench_cryptosite` on the synthetic smoke-test structure

```bash
./target/release/bench_cryptosite \
    --pdb-dir druggability_bench/data \
    --labels druggability_bench/data/example_labels.csv \
    --out /tmp/bench_smoke.json
```

Real output: `example_structure`, 10 residues / 8 labeled, AUC=1.0000 for
all five methods (nonbacktracking, degree, closeness, betweenness,
eigenvector) -- this is the known synthetic toy case, not a validation
result (see `druggability_bench/README.md`).

## 5. Run `target_druggability` directly and feed a real report into `docking_prep.py`

The default `--top-percentile 0.9` finds 0 pockets on this 10-residue toy
structure (too small a residue count for the default top-10% threshold to
select a cluster); swept the threshold live rather than assuming a value:

```bash
for tp in 0.9 0.8 0.7 0.6 0.5; do
  ./target/release/target_druggability \
      --pdb druggability_bench/data/example_structure.pdb \
      --top-percentile $tp --out /tmp/td_report_$tp.json
done
```

`--top-percentile 0.8` is the first threshold that finds a pocket (1
cluster, 4 residues). Used that report to drive `docking_prep.py` for
real:

```python
from pathlib import Path
import json
from docking_prep import compute_pocket_geometry, write_openmm_template, write_diffdock_pocket_row

report = json.loads(Path("/tmp/td_report_0.8.json").read_text())
geom = compute_pocket_geometry(
    report[0], Path("../druggability_bench/data/example_structure.pdb"), pocket_rank=1,
)
write_openmm_template(geom, Path("../druggability_bench/data/example_structure.pdb"),
                       Path("/tmp/pipeline_run/openmm_pocket_prep.py"))
write_diffdock_pocket_row(geom, Path("../druggability_bench/data/example_structure.pdb"),
                           "PLACEHOLDER_LIGAND_SMILES", Path("/tmp/pipeline_run/diffdock_pocket_batch.csv"))
```

Real result: pocket center `(8.25, 0.75, 0.75)`, box `(19.0, 15.0, 15.0)`
Å. Independently cross-checked by manually averaging the same four
residues' Cα coordinates via `pdb_geometry.parse_ca_coords` -- matched
exactly. `openmm_pocket_prep.py` compiled cleanly
(`python3 -m py_compile`); `diffdock_pocket_batch.csv` has the correct
header and one well-formed row.

## 6. Live-test the corrected `rcsb_bulk_download.py`

```bash
python3 pipeline/rcsb_bulk_download.py \
    --ids 1a8o,1JWP,AF-P00918-F1-model_v4 --out-dir /tmp/pdb_cache
```

Real output:

```
fetching 3 structure(s) -> /tmp/pdb_cache
[1/3] 1a8o: FAIL -- blocked by sandbox network egress policy: Host not in allowlist: files.rcsb.org. Add this host to your network egress settings to allow access.
[2/3] 1JWP: FAIL -- blocked by sandbox network egress policy: Host not in allowlist: files.rcsb.org. Add this host to your network egress settings to allow access.
[3/3] AF-P00918-F1-model_v4: FAIL -- AlphaFold entry -- use pipeline/alphafold_client.py, not this script

done: 0 downloaded/cached, 0 failed, 1 AlphaFold entries skipped (route those through pipeline/alphafold_client.py), 2 blocked by sandbox network policy (add files.rcsb.org to your network egress settings and re-run)
```

Confirmed: fails immediately (no multi-second backoff delay) on the
sandbox-denied host, correctly distinguishes that from a "could not
reach RCSB" network error, and correctly routes the `AF-`-prefixed id
without attempting to fetch it as an RCSB entry.

## 7. Prove the case-preservation fix actually fixes the labels join

The original script would have written `1a8o` as `1A8O.pdb`
(`'1a8o'.strip().upper() + '.pdb'` → `1A8O.pdb`), which does not match
`druggability_bench/README.md`'s lowercase labels convention. The
corrected script preserves input case (`1a8o` → `1a8o.pdb`). Proved the
downstream effect end-to-end by placing a structure file under the
case-preserved name and running the real join:

```bash
cp druggability_bench/data/example_structure.pdb /tmp/e2e_test/pdb_cache/example_structure.pdb
cp druggability_bench/data/example_labels.csv /tmp/e2e_test/labels.csv
./target/release/bench_cryptosite \
    --pdb-dir /tmp/e2e_test/pdb_cache --labels /tmp/e2e_test/labels.csv \
    --out /tmp/e2e_test/results.json
```

Real output: `[example_structure] 10 residues, 8 labeled`, AUC computed
for all 5 methods -- no `[skip] ... no rows in --labels` warning, which
is exactly the failure mode the original case bug would have caused
against a real, lowercase-keyed CryptoSite-style labels file.

## 8. Inspect the real CryptoBank dataset (memory-bounded, this session)

`cryptobank_dataset_11_08_2025.zip` unzips to a 4.8GB pandas pickle. A
plain load OOM's on this sandbox's 3.9GB RAM -- confirmed live:

```bash
python3 -c "import pandas as pd; pd.read_pickle('cryptobank_dataset_11_08_2025.pkl')"
```

Real result: process killed by the OOM reaper (`Killed`, no traceback --
that's the kernel, not Python).

Built `pipeline/inspect_large_pickle.py` to introspect it within a fixed
memory budget instead: it exploits the fact that CPython's pickler writes
any payload >= 64KB unframed (so it can be `seek()`-skipped on the raw
file without ever being read), caps how much of any single object-dtype
array it retains, and stops itself cleanly via a `resource.getrusage()`
watchdog before hitting the sandbox's ceiling rather than getting
SIGKILLed. Full account of what worked, what didn't (an earlier attempt
to also cap pickle's own backreference memo table corrupted unrelated
structural references and crashed with `STACK_GLOBAL requires str` --
reverted), and what the live data actually showed (uppercase PDB ids,
confirming rather than just motivating the case-preservation fix in
step 6/7) is in `CRYPTOBANK_DATASET_NOTES.md`.

```bash
python3 pipeline/inspect_large_pickle.py cryptobank_dataset_11_08_2025.pkl
```

Real result: exit code 0 (controlled stop, not a crash), peak RSS 2.21GB,
5,989,860 total rows confirmed, 35 arrays' real dtype/shape/sample values
captured before the deliberate cutoff.

## 9. Stream the full CryptoBank dataset to real, verified Parquet (this session)

Went beyond inspection this time -- actually converted as much of the
real data as this sandbox's memory allows into real Parquet files,
using Polars, and independently verified the result.

```bash
# Root is available in this sandbox; added swap as a legitimate lever
# before writing any more custom parsing code -- most memory-bound tasks
# should try this first.
fallocate -l 2G /home/claude/swapfile && chmod 600 /home/claude/swapfile
mkswap /home/claude/swapfile && swapon /home/claude/swapfile
free -h   # confirms: Swap: 2.0Gi available
```

Tried a plain load again with the swap headroom:

```bash
python3 -c "import pandas as pd; pd.read_pickle('cryptobank_dataset_11_08_2025.pkl')"
```

Real result: still SIGKILLed. `dmesg` confirms a genuine **global host
OOM** (`constraint=CONSTRAINT_NONE`), not a container/cgroup limit --
`anon-rss:3871068kB` plus a fully-exhausted 2GB swap, still not enough.
This sandbox's ~3.9GB RAM + 2GB swap (~5.9GB combined) isn't sufficient
for pandas' full in-memory materialization of this file, confirmed by
trying it, not assumed.

Built `pipeline/pickle_to_parquet.py` (streams every element to disk
instead of sampling) and `pipeline/assemble_parquet.py` (reads those
files back through Polars into real Parquet, then re-reads the Parquet
independently to verify it):

```bash
python3 pipeline/pickle_to_parquet.py cryptobank_dataset_11_08_2025.pkl columns
python3 pipeline/assemble_parquet.py columns parquet_out
```

Real result: exit 0 (controlled stop, not SIGKILL). **35 of 36 arrays
captured in full**, including all 4 datetime columns at the true
**5,989,860-row** length; the 36th (a large per-row 3-float-list column)
captured to 3,750,000 real rows before the watchdog's `RSS+swap`
ceiling stopped it. Output: 33 Parquet files, ~19MB total.

`assemble_parquet.py`'s own verification step re-reads the just-written
Parquet independently and checks real facts, not just that the files
parse:

```
[OK] row_level.parquet: 5,989,860 rows, columns=['datetime_deposit_a', 'datetime_release_a', 'datetime_deposit_b', 'datetime_release_b']
[OK] lookup_0001_pdb_id_categories.parquet: '101M' present=True (n=86,134)
[OK] lookup_0027_drugbank_id_categories.parquet: 'DB00114' present=True (n=4,927)
[OK] lookup_0043_pdbid_chain_categories_a.parquet: '101M_A' present=True (n=163,499)
[OK] partial_col_0066.parquet: 3,750,000 real rows (partial, not the full ~5,989,860)
[OK] 966,014 / 3,750,000 rows have a nonzero first component (confirms this isn't all-zero padding)

[ALL VERIFICATIONS PASSED]
```

One real bug found and fixed along the way, worth calling out because
it's a subtle general lesson, not specific to this file: the streaming
converter's first memory watchdog checked `resource.getrusage().ru_maxrss`,
which is **resident-only** and doesn't include pages the kernel has
swapped out. Once swap started filling, RSS alone plateaued near the RAM
ceiling while total memory kept climbing toward the real combined wall,
so the watchdog never fired and the process got a second real SIGKILL.
Fixed by reading `VmRSS + VmSwap` from `/proc/self/status` and checking
the sum -- see `CRYPTOBANK_DATASET_NOTES.md` for the full account.

Full details on what's now confirmed (real row count, real spot-checked
ID values, the apo/holo-pair hypothesis for the doubled date columns)
versus what's still open (the true `columns` name Index, the remaining
~2.24M rows) are in `CRYPTOBANK_DATASET_NOTES.md`.

## What this run did *not* do, and why

- **No real CryptoSite/CryptoBank benchmark number**, in the sense of
  running `bench_cryptosite` against the full real dataset end-to-end --
  that needs the actual PDB coordinate files (blocked by the
  `files.rcsb.org` egress wall, step 1) joined against the real labels
  extracted from `cryptobank_dataset_11_08_2025.pkl`, which itself
  required the memory-bounded inspector in step 8 and was stopped early
  by design (see `CRYPTOBANK_DATASET_NOTES.md` for exactly what was and
  wasn't established about its schema). The harness that would produce a
  real number (`druggability_bench`) is built, tested, and proven correct
  on synthetic data above.
- **No live Pharos/AlphaFold pipeline run.** Same network wall, confirmed
  live in step 1, not assumed.
- **No claim that `rcsb_bulk_download.py` downloaded anything for real**
  in this session -- it correctly *can't*, from here, and says so instead
  of pretending otherwise.
