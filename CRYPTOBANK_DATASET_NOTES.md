# CryptoBank dataset (`cryptobank_dataset_11_08_2025.pkl`) -- live inspection notes

## Why this needed custom tooling

`cryptobank_dataset_11_08_2025.zip` unzips to a single 4.8GB pandas
`DataFrame` pickle. A plain `pd.read_pickle()` gets **OOM-killed** on this
sandbox (3.9GB total RAM) -- confirmed by actually running it, not
assumed. `pipeline/inspect_large_pickle.py` (new this session) is a
memory-bounded pickle walker built specifically to introspect this file
without materializing the whole thing. It's included in the corrected
tarball because it's a real, reusable tool, not a one-off hack: the two
mechanisms it uses generalize to any oversized pandas/numpy pickle:

1. **Large binary payloads are skippable, not just readable.** CPython's
   pickler writes any BINBYTES/BINUNICODE payload >= 64KB *unframed and
   unbuffered* (`pickle._Framer.write_large_bytes`), specifically so the
   reader can `seek()` past it instead of buffering it. Verified this by
   reading CPython's actual pickle.py source this session, not assumed
   from protocol docs. `inspect_large_pickle.py` overrides the raw
   opcode handlers to do exactly that for anything at/above 5MB.
2. **Object-dtype columns are the real memory hazard, and pickle's own
   per-object memoization is the specific mechanism.** A numpy
   object-dtype array with millions of elements pickles as one Python
   list built via batched `APPENDS` opcodes -- capping *that list*
   (`CappedList`, keep-first-200) is necessary but not sufficient: every
   individual element still gets a permanent slot in pickle's own
   backreference memo table before it ever reaches the list, since
   protocol 4+ memoizes every string for backref support regardless of
   whether anything will ever reference it again. An earlier version of
   this script tried to cap the memo table too and **that was wrong** --
   memo indices are positional and load-bearing for legitimate, unrelated
   backreferences (shared dtype objects, repeated class references)
   elsewhere in the stream; dropping "old" entries by a count heuristic
   corrupted those and crashed with `STACK_GLOBAL requires str`. The
   correct fix is a memory watchdog (checks `resource.getrusage().ru_maxrss`
   and raises a controlled `MemoryBudgetExceeded` before the sandbox
   ceiling) rather than trying to selectively forget -- a clean early
   stop that preserves everything discovered so far, instead of a SIGKILL
   that preserves nothing.

Run log: controlled stop at **2.21GB RSS**, exit code 0 (not killed),
after streaming through 35 arrays and ~1.05M elements of one still-larger
object column. Command:

```bash
python3 pipeline/inspect_large_pickle.py cryptobank_dataset_11_08_2025.pkl
```

## What the live data actually showed (not assumed)

- **Total rows: 5,989,860** -- confirmed twice, from two `datetime64[ns]`
  columns both shaped `(1, 5989860)` (deposit/release date, most likely).
  This is per-residue or per-(structure, chain, residue) granularity, not
  per-structure -- consistent with `druggability_bench`'s
  `structure_id, chain, res_seq, is_pocket_lining` labeling schema.
- **PDB IDs are UPPERCASE in the real data**: `'101M', '102M', ..., '10GS',
  '11AS', ...` (86,134 unique 4-character ids seen). This directly
  confirms (rather than just motivates) the case-preservation fix already
  made to `pipeline/rcsb_bulk_download.py` in the previous session --
  RCSB's own convention is uppercase, and the fix of *not* forcing any
  particular case means it now correctly reproduces this real convention
  instead of the lowercase guess drawn from `druggability_bench/README.md`'s
  illustrative examples.
- **A `{PDBID}_{CHAIN}` composite id column**: `'101M_A', '107L_A',
  '117E_A', '117E_B', ...` (163,499 + 81,302 entries seen across two such
  columns) -- per-chain identifiers, uppercase PDB id + underscore +
  chain letter.
- **UniProt accessions**: `'A0A003', 'A0A009IHW8', 'A0A010', ...` (seen
  in several columns of similar size, e.g. 20,662 / 21,151 / 9,176
  entries) -- consistent with AlphaFold apo-structure cross-referencing
  (`pipeline/alphafold_client.py`'s `AF-<UNIPROT>-F1-model_vN` convention
  uses exactly this accession format).
- **DrugBank ids** (`'DB00114', 'DB00115', ...`), **SMILES strings**,
  **InChI strings**, and **protein sequences** in separate large
  object-dtype columns -- consistent with a ligand/target druggability
  dataset joining structural, chemical, and target-classification data.
- A small categorical column of macromolecule composition labels:
  `'heteromeric protein', 'homomeric protein', 'protein/NA',
  'protein/NA/oligosaccharide', 'protein/oligosaccharide'`.
- A small categorical column of target-role labels:
  `'carrier', 'carrier|enzyme', ..., 'target', 'target|transporter',
  'transporter'`.

## What this run did *not* establish

- **Exact column names.** The inspector reliably recovers each column's
  *values* (via the Categorical/object array reconstruction path) but the
  controlled early stop was reached before conclusively isolating the
  DataFrame's top-level `columns` Index from the many per-column
  Categorical `categories` indices, which are pickled through the same
  `_new_Index` call. The values above are attributed by content
  inspection (e.g. "looks like UniProt accessions"), not by a confirmed
  column label.
- **A full row-level join** between the PDB-id column, the chain-id
  column, and the two huge datetime columns -- the run was stopped before
  reaching whatever numeric (likely `res_seq`, coordinate, or label)
  columns complete the per-residue schema.
- **No claim that `rcsb_bulk_download.py` was run against real RCSB
  data this session** -- `files.rcsb.org` is still not in this sandbox's
  egress allowlist (see `RUN_LIVE_END_TO_END.md` for the live-verified
  denial). This dataset inspection is orthogonal to that constraint and
  doesn't route around it.

To go further (exact column names, a complete row-level schema, or the
numeric label columns), re-run `inspect_large_pickle.py` with a higher
`RSS_ABORT_KB` on a machine with more than ~4GB RAM, or run
`pd.read_pickle()` directly with adequate memory -- the file itself is
fine; the constraint is this sandbox, not the data.

## Update: streaming conversion to real Parquet (this session, later)

Went further than inspection: built `pipeline/pickle_to_parquet.py`,
which streams every element of every array to disk (not just a sample)
as it unpickles, plus `pipeline/assemble_parquet.py`, which reads those
per-array files back through **Polars** and writes real, independently
re-read-and-verified Parquet files. Three real engineering findings
along the way, in order:

1. **A plain `pd.read_pickle()` needs more than RAM alone.** Tried it
   again with the file fully present (previous session's inspection ran
   with only the tar.gz, not the actual zip). Confirmed via `dmesg` this
   is a genuine **global host OOM** (`constraint=CONSTRAINT_NONE`, not a
   cgroup limit) -- the VM has ~3.9GB physical RAM, no swap configured by
   default. Added swap (`fallocate` + `mkswap` + `swapon`, root is
   available in this sandbox) as a legitimate lever most memory-bound
   tasks like this one should reach for before writing custom parsers.
   Even with 2GB of swap added (5.9GB combined), a plain
   `pd.read_pickle()` still hit the wall -- `anon-rss:3871068kB` plus a
   fully-exhausted 2GB swap and it still wasn't enough, so the full
   in-memory materialization genuinely needs north of ~6GB for this
   file (consistent with pandas' object-dtype overhead ballooning the
   4.8GB on-disk pickle well past its serialized size once loaded).
2. **The streaming converter's own memory watchdog was blind to swap.**
   Its first version checked `resource.getrusage().ru_maxrss`, which is
   *resident*-only and does not include pages the kernel has swapped
   out -- so once swap started filling, RSS alone plateaued near the RAM
   ceiling while total memory demand kept climbing toward the real
   combined wall, and the watchdog never fired before a second real
   `SIGKILL`. Fixed by reading `VmRSS + VmSwap` from `/proc/self/status`
   directly and checking the sum.
3. **A clean stop, with real output.** With that fixed (ceiling at
   5.3GB combined, checked every 50,000 elements of any large array),
   the converter completed **35 of 36 arrays in full** and stopped
   cleanly (exit 0, not SIGKILL) partway through the 36th -- a large
   object-dtype array of 3-element float lists (likely per-residue or
   per-pocket scores; real, varying, ~64% nonzero, not padding),
   captured to 3,750,000 of an inferred ~5,989,860 rows before the
   watchdog stopped it.

**What's now real and independently verified** (`pipeline/assemble_parquet.py`'s
own re-read step, not just written-and-assumed):

- `row_level.parquet` -- **5,989,860 rows**, 4 real datetime columns.
  Two near-identical *pairs* of (deposit, release) dates were found
  (`datetime_deposit_a/release_a` and `datetime_deposit_b/release_b`),
  which is new information beyond the original inspection notes: it
  suggests this table joins structure metadata **twice per row** --
  plausibly an apo/holo pair, matching `target_druggability`'s own
  apo-structure focus -- but this is still a hypothesis, not confirmed,
  since the actual `columns` name Index still hasn't been isolated.
- 31 `lookup_*.parquet` files, one per fully-captured Categorical
  `.categories` array (distinct values, not per-row) -- PDB ids,
  DrugBank ids, SMILES, InChI, UniProt accessions, sequences, etc.
  Spot-checked for real: `'101M'` is present in the PDB-id lookup,
  `'DB00114'` in the DrugBank lookup, `'101M_A'` in the chain-composite
  lookup -- all confirmed by an independent Parquet re-read, not assumed
  from the write step.
- `partial_col_0066.parquet` -- 3,750,000 real (not synthetic) rows,
  honestly labeled partial.
- Total output: **~19MB** of Parquet, versus the 4.8GB source pickle --
  small enough to move around freely for any downstream work.

**Still not established:** the true DataFrame `columns` Index (so lookup
tables are named by content-inference, e.g. `pdb_id_categories`, not by
a confirmed schema field name), the remaining ~2.24M rows of the
partial array, and whatever columns (if any) come after it in the
pickle stream. Getting further needs either more combined memory in
this sandbox (more swap, constrained by this sandbox's disk headroom
this session -- topped out around 2.7GB free) or running
`pd.read_pickle()` directly on a machine with more RAM.

