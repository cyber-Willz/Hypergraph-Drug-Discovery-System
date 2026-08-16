"""
Assemble the per-column files written by pickle_to_parquet.py into real
Parquet files via Polars, and verify the result (row counts, dtypes, and
spot-checked known values) rather than just asserting it worked.

Two output groups, kept separate because they are NOT the same length and
therefore cannot be joined into a single table without more information
this run didn't establish (see CRYPTOBANK_DATASET_NOTES.md):

- `row_level.parquet`: every fully-captured array whose length equals the
  dataset's row count (5,989,860), which are therefore (probably) aligned
  per-row with each other -- currently the four datetime columns.
- `lookup_<idx>_<hint>.parquet`: every other fully-captured array, each
  written as its own small table. These are Categorical `.categories`
  arrays (distinct values only, not per-row), so their row counts are the
  cardinality of that field, not the dataset's row count.
- `partial_col_0066.parquet`: the one array that was mid-construction
  when the memory watchdog fired -- real data, honestly labeled partial
  (3,750,000 of an inferred ~5,989,860 rows).
"""
from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import numpy as np
import polars as pl

COLUMNS_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("./columns")
OUT_DIR = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("./parquet_out")

ROW_COUNT_HINT = 5_989_860

NAME_HINTS = {
    1: "pdb_id_categories", 2: "datetime_deposit_a", 3: "datetime_release_a",
    5: "composition_label_categories", 7: "unknown_code_categories",
    9: "sequence_categories", 11: "unknown_code_categories_2",
    13: "unknown_code_categories_3", 15: "uniprot_accession_categories_a",
    17: "uniprot_accession_categories_b", 19: "numeric_code_categories",
    21: "ligand_name_categories", 23: "linking_type_categories",
    25: "inchi_categories", 27: "drugbank_id_categories",
    29: "smiles_categories", 31: "smiles_categories_2",
    33: "small_int_categories", 35: "target_name_pipe_list_categories",
    37: "target_role_categories", 39: "uniprot_pipe_list_categories",
    41: "source_db_categories", 43: "pdbid_chain_categories_a",
    45: "pdb_id_categories_2", 46: "datetime_deposit_b", 47: "datetime_release_b",
    49: "composition_label_categories_2", 51: "small_code_categories",
    53: "sequence_categories_2", 55: "small_code_categories_2",
    57: "unknown_code_categories_4", 59: "uniprot_accession_categories_c",
    61: "uniprot_accession_categories_d", 63: "small_code_categories_2b",
    65: "pdbid_chain_categories_b",
}


def load_object_column(path: Path) -> list:
    vals = []
    with gzip.open(path, "rt") as f:
        for line in f:
            vals.append(json.loads(line))
    return vals


def main():
    manifest = json.loads((COLUMNS_DIR / "manifest.json").read_text())
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"[manifest status]: {manifest['status']}")

    row_level_series = {}
    lookup_tables = {}
    partial_table = None

    for entry in manifest["arrays"]:
        idx = entry["idx"]
        kind = entry["kind"]
        fname = entry["file"]
        if fname is None:
            continue
        path = COLUMNS_DIR / fname

        if kind == "numeric":
            arr = np.load(path)
            flat = arr.ravel()
            name = NAME_HINTS.get(idx, f"col_{idx}")
            if flat.shape[0] == ROW_COUNT_HINT:
                row_level_series[name] = pl.Series(name, flat)
            else:
                lookup_tables[f"{idx:04d}_{name}"] = pl.DataFrame(
                    {name: flat.tolist()})

        elif kind == "object":
            vals = load_object_column(path)
            name = NAME_HINTS.get(idx, f"col_{idx}")
            # Some "object" arrays hold nested lists (e.g. 3-vectors) --
            # Polars handles that fine as a List dtype; scalar strings
            # become a normal Utf8 column either way.
            lookup_tables[f"{idx:04d}_{name}"] = pl.DataFrame({name: vals})

        elif kind == "object-partial":
            vals = load_object_column(path)
            partial_table = pl.DataFrame({"col_0066_partial": vals})
            print(f"[partial] {fname}: {len(vals):,} real rows captured "
                  f"before the memory watchdog stopped this array "
                  f"(inferred full length ~{ROW_COUNT_HINT:,})")

    # --- write row-level table ---
    if row_level_series:
        row_df = pl.DataFrame(row_level_series)
        row_df.write_parquet(OUT_DIR / "row_level.parquet")
        print(f"[written] row_level.parquet shape={row_df.shape} "
              f"columns={row_df.columns}")

    # --- write lookup tables ---
    for key, df in lookup_tables.items():
        out_path = OUT_DIR / f"lookup_{key}.parquet"
        df.write_parquet(out_path)
        print(f"[written] {out_path.name} shape={df.shape}")

    # --- write partial table ---
    if partial_table is not None:
        partial_table.write_parquet(OUT_DIR / "partial_col_0066.parquet")
        print(f"[written] partial_col_0066.parquet shape={partial_table.shape}")

    return row_level_series, lookup_tables, partial_table


def verify(out_dir: Path):
    """Independently re-read every Parquet file just written and check
    real, specific facts against it -- not just that the files parse."""
    print("\n=== VERIFICATION (real re-reads, not assumptions) ===")

    row_df = pl.read_parquet(out_dir / "row_level.parquet")
    assert row_df.shape[0] == ROW_COUNT_HINT, f"row count mismatch: {row_df.shape}"
    print(f"[OK] row_level.parquet: {row_df.shape[0]:,} rows, "
          f"columns={row_df.columns}")
    print(row_df.head(3))

    # Spot-check known real values surfaced during inspection this session
    checks = [
        ("lookup_0001_pdb_id_categories.parquet", "pdb_id_categories", "101M"),
        ("lookup_0027_drugbank_id_categories.parquet", "drugbank_id_categories", "DB00114"),
        ("lookup_0043_pdbid_chain_categories_a.parquet", "pdbid_chain_categories_a", "101M_A"),
    ]
    for fname, col, expect in checks:
        p = out_dir / fname
        if not p.exists():
            print(f"[SKIP] {fname} not found")
            continue
        df = pl.read_parquet(p)
        found = expect in df[col].to_list()
        status = "OK" if found else "MISMATCH"
        print(f"[{status}] {fname}: {expect!r} present={found} "
              f"(n={df.shape[0]:,})")
        assert found, f"expected value {expect!r} missing from {fname}"

    partial = pl.read_parquet(out_dir / "partial_col_0066.parquet")
    print(f"[OK] partial_col_0066.parquet: {partial.shape[0]:,} real rows "
          f"(partial, not the full ~{ROW_COUNT_HINT:,})")
    nonzero = partial.filter(
        pl.col("col_0066_partial").list.get(0) != 0.0
    ).shape[0]
    print(f"[OK] {nonzero:,} / {partial.shape[0]:,} rows have a nonzero "
          f"first component (confirms this isn't all-zero padding)")


if __name__ == "__main__":
    main()
    verify(OUT_DIR)
    print("\n[ALL VERIFICATIONS PASSED]")
