"""
End-to-end orchestrator implementing the roadmap pipeline:

    Pharos GraphQL API -> AlphaFold PDB -> target_druggability -> OpenMM / DiffDock-Pocket prep

Two modes:

  Single target:
    python pipeline.py --symbol CA2 --out-dir runs/ca2

  Disease-driven discovery + batch triage:
    python pipeline.py --disease "pancreatic cancer" --top 20 --out-dir runs/panc \\
        --difficult-only

`--difficult-only` restricts to Pharos Tdark/Tbio targets -- the actual
IDG "understudied, worth illuminating" tier, which is the whole point of
cross-referencing TCRD/Pharos in the first place (see
target_druggability's own report ranking: it already prioritizes
Tdark/Tbio + strong pocket signal, this flag just avoids spending
AlphaFold/compute budget on well-studied Tclin/Tchem targets you'd filter
out later anyway).

Requires: `pip install requests`, a built `target_druggability` release
binary (built automatically via cargo if missing), and network access to
pharos-api.ncats.io / alphafold.ebi.ac.uk (neither of which is reachable
from a sandboxed environment without egress to those hosts -- run this on
a machine that has it).
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import alphafold_client
import docking_prep
import druggability_runner
import pharos_client


def process_one_target(
    target: pharos_client.PharosTarget,
    workspace_root: Path,
    out_dir: Path,
    binary_path: Path,
    cutoff: float,
    min_seq_sep: int,
    top_percentile: float,
) -> dict | None:
    """Run the full pipeline for one Pharos target. Returns a summary dict,
    or None if it couldn't be processed (no AlphaFold model, parse
    failure, etc -- logged to stderr, not fatal for batch runs)."""
    target_dir = out_dir / target.sym
    target_dir.mkdir(parents=True, exist_ok=True)

    print(f"[{target.sym}] fetching AlphaFold model for {target.uniprot}...", file=sys.stderr)
    pred = alphafold_client.get_prediction(target.uniprot)
    if pred is None:
        print(f"[{target.sym}] no AlphaFold model for {target.uniprot}, skipping", file=sys.stderr)
        return None
    pdb_path = alphafold_client.download_pdb(pred, target_dir)
    print(f"[{target.sym}] model {pred.entry_id} -> {pdb_path}", file=sys.stderr)

    tcrd_csv_path = target_dir / "tcrd_row.csv"
    tcrd_csv_path.write_text(pharos_client.TCRD_CSV_HEADER + "\n" + target.to_tcrd_csv_row() + "\n")

    config = druggability_runner.DruggabilityRunConfig(
        binary_path=binary_path, cutoff=cutoff, min_seq_sep=min_seq_sep, top_percentile=top_percentile,
    )
    report_path = target_dir / "report.json"
    print(f"[{target.sym}] running target_druggability...", file=sys.stderr)
    report = druggability_runner.run(
        config, pdb_path, report_path, tcrd_csv_path=tcrd_csv_path,
        symbol=target.sym, uniprot=target.uniprot,
    )
    entry = report[0]
    n_pockets = len(entry.get("pockets", []))
    print(f"[{target.sym}] {n_pockets} candidate pocket(s), rho_B={entry['global_coupling_rho_b']:.3f}", file=sys.stderr)

    docking_summary = None
    if n_pockets > 0:
        geom = docking_prep.compute_pocket_geometry(entry, pdb_path, pocket_rank=1)
        if geom is not None:
            docking_prep.write_pocket_grid_json(geom, target_dir / "pocket_grid.json")
            docking_prep.write_openmm_template(geom, pdb_path, target_dir / "openmm_pocket_prep.py")
            docking_prep.write_diffdock_pocket_row(
                geom, pdb_path, ligand_description="REPLACE_WITH_SMILES_OR_SDF_PATH",
                out_path=out_dir / "diffdock_pocket_batch.csv",
            )
            docking_summary = {"center": geom.center, "box_size": geom.box_size, "n_residues": len(geom.residues)}
            print(f"[{target.sym}] pocket center {geom.center} written to {target_dir/'pocket_grid.json'}", file=sys.stderr)

    return {
        "symbol": target.sym,
        "uniprot": target.uniprot,
        "tdl": target.tdl,
        "family": target.fam,
        "alphafold_entry": pred.entry_id,
        "n_pockets": n_pockets,
        "global_coupling_rho_b": entry["global_coupling_rho_b"],
        "top_pocket": docking_summary,
        "report_path": str(report_path),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--symbol", help="single-target mode: gene symbol")
    ap.add_argument("--uniprot", help="single-target mode: UniProt accession (used if --symbol not given)")
    ap.add_argument("--disease", help="batch mode: discover targets associated with this disease via Pharos")
    ap.add_argument("--top", type=int, default=25, help="batch mode: max targets to discover (default 25)")
    ap.add_argument("--difficult-only", action="store_true", help="batch mode: keep only Tdark/Tbio targets")
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--workspace-root", type=Path, default=Path(__file__).resolve().parent.parent,
                     help="path to the target_druggability cargo workspace (default: parent of this pipeline/ dir)")
    ap.add_argument("--cutoff", type=float, default=8.0)
    ap.add_argument("--min-seq-sep", type=int, default=3)
    ap.add_argument("--top-percentile", type=float, default=0.9)
    args = ap.parse_args()

    if not args.symbol and not args.uniprot and not args.disease:
        ap.error("provide --symbol/--uniprot (single target) or --disease (batch discovery)")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    binary_path = druggability_runner.build_binary(args.workspace_root)

    if args.disease:
        targets = pharos_client.discover_targets_for_disease(args.disease, top=args.top)
        if args.difficult_only:
            targets = [t for t in targets if t.is_difficult]
        print(f"discovered {len(targets)} target(s) for {args.disease!r}"
              f"{' (Tdark/Tbio only)' if args.difficult_only else ''}", file=sys.stderr)
    else:
        t = pharos_client.get_target(symbol=args.symbol, uniprot=args.uniprot)
        if t is None:
            print("no Pharos record found for that target", file=sys.stderr)
            sys.exit(1)
        targets = [t]

    results = []
    for t in targets:
        try:
            r = process_one_target(
                t, args.workspace_root, args.out_dir, binary_path,
                args.cutoff, args.min_seq_sep, args.top_percentile,
            )
            if r is not None:
                results.append(r)
        except Exception as e:  # noqa: BLE001 -- batch runs must not die on one bad target
            print(f"[{t.sym}] FAILED: {e}", file=sys.stderr)

    results.sort(key=lambda r: (r["tdl"] not in ("Tdark", "Tbio"), -(r["n_pockets"])))
    summary_path = args.out_dir / "summary.json"
    summary_path.write_text(json.dumps(results, indent=2))
    print(f"\n{len(results)} target(s) processed. Summary: {summary_path}", file=sys.stderr)
    if any(r["top_pocket"] for r in results):
        print(f"Batch DiffDock-Pocket rows (fill in ligand SMILES/SDF before running): "
              f"{args.out_dir / 'diffdock_pocket_batch.csv'}", file=sys.stderr)


if __name__ == "__main__":
    main()
