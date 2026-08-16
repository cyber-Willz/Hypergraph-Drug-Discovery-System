"""
Minimal client for the AlphaFold Protein Structure Database REST API.

All endpoints are keyed on UniProt accession, confirmed against the
AlphaFold DB paper (Varadi et al., Nucleic Acids Research) and current
API docs at https://alphafold.ebi.ac.uk/api-docs:

    GET https://alphafold.ebi.ac.uk/api/prediction/{uniprot_accession}
    -> JSON list (usually one entry) with entryId, pdbUrl, cifUrl,
       paeImageUrl, paeDocUrl, latestVersion, ...

Confidence caveat this module deliberately does NOT paper over: AlphaFold
models are predictions, not experimental structures. A model's per-residue
pLDDT (fetchable via the PAE/confidence files, not exposed directly by this
thin client) matters a lot for how much to trust a predicted cryptic pocket
in a low-confidence region -- low-pLDDT regions are often just disordered,
and non-backtracking centrality has no way to know that. If you're
triaging Tdark targets that only have AlphaFold models (no experimental
structure), treat pocket calls in low-confidence regions with real
skepticism, not as equivalent evidence to a pocket found in a
crystallographic apo structure.
"""
from __future__ import annotations

import dataclasses
from pathlib import Path
from typing import Optional

import requests

API_BASE = "https://alphafold.ebi.ac.uk/api/prediction"


@dataclasses.dataclass
class AlphaFoldPrediction:
    uniprot: str
    entry_id: str
    pdb_url: str
    cif_url: str
    model_created_date: Optional[str] = None
    latest_version: Optional[int] = None


def get_prediction(uniprot: str, timeout: float = 30.0) -> Optional[AlphaFoldPrediction]:
    """Fetch prediction metadata for a UniProt accession. Returns None if
    AlphaFold DB has no model (common for very short peptides, some
    multi-pass membrane proteins historically, or accessions AlphaFold
    hasn't covered yet -- check https://alphafold.ebi.ac.uk/ manually
    before assuming this is a bug)."""
    resp = requests.get(f"{API_BASE}/{uniprot}", timeout=timeout)
    if resp.status_code == 404:
        return None
    resp.raise_for_status()
    entries = resp.json()
    if not entries:
        return None
    e = entries[0]  # AlphaFold DB returns the latest model version first
    return AlphaFoldPrediction(
        uniprot=uniprot,
        entry_id=e["entryId"],
        pdb_url=e["pdbUrl"],
        cif_url=e.get("cifUrl", ""),
        model_created_date=e.get("modelCreatedDate"),
        latest_version=e.get("latestVersion"),
    )


def download_pdb(prediction: AlphaFoldPrediction, out_dir: Path, timeout: float = 60.0) -> Path:
    """Download the model's PDB coordinate file to `out_dir/<uniprot>.pdb`."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{prediction.uniprot}.pdb"
    resp = requests.get(prediction.pdb_url, timeout=timeout)
    resp.raise_for_status()
    out_path.write_bytes(resp.content)
    return out_path


def fetch_structure(uniprot: str, out_dir: Path) -> Optional[Path]:
    """Convenience: metadata fetch + download in one call. Returns the
    local PDB path, or None if no AlphaFold model exists for this
    accession."""
    pred = get_prediction(uniprot)
    if pred is None:
        return None
    return download_pdb(pred, out_dir)


if __name__ == "__main__":
    import argparse
    import sys

    ap = argparse.ArgumentParser(description="Download an AlphaFold model PDB by UniProt accession.")
    ap.add_argument("uniprot")
    ap.add_argument("--out-dir", default="./structures")
    args = ap.parse_args()

    path = fetch_structure(args.uniprot, Path(args.out_dir))
    if path is None:
        print(f"no AlphaFold model found for {args.uniprot}", file=sys.stderr)
        sys.exit(1)
    print(path)
