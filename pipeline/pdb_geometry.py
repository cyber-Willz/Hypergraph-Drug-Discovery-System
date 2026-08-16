"""
Minimal C-alpha coordinate extraction from a PDB file, using the exact
same fixed-width column layout as ../target_druggability/src/pdb.rs
(atom name 13-16, resName 18-20, chainID 22, resSeq 23-26, coords 31-54)
so `(chain, res_seq)` keys line up 1:1 with what the Rust report's
residue list uses. This exists because target_druggability's JSON report
carries residue identity (chain/res_seq/res_name) and score, not
coordinates -- docking prep needs the coordinates back to compute a
pocket center.
"""
from __future__ import annotations

from pathlib import Path


def parse_ca_coords(pdb_path: Path) -> dict[tuple[str, int], tuple[float, float, float]]:
    coords: dict[tuple[str, int], tuple[float, float, float]] = {}
    text = Path(pdb_path).read_text(errors="replace")
    for line in text.splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        if len(line) < 54:
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        chain = line[21:22].strip() or " "
        try:
            res_seq = int(line[22:26].strip())
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        key = (chain, res_seq)
        if key not in coords:  # first altloc wins, matching the Rust parser
            coords[key] = (x, y, z)
    return coords
