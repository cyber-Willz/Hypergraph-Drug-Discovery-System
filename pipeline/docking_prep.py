"""
Converts a target_druggability JSON report's top-ranked pocket into:

1. A pocket center (mean C-alpha coordinate of the pocket's residues) and
   a bounding-box size -- the input a pocket-aware docking tool needs to
   constrain its search instead of docking blind.
2. A flexible-residue list (chain, res_seq) for tools that model
   side-chain flexibility at the pocket (e.g. DiffDock-Pocket).
3. An OpenMM restraint-selection template: which residues to leave
   unrestrained (the pocket) versus positionally restrained (everything
   else) during an equilibration/refinement run before docking, plus the
   actual OpenMM snippet that applies it.

Honest limits on step 1: "pocket center" here is the geometric mean of
Cα positions, a coarse stand-in for a real binding-site centroid (which
would use a probe-accessible cavity detection, e.g. fpocket/CASTp, not
just residue Cα averaging). It is what non-backtracking centrality's
output can support without extra tooling; treat the resulting grid/box as
a starting point to inspect visually (PyMOL/ChimeraX), not a final
docking-ready site definition.

Docking-tool CSV schema note: vanilla DiffDock (gcorso/DiffDock) does
*blind* docking and has no pocket-center input at all -- the pocket
information here is wasted on it. DiffDock-Pocket
(github.com/plainerman/DiffDock-Pocket) is the pocket-aware variant this
is meant for, but its exact CSV column names for pocket center / flexible
residues aren't stable enough across versions to hardcode here reliably
-- check `data/protein_ligand_example_csv.csv` in your DiffDock-Pocket
checkout and map the fields below onto it, rather than trusting an
assumed column name silently.
"""
from __future__ import annotations

import dataclasses
import json
from pathlib import Path
from typing import Optional

from pdb_geometry import parse_ca_coords


@dataclasses.dataclass
class PocketGeometry:
    structure_id: str
    rank: int
    center: tuple[float, float, float]
    box_size: tuple[float, float, float]  # padded bounding box, Angstroms
    residues: list[dict]  # [{chain, res_seq, res_name, score}]


def compute_pocket_geometry(
    report_entry: dict, pdb_path: Path, pocket_rank: int = 1, box_padding: float = 6.0,
) -> Optional[PocketGeometry]:
    """`report_entry` is one element of target_druggability's JSON output
    (one structure's TargetReport). Returns None if that structure has no
    pocket at the requested rank."""
    pockets = report_entry.get("pockets", [])
    pocket = next((p for p in pockets if p["rank"] == pocket_rank), None)
    if pocket is None:
        return None

    ca = parse_ca_coords(pdb_path)
    coords = []
    missing = []
    for r in pocket["residues"]:
        key = (r["chain"], r["res_seq"])
        if key in ca:
            coords.append(ca[key])
        else:
            missing.append(key)
    if missing:
        raise ValueError(
            f"{len(missing)} pocket residue(s) not found in {pdb_path} by (chain, res_seq): {missing[:5]}..."
            " -- is this the same structure file the report was generated from?"
        )
    if not coords:
        return None

    n = len(coords)
    cx = sum(c[0] for c in coords) / n
    cy = sum(c[1] for c in coords) / n
    cz = sum(c[2] for c in coords) / n
    dx = (max(c[0] for c in coords) - min(c[0] for c in coords)) + 2 * box_padding
    dy = (max(c[1] for c in coords) - min(c[1] for c in coords)) + 2 * box_padding
    dz = (max(c[2] for c in coords) - min(c[2] for c in coords)) + 2 * box_padding

    return PocketGeometry(
        structure_id=report_entry["structure_id"],
        rank=pocket_rank,
        center=(cx, cy, cz),
        box_size=(dx, dy, dz),
        residues=pocket["residues"],
    )


def write_pocket_grid_json(geom: PocketGeometry, out_path: Path) -> None:
    out_path.write_text(json.dumps(dataclasses.asdict(geom), indent=2))


OPENMM_RESTRAINT_TEMPLATE = '''\
"""
Auto-generated OpenMM pocket-refinement restraint template for
{structure_id}, pocket rank {rank}.

Strategy: restrain the whole protein to its input coordinates EXCEPT the
candidate pocket residues (and a small penumbra you may want to widen),
then run a short NVT/NPT equilibration so the pocket region can relax and
partially open before docking -- this is the standard "let the cryptic
pocket breathe" prep step; non-backtracking centrality only told you
*where* to look, it did not simulate pocket opening. Combine this with
DiffDock-Pocket / AutoDock or with an unrestrained follow-up production
run if you want actual pocket-opening sampling (e.g. mixed-solvent MD,
adaptive sampling) rather than just a relaxed static snapshot.

Fill in FORCEFIELD/WATER MODEL choices for your system; this template
assumes a fixed protein-only or already-solvated PDB and existing OpenMM
+ PDBFixer installation (`pip install openmm pdbfixer`).
"""
from openmm.app import *
from openmm import *
from openmm.unit import *

PDB_PATH = "{pdb_path}"
PROTEIN_RESTRAINT_K = 5.0 * kilocalories_per_mole / angstrom**2

# (chain, res_seq) pairs to leave UNRESTRAINED -- the candidate pocket.
POCKET_RESIDUES = {pocket_residues!r}

pdb = PDBFile(PDB_PATH)
forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff, constraints=HBonds)

restraint = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
restraint.addGlobalParameter("k", PROTEIN_RESTRAINT_K)
restraint.addPerParticleParameter("x0")
restraint.addPerParticleParameter("y0")
restraint.addPerParticleParameter("z0")

pocket_set = set(POCKET_RESIDUES)
for atom in pdb.topology.atoms():
    chain_id = atom.residue.chain.id
    res_seq = int(atom.residue.id)
    if (chain_id, res_seq) in pocket_set:
        continue  # leave pocket residues free to move
    if atom.element is not None and atom.element.symbol == "H":
        continue  # don't restrain hydrogens
    pos = pdb.positions[atom.index].value_in_unit(nanometers)
    restraint.addParticle(atom.index, list(pos))
system.addForce(restraint)

integrator = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 0.002 * picoseconds)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()
simulation.step(50000)  # 100 ps NVT relaxation of the pocket region; extend as needed

simulation.reporters.append(PDBReporter("{structure_id}_pocket_relaxed.pdb", 50000))
simulation.step(1)
'''


def write_openmm_template(geom: PocketGeometry, pdb_path: Path, out_path: Path) -> None:
    pocket_residues = [(r["chain"], r["res_seq"]) for r in geom.residues]
    out_path.write_text(OPENMM_RESTRAINT_TEMPLATE.format(
        structure_id=geom.structure_id, rank=geom.rank, pdb_path=str(pdb_path),
        pocket_residues=pocket_residues,
    ))


def write_diffdock_pocket_row(geom: PocketGeometry, protein_path: Path, ligand_description: str, out_path: Path) -> None:
    """Appends (or creates) a batch CSV row with the fields DiffDock-Pocket
    needs conceptually -- complex name, protein path, ligand, pocket
    center. Verify the exact header your DiffDock-Pocket checkout expects
    (see module docstring) before relying on this being drop-in."""
    header = "complex_name,protein_path,ligand_description,pocket_center_x,pocket_center_y,pocket_center_z\n"
    row = (
        f"{geom.structure_id}_pocket{geom.rank},{protein_path},{ligand_description},"
        f"{geom.center[0]:.3f},{geom.center[1]:.3f},{geom.center[2]:.3f}\n"
    )
    write_header = not out_path.exists()
    with open(out_path, "a") as f:
        if write_header:
            f.write(header)
        f.write(row)
