"""
Thin wrapper around the compiled `target_druggability` CLI binary (see
../target_druggability). This deliberately shells out rather than
reimplementing the spectral analysis in Python: the Rust crate is the
validated implementation (see ../druggability_bench for the benchmark
harness), and re-deriving it here would just create a second
implementation to keep in sync.
"""
from __future__ import annotations

import dataclasses
import json
import subprocess
from pathlib import Path
from typing import Optional


@dataclasses.dataclass
class DruggabilityRunConfig:
    binary_path: Path
    cutoff: float = 8.0
    min_seq_sep: int = 3
    top_percentile: float = 0.9


def run(
    config: DruggabilityRunConfig,
    pdb_path: Path,
    out_json_path: Path,
    tcrd_csv_path: Optional[Path] = None,
    symbol: Optional[str] = None,
    uniprot: Optional[str] = None,
) -> list[dict]:
    """Run target_druggability on one structure and return the parsed
    report (a list with one entry per --pdb given; always length 1 here
    since we call it per-structure so batch runs can fail independently)."""
    cmd = [
        str(config.binary_path),
        "--pdb", str(pdb_path),
        "--cutoff", str(config.cutoff),
        "--min-seq-sep", str(config.min_seq_sep),
        "--top-percentile", str(config.top_percentile),
        "--out", str(out_json_path),
    ]
    if tcrd_csv_path is not None:
        cmd += ["--tcrd", str(tcrd_csv_path)]
    if symbol is not None:
        cmd += ["--symbol", symbol]
    if uniprot is not None:
        cmd += ["--uniprot", uniprot]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"target_druggability failed on {pdb_path} (exit {result.returncode}):\n{result.stderr}"
        )
    return json.loads(out_json_path.read_text())


def build_binary(workspace_root: Path) -> Path:
    """Build the release binary if it isn't already built. Returns its path."""
    binary = workspace_root / "target" / "release" / "target_druggability"
    if binary.exists():
        return binary
    result = subprocess.run(
        ["cargo", "build", "--release", "-p", "target_druggability"],
        cwd=str(workspace_root), capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"cargo build failed:\n{result.stderr}")
    if not binary.exists():
        raise RuntimeError(f"cargo build succeeded but {binary} is still missing")
    return binary
