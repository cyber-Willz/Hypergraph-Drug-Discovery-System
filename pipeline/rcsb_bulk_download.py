"""
Bulk-download real PDB coordinate files for a list of structure IDs, keyed
on the target_druggability/druggability_bench naming convention used
throughout this repo:

  - Ordinary RCSB IDs (e.g. "1a8o") -> fetched from
    https://files.rcsb.org/download/{ID}.pdb -- this is RCSB's documented
    File Download Service endpoint (files.rcsb.org).
  - AlphaFold-derived apo entries, conventionally written
    "AF-<UNIPROT>-F1-model_vN" -> NOT handled here. Route these through
    pipeline/alphafold_client.py instead, which already implements the
    correct AlphaFold DB endpoint (alphafold.ebi.ac.uk) and download path
    -- duplicating that logic here would just be a second place for the
    same bug to hide.

CASE MATTERS. `druggability_bench` derives each structure's `structure_id`
from the PDB file's stem, byte-for-byte, with no case-folding (see
`druggability_bench/src/pdb.rs::parse_ca_atoms_file` and the join logic in
`druggability_bench/src/labels.rs`). A live inspection of the real
CryptoBank dataset this session (see `CRYPTOBANK_DATASET_NOTES.md`)
confirmed its PDB ids are UPPERCASE (`'101M', '10GS', '11AS', ...`, and
composite `{PDBID}_{CHAIN}` ids like `'101M_A'`) -- the opposite case
convention from `druggability_bench/README.md`'s lowercase illustrative
examples (`'1jwp'`, `'1a8o'`). This script therefore writes each file
using EXACTLY the case you passed in (`--ids 101M` -> `101M.pdb`;
`--ids 1a8o` -> `1a8o.pdb`) and does not uppercase/lowercase anything --
match whatever case your `--labels` CSV actually uses, or
`bench_cryptosite` will silently skip every structure with "no rows in
--labels for structure_id ...". (RCSB's download endpoint itself is
case-insensitive on the ID segment of the URL, so this choice costs
nothing on the fetch side.)

Usage:
    python rcsb_bulk_download.py --ids 1a8o,1jwp,4hhb --out-dir pdb_cache
    python rcsb_bulk_download.py --ids-file structure_ids.txt --out-dir pdb_cache

`--ids-file` expects one PDB ID per line (blank lines and lines starting
with '#' are skipped). IDs starting with "AF-" are detected and skipped
with a pointer to alphafold_client.py rather than silently mishandled.

Network requirement: files.rcsb.org. Confirmed NOT in this sandbox's
egress allowlist by a live request this session -- the proxy returns
`403` with header `x-deny-reason: host_not_allowed` and a body telling
you to add the host to your network egress settings, rather than a
connection failure. That distinction matters: it's a policy wall, not a
flaky network, so this script detects it and fails fast on the first
attempt instead of burning retries with backoff on something that will
never succeed from here. Run this on a machine/session with real egress
to files.rcsb.org (e.g. via your network settings, or a machine outside
this sandbox).
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import requests

DOWNLOAD_BASE = "https://files.rcsb.org/download"


class SandboxDenied(Exception):
    """Raised when the egress proxy itself refuses the host -- a policy
    wall, not a transient network problem. Never worth retrying."""


def _check_sandbox_denial(resp: requests.Response) -> None:
    if resp.status_code == 403 and resp.headers.get("x-deny-reason") == "host_not_allowed":
        raise SandboxDenied(resp.text.strip() or "host not in network egress allowlist")


def fetch_one(
    pdb_id: str, out_dir: Path, session: requests.Session, timeout: float = 30.0, retries: int = 3,
) -> tuple[str, bool, str]:
    """Returns (pdb_id, success, message). `pdb_id`'s case is preserved
    exactly as given -- see the case-matters note in the module docstring."""
    pdb_id = pdb_id.strip()
    if not pdb_id:
        return pdb_id, False, "empty id"
    if pdb_id.upper().startswith("AF-"):
        return pdb_id, False, "AlphaFold entry -- use pipeline/alphafold_client.py, not this script"

    out_path = out_dir / f"{pdb_id}.pdb"
    if out_path.exists() and out_path.stat().st_size > 0:
        return pdb_id, True, "already cached"

    url = f"{DOWNLOAD_BASE}/{pdb_id}.pdb"
    last_err = None
    for attempt in range(1, retries + 1):
        try:
            resp = session.get(url, timeout=timeout)
            _check_sandbox_denial(resp)
            if resp.status_code == 404:
                return pdb_id, False, f"not found at {url} (obsolete/superseded entry, or CIF-only structure)"
            resp.raise_for_status()
            out_path.write_bytes(resp.content)
            return pdb_id, True, f"downloaded {len(resp.content):,} bytes"
        except SandboxDenied as e:
            # Deterministic policy denial from our own egress proxy, not a
            # flaky remote -- retrying it is pointless, so don't.
            return pdb_id, False, f"blocked by sandbox network egress policy: {e}"
        except requests.exceptions.RequestException as e:
            last_err = e
            if attempt < retries:
                time.sleep(1.5 * attempt)  # simple backoff, be polite to a shared public API
    return pdb_id, False, f"could not reach RCSB ({url}) after {retries} attempts: {last_err}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--ids", help="comma-separated PDB IDs")
    src.add_argument("--ids-file", type=Path, help="file with one PDB ID per line")
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--rate-limit-sleep", type=float, default=0.2,
                     help="seconds to sleep between requests (be a good API citizen)")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    if args.ids:
        ids = [x.strip() for x in args.ids.split(",") if x.strip()]
    else:
        ids = [
            line.strip() for line in args.ids_file.read_text().splitlines()
            if line.strip() and not line.strip().startswith("#")
        ]

    print(f"fetching {len(ids)} structure(s) -> {args.out_dir}", file=sys.stderr)

    ok, failed, af_skipped, denied = 0, 0, 0, 0
    with requests.Session() as session:
        for i, pdb_id in enumerate(ids, 1):
            sid, success, msg = fetch_one(pdb_id, args.out_dir, session)
            status = "OK" if success else "FAIL"
            print(f"[{i}/{len(ids)}] {sid}: {status} -- {msg}", file=sys.stderr)
            if success:
                ok += 1
            elif "AlphaFold" in msg:
                af_skipped += 1
            elif "sandbox network egress policy" in msg:
                denied += 1
            else:
                failed += 1
            time.sleep(args.rate_limit_sleep)

    print(
        f"\ndone: {ok} downloaded/cached, {failed} failed, {af_skipped} AlphaFold entries skipped "
        f"(route those through pipeline/alphafold_client.py), {denied} blocked by sandbox network policy "
        f"(add files.rcsb.org to your network egress settings and re-run)",
        file=sys.stderr,
    )
    if failed or denied:
        sys.exit(1)


if __name__ == "__main__":
    main()
