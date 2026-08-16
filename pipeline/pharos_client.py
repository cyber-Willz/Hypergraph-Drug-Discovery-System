"""
Minimal client for the Pharos (NIH Illuminating the Druggable Genome / IDG)
GraphQL API.

Endpoint and schema confirmed against Pharos' own published example queries
at https://pharos.nih.gov/api (GraphQL Playground) as of this writing:

    query targetDetails{
      target(q:{sym:"ACE2"}) { name tdl fam sym description novelty }
    }
    query associatedTargets{
      targets(filter: { associatedDisease: "asthma" }) {
        targets(top: 5) { name sym diseaseAssociationDetails { name dataType evidence } }
      }
    }

Pharos' schema does evolve, and the interactive playground at
https://pharos.nih.gov/api is the authoritative source if a query here
starts failing -- this module fails loudly (raises) rather than silently
returning partial/wrong data on a schema mismatch.
"""
from __future__ import annotations

import dataclasses
from typing import Optional

import requests

GRAPHQL_URL = "https://pharos-api.ncats.io/graphql"
DIFFICULT_TIERS = {"TDARK", "TBIO"}


@dataclasses.dataclass
class PharosTarget:
    name: str
    sym: str
    uniprot: str
    tdl: str
    fam: str
    novelty: Optional[float] = None
    description: Optional[str] = None

    @property
    def is_difficult(self) -> bool:
        """Tdark/Tbio -- the IDG 'understudied' tiers this whole roadmap exists to prioritize."""
        return self.tdl.upper() in DIFFICULT_TIERS

    def to_tcrd_csv_row(self) -> str:
        """One row in the schema `target_druggability`'s --tcrd loader expects:
        symbol,uniprot,name,tdl,family,novelty_score
        """
        novelty = "" if self.novelty is None else str(self.novelty)
        name = self.name.replace(",", ";")
        return f"{self.sym},{self.uniprot},{name},{self.tdl},{self.fam},{novelty}"


TCRD_CSV_HEADER = "symbol,uniprot,name,tdl,family,novelty_score"


def _post(query: str, variables: dict, timeout: float = 30.0) -> dict:
    resp = requests.post(GRAPHQL_URL, json={"query": query, "variables": variables}, timeout=timeout)
    resp.raise_for_status()
    payload = resp.json()
    if "errors" in payload and payload["errors"]:
        raise RuntimeError(f"Pharos GraphQL error: {payload['errors']}")
    return payload["data"]


_TARGET_DETAIL_QUERY = """
query GetTarget($q: ITarget!) {
  target(q: $q) {
    name
    sym
    uniprot
    tdl
    fam
    novelty
    description
  }
}
"""


def get_target(symbol: Optional[str] = None, uniprot: Optional[str] = None) -> Optional[PharosTarget]:
    """Look up one target by gene symbol or UniProt accession. Returns
    None if Pharos has no record (this happens for real -- not every
    UniProt accession has a TCRD/Pharos entry)."""
    if not symbol and not uniprot:
        raise ValueError("provide symbol or uniprot")
    q = {"uniprot": uniprot} if uniprot else {"sym": symbol}
    data = _post(_TARGET_DETAIL_QUERY, {"q": q})
    t = data.get("target")
    if not t:
        return None
    return PharosTarget(
        name=t["name"], sym=t["sym"], uniprot=t["uniprot"], tdl=t["tdl"],
        fam=t.get("fam") or "Other", novelty=t.get("novelty"), description=t.get("description"),
    )


_DISEASE_TARGETS_QUERY = """
query DiseaseTargets($disease: String!, $top: Int!) {
  targets(filter: { associatedDisease: $disease }) {
    targets(top: $top) {
      name
      sym
      uniprot
      tdl
      fam
      novelty
    }
  }
}
"""


def discover_targets_for_disease(disease: str, top: int = 25) -> list[PharosTarget]:
    """Batch discovery: targets associated with a disease, for the
    'Pharos GraphQL API -> ...' pipeline entry point in the roadmap.
    Ranked by however Pharos' own disease-association relevance ranking
    orders them; callers filter down to `.is_difficult` for Tdark/Tbio
    triage.
    """
    data = _post(_DISEASE_TARGETS_QUERY, {"disease": disease, "top": top})
    raw = data.get("targets", {}).get("targets", []) or []
    out = []
    for t in raw:
        if not t.get("uniprot"):
            continue  # can't proceed to AlphaFold without a UniProt accession
        out.append(PharosTarget(
            name=t["name"], sym=t["sym"], uniprot=t["uniprot"], tdl=t.get("tdl") or "Tdark",
            fam=t.get("fam") or "Other", novelty=t.get("novelty"),
        ))
    return out


if __name__ == "__main__":
    import argparse
    import sys

    ap = argparse.ArgumentParser(description="Query Pharos for one target or a disease's target list.")
    ap.add_argument("--symbol")
    ap.add_argument("--uniprot")
    ap.add_argument("--disease")
    ap.add_argument("--top", type=int, default=25)
    args = ap.parse_args()

    if args.disease:
        for t in discover_targets_for_disease(args.disease, args.top):
            flag = "  [DIFFICULT]" if t.is_difficult else ""
            print(f"{t.sym}\t{t.uniprot}\t{t.tdl}\t{t.fam}{flag}")
    elif args.symbol or args.uniprot:
        t = get_target(symbol=args.symbol, uniprot=args.uniprot)
        if t is None:
            print("no Pharos record found", file=sys.stderr)
            sys.exit(1)
        print(t)
    else:
        ap.print_help()
        sys.exit(2)
