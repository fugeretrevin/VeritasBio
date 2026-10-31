"""Is CR9114's 25% missing data depth-correlated?

Phillips et al. report 65,091 measured CR9114 variants against H1. Our prepared file
holds 48,840 - a 25% attrition the paper does not account for. That matters because
CR9114 is one of the three proteins in the pre-registered pooled aggregate, and the
quantity we claim about it is a PER-DEPTH statistic (shallow-shell enrichment, 0.622x).

If the missing rows are spread evenly across mutation depth, the shell statistics are
unbiased and only the counts shrink. If the attrition is concentrated at particular
depths, then every per-depth number we report for CR9114 is biased by the filter rather
than by biology, and the H1 ranking that depends on it is unsafe.

The library is 16 binary sites, so the complete design space is 2^16 = 65,536 and the
expected count at depth k is exactly C(16, k). Comparing observed counts against that
binomial is a direct test: flat retention across depths means the filter is depth-blind.
"""
from __future__ import annotations

import csv
import json
from math import comb
from pathlib import Path
from typing import Dict, Optional, Tuple

DATA = Path(r"C:\Users\stoof\OneDrive\Documents\GitHub\VeritasBio\scripts\data")
FILES = {
    "CR9114_H1": ("CR9114_H1_Phillips_2021.csv", 16, 65_091),
    "CR6261_H1": ("CR6261_H1_Phillips_2021.csv", 11, 1_812),
    "CR6261_H9": ("CR6261_H9_Phillips_2021.csv", 11, 1_812),
}

wt_seqs = json.loads((DATA / "library_wt_sequences.json").read_text(encoding="utf-8"))


def depth_of(mutant: str) -> Optional[int]:
    """Count substitutions. Binary sites, so depth is just the token count."""
    s = mutant.strip()
    if not s:
        return None
    if s.upper() == "WT":
        return 0
    toks = [t for t in s.split(":") if t]
    for t in toks:
        digits = "".join(c for c in t if c.isdigit())
        if not digits or len(t) < 3:
            return None
    return len(toks)


print("=" * 96)
print("IS THE MISSING DATA DEPTH-CORRELATED?")
print("=" * 96)

for label, (fname, n_sites, paper_n) in FILES.items():
    path = DATA / fname
    if not path.exists():
        print(f"  ! missing {fname}")
        continue
    obs: Dict[int, int] = {}
    floored = 0
    total = 0
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            d = depth_of(str(r.get("mutant", "")))
            if d is None:
                continue
            obs[d] = obs.get(d, 0) + 1
            total += 1
            try:
                if abs(float(r["DMS_score"]) - 7.0) < 1e-9:
                    floored += 1
            except (TypeError, ValueError, KeyError):
                pass

    space = 2 ** n_sites
    print()
    print(f"  {label}   {n_sites} binary sites, design space {space:,}")
    print(f"    our rows {total:,}   paper reports {paper_n:,}   "
          f"complete space {space:,}   at floor 7.0: {floored:,}")
    print(f"    {'depth':<7}{'expected':>10}{'observed':>10}{'retention':>11}")
    retentions = []
    for d in range(n_sites + 1):
        exp = comb(n_sites, d)
        got = obs.get(d, 0)
        ret = got / exp if exp else float("nan")
        retentions.append((d, exp, got, ret))
        flag = ""
        if exp >= 10:
            if ret < 0.5:
                flag = "  <-- heavily depleted"
            elif ret > 0.995:
                flag = "  <-- complete"
        print(f"    {d:<7}{exp:>10,}{got:>10,}{ret:>10.1%}{flag}")

    big = [r for r in retentions if r[1] >= 10]
    if big:
        lo = min(big, key=lambda r: r[3])
        hi = max(big, key=lambda r: r[3])
        print(f"    retention range over depths with >=10 expected: "
              f"{lo[3]:.1%} (depth {lo[0]}) to {hi[3]:.1%} (depth {hi[0]})")
        spread = hi[3] - lo[3]
        verdict = ("DEPTH-CORRELATED - per-depth statistics are biased by the filter"
                   if spread > 0.15 else
                   "approximately depth-blind - shell statistics unbiased, counts only shrink")
        print(f"    => {verdict}")
