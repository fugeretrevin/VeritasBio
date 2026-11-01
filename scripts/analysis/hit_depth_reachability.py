"""Where do the top-1% hits live, by mutation depth, versus the generator's depth-4 ceiling?

Mutation strings omit wild-type positions (HANDOFF rule), so depth is the count of
listed substitutions. Asserted below rather than assumed: no listed substitution may be
a synonymous wt->wt, and max depth may not exceed the library's position count.
"""
import csv, collections, sys
from pathlib import Path

DATA = Path("scripts/data")
CEILING = 4

FILES = [
    ("TrpB",         "TRPB_THEMA_Johnston_2024.csv",   "functional"),
    ("GB1",          "GB1_STRSG_Wu_2016.csv",          "functional"),
    ("CreiLOV",      "CREILOV_CHLRE_Chen_2023.csv",    "functional"),
    ("eqFP611-blue", "EQFP611_BLUE_Poelwijk_2019.csv", "functional"),
    ("eqFP611-red",  "EQFP611_RED_Poelwijk_2019.csv",  "NON-functional"),
    ("CR9114-H1",    "CR9114_H1_Phillips_2021.csv",    "NON-functional"),
    ("CR6261-H1",    "CR6261_H1_Phillips_2021.csv",    "NON-functional"),
    ("CR6261-H9",    "CR6261_H9_Phillips_2021.csv",    "NON-functional"),
    ("HER2",         "HER2_TRASTUZUMAB_Mason_2021.csv", "functional"),
]

print(f"{'landscape':<14} {'parent':<15} {'N':>7} {'pos':>4} {'top1%':>6} "
      f"{'reachable(<=4)':>14} {'hits<=4':>8} {'hits>4':>7} {'% reachable':>12}")
print("-" * 104)

for label, fname, parent in FILES:
    p = DATA / fname
    if not p.exists():
        print(f"{label:<14} MISSING {p}"); continue
    rows = []
    positions = set()
    for r in csv.DictReader(open(p)):
        mut = r["mutant"].strip()
        try:
            score = float(r["DMS_score"])
        except (TypeError, ValueError):
            continue
        subs = [s for s in mut.split(":") if s]
        # Guard the depth definition: reject a synonymous entry, which would mean
        # wild-type positions are NOT omitted and depth != len(subs).
        for s in subs:
            if len(s) >= 3 and s[0] == s[-1]:
                raise SystemExit(f"{label}: synonymous substitution {s!r} in {mut!r} - "
                                 "wild-type positions are not omitted; depth is not the count")
            if len(s) >= 3:
                positions.add(s[1:-1])
        rows.append((len(subs), score))
    if not rows:
        print(f"{label:<14} no usable rows"); continue

    n = len(rows)
    k = max(1, round(n * 0.01))
    thresh = sorted((s for _, s in rows), reverse=True)[k - 1]
    hits = [d for d, s in rows if s >= thresh]
    npos = len(positions)
    maxd = max(d for d, _ in rows)
    if maxd > npos:
        raise SystemExit(f"{label}: max depth {maxd} exceeds {npos} positions")

    reach = sum(1 for d, _ in rows if d <= CEILING)
    h_in = sum(1 for d in hits if d <= CEILING)
    h_out = len(hits) - h_in
    pct = 100.0 * h_in / len(hits) if hits else 0.0
    print(f"{label:<14} {parent:<15} {n:>7} {npos:>4} {len(hits):>6} "
          f"{reach:>14} {h_in:>8} {h_out:>7} {pct:>11.1f}%")

print()
print("Hit-depth histogram for the deep landscapes (top-1% only):")
for label, fname, parent in FILES:
    p = DATA / fname
    if not p.exists(): continue
    rows = []
    for r in csv.DictReader(open(p)):
        try: score = float(r["DMS_score"])
        except (TypeError, ValueError): continue
        rows.append((len([s for s in r["mutant"].split(":") if s]), score))
    if max(d for d, _ in rows) <= CEILING:
        continue
    n = len(rows); k = max(1, round(n * 0.01))
    thresh = sorted((s for _, s in rows), reverse=True)[k - 1]
    hist = collections.Counter(d for d, s in rows if s >= thresh)
    tot = sum(hist.values())
    bars = " ".join(f"d{d}:{hist[d]}" for d in sorted(hist))
    print(f"  {label:<14} ({parent:<14}) n_hits={tot:<5} {bars}")
