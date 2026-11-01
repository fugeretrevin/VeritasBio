"""The Day-1 dataset facts the expert asked for, on CR9114-H1 specifically.

Answers: full depth-by-depth hit table, depth of the global optimum, the rank of the
fully-mutated (matured) reference - which decides whether seeding a reverse search from
it is leakage - and the expected depth profile of a uniform 25-variant initial sample.
"""
import csv, collections, statistics
from pathlib import Path

p = Path("scripts/data/CR9114_H1_Phillips_2021.csv")
rows = []
for r in csv.DictReader(open(p)):
    try: s = float(r["DMS_score"])
    except (TypeError, ValueError): continue
    subs = [x for x in r["mutant"].split(":") if x]
    rows.append((len(subs), s, r["mutant"]))

n = len(rows)
by_score = sorted(rows, key=lambda t: -t[1])
k = max(1, round(n * 0.01))
thresh = by_score[k - 1][1]
hits = [t for t in rows if t[1] >= thresh]

print(f"CR9114-H1: {n} variants, top-1% = {len(hits)} (threshold score {thresh:.4f})")
print(f"assay floor count (score==7.0): {sum(1 for d,s,_ in rows if s == 7.0)}")
print()

print("FULL depth table (all depths 0-16):")
print(f"{'depth':>5} {'variants':>9} {'top1 hits':>10} {'hit rate':>9} {'enrichment':>11} {'max score':>10}")
allc = collections.Counter(d for d, _, _ in rows)
hitc = collections.Counter(d for d, _, _ in hits)
maxs = collections.defaultdict(lambda: float("-inf"))
for d, s, _ in rows: maxs[d] = max(maxs[d], s)
base = len(hits) / n
cum_reach = 0
for d in range(0, 17):
    if allc[d] == 0 and hitc[d] == 0: continue
    rate = hitc[d] / allc[d] if allc[d] else 0.0
    print(f"{d:>5} {allc[d]:>9} {hitc[d]:>10} {100*rate:>8.3f}% "
          f"{rate/base if base else 0:>10.2f}x {maxs[d]:>10.4f}")
print()

gd, gs, gm = by_score[0]
print(f"GLOBAL OPTIMUM: depth {gd}, score {gs:.4f}")
print(f"  mutant: {gm}")
print()

# The fully-mutated reference (all 16 substitutions) - the matured antibody.
full = [t for t in rows if t[0] == max(allc)]
deepest = max(allc)
print(f"deepest depth present in table: {deepest} ({allc[deepest]} variant(s))")
for d, s, m in sorted(full, key=lambda t: -t[1])[:3]:
    rank = next(i for i, t in enumerate(by_score, 1) if t[2] == m)
    pct = 100.0 * rank / n
    print(f"  depth-{d} reference score {s:.4f} -> rank {rank}/{n} (top {pct:.3f}%)"
          f"  {'IN top-1%' if s >= thresh else 'NOT in top-1%'}")
print()

CEIL = 4
print(f"Reachable at the generator's measured ceiling (depth <= {CEIL}):")
r_all = sum(allc[d] for d in range(CEIL + 1))
r_hit = sum(hitc[d] for d in range(CEIL + 1))
print(f"  {r_all} of {n} variants ({100*r_all/n:.2f}%)")
print(f"  {r_hit} of {len(hits)} top-1% hits ({100*r_hit/len(hits):.2f}%)")
print(f"  best score reachable: {max(maxs[d] for d in range(CEIL+1)):.4f} "
      f"vs global {gs:.4f}")
print()

# Expected initial-sample depth profile: uniform 25 from the table.
print("Uniform 25-variant initial sample, expected depth profile:")
exp = {d: 25.0 * allc[d] / n for d in sorted(allc)}
print("  " + "  ".join(f"d{d}:{v:.2f}" for d, v in exp.items() if v >= 0.01))
p_deep = 1.0 - (1.0 - sum(allc[d] for d in range(9, 17)) / n) ** 25
p_hit = 1.0 - (1.0 - len(hits) / n) ** 25
print(f"  P(at least one initial variant at depth>=9) = {100*p_deep:.3f}%")
print(f"  P(at least one initial variant is top-1%)   = {100*p_hit:.3f}%")
