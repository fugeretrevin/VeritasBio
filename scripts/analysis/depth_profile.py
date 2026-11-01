"""Per-round depth profile of what the generator PROPOSES vs what it SELECTS.

The question is whether the engine's reachable depth intersects the depth shells where
the top-1% hits actually live. depth_returns.csv already records this per seed, per
round, per depth; nothing here instruments anything, it only aggregates.
"""
import csv, sys, collections
from pathlib import Path

def profile(run_dir: Path, label: str) -> None:
    rows = list(csv.DictReader(open(run_dir / "depth_returns.csv")))
    if not rows:
        print(f"{label}: no depth rows"); return
    seeds = {r["seed"] for r in rows}
    # per (round, depth): totals across seeds
    agg = collections.defaultdict(lambda: [0, 0, 0])
    maxprop = collections.defaultdict(set)   # round -> set of max proposed depth per seed
    maxsel = collections.defaultdict(set)
    per_seed_round = collections.defaultdict(lambda: ([], []))
    for r in rows:
        rd, d = int(r["round_index"]), int(r["depth"])
        p, s, h = int(r["proposed"]), int(r["selected"]), int(r["top1_hits"])
        a = agg[(rd, d)]; a[0] += p; a[1] += s; a[2] += h
        if p > 0: per_seed_round[(r["seed"], rd)][0].append(d)
        if s > 0: per_seed_round[(r["seed"], rd)][1].append(d)
    for (seed, rd), (props, sels) in per_seed_round.items():
        if props: maxprop[rd].add((seed, max(props)))
        if sels: maxsel[rd].add((seed, max(sels)))

    print(f"\n{'='*78}\n{label}   ({len(seeds)} seeds)\n{'='*78}")
    print(f"{'round':>5} {'depth':>6} {'proposed':>10} {'selected':>9} {'top1 hits':>10} {'hit rate':>9}")
    for (rd, d) in sorted(agg):
        p, s, h = agg[(rd, d)]
        rate = f"{100.0*h/s:7.1f}%" if s else "       -"
        print(f"{rd:>5} {d:>6} {p:>10} {s:>9} {h:>10} {rate:>9}")

    print(f"\n  max depth REACHED per round (mean over seeds):")
    for rd in sorted(maxprop):
        pv = [v for _, v in maxprop[rd]]
        sv = [v for _, v in maxsel.get(rd, [])]
        mp = sum(pv)/len(pv) if pv else 0
        ms = sum(sv)/len(sv) if sv else 0
        print(f"    round {rd}: proposed max {mp:.2f} (abs {max(pv) if pv else 0})"
              f"   selected max {ms:.2f} (abs {max(sv) if sv else 0})")

for arg in sys.argv[1:]:
    d, label = arg.split("=", 1)
    profile(Path(d), label)
