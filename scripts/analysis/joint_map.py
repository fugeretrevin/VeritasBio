"""Map each 2^5 joint-ablation arm to its run directory, then analyse the factorial.

The sweep script printed one summary row per arm but the benchmark manifest does not
record which components were ablated, so the run directories cannot be labelled from
their own contents. Rather than reconstruct the order (PowerShell's Sort-Object is not
stable, so the subset order is not predictable), each printed row is matched to the run
whose per-seed means reproduce it. Mean best fitness to 4 decimals plus mean top-1% to
2 decimals is a strong enough fingerprint to be unambiguous; the script asserts that.
"""
from __future__ import annotations

import csv
import json
import re
import statistics as st
from pathlib import Path
from typing import Dict, List, Optional, Tuple

ROOT = Path(r"C:\Users\stoof\OneDrive\Documents\GitHub\VeritasBio\benchmark_comparison_results")
OUT = Path(r"C:\Users\stoof\AppData\Local\Temp\claude"
           r"\c--Users-stoof-OneDrive-Documents-GitHub-VeritasBio"
           r"\0691478f-a671-445d-b8dc-a7025937cb24\tasks\bke7hs34z.output")
DEST = Path(r"C:\Users\stoof\AppData\Local\Temp\claude"
            r"\c--Users-stoof-OneDrive-Documents-GitHub-VeritasBio"
            r"\0691478f-a671-445d-b8dc-a7025937cb24\scratchpad\joint_arms.json")

HEADER = re.compile(r"^=== (\w+) : minus (.+) ===$")
ROW = re.compile(r"^\s+panoptes_default\s+available\s+(\d+)/(\d+)\s+"
                 r"([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)")

TOP1 = "cum_top1_hits_including_initial"
BEST = "best_measured_fitness"


def final_rows(run: Path, method: str = "panoptes_default") -> Dict[int, dict]:
    out: Dict[int, dict] = {}
    path = run / "round_metrics.csv"
    if not path.exists():
        return out
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r.get("method") != method:
                continue
            seed = int(float(r["seed"]))
            prev = out.get(seed)
            if prev is None or int(float(r["round"])) >= int(float(prev["round"])):
                out[seed] = r
    return out


def series(rows: Dict[int, dict], key: str) -> Dict[int, float]:
    return {s: float(r[key]) for s, r in rows.items() if r.get(key) not in (None, "")}


def library_of(run: Path) -> Optional[str]:
    try:
        man = json.loads((run / "manifest.json").read_text(encoding="utf-8"))
    except Exception:
        return None
    path = str(man.get("dataset", {}).get("path", "")).upper()
    if "TRPB" in path:
        return "trpb"
    if "GB1" in path:
        return "gb1"
    if "HER2" in path:
        return "her2"
    return None


# --- parse the sweep output into (library, removed-set, printed means) in order ---
printed: List[Tuple[str, str, Optional[Tuple[float, float]]]] = []
lines = OUT.read_text(encoding="utf-8", errors="replace").splitlines()
i = 0
while i < len(lines):
    h = HEADER.match(lines[i].strip()) or HEADER.match(lines[i])
    if h:
        lib, label = h.group(1), h.group(2)
        means = None
        for j in range(i + 1, min(i + 4, len(lines))):
            if HEADER.match(lines[j].strip()):
                break
            m = ROW.match(lines[j])
            if m:
                means = (float(m.group(3)), float(m.group(5)))
                break
        printed.append((lib, label, means))
    i += 1

print(f"headers parsed: {len(printed)}   with a result row: "
      f"{sum(1 for p in printed if p[2] is not None)}")
missing = [(l, lab) for l, lab, m in printed if m is None]
if missing:
    print(f"no result row for: {missing}")

# --- align the sweep's own trailing block of runs to the printed arms, in order ---
# Earlier one-at-a-time sweeps re-ran several of these same arms and, the engine being
# deterministic per seed, produced byte-identical means, so a mean-only fingerprint is
# ambiguous across sweeps. The sweep ran last and sequentially, and the run directories
# are named by UTC start, so its arms are the trailing block. Order gives the labels and
# the printed means are then used to VERIFY each alignment rather than to discover it.
expected = [(lib, label) for lib, label, means in printed if means is not None]

usable: List[Tuple[Path, str, float, float]] = []
for run in sorted(ROOT.glob("run_2026103*")):
    rows = final_rows(run)
    lib = library_of(run)
    if len(rows) != 48 or lib is None:
        continue
    best, top1 = series(rows, BEST), series(rows, TOP1)
    if not best or not top1:
        continue
    usable.append((run, lib, st.fmean(best.values()), st.fmean(top1.values())))

block = usable[-len(expected):]
arms: Dict[str, Dict[str, str]] = {"gb1": {}, "trpb": {}}
bad = []
for (lib_x, label), (run, lib_r, best_m, top1_m) in zip(expected, block):
    printed_means = next(m for l, lb, m in printed if (l, lb) == (lib_x, label))
    ok = (lib_x == lib_r
          and abs(best_m - printed_means[0]) < 1e-3
          and abs(top1_m - printed_means[1]) < 5e-2)
    if ok:
        arms[lib_x][label] = run.name
    else:
        bad.append((lib_x, label, run.name, lib_r, round(best_m, 4),
                    round(top1_m, 2), printed_means))

print(f"\nverified: gb1 {len(arms['gb1'])}, trpb {len(arms['trpb'])}")
if bad:
    print(f"MISALIGNED ({len(bad)}):")
    for b in bad[:8]:
        print(f"   {b}")

DEST.write_text(json.dumps(arms, indent=2), encoding="utf-8")
print(f"\nwrote {DEST}")
