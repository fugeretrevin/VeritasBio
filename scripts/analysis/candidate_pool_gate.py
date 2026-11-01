"""Why does the engine find 0 of 82 top-1% hits on eqFP611-red while ridge finds 15.35?

It spent 80.3 of 81 measurements, so this is not the CR9114 starvation failure. Two
mutually exclusive causes, needing opposite fixes:

  GENERATION  the hits never enter the candidate pool  -> the geometry fix is incomplete
  RANKING     they enter and are ranked down           -> prior / acquisition / Pareto

`candidate_trace.csv` is the engine's real ballot and carries, per candidate, the prior
contribution, the acquisition score, the Pareto rank and the final rank. So the question
is answerable from disk. Read-only.
"""
import csv
import collections
import json
import statistics
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
import benchmark_panoptes_isolated as bench  # noqa: E402

if len(sys.argv) < 3:
    raise SystemExit(
        "usage: candidate_pool_gate.py <run-dir> <library-key>\n"
        f"  library keys: {', '.join(sorted(bench.LIBRARIES))}"
    )
RUN = Path(sys.argv[1])
LIBRARY = sys.argv[2]
if LIBRARY not in bench.LIBRARIES:
    raise SystemExit(
        f"unknown library {LIBRARY!r}; known: "
        f"{', '.join(sorted(bench.LIBRARIES))}"
    )
DATASET = bench.LIBRARIES[LIBRARY].dataset

manifest = json.load(open(RUN / "manifest.json"))
threshold = manifest["results"]["top1"]["threshold_value"]
print(f"top-1% threshold from the run's own manifest: {threshold}")

truth = {}
for row in csv.DictReader(open(DATASET)):
    try:
        truth[row["mutant"].strip()] = float(row["DMS_score"])
    except (TypeError, ValueError):
        continue
hits = {v for v, s in truth.items() if s >= threshold}
print(f"universe {len(truth)}, top-1% hits {len(hits)}")

rows = [r for r in csv.DictReader(open(RUN / "candidate_trace.csv"))]
print(f"ballot rows: {len(rows)}")
print()

# ---- 1. Generation: did any hit ever get proposed? --------------------------------
proposed = {r["variant_id"] for r in rows}
proposed_hits = proposed & hits
print("=" * 70)
print("1. GENERATION")
print(f"   distinct candidates ever proposed : {len(proposed)}")
print(f"   of which are top-1% hits          : {len(proposed_hits)} of {len(hits)} "
      f"({100.0 * len(proposed_hits) / len(hits):.1f}%)")
by_source = collections.Counter(r["generation_source"] for r in rows)
print(f"   generation_source                 : {dict(by_source)}")
print()

if not proposed_hits:
    print("   VERDICT: GENERATION. No hit ever entered the pool, so no amount of")
    print("   re-ranking could have selected one. The geometry fix is incomplete here.")
    raise SystemExit

# ---- 2. Ranking: where did the proposed hits land? ---------------------------------
def nums(values):
    out = []
    for v in values:
        try:
            out.append(float(v))
        except (TypeError, ValueError):
            pass
    return out

hit_rows = [r for r in rows if r["variant_id"] in hits]
miss_rows = [r for r in rows if r["variant_id"] not in hits]

print("=" * 70)
print("2. RANKING - hits that WERE proposed, versus everything else")
print(f"   hit ballot rows {len(hit_rows)}, non-hit ballot rows {len(miss_rows)}")
print()
print(f"   {'column':<24}{'hits (mean)':>14}{'non-hits (mean)':>17}{'  signal'}")
print("   " + "-" * 62)
for column in ("predicted_fitness", "posterior_uncertainty", "prior_contribution",
               "prior_structural", "prior_raw_total", "acquisition_score",
               "pareto_rank", "quality_score", "diversity_adjustment", "final_rank"):
    h, m = nums(r[column] for r in hit_rows), nums(r[column] for r in miss_rows)
    if not h or not m:
        continue
    hm, mm = statistics.fmean(h), statistics.fmean(m)
    # For ranks, lower is better; for scores, higher is better.
    lower_better = "rank" in column
    better = (hm < mm) if lower_better else (hm > mm)
    note = "hits favoured" if better else "** hits PENALISED **"
    print(f"   {column:<24}{hm:>14.5f}{mm:>17.5f}  {note}")

print()
status = collections.Counter(r["status"] for r in hit_rows)
print(f"   status of hit rows: {dict(status)}")
selected = [r for r in hit_rows if r["status"] not in ("rejected",)]
print(f"   hit rows not rejected: {len(selected)}")

# How deep in the ranked list did hits sit? If the engine proposed a hit but ranked it
# 900th, the batch of 16 never reached it.
ranks = nums(r["final_rank"] for r in hit_rows)
if ranks:
    ranks.sort()
    print(f"   hit final_rank: min {ranks[0]:.0f}, median {statistics.median(ranks):.0f}, "
          f"max {ranks[-1]:.0f}")
    print(f"   hits ranked inside the top 16 of their round: "
          f"{sum(1 for r in ranks if r <= 16)} of {len(ranks)}")
print()
print("   VERDICT: RANKING. The hits were proposed. Whichever column above shows")
print("   'hits PENALISED' is where selection is going wrong, and each of prior,")
print("   acquisition and pareto is separately ablatable.")
