# Evidence — the run artifacts behind every number we quote

`benchmark_comparison_results/` is gitignored, because a benchmark that writes hundreds of
run directories should not push them. The consequence is that a public clone carries no
artifacts at all and every figure becomes unverifiable without re-running. This directory
is the fix: the runs that the README and the submission actually cite, copied verbatim.

Nothing here is regenerated or reformatted. Each subdirectory keeps its original run id,
so it matches the id quoted in the text, and holds the five files that let a reader check
a claim without a GPU or an API key:

| file | what it is for |
| --- | --- |
| `manifest.json` | the authoritative record: library, budget, seed block, component states, commit, dirty-tree flag, environment, label provenance, and the configuration fingerprint. Contains no outcomes, deliberately. |
| `summary.csv` | one row per arm: status, seeds completed, mean measurements used, best fitness, top-1% yield (including and excluding the initial sample), token usage |
| `round_metrics.csv` | per arm, per seed, per round. This is the file to recompute paired statistics from — seed is the unit of independence |
| `report.txt` | the console transcript, including the limitations and cost sections |
| `benchmark.log` | per-seed progress and any error, verbatim |

## Which run backs which claim

| claim | run id | library |
| --- | --- | --- |
| Live Meta agent beats one-hot ridge on top-1% yield; the 10x variance result; 568,189 tokens over 96 calls | `run_20261101_043406_8662bdf0` | trpb |
| Engine beats ridge on 12 of 12 seeds; one seed reaches the global optimum 9.83451; ridge not separable from random | `run_20261101_051322_4bb22522` | cr9114_h1 |
| Six-landscape sweep — the three losses, the tie, and the development-data win | the six runs below | — |
| The content-policy refusal, verbatim 400 body | `run_20261101_052018_8755fea6` | cr9114_h1 |

The sweep, all at commit `406ef58`, 48 paired seeds each:

| run id | library |
| --- | --- |
| `run_20261101_024225_8fa1d86d` | eqfp611_red |
| `run_20261101_024403_fe1ae236` | eqfp611_blue |
| `run_20261101_024544_e82e6214` | creilov |
| `run_20261101_025014_5e79e3ae` | gb1 |
| `run_20261101_025156_e20ba927` | cr9114_h1 |
| `run_20261101_025436_dfa910c3` | her2 |

## Three traps in this data, named rather than left to be discovered

**1. `run_20261101_051322_4bb22522` has run state `PARTIAL`, and its `panoptes_muse` rows
are not a result.** The live arm was requested and failed on all 12 seeds with an
`OrchestratorConfig` validation error (`target_accession=None`, fixed in `c97b4fb`). Its
rows are still in `round_metrics.csv` with `status=completed` **for round 0 only** — and
round 0 is the shared random starting sample, identical across arms by construction. Any
aggregation that does not filter on completed campaigns will produce a plausible-looking
`panoptes_muse` figure for CR9114 that is just the random start. The engine figures we
quote for this landscape come from `panoptes_default`, which completed 12 of 12.

**2. `her2` must not be quoted from the sweep.** Its pre-flight gate FAILED on
`endpoint.her2`: the d>=2 diversity ceiling equals the raw hit count, so 0% of the hit set
is removed against a 5% informativeness floor. That is a property of the library and no
budget change fixes it. The run is included for completeness and because the failing gate
is itself worth seeing.

**3. `panoptes_default` changed meaning on 2026-10-31.** `DPP_DIVERSITY_ENABLED` and
`TWO_STAGE_GP_ENABLED` both default `False` from that date, because both were measured to
hurt top-1% yield. So a bare arm name is not a sufficient description of a configuration:
read `configuration.code.commit` from the manifest, not the arm label. Note also that
manifests written before `e083333` report `dpp` and `dkl` as *enabled* regardless, because
`component_states()` read only the CLI flag and could not see the module defaults. The
commit is authoritative; the component field, in those runs, is not.

## Reproducing a figure

Paired statistics, from any run's `round_metrics.csv`. Take the LAST round per
(method, seed) — that is the campaign's final cumulative count — then difference the two
arms within each seed:

```python
import csv, statistics
from collections import defaultdict

final = defaultdict(dict)          # method -> seed -> value
with open("evidence/run_20261101_043406_8662bdf0/round_metrics.csv",
          encoding="utf-8", newline="") as fh:
    for row in csv.DictReader(fh):
        m, s, r = row["method"], int(row["seed"]), int(row["round"])
        prev = final[m].get(s, (-1, None))
        if r >= prev[0]:
            final[m][s] = (r, float(row["cum_top1_hits_excluding_initial"]))

a = {s: v for s, (_r, v) in final["panoptes_muse"].items()}
b = {s: v for s, (_r, v) in final["ridge"].items()}
seeds = sorted(set(a) & set(b))
diffs = [a[s] - b[s] for s in seeds]
print(len(seeds), statistics.fmean(diffs),
      sum(1 for d in diffs if d > 0), sum(1 for d in diffs if d < 0))
```

Two rules that are load-bearing, both learned the hard way and recorded in
`docs/HANDOFF.md` §2:

- **Seed is the unit of independence.** Do not pool per-round rows as independent
  observations; it inflated t-statistics here by 1.3-2.6x once.
- **Check `mean_measurements_used` is equal across arms before comparing yields.** Engine
  arms once measured 120-124 against baselines' 125, which voided a set of headline
  numbers.
