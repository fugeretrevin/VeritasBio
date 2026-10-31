# Benchmark analysis scripts

Every number in the benchmark review package is derived by one of these scripts from the
runs under `benchmark_comparison_results/`. They are kept so the findings — including the
negative ones and the two analysis errors that were caught and corrected — are
reproducible rather than asserted.

Each script writes a `*_report.txt` beside itself; those outputs are committed too, so the
results can be read without re-running anything.

## Running them

Two constants at the top of each script are absolute paths from the machine the analysis
was done on:

- `ROOT` — the `benchmark_comparison_results/` directory
- `SCRATCH` — where the report is written

Point both somewhere local before running. Then:

```
.venv\Scripts\python.exe -X utf8 -I scripts\analysis\<script>.py
```

`-X utf8` matters: the reports contain non-ASCII characters and the Windows console
default (cp1252) will raise `UnicodeEncodeError` without it. `-I` is used because these
read data files, and isolated mode keeps the script's directory off `sys.path`.

Several scripts hard-code the run ids of specific arms (held-out seeds 48–95). Those ids
are listed in `joint_arms.json` and in the `HELDOUT` / `ARMS` tables inside each script.

## What each one does

| script | question it answers |
| --- | --- |
| `joint_map.py` | which of the 62 completed runs corresponds to which configuration |
| `joint_analyse.py` | the 2^5 component factorial: main effects and two-way interactions |
| `joint_choice.py` | ranks configurations on the development seeds and picks one to confirm |
| `heldout.py` | confirms the selected configuration on held-out seeds 48–95 |
| `clustered.py` | **re-does the factorial with seed-clustered inference** (see below) |
| `budget.py` | performance at 25/50/75/100/125 measurements, from existing runs |
| `diversity.py` | whether the advantage is distinct solutions or one dense basin |
| `responsive.py` | whether the degenerate TrpB configuration ignores its measurements |
| `reachable.py` | what the degenerate configuration's candidate generator can actually reach |
| `gb1_mechanism.py` | why reliability wins on TrpB and loses on GB1 |

## Two corrections worth knowing about

These are the reason the scripts are committed rather than discarded.

**Seed-clustered inference (`clustered.py` supersedes part of `joint_analyse.py`).** The
original factorial pooled 768 background-pair differences and ran a t-test over them as if
independent, when each of the 48 seeds contributed 16 of them. That understates the
standard error. `clustered.py` redoes it with the seed as the unit of independence —
per-seed contrasts, a paired t over seeds, and a bootstrap resampling whole seeds.
t-statistics were inflated 1.3x–2.6x. No headline conclusion reversed, but three marginal
calls changed. `joint_report.txt` is retained as the archived original.

**Genotype canonicalisation (`diversity.py`).** Mutation strings omit wild-type positions,
so `V183A:F184A:V227A` is a four-site genotype whose fourth site is wild-type. An early
version compared such strings directly, and `zip()` silently truncated the 4-tuple to 3
and understated Hamming distances — which produced the *opposite* conclusion about
diversity. The committed version fills wild-type at omitted positions and asserts equal
length rather than truncating. It validates by recovering exactly 1,600 of 1,600 TrpB and
1,494 of 1,494 GB1 top-1% library members, matching each manifest's independently recorded
count.

Any future script comparing genotypes must canonicalise the same way. Never zip or
truncate genotypes of unequal length.

## A note on thresholds

Top-1% thresholds are read from each run's own `manifest.json` (`top1.threshold`, the 99th
percentile over the eligible set) so that analysis scores variants with exactly the
threshold the benchmark used. Do not recompute them independently — the eligible universe
depends on library-specific eligibility rules, and HER2 in particular excludes an
unlabelled split whose rows carry `DMS_score` of exactly 0.0.
